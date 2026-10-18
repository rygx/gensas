//
// Copyright (C) 2022, UChicago Argonne, LLC. All rights reserved.
//
// Software Name: Generic Semi-Analytical Simulation Tool (GenSAS)
// By: Argonne National Laboratory
// OPEN SOURCE LICENSE
//
// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
// 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
//
//
// ******************************************************************************************************
// DISCLAIMER
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
// PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// ***************************************************************************************************
//
#include "sas/HybridSimulator.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>

using namespace std;

namespace che {
	namespace core {

		// Added to the errors for events that pile up. The cause is in the model, so the message says what to do about it.
		static const char *const kZenoAdvice =
			": the model switches infinitely often, or without progress, in a finite time (for example something that loses a part of its energy at every impact, like a bouncing ball). "
			"Stop the run before this time, or let the model settle explicitly: in a when equation, end the run with terminate() or set the state to its resting value once it is small.";

		HybridSimulator::HybridSimulator(shared_ptr<SasModel> model) : model(model) {}

		bool HybridSimulator::fail(const string &msg) {
			err = msg;
			return false;
		}

		SasSolutionSet *HybridSimulator::run(const SasComputationOptions &options) {
			err.clear();
			wasCancelled = false;
			firedEvents.clear();
			finalVals.clear();
			warningList.clear();
			wasTerminated = false;
			terminateTime = 0.0;
			terminateMessage.clear();
			nSegments = 0;
			if (model == nullptr) {
				fail("no model");
				return NULL;
			}

			HybridEnv env = makeHybridEnv(*model);
			vector<Guard> guards = collectGuards(*model, env); // before the continuous variables enter env
			bool hasDynamic = false;
			for (const Guard &g : guards) {
				hasDynamic = hasDynamic || g.dynamic;
			}
			// Start values of the continuous variables, for event expressions at the start time.
			for (const shared_ptr<IdPearl> &id : model->idTable) {
				if ((id->idTy == ID_VAR || id->idTy == ID_IVAR) && id->str != NULL) {
					env.values[(const char *)id->str] = (id->initVal && id->ty != TY_INT && id->ty != TY_BOOL) ? id->value.d : 0.0;
				}
			}
			vector<WhenClause> clauses = collectWhenClauses(*model);
			unique_ptr<SasSolutionSet> total(new SasSolutionSet());
			map<string, double> state; // continuous variables by name, from the previous segment and reinits
			// The continuous variables at the end of the last segment of each mode (the values of the discrete variables),
			// to start the algebraic variables from when a mode is entered again (the voltages after a fault is cleared are
			// close to those before it, not to those during it).
			map<string, map<string, double>> modeSnapshot; // names of the algebraic variables to values
			bool allowStartHint = true;
			string previousModeKey;
			auto modeKey = [&]() {
				ostringstream os;
				os.precision(17);
				for (const string &n : env.discrete) {
					const auto it = env.values.find(n);
					os << n << '=' << (it != env.values.end() ? it->second : 0.0) << ';';
				}
				return os.str();
			};
			vector<string> stateNames; // names of the continuous variables, to detect structure changes
			map<int, bool> cur;		   // truth of the dynamic guards
			int closeEvents = 0;	   // consecutive events within 1e-8 of each other, to detect chattering
			double lastEventT = -1.0;  // times of the last two handled events (< 0: none), to size the
			double prevEventT = -1.0;  // detection sampling after an event

			const double endTol = 1e-9 * max(1.0, fabs(options.endTime));
			double t = options.startTime;
			env.time = t;
			if (!evalGuards(guards, env, cur, err)) {
				return NULL;
			}
			// An assertion that is already violated at the start time.
			auto reportAssertions = [&](double T, const map<int, bool> &now, const map<int, bool> &before) -> bool {
				for (const Guard &g : guards) {
					if (g.kind != G_ASSERT || !now.at(g.id) || (before.count(g.id) > 0 && before.at(g.id))) {
						continue;
					}
					ostringstream os;
					os << "assertion violated at t=" << T << ": " << g.message;
					if (!g.warning) {
						fail(os.str());
						return false;
					}
					warningList.push_back(os.str());
				}
				return true;
			};
			if (!reportAssertions(t, cur, map<int, bool>())) {
				return NULL;
			}

			auto record = [&](const HybridEvent &ev) {
				if (ev.assignments.empty() && ev.reinits.empty()) {
					return;
				}
				for (const pair<string, double> &r : ev.reinits) {
					state[r.first] = r.second;
					env.values[r.first] = r.second;
				}
				firedEvents.push_back(ev);
			};

			// Runs the events at time T: time events (if `timeEvents`), sample() events at the start time
			// (if `initial`) and the when-conditions that became true. Repeats while the event's own
			// changes make further conditions true.
			auto handleEvents = [&](double T, bool timeEvents, bool initial) -> bool {
				env.time = T;
				closeEvents = (lastEventT >= 0.0 && T - lastEventT < 1e-8 * (1.0 + fabs(T))) ? closeEvents + 1 : 0;
				if (closeEvents >= 10) {
					ostringstream os;
					os << "events accumulate at t=" << T << " (Zeno behavior or chattering)" << kZenoAdvice;
					fail(os.str());
					return false;
				}
				prevEventT = lastEventT;
				lastEventT = T;
				for (const pair<const string, double> &s : state) {
					env.values[s.first] = s.second;
				}
				// edge() and change() compare with the values before this event instant.
				env.eventPre.clear();
				for (const string &n : env.discrete) {
					env.eventPre[n] = env.values[n];
				}
				struct AtEvent {
					HybridEnv &e;
					explicit AtEvent(HybridEnv &env) : e(env) {
						e.atEvent = true;
					}
					~AtEvent() {
						e.atEvent = false;
					}
				} atEvent(env);
				set<int> fired;
				for (int iter = 0; iter < 50; iter++) {
					map<int, bool> now;
					if (!evalGuards(guards, env, now, err)) {
						return false;
					}
					fired.clear();
					for (const Guard &g : guards) {
						if (g.dynamic && g.kind == G_WHEN && now[g.id] && !cur[g.id]) {
							fired.insert(g.id);
						}
					}
					if (!reportAssertions(T, now, cur)) {
						return false;
					}
					cur = now;
					bool timeNow = (iter == 0) && (timeEvents || initial);
					if (fired.empty() && !timeNow) {
						return true;
					}
					HybridEvent ev;
					if (!fireEvents(clauses, env, T, initial && iter == 0, ev, err, &fired, !timeNow)) {
						return false;
					}
					record(ev);
					for (const string &w : ev.warnings) {
						warningList.push_back(w);
					}
					if (ev.terminate) {
						wasTerminated = true;
						terminateTime = T;
						terminateMessage = ev.message;
						return true;
					}
				}
				ostringstream os;
				os << "event iteration does not converge at t=" << T;
				fail(os.str());
				return false;
			};

			// Discrete variables over time, for the result.
			vector<string> discreteNames(env.discrete.begin(), env.discrete.end());
			total->discrete.names = discreteNames;
			auto snapshot = [&](double T) {
				vector<double> vals;
				for (const string &n : discreteNames) {
					vals.push_back(env.values[n]);
				}
				total->discrete.times.push_back(T);
				total->discrete.values.push_back(vals);
			};

			if (!handleEvents(t, false, true)) { // sample() and initial() events at the start time
				return NULL;
			}
			snapshot(t);
			if (wasTerminated) {
				fail("terminate() ran at the start time (" + terminateMessage + "); there is nothing to simulate");
				return NULL;
			}

			const size_t maxEvents = 1000000;
			int stuck = 0; // consecutive iterations that hardly advance the time
			while (true) {
				double te = 0.0;
				bool hasEvent = nextEventTime(clauses, env, t, te, err);
				if (!err.empty()) {
					return NULL;
				}
				bool timeEvent = hasEvent && te < options.endTime - endTol;
				double segEnd = timeEvent ? te : options.endTime;
				double dynT = -1.0; // time of a detected state event, if any

				if (segEnd - t > options.alphaTol) {
					shared_ptr<SasComputationModel> comp = make_shared<SasComputationModel>();
					comp->log = options.log;
					comp->generateCompModel(model);
					SasModel &cm = *comp->topModels.front();
					if (!specializeModel(cm, env, err, &cur)) {
						return NULL;
					}
					for (const shared_ptr<IdPearl> &id : cm.idTable) {
						if ((id->idTy == ID_VAR || id->idTy == ID_IVAR) && id->str != NULL) {
							map<string, double>::const_iterator it = state.find((const char *)id->str);
							if (it != state.end()) {
								id->value.d = it->second;
								id->initVal = true;
							}
						}
					}
					// Returning to a mode that was active before (a fault is cleared, a line is reconnected): the algebraic
					// variables start from the values at the end of the last time in that mode, not from the values of
					// the mode just left (the voltages are about 1 again, not close to zero as during the fault).
					const string segmentModeKey = modeKey();
					vector<shared_ptr<IdPearl>> hintedIds;
					bool hinted = false;
					if (allowStartHint && !previousModeKey.empty() && segmentModeKey != previousModeKey) {
						map<string, map<string, double>>::const_iterator snap = modeSnapshot.find(segmentModeKey);
						if (snap != modeSnapshot.end()) {
							for (const shared_ptr<IdPearl> &id : cm.idTable) {
								if ((id->idTy == ID_VAR || id->idTy == ID_IVAR) && id->str != NULL) {
									map<string, double>::const_iterator v = snap->second.find((const char *)id->str);
									if (v != snap->second.end()) {
										id->value.d = v->second;
										id->initVal = true;
										hintedIds.push_back(id);
										hinted = true;
									}
								}
							}
						}
					}
					comp->preprocessModels();
					comp->generateDAEs();
					// An override that turned out to be for a state in this mode is undone (continuity of the states).
					for (const shared_ptr<IdPearl> &id : hintedIds) {
						if (id->index >= 0 && id->index < comp->nX) {
							map<string, double>::const_iterator it = state.find((const char *)id->str);
							if (it != state.end()) {
								id->value.d = it->second;
							}
						}
					}
					if (comp->nX != comp->nDE || comp->nY != comp->nAE) {
						fail("the DAE is not well defined: states/differential equations and algebraic variables/equations do not match");
						return NULL;
					}
					vector<string> names;
					vector<pair<string, int>> nameIdx;
					for (const shared_ptr<IdPearl> &id : comp->idTable) {
						if (id->idTy == ID_VAR || id->idTy == ID_IVAR) {
							names.push_back((const char *)id->str);
							nameIdx.push_back(make_pair(string((const char *)id->str), id->index));
						}
					}
					// A mode may reorder the variables or turn a state into an algebraic variable (and
					// back); the results are matched by name. The set of variables must stay the same.
					sort(names.begin(), names.end());
					if (!stateNames.empty() && names != stateNames) {
						ostringstream os;
						os << "the set of continuous variables changes at t=" << t << "; models that add or remove variables in a mode are not supported yet";
						fail(os.str());
						return NULL;
					}
					stateNames = names;

					// State-event detection: samples the guards along each new solution segment.
					HybridEnv scanEnv = env;
					map<int, bool> scanCur = cur;
					bool scanFailed = false;
					// Only the variables the guards refer to are copied into the environment (a grid model has
					// hundreds of variables and a few guard variables; a string-keyed map update per variable
					// and sample dominated the run time), and a solution is evaluated once per sample time.
					const set<string> guardNames = collectGuardNames(guards);
					vector<pair<string, int>> scanIdx;
					for (const pair<string, int> &ni : nameIdx) {
						if (guardNames.count(ni.first) > 0) {
							scanIdx.push_back(ni);
						}
					}
					const SasSolution *scanSol = NULL;
					double scanAlpha = 0.0;
					auto setEnvAt = [&](const shared_ptr<SasSolution> &sol, double alpha) {
						if (scanSol == sol.get() && scanAlpha == alpha) {
							return;
						}
						if (!scanIdx.empty()) {
							arma::vec v = sol->solution->getSolValue(alpha - sol->absStart);
							for (const pair<string, int> &ni : scanIdx) {
								if (ni.second >= 0 && ni.second < (int)v.n_rows) {
									scanEnv.values[ni.first] = v(ni.second);
								}
							}
						}
						scanEnv.time = alpha;
						scanSol = sol.get();
						scanAlpha = alpha;
					};
					auto truthAt = [&](const shared_ptr<SasSolution> &sol, double alpha, map<int, bool> &out) -> bool {
						setEnvAt(sol, alpha);
						if (!evalGuards(guards, scanEnv, out, err)) {
							scanFailed = true;
							return false;
						}
						return true;
					};
					// Margin (left minus right side) of relation k at alpha, k = -1 for all of them.
					const vector<const AstNode *> relations = collectRelations(guards);
					auto marginsAt = [&](const shared_ptr<SasSolution> &sol, double alpha, vector<double> &out) -> bool {
						setEnvAt(sol, alpha);
						out.assign(relations.size(), 0.0);
						for (size_t k = 0; k < relations.size(); k++) {
							if (!relationMargin(relations[k], scanEnv, out[k], err)) {
								scanFailed = true;
								return false;
							}
						}
						return true;
					};
					auto marginOf = [&](const shared_ptr<SasSolution> &sol, double alpha, size_t k, double &m) -> bool {
						setEnvAt(sol, alpha);
						if (!relationMargin(relations[k], scanEnv, m, err)) {
							scanFailed = true;
							return false;
						}
						return true;
					};
					// A relation whose margin has the same sign at both ends of (lo, hi) can still cross zero
					// twice inside (a short pulse between two samples). Looks for the extremum of the margin
					// and returns its time if the margin has the opposite sign there, else -1.
					auto hiddenCrossing = [&](const shared_ptr<SasSolution> &sol, double lo, double hi, const vector<double> &m0, const vector<double> &m1) -> double {
						for (size_t k = 0; k < relations.size(); k++) {
							if (m0[k] == 0.0 || m1[k] == 0.0 || (m0[k] > 0.0) != (m1[k] > 0.0)) {
								continue; // a sign change at the ends is handled by the guard test
							}
							const double sgn = m0[k] > 0.0 ? 1.0 : -1.0;
							double fm;
							double mid = 0.5 * (lo + hi);
							if (!marginOf(sol, mid, k, fm)) {
								return -1.0;
							}
							fm *= sgn;
							if (fm < 0.0) {
								return mid;
							}
							// Parabola through the three samples: f(x) = f0 + b x + c x^2 on [0, 1].
							double f0 = sgn * m0[k];
							double f1 = sgn * m1[k];
							double c = 2.0 * f0 + 2.0 * f1 - 4.0 * fm;
							double b = -3.0 * f0 - f1 + 4.0 * fm;
							if (!(c > 0.0)) {
								continue; // no interior minimum of sgn * margin
							}
							double x = -b / (2.0 * c);
							if (x <= 0.0 || x >= 1.0) {
								continue;
							}
							double fv = f0 + b * x + c * x * x;
							if (fv > 0.25 * min(f0, f1)) {
								continue; // the margin stays well away from zero
							}
							double l = lo;
							double h = hi;
							for (int it = 0; it < 60; it++) {
								double a1 = l + (h - l) / 3.0;
								double a2 = h - (h - l) / 3.0;
								double v1, v2;
								if (!marginOf(sol, a1, k, v1) || !marginOf(sol, a2, k, v2)) {
									return -1.0;
								}
								if (sgn * v1 < sgn * v2) {
									h = a2;
								} else {
									l = a1;
								}
							}
							double xm = 0.5 * (l + h);
							double vm;
							if (!marginOf(sol, xm, k, vm)) {
								return -1.0;
							}
							if (sgn * vm < 0.0 && xm > lo && xm < hi) {
								return xm;
							}
						}
						return -1.0;
					};
					// A guard has flipped if its truth differs from the tracked one, except when-conditions
					// going true -> false, which only re-arm.
					auto flipped = [&](const map<int, bool> &now) {
						for (const Guard &g : guards) {
							if (g.dynamic && now.at(g.id) != scanCur.at(g.id) && !((g.kind == G_WHEN || g.kind == G_ASSERT) && !now.at(g.id))) {
								return true;
							}
						}
						return false;
					};

					SasComputationOptions opts = options;
					opts.startTime = t;
					opts.endTime = segEnd;
					opts.tolerateNonFinite = true; // the part after an event is cut off; the rest is validated below
					if (options.progress) {
						// Report against the end of the whole run, not of this inter-event segment.
						opts.progress = [&](double tm, double) { return options.progress(tm, options.endTime); };
					}
					if (hasDynamic) {
						opts.onSolution = [&](const shared_ptr<SasSolution> &sol) -> bool {
							double a0 = sol->absStart;
							double a1 = sol->absEnd;
							if (!(a1 > a0)) {
								return false;
							}
							// Sample at least guardSamples times per segment. Right after an event the
							// spacing is a fraction of the previous gap between events (bouncing and
							// chattering systems have geometrically shrinking gaps) and then grows with
							// the time since the event, so that quiet stretches stay cheap.
							double span = a1 - a0;
							double base = span / guardSamples;
							double floorStep = base;
							if (prevEventT >= 0.0 && lastEventT > prevEventT) {
								floorStep = min(base, (lastEventT - prevEventT) / 8.0);
							}
							double prevA = a0;
							vector<double> prevM;
							vector<double> curM;
							if (!relations.empty() && !marginsAt(sol, a0, prevM)) {
								return true;
							}
							auto rearm = [&](map<int, bool> &now) {
								for (const Guard &g : guards) {
									if (g.dynamic && (g.kind == G_WHEN || g.kind == G_ASSERT) && !now[g.id]) {
										scanCur[g.id] = false;
									}
								}
							};
							for (int iter = 0; prevA < a1 && iter < 5000000; iter++) {
								// The step must stay representable next to the time (events that follow each
								// other at 1e-16 intervals would otherwise give a step of zero).
								double step = min(base, max(max(floorStep, (prevA - a0) / 32.0), 1e-13 * (1.0 + fabs(prevA))));
								double a = prevA + step;
								if (a1 - a < 0.25 * step) {
									a = a1;
								}
								a = min(a, a1);
								map<int, bool> now;
								if (!truthAt(sol, a, now)) {
									return true;
								}
								if (!relations.empty() && !marginsAt(sol, a, curM)) {
									return true;
								}
								double hi = a;
								bool hit = flipped(now);
								if (!hit && !relations.empty()) {
									double pt = hiddenCrossing(sol, prevA, a, prevM, curM);
									if (scanFailed) {
										return true;
									}
									if (pt > prevA) {
										map<int, bool> pm;
										if (!truthAt(sol, pt, pm)) {
											return true;
										}
										if (flipped(pm)) {
											hit = true;
											hi = pt;
										} else {
											// The pulse does not change a tracked guard (a when-condition that
											// went false is re-armed); continue from the extremum.
											rearm(pm);
											prevA = pt;
											if (!marginsAt(sol, pt, prevM)) {
												return true;
											}
											continue;
										}
									}
								}
								if (!hit) {
									rearm(now);
									prevA = a;
									prevM = curM;
									continue;
								}
								double lo = prevA;
								for (int it = 0; it < 100 && hi - lo > 1e-13 * (1.0 + fabs(hi)); it++) {
									double mid = 0.5 * (lo + hi);
									map<int, bool> m;
									if (!truthAt(sol, mid, m)) {
										return true;
									}
									if (flipped(m)) {
										hi = mid;
									} else {
										lo = mid;
									}
								}
								dynT = hi;
								return true;
							}
							if (prevA < a1) {
								err = "the state-event scan does not advance at t=" + to_string(prevA);
								scanFailed = true;
								return true;
							}
							return false;
						};
					}
					unique_ptr<SasSolutionSet> sols(comp->solve(opts));
					if (hinted && sols->solutionLink.empty() && sols->stalled && !sols->cancelled) {
						// The values of the earlier mode did not lead to a solution: once more from the values
						// of the mode just left.
						allowStartHint = false;
						continue;
					}
					allowStartHint = true;
					if (scanFailed) {
						return NULL;
					}
					if (sols->cancelled) {
						wasCancelled = true;
						fail("cancelled");
						return NULL;
					}
					for (const string &w : sols->warnings) {
						if (find(warningList.begin(), warningList.end(), w) == warningList.end()) {
							warningList.push_back(w);
						}
					}
					if (dynT < 0.0 && sols->stalled) {
						fail(sols->stallMessage);
						return NULL;
					}
					if (dynT >= 0.0) {
						segEnd = dynT;
						while (!sols->solutionLink.empty() && sols->solutionLink.back()->absStart >= dynT) {
							sols->solutionLink.pop_back();
						}
						if (!sols->solutionLink.empty()) {
							sols->solutionLink.back()->absEnd = dynT;
						}
					} else if (sols->solutionLink.empty()) {
						ostringstream os;
						os << "the solver made no progress at t=" << t;
						fail(os.str());
						return NULL;
					} else if (segEnd - sols->solutionLink.back()->absEnd > 100.0 * options.alphaTol) {
						ostringstream os;
						os << "the solver stalled at t=" << sols->solutionLink.back()->absEnd << " before the segment end t=" << segEnd;
						fail(os.str());
						return NULL;
					} else {
						shared_ptr<SasSolution> last = sols->solutionLink.back();
						last->absEnd = max(last->absEnd, segEnd); // close the sub-tolerance gap
					}
					for (const shared_ptr<SasSolution> &kept : sols->solutionLink) {
						string why;
						if (!comp->validateDomain(*kept, &why)) {
							fail("the solution leaves the domain of a function: " + why);
							return NULL;
						}
					}
					if (!sols->solutionLink.empty()) {
						shared_ptr<SasSolution> last = sols->solutionLink.back();
						arma::vec v = last->solution->getSolValue(max(0.0, min(segEnd, last->absEnd) - last->absStart));
						map<string, double> algebraic;
						for (const shared_ptr<IdPearl> &id : comp->idTable) {
							if ((id->idTy == ID_VAR || id->idTy == ID_IVAR) && id->str != NULL && id->index >= comp->nX && id->index < (int)v.n_rows) {
								algebraic[(const char *)id->str] = v(id->index);
							}
						}
						modeSnapshot[segmentModeKey] = algebraic;
						previousModeKey = segmentModeKey;
						for (const shared_ptr<IdPearl> &id : comp->idTable) {
							if ((id->idTy == ID_VAR || id->idTy == ID_IVAR) && id->index >= 0 && id->index < (int)v.n_rows) {
								state[(const char *)id->str] = v(id->index);
							}
						}
					}
					for (const shared_ptr<SasSolution> &s : sols->solutionLink) {
						total->solutionLink.push_back(s);
					}
					total->keepAlive.push_back(comp);
					nSegments++;
					cur = scanCur; // truth at the end of the scanned part of the segment
				}

				bool dynEvent = dynT >= 0.0;
				if (!timeEvent && !dynEvent) {
					t = segEnd;
					break;
				}
				stuck = (segEnd - t < 1e-9 * (1.0 + fabs(t))) ? stuck + 1 : 0;
				if (stuck > 1000) {
					ostringstream os;
					os << "the simulation does not advance at t=" << t << " (Zeno behavior or chattering)" << kZenoAdvice;
					fail(os.str());
					return NULL;
				}
				t = segEnd;
				bool timeNow = timeEvent && (!dynEvent || fabs(dynT - te) <= 1e-9 * max(1.0, fabs(te)));
				size_t nBefore = firedEvents.size();
				if (!handleEvents(t, timeNow, false)) {
					return NULL;
				}
				if (firedEvents.size() > nBefore) {
					snapshot(t);
				}
				if (wasTerminated) {
					break;
				}
				if (firedEvents.size() > maxEvents) {
					fail("too many events");
					return NULL;
				}
			}

			for (const HybridEvent &ev : firedEvents) {
				SasEventRecord rec;
				rec.time = ev.time;
				rec.assignments = ev.assignments;
				rec.reinits = ev.reinits;
				total->events.push_back(rec);
			}
			total->terminated = wasTerminated;
			total->terminateTime = terminateTime;
			total->terminateMessage = terminateMessage;
			total->warnings = warningList;
			finalVals = env.values;
			for (const pair<const string, double> &s : state) {
				finalVals[s.first] = s.second;
			}
			return total.release();
		}

	} // namespace core
} // namespace che
