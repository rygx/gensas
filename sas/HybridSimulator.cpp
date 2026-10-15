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

		HybridSimulator::HybridSimulator(shared_ptr<SasModel> model) : model(model) {}

		bool HybridSimulator::fail(const string &msg) {
			err = msg;
			return false;
		}

		SasSolutionSet *HybridSimulator::run(const SasComputationOptions &options) {
			err.clear();
			firedEvents.clear();
			finalVals.clear();
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
					os << "events accumulate at t=" << T << " (chattering or Zeno behavior)";
					fail(os.str());
					return false;
				}
				prevEventT = lastEventT;
				lastEventT = T;
				for (const pair<const string, double> &s : state) {
					env.values[s.first] = s.second;
				}
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

			if (!handleEvents(t, false, true)) { // sample() events at the start time
				return NULL;
			}
			snapshot(t);

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
					comp->preprocessModels();
					comp->generateDAEs();
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
					if (!stateNames.empty() && names != stateNames) {
						ostringstream os;
						os << "the set or order of the continuous variables changes at t=" << t << "; models whose structure changes are not supported yet";
						fail(os.str());
						return NULL;
					}
					stateNames = names;

					// State-event detection: samples the guards along each new solution segment.
					HybridEnv scanEnv = env;
					map<int, bool> scanCur = cur;
					bool scanFailed = false;
					auto truthAt = [&](const shared_ptr<SasSolution> &sol, double alpha, map<int, bool> &out) -> bool {
						arma::vec v = sol->solution->getSolValue(alpha - sol->absStart);
						for (const pair<string, int> &ni : nameIdx) {
							if (ni.second >= 0 && ni.second < (int)v.n_rows) {
								scanEnv.values[ni.first] = v(ni.second);
							}
						}
						scanEnv.time = alpha;
						if (!evalGuards(guards, scanEnv, out, err)) {
							scanFailed = true;
							return false;
						}
						return true;
					};
					// A guard has flipped if its truth differs from the tracked one, except when-conditions
					// going true -> false, which only re-arm.
					auto flipped = [&](const map<int, bool> &now) {
						for (const Guard &g : guards) {
							if (g.dynamic && now.at(g.id) != scanCur.at(g.id) && !(g.kind == G_WHEN && !now.at(g.id))) {
								return true;
							}
						}
						return false;
					};

					SasComputationOptions opts = options;
					opts.startTime = t;
					opts.endTime = segEnd;
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
							for (int iter = 0; prevA < a1 && iter < 5000000; iter++) {
								double step = min(base, max(floorStep, (prevA - a0) / 32.0));
								double a = prevA + step;
								if (a1 - a < 0.25 * step) {
									a = a1;
								}
								a = min(a, a1);
								map<int, bool> now;
								if (!truthAt(sol, a, now)) {
									return true;
								}
								if (!flipped(now)) {
									for (const Guard &g : guards) { // re-arm
										if (g.dynamic && g.kind == G_WHEN && !now[g.id]) {
											scanCur[g.id] = false;
										}
									}
									prevA = a;
									continue;
								}
								double lo = prevA;
								double hi = a;
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
							return false;
						};
					}
					unique_ptr<SasSolutionSet> sols(comp->solve(opts));
					if (scanFailed) {
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
					if (!sols->solutionLink.empty()) {
						shared_ptr<SasSolution> last = sols->solutionLink.back();
						arma::vec v = last->solution->getSolValue(max(0.0, min(segEnd, last->absEnd) - last->absStart));
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
					os << "the simulation does not advance at t=" << t << " (chattering or Zeno behavior)";
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
			finalVals = env.values;
			for (const pair<const string, double> &s : state) {
				finalVals[s.first] = s.second;
			}
			return total.release();
		}

	} // namespace core
} // namespace che
