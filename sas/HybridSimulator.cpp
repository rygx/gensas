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
			const double endTol = 1e-9 * max(1.0, fabs(options.endTime));
			double t = options.startTime;

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

			// Events at the start time (sample() with start == startTime).
			{
				HybridEvent ev;
				if (!fireEvents(clauses, env, t, true, ev, err)) {
					return NULL;
				}
				record(ev);
			}

			const size_t maxEvents = 1000000;
			while (true) {
				double te = 0.0;
				bool hasEvent = nextEventTime(clauses, env, t, te, err);
				if (!err.empty()) {
					return NULL;
				}
				bool atEvent = hasEvent && te < options.endTime - endTol;
				double segEnd = atEvent ? te : options.endTime;

				if (segEnd - t > options.alphaTol) {
					shared_ptr<SasComputationModel> comp = make_shared<SasComputationModel>();
					comp->generateCompModel(model);
					SasModel &cm = *comp->topModels.front();
					if (!specializeModel(cm, env, err)) {
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
					for (const shared_ptr<IdPearl> &id : comp->idTable) {
						if (id->idTy == ID_VAR || id->idTy == ID_IVAR) {
							names.push_back((const char *)id->str);
						}
					}
					if (!stateNames.empty() && names != stateNames) {
						ostringstream os;
						os << "the set or order of the continuous variables changes at t=" << t << "; models whose structure changes are not supported yet";
						fail(os.str());
						return NULL;
					}
					stateNames = names;

					SasComputationOptions opts = options;
					opts.startTime = t;
					opts.endTime = segEnd;
					unique_ptr<SasSolutionSet> sols(comp->solve(opts));
					if (sols->solutionLink.empty()) {
						ostringstream os;
						os << "the solver made no progress at t=" << t;
						fail(os.str());
						return NULL;
					}
					shared_ptr<SasSolution> last = sols->solutionLink.back();
					if (segEnd - last->absEnd > 100.0 * options.alphaTol) {
						ostringstream os;
						os << "the solver stalled at t=" << last->absEnd << " before the segment end t=" << segEnd;
						fail(os.str());
						return NULL;
					}
					last->absEnd = max(last->absEnd, segEnd); // close the sub-tolerance gap
					arma::vec v = last->solution->getSolValue(segEnd - last->absStart);
					for (const shared_ptr<IdPearl> &id : comp->idTable) {
						if ((id->idTy == ID_VAR || id->idTy == ID_IVAR) && id->index >= 0 && id->index < (int)v.n_rows) {
							state[(const char *)id->str] = v(id->index);
						}
					}
					for (const shared_ptr<SasSolution> &s : sols->solutionLink) {
						total->solutionLink.push_back(s);
					}
					total->keepAlive.push_back(comp);
					nSegments++;
				}

				if (!atEvent) {
					t = segEnd;
					break;
				}
				t = segEnd;
				for (const pair<const string, double> &s : state) {
					env.values[s.first] = s.second;
				}
				HybridEvent ev;
				if (!fireEvents(clauses, env, t, false, ev, err)) {
					return NULL;
				}
				record(ev);
				if (firedEvents.size() > maxEvents) {
					fail("too many events");
					return NULL;
				}
			}

			finalVals = env.values;
			for (const pair<const string, double> &s : state) {
				finalVals[s.first] = s.second;
			}
			return total.release();
		}

	} // namespace core
} // namespace che
