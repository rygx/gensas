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
#include "sas/SasResult.h"
#include "sas/SasComputation.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <json/json.h>

using namespace std;
using namespace arma;

namespace che {
	namespace core {

		namespace {
			// Index of the last discrete trace entry at or before time t (a sample takes the values that
			// hold from the start of its solution segment, so the two samples at an event time show the
			// values before and after the event).
			size_t traceIndexAt(const SasDiscreteTrace &trace, double t) {
				size_t k = 0;
				for (size_t m = 0; m < trace.times.size(); m++) {
					if (trace.times[m] <= t + 1e-12 * (1.0 + fabs(t))) {
						k = m;
					}
				}
				return k;
			}
		} // namespace

		SasResult buildResult(SasSolutionSet &set, double interval) {
			SasResult res;
			res.interval = interval;
			if (set.solutionLink.empty() || !(interval > 0.0)) {
				return res;
			}
			SasComputationModel *lastModel = set.solutionLink.back()->linkCompModel;
			for (size_t i = 0; i < lastModel->idTable.size(); i++) {
				const shared_ptr<IdPearl> &id = lastModel->idTable[i];
				if (id->idTy == ID_VAR || id->idTy == ID_IVAR) {
					res.names.push_back((const char *)id->str);
				}
			}
			// Per computation model: row of each name in its solution vector, or -1.
			map<SasComputationModel *, vector<int>> rowOf;
			auto rowsFor = [&](SasComputationModel *cm) -> const vector<int> & {
				map<SasComputationModel *, vector<int>>::iterator it = rowOf.find(cm);
				if (it != rowOf.end()) {
					return it->second;
				}
				map<string, int> pos;
				for (size_t i = 0; i < cm->idTable.size(); i++) {
					const shared_ptr<IdPearl> &id = cm->idTable[i];
					if (id->idTy == ID_VAR || id->idTy == ID_IVAR) {
						pos[(const char *)id->str] = (int)i;
					}
				}
				vector<int> rows;
				for (const string &n : res.names) {
					map<string, int>::const_iterator f = pos.find(n);
					rows.push_back(f == pos.end() ? -1 : f->second);
				}
				return rowOf.insert(make_pair(cm, rows)).first->second;
			};

			res.startTime = set.solutionLink.front()->absStart;
			double maxAlpha = set.solutionLink.back()->absEnd;
			res.endTime = maxAlpha;
			vector<vec> samples;
			vector<SasComputationModel *> sampleModel;
			vector<double> segStart; // start time of the solution segment each sample was taken from
			double t = res.startTime;
			// The grid starts at the start time and always ends with a sample at the end time.
			for (;;) {
				list<shared_ptr<SasSolution>> found = set.findSasSolution(t);
				for (const shared_ptr<SasSolution> &sol : found) {
					samples.push_back(sol->solution->getSolValue(t - sol->absStart));
					sampleModel.push_back(sol->linkCompModel);
					res.time.push_back(t);
					segStart.push_back(sol->absStart);
				}
				if (t >= maxAlpha) {
					break;
				}
				t += interval;
				if (t > maxAlpha - interval / 100) {
					t = maxAlpha;
				}
			}
			res.values.assign(res.names.size(), vector<double>(samples.size(), 0.0));
			for (size_t j = 0; j < samples.size(); j++) {
				const vector<int> &rows = rowsFor(sampleModel[j]);
				for (size_t r = 0; r < res.names.size(); r++) {
					bool ok = rows[r] >= 0 && rows[r] < (int)samples[j].n_rows;
					res.values[r][j] = ok ? samples[j](rows[r]) : datum::nan;
				}
			}

			if (!set.discrete.names.empty() && !set.discrete.times.empty()) {
				res.discreteNames = set.discrete.names;
				res.discrete.assign(res.discreteNames.size(), vector<double>(segStart.size(), 0.0));
				for (size_t j = 0; j < segStart.size(); j++) {
					size_t k = traceIndexAt(set.discrete, segStart[j]);
					for (size_t n = 0; n < res.discreteNames.size(); n++) {
						res.discrete[n][j] = set.discrete.values[k][n];
					}
				}
			}
			res.events = set.events;
			return res;
		}

		namespace {
			class JsonResultWriter : public ResultWriter {
			public:
				std::string format() const override {
					return "json";
				}

				bool write(const SasResult &result, const std::string &path, std::string &error) const override {
					Json::Value res(Json::objectValue);
					Json::Value meta(Json::objectValue);
					meta["format"] = "result";
					meta["version"] = kResultFormatVersion;
					res["gensas"] = meta;

					Json::Value tVals(Json::arrayValue);
					for (double t : result.time) {
						tVals.append(Json::Value(t));
					}
					res["time"] = tVals;
					for (size_t r = 0; r < result.names.size(); r++) {
						Json::Value vVals(Json::arrayValue);
						for (double v : result.values[r]) {
							vVals.append(Json::Value(v));
						}
						res[result.names[r]] = vVals;
					}
					if (!result.discreteNames.empty()) {
						Json::Value dVals(Json::objectValue);
						for (size_t n = 0; n < result.discreteNames.size(); n++) {
							Json::Value series(Json::arrayValue);
							for (double v : result.discrete[n]) {
								series.append(Json::Value(v));
							}
							dVals[result.discreteNames[n]] = series;
						}
						res["discrete"] = dVals;
					}
					if (!result.events.empty()) {
						Json::Value evVals(Json::arrayValue);
						for (const SasEventRecord &ev : result.events) {
							Json::Value e(Json::objectValue);
							e["time"] = Json::Value(ev.time);
							Json::Value as(Json::objectValue);
							for (const pair<string, double> &a : ev.assignments) {
								as[a.first] = Json::Value(a.second);
							}
							e["assignments"] = as;
							Json::Value rs(Json::objectValue);
							for (const pair<string, double> &r : ev.reinits) {
								rs[r.first] = Json::Value(r.second);
							}
							e["reinits"] = rs;
							evVals.append(e);
						}
						res["events"] = evVals;
					}

					std::ofstream out(path.c_str());
					if (!out) {
						error = "cannot open \"" + path + "\" for writing";
						return false;
					}
					Json::StreamWriterBuilder builder;
					builder["commentStyle"] = "None";
					builder["indentation"] = "";
					std::unique_ptr<Json::StreamWriter> writer(builder.newStreamWriter());
					writer->write(res, &out);
					out << std::endl;
					out.close();
					if (!out) {
						error = "error while writing \"" + path + "\"";
						return false;
					}
					return true;
				}
			};

			// CSV needs quoting only for names with a comma, quote or line break.
			string csvField(const string &s) {
				if (s.find_first_of(",\"\r\n") == string::npos) {
					return s;
				}
				string q = "\"";
				for (char c : s) {
					if (c == '"') {
						q += '"';
					}
					q += c;
				}
				return q + "\"";
			}

			class CsvResultWriter : public ResultWriter {
			public:
				std::string format() const override {
					return "csv";
				}

				bool write(const SasResult &result, const std::string &path, std::string &error) const override {
					std::ofstream out(path.c_str());
					if (!out) {
						error = "cannot open \"" + path + "\" for writing";
						return false;
					}
					out << setprecision(numeric_limits<double>::max_digits10);
					out << "time";
					for (const string &n : result.names) {
						out << ',' << csvField(n);
					}
					for (const string &n : result.discreteNames) {
						out << ',' << csvField(n);
					}
					out << '\n';
					for (size_t j = 0; j < result.time.size(); j++) {
						out << result.time[j];
						for (size_t r = 0; r < result.names.size(); r++) {
							out << ',' << result.values[r][j];
						}
						for (size_t r = 0; r < result.discreteNames.size(); r++) {
							out << ',' << result.discrete[r][j];
						}
						out << '\n';
					}
					out.close();
					if (!out) {
						error = "error while writing \"" + path + "\"";
						return false;
					}
					if (result.events.empty()) {
						return true;
					}

					// Events go to a second file: <path without .csv>.events.csv
					string evPath = path;
					if (evPath.size() >= 4 && evPath.compare(evPath.size() - 4, 4, ".csv") == 0) {
						evPath.resize(evPath.size() - 4);
					}
					evPath += ".events.csv";
					std::ofstream ev(evPath.c_str());
					if (!ev) {
						error = "cannot open \"" + evPath + "\" for writing";
						return false;
					}
					ev << setprecision(numeric_limits<double>::max_digits10);
					ev << "time,kind,name,value\n";
					for (const SasEventRecord &e : result.events) {
						for (const pair<string, double> &a : e.assignments) {
							ev << e.time << ",assign," << csvField(a.first) << ',' << a.second << '\n';
						}
						for (const pair<string, double> &r : e.reinits) {
							ev << e.time << ",reinit," << csvField(r.first) << ',' << r.second << '\n';
						}
					}
					ev.close();
					if (!ev) {
						error = "error while writing \"" + evPath + "\"";
						return false;
					}
					return true;
				}
			};
		} // namespace

		std::unique_ptr<ResultWriter> makeResultWriter(const std::string &format) {
			if (format == "json") {
				return std::unique_ptr<ResultWriter>(new JsonResultWriter());
			}
			if (format == "csv") {
				return std::unique_ptr<ResultWriter>(new CsvResultWriter());
			}
			return nullptr;
		}

	} // namespace core
} // namespace che
