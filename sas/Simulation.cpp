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
#include "sas/Simulation.h"

#include "sas/HybridSimulator.h"
#include "sas/ModelChecker.h"
#include "sas/SasExpr.h"
#include "sas/SasLexico.h"
#include "sas/SasResultMat.h"

#include "json/json.h"

#include <chrono>
#include <cmath>
#include <iostream>
#include <set>
#include <sstream>

using namespace std;

namespace che {
	namespace core {

		ModelicaJob::ModelicaJob() {
			options.alphaTol = 1e-3;
			options.errorTol = 1e-5;
			options.segment = 1.0;
		}

		namespace {

			bool isAbsolutePath(const string &p) {
				if (p.empty()) {
					return false;
				}
				if (p[0] == '/' || p[0] == '\\') {
					return true;
				}
				return p.size() >= 2 && p[1] == ':' && isalpha((unsigned char)p[0]);
			}

			string resolvePath(const string &baseDir, const string &p) {
				if (baseDir.empty() || isAbsolutePath(p)) {
					return p;
				}
				char last = baseDir.back();
				return (last == '/' || last == '\\') ? baseDir + p : baseDir + "/" + p;
			}

			// Checks that `obj` is an object whose keys are all in `allowed`.
			bool checkKeys(const Json::Value &obj, const string &path, const set<string> &allowed, string &error) {
				if (!obj.isObject()) {
					error = (path.empty() ? string("the job") : path) + ": expected an object";
					return false;
				}
				for (const string &k : obj.getMemberNames()) {
					if (allowed.count(k) == 0) {
						error = "unknown key " + (path.empty() ? k : path + "." + k);
						return false;
					}
				}
				return true;
			}

			bool getNumber(const Json::Value &obj, const string &path, const char *key, double &out, string &error) {
				if (!obj.isMember(key)) {
					return true;
				}
				const Json::Value &v = obj[key];
				if (!v.isNumeric() || !std::isfinite(v.asDouble())) {
					error = path + "." + key + ": expected a number";
					return false;
				}
				out = v.asDouble();
				return true;
			}

			// An object of name -> number (or true/false).
			bool getValueMap(const Json::Value &obj, const string &path, const char *key, map<string, double> &out, string &error) {
				if (!obj.isMember(key)) {
					return true;
				}
				const Json::Value &m = obj[key];
				if (!m.isObject()) {
					error = path + "." + key + ": expected an object that maps names to numbers";
					return false;
				}
				for (const string &name : m.getMemberNames()) {
					const Json::Value &v = m[name];
					if (v.isBool()) {
						out[name] = v.asBool() ? 1.0 : 0.0;
					} else if (v.isNumeric() && std::isfinite(v.asDouble())) {
						out[name] = v.asDouble();
					} else {
						error = path + "." + key + "." + name + ": expected a number";
						return false;
					}
				}
				return true;
			}

			bool getString(const Json::Value &obj, const string &path, const char *key, string &out, string &error) {
				if (!obj.isMember(key)) {
					return true;
				}
				if (!obj[key].isString()) {
					error = path + "." + key + ": expected a string";
					return false;
				}
				out = obj[key].asString();
				return true;
			}

			// Clamps `value` into [lo, hi]; appends a warning when it had to move.
			void clampValue(double &value, double lo, double hi, const char *what, vector<string> &warnings) {
				if (value > hi) {
					ostringstream os;
					os << what << " is limited to at most " << hi << " (" << value << " given); using " << hi << ".";
					warnings.push_back(os.str());
					value = hi;
				} else if (value < lo) {
					ostringstream os;
					os << what << " must be at least " << lo << " (" << value << " given); using " << lo << ".";
					warnings.push_back(os.str());
					value = lo;
				}
			}

		} // namespace

		static bool parseModelicaJobImpl(const string &json, const string &baseDir, ModelicaJob &job, string &error);

		bool parseModelicaJob(const string &json, const string &baseDir, ModelicaJob &job, string &error) {
			try {
				return parseModelicaJobImpl(json, baseDir, job, error);
			} catch (const std::exception &e) {
				error = string("invalid job: ") + e.what();
				return false;
			}
		}

		static bool parseModelicaJobImpl(const string &json, const string &baseDir, ModelicaJob &job, string &error) {
			Json::Value root;
			Json::CharReaderBuilder builder;
			builder["allowComments"] = false;
			string parseErrors;
			istringstream in(json);
			if (!Json::parseFromStream(builder, in, &root, &parseErrors)) {
				error = "invalid JSON: " + parseErrors;
				return false;
			}
			if (!checkKeys(root, "", {"gensas", "kind", "model", "solver", "output"}, error)) {
				return false;
			}
			const Json::Value &head = root["gensas"];
			if (!head.isObject() || !head["format"].isString() || head["format"].asString() != "job") {
				error = "gensas.format: expected \"job\"";
				return false;
			}
			if (!head["version"].isInt() || head["version"].asInt() != 1) {
				error = "gensas.version: unsupported job version (this build reads version 1)";
				return false;
			}
			if (!root["kind"].isString() || root["kind"].asString() != "modelica") {
				error = "kind: only \"modelica\" jobs are supported";
				return false;
			}

			ModelicaJob out;
			// model
			if (!root["model"].isObject()) {
				error = "model: expected an object with \"file\" or \"source\"";
				return false;
			}
			const Json::Value &model = root["model"];
			if (!checkKeys(model, "model", {"file", "source", "className", "piecewise", "smoothEps", "parameters", "start"}, error)) {
				return false;
			}
			if (!getValueMap(model, "model", "parameters", out.parameters, error) || !getValueMap(model, "model", "start", out.start, error)) {
				return false;
			}
			string file, source, piecewise;
			if (!getString(model, "model", "className", out.className, error) || !getString(model, "model", "file", file, error) || !getString(model, "model", "source", source, error) || !getString(model, "model", "piecewise", piecewise, error) || !getNumber(model, "model", "smoothEps", out.smoothEps, error)) {
				return false;
			}
			if (file.empty() == source.empty()) {
				error = "model: give exactly one of \"file\" and \"source\"";
				return false;
			}
			out.modelFile = file.empty() ? string() : resolvePath(baseDir, file);
			out.modelSource = source;
			if (!piecewise.empty()) {
				if (piecewise != "events" && piecewise != "smooth") {
					error = "model.piecewise: expected \"events\" or \"smooth\"";
					return false;
				}
				out.piecewise = piecewise == "smooth" ? PW_SMOOTH : PW_EVENTS;
			}
			if (!(out.smoothEps > 0.0)) {
				error = "model.smoothEps: expected a positive number";
				return false;
			}

			// solver
			if (root.isMember("solver")) {
				const Json::Value &solver = root["solver"];
				if (!checkKeys(solver, "solver", {"endTime", "segment", "level", "alphaTol", "errorTol"}, error)) {
					return false;
				}
				double level = out.options.nLvl;
				if (!getNumber(solver, "solver", "endTime", out.options.endTime, error) || !getNumber(solver, "solver", "segment", out.options.segment, error) || !getNumber(solver, "solver", "level", level, error) || !getNumber(solver, "solver", "alphaTol", out.options.alphaTol, error) || !getNumber(solver, "solver", "errorTol", out.options.errorTol, error)) {
					return false;
				}
				if (level != std::floor(level) || level < 1 || level > 1000) {
					error = "solver.level: expected an integer between 1 and 1000";
					return false;
				}
				out.options.nLvl = (int)level;
			}

			// output
			if (root.isMember("output")) {
				const Json::Value &output = root["output"];
				if (!checkKeys(output, "output", {"step", "json", "csv", "mat"}, error)) {
					return false;
				}
				string json, csv, mat;
				if (!getNumber(output, "output", "step", out.outputStep, error) || !getString(output, "output", "json", json, error) || !getString(output, "output", "csv", csv, error) || !getString(output, "output", "mat", mat, error)) {
					return false;
				}
				out.jsonOutput = json.empty() ? json : resolvePath(baseDir, json);
				out.csvOutput = csv.empty() ? csv : resolvePath(baseDir, csv);
				out.matOutput = mat.empty() ? mat : resolvePath(baseDir, mat);
			}
			job = out;
			return true;
		}

		bool applyJobOverrides(string &json, const vector<string> &assignments, string &error) {
			Json::Value root;
			Json::CharReaderBuilder builder;
			string parseErrors;
			istringstream in(json);
			if (!Json::parseFromStream(builder, in, &root, &parseErrors)) {
				error = "invalid JSON: " + parseErrors;
				return false;
			}
			for (const string &a : assignments) {
				size_t eq = a.find('=');
				if (eq == string::npos || eq == 0) {
					error = "expected key=value, got \"" + a + "\"";
					return false;
				}
				string key = a.substr(0, eq);
				string text = a.substr(eq + 1);
				Json::Value value;
				istringstream vin(text);
				string verr;
				if (!Json::parseFromStream(builder, vin, &value, &verr)) {
					value = text;
				}
				Json::Value *node = &root;
				size_t pos = 0;
				while (true) {
					size_t dot = key.find('.', pos);
					string part = key.substr(pos, dot == string::npos ? string::npos : dot - pos);
					if (part.empty()) {
						error = "bad key in \"" + a + "\"";
						return false;
					}
					if (!node->isObject() && !node->isNull()) {
						error = "cannot set " + key + ": " + key.substr(0, pos > 0 ? pos - 1 : 0) + " is not an object";
						return false;
					}
					node = &(*node)[part];
					if (dot == string::npos) {
						break;
					}
					pos = dot + 1;
				}
				*node = value;
			}
			Json::StreamWriterBuilder writer;
			writer["indentation"] = "";
			json = Json::writeString(writer, root);
			return true;
		}

		void applyJobLimits(ModelicaJob &job, vector<string> &warnings) {
			const double MAX_TIME = 1000.0;
			const double MIN_TIME = 0.0;
			const double MAX_SEGMENT = 10.0;
			const double MIN_SEGMENT = 1e-4;
			const double MAX_ALPHA_IN_SEG = 1;
			const double MIN_ALPHA_IN_SEG = 1e-9;
			const double MAX_ERRTOL = 1e-2;
			const double MIN_ERRTOL = 1e-10;
			const double MAX_OUTSTEP = 1;
			const double MIN_OUTSTEP = 1e-3;
			SasComputationOptions &o = job.options;
			clampValue(o.endTime, MIN_TIME, MAX_TIME, "the simulation length (endTime)", warnings);
			clampValue(o.segment, MIN_SEGMENT, MAX_SEGMENT, "the segment length", warnings);
			clampValue(o.alphaTol, MIN_ALPHA_IN_SEG * o.segment, MAX_ALPHA_IN_SEG * o.segment, "the alpha tolerance", warnings);
			clampValue(o.errorTol, MIN_ERRTOL, MAX_ERRTOL, "the error tolerance", warnings);
			clampValue(job.outputStep, MIN_OUTSTEP, MAX_OUTSTEP, "the output step", warnings);
		}

		Simulation::Simulation(const ModelicaJob &job) : jobData(job) {
			applyJobLimits(jobData, warns);
			nJobWarnings = warns.size();
		}

		bool Simulation::run(const ProgressFn &progress) {
			// The solver reports some numerical failures by throwing; a library call must not.
			try {
				return runImpl(progress);
			} catch (const std::exception &e) {
				err = string("the solver failed: ") + e.what();
			} catch (...) {
				err = "the solver failed with an unknown error";
			}
			res = SasResult();
			return false;
		}

		namespace {

			// Names of the variables that appear as der(x).
			void collectDerivedNames(const AstNode *p, set<string> &names) {
				for (; p != NULL; p = p->next) {
					if (p->op == OP_CALL && p->subs[0] != NULL && p->subs[0]->op == OP_ID && p->subs[0]->value.p != NULL &&
						string((const char *)p->subs[0]->value.p) == "der" && p->subs[1] != NULL && p->subs[1]->op == OP_ID && p->subs[1]->value.p != NULL) {
						names.insert((const char *)p->subs[1]->value.p);
					}
					collectDerivedNames(p->subs[0], names);
					collectDerivedNames(p->subs[1], names);
				}
			}

			string typeName(Type ty) {
				return ty == TY_INT ? "Integer" : ty == TY_BOOL ? "Boolean"
																: "Real";
			}

			string textOf(const unsigned char *s) {
				return s != NULL ? string((const char *)s) : string();
			}

			double valueOf(const IdPearl &id) {
				return (id.ty == TY_INT || id.ty == TY_BOOL) ? (double)id.value.i : id.value.d;
			}

			ModelInfo buildModelInfo(const SasModel &model, bool hybrid) {
				ModelInfo info;
				if (model.modelInfo != nullptr) {
					info.name = textOf(model.modelInfo->str);
					info.description = textOf(model.modelInfo->annotation);
				}
				info.hybrid = hybrid;
				set<string> derived;
				for (const shared_ptr<AstTree> &tree : model.modelTrees) {
					collectDerivedNames(tree->pHead, derived);
				}
				for (const shared_ptr<IdPearl> &id : model.idTable) {
					if (id->str == NULL) {
						continue;
					}
					const string name = (const char *)id->str;
					if (id->idTy == ID_PARAM) {
						ModelParameter p;
						p.name = name;
						p.type = typeName(id->ty);
						p.unit = textOf(id->unit);
						p.description = textOf(id->annotation);
						p.value = valueOf(*id);
						p.constant = id->isConstant;
						p.line = id->line;
						info.parameters.push_back(p);
					} else if (id->idTy == ID_VAR || id->idTy == ID_DISC) {
						ModelVariable v;
						v.name = name;
						v.kind = id->idTy == ID_DISC ? "discrete" : derived.count(name) > 0 ? "state"
																							: "algebraic";
						v.type = typeName(id->ty);
						v.unit = textOf(id->unit);
						v.description = textOf(id->annotation);
						v.hasStart = id->initVal;
						v.start = id->initVal ? valueOf(*id) : 0.0;
						v.input = id->isInput;
						v.line = id->line;
						info.variables.push_back(v);
					}
				}
				return info;
			}

		} // namespace

		string modelInfoJson(const ModelInfo &info) {
			Json::Value root(Json::objectValue);
			root["name"] = info.name;
			root["description"] = info.description;
			Json::Value models(Json::arrayValue);
			for (const string &m : info.models) {
				models.append(m);
			}
			root["models"] = models;
			root["complete"] = info.complete;
			root["hybrid"] = info.hybrid;
			Json::Value vars(Json::arrayValue);
			for (const ModelVariable &v : info.variables) {
				Json::Value j(Json::objectValue);
				j["name"] = v.name;
				j["kind"] = v.kind;
				j["type"] = v.type;
				j["unit"] = v.unit;
				j["description"] = v.description;
				j["start"] = v.hasStart ? Json::Value(v.start) : Json::Value();
				j["input"] = v.input;
				j["line"] = v.line;
				vars.append(j);
			}
			root["variables"] = vars;
			Json::Value params(Json::arrayValue);
			for (const ModelParameter &p : info.parameters) {
				Json::Value j(Json::objectValue);
				j["name"] = p.name;
				j["type"] = p.type;
				j["unit"] = p.unit;
				j["description"] = p.description;
				j["value"] = p.value;
				j["constant"] = p.constant;
				j["line"] = p.line;
				params.append(j);
			}
			root["parameters"] = params;
			Json::StreamWriterBuilder builder;
			builder["indentation"] = "";
			return Json::writeString(builder, root);
		}

		bool Simulation::checkDaeShape(const SasComputationModel &comp) {
			string msg;
			if (comp.nX != comp.nDE) {
				msg = "the number of state variables (" + to_string(comp.nX) + ") differs from the number of differential equations (" + to_string(comp.nDE) + "); the DAE is not well defined";
			} else if (comp.nY != comp.nAE) {
				msg = "the number of algebraic variables (" + to_string(comp.nY) + ") differs from the number of algebraic equations (" + to_string(comp.nAE) + "); the DAE is not well defined";
			}
			if (msg.empty()) {
				return true;
			}
			err = msg;
			Diagnostic d;
			d.message = msg;
			d.source = "checker";
			diags.push_back(d);
			return false;
		}

		bool Simulation::frontEnd(SasModelParser &parser, bool &hybrid, SolverFeatures &solverFeatures) {
			diags.clear();
			info = ModelInfo();
			infoValid = false;
			ostringstream parseMessages; // errors and warnings of the parser go to the result, not to stderr
			parser.diag = &parseMessages;
			if (!jobData.modelFile.empty()) {
				parser.parseFromFile(jobData.modelFile.c_str());
			} else {
				parser.parseString(jobData.modelSource.c_str());
			}
			if (parser.reader.buffer == NULL || parser.reader.size <= 1) {
				err = jobData.modelFile.empty() ? "the model text is empty" : "cannot read the model file " + jobData.modelFile;
				Diagnostic d;
				d.message = err;
				d.source = "parser";
				diags.push_back(d);
				return false;
			}
			parser.currentModel = parser.pool.createNewTopModel();
			parser.currentModel->diag = &parseMessages;
			parser.className = jobData.className;
			parser.parameterOverrides = jobData.parameters;
			parser.startOverrides = jobData.start;
			parseSingleSasModel(&parser);
			if (parser.errCount == 0) {
				// every replacement must have found its declaration
				for (int pass = 0; pass < 2; pass++) {
					const string what = pass == 0 ? "parameters" : "start";
					for (const auto &kv : pass == 0 ? jobData.parameters : jobData.start) {
						if (parser.usedOverrides.count(what + "." + kv.first) == 0) {
							const string msg = what + "." + kv.first + ": the model has no such " + (pass == 0 ? "parameter" : "variable");
							Diagnostic d;
							d.message = msg;
							d.source = "job";
							parser.diagnostics.push_back(d);
							parser.errCount++;
							parseMessages << msg << "\n";
						}
					}
				}
			}
			parser.diag = &cerr; // parseMessages ends with this scope, the model may live longer
			diags = parser.diagnostics;
			if (parser.currentModel != nullptr) {
				parser.currentModel->diag = &cerr;
			}
			if (parser.errCount > 0) {
				info = ModelInfo(); // the model could not be read, but the editor can still list the models of the source
				info.models = parser.modelNames;
				info.complete = false;
				infoValid = true;
				err = "parsing the model failed:\n" + parseMessages.str();
				return false;
			}

			// The default SolverFeatures describe what the plain solver can run today. Models with discrete
			// variables, if/when statements or `time` go to the time-event simulator.
			// abs/min/max/sign in continuous equations become if-expressions (their switching points are state
			// events) or, with PW_SMOOTH, smooth approximations.
			rewritePiecewiseFunctions(*parser.currentModel, jobData.piecewise, jobData.smoothEps);
			hybrid = isHybridModel(*parser.currentModel);
			info = buildModelInfo(*parser.currentModel, hybrid);
			info.models = parser.modelNames;
			infoValid = true;
			if (hybrid) {
				solverFeatures.discreteVariables = solverFeatures.ifExpressions = solverFeatures.ifEquations = true;
				solverFeatures.whenEquations = solverFeatures.piecewiseFunctions = solverFeatures.eventFunctions = true;
				solverFeatures.simulationTime = true;
				solverFeatures.terminateAndAssert = true;
				solverFeatures.timeEventsOnly = true;
				solverFeatures.stateEvents = true;
			}
			ostringstream diagnostics;
			if (checkModel(*parser.currentModel, solverFeatures, jobData.modelFile.empty() ? string("<string>") : jobData.modelFile, diagnostics, &diags) > 0) {
				err = "the model uses constructs that are not supported yet:\n" + diagnostics.str();
				return false;
			}
			return true;
		}

		bool Simulation::check() {
			err.clear();
			SasModelParser parser;
			bool hybrid = false;
			SolverFeatures features;
			if (!frontEnd(parser, hybrid, features)) {
				return false;
			}
			if (hybrid) {
				return true; // the equations of a hybrid model are only known mode by mode
			}
			try {
				SasComputationModel comp;
				comp.generateCompModel(parser.pool);
				comp.preprocessModels();
				comp.generateDAEs();
				return checkDaeShape(comp);
			} catch (const std::exception &e) {
				err = e.what();
				Diagnostic d;
				d.message = err;
				d.source = "checker";
				diags.push_back(d);
				return false;
			}
		}

		bool Simulation::runImpl(const ProgressFn &progress) {
			err.clear();
			warns.resize(nJobWarnings); // warnings of an earlier run
			wasCancelled = false;
			res = SasResult();
			const auto started = chrono::steady_clock::now();

			SasModelParser parser;
			bool hybrid = false;
			SolverFeatures solverFeatures;
			if (!frontEnd(parser, hybrid, solverFeatures)) {
				return false;
			}

			SasComputationOptions options = jobData.options;
			options.progress = progress;
			options.log = jobData.verbose ? &cout : nullptr;
			SasComputationModel compModel; // must outlive the solutions
			compModel.log = options.log;
			unique_ptr<SasSolutionSet> solutions;
			HybridSimulator hybridSim(parser.currentModel);
			if (hybrid) {
				solutions.reset(hybridSim.run(options));
				if (solutions == nullptr) {
					wasCancelled = hybridSim.cancelled();
					err = hybridSim.error();
					warns.insert(warns.end(), hybridSim.warnings().begin(), hybridSim.warnings().end());
					return false;
				}
				warns.insert(warns.end(), hybridSim.warnings().begin(), hybridSim.warnings().end());
			} else {
				compModel.generateCompModel(parser.pool);
				compModel.preprocessModels();
				compModel.generateDAEs();
				// A state without a start value begins at 0, which is easy to overlook.
				// Algebraic variables only need a guess for the solver and are not reported.
				vector<string> noStart;
				for (const shared_ptr<SasModel> &m : compModel.topModels) {
					for (const shared_ptr<IdPearl> &id : m->idTable) {
						if (id->idTy == ID_VAR && id->index >= 0 && id->index < compModel.nX && !id->initVal && id->str != NULL) {
							noStart.push_back((const char *)id->str);
						}
					}
				}
				if (!noStart.empty()) {
					const size_t shown = std::min<size_t>(noStart.size(), 5);
					string names;
					for (size_t i = 0; i < shown; i++) {
						names += (i == 0 ? "" : ", ") + noStart[i];
					}
					if (noStart.size() > shown) {
						names += ", ...";
					}
					warns.push_back("state" + string(noStart.size() == 1 ? "" : "s") + " without a start value begin" + (noStart.size() == 1 ? "s" : "") + " at 0: " + names);
				}
				if (jobData.verbose) {
					cout << "Original model:" << endl;
					parser.currentModel->printAllModelTrees();
					cout << "Copied model:" << endl;
					for (const shared_ptr<SasModel> &m : compModel.topModels) {
						m->printAllModelTrees();
						m->printIdTable();
					}
				}
				if (!checkDaeShape(compModel)) {
					return false;
				}
				solutions.reset(compModel.solve(options));
				if (solutions != nullptr && solutions->cancelled) {
					wasCancelled = true;
					err = "cancelled";
					return false;
				}
				if (solutions != nullptr && solutions->stalled) {
					err = solutions->stallMessage;
					return false;
				}
				if (solutions != nullptr) {
					warns.insert(warns.end(), solutions->warnings.begin(), solutions->warnings.end());
				}
			}

			res = buildResult(*solutions, jobData.outputStep);
			elapsed = chrono::duration<double>(chrono::steady_clock::now() - started).count();
			if (res.empty()) {
				err = "the simulation produced no results";
				return false;
			}
			return true;
		}

		bool Simulation::writeResult(const string &format, const string &path, string &error) const {
			unique_ptr<ResultWriter> writer = format == "mat" ? makeMatResultWriter() : makeResultWriter(format);
			if (writer == nullptr) {
				error = "unknown output format " + format;
				return false;
			}
			return writer->write(res, path, error);
		}

	} // namespace core
} // namespace che
