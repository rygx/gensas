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
#include "util/AbstractCheCalculator.h"
#include "pf/ChePFCalculator.h"
#include "io/MatPsatDataRW.h"
#include "util/SafeArmadillo.h"
#include "util/CheCompUtil.h"
#include "nvwa/pctimer.h"

#include "sas/SasInput.h"
#include "sas/SasLexico.h"
#include "sas/SasExpr.h"
#include "sas/SasComputation.h"
#include "sas/HybridSimulator.h"
#include "sas/SasResult.h"
#include "sas/SasResultMat.h"
#include "sas/ModelChecker.h"
#include "sas/Simulation.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <memory>
// #include <memory>
// #include "debug_new.h"

using namespace arma;
using namespace che::core;
using namespace che::util;
using namespace nvwa;

#ifdef __DIR__
#undef __DIR__
#endif

#ifdef WIN32
#define __DIR__ std::string(__FILE__).substr(0, std::string(__FILE__).find_last_of("\\"))
#elif defined linux
#define __DIR__ (std::getenv("PWD"))
#else
#define __DIR__ ""
#endif

#ifdef WIN32
#include <direct.h>
#define GetCurrentDir _getcwd
#else
#include <unistd.h>
#define GetCurrentDir getcwd
#endif
#include <iomanip>
#include <iostream>

// #define new DEBUG_NEW

// using namespace nvwa;

std::string GetCurrentWorkingDir(void) {
	char buff[FILENAME_MAX];
	GetCurrentDir(buff, FILENAME_MAX);
	std::string current_working_dir(buff);
	return current_working_dir;
}

static void printUsage() {
	cout << "Usage:\n"
			"  app run JOB.json [--set key=value ...] [--validate]       run a Modelica job file\n"
			"  app -g -m file|string -i MODEL [options]                   Modelica simulation from flags\n"
			"  app -p -f CASE [options]                                   AC power flow\n"
			"(app is the executable; app.exe on Windows.) See doc/Usage.md and doc/JobFormat.md."
		 << endl;
}

// Parses a number argument; false (with a message) if it is missing or not a number.
static bool numberArg(int argc, char **argv, int &i, const char *flag, double &out) {
	if (++i >= argc) {
		cerr << "Error: " << flag << " needs a value." << endl;
		return false;
	}
	try {
		size_t used = 0;
		out = stod(argv[i], &used);
		if (used != string(argv[i]).size() || !std::isfinite(out)) {
			throw std::invalid_argument("trailing characters");
		}
	} catch (const std::exception &) {
		cerr << "Error: " << flag << " expects a number, got \"" << argv[i] << "\"." << endl;
		return false;
	}
	return true;
}

static bool stringArg(int argc, char **argv, int &i, const char *flag, string &out) {
	if (++i >= argc) {
		cerr << "Error: " << flag << " needs a value." << endl;
		return false;
	}
	out = argv[i];
	return true;
}

// The flags of `app -g`; they build the same job as a job file.
static bool parseGeneralArgs(int argc, char **argv, ModelicaJob &job, bool &verbose) {
	bool fileMode = false;
	string input;
	for (int i = 2; i < argc; i++) {
		string arg = argv[i];
		double v = 0.0;
		string s;
		if (arg == "--mode" || arg == "-m") {
			if (!stringArg(argc, argv, i, "--mode", s)) {
				return false;
			}
			if (s != "file" && s != "string") {
				cerr << "Error: --mode expects file or string." << endl;
				return false;
			}
			fileMode = s == "file";
		} else if (arg == "--input" || arg == "-i") {
			if (!stringArg(argc, argv, i, "--input", input)) {
				return false;
			}
		} else if (arg == "--output" || arg == "-o") {
			if (!stringArg(argc, argv, i, "--output", job.matOutput)) {
				return false;
			}
		} else if (arg == "--json" || arg == "-j") {
			if (!stringArg(argc, argv, i, "--json", job.jsonOutput)) {
				return false;
			}
		} else if (arg == "--csv" || arg == "-c") {
			if (!stringArg(argc, argv, i, "--csv", job.csvOutput)) {
				return false;
			}
		} else if (arg == "--piecewise") {
			if (!stringArg(argc, argv, i, "--piecewise", s) || (s != "events" && s != "smooth")) {
				cerr << "Error: --piecewise expects events or smooth." << endl;
				return false;
			}
			job.piecewise = s == "smooth" ? PW_SMOOTH : PW_EVENTS;
		} else if (arg == "--smooth-eps") {
			if (!numberArg(argc, argv, i, "--smooth-eps", v)) {
				return false;
			}
			if (!(v > 0.0)) {
				cerr << "Error: --smooth-eps expects a positive number." << endl;
				return false;
			}
			job.smoothEps = v;
		} else if (arg == "--time" || arg == "-t") {
			if (!numberArg(argc, argv, i, "--time", job.options.endTime)) {
				return false;
			}
		} else if (arg == "--level" || arg == "-l") {
			if (!numberArg(argc, argv, i, "--level", v)) {
				return false;
			}
			if (v != std::floor(v) || v < 1 || v > 1000) {
				cerr << "Error: --level expects an integer between 1 and 1000." << endl;
				return false;
			}
			job.options.nLvl = (int)v;
		} else if (arg == "--atol" || arg == "-a") {
			if (!numberArg(argc, argv, i, "--atol", job.options.alphaTol)) {
				return false;
			}
		} else if (arg == "--segment" || arg == "-s") {
			if (!numberArg(argc, argv, i, "--segment", job.options.segment)) {
				return false;
			}
		} else if (arg == "--etol" || arg == "-e") {
			if (!numberArg(argc, argv, i, "--etol", job.options.errorTol)) {
				return false;
			}
		} else if (arg == "--step" || arg == "-k") {
			if (!numberArg(argc, argv, i, "--step", job.outputStep)) {
				return false;
			}
		} else if (arg == "--verbose" || arg == "-v") {
			verbose = true;
		} else {
			cerr << "Error: unknown option " << arg << "." << endl;
			return false;
		}
	}
	if (input.empty()) {
		cerr << "Error: --input is required." << endl;
		return false;
	}
	if (fileMode) {
		job.modelFile = input;
	} else {
		// The string is usually passed in quotes; strip them.
		job.modelSource = input.size() >= 2 && (input.front() == '"' || input.front() == '\'') && input.back() == input.front() ? input.substr(1, input.size() - 2) : input;
		if (!job.modelSource.empty()) {
			cout << job.modelSource << endl;
		}
	}
	return true;
}

// Runs a job, prints events and warnings, and writes the requested outputs. Returns the exit code.
static int runModelica(const ModelicaJob &job, bool validateOnly) {
	Simulation sim(job);
	for (const string &w : sim.warnings()) {
		cout << "Warning: " << w << endl;
	}
	if (validateOnly) {
		const ModelicaJob &j = sim.job();
		cout << "Job is valid: " << (j.modelFile.empty() ? string("<inline model>") : j.modelFile) << ", endTime=" << j.options.endTime << ", segment=" << j.options.segment << ", level=" << j.options.nLvl << ", alphaTol=" << j.options.alphaTol << ", errorTol=" << j.options.errorTol << ", output step=" << j.outputStep << endl;
		return 0;
	}
	if (!sim.run()) {
		cerr << "Error: " << sim.error() << endl;
		return 1;
	}
	for (const SasEventRecord &ev : sim.result().events) {
		cout << "Event at t=" << ev.time << ":";
		for (const auto &a : ev.assignments) {
			cout << " " << a.first << "=" << a.second;
		}
		for (const auto &r : ev.reinits) {
			cout << " reinit(" << r.first << ")=" << r.second;
		}
		cout << endl;
	}
	cout << "Computation time: " << sim.elapsedSeconds() << " s." << endl;
	int status = 0;
	const struct {
		const char *format;
		const string &path;
	} targets[] = {{"json", sim.job().jsonOutput}, {"csv", sim.job().csvOutput}, {"mat", sim.job().matOutput}};
	for (const auto &target : targets) {
		if (target.path.empty()) {
			continue;
		}
		string err;
		if (sim.writeResult(target.format, target.path, err)) {
			cout << target.format << " output written to file " << target.path << "." << endl;
		} else {
			cerr << "Error: " << err << "." << endl;
			status = 1;
		}
	}
	return status;
}

// app run JOB.json [--set key=value ...] [--validate]
static int runJobFile(int argc, char **argv) {
	string path;
	vector<string> sets;
	bool validateOnly = false;
	for (int i = 2; i < argc; i++) {
		string arg = argv[i];
		if (arg == "--set") {
			string s;
			if (!stringArg(argc, argv, i, "--set", s)) {
				return 1;
			}
			sets.push_back(s);
		} else if (arg == "--validate") {
			validateOnly = true;
		} else if (path.empty() && !arg.empty() && arg[0] != '-') {
			path = arg;
		} else {
			cerr << "Error: unexpected argument " << arg << "." << endl;
			return 1;
		}
	}
	if (path.empty()) {
		cerr << "Error: app run needs a job file." << endl;
		return 1;
	}
	ifstream in(path, ios::binary);
	if (!in) {
		cerr << "Error: cannot read the job file " << path << "." << endl;
		return 1;
	}
	stringstream buffer;
	buffer << in.rdbuf();
	string text = buffer.str();
	string err;
	if (!applyJobOverrides(text, sets, err)) {
		cerr << "Error: " << err << endl;
		return 1;
	}
	size_t slash = path.find_last_of("/\\");
	string baseDir = slash == string::npos ? string() : path.substr(0, slash);
	ModelicaJob job;
	if (!parseModelicaJob(text, baseDir, job, err)) {
		cerr << "Error: " << path << ": " << err << endl;
		return 1;
	}
	return runModelica(job, validateOnly);
}

static int realMain(int argc, char **argv) {

	if (argc < 2) {
		printUsage();
		return 1;
	}
	string compMode = argv[1];
	if (compMode == "--help" || compMode == "-h") {
		printUsage();
		return 0;
	}
	cout << "compMode=" << compMode << endl;

	if (compMode == "run") {
		return runJobFile(argc, argv);
	} else if (compMode == "-g") {
		ModelicaJob job;
		bool verbose = false;
		if (!parseGeneralArgs(argc, argv, job, verbose)) {
			return 1;
		}
		job.verbose = verbose;
		return runModelica(job, false);

	} else if (compMode == "-p") {

		string filePath = "";
		string outputPath = "detault.mat";
		int nlvl = 15;
		double segment = 1.0;
		double alphaTol = 1e-4;
		double diffTol = 1e-6;
		double diffTolMax = 1e-2;
		int repeat = 1;
		for (int iArg = 2; iArg < argc; iArg++) {
			string arg = argv[iArg];
			if (arg == "--file" || arg == "-f") {
				if (++iArg < argc) {
					filePath = argv[iArg];
				} else {
					cerr << "File name should be specified after --file or -f." << endl;
				}
			} else if (arg == "--level" || arg == "-l") {
				if (++iArg < argc) {
					nlvl = stoi(argv[iArg]);
				} else {
					cerr << "nlvl should be specified after --level or -l. Using nlvl=" << nlvl << " as default." << endl;
				}
			} else if (arg == "--repeat" || arg == "-r") {
				if (++iArg < argc) {
					repeat = stoi(argv[iArg]);
				} else {
					cerr << "repeat should be specified after --repeat or -r. Using repeat=" << repeat << " as default." << endl;
				}
			} else if (arg == "--segment" || arg == "-s") {
				if (++iArg < argc) {
					segment = stod(argv[iArg]);
				} else {
					cerr << "segment should be specified after --segment or -s. Using segment=" << segment << " as default." << endl;
				}
			} else if (arg == "--alphatol" || arg == "-a") {
				if (++iArg < argc) {
					alphaTol = stod(argv[iArg]);
				} else {
					cerr << "alphaTol should be specified after --alphatol or -a. Using alphaTol=" << alphaTol << " as default." << endl;
				}
			} else if (arg == "--difftol" || arg == "-d") {
				if (++iArg < argc) {
					diffTol = stod(argv[iArg]);
				} else {
					cerr << "diffTol should be specified after --difftol or -d. Using diffTol=" << diffTol << " as default." << endl;
				}
			} else if (arg == "--output" || arg == "-o") {
				if (++iArg < argc) {
					string subArg = argv[iArg];
					outputPath = GetCurrentWorkingDir() + "/" + subArg;
				} else {
					cerr << "Output file name should be specified after --file or -f." << endl;
				}
			} else if (arg == "--difftolmax" || arg == "-e") {
				if (++iArg < argc) {
					diffTolMax = stod(argv[iArg]);
				} else {
					cerr << "diffTolMax should be specified after --difftolmax or -e. Using diffTolMax=" << diffTolMax << " as default." << endl;
				}
			}
		}

		if (nlvl < 3) {
			nlvl = 3;
		}
		if (nlvl > 100) {
			nlvl = 100;
		}
		if (repeat > 100) {
			repeat = 100;
		}
		if (segment < 0.01) {
			segment = 0.01;
		}
		if (segment > 1.0) {
			segment = 1.0;
		}
		if (alphaTol < 1e-8) {
			alphaTol = 1e-8;
		}
		if (alphaTol > 1e-2) {
			alphaTol = 1e-2;
		}
		if (diffTol < 1e-10) {
			diffTol = 1e-10;
		}
		if (diffTol > 1e-2) {
			diffTol = 1e-2;
		}
		if (diffTolMax < 1e-4) {
			diffTolMax = 1e-4;
		}
		if (diffTolMax < 100.0 * diffTol) {
			diffTolMax = 100.0 * diffTol;
		}

		chedata::MatPsatReader matReader;
		chedata::PsatDataSet psatData;
		//  string filePath = GetCurrentWorkingDir() + "/resources/psat_mat/d_dcase2383wp_mod2_ind_zip3.mat";
		// string filePath = GetCurrentWorkingDir() + "/resources/psat_mat/d_014_ind_zip1.mat";
		// string filePath = GetCurrentWorkingDir() + "/resources/psat_mat/d_ei_458_100.mat";
		// string filePath = GetCurrentWorkingDir() + "/resources/psat_mat/d_dcase2383wp_mod2_zip9x.mat";
		// string filePath = GetCurrentWorkingDir() + "/resources/psat_mat/d_70k.mat";
		if (filePath.size() >= 4 && filePath.compare(filePath.size() - 4, 4, ".mat") == 0) {
			cerr << "Error: PSAT .mat files are not read directly. Convert \"" << filePath << "\" with scripts/psat_to_case.py first (see doc/CaseFormat.md)." << endl;
			return -1;
		}
		int flag = matReader.parse(filePath.c_str(), &psatData);
		if (flag != CHE_IO_SUCCESS) {
			cerr << "Error: cannot read the input file." << endl;
			return -1;
		}
		psatData.renumberBuses();
		uvec islands = CheCompUtil::searchIslands(psatData);
		// islands.print("islands");

		CheCompOptions compOpt(nlvl, 1.0, alphaTol, segment, diffTol, diffTolMax);

		pctimer_t totalTime = 0.;

		for (int i = 0; i < repeat; i++) {
			AbstractCheCalculator *pCalculator = new ChePfCalculator(psatData, compOpt, islands);
			pctimer_t stTime = pctimer();
			int pfFlag = pCalculator->calc();
			pctimer_t endTime = pctimer();
			{
				// A compact fingerprint of the converged state, for regression checks.
				CheState finalState = pCalculator->exportResult();
				cout << "Final state: n=" << finalState.state.n_rows << ", norm=" << setprecision(10) << norm(finalState.state, 2) << ", sum=" << accu(finalState.state) << endl;
			}
			pCalculator->writeMatFile(outputPath.c_str());
			delete pCalculator;
			totalTime += endTime - stTime;

			cout << "Computation No." << i + 1 << endl;
		}

		cout << "Computation time: " << totalTime << " s." << endl;

		return 0;

	} else {

		cerr << "The first arg should be run (job file), -g (general SAS) or -p (power flow)." << endl;
		printUsage();
		return 1;
	}
}

int main(int argc, char **argv) {
	try {
		return realMain(argc, argv);
	} catch (const std::exception &e) {
		cerr << "Error: " << e.what() << endl;
	} catch (...) {
		cerr << "Error: unknown failure." << endl;
	}
	return 1;
}
