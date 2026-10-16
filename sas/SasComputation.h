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
#ifndef SAS_COMPUTATION_H
#define SAS_COMPUTATION_H

#include "sas/SasConfig.h"
#include "sas/SasInput.h"
#include "sas/SasLexico.h"
#include "sas/SasResult.h"
#include <functional>
#include <list>
#include <vector>
#include "util/AbstractCheCalculator.h"

namespace che {
	namespace core {

		class IdReplaceTable {
		public:
			AstNode *ori = NULL;
			AstNode *rep = NULL;
			IdReplaceTable(AstNode *ori, AstNode *rep);
		};

		class SasComputationModel;
		class SasSolution;

		class SasComputationOptions {
		public:
			double startTime = 0.0; // absolute time at which the initial values apply
			double endTime = 10.0;	// absolute end time
			double segment = 1.0;
			double alphaTol = 1e-4;
			double errorTol = 1e-6;
			int nLvl = 15;
			// Called with each solution segment right after it was computed. If it returns true, solve()
			// stops after this segment (used to cut the simulation at a detected event).
			std::function<bool(const std::shared_ptr<SasSolution> &)> onSolution;
			// Called after each solution segment with the absolute time reached and the end time. If it
			// returns true, the run is cancelled (the solution set is marked `cancelled`).
			std::function<bool(double time, double endTime)> progress;
			// Where diagnostic messages (stage lengths, equation counts) go; null (default) is silent.
			std::ostream *log = nullptr;
		};

		class SasSolution {
		public:
			double absStart;
			double absEnd;
			CheSolution *solution = NULL;
			SasComputationModel *linkCompModel = NULL;

			SasSolution(SasComputationModel *linkCompModel, int type, int nState, int nLvl);

			virtual ~SasSolution();
		};

		class SasSolutionSet {
		public:
			SasSolutionSet();

			// Hybrid simulations only: the events that fired and the discrete variables over time.
			std::vector<SasEventRecord> events;
			SasDiscreteTrace discrete;

			std::list<std::shared_ptr<SasSolution>> solutionLink;

			// Set when options.progress asked to stop; the segments computed so far are kept.
			bool cancelled = false;

			// Hybrid simulations: terminate() ended the run at terminateTime, and the warnings of failed
			// assert() calls with AssertionLevel.warning.
			bool terminated = false;
			double terminateTime = 0.0;
			std::string terminateMessage;
			std::vector<std::string> warnings;

			// Set when the solver could not advance (singular solution or a function evaluated outside its
			// domain); the segments computed up to stallTime are kept.
			bool stalled = false;
			double stallTime = 0.0;
			std::string stallMessage;

			// Computation models referenced by the solutions (SasSolution::linkCompModel is a raw pointer).
			std::vector<std::shared_ptr<SasComputationModel>> keepAlive;

			std::list<std::shared_ptr<SasSolution>> findSasSolution(double time);

			std::shared_ptr<SasSolution> findSasSolution(std::shared_ptr<SasComputationModel> &, double time);

			virtual ~SasSolutionSet();
		};

		// An operand with a restricted range (argument of log, sqrt, asin, acos, denominator).
		struct DomainCheck {
			enum Kind { POSITIVE = 0,
						NONNEGATIVE = 1,
						UNIT = 2,
						NONZERO = 3 };
			AstNode *arg;
			Kind kind;
			const char *label;
		};

		class SasComputationModel {
		public:
			// Model storage
			std::list<std::shared_ptr<StrPearl>> strTable;
			std::list<std::shared_ptr<SasModel>> topModels;
			// DAE
			std::vector<shared_ptr<IdPearl>> idTable;
			std::vector<shared_ptr<AstTree>> eqnTable;
			int nDE = 0;
			int nAE = 0;
			int nX = 0;
			int nY = 0;
			// Where diagnostic messages go; null (default) is silent. Set before generateDAEs().
			std::ostream *log = nullptr;

			SasComputationModel();

			void generateCompModel(ModelPool &);

			// Copies a single model (and its equations) from the pool into this computation model.
			void generateCompModel(const std::shared_ptr<SasModel> &);

			void preprocessModels();

			void generateDAEs();

			SasSolutionSet *solve(const SasComputationOptions &options);

			virtual ~SasComputationModel();

		private:
			int tmpVarCnt = 1;

			AstNode *copyAstNodeItself(std::shared_ptr<SasModel> &pModel, AstNode *ori);

			AstNode *copyAstTree(std::shared_ptr<SasModel> &pModel, AstNode *ori);

			std::shared_ptr<StrPearl> generateNewTempVar(std::shared_ptr<SasModel> &pModel);

			bool compareTwoSubTrees(AstNode *aTree, AstNode *bTree);

			void mergeConst(AstNode *aTree);

			void modelFlatten();

			void doReplacement(AstNode *p, bool waiveNext, int depth, std::shared_ptr<SasModel> &currentModel,
							   std::vector<std::list<std::shared_ptr<IdReplaceTable>>> &repTable);

			void processInitialValues();

			vec getStartStateIVP();

			vec calcDiff(CheSolution *solution, double alpha);

			// Why the equations are undefined at a point of a solution (a function argument outside its
			// domain); empty if nothing is found.
			std::string describeDomainViolation(CheSolution *solution, double alpha);

			// True if the operands with restricted ranges stay in their range along [0, alpha] of a
			// solution; otherwise the reason goes to *why.
			bool domainOk(CheSolution *solution, double alpha, std::string *why);
			std::vector<DomainCheck> domainChecks;

			// Messages of the last solve() about steps that were accepted with a relaxed tolerance.
			std::vector<std::string> solverWarnings;
			static const size_t MAX_SOLVER_WARNINGS = 10;
			// Why the last segment could not advance (see SasSolutionSet::stalled).
			std::string stallReason;

			shared_ptr<SasSolution> solveSegment(const vec &init, double curAlpha, double seg, const SasComputationOptions &options);
		};

		enum IdLocatorType { LOCATOR_INIT,
							 LOCATOR_STATE };

		class IdLocator {
		public:
			int type = LOCATOR_INIT;
			std::list<std::shared_ptr<SasModel>> searchDomain;
			arma::vec state;
			arma::vec der;

			IdLocator(const std::list<std::shared_ptr<SasModel>> &);

			void setStates(arma::vec &, arma::vec &);

			void unsetStates();

			static std::shared_ptr<IdPearl> findId(std::list<std::shared_ptr<SasModel>> &searchDomain, const char *iden, int len);

			double findVarValue(const char *iden, int len);

			double findDerValue(const char *iden, int len);
		};
	} // namespace core
} // namespace che

#endif