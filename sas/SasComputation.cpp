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
#include "sas/SasComputation.h"
#include <map>
#include <memory>
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <limits>
#include <algorithm>

using namespace std;

namespace che {
	namespace core {

		static int tokenPrec[]{
#define OPINFO(op, prec, name, func, opcode) prec,
#include "sas/opinfo.txt"
#undef OPINFO
		};

		static int numOperand[]{
#define OPINFO(op, prec, name, func, opcode) opcode,
#include "sas/opinfo.txt"
#undef OPINFO
		};

		IdReplaceTable::IdReplaceTable(AstNode *ori, AstNode *rep) {
			this->ori = ori;
			this->rep = rep;
		}

		AstNode *SasComputationModel::copyAstNodeItself(shared_ptr<SasModel> &pModel, AstNode *ori) {
			if (ori == NULL) {
				return NULL;
			}
			AstNode *pNewNode = new AstNode();
			pNewNode->isfunc = ori->isfunc;
			pNewNode->nodeKind = ori->nodeKind;
			pNewNode->op = ori->op;
			pNewNode->ty = ori->ty;
			pNewNode->index = ori->index;
			pNewNode->pairIdx = ori->pairIdx;
			if (pNewNode->op == OP_ID) {
				if (pModel != nullptr) {
					char *iden = (char *)ori->value.p;
					shared_ptr<IdPearl> pId = pModel->lookupIdFromModel(iden, strlen(iden));
					if (pId == nullptr) {
						// error: undefined identifier
						pModel->error("Undefined identifier: %s", iden);
						pId = pModel->injectNewId(iden, strlen(iden));
					}
					pNewNode->value.p = pId->str;
				}
			} else if (pNewNode->op == OP_STR) {
				if (pModel != nullptr) {
					char *str = (char *)ori->value.p;
					this->strTable.push_back(make_shared<IdPearl>((uchar *)str, strlen(str)));
				}
			} else if (pNewNode->op == OP_CONST) {
				if (pNewNode->ty == TY_REAL) {
					pNewNode->value.d = ori->value.d;
				} else {
					pNewNode->value.i[0] = ori->value.i[0];
					pNewNode->value.i[1] = ori->value.i[1];
				}
			} else {
				pNewNode->value = ori->value;
			}

			pNewNode->line = ori->line;
			return pNewNode;
		}

		AstNode *SasComputationModel::copyAstTree(shared_ptr<SasModel> &pModel, AstNode *ori) {
			AstNode *pNewNode = copyAstNodeItself(pModel, ori);
			if (ori->subs[0] != NULL) {
				pNewNode->subs[0] = copyAstTree(pModel, ori->subs[0]);
			}
			if (ori->subs[1] != NULL) {
				pNewNode->subs[1] = copyAstTree(pModel, ori->subs[1]);
			}
			// The `next` chain (call arguments, statement lists) is copied with its subtrees.
			if (ori->next != NULL) {
				pNewNode->next = copyAstTree(pModel, ori->next);
			}
			return pNewNode;
		}

		SasComputationModel::SasComputationModel() : strTable(), topModels() {
		}

		void SasComputationModel::generateCompModel(ModelPool &pool) {
			// Only considered top-level models currently
			list<shared_ptr<SasModel>>::iterator oldModelIt;
			for (oldModelIt = pool.topModels.begin(); oldModelIt != pool.topModels.end(); oldModelIt++) {
				generateCompModel(*oldModelIt);
			}
		}

		void SasComputationModel::generateCompModel(const shared_ptr<SasModel> &oldModel) {
			// currently do not consider the crossing-reference of the IDs in models
			// Inside a model
			// First copy the idTables and necklaces
			shared_ptr<SasModel> pNewModel = make_shared<SasModel>();
			pNewModel->modelInfo = make_shared<IdPearl>(*oldModel->modelInfo); // Copy model information
			// Copy the id table
			for (int idxNck = 0; idxNck < pNewModel->idNecklaces.size(); idxNck++) {
				list<shared_ptr<IdPearl>>::iterator idIt;
				for (idIt = oldModel->idNecklaces[idxNck].begin();
					 idIt != oldModel->idNecklaces[idxNck].end();
					 idIt++) {
					shared_ptr<IdPearl> pNewIdPearl = make_shared<IdPearl>(*(*idIt));
					pNewModel->idNecklaces[idxNck].push_back(pNewIdPearl);
					pNewModel->idTable.push_back(pNewIdPearl);
				}
			}
			this->topModels.push_back(pNewModel);

			// Copy the Lexico tree (and check undefined identifiers)
			list<shared_ptr<AstTree>>::iterator oldTreeIt;
			for (oldTreeIt = oldModel->modelTrees.begin(); oldTreeIt != oldModel->modelTrees.end(); oldTreeIt++) {
				shared_ptr<AstTree> pNewTree = make_shared<AstTree>();
				pNewTree->pHead = copyAstTree(pNewModel, (*oldTreeIt)->pHead);
				pNewTree->line = (*oldTreeIt)->line;
				pNewModel->modelTrees.push_back(pNewTree);
			}
		}

		shared_ptr<IdPearl> SasComputationModel::generateNewTempVar(shared_ptr<SasModel> &pModel) {
			char buffer[MAX_ID_LEN];
			memset(buffer, 0, MAX_ID_LEN * sizeof(char));
			buffer[0] = '1';
			int res = tmpVarCnt;
			int div = tmpVarCnt;
			int idx = 1;
			do {
				div = div / 26;
				res = res - div * 26;
				buffer[idx++] = 'A' + res;
				res = div;
			} while (div > 0);
			shared_ptr<IdPearl> pId = pModel->injectNewId(buffer, idx);
			pId->idTy = ID_TMPVAR;
			tmpVarCnt++;
			return pId;
		}

		bool SasComputationModel::compareTwoSubTrees(AstNode *aTree, AstNode *bTree) {
			if (aTree->op != bTree->op)
				return false;
			if (aTree->ty != bTree->ty)
				return false;
			if ((aTree->next == NULL && bTree->next != NULL) || (aTree->next != NULL && bTree->next == NULL))
				return false;
			if ((aTree->subs[0] == NULL && bTree->subs[0] != NULL) || (aTree->subs[0] != NULL && bTree->subs[0] == NULL))
				return false;
			if ((aTree->subs[1] == NULL && bTree->subs[1] != NULL) || (aTree->subs[1] != NULL && bTree->subs[1] == NULL))
				return false;
			if (aTree->op == OP_CONST) {
				if (aTree->ty == TY_INT) {
					if (aTree->value.i[0] != bTree->value.i[0])
						return false;
					if (aTree->value.i[1] != bTree->value.i[1])
						return false;
				} else {
					if (aTree->value.d != bTree->value.d)
						return false;
				}
			}
			if (aTree->op == OP_STR) {
				if (strcmp((char *)aTree->value.p, (char *)bTree->value.p) != 0)
					return false;
			}
			if (aTree->op == OP_ID) {
				if (aTree->value.p != bTree->value.p)
					return false;
			}
			if (aTree->next != NULL) {
				if (!compareTwoSubTrees(aTree->next, bTree->next))
					return false;
			}
			if (aTree->subs[0] != NULL) {
				if (!compareTwoSubTrees(aTree->subs[0], bTree->subs[0]))
					return false;
			}
			if (aTree->subs[1] != NULL) {
				if (!compareTwoSubTrees(aTree->subs[1], bTree->subs[1]))
					return false;
			}
			return true;
		}

		void freeAstTreeRecursive(AstNode *p);

		double calcUnaryOperation(double x, Op op) {
			double d = NAN;
			if (op == OP_POS) {
				d = x;
			} else if (op == OP_MINUS) {
				d = -x;
			}
			return d;
		}

		double calcBinaryOperation(double left, double right, Op op) {
			double d = NAN;
			if (op == OP_ADD) {
				d = left + right;
			} else if (op == OP_SUB || op == OP_EQN) {
				d = left - right;
			} else if (op == OP_MUL) {
				d = left * right;
			} else if (op == OP_DIV) {
				d = left / right;
			} else if (op == OP_POW) {
				d = pow(left, right);
			}
			return d;
		}

		static func defaultFunctions[] =
			{
				{"der", 3, TK_F_DER},
				{"sin", 3, TK_F_SIN},
				{"cos", 3, TK_F_COS},
				{"exp", 3, TK_F_EXP},
				{"sqrt", 4, TK_F_SQRT},
				{"inv", 3, TK_F_INV},
				{"log", 3, TK_F_LOG},
				{"atan", 4, TK_F_ATAN},
				{"asin", 4, TK_F_ASIN},
				{"acos", 4, TK_F_ACOS},
				{NULL, 0, TK_ID}};

		static int findDefaultFunctions(char *str, int len) {
			func *p = defaultFunctions;
			while (p->name) {
				if (p->len == len && strncmp(str, p->name, len) == 0) {
					break;
				}
				p++;
			}
			return p->tok;
		}

		double calcFunctionSingleArg(double x, TokenKind tk) {
			if (tk == TK_F_SIN) {
				return sin(x);
			} else if (tk == TK_F_COS) {
				return cos(x);
			} else if (tk == TK_F_EXP) {
				return exp(x);
			} else if (tk == TK_F_SQRT) {
				return sqrt(x);
			} else if (tk == TK_F_INV) {
				return 1.0 / x;
			} else if (tk == TK_F_LOG) {
				return log(x);
			} else if (tk == TK_F_ATAN) {
				return atan(x);
			} else if (tk == TK_F_ASIN) {
				return asin(x);
			} else if (tk == TK_F_ACOS) {
				return acos(x);
			}
			return NAN;
		}

		double calcTreeValue(AstNode *p, IdLocator *locator) {
			if (p == NULL) {
				return NAN;
			}
			if (p->op == OP_CONST) {
				return p->ty == TY_INT ? (double)p->value.i[0] : p->value.d;
			} else if (p->op == OP_CALL) {
				if (p->subs[0] != NULL && p->subs[1] != NULL) {
					char *funcName = (char *)p->subs[0]->value.p;
					TokenKind tk = (TokenKind)findDefaultFunctions(funcName, strlen(funcName));
					if (tk == TK_F_DER) {
						if (p->subs[1]->op == OP_ID && p->subs[1]->ty != TY_FUNCTION) {
							if (locator->type == LOCATOR_STATE && p->subs[1]->index != -1) {
								return locator->der(p->subs[1]->index);
							} else {
								char *b = (char *)p->subs[1]->value.p;
								return locator->findDerValue(b, strlen(b));
							}
						} else if (p->subs[1]->op == OP_CONST) {
							return 0.0;
						}
					} else {
						return calcFunctionSingleArg(calcTreeValue(p->subs[1], locator), tk);
					}
				} else {
					return NAN;
				}
			} else if (p->op == OP_ID) {
				if (locator != NULL) {
					if (locator->type == LOCATOR_STATE && p->index != -1) {
						return locator->state(p->index);
					} else {
						char *b = (char *)p->value.p;
						return locator->findVarValue(b, strlen(b));
					}
				} else {
					return NAN;
				}
			} else {
				if (numOperand[p->op] == 1) {
					return calcUnaryOperation(calcTreeValue(p->subs[0], locator), p->op);
				} else if (numOperand[p->op] == 2) {
					return calcBinaryOperation(calcTreeValue(p->subs[0], locator), calcTreeValue(p->subs[1], locator), p->op);
				} else {
					return NAN;
				}
			}
			return NAN;
		}

		void tempConvertIntToDouble(AstNode *node) {
			if (node == NULL) {
				return;
			}
			// Temporary: convert all integer variables to double
			if (node->op == OP_CONST && node->ty == TY_INT) {
				int v = node->value.i[0];
				node->value.d = (double)v;
				node->ty = TY_REAL;
			}
			tempConvertIntToDouble(node->subs[0]);
			tempConvertIntToDouble(node->subs[1]);
			tempConvertIntToDouble(node->next);
		}

		void SasComputationModel::mergeConst(AstNode *node) {
			if (node == NULL) {
				return;
			}
			if (node->op == OP_CONST) {
				return;
			}

			if (node->op == OP_ADD || node->op == OP_MUL) {
				// Exchange constant to the left
				if (node->subs[0] != NULL && node->subs[0]->op != OP_CONST && node->subs[1] != NULL && node->subs[1]->op == OP_CONST) {
					AstNode *tempNode = node->subs[0];
					node->subs[0] = node->subs[1];
					node->subs[1] = tempNode;
				}
			}
			if (node->subs[0] != NULL && node->subs[0]->op == OP_POS) {
				AstNode *tempNode = node->subs[0];
				node->subs[0] = tempNode->subs[0];
				tempNode->subs[0] = NULL;
				freeAstTreeRecursive(tempNode);
			}

			mergeConst(node->subs[0]);
			mergeConst(node->subs[1]);

			if (node->op == OP_ADD || node->op == OP_SUB || node->op == OP_MUL || node->op == OP_DIV || node->op == OP_POW) {
				/*
					  (op)
					  /  \      => const
				 const   const
				*/
				// directly calculate value
				if (node->subs[0] != NULL && node->subs[0]->op == OP_CONST && node->subs[1] != NULL && node->subs[1]->op == OP_CONST) {

					node->value.d = calcBinaryOperation(node->subs[0]->value.d, node->subs[1]->value.d, node->op);
					node->op = OP_CONST;
					node->ty = TY_REAL;
					freeAstTreeRecursive(node->subs[0]);
					node->subs[0] = NULL;
					freeAstTreeRecursive(node->subs[1]);
					node->subs[1] = NULL;
					if (node->next != NULL) {
						freeAstTreeRecursive(node->next);
						node->next = NULL;
					}
				}
				/*
					  (op1)
					  /   \      =>     (op2)
				 const     (op2)        /   \
						   /   \     const   ANY
					   const   ANY

					  (op1)
					  /   \      =>     (op1)
				 const     (op2)        /   \
						   /   \     const   ANY
						 ANY   const
				*/
				// merge constants
				else if (node->subs[0] != NULL && node->subs[0]->op == OP_CONST) {
					if (node->subs[1] != NULL && node->subs[1]->subs[1] != NULL && tokenPrec[node->op] == tokenPrec[node->subs[1]->op]) {
						if (node->subs[1]->subs[0] != NULL && node->subs[1]->subs[0]->op == OP_CONST) {
							node->subs[0]->value.d = calcBinaryOperation(node->subs[0]->value.d, node->subs[1]->subs[0]->value.d, node->op);
							node->op = node->subs[1]->op;
							AstNode *tempNode = node->subs[1];
							node->subs[1] = tempNode->subs[1];
							tempNode->subs[1] = NULL;
							freeAstTreeRecursive(tempNode);
						} else if (node->subs[1]->subs[1] != NULL && node->subs[1]->subs[1]->op == OP_CONST) {
							node->subs[0]->value.d = calcBinaryOperation(node->subs[0]->value.d, node->subs[1]->subs[0]->value.d, node->subs[1]->op);
							AstNode *tempNode = node->subs[1];
							node->subs[1] = tempNode->subs[0];
							tempNode->subs[0] = NULL;
							freeAstTreeRecursive(tempNode);
						}
					}
				}
			}
			if (node->op == OP_POS || node->op == OP_MINUS) {
				if (node->subs[0] != NULL && node->subs[0]->op == OP_CONST) {
					double d = node->subs[0]->value.d;
					node->value.d = (node->op == OP_MINUS) ? -d : d;
					node->op = OP_CONST;
					node->ty = TY_REAL;
					node->nodeKind = node->subs[0]->nodeKind;
					node->isfunc = 0;
					freeAstTreeRecursive(node->subs[0]);
					node->subs[0] = NULL;
				}
			}
			if (node->subs[0] != NULL && node->subs[0]->op == OP_MINUS) {
				if (node->subs[0]->subs[0] != NULL && node->subs[0]->subs[0]->op == OP_MINUS) {
					AstNode *tempNode = node->subs[0];
					node->subs[0] = tempNode->subs[0]->subs[0];
					tempNode->subs[0]->subs[0] = NULL;
					freeAstTreeRecursive(tempNode);
				}
			}
		}

		unsigned int calcNodeHash(AstNode *p) {
			unsigned int sHash = 0;
			if (p == NULL) {
				return 1;
			} else if (p->op == OP_ID || p->op == OP_STR) {
				char *b = (char *)p->value.p;
				sHash = SasCompUtil::ELFHashOfTwo(
					SasCompUtil::ELFHashInt(p->op), SasCompUtil::ELFHashStr(b, strlen(b)));
			} else if (p->op == OP_CONST) {
				if (p->ty == TY_INT) {
					sHash = SasCompUtil::ELFHashInt(p->value.i[0]);
				} else {
					sHash = SasCompUtil::ELFHashDouble(p->value.d);
				}
			} else {
				sHash = SasCompUtil::ELFHashInt(p->op);
			}

			sHash = SasCompUtil::ELFHashOfTwo(sHash, calcNodeHash(p->subs[0]));
			sHash = SasCompUtil::ELFHashOfTwo(sHash, calcNodeHash(p->subs[1]));
			sHash = SasCompUtil::ELFHashOfTwo(sHash, calcNodeHash(p->next));
			return sHash;
		}

		void SasComputationModel::doReplacement(AstNode *p, bool waive, int depth, shared_ptr<SasModel> &currentModel,
												vector<list<shared_ptr<IdReplaceTable>>> &repTable) {
			if (p == NULL) {
				return;
			}
			int sDepth = depth;
			if (p->op == OP_MUL || p->op == OP_DIV || p->op == OP_POW || p->op == OP_CALL) {
				sDepth++;
			}
			bool waiveNext = p->op == OP_EQN && depth == 0 &&
							 (p->subs[0]->op == OP_MUL || p->subs[0]->op == OP_DIV || p->subs[0]->op == OP_POW || p->subs[0]->op == OP_CALL);
			doReplacement(p->subs[0], waiveNext, sDepth, currentModel, repTable);
			doReplacement(p->subs[1], false, sDepth, currentModel, repTable);
			doReplacement(p->next, false, sDepth, currentModel, repTable);

			if (!waive) {
				if ((p->op == OP_MUL && depth > 0) || p->op == OP_DIV || p->op == OP_POW || p->op == OP_CALL ||
					(p->op == OP_ADD && depth > 0) || (p->op == OP_SUB && depth > 0) ||
					(p->op == OP_POS && depth > 0) || (p->op == OP_MINUS && depth > 0)) {
					unsigned int sHash = calcNodeHash(p) & MAX_ID_PEARL;
					list<shared_ptr<IdReplaceTable>>::iterator it;
					AstNode *pMatch = NULL;
					for (it = repTable[sHash].begin(); it != repTable[sHash].end(); it++) {
						if (compareTwoSubTrees(p, (*it)->ori) && p != (*it)->ori) {
							pMatch = (*it)->rep;
						}
					}
					if (pMatch != NULL) {
						p->op = pMatch->op;
						p->ty = pMatch->ty;
						p->value = pMatch->value;
						freeAstTreeRecursive(p->subs[0]);
						freeAstTreeRecursive(p->subs[1]);
						freeAstTreeRecursive(p->next);
						p->subs[0] = NULL;
						p->subs[1] = NULL;
						p->next = NULL;
					} else {
						shared_ptr<IdPearl> pId = generateNewTempVar(currentModel);

						AstNode *pNewTreeHead = new AstNode();
						pNewTreeHead->op = OP_EQN;
						pNewTreeHead->subs[0] = new AstNode();
						pNewTreeHead->subs[0]->op = p->op;
						pNewTreeHead->subs[0]->subs[0] = p->subs[0];
						pNewTreeHead->subs[0]->subs[1] = p->subs[1];
						pNewTreeHead->subs[0]->next = p->next;
						pNewTreeHead->subs[1] = new AstNode();
						pNewTreeHead->subs[1]->op = OP_ID;
						pNewTreeHead->subs[1]->value.p = pId->str;

						p->op = OP_ID;
						p->subs[0] = NULL;
						p->subs[1] = NULL;
						p->next = NULL;
						p->value.p = pId->str;

						shared_ptr<AstTree> tree = make_shared<AstTree>();
						tree->pHead = pNewTreeHead;
						// Note the order of inserting! For determining initial values.
						currentModel->modelTrees.push_back(tree);

						sHash = calcNodeHash(pNewTreeHead->subs[0]) & MAX_ID_PEARL;
						repTable[sHash].push_back(make_shared<IdReplaceTable>(pNewTreeHead->subs[0], pNewTreeHead->subs[1]));

						// some functions need pairing functions
						if (pNewTreeHead->subs[0]->op == OP_CALL) {
							if (pNewTreeHead->subs[0]->subs[0]->op == OP_ID && pNewTreeHead->subs[0]->subs[0]->ty == TY_FUNCTION) {
								char *funcName = (char *)pNewTreeHead->subs[0]->subs[0]->value.p;
								int tok = findDefaultFunctions(funcName, strlen(funcName));
								if (tok == TK_F_SIN) {
									shared_ptr<IdPearl> pairId = generateNewTempVar(currentModel);

									AstNode *pairTreeHead = new AstNode();
									pairTreeHead->op = OP_EQN;
									pairTreeHead->subs[0] = new AstNode();
									pairTreeHead->subs[0]->op = OP_CALL;
									pairTreeHead->subs[0]->subs[0] = new AstNode();
									pairTreeHead->subs[0]->subs[0]->op = OP_ID;
									pairTreeHead->subs[0]->subs[0]->ty = TY_FUNCTION;
									shared_ptr<IdPearl> funcId = currentModel->lookupIdentifier("cos", 3);
									pairTreeHead->subs[0]->subs[0]->value.p = funcId->str;
									pairTreeHead->subs[0]->subs[1] = copyAstTree(currentModel, pNewTreeHead->subs[0]->subs[1]); // copy the tree HERE
									pairTreeHead->subs[1] = new AstNode();
									pairTreeHead->subs[1]->op = OP_ID;
									pairTreeHead->subs[1]->value.p = pairId->str;

									shared_ptr<AstTree> pairTree = make_shared<AstTree>();
									pairTree->pHead = pairTreeHead;
									// Note the order of inserting! For determining initial values.
									currentModel->modelTrees.push_back(pairTree);

									sHash = calcNodeHash(pairTreeHead->subs[0]) & MAX_ID_PEARL;
									repTable[sHash].push_back(make_shared<IdReplaceTable>(pairTreeHead->subs[0], pairTreeHead->subs[1]));

									pairId->pair = pId.get();
									pId->pair = pairId.get();
								} else if (tok == TK_F_COS) {
									shared_ptr<IdPearl> pairId = generateNewTempVar(currentModel);

									AstNode *pairTreeHead = new AstNode();
									pairTreeHead->op = OP_EQN;
									pairTreeHead->subs[0] = new AstNode();
									pairTreeHead->subs[0]->op = OP_CALL;
									pairTreeHead->subs[0]->subs[0] = new AstNode();
									pairTreeHead->subs[0]->subs[0]->op = OP_ID;
									pairTreeHead->subs[0]->subs[0]->ty = TY_FUNCTION;
									shared_ptr<IdPearl> funcId = currentModel->lookupIdentifier("sin", 3);
									pairTreeHead->subs[0]->subs[0]->value.p = funcId->str;
									pairTreeHead->subs[0]->subs[1] = copyAstTree(currentModel, pNewTreeHead->subs[0]->subs[1]); // copy the tree HERE
									pairTreeHead->subs[1] = new AstNode();
									pairTreeHead->subs[1]->op = OP_ID;
									pairTreeHead->subs[1]->value.p = pairId->str;

									shared_ptr<AstTree> pairTree = make_shared<AstTree>();
									pairTree->pHead = pairTreeHead;
									// Note the order of inserting! For determining initial values.
									currentModel->modelTrees.push_back(pairTree);

									sHash = calcNodeHash(pairTreeHead->subs[0]) & MAX_ID_PEARL;
									repTable[sHash].push_back(make_shared<IdReplaceTable>(pairTreeHead->subs[0], pairTreeHead->subs[1]));

									pairId->pair = pId.get();
									pId->pair = pairId.get();
								}
							}
						}
					}
				}
			}
		}

		// Replace parameters with concrete values, this will happen before pre-process the model
		void regulateIds(AstNode *p, shared_ptr<SasModel> &currentModel) {
			if (p == NULL) {
				return;
			}
			if (p->op == OP_ID) {
				char *b = (char *)p->value.p;
				shared_ptr<IdPearl> idp = currentModel->lookupIdFromModel(b, strlen(b));
				if (idp == nullptr) {
					currentModel->error("ID \"%s\" not found.", b);
				} else if (idp->idTy == ID_PARAM) {
					p->op = OP_CONST;
					p->ty = TY_REAL;
					p->value.d = idp->ty == TY_INT ? (double)idp->value.i : idp->value.d;
				}
			}
			regulateIds(p->subs[0], currentModel);
			regulateIds(p->subs[1], currentModel);
			regulateIds(p->next, currentModel);
		}

		void SasComputationModel::modelFlatten() {
			vector<list<shared_ptr<IdReplaceTable>>> replaceTable(ID_HASH_BUCKETS);

			list<shared_ptr<SasModel>>::iterator modelIt;
			for (modelIt = topModels.begin(); modelIt != topModels.end(); modelIt++) {
				list<shared_ptr<AstTree>>::iterator treeIt;
				for (treeIt = (*modelIt)->modelTrees.begin(); treeIt != (*modelIt)->modelTrees.end(); treeIt++) {
					doReplacement((*treeIt)->pHead, false, 0, *modelIt, replaceTable);
				}
			}
		}

		void markEqnTypes(shared_ptr<AstTree> pTree) {
			AstNode *pHead = pTree->pHead;
			if (pHead == NULL || pHead->subs[0] == NULL || pHead->subs[1] == NULL || pHead->op != OP_EQN) {
				pTree->eqnType = EQN_NONE;
				return;
			}
			if (pHead->subs[0]->op == OP_CALL) {
				if (pHead->subs[0]->subs[0] != NULL && pHead->subs[0]->subs[0]->op == OP_ID) {
					char *b = (char *)pHead->subs[0]->subs[0]->value.p;
					TokenKind tk = (TokenKind)findDefaultFunctions(b, strlen(b));
					if (tk == TK_F_DER) {
						pTree->eqnType = EQN_DIFF;
						return;
					} else {
						pTree->eqnType = EQN_ALGE; // incomplete, there may be der terms, which makes eqn EQN_COMP
						return;
					}
				}
			}
			pTree->eqnType = EQN_ALGE;
			return;
		}

		// True if every variable in the expression has an initial value.
		static bool allIdsInitialized(AstNode *p, list<shared_ptr<SasModel>> &searchDomain) {
			if (p == NULL) {
				return true;
			}
			if (p->op == OP_ID && p->ty != TY_FUNCTION) {
				char *b = (char *)p->value.p;
				shared_ptr<IdPearl> id = IdLocator::findId(searchDomain, b, strlen(b));
				return id != nullptr && id->initVal && allIdsInitialized(p->next, searchDomain);
			}
			if (p->op == OP_CALL) {
				return allIdsInitialized(p->subs[1], searchDomain) && allIdsInitialized(p->next, searchDomain);
			}
			return allIdsInitialized(p->subs[0], searchDomain) && allIdsInitialized(p->subs[1], searchDomain) && allIdsInitialized(p->next, searchDomain);
		}

		// Initial values of the variables that no start attribute covers: the value of a bare identifier on one
		// side of an equation is taken from the other side. Definitions of temporaries come first
		// (sqrt(u) = T, T temporary): an alias such as z = T must not give T a value before its definition does.
		void SasComputationModel::processInitialValues() {
			list<shared_ptr<SasModel>>::iterator modelIt;
			for (modelIt = topModels.begin(); modelIt != topModels.end(); modelIt++) {
				list<shared_ptr<SasModel>> searchDomain;
				searchDomain.push_back(*modelIt);
				IdLocator locator(searchDomain);
				for (int pass = 0; pass < 2; pass++) {
					list<shared_ptr<AstTree>>::iterator treeIt;
					for (treeIt = (*modelIt)->modelTrees.begin(); treeIt != (*modelIt)->modelTrees.end(); treeIt++) {
						AstNode *head = (*treeIt)->pHead;
						if (head == NULL || head->subs[0] == NULL || head->subs[1] == NULL || head->op != OP_EQN) {
							continue;
						}
						// "der(x) = T" says nothing about the value of T at the start; only algebraic equations do.
						if ((*treeIt)->eqnType == EQN_DIFF) {
							continue;
						}
						const bool definition = head->subs[0]->op != OP_ID;
						if (definition != (pass == 0)) {
							continue;
						}
						if (head->subs[1]->op == OP_ID) { // must be an ID
							char *b = (char *)head->subs[1]->value.p;
							shared_ptr<IdPearl> id = IdLocator::findId(searchDomain, b, strlen(b));
							if (id != nullptr && !id->initVal) {
								id->value.d = calcTreeValue(head->subs[0], &locator);
								id->initVal = true;
							}
						}
						if (!definition) {
							// a = expression: the left side takes the value of the expression if it has none and
							// everything the expression uses has a value.
							char *a = (char *)head->subs[0]->value.p;
							shared_ptr<IdPearl> left = IdLocator::findId(searchDomain, a, strlen(a));
							if (left != nullptr && !left->initVal && allIdsInitialized(head->subs[1], searchDomain)) {
								double v = calcTreeValue(head->subs[1], &locator);
								if (std::isfinite(v)) {
									left->value.d = v;
									left->initVal = true;
								}
							}
						}
					}
				}
			}
		}

		void SasComputationModel::preprocessModels() {
			list<shared_ptr<SasModel>>::iterator modelIt;
			for (modelIt = topModels.begin(); modelIt != topModels.end(); modelIt++) {
				list<shared_ptr<AstTree>>::iterator treeIt;
				for (treeIt = (*modelIt)->modelTrees.begin(); treeIt != (*modelIt)->modelTrees.end(); treeIt++) {
					regulateIds((*treeIt)->pHead, (*modelIt));
					tempConvertIntToDouble((*treeIt)->pHead);
					mergeConst((*treeIt)->pHead);
				}
			}
			modelFlatten();

			for (modelIt = topModels.begin(); modelIt != topModels.end(); modelIt++) {
				list<shared_ptr<AstTree>>::iterator treeIt;
				for (treeIt = (*modelIt)->modelTrees.begin(); treeIt != (*modelIt)->modelTrees.end(); treeIt++) {
					markEqnTypes(*treeIt);
				}
			}

			processInitialValues();
		}

		void numberNodesRecursive(AstNode *p, IdLocator *locator) {
			if (p == NULL) {
				return;
			}
			if (p->op == OP_ID && p->ty != TY_FUNCTION) {
				char *b = (char *)p->value.p;
				shared_ptr<IdPearl> pId = IdLocator::findId(locator->searchDomain, b, strlen(b));
				if (pId != nullptr) {
					p->index = pId->index;
					if (pId->pair != NULL) {
						p->pairIdx = pId->pair->index;
					}
				} else {
					cerr << "Cannot find ID: " << b << endl;
				}
			}
			numberNodesRecursive(p->subs[0], locator);
			numberNodesRecursive(p->subs[1], locator);
			numberNodesRecursive(p->next, locator);
		}

		void SasComputationModel::generateDAEs() {
			idTable.clear();
			eqnTable.clear();
			vector<shared_ptr<IdPearl>> yTable;
			vector<shared_ptr<AstTree>> aeTable;
			nDE = 0;
			nAE = 0;
			nX = 0;
			nY = 0;

			list<shared_ptr<SasModel>>::iterator modelIt;
			for (modelIt = topModels.begin(); modelIt != topModels.end(); modelIt++) {
				list<shared_ptr<AstTree>>::iterator treeIt;
				list<shared_ptr<SasModel>> searchDomain;
				searchDomain.push_back(*modelIt);
				IdLocator locator(searchDomain);
				for (treeIt = (*modelIt)->modelTrees.begin(); treeIt != (*modelIt)->modelTrees.end(); treeIt++) {
					if ((*treeIt)->eqnType == EQN_DIFF) {
						(*treeIt)->eqnIdx = nDE;
						eqnTable.push_back(*treeIt);
						nDE++;
						if ((*treeIt)->pHead != NULL && (*treeIt)->pHead->subs[0] != NULL && (*treeIt)->pHead->subs[1] != NULL && (*treeIt)->pHead->op == OP_EQN && (*treeIt)->pHead->subs[0]->op == OP_CALL && (*treeIt)->pHead->subs[0]->subs[0] != NULL && (*treeIt)->pHead->subs[0]->subs[1] != NULL && (*treeIt)->pHead->subs[0]->subs[1]->op == OP_ID) {
							// retrieve the x from der(x)=*;
							char *b = (char *)(*treeIt)->pHead->subs[0]->subs[1]->value.p;
							shared_ptr<IdPearl> id = IdLocator::findId(searchDomain, b, strlen(b));
							if (id != nullptr && id->index == -1) {
								id->index = nX;
								idTable.push_back(id);
								nX++;
							} else if (id == nullptr) {
								(*modelIt)->error("Unknown ID: %s.", b);
							} else {
								(*modelIt)->error("conflicting differential equations for der(%s).", b);
							}
						}
					} else if ((*treeIt)->eqnType == EQN_ALGE) {
						aeTable.push_back(*treeIt);
					}
				}

				list<shared_ptr<IdPearl>>::iterator idIt;
				for (idIt = (*modelIt)->idTable.begin(); idIt != (*modelIt)->idTable.end(); idIt++) {
					if (((*idIt)->idTy == ID_VAR || (*idIt)->idTy == ID_TMPVAR) && ((*idIt)->index == -1)) {
						yTable.push_back(*idIt);
					}
				}
			}

			vector<shared_ptr<IdPearl>>::iterator yIt;
			for (yIt = yTable.begin(); yIt != yTable.end(); yIt++) {
				(*yIt)->index = nX + nY;
				nY++;
			}
			vector<shared_ptr<AstTree>>::iterator aeIt;
			for (aeIt = aeTable.begin(); aeIt != aeTable.end(); aeIt++) {
				(*aeIt)->eqnIdx = nDE + nAE;
				nAE++;
			}
			idTable.insert(idTable.end(), yTable.begin(), yTable.end());
			eqnTable.insert(eqnTable.end(), aeTable.begin(), aeTable.end());

			if (log != nullptr) {
				*log << "Number of DE: " << nDE << endl;
			}
			if (log != nullptr) {
				*log << "Number of AE: " << nAE << endl;
			}
			if (log != nullptr) {
				*log << "Number of X: " << nX << endl;
			}
			if (log != nullptr) {
				*log << "Number of Y: " << nY << endl;
			}
			if (nDE != nX) {
				cerr << "Number of DE not equal number of X" << endl;
			}
			if (nAE != nY) {
				cerr << "Number of AE not equal number of Y" << endl;
			}

			for (modelIt = topModels.begin(); modelIt != topModels.end(); modelIt++) {
				list<shared_ptr<AstTree>>::iterator treeIt;
				list<shared_ptr<SasModel>> searchDomain;
				searchDomain.push_back(*modelIt);
				IdLocator locator(searchDomain);
				for (treeIt = (*modelIt)->modelTrees.begin(); treeIt != (*modelIt)->modelTrees.end(); treeIt++) {
					numberNodesRecursive((*treeIt)->pHead, &locator);
				}
			}
		}

		vec SasComputationModel::getStartStateIVP() {
			int nState = idTable.size();
			vec state(nState, fill::zeros);
			vector<shared_ptr<IdPearl>>::iterator idIt;
			for (idIt = idTable.begin(); idIt != idTable.end(); idIt++) {
				state((*idIt)->index) = (*idIt)->value.d;
			}
			return state;
		}

		// The inverse functions y = log(x), atan(x), asin(x) and acos(x) satisfy y' = sign * w * x' with
		// w = g(x)^alpha, where g = x, 1+x^2, 1-x^2 and 1-x^2 and alpha = -1, -1, -1/2 and -1/2.
		static bool inverseFunction(int tok, double &sign, double &alpha) {
			switch (tok) {
			case TK_F_LOG:
			case TK_F_ATAN:
				sign = 1.0;
				alpha = -1.0;
				return true;
			case TK_F_ASIN:
				sign = 1.0;
				alpha = -0.5;
				return true;
			case TK_F_ACOS:
				sign = -1.0;
				alpha = -0.5;
				return true;
			}
			return false;
		}

		// Series w_0..w_{count-1} of w = g(x)^alpha, given the series x_0..x_{count-1} of the argument.
		static vector<double> inverseFunctionFactor(int tok, double alpha, const vector<double> &x, int count) {
			vector<double> g(count), w(count);
			for (int k = 0; k < count; k++) {
				if (tok == TK_F_LOG) {
					g[k] = x[k];
					continue;
				}
				double square = 0.0;
				for (int i = 0; i <= k; i++) {
					square += x[i] * x[k - i];
				}
				g[k] = (tok == TK_F_ATAN) ? square : -square;
				if (k == 0) {
					g[0] += 1.0;
				}
			}
			w[0] = pow(g[0], alpha);
			for (int m = 1; m < count; m++) {
				double s = 0.0;
				for (int k = 1; k <= m; k++) {
					s += (alpha * k - (m - k)) * g[k] * w[m - k];
				}
				w[m] = s / (m * g[0]);
			}
			return w;
		}

		// Prerequisites: (1) the equation has been flattened; (2) the equations and IDs have been numbered.
		void insertLHS(AstNode *p, const vec &init, vector<int> &locRow, vector<int> &locCol, vector<double> &value, int row, int colOffset, double coeff) {
			if (p == NULL) {
				return;
			}
			if (p->op == OP_ID) {
				if (p->index - colOffset >= 0) {
					locRow.push_back(row);
					locCol.push_back(p->index - colOffset);
					value.push_back(coeff);
				}
				return;
			} else if (p->op == OP_EQN) {
				if (p->subs[0] != NULL && p->subs[1] != NULL && p->subs[0]->op == OP_CALL) {
					if (p->subs[1]->op == OP_ID && p->subs[0]->subs[1]->op == OP_ID) {
						char *funcName = (char *)p->subs[0]->subs[0]->value.p;
						int tok = findDefaultFunctions(funcName, strlen(funcName));
						if (tok == TK_F_SIN) {																								  // sin(x)=y;
							insertLHS(p->subs[0]->subs[1], init, locRow, locCol, value, row, colOffset, coeff * init(p->subs[1]->pairIdx));	  // x coefficient
							insertLHS(p->subs[1], init, locRow, locCol, value, row, colOffset, -coeff);										  // y coefficient
						} else if (tok == TK_F_COS) {																						  // cos(x)=y;
							insertLHS(p->subs[0]->subs[1], init, locRow, locCol, value, row, colOffset, -coeff * init(p->subs[1]->pairIdx));  // x coefficient
							insertLHS(p->subs[1], init, locRow, locCol, value, row, colOffset, -coeff);										  // y coefficient
						} else if (tok == TK_F_EXP) {																						  // exp(x)=y;
							insertLHS(p->subs[0]->subs[1], init, locRow, locCol, value, row, colOffset, coeff * init(p->subs[1]->index));	  // x coefficient
							insertLHS(p->subs[1], init, locRow, locCol, value, row, colOffset, -coeff);										  // y coefficient
						} else if (tok == TK_F_SQRT) {																						  // sqrt(x)=y;
							insertLHS(p->subs[0]->subs[1], init, locRow, locCol, value, row, colOffset, coeff / 2 / init(p->subs[1]->index)); // x coefficient
							insertLHS(p->subs[1], init, locRow, locCol, value, row, colOffset, -coeff);										  // y coefficient
						} else if (tok == TK_F_LOG || tok == TK_F_ATAN || tok == TK_F_ASIN || tok == TK_F_ACOS) {							  // y=f(x), y'=sign*w*x'
							double sign = 1.0, alpha = 1.0;
							inverseFunction(tok, sign, alpha);
							const double x0 = init(p->subs[0]->subs[1]->index);
							const double g0 = (tok == TK_F_LOG) ? x0 : (tok == TK_F_ATAN ? 1.0 + x0 * x0 : 1.0 - x0 * x0);
							insertLHS(p->subs[0]->subs[1], init, locRow, locCol, value, row, colOffset, coeff * sign * pow(g0, alpha)); // x coefficient
							insertLHS(p->subs[1], init, locRow, locCol, value, row, colOffset, -coeff);									// y coefficient
						}
						// There should not be der(x) because only AEs are processed.
					} // TODO: There may be cases such as sin(x)=0;
				} else if (p->subs[0] != NULL && p->subs[1] != NULL && p->subs[0]->op == OP_DIV) { // x/y=z;
					if (p->subs[0]->subs[0] != NULL && p->subs[0]->subs[1] != NULL) {
						double y0 = p->subs[0]->subs[1]->op == OP_ID ? init(p->subs[0]->subs[1]->index) : p->subs[0]->subs[1]->value.d;
						double z0 = p->subs[1]->op == OP_ID ? init(p->subs[1]->index) : p->subs[1]->value.d;

						insertLHS(p->subs[0]->subs[0], init, locRow, locCol, value, row, colOffset, coeff / y0);	   // x coefficient
						insertLHS(p->subs[0]->subs[1], init, locRow, locCol, value, row, colOffset, -coeff * z0 / y0); // y coefficient
						insertLHS(p->subs[1], init, locRow, locCol, value, row, colOffset, -coeff);					   // y coefficient
					}
				} else {
					insertLHS(p->subs[0], init, locRow, locCol, value, row, colOffset, coeff);	// LHS
					insertLHS(p->subs[1], init, locRow, locCol, value, row, colOffset, -coeff); // RHS
				}
			} else if (p->op == OP_ADD) {
				insertLHS(p->subs[0], init, locRow, locCol, value, row, colOffset, coeff);
				insertLHS(p->subs[1], init, locRow, locCol, value, row, colOffset, coeff);
			} else if (p->op == OP_SUB) {
				insertLHS(p->subs[0], init, locRow, locCol, value, row, colOffset, coeff);
				insertLHS(p->subs[1], init, locRow, locCol, value, row, colOffset, -coeff);
			} else if (p->op == OP_MUL) {
				if (p->subs[0] != NULL && p->subs[1] != NULL) {
					insertLHS(p->subs[0], init, locRow, locCol, value, row, colOffset, coeff * (p->subs[1]->op == OP_ID ? init(p->subs[1]->index) : p->subs[1]->value.d));
					insertLHS(p->subs[1], init, locRow, locCol, value, row, colOffset, coeff * (p->subs[0]->op == OP_ID ? init(p->subs[0]->index) : p->subs[0]->value.d));
				}
			} else if (p->op == OP_POS) {
				insertLHS(p->subs[0], init, locRow, locCol, value, row, colOffset, coeff);
			} else if (p->op == OP_MINUS) {
				insertLHS(p->subs[0], init, locRow, locCol, value, row, colOffset, -coeff);
			}
		}

		void updateDEcoeff(int lvl, int row, mat &sol, AstNode *p, double coeff) {
			if (p == NULL) {
				return;
			}
			if (p->op == OP_ID) {
				sol(row, lvl) += 1.0 / lvl * coeff * sol(p->index, lvl - 1);
			} else if (p->op == OP_ADD) {
				updateDEcoeff(lvl, row, sol, p->subs[0], coeff);
				updateDEcoeff(lvl, row, sol, p->subs[1], coeff);
			} else if (p->op == OP_SUB) {
				updateDEcoeff(lvl, row, sol, p->subs[0], coeff);
				updateDEcoeff(lvl, row, sol, p->subs[1], -coeff);
			} else if (p->op == OP_POS) {
				updateDEcoeff(lvl, row, sol, p->subs[0], coeff);
			} else if (p->op == OP_MINUS) {
				updateDEcoeff(lvl, row, sol, p->subs[0], -coeff);
			} else if (p->op == OP_MUL) {
				if (p->subs[0] != NULL && p->subs[1] != NULL) {
					if (p->subs[0]->op == OP_ID && p->subs[1]->op == OP_ID) {
						double d = 0.0;
						for (int i = 0; i < lvl; i++) {
							d += sol(p->subs[0]->index, i) * sol(p->subs[1]->index, lvl - 1 - i);
						}
						sol(row, lvl) += 1.0 / lvl * coeff * d;
					} else if (p->subs[0]->op == OP_ID && p->subs[1]->op == OP_CONST) {
						sol(row, lvl) += 1.0 / lvl * coeff * (p->subs[1]->value.d) * sol(p->subs[0]->index, lvl - 1);
					} else if (p->subs[0]->op == OP_CONST && p->subs[1]->op == OP_ID) {
						sol(row, lvl) += 1.0 / lvl * coeff * (p->subs[0]->value.d) * sol(p->subs[1]->index, lvl - 1);
					}
				}
			} else if (p->op == OP_CONST) {
				if (lvl == 1) {
					sol(row, lvl) += 1.0 / lvl * coeff * (p->value.d);
				}
			}
		}

		void updateAERHS(int lvl, int row, mat &sol, vec &rhs, AstNode *p, double coeff) {
			if (p == NULL) {
				return;
			}
			if (p->op == OP_ID) {
			} else if (p->op == OP_EQN) {
				if (p->subs[0] != NULL && p->subs[1] != NULL && p->subs[0]->op == OP_CALL) {
					if (p->subs[1]->op == OP_ID && p->subs[0]->subs[1]->op == OP_ID) {
						char *funcName = (char *)p->subs[0]->subs[0]->value.p;
						int tok = findDefaultFunctions(funcName, strlen(funcName));
						if (tok == TK_F_SIN) { // sin(x)=y;
							double d = 0.0;
							for (int k = 1; k <= lvl - 1; k++) {
								d += (lvl - k) * sol(p->subs[1]->pairIdx, k) * sol(p->subs[0]->subs[1]->index, lvl - k);
							}
							rhs(row) = -coeff * d / lvl;
						} else if (tok == TK_F_COS) { // cos(x)=y;
							double d = 0.0;
							for (int k = 1; k <= lvl - 1; k++) {
								d += (lvl - k) * sol(p->subs[1]->pairIdx, k) * sol(p->subs[0]->subs[1]->index, lvl - k);
							}
							rhs(row) = coeff * d / lvl;
						} else if (tok == TK_F_EXP) { // exp(x)=y;
							double d = 0.0;
							for (int k = 1; k <= lvl - 1; k++) {
								d += (lvl - k) * sol(p->subs[1]->index, k) * sol(p->subs[0]->subs[1]->index, lvl - k);
							}
							rhs(row) = -coeff * d / lvl;
						} else if (tok == TK_F_SQRT) { // sqrt(x)=y;
							double d = 0.0;
							for (int k = 1; k <= lvl - 1; k++) {
								d += sol(p->subs[1]->index, k) * sol(p->subs[1]->index, lvl - k);
							}
							rhs(row) = coeff * d / 2 / sol(p->subs[1]->index, 0);
						} else if (tok == TK_F_LOG || tok == TK_F_ATAN || tok == TK_F_ASIN || tok == TK_F_ACOS) {
							double sign = 1.0, alpha = 1.0;
							inverseFunction(tok, sign, alpha);
							const int xIdx = p->subs[0]->subs[1]->index;
							vector<double> x(lvl);
							for (int k = 0; k < lvl; k++) {
								x[k] = sol(xIdx, k);
							}
							vector<double> w = inverseFunctionFactor(tok, alpha, x, lvl);
							double d = 0.0;
							for (int k = 1; k <= lvl - 1; k++) {
								d += k * x[k] * w[lvl - k];
							}
							rhs(row) = -coeff * sign * d / lvl;
						}
						// There should not be der(x) because only AEs are processed.
					} // TODO: There may be cases such as sin(x)=0;
				} else if (p->subs[0] != NULL && p->subs[1] != NULL && p->subs[0]->op == OP_DIV) { // x/y=z;
					if (p->subs[0]->subs[0] != NULL && p->subs[0]->subs[1] != NULL) {
						if (p->subs[0]->subs[1]->op == OP_ID && p->subs[1]->op == OP_ID) {
							double d = 0.0;
							for (int k = 1; k <= lvl - 1; k++) {
								d += sol(p->subs[0]->subs[1]->index, k) * sol(p->subs[1]->index, lvl - k);
							}
							rhs(row) = coeff * d / sol(p->subs[0]->subs[1]->index, 0);
						}
					}
				} else {
					updateAERHS(lvl, row, sol, rhs, p->subs[0], -coeff);
					updateAERHS(lvl, row, sol, rhs, p->subs[1], coeff);
				}
			} else if (p->op == OP_ADD) {
				updateAERHS(lvl, row, sol, rhs, p->subs[0], coeff);
				updateAERHS(lvl, row, sol, rhs, p->subs[1], coeff);
			} else if (p->op == OP_SUB) {
				updateAERHS(lvl, row, sol, rhs, p->subs[0], coeff);
				updateAERHS(lvl, row, sol, rhs, p->subs[1], -coeff);
			} else if (p->op == OP_MUL) {
				if (p->subs[0] != NULL && p->subs[1] != NULL && p->subs[0]->op == OP_ID && p->subs[1]->op == OP_ID) {
					double d = 0.0;
					for (int k = 1; k <= lvl - 1; k++) {
						d += sol(p->subs[0]->index, k) * sol(p->subs[1]->index, lvl - k);
					}
					rhs(row) += coeff * d;
				}
			} else if (p->op == OP_POS) {
				updateAERHS(lvl, row, sol, rhs, p->subs[0], coeff);
			} else if (p->op == OP_MINUS) {
				updateAERHS(lvl, row, sol, rhs, p->subs[0], -coeff);
			} else if (p->op == OP_CONST) {
				// you do not have to deal with standalone const on AEs
			}
		}

		// Largest absolute entry; infinity when any entry is not finite (arma's "inf" norm skips NaN, which
		// would let a step with an undefined residual pass the tolerance check).
		static double maxAbsEntry(const vec &v, bool skipNonFinite = false) {
			double m = 0.0;
			for (uword i = 0; i < v.n_elem; i++) {
				if (!std::isfinite(v(i))) {
					if (skipNonFinite) {
						continue;
					}
					return std::numeric_limits<double>::infinity();
				}
				m = std::max(m, std::fabs(v(i)));
			}
			return m;
		}

		// Collects the places in a tree where an operand has a restricted range: arguments of log, sqrt,
		// asin, acos, and the denominators of divisions and of 1/x.
		static void collectDomainChecks(AstNode *p, vector<DomainCheck> &out) {
			if (p == NULL) {
				return;
			}
			if (p->op == OP_CALL && p->subs[0] != NULL && p->subs[1] != NULL && p->subs[0]->op == OP_ID) {
				char *funcName = (char *)p->subs[0]->value.p;
				TokenKind tk = (TokenKind)findDefaultFunctions(funcName, strlen(funcName));
				if (tk == TK_F_LOG) {
					out.push_back({p->subs[1], DomainCheck::POSITIVE, "log"});
				} else if (tk == TK_F_SQRT) {
					out.push_back({p->subs[1], DomainCheck::NONNEGATIVE, "sqrt"});
				} else if (tk == TK_F_ASIN || tk == TK_F_ACOS) {
					out.push_back({p->subs[1], DomainCheck::UNIT, tk == TK_F_ASIN ? "asin" : "acos"});
				} else if (tk == TK_F_INV) {
					out.push_back({p->subs[1], DomainCheck::NONZERO, "1/x"});
				}
			} else if (p->op == OP_DIV && p->subs[1] != NULL) {
				out.push_back({p->subs[1], DomainCheck::NONZERO, "division"});
			}
			collectDomainChecks(p->subs[0], out);
			collectDomainChecks(p->subs[1], out);
			collectDomainChecks(p->next, out);
		}

		static bool domainValueBad(const DomainCheck &c, double x) {
			switch (c.kind) {
			case DomainCheck::POSITIVE:
				return !(x > 0.0);
			case DomainCheck::NONNEGATIVE:
				return !(x >= 0.0);
			case DomainCheck::UNIT:
				return !(x >= -1.0 && x <= 1.0);
			default:
				return !(x != 0.0) || !std::isfinite(x);
			}
		}

		static string domainMessage(const DomainCheck &c, double x) {
			static const char *need[] = {"must be positive", "must not be negative", "must be in [-1, 1]", "must not be zero"};
			ostringstream os;
			os << c.label << ": the " << (strcmp(c.label, "division") == 0 ? "denominator" : "argument") << " is " << x << " but it " << need[c.kind];
			return os.str();
		}

		// Checks the domain conditions along [0, alpha] of a solution at a few points. A denominator that
		// changes sign between two points has a pole in between, which a series can pass without any
		// residual (a rational function is represented exactly).
		bool SasComputationModel::domainOk(CheSolution *solution, double alpha, string *why, double *at) {
			if (domainChecks.empty()) {
				return true;
			}
			const int samples = 8;
			vector<double> prev(domainChecks.size(), 0.0);
			list<shared_ptr<SasModel>> tempList;
			IdLocator locator(tempList);
			for (int k = 0; k <= samples; k++) {
				vec state = solution->getSolValue(alpha * k / samples);
				vec der(state.n_rows, fill::zeros);
				locator.setStates(state, der);
				for (size_t i = 0; i < domainChecks.size(); i++) {
					const DomainCheck &c = domainChecks[i];
					double x = calcTreeValue(c.arg, &locator);
					// k == 0 is the start point (checked by the caller), used as the reference for sign changes
					bool bad = k > 0 && domainValueBad(c, x);
					bool crossed = k > 0 && c.kind == DomainCheck::NONZERO && std::isfinite(x) && ((x < 0.0) != (prev[i] < 0.0));
					if (bad || crossed) {
						if (at != NULL) {
							*at = alpha * k / samples;
						}
						if (why != NULL) {
							*why = bad ? domainMessage(c, x) : string(c.label) + ": the denominator changes sign (from " + to_string(prev[i]) + " to " + to_string(x) + "), the solution has a pole";
						}
						return false;
					}
					prev[i] = x;
				}
			}
			return true;
		}

		bool SasComputationModel::validateDomain(const SasSolution &sol, string *why) {
			const double len = sol.absEnd - sol.absStart;
			if (len <= 0.0) {
				return true;
			}
			if (domainChecks.empty()) {
				for (const shared_ptr<AstTree> &t : eqnTable) {
					collectDomainChecks(t->pHead, domainChecks);
				}
			}
			string msg;
			double at = 0.0;
			if (domainOk(sol.solution, 0.98 * len, &msg, &at)) {
				return true;
			}
			if (why != NULL) {
				*why = msg + " at t = " + to_string(sol.absStart + at);
			}
			return false;
		}

		// Describes why the equations cannot be evaluated at the given point of a solution ("" if no
		// function argument is outside its domain).
		string SasComputationModel::describeDomainViolation(CheSolution *solution, double alpha) {
			vec state = solution->getSolValue(alpha);
			vec der(state.n_rows, fill::zeros);
			list<shared_ptr<SasModel>> tempList;
			IdLocator locator(tempList);
			locator.setStates(state, der);
			for (const DomainCheck &c : domainChecks) {
				double x = calcTreeValue(c.arg, &locator);
				if (domainValueBad(c, x)) {
					return domainMessage(c, x);
				}
			}
			return "";
		}

		vec SasComputationModel::calcDiff(CheSolution *solution, double alpha, CheSolution *derSol) {
			vec state = solution->getSolValue(alpha);
			std::unique_ptr<CheSolution> ownDerSol;
			if (derSol == nullptr) {
				ownDerSol.reset(CheSolutionFactory::makeDerivativeSol(solution));
				derSol = ownDerSol.get();
			}
			vec der = derSol->getSolValue(alpha);
			return calcDiffAt(state, der);
		}

		vec SasComputationModel::calcDiffAt(vec &state, vec &der) {
			vec diff(nDE + nAE, fill::zeros);
			list<shared_ptr<SasModel>> tempList;
			IdLocator locator(tempList);
			/*state.print("state");
			der.print("der");*/
			locator.setStates(state, der);
			for (int i = 0; i < nDE + nAE; i++) {
				diff(i) = calcTreeValue(eqnTable[i]->pHead, &locator);
			}

			return diff;
		}

		shared_ptr<SasSolution> SasComputationModel::solveSegment(const vec &init, double curAlpha, double seg, const SasComputationOptions &options) {
			double startAbsAlpha = curAlpha;
			double maxAlpha = seg;
			double tol = options.errorTol;

			int nState = init.n_rows;

			shared_ptr<SasSolution> sasSol = make_shared<SasSolution>(this, CHESOL_PADE, nState, options.nLvl);
			sasSol->solution->solution.col(0) = init;

			// Jacobian of the algebraic equations at the given state, split into the parts for the
			// differential (X) and algebraic (Y) variables.
			sp_mat LHSY, LHSX;
			auto buildLHS = [&](const vec &state) {
				vector<int> locRow;
				vector<int> locCol;
				vector<double> value;
				for (int iAE = 0; iAE < nAE; iAE++) {
					insertLHS(eqnTable[nDE + iAE]->pHead, state, locRow, locCol, value, iAE, 0, 1.0);
				}
				uvec locRowVec = conv_to<uvec>::from(locRow);
				uvec locColVec = conv_to<uvec>::from(locCol);
				vec valueVec = conv_to<vec>::from(value);
				sp_mat LHStotal(true, join_cols(locRowVec.t(), locColVec.t()), valueVec, nAE, nX + nY);
				LHSY = LHStotal.tail_cols(nY);
				LHSX = LHStotal.head_cols(nX);
			};
			buildLHS(init);

			// A Pade solution caches its coefficients on the first evaluation. The series is changed after
			// that (AE calibration, filling the higher orders), so the cache must be dropped each time.
			auto resetPade = [&sasSol]() {
				if (CheSolutionPade *pade = dynamic_cast<CheSolutionPade *>(sasSol->solution)) {
					pade->ready = false;
				}
			};

			// The residual of the equations at the initial state. Only that state is in the series so far, so
			// the solution and its derivative at the start are the initial state and zero, with no need to
			// build a Pade approximant of either.
			auto residualAtStart = [&]() {
				vec state = sasSol->solution->solution.col(0);
				vec der(state.n_rows, fill::zeros);
				return calcDiffAt(state, der);
			};

			if (nAE > 0) {
				// Check and calibrate AE imbalances.
				vec diff = residualAtStart();
				vec initY = sasSol->solution->solution.col(0).tail(nY);
				vec diffAE = diff.tail(nAE);
				double diffMax = maxAbsEntry(diffAE);
				if (diffMax > tol / 2) {
					// Newton's method; the Jacobian is rebuilt at each iterate.
					int iter = 0;
					int maxIter = 10;
					bool singular = false;
					while (diffMax > tol / 10 && iter < maxIter) {
						try {
							buildLHS(sasSol->solution->solution.col(0));
							sasSol->solution->solution.col(0).tail(nY) -= spsolve(LHSY, diffAE);
						} catch (const std::runtime_error &) {
							singular = true;
							break;
						}
						diff = residualAtStart();
						diffAE = diff.tail(nAE);
						diffMax = maxAbsEntry(diffAE);
						iter++;
					}
					if (singular || (diffMax > tol / 10 && iter >= maxIter) || !std::isfinite(diffMax)) { // if not good, then restore original solution.
						sasSol->solution->solution.col(0).tail(nY) = initY;
						const double before = maxAbsEntry(residualAtStart().tail(nAE));
						if (!(before < 1e-4)) {
							string why = describeDomainViolation(sasSol->solution, 0.0);
							if (!why.empty()) {
								throw std::runtime_error("the equations cannot be evaluated at t = " + std::to_string(startAbsAlpha) + ": " + why);
							}
							throw std::runtime_error("the algebraic equations cannot be solved consistently at t = " + std::to_string(startAbsAlpha) + " (residual " + std::to_string(before) + "); give the algebraic variables start values that are close to the solution");
						}
						if (options.log != nullptr) {
							*options.log << "Refinement of AE solution fails, restoring original solution" << endl;
						}
					}
					// the series is expanded around the refined point
					buildLHS(sasSol->solution->solution.col(0));
				}
			}

			// The matrix of the algebraic variables is the same at every order. For a few unknowns SuperLU's
			// set-up costs more than the solve, so a small system is factored once, densely, and each order is
			// two triangular solves.
			const bool denseAE = nAE > 0 && nAE == nY && nY <= 64;
			mat luL, luU, luP;
			if (denseAE) {
				bool singular = !lu(luL, luU, luP, mat(LHSY));
				for (int i = 0; i < nY && !singular; i++) {
					singular = !(std::fabs(luU(i, i)) > 0.0) || !std::isfinite(luU(i, i));
				}
				if (singular) {
					throw std::runtime_error("the matrix of the algebraic equations is singular");
				}
			}

			for (int lvl = 1; lvl < options.nLvl; lvl++) {
				for (int iDE = 0; iDE < nDE; iDE++) {
					updateDEcoeff(lvl, iDE, sasSol->solution->solution, eqnTable[iDE]->pHead->subs[1], 1.0);
				}

				if (nAE > 0) {
					vec rhs(nAE, fill::zeros);
					for (int iAE = 0; iAE < nAE; iAE++) {
						updateAERHS(lvl, iAE, sasSol->solution->solution, rhs, eqnTable[nDE + iAE]->pHead, 1.0);
					}

					rhs -= LHSX * sasSol->solution->solution.col(lvl).head(nX);

					if (denseAE) {
						// P*A = L*U: forward substitution with the unit lower factor, then back substitution
						const vec pr = luP * rhs;
						vec z(nY);
						for (int i = 0; i < nY; i++) {
							double acc = pr(i);
							for (int j = 0; j < i; j++) {
								acc -= luL(i, j) * z(j);
							}
							z(i) = acc;
						}
						for (int i = nY - 1; i >= 0; i--) {
							double acc = z(i);
							for (int j = i + 1; j < nY; j++) {
								acc -= luU(i, j) * z(j);
							}
							z(i) = acc / luU(i, i);
						}
						sasSol->solution->solution.col(lvl).tail(nY) = z;
					} else {
						sasSol->solution->solution.col(lvl).tail(nY) = spsolve(LHSY, rhs, "superlu");
					}
				}
			}

			resetPade();

			// sasSol->solution->solution.print("sol");
			// The series is final from here on: the derivative, and its Pade approximant, are built once
			// and used by every trial of the step size search.
			std::unique_ptr<CheSolution> derSol(CheSolutionFactory::makeDerivativeSol(sasSol->solution));
			vec diff = calcDiff(sasSol->solution, 0.0, derSol.get());
			double diffMax = maxAbsEntry(diff);
			if (!std::isfinite(diffMax)) {
				string why = describeDomainViolation(sasSol->solution, 0.0);
				throw std::runtime_error("the equations cannot be evaluated at t = " + std::to_string(startAbsAlpha) + (why.empty() ? string("") : ": " + why));
			}
			if (diffMax > tol) {
				if (solverWarnings.size() < MAX_SOLVER_WARNINGS) {
					ostringstream os;
					os << "at t = " << startAbsAlpha << " the residual at the start of a step (" << diffMax << ") exceeds the error tolerance (" << tol << "); the tolerance was relaxed for that step";
					solverWarnings.push_back(os.str());
				}
				tol = diffMax;
			}
			double alphaLeft = 0.0;
			double alphaRight = maxAlpha;
			double alpha = alphaRight;
			// The search resolution is alphaTol, but a short segment (the last one of a run, or one up to
			// an event) is resolved finer so that a run can end exactly at the requested time.
			const double resolution = std::min(options.alphaTol, maxAlpha / 64.0);
			double rejectAlpha = -1.0; // the smallest step rejected because the residual is not finite
			string rejectWhy;		   // the reason of the smallest step rejected because of a domain violation
			while (alphaRight - alphaLeft > resolution) {
				diff = calcDiff(sasSol->solution, alpha, derSol.get());
				// diff.print("diff");
				diffMax = maxAbsEntry(diff, options.tolerateNonFinite);
				string why;
				if (diffMax < tol && (options.tolerateNonFinite || domainOk(sasSol->solution, alpha, &why))) {
					alphaLeft = alpha;
					alpha = (alphaLeft + alphaRight) / 2;
				} else {
					if (!why.empty()) {
						rejectWhy = why;
					} else if (!std::isfinite(diffMax)) {
						rejectAlpha = alpha;
					}
					alphaRight = alpha;
					alpha = (alphaLeft + alphaRight) / 2;
				}
			}
			alpha = alphaLeft;
			if (alpha <= 0.0) {
				// The step cannot advance because the equations are undefined right after this point.
				if (!rejectWhy.empty()) {
					stallReason = rejectWhy;
				} else if (rejectAlpha >= 0.0) {
					stallReason = describeDomainViolation(sasSol->solution, rejectAlpha);
				}
			}
			if (options.log != nullptr) {
				*options.log << "Stage start " << startAbsAlpha << " end " << startAbsAlpha + alpha << ", length " << alpha << endl;
			}

			sasSol->absStart = startAbsAlpha;
			sasSol->absEnd = startAbsAlpha + alpha;

			return sasSol;
		}

		SasSolutionSet *SasComputationModel::solve(const SasComputationOptions &options) {
			vec init = getStartStateIVP();
			if (options.log != nullptr) {
				*options.log << "Init values generated." << endl;
			}
			SasSolutionSet *solLink = new SasSolutionSet();
			double alpha = options.startTime;
			int noMoveCnt = 0;
			// Runs until the end time is reached (up to rounding), not just to within alphaTol.
			const double timeEps = 1e-12 * std::max(1.0, std::fabs(options.endTime));
			solverWarnings.clear();
			stallReason.clear();
			domainChecks.clear();
			for (const shared_ptr<AstTree> &t : eqnTable) {
				collectDomainChecks(t->pHead, domainChecks);
			}
			while (alpha < options.endTime - timeEps) {
				double seg = options.endTime - alpha;
				if (seg > options.segment) {
					seg = options.segment;
				}
				shared_ptr<SasSolution> sasSol = solveSegment(init, alpha, seg, options);
				if (sasSol->absEnd >= options.endTime - timeEps) {
					sasSol->absEnd = options.endTime; // land exactly on the end time
				}
				if (sasSol->absEnd <= alpha) {
					noMoveCnt++;
				} else {
					noMoveCnt = 0;
				}
				if (noMoveCnt >= 2) {
					solLink->stalled = true;
					solLink->stallTime = alpha;
					ostringstream os;
					os << "the solver cannot advance beyond t = " << alpha;
					if (!stallReason.empty()) {
						os << ": " << stallReason;
					} else {
						os << " within the error tolerance; the solution may be singular here";
					}
					solLink->stallMessage = os.str();
					break;
				}
				alpha = sasSol->absEnd;
				init = sasSol->solution->getSolValue(sasSol->absEnd - sasSol->absStart);
				solLink->solutionLink.push_back(sasSol);
				if (options.progress && options.progress(alpha, options.endTime)) {
					solLink->cancelled = true;
					break;
				}
				if (options.onSolution && options.onSolution(sasSol)) {
					break;
				}
			}

			solLink->warnings.insert(solLink->warnings.end(), solverWarnings.begin(), solverWarnings.end());
			return solLink;
		}

		SasComputationModel::~SasComputationModel() {
			strTable.clear();
			topModels.clear();
			eqnTable.clear();
			idTable.clear();
		}

		shared_ptr<IdPearl> IdLocator::findId(list<shared_ptr<SasModel>> &searchDomain, const char *iden, int len) {
			shared_ptr<IdPearl> id = nullptr;
			list<shared_ptr<SasModel>>::iterator searchIt;
			for (searchIt = searchDomain.begin(); searchIt != searchDomain.end(); searchIt++) {
				if ((id = (*searchIt)->lookupIdFromModel(iden, len)) != nullptr) {
					return id;
				}
			}
			return id;
		}

		SasSolution::SasSolution(SasComputationModel *linkCompModel, int type, int nState, int nLvl) {
			this->linkCompModel = linkCompModel;
			this->solution = CheSolutionFactory::makeInitCheSol(type, nState, nLvl);
		}

		SasSolution::~SasSolution() {
			if (solution != NULL) {
				delete solution;
				solution = NULL;
			}
		}

		SasSolutionSet::SasSolutionSet() : solutionLink() {
		}

		list<shared_ptr<SasSolution>> SasSolutionSet::findSasSolution(double time) {
			list<shared_ptr<SasSolution>> solList;
			list<shared_ptr<SasSolution>>::iterator solIt;
			for (solIt = solutionLink.begin(); solIt != solutionLink.end(); solIt++) {
				if ((*solIt)->absStart <= time && (*solIt)->absEnd >= time) {
					solList.push_back(*solIt);
				}
			}
			return solList;
		}

		shared_ptr<SasSolution> SasSolutionSet::findSasSolution(shared_ptr<SasComputationModel> &pCompModel, double time) {
			list<shared_ptr<SasSolution>>::iterator solIt;
			for (solIt = solutionLink.begin(); solIt != solutionLink.end(); solIt++) {
				if ((*solIt)->absStart <= time && (*solIt)->absEnd >= time && pCompModel.get() == (*solIt)->linkCompModel) {
					return (*solIt);
				}
			}
			return nullptr;
		}

		SasSolutionSet::~SasSolutionSet() {
			solutionLink.clear();
		}

		IdLocator::IdLocator(const std::list<std::shared_ptr<SasModel>> &searchDomain) : searchDomain(searchDomain) {}

		void IdLocator::setStates(arma::vec &state, arma::vec &der) {
			type = LOCATOR_STATE;
			this->state = state;
			this->der = der;
		}

		void IdLocator::unsetStates() {
			type = LOCATOR_INIT;
		}

		double IdLocator::findVarValue(const char *iden, int len) {
			shared_ptr<IdPearl> id = IdLocator::findId(this->searchDomain, iden, len);
			if (id == nullptr) {
				return NAN;
			}
			if (this->type == LOCATOR_INIT) {
				return id->value.d;
			} else if (this->type == LOCATOR_STATE) {
				return state(id->index);
			} else {
				return NAN;
			}
		}

		double IdLocator::findDerValue(const char *iden, int len) {
			shared_ptr<IdPearl> id = IdLocator::findId(this->searchDomain, iden, len);
			if (id == nullptr) {
				return NAN;
			}
			if (this->type == LOCATOR_INIT) {
				return id->derVal.d;
			} else if (this->type == LOCATOR_STATE) {
				return der(id->index);
			} else {
				return NAN;
			}
		}
	} // namespace core
} // namespace che