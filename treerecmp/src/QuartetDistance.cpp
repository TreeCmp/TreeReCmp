
//
// File: QPartitionList.h
// Created by: Anna Pawelczyk
// Created on: 15 Aug 2011, 11:10
//

/*
This software is a computer program whose purpose is to provide classes
for phylogenetic data analysis.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "QuartetDistance.h"
using namespace bpp;
using namespace std;


namespace tools {

	TreeParams2::TreeParams2(Node* root, int n, int l)
	{
		inSize = n - l;
		lSize = l;
		rootId = root->getId() - l;

		subTr = new int[inSize * 2];
		countSubTrSize(root);

		inSons = new list<int>[inSize];
		setInSons(root, rootId);
	}
	TreeParams2::~TreeParams2()
	{
		delete[] subTr;
		delete[] inSons;
	}

	int TreeParams2::choose2(int a)
	{
		return a * (a - 1) / 2;
	}
	int TreeParams2::choose3(int a)
	{
		return a * (a - 1) * (a - 2) / (2 * 3);
	}

	int TreeParams2::countSubTrSize(Node* r)
	{
		vector<Node*> s = r->getSons();
		if (s.size() == 0) {
			return 1;           // a leaf
		}
		int id = r->getId() - lSize;
		subTr[id] = 0;
		for (vector<Node*>::iterator it = s.begin(); it != s.end(); it++) {
			subTr[id] += countSubTrSize(*it);
		}
		return subTr[id];
	}

	void TreeParams2::setInSons(Node* root, int rootNewId)
	{
		vector<Node*> sons = root->getSons();
		for (vector<Node*>::iterator sonIt = sons.begin(); sonIt != sons.end(); sonIt++) {
			int sonId = (*sonIt)->getId() - lSize;
			if (sonId > 0) {
				inSons[rootNewId].push_back(sonId);
				setInSons(*sonIt, sonId);
			}
		}
	}


	QuartetDistance::QuartetDistance(const TreeTemplate<Node>& tr1In, const TreeTemplate<Node>& tr2In)
	{
		lSize = tr1In.getNumberOfLeaves();
		r1 = const_cast<Node*> (tr1In.getRootNode());
		r2 = const_cast<Node*> (tr2In.getRootNode());
		trP1 = new TreeParams2(r1, tr1In.getNumberOfNodes(), lSize);
		trP2 = new TreeParams2(r2, tr2In.getNumberOfNodes(), lSize);

		intersection = new int*[trP1->inSize];
		for (int i = 0; i < trP1->inSize; i++) {
			intersection[i] = new int[trP2->inSize];
			for (int j = 0; j < trP2->inSize; j++) {
				intersection[i][j] = -1;
			}
		}
		countIntersection(r1, r2);
	}
	QuartetDistance::~QuartetDistance()
	{
		for (int i = 0; i < trP1->inSize; i++) {
			delete[] intersection[i];
		}
		delete[] intersection;
		delete trP1;
		delete trP2;
	}

	int QuartetDistance::getDistance()
	{
		int B1 = getSingleTrQuartetsSize(trP1->inSize, trP1->subTr, trP1->inSons);
		int B2 = getSingleTrQuartetsSize(trP2->inSize, trP2->subTr, trP2->inSons);
		int S = getShared();
		int N = getNonshared();
		return B1 + B2 - 2 * S - N;
	}

	int QuartetDistance::getNonshared()
	{
		int nonshared = 0;

		for (int v1 = 0; v1 < trP1->inSize; v1++) {
			list<int> I = trP1->inSons[v1];
			for (int v2 = 0; v2 < trP2->inSize; v2++) {
				list<int> J = trP2->inSons[v2];

				for (list<int>::iterator it1 = I.begin(); it1 != I.end(); it1++) {
					for (list<int>::iterator it2 = J.begin(); it2 != J.end(); it2++) {
						int i = *it1, j = *it2;
						nonshared += in_a_b(i, j) * in_a_bneg(i, j) * in_aneg_b(i, j) * in_aneg_bneg(i, j);

						for (list<int>::iterator it1_k = I.begin(); it1_k != I.end(); it1_k++) {
							if (it1 != it1_k) {
								nonshared -= in_a_b(i, j) * in_a_bneg(i, j) * in_a_b(*it1_k, j) * in_a_bneg(*it1_k, j);
							}
						}
						nonshared -= in_a_b(i, j) * in_a_bneg(i, j) * in_aneg_b(v1, j) * in_aneg_bneg(v1, j);

						for (list<int>::iterator it2_l = J.begin(); it2_l != J.end(); it2_l++) {
							if (it2 != it2_l) {
								nonshared -= in_a_b(i, j) * in_aneg_b(i, j) * in_a_b(i, *it2_l) * in_aneg_b(i, *it2_l);
							}
						}
						nonshared -= in_a_b(i, j) * in_aneg_b(i, j) * in_a_bneg(i, v2) * in_aneg_bneg(i, v2);

						for (list<int>::iterator it1_k = I.begin(); it1_k != I.end(); it1_k++) {
							for (list<int>::iterator it2_l = J.begin(); it2_l != J.end(); it2_l++) {
								if (it1 != it1_k && it2 != it2_l) {
									nonshared += in_a_b(i, j) * in_a_b(i, *it2_l) * in_a_b(*it1_k, j) * in_a_b(*it1_k, *it2_l);
								}
							}
							nonshared += in_a_b(i, j) * in_a_bneg(i, v2) * in_a_b(*it1_k, j) * in_a_bneg(*it1_k, v2);
						}
						for (list<int>::iterator it2_l = J.begin(); it2_l != J.end(); it2_l++) {
							nonshared += in_a_b(i, j) * in_a_b(i, *it2_l) * in_aneg_b(v1, j) * in_aneg_b(v1, *it2_l);
						}
						nonshared += in_a_b(i, j) * in_a_bneg(i, v2) * in_aneg_b(v1, j) * in_aneg_bneg(v1, v2);
					}
				}

				for (list<int>::iterator it2 = J.begin(); it2 != J.end(); it2++) {
					int j = *it2;
					nonshared += in_aneg_b(v1, j) * in_aneg_bneg(v1, j) * in_a_b(v1, j) * in_a_bneg(v1, j);

					for (list<int>::iterator it1_k = I.begin(); it1_k != I.end(); it1_k++) {
						nonshared -= in_aneg_b(v1, j) * in_aneg_bneg(v1, j) * in_a_b(*it1_k, j) * in_a_bneg(*it1_k, j);
					}
					for (list<int>::iterator it2_l = J.begin(); it2_l != J.end(); it2_l++) {
						if (it2 != it2_l) {
							nonshared -= in_aneg_b(v1, j) * in_a_b(v1, j) * in_aneg_b(v1, *it2_l) * in_a_b(v1, *it2_l);
						}
					}
					nonshared -= in_aneg_b(v1, j) * in_a_b(v1, j) * in_aneg_bneg(v1, v2) * in_a_bneg(v1, v2);

					for (list<int>::iterator it1_k = I.begin(); it1_k != I.end(); it1_k++) {
						for (list<int>::iterator it2_l = J.begin(); it2_l != J.end(); it2_l++) {
							if (it2 != it2_l) {
								nonshared += in_aneg_b(v1, j) * in_aneg_b(v1, *it2_l) * in_a_b(*it1_k, j) * in_a_b(*it1_k, *it2_l);
							}
						}
						nonshared += in_aneg_b(v1, j) * in_aneg_bneg(v1, v2) * in_a_b(*it1_k, j) * in_a_bneg(*it1_k, v2);
					}
				}

				for (list<int>::iterator it1 = I.begin(); it1 != I.end(); it1++) {
					int i = *it1;
					nonshared += in_a_bneg(i, v2) * in_a_b(i, v2) * in_aneg_bneg(i, v2) * in_aneg_b(i, v2);

					for (list<int>::iterator it1_k = I.begin(); it1_k != I.end(); it1_k++) {
						if (it1 != it1_k) {
							nonshared -= in_a_bneg(i, v2) * in_a_b(i, v2) * in_a_bneg(*it1_k, v2) * in_a_b(*it1_k, v2);
						}
					}
					nonshared -= in_a_bneg(i, v2) * in_a_b(i, v2) * in_aneg_bneg(v1, v2) * in_aneg_b(v1, v2);

					for (list<int>::iterator it2_l = J.begin(); it2_l != J.end(); it2_l++) {
						nonshared -= in_a_bneg(i, v2) * in_aneg_bneg(i, v2) * in_a_b(i, *it2_l) * in_aneg_b(i, *it2_l);
					}

					for (list<int>::iterator it2_l = J.begin(); it2_l != J.end(); it2_l++) {
						for (list<int>::iterator it1_k = I.begin(); it1_k != I.end(); it1_k++) {
							if (it1 != it1_k) {
								nonshared += in_a_bneg(i, v2) * in_a_b(i, *it2_l) * in_a_bneg(*it1_k, v2) * in_a_b(*it1_k, *it2_l);
							}
						}
						nonshared += in_a_bneg(i, v2) * in_a_b(i, *it2_l) * in_aneg_bneg(v1, v2) * in_aneg_b(v1, *it2_l);
					}
				}


				nonshared += in_a_b(v1, v2) * in_a_bneg(v1, v2) * in_aneg_b(v1, v2) * in_aneg_bneg(v1, v2);

				for (list<int>::iterator it1_k = I.begin(); it1_k != I.end(); it1_k++) {
					nonshared -= in_aneg_bneg(v1, v2) * in_aneg_b(v1, v2) * in_a_bneg(*it1_k, v2) * in_a_b(*it1_k, v2);
				}

				for (list<int>::iterator it2_l = J.begin(); it2_l != J.end(); it2_l++) {
					nonshared -= in_aneg_bneg(v1, v2) * in_a_bneg(v1, v2) * in_aneg_b(v1, *it2_l) * in_a_b(v1, *it2_l);
				}

				for (list<int>::iterator it1_k = I.begin(); it1_k != I.end(); it1_k++) {
					for (list<int>::iterator it2_l = J.begin(); it2_l != J.end(); it2_l++) {
						nonshared += in_aneg_bneg(v1, v2) * in_aneg_b(v1, *it2_l) * in_a_bneg(*it1_k, v2) * in_a_b(*it1_k, *it2_l);
					}
				}
			}
		}

		return nonshared / 4;
	}

	int QuartetDistance::getShared()
	{
		int shared = 0;

		for (int v1 = 0; v1 < trP1->inSize; v1++) {
			list<int> I = trP1->inSons[v1];                 // SPR CZY 0
			for (int v2 = 0; v2 < trP2->inSize; v2++) {
				list<int> J = trP2->inSons[v2];

				list<int>::iterator it1;
				list<int>::iterator it2;
				//todo: W Visualu deklaracja tablic wymaga sta³ych wartoœci
				int* S1 = new int[trP1->inSize];
				int* S1_neg = new int[trP1->inSize];
				int* S2 = new int[trP2->inSize];
				int* S2_neg = new int[trP2->inSize];
				//int S1[trP1->inSize], S1_neg[trP1->inSize];
				//int S2[trP2->inSize], S2_neg[trP2->inSize];
				int minSize = min(trP1->inSize, trP2->inSize);
				for (int i = 0; i < minSize; i++) { S1[i] = 0; S1_neg[i] = 0; S2[i] = 0; S2_neg[i] = 0; }
				for (int i = minSize; i < trP1->inSize; i++) { S1[i] = 0; S1_neg[i] = 0; }
				for (int i = minSize; i < trP2->inSize; i++) { S2[i] = 0; S2_neg[i] = 0; }
				int S = 0;

				//******** S *****************************            
				for (it1 = I.begin(); it1 != I.end(); it1++) {
					int i = *it1;
					for (it2 = J.begin(); it2 != J.end(); it2++) {
						int j = *it2;
						S1[i] += a_b(i, j);
						S1_neg[i] += aneg_b(i, j);

						S += a_b(i, j);
					}
					S1[i] += a_bneg(i, v2);
					S1_neg[i] += aneg_bneg(i, v2);

					S2[v2] += a_bneg(i, v2);
					S2_neg[v2] += a_b(i, v2);

					S += a_bneg(i, v2);
				}


				for (it2 = J.begin(); it2 != J.end(); it2++) {
					int j = *it2;
					for (it1 = I.begin(); it1 != I.end(); it1++) {
						int i = *it1;
						S2[j] += a_b(i, j);
						S2_neg[j] += a_bneg(i, j);
					}
					S2[j] += aneg_b(v1, j);
					S2_neg[j] += aneg_bneg(v1, j);

					S1[v1] += aneg_b(v1, j);
					S1_neg[v1] += a_b(v1, j);

					S += aneg_b(v1, j);
				}

				S1[v1] += aneg_bneg(v1, v2);
				S1_neg[v1] += a_bneg(v1, v2);
				S2[v2] += aneg_bneg(v1, v2);
				S2_neg[v2] += aneg_b(v1, v2);

				S += aneg_bneg(v1, v2);

				//******** EQUATION *****************************
				for (it1 = I.begin(); it1 != I.end(); it1++) {
					int i = *it1;
					for (it2 = J.begin(); it2 != J.end(); it2++) {
						int j = *it2;
						int val =
							a_b(i, j) * (
								aneg_bneg(i, j)
								+ (a_bneg(i, j) - S2_neg[j])
								+ (aneg_b(i, j) - S1_neg[i])
								+ (S - S1[i] - S2[j] + a_b(i, j))
								);
						shared += val;
					}
				}

				if (v2 != r2->getId() - lSize) {
					for (it1 = I.begin(); it1 != I.end(); it1++) {
						int i = *it1;
						int val =
							a_bneg(i, v2) * (
								aneg_b(i, v2)
								+ (a_b(i, v2) - S2_neg[v2])
								+ (aneg_bneg(i, v2) - S1_neg[i])
								+ (S - S1[i] - S2[v2] + a_bneg(i, v2))
								);
						shared += val;
					}
				}
				if (v1 != r1->getId() - lSize) {
					for (it1 = J.begin(); it1 != J.end(); it1++) {
						int j = *it1;
						int val =
							aneg_b(v1, j) * (
								a_bneg(v1, j)
								+ (aneg_bneg(v1, j) - S2_neg[j])
								+ (a_b(v1, j) - S1_neg[v1])
								+ (S - S1[v1] - S2[j] + aneg_b(v1, j))
								);
						shared += val;
					}
				}
				if (v1 != r1->getId() - lSize && v2 != r2->getId() - lSize) {
					int val =
						aneg_bneg(v1, v2) * (
							a_b(v1, v2)
							+ (aneg_b(v1, v2) - S2_neg[v2])
							+ (a_bneg(v1, v2) - S1_neg[v1])
							+ (S - S1[v1] - S2[v2] + aneg_bneg(v1, v2))
							);
					shared += val;
				}
				// todo: tutaj musimy skasowaæ tablice
				delete[] S1;
				delete[] S1_neg;
				delete[] S2;
				delete[] S2_neg;
			}
		}

		return shared / 2;
	}

	int QuartetDistance::getSingleTrQuartetsSize(int size, int* subTrsize, list<int>* inSons)
	{
		int single = 0;
		for (int v = 0; v < size; v++) {
			int Ss = 0;
			list<int> &I = inSons[v];
			for (list<int>::iterator it = I.begin(); it != I.end(); it++) {
				Ss += binomCoef_bin(subTrsize[*it]);
			}
			Ss += binomCoef_bin(lSize - subTrsize[v]);


			for (list<int>::iterator it = I.begin(); it != I.end(); it++) {
				single += binomCoef_bin(subTrsize[*it]) * (binomCoef_bin(lSize - subTrsize[*it]) - Ss + binomCoef_bin(subTrsize[*it]));
			}
			single += binomCoef_bin(lSize - subTrsize[v]) * (binomCoef_bin(subTrsize[v]) - Ss + binomCoef_bin(lSize - subTrsize[v]));
		}
		return single / 2;
	}

	int QuartetDistance::binomCoef_bin(int a)
	{
		return a * (a - 1) / 2;
	}


	int QuartetDistance::non_product(int i, int j)
	{
		return a_b(i, j) * a_bneg(i, j) * aneg_b(i, j) * aneg_bneg(i, j);

	}

	int QuartetDistance::a_b(int i, int j)
	{
		return binomCoef_n_2(intersection[i][j]);
	}

	int QuartetDistance::aneg_bneg(int i, int j)
	{
		return binomCoef_n_2(lSize - (trP1->subTr[i] + trP2->subTr[j] - intersection[i][j]));
	}

	int QuartetDistance::a_bneg(int i, int j)
	{
		return binomCoef_n_2(trP1->subTr[i] - intersection[i][j]);
	}

	int QuartetDistance::aneg_b(int i, int j)
	{
		return binomCoef_n_2(trP2->subTr[j] - intersection[i][j]);
	}

	int QuartetDistance::in_a_b(int i, int j)
	{
		return (intersection[i][j]);
	}

	int QuartetDistance::in_aneg_bneg(int i, int j)
	{
		return (lSize - (trP1->subTr[i] + trP2->subTr[j] - intersection[i][j]));
	}

	int QuartetDistance::in_a_bneg(int i, int j)
	{
		return (trP1->subTr[i] - intersection[i][j]);
	}

	int QuartetDistance::in_aneg_b(int i, int j)
	{
		return (trP2->subTr[j] - intersection[i][j]);
	}

	int QuartetDistance::countIntersection(Node* root1, Node* root2)
	{
		int id1 = root1->getId() - lSize;
		int id2 = root2->getId() - lSize;
		bool internal = (min(id1, id2) >= 0) ? true : false;
		if (internal && intersection[id1][id2] != -1) return intersection[id1][id2];

		int val = 0;
		vector<Node*> s1 = root1->getSons();
		vector<Node*> s2 = root2->getSons();

		if (s1.size() == 0 && s2.size() == 0) {
			return id1 == id2;
		}

		for (vector<Node*>::iterator it1 = s1.begin(); it1 != s1.end(); it1++) {
			val += countIntersection(*it1, root2);
			for (vector<Node*>::iterator it2 = s2.begin(); it2 != s2.end(); it2++) {
				val -= countIntersection(*it1, *it2);
			}
		}
		for (vector<Node*>::iterator it2 = s2.begin(); it2 != s2.end(); it2++) {
			val += countIntersection(root1, *it2);
		}

		if (internal) {
			intersection[id1][id2] = val;
		}
		return val;
	}

	int QuartetDistance::binomCoef_n_2(int x)
	{
		return x * (x - 1) / 2;
	}

} // end of namespace
