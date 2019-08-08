//
// File: Partitioning.cpp
// Created by: Anna Pawelczyk
// Created on: 23 May 2011, 10:36
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

#include <limits.h>
#include "Partitioning.h"

namespace tools
{
LeavesSets::LeavesList* PairLeavesSets::browseTree(int rootId, const TreeTemplate<Node>& tr, int& rootIndex, vector<vector<int> >& splitNode)
{
    vector<int> sons = tr.getSonsId(rootId);
    if (sons.size() == 0) {			// a leaf
        return new LeavesList(rootId);		
    } else {				//2 sons in bifurcating tree
        LeavesList* l1 = browseTree(sons.at(0), tr, rootIndex, splitNode);
        LeavesList* l2 = browseTree(sons.at(1), tr, rootIndex, splitNode);
        for (l1->resetIteration(); !l1->iterationFinished();) {
            int leafIdTr1 = l1->nextId();
            for (l2->resetIteration(); !l2->iterationFinished();) {
                int leafIdTr2 = l2->nextId();
                splitNode[leafIdTr1][leafIdTr2] = rootIndex;
                splitNode[leafIdTr2][leafIdTr1] = rootIndex;
            }
        }
        rootIndex++;
        l1->concat(l2);
        return l1;
    }
}


PairLeavesSets::PairLeavesSets(const TreeTemplate<Node> &tr1, const TreeTemplate<Node> &tr2)
{       
    leavesNum = tr1.getNumberOfLeaves();
    vector<int> tmp;
    tmp.assign(leavesNum , 0);
    splitNode1.assign(leavesNum, tmp);
    splitNode2.assign(leavesNum, tmp);

    int rootIndex = 0;
    browseTree(tr1.getRootId(), tr1, rootIndex, splitNode1);
    rootIndex = 0;
    browseTree(tr2.getRootId(), tr2, rootIndex, splitNode2);
    inNodesNum = rootIndex;   // both tree, as thet are bifurcating and have the same leaves sets - have the equal number of internal nodes           
}
int PairLeavesSets::getSize() { return inNodesNum; }

int PairLeavesSets::getCostMatrix(int** costMatrix)
{
// internNodes1Sizes[i] stores the number of pairs that node i is the ascentor of 
// acording t the splitNode1 sense
    int* internNodes1Sizes = new int[inNodesNum];        
    int* internNodes2Sizes = new int[inNodesNum];
    for (int a = 0; a < inNodesNum; a++) {
        for (int b = 0; b < inNodesNum; b++) {
            costMatrix[a][b] = 0;
        }
        internNodes1Sizes[a] = 0;
        internNodes2Sizes[a] = 0;
    }

    for (int i = 0; i < leavesNum; i++) {
        for (int j = i+1; j < leavesNum; j++)  {
            int internNode1 = splitNode1[i][j];
            int internNode2 = splitNode2[i][j];
            costMatrix[internNode1][internNode2]++;
            internNodes1Sizes[internNode1]++;
            internNodes2Sizes[internNode2]++;
        }
    }
    
    // count costs and reverse values
    for (int i = 0; i < inNodesNum; i++) {
        for (int j = 0; j < inNodesNum; j++)  {
            costMatrix[i][j] = internNodes1Sizes[i] + internNodes2Sizes[j] - 2*costMatrix[i][j];
        }
    }
    delete[] internNodes1Sizes;
    delete[] internNodes2Sizes;
	return 0;
}



/*******************PARTITIONING************************
 *******************************************************/
    
const int Partitioning::BITS_IN_INT = 32;
int taxonNum_or_MaxInt;

int Partitioning::getSize() {return largerBitList.size();}

void Partitioning::getCostReversedMatrix(vector<vector<int> >& costMatrix)
{
    int taxonNum_or_MaxInt = taxonsNumber;
    
    int smallerSize = smallerBitList.size();
    int largerSize = largerBitList.size();
    
    costMatrix.resize(largerSize);
        
    vector<int*>::iterator lIt;
    vector<int*>::iterator sIt;

    lIt = largerBitList.begin();
    for (int a = 0; a < largerSize; a++) {        
            costMatrix[a].resize(largerSize);
        vector<int*>::iterator sIt = smallerBitList.begin();
        for (int b = 0; b < smallerSize; b++) {
            int a1b1 = 0;
            int iterationsNum = taxonsNumber / (BITS_IN_INT*8) + 1;
            for (int i = 0; i < iterationsNum; i++) {
                a1b1 += getNumberOfSetBits( (*lIt)[i] ^ (*sIt)[i] );
            }
            a1b1 = taxonsNumber - a1b1;
            costMatrix[a][b] = min(a1b1, taxonNum_or_MaxInt - a1b1);
            sIt++;
        }
        lIt++;
    }   

    if (largerSize != smallerSize) {	// Usig dGMS1 for dummy v
        lIt = largerBitList.begin();
        for (int a = 0; a < largerSize; a++) {
            for (int b = smallerSize; b < largerSize; b++) {
                costMatrix[a][b] = (this->*dummyFun)(largerBitList.at(a));
                sIt++;
            }
            lIt++;
        }
    }
}


int Partitioning::getCostMatrix(int** costMatrix)
{    
    int smallerSize = smallerBitList.size();
    int largerSize = largerBitList.size();

    vector<int*>::iterator lIt;
    vector<int*>::iterator sIt;

    lIt = largerBitList.begin();
    for (int a = 0; a < largerSize; a++) {
        vector<int*>::iterator sIt = smallerBitList.begin();
        for (int b = 0; b < smallerSize; b++) {
            int a1b1 = 0;
            int iterationsNum = taxonsNumber / (BITS_IN_INT*8) + 1;
            for (int i = 0; i < iterationsNum; i++) {
                a1b1 += getNumberOfSetBits( (*lIt)[i] ^ (*sIt)[i] );
            }                    
            costMatrix[a][b] = min(a1b1, taxonNum_or_MaxInt - a1b1);
            sIt++;
        }
        lIt++;
    }   

    if (largerSize != smallerSize) {	// Usig dGMS1 for dummy v
        lIt = largerBitList.begin();
        for (int a = 0; a < largerSize; a++) {
            for (int b = smallerSize; b < largerSize; b++) {
                costMatrix[a][b] = (this->*dummyFun)(largerBitList.at(a));
                sIt++;
            }
            lIt++;
        }
    }
    return largerSize;
}

int Partitioning::dummyFunction1(int* unused)
{
    return ceil( floor((double)taxonsNumber / 2) / 2);
}
int Partitioning::dummyFunction2(int* bitBipart)
{
    int count = 0;
    int iterationsNum = taxonsNumber / BITS_IN_INT + 1;
    for (int i = 0; i < iterationsNum; i++) {
        count += getNumberOfSetBits(bitBipart[i]);
    }
    return min (count, taxonsNumber - count);
}

Partitioning::~Partitioning()
{
    delete pl1;
    delete pl2;
}

void Partitioning::setInitFields(const TreeTemplate<Node>& tr1, dummyFunType d)
{        
    if (pl1->size() > pl2->size()) {
        smallerBitList = pl2->getBitList();
        largerBitList = pl1->getBitList();
    } else {
        smallerBitList = pl1->getBitList();
        largerBitList = pl2->getBitList();
    }
    taxonsNumber = tr1.getNumberOfLeaves();
    switch (d) {
        case GMS1 : {dummyFun = &Partitioning::dummyFunction1; break;}
        case GMS2 : {dummyFun = &Partitioning::dummyFunction2; break;}
    }
}

Splitting::Splitting(const TreeTemplate<Node>& tr1, const TreeTemplate<Node>& tr2, dummyFunType d) 
{        
    pl1 = new BipartitionList(tr1); //BipartitionList bpl1(tr1, true);
    pl2 = new BipartitionList(tr2);
    setInitFields(tr1, d);
    taxonNum_or_MaxInt = taxonsNumber;
}


Clustering::Clustering(const TreeTemplate<Node>& tr1, const TreeTemplate<Node>& tr2, dummyFunType d)
{	
    pl1 = new ClusterList(tr1);
    pl2 = new ClusterList(tr2);
    setInitFields(tr1, d);
    taxonNum_or_MaxInt = INT_MAX;
}


//****************Deprecated*******************8
int Clustering::countDistance(int *bitClust1, int *bitClust2)
{
    int a1b1 = 0;
    int iterationsNum = taxonsNumber / BITS_IN_INT + 1;
    for (int i = 0; i < iterationsNum; i++) {
        a1b1 += getNumberOfSetBits( bitClust1[i] ^ bitClust2[i] );
    }
    return a1b1;
}

int Splitting::countDistance(int *bitBipart1, int *bitBipart2)
{
    int a1b1 = 0;
    int iterationsNum = taxonsNumber / BITS_IN_INT + 1;
    for (int i = 0; i < iterationsNum; i++) {
        a1b1 += getNumberOfSetBits( bitBipart1[i] ^ bitBipart2[i] );
    }
    return min(a1b1, taxonsNumber - a1b1);
}

int Partitioning::fixReversedWeightsSum(int sum)
{
    return taxonsNumber * largerBitList.size() - sum;
}  

void PairLeavesSets::getCostReversedMatrix(vector<vector<int> >& costMatrix)
{
    // internNodes1Sizes[i] stores the number of pairs that node i is the ascentor of 
    // acording t the splitNode1 sense
    int* internNodes1Sizes = new int[inNodesNum];        
    int* internNodes2Sizes = new int[inNodesNum];
    memset(internNodes1Sizes, 0, inNodesNum * sizeof(int));
    memset(internNodes2Sizes, 0, inNodesNum * sizeof(int));

    vector<int> emptyVec;
    emptyVec.resize(inNodesNum, 0);
    costMatrix.resize(inNodesNum, emptyVec);

    for (int i = 0; i < leavesNum; i++) {
        for (int j = i+1; j < leavesNum; j++)  {
            int internNode1 = splitNode1[i][j];
            int internNode2 = splitNode2[i][j];
            costMatrix[internNode1][internNode2]++;
            internNodes1Sizes[internNode1]++;
            internNodes2Sizes[internNode2]++;
        }
    }

    // count costs and reverse values
    for (int i = 0; i < inNodesNum; i++) {
        for (int j = 0; j < inNodesNum; j++)  {
            costMatrix[i][j] = internNodes1Sizes[i] + internNodes2Sizes[j] - 2*costMatrix[i][j];
            costMatrix[i][j] = leavesNum * leavesNum - 2*costMatrix[i][j];
        }
    }
}
int PairLeavesSets::fixReversedWeightsSum(int sum)
{
    return leavesNum * leavesNum * inNodesNum - sum ;// * (this->size * this->size) - dist;  
}

} // end of namespace

