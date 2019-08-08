//
// File: Partitioning.h
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

#ifndef PARTITIONING_H
#define	PARTITIONING_H

#include <Bpp/Phyl/TreeTemplate.h>
#include <cmath>
#include <cstring>
#include <Bpp/Phyl/BipartitionList.h>
#include "PartitionList.h"
using namespace std;
using namespace bpp;

namespace tools {
/**
 * @brief Interface for operationg on two description elements.\n
 * A description element is a piece of information about topology of a tree,
 * for example, it can be a branch, split or cluster in some phylogenetic
 * tree description.
 */    
class ITwoTreesDescriptionElements
{    
public:
    /**
     * @brief Counts a cost matrix and reverses its values:
     * - cost is the distance between a pair of splits where the pair incldes
     * one split from T1 and one from T2.
     * - each finally returned in the cost matrix value equals (maxCost - value).
     * This is done to find minimum weight perfect matching by a Hungarian (Maximum
     * Weight Perfect Matching) algorithm.
     *
     * O(n^3) == O(n^2*countDistance complexity)
     * @param costMatrix - a matrix that would store the result.
     */
    virtual int getCostMatrix(int** costMatrix) = 0;
    virtual int getSize() = 0;

    virtual int fixReversedWeightsSum(int sum) = 0;
    virtual void getCostReversedMatrix(vector<vector<int> >& costMatrix) = 0;
};
/****************************k-leaves Subsets*************************/
/**
 * @brief Abstract class for operating on two description elements, where a description
 * element is a subset of leaves with specified size.\n
 */  
class LeavesSets : public ITwoTreesDescriptionElements
{
protected:   
    class LeavesList
    {
    private:
        struct Element
        {
        public:
            Element(int valIn)
            {
                val = valIn;
                next = NULL;
            }
            int val;
            Element* next;
        };

        Element* first;
        Element* last;
        int listSize;
        Element* iterator;
    public:
        LeavesList(int val)
        {
            first = new Element(val);
            last = first;
            listSize = 1;
        }
        ~LeavesList()
        {
            Element* el = first;
            while (el != NULL) {
                Element* elNext = el->next;
                delete el;
                el = elNext;
            }
            delete el;
        }
        int size()
        {
            return listSize;
        }
        void concat(LeavesList* other)
        {
            last->next = other->first;
            last = other->last;
            listSize += other->size();
        }
        void resetIteration()
        {
            iterator = first;
        }
        int nextId()
        {
            if (iterator == NULL) throw Exception("LeavesList's iterator out of range");
            int val = iterator->val;
            iterator = iterator->next;
            return val;
        }
        bool iterationFinished()
        {
            return iterator == NULL;
        }
    };
    // Matrices of size: number of leaves. They keep the internal node's generational indice.
    // The splitNode1[i][j] informs what is the generational index and level in the tree1
    // at which two leaves of id i and j split. 
    // A generational index is a value 0..inNodesNum-1 assigned to n internal node during browseTree.
    // It has no other usage than to be the smallest possible distinctive internal node description.
    vector<vector<int> > splitNode1;        
    vector<vector<int> > splitNode2;        
    int inNodesNum;                 // the number internal nodes
    int leavesNum;
public:
    virtual int getCostMatrix(int** costMatrix) = 0; 

    virtual void getCostReversedMatrix(vector<vector<int> >& costMatrix) = 0;
    virtual int fixReversedWeightsSum(int sum) = 0;
private:
    virtual LeavesList* browseTree(int rootId, const TreeTemplate<Node>& tr, int& rootIndex, vector<vector<int> >& splitNode) = 0;
};

/**
 * @brief Two description elements operations, where a description element is 
 * a two-elements subset of leaves with specified size.\n
 */
class PairLeavesSets : public LeavesSets
{
public:
    PairLeavesSets(const TreeTemplate<Node> &tr1, const TreeTemplate<Node> &tr2);        
    int getSize();        
    int getCostMatrix(int** costMatrix); 
    void getCostReversedMatrix(vector<vector<int> >& costMatrix);        
    int fixReversedWeightsSum(int sum);        
private:
    LeavesList* browseTree(int rootId, const TreeTemplate<Node>& tr, int& rootIndex, vector<vector<int> >& splitNode);
};


/****************************Partitioning*************************/
/**
 * @brief Abstract class for operating on two description elements, where a description
 * element is a subset of leaves that share a particular property.\n
 */
class Partitioning : public ITwoTreesDescriptionElements
{
public:
    enum dummyFunType {
        GMS1,
        GMS2
    };
private:
    typedef int (Partitioning::*methodPtr)(int* bitBipart);
    methodPtr dummyFun;
protected:
    PartitionList *pl1;
    PartitionList *pl2;
    static const int BITS_IN_INT;
    vector<int *> largerBitList;
    vector<int *> smallerBitList;
    int taxonNum_or_MaxInt;
    int taxonsNumber;
    /**
     * @brief Each bipartition is a sequence of bits corresponding to the trees
     * the bipartitions come from. That trees have the same leaves set.
     * Bit position in both bipartitions corresponds to the same leaf.
     * In this way XORing both bipartitions data - we obtain the value which
     * number of bits is the number of leaves.
     * Lets A1|A2 denote a bipartition in treeA and B1|B2 denote a bipartition in treeB
     * According to Bogdanowicz "Comparing Phylometric Trees(...)"
     * the dist is 0.5 min{ |A1xorB1|+|A2xorB2|, |A1xorB2|+|A2xorB1| } which equals to
     *  min {|A1|+|B1|-2|A1&B1|, leavesSize - |A1|+|B1|-2|A1&B1|}
     *
     * @param bipA - the array of bits of (leavesSize / sizeof(int) + 1) size. If a bit
     * is set it means that a leaf that corresponds to the bit belongs to the bipartition.
     * @param bipB - the array of bits of (leavesSize / sizeof(int) + 1) size. If a bit
     * is set it means that a leaf that corresponds to the bit belongs to the bipartition.
     * @param leavesSize - the number of leaves of the tree bipartition bipA comes from which
     * is equal to the number leaves of the tree bipartition bipb comes from.
     *
     * @return The distance between bipA and bipB
     */
    virtual int countDistance(int *bitDescrEl1, int *bitDescrEl2) = 0;
    int dummyFunction1(int* unused);
    int dummyFunction2(int* bitBipart);
    void setInitFields(const TreeTemplate<Node>& tr1, dummyFunType d);
    /**
     * @brief Counting bits of a 32bit value. The variable-precision SWAR algorithm.
     * For delatis see http://graphics.stanford.edu/~seander/bithacks.html#CountBitsSetParallel
     * @param val - the 32bit value of which the number of set bits is count.
     * @return The number of set bits
     */
    static int getNumberOfSetBits(int val)
    {
        val = val - ((val >> 1) & 0x55555555);
        val = (val & 0x33333333) + ((val >> 2) & 0x33333333);
        return ((val + (val >> 4) & 0xF0F0F0F) * 0x1010101) >> 24;
    }
public:
    int fixReversedWeightsSum(int sum);
    void getCostReversedMatrix(vector<vector<int> >& costMatrix);        
    int getCostMatrix(int** costMatrix);
    int getSize();
    ~Partitioning();
};

/**
 * @brief This class defines two description elements operations. Here the description
 * element is a subset of leaves called split (or bipartition): if we remove an edge 
 * from a tree then we divide it up into two components. A split is a set of 
 * leaves of one of those components.
 */
class Splitting : public Partitioning
{
public:
    Splitting(const TreeTemplate<Node>& tr1, const TreeTemplate<Node>& tr2, dummyFunType d = GMS2);

private:
    int countDistance(int *bipA, int *bipB);
};

/**
 * @brief This class defines two description elements operations. Here the description
 * element is a subset of leaves called cluster: Given a rooted tree T, if we 
 * choose a one internal node v, then the set of leaves in T that are 
 * descendants of v (including v if it is a leaf) is called a cluster.
 */
class Clustering : public Partitioning
{
public:
    Clustering(const TreeTemplate<Node>& tr1, const TreeTemplate<Node>& tr2, dummyFunType d = GMS1);
private:
    int countDistance(int *bipA, int *bipB);	
};


} // end of namespace
#endif	/* PARTITIONING_H */

