//
// File: ClusterTable.h
// Created by: Anna Pawelczyk
// Created on: 23 May 2011, 10:41
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



#ifndef PARTITIONLIST_H
#define	PARTITIONLIST_H
#include <vector>
#include <Bpp/Phyl/TreeTemplate.h>
#include <cmath>
using namespace std;
using namespace bpp;

namespace tools
{
/**
 * @brief Interface for description elements container.\n
 * A description element is a piece of information about topology of a tree,
 * for example, it can be a branch, split or cluster in some phylogenetic
 * tree description. A partition is a subset of leaves that share a particular property.
 */ 
class PartitionList
{
protected:
    vector<int *> bitList;
    int intCount;
    static const int BITS_IN_INT;
    int* browseTree(Node* root);
public:
    PartitionList();
    ~PartitionList();
    PartitionList(int count);
    /**
     *  O(n^2 logn)
     **/        
    const vector<int *>& getBitList();
    int size() const;
private:
    /*
     * This bases on the fact that the trees have the same leaves ids
     */
    void setLeafBit(int leafId, int* clusterBitList);
    void joinBits(int *cl1, const int *cl2);

    int* browseTree_do(Node* root);
    vector<int *> leavesBitList;
}; 


/**
 * @brief Description elements container. Here a description
 * element is a subset of leaves called cluster: Given a rooted tree T, if we 
 * choose a one internal node v, then the set of leaves in T that are 
 * descendants of v (including v if it is a leaf) is called a cluster.
 */
class ClusterList : public PartitionList
{     
public:   
    ClusterList(const TreeTemplate<Node>& tr);
};

/**
 * @brief Description elements container. Here a description
 * element is a subset of leaves called split (or bipartition): if we remove an edge 
 * from a tree then we divide it up into two components. A split is a set of 
 * leaves of one of those components.
 */
class BipartitionList : public PartitionList
{
public:
    BipartitionList(const TreeTemplate<Node>& tr);
};

}	// end of namespace

#endif	/* PARTITIONLIST_H */

