//
// File: ClusterTable.h
// Created by: Anna Pawelczyk
// Created on: 20 May 2011, 12:21
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

#ifndef CLUSTERTABLE_H
#define	CLUSTERTABLE_H
#include <limits.h>
#include <stack>
#include <list>
#include <vector>
using namespace std;
#include "PostorderTree.h"

namespace tools {
/**
 * @brief structure that stores information on the clusters, used for 
 * the O(n) Robinson-Foulds algorithm. 
 * For each cluster in the examined tree the algorithm checks 
 * whether the cluster is already in the ClustersStructure. Finally 
 * the cluster table stores just the clusters of the resultant consensus tree.
 */
class ClusterTable {
public:
    /**
    * @brief Stores information about an element of ClusterTable.
    */
    class ClusterElement {
    public:
        int lowestElementPos;
        int highestElementPos;
        int postorderLeafPos;
        double branchW;
        bool valid;
        ClusterElement()
        {
                valid = false;
                lowestElementPos = INT_MIN;
                highestElementPos = INT_MIN;
                postorderLeafPos = INT_MIN;
        }
        ClusterElement(const ClusterElement& oryginal)
        {
                lowestElementPos = oryginal.lowestElementPos;
                highestElementPos = oryginal.highestElementPos;
                postorderLeafPos = oryginal.postorderLeafPos;
                valid = oryginal.valid;
        }
        bool isElement(int lowest, int highest)
        {
                return (lowestElementPos == lowest && highestElementPos == highest);
        }
    };
    class Listing
    {
    public:
        int lowestLeafPos;
        int highestLeafPos;
        int n;
        double branchW;
        int subNodesSize;
        Listing(int lPosIn, int hPosIn, int nIn, int subNodesIn)
        {
                lowestLeafPos = lPosIn;
                highestLeafPos = hPosIn;
                n = nIn;
                subNodesSize = subNodesIn;
        }
        Listing()
        {
                lowestLeafPos = INT_MAX;
                highestLeafPos = INT_MIN;
                n = 0;
                subNodesSize = 1;
        }
        void updateWithNewElement(Listing *other)
        {
                lowestLeafPos = (other->lowestLeafPos < lowestLeafPos) ? other->lowestLeafPos : lowestLeafPos ;
                highestLeafPos = (other->highestLeafPos > highestLeafPos) ? other->highestLeafPos : highestLeafPos;
                n += other->n;
                subNodesSize += other->subNodesSize;
        }
        bool isCompatible()
        {
                return (n == highestLeafPos - lowestLeafPos + 1);
        }
    };
private:
    vector<ClusterElement *> clusterArray;
    int internalNodesNum;
    double weight;

public:
    ClusterTable(int leavesSize, PostorderTree& tr);
    ClusterTable(const ClusterTable& orig);
    virtual ~ClusterTable();
    ClusterElement * operator[](int pos) { return clusterArray.at(pos); }
    void removeUncommonElements(PostorderTree& otherTr);
    int getNumberOfInternalNodes() { return internalNodesNum; }
    int getWeight() { return weight; }
private:
    int getPostorderPosForNode(int nodeId);
    void markValid(ClusterTable::Listing *listing);
    int getPositionIfContains(ClusterTable::Listing *listing);

};
} // end of namespace

#endif	/* CLUSTERTABLE_H */

