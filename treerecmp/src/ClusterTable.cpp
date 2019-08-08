//
// File: ClusterTable.cpp
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

#include "ClusterTable.h"
namespace tools
{
ClusterTable::ClusterTable(int leavesSize, PostorderTree& tr)
{
    internalNodesNum = 0;
    for (int i = 0; i < leavesSize; i++) {
            clusterArray.push_back(new ClusterElement());
    }
    int postorderLeafPosition = 0;
    int top;
    int bottom;

    int sizeTr = tr.getNumberOfNodes();
    for (int i = 0; i < sizeTr; i++) {
        PostorderTree::NodeAgent* agent = tr[i];
        if (agent->subNodesSize == 0) {
            clusterArray.at(agent->nodeId)->postorderLeafPos = postorderLeafPosition;
            top = postorderLeafPosition;
            postorderLeafPosition++;
        } else {
            int loc;			// Location of cluster - the position in cllusterArray
            int leftLeaf = tr[i - agent->subNodesSize]->nodeId;
            bottom = clusterArray.at(leftLeaf)->postorderLeafPos;
            PostorderTree::NodeAgent *agentNext = tr[i + 1];
            loc = (agentNext->subNodesSize == 0) ? top : bottom;
            clusterArray.at(loc)->lowestElementPos = bottom;
            clusterArray.at(loc)->highestElementPos = top;
            clusterArray.at(loc)->branchW = agent->branchW;
            weight += clusterArray.at(loc)->branchW;
            internalNodesNum++;
        }
    }
}

ClusterTable::ClusterTable(const ClusterTable& orig) {
    clusterArray = orig.clusterArray;
}

ClusterTable::~ClusterTable()
{
    for (vector<ClusterElement*>::iterator it = clusterArray.begin(); it != clusterArray.end(); it++) {
        delete *it;
    }
}

void ClusterTable::removeUncommonElements(PostorderTree& otherTr)
{
    internalNodesNum = 0;
    stack<ClusterTable::Listing*> otherTrListingStack;
    for (PostorderTree::const_iterator naIt = otherTr.begin(); naIt != otherTr.end(); naIt++) {
        PostorderTree::NodeAgent* nAgent = *naIt;
        if (nAgent->subNodesSize == 0) {	// agent of the leaf node
            int pos = getPostorderPosForNode(nAgent->nodeId);
            Listing *leafListing = new Listing(pos, pos, 1, 1);
            otherTrListingStack.push(leafListing);
        } else {				// agent of the internal node
            Listing *clusterListing = new Listing();
            int subSize = nAgent->subNodesSize;
            do {
                Listing *childListing = otherTrListingStack.top();
                clusterListing->updateWithNewElement(childListing);
                subSize -= childListing->subNodesSize;
                otherTrListingStack.pop();
                delete childListing;
            } while (subSize != 0);
            clusterListing->branchW = nAgent->branchW;
            otherTrListingStack.push(clusterListing);

            /* Here the clusterListing element (element that describes a cluster 
             * from the otherTr tree) is checked whether it exists in
             * this clusterTable. If so, the clusterElement that corresponds to the
             * clusterListing is marked valid.
             *
             * The this clusterTable instance contains clusterElements that are described as so:
             * if clusterElement.lowestElementPos = a and clusterElement.highesttElementPos = b,
             * then cluster contains the leaves with a, a+1,..,b-1, b postorderLeafPos.
             *
             * The clusterListing item describing  cluster of the tree that is compared with this clusterTable
             * is here checked whether it contains the same leaves as one of those in this clusterTable.
             * To be so it must first of all contain n nodes that in this clusterTable have succesive
             * a, a+1,..,b-1, b postorderLeafPos values (retrived in getPostorderPosForNode()).
             * This is cheked by the clusterListing->isCompatible() method.
             * Secondly, the algorithm checks whether the compatible clusterListing exists
             * in the clusterTable so whether it is a clusterElement. This is performed by this.getPositionIfContains.
             * The clusterListing may occur in the this.clusterArray either on the position of
             * the lowest (a) or highets (b) postorderLeafPosition of a clusterElement.

             */
            weight += clusterListing->branchW;
            if (clusterListing->isCompatible()) {
                int pos = this->getPositionIfContains(clusterListing);
                if (pos != -1) {
                    clusterArray[pos]->valid = true;
                    internalNodesNum++;
                    double w = clusterListing->branchW + clusterArray[pos]->branchW
                        - abs(clusterListing->branchW - clusterArray[pos]->branchW);
                    weight -= w;
                }
            }
        }
    }
}

int ClusterTable::getPostorderPosForNode(int nodeId)
{
    return clusterArray.at(nodeId)->postorderLeafPos;
}

int ClusterTable::getPositionIfContains(ClusterTable::Listing *listing)
{
    int l = listing->lowestLeafPos;
    int h = listing->highestLeafPos;
    if (clusterArray[l]->isElement(l, h)) return l;
    if (clusterArray[h]->isElement(l, h)) return h;
    return -1;
}


} // end of namespace

