//
// File: TreesManip.cpp
// Created by: Anna Pawelczyk
// Created on: 19 May 2011, 20:42
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

#include "TreesManip.h"
namespace tools {    

TreeTemplate<Node>* TreesManip::createOrderedTrees(const TreeTemplate<Node>& trIn)
{
    TreeTemplate<Node>* trOut = trIn.clone();
    //Internally, map containers keep their elements ordered by their keys from lower to higher
    map<string, Node*> sortedLeaves;
    vector<Node*> nodesList = trOut->getNodes();
    int internalNodeId = nodesList.size() - 1;
    for(vector<Node*>::iterator nIt = nodesList.begin(); nIt != nodesList.end(); nIt++) {
        Node *n = *nIt;
        if (n->isLeaf()) {
            sortedLeaves[n->getName()] = n;
        } else {
            n->setId(internalNodeId);
            internalNodeId--;
        }
    }
    int leafId = 0;
    for(map<string, Node*>::iterator pairIt = sortedLeaves.begin(); pairIt != sortedLeaves.end(); pairIt++) {
        (*pairIt).second->setId(leafId);
        leafId++;
    }
    if (leafId -1 != internalNodeId) {
        throw Exception("Unknown error - No1");
    }
    return trOut;	
}


void TreesManip::getLeavesLevels(Node* root, int level, vector<int>& leavesIdLevel)
{
    vector<Node*> sons = root->getSons();
    if (sons.size() == 0) {			// a leaf
            leavesIdLevel[root->getId()] = level;		
    } else {				//2 sons in bifurcating tree
        for (int i = 0; i < sons.size(); i++) {
            getLeavesLevels(sons.at(i), level + 1, leavesIdLevel);
        }
    }
}

} // end of namespace

