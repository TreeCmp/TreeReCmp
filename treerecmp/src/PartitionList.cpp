//
// File: ClusterList.cpp
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

#include "PartitionList.h"
#include <string.h>     // strcat

namespace tools
{
const int PartitionList::BITS_IN_INT = 32;

PartitionList::PartitionList() {}        

PartitionList::~PartitionList() 
{
    for(vector<int*>::iterator it = bitList.begin(); it != bitList.end(); it++) {
        delete[] *it;
    }
}

PartitionList::PartitionList(int count)
{
    intCount = ceil((float)count / BITS_IN_INT );  
}


const vector<int *>& PartitionList::getBitList()
{
    return bitList;
}

int PartitionList::size() const
{
    return bitList.size();
}


int* PartitionList::browseTree(Node* root)
{
    browseTree_do(root);

    //clear undeleted clusters
    for(vector<int*>::iterator it = leavesBitList.begin(); it !=leavesBitList.end(); it++) {
        delete[] *it;
    }
    leavesBitList.clear();
	return 0;
}

int* PartitionList::browseTree_do(Node* root)
{
    int* cluster = new int[intCount];
    for (int i = 0; i < intCount; i++) cluster[i] = 0;
    vector<Node*> sons = root->getSons();
    if (sons.size() == 0) {		// a leaf
        setLeafBit(root->getId(), cluster);
        leavesBitList.push_back(cluster);
    } else {
        for (vector<Node*>::iterator it = sons.begin(); it != sons.end(); it++) {
            joinBits(cluster, browseTree_do(*it));
        }
        bitList.push_back(cluster);
    }
    return cluster;
}

/*
* This bases on the fact that the trees have the same leaves ids
*/
void PartitionList::setLeafBit(int leafId, int* clusterBitList)
{
    int listIntPosition = floor((float)leafId / BITS_IN_INT);
    int listPosition = leafId % BITS_IN_INT;
    clusterBitList[listIntPosition] = 1 << listPosition;
}
void PartitionList::joinBits(int *cl1, const int *cl2)
{
    for (int i = 0; i < intCount; i++) {
        cl1[i] = cl1[i] | cl2[i];
    }
}


ClusterList::ClusterList(const TreeTemplate<Node>& tr) : PartitionList(tr.getNumberOfLeaves())
{
    //intCount = ceil(43 / BITS_IN_INT );
    Node *root = ((const_cast<TreeTemplate<Node>& >(tr)).getRootNode());

    vector<Node*> sons = root->getSons();
    for (vector<Node*>::iterator sIt = sons.begin(); sIt != sons.end(); sIt++) {
        browseTree(*sIt);       
    }
}

BipartitionList::BipartitionList(const TreeTemplate<Node>& tr) : PartitionList(tr.getNumberOfLeaves())
{
    //intCount = ceil(43 / BITS_IN_INT );
    Node *root = ((const_cast<TreeTemplate<Node>& >(tr)).getRootNode());

    vector<Node*> sons = root->getSons();
    for (vector<Node*>::iterator sIt = sons.begin(); sIt != sons.end(); sIt++) {
        browseTree(*sIt);       
    }

    //remove one of repeated fake-root bipartitions
    if (root->getNumberOfSons() == 2) {
        delete[] bitList.back();
        bitList.pop_back();
    }

}

} // end of namespace
