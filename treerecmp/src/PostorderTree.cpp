//
// File: PostorderTree.cpp
// Created by: Anna Pawelczyk
// Created on: 20 May 2011, 11:51
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

#include "PostorderTree.h"
namespace tools {

PostorderTree::PostorderTree(const Node* rootIn, bool isWeightedIn)
{
    isWeighted = isWeightedIn;
    Node root(*rootIn);
    setPostorderList(root);
    // Fake last node for the purposes of Day's R-F algorithm
    NodeAgent* fakeAgent = new NodeAgent(-1, 0, -1000);
    postorderNodes.push_back(fakeAgent);
}

PostorderTree::PostorderTree(const PostorderTree& orig)
{
    postorderNodes = orig.postorderNodes;
}

PostorderTree::~PostorderTree()
{
    for (vector<NodeAgent*>::iterator it = postorderNodes.begin(); it != postorderNodes.end(); it++) {
        delete *it;
    }
}

void PostorderTree::rootTree(TreeTemplate<Node>& tr)
{
    int nullId = tr.getNumberOfLeaves() - 1;
    Node* nullNode = tr.getNode(nullId);
    Node* newRootNode = nullNode->getFather();
    newRootNode->removeSon(nullNode);
    delete nullNode;
    tr.rootAt(newRootNode);
    newRootNode->setDistanceToFather(1);
}

PostorderTree::NodeAgent* PostorderTree::operator[](int pos)
{
    return getNodeAgent(pos);
}

PostorderTree::NodeAgent* PostorderTree::getNodeAgent(int pos)
{
    return postorderNodes.at(pos);
}
int PostorderTree::getNumberOfNodes()
{
    // The number of elements without fakeAgent
    return postorderNodes.size() - 1 ;
}

int PostorderTree::setPostorderList(Node& root)
{
    if (root.isLeaf()) {
        double w = isWeighted ? root.getDistanceToFather() : 1;
        NodeAgent* agent = new NodeAgent(root.getId(), 0, w);
        postorderNodes.push_back(agent);
        return 1;
    }

    vector<Node*> sons = root.getSons() ;
    int subtreeSize = 0;
    for (vector<Node*>::iterator it = sons.begin(); it != sons.end(); it++) {
        Node* n = *it;
        subtreeSize += setPostorderList(*n);
    }
    double w = isWeighted ? root.getDistanceToFather() : 1;
    NodeAgent* agent = new NodeAgent(root.getId(), subtreeSize, w);
    postorderNodes.push_back(agent);
    return subtreeSize + 1;
}


} // end of namespace
