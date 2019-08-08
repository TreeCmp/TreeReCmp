//
// File: PostorderTree.h
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

#ifndef POSTORDERTREE_H
#define	POSTORDERTREE_H
#include <vector>
using namespace std;
#include <Bpp/Phyl/Node.h>
#include <Bpp/Phyl/TreeTemplate.h>
using namespace bpp;

namespace tools {
    /**
     * @brief The way of storing phylogenetic tree's topology, convenient e.i. for 
     * O(n) Robinson-Foulds algorithm.
     * This structure allows to scan the tree in a way that all the leaves of each
     * cluster are examined first, and then the internal node that is the root
     * of the cluster is examined. Its subNodesSize value allows to determine how many
     * of preceding elements in a PSW array are the cluster’s elements (leaves
     * and internal nodes).
     */
class PostorderTree {
public:
    /**
     * @brief Structure of a node of PostorderTree that stores information 
     * about the real id of a node and the number of its descendants.
     */
    class NodeAgent {
    public:
        int nodeId;
        int subNodesSize;       // number of node's descendants
        double branchW;         // branch weight
        NodeAgent(int nodeId, int leavesSize, double branchWIn)
        {
            this->nodeId = nodeId;
            this->subNodesSize = leavesSize;
            this->branchW = branchWIn;
        }
    };
    typedef vector<NodeAgent*>::iterator iterator;
    typedef vector<NodeAgent*>::const_iterator const_iterator;
private:
    vector<NodeAgent*> postorderNodes;
    bool isWeighted;

public:
    /**
     * @brief During a DFS search of a tree that root is given as parameter,
     * a node x on which searching is finished (there is no more child to
     * examine) obtains a successive position in a PSW array and this position
     * is its pswId. Under this position, the array stores data about the
     * node in NodeAgent class: node’s original label and the number of node’s descendants.
     * @param[in]   root    root to the bpp:Tree that will be the base for the new PostorderTree.
     * @param[in]   isWeightedIn    informs whether input tree is weighted.
     */
    PostorderTree(const Node* root, bool isWeightedIn = false);
    PostorderTree(const PostorderTree& orig);
    virtual ~PostorderTree();
    static void rootTree(TreeTemplate<Node>& tr);
    NodeAgent* operator[](int pos);
    iterator begin() { return postorderNodes.begin(); }
    const_iterator begin() const { return postorderNodes.begin(); }
    iterator end() { return postorderNodes.end() - 1; }
    const_iterator end() const { return postorderNodes.end() - 1; }
    NodeAgent* getNodeAgent(int pos);
    int getNumberOfNodes();
private:
    int setPostorderList(Node& root);

};
} // end of namespace

#endif	/* POSTORDERTREE_H */

