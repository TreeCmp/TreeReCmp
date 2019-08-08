//
// File: TreesManip.h
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

#ifndef TREESMANIP_H
#define	TREESMANIP_H

#include <vector>
using namespace std;
#include <Bpp/Phyl/TreeTemplate.h>
using namespace bpp;

namespace tools {
/**
 * @brief General methods used across the PhylotreeDist package
 * For manipulating the information about trees.
 */
class TreesManip {
public:
    /**
     * @author Ania Pawelczyk
     * @brief The method returns a new tree with the topology of input tree, but with
     * nodes ids ordered in the following way (lNum = number_of_leaves) :
     * \n  - the ids [ 0 .. lNum-1 ] are given to leaves in the alphabetical order of their names
     * \n  - the ids [ lNum .. tree.numberOfNodes-1 ] are given arbitrary to internal nodes
     * 
     * \n\nThe ordering method assures that if two trees with the same leaves sets are ordered,
     * their leaves have the same ids.
     * \n\n time complexity: O(n logn)
     * @param trIn - the tree which topology will have the new ordered tree.
     * @returns a pointer to new tree with the same topology as the input has and with ordered leaves.
     */
    static TreeTemplate<Node>* createOrderedTrees(const TreeTemplate<Node>& trIn);
    static void getLeavesLevels(Node* root, int level, vector<int>& leavesIdLevel);
private:

};
} // end of namespace

#endif	/* TREESMANIP_H */

