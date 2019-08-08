//
// File: Subtree.h
// Created by: Anna Pawelczyk
// Created on: 19 May 2011, 12:15
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

#ifndef TRIPLETS_H
#define	TRIPLETS_H

#include <Bpp/Phyl/TreeTemplate.h>
#include <list>
#include <cstring>
using namespace bpp;
using namespace std;
namespace tools {
/**
 * @brief Helper class for Triplet distance - stores helpful data about trees.
 */    
class TreeParams
{
    friend class Triplets;
private:
    int inSize;
    int lSize;
    list<int>* inSons;    
    int* subTr;
    int rootId;
    
public:
    TreeParams(Node* root, int n, int l);
    ~TreeParams();
    int getResolved();
    int getUnresolved(int R);
    
private:
    int getResolved_in(int id);
    static int choose2(int a);
    static int choose3(int a);
    int countSubTrSize(Node* r);
    //collecting sons being internal nodes. Attantion! Internal nodes have ids >= l
    void setInSons(Node* root, int rootNewId);
};

    

/**
 * @brief Computing the Triplets distance.
 * \n Adapted from
 * M. Bansal, J. Dong, and D. Fernández-Baca. Comparing and aggregating partially 
 * resolved trees. LATIN 2008: Theoretical Informatics, 17:72–83, 2008.
 */
class Triplets
{
private:
    TreeParams *trP1;
    TreeParams *trP2;
    int ** intersection;
    int lSize;
    
public:
    Triplets(const TreeTemplate<Node>& tr1In, const TreeTemplate<Node>& tr2In);
    ~Triplets();
    int getDistance();
    
private:    
    int countIntersection(Node* root1, Node* root2, int** intersection);
    int getSameResolved();
    int getResolvedT1();
    int ro(int u2, int u, int v);
    int b_nega(int b, int a);
};


} // end of namespace
#endif	/* SUBTREE_H */
