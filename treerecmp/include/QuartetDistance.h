//
// File: QPartitionList.h
// Created by: Anna Pawelczyk
// Created on: 15 Aug 2011, 11:10
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

#ifndef QPARTETS_H
#define	QPARTETS_H
#include <Bpp/Phyl/TreeTemplate.h>
#include <list>
using namespace bpp;
using namespace std;


namespace tools {
/**
 * @brief Helper class for Quartet distance - stores helpful data about trees.
 */
class TreeParams2
{
    friend class QuartetDistance;
private:
    int inSize;
    int lSize;
    list<int>* inSons;    
    int* subTr;
    int rootId;
    
public:
    TreeParams2(Node* root, int n, int l);
    ~TreeParams2();
    
private:    
    static int choose2(int a);
    static int choose3(int a);    
    int countSubTrSize(Node* r);
    //collecting sons being internal nodes. Attantion! Internal nodes have ids >= l
    void setInSons(Node* root, int rootNewId);
};


/**
 * @brief Computing the Quartet distance. 
 * The method bases on finding in the two input trees the butterfly topology
 * quartets that are and are not shared by the trees. A butterfly topology quartet
 * is a topology of four leaves where one pair of them is separated from the other
 * pair by an edge.
 * \n Adapted from
 * C. Christiansen, T. Mailund, C. Pedersen, and M. Randers. Computing
 * the quartet distance between trees of arbitrary degree. Algorithms in
 * Bioinformatics, 3692:77–88, 2005.
 */
class QuartetDistance
{
private:    
    TreeParams2 *trP1;
    TreeParams2 *trP2;
    int ** intersection;
    int lSize;
    Node* r1;
    Node* r2;

public:
    QuartetDistance(const TreeTemplate<Node>& tr1In, const TreeTemplate<Node>& tr2In);
    ~QuartetDistance();    
    int getDistance();    
    int getNonshared();
    int getShared();

private:
    int getSingleTrQuartetsSize(int size, int* subTrsize, list<int>* inSons);
    int binomCoef_bin(int a);
    int non_product(int i, int j);
    int a_b(int i, int j);
    int aneg_bneg(int i, int j);
    int a_bneg(int i, int j);
    int aneg_b(int i, int j);
    int in_a_b(int i, int j);
    int in_aneg_bneg(int i, int j);
    int in_a_bneg(int i, int j);
    int in_aneg_b(int i, int j);    
    int countIntersection(Node* root1, Node* root2);    
    int binomCoef_n_2(int x);
};

} // end of namespace
#endif	/* QPARTITIONLIST_H */

