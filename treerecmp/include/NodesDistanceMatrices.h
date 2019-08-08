//
// File: NodesDistanceMatrices.h
// Created by: Anna Pawelczyk
// Created on: 25 May 2011, 10:43
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

#ifndef NODESDISTANCEMATRICES_H
#define	NODESDISTANCEMATRICES_H
#include "TreesManip.h"
#include <Bpp/Phyl/TreeTemplate.h>
#include <vector>
#include <list>
using namespace bpp;
using namespace std;

namespace tools {
/**
 * @brief The abstract class for nodal distances algorithms.
 */
class INodesDist {
protected:
    int k;
    vector <vector <double> > tr1NDists;
    vector <vector <double> > tr2NDists;
public:
    virtual void init(const TreeTemplate<Node> &tr1, const TreeTemplate<Node> &tr2) = 0;
    double getDistance();
};

/**
 * @brief Nodal distance algorithm for weighted trees.
 */
class WeigthedNodesDist : public INodesDist  {
private:
    void getTreeNodesDists(TreeTemplate<Node>& tr, vector<vector <double> >& trNDists);
public:
    WeigthedNodesDist(int kIn);
    void init(const TreeTemplate<Node> &tr1, const TreeTemplate<Node> &tr2);
};

/**
 * @brief Nodal distance algorithm for unweighted trees.
 */
class UnWeigthedNodesDist : public INodesDist  {
private:
    void getTreeNodesDists(TreeTemplate<Node>& tr, vector<vector <double> >& trNDists);
public:
    UnWeigthedNodesDist(int kIn);
    void init(const TreeTemplate<Node> &tr1, const TreeTemplate<Node> &tr2);
};

}
#endif	/* NODESDISTANCEMATRICES_H */

