//
// File: PhylotreeDist.cpp
// Created by: Anna Pawelczyk
// Created on: 19 May 2011, 12:10
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

#include "PhylotreeDist.h"
namespace dist {

bool PhylotreeDist::checkLeavesNames(const TreeTemplate<Node>& trIn1, const TreeTemplate<Node>& trIn2)
            throw (bpp::Exception)
{
    if(!VectorTools::haveSameElements(trIn1.getLeavesNames(), trIn2.getLeavesNames()))
            throw Exception("Trees have different sets of leaves.\n");
    return true;
}
bool PhylotreeDist::checkRooted(bool condition, const TreeTemplate<Node>& trIn1, const TreeTemplate<Node>& trIn2)
            throw (bpp::Exception)
{
    if(condition) {
        if (!trIn1.isRooted()) throw Exception("Unrooted tree trIn1. Trees must be rooted.");
        if (!trIn2.isRooted()) throw Exception("Unrooted tree trIn2. Trees must be rooted.");            
    } else {
        if (trIn1.isRooted()) throw Exception("Rooted tree trIn1. Trees must be unrooted.");
        if (trIn2.isRooted()) throw Exception("Rooted tree trIn2. Trees must be unrooted.");
    }
    return true;
}

int PhylotreeDist::robinsonFoulds(const TreeTemplate<Node>& trIn1, const TreeTemplate<Node>& trIn2, bool setNodesId, bool checkNames) 
            throw (bpp::Exception)
{
    //Constraints assurance  
    if (trIn1.isRooted()  ^ trIn2.isRooted()) throw Exception("Bad input trees. Both tree must be either rooted or unrooted.");  
    if (checkNames) {
        checkLeavesNames(trIn1, trIn2); 
    }
    const TreeTemplate<Node> *tr1 = setNodesId ? tools::TreesManip::createOrderedTrees(trIn1) : &trIn1;
    const TreeTemplate<Node> *tr2 = setNodesId ? tools::TreesManip::createOrderedTrees(trIn2) : &trIn2;

    //The algorithm  
    if (!trIn2.isRooted()) {
        if (checkNames) {
        }
            tr1 = tr1->clone();
            tr2 = tr2->clone();
        // Pay attention - changing tree to speed the algorithm
        if (!trIn1.isRooted()) PostorderTree::rootTree(const_cast<TreeTemplate<Node>&>(*tr1));
        if (!trIn2.isRooted()) PostorderTree::rootTree(const_cast<TreeTemplate<Node>&>(*tr2));
    }      
    PostorderTree pTrA(tr1->getRootNode());
    PostorderTree pTrB(tr2->getRootNode());
    ClusterTable *clusters = new ClusterTable(tr1->getNumberOfLeaves(), pTrA);
    clusters->removeUncommonElements(pTrB);
    return tr1->getNumberOfNodes() + tr2->getNumberOfNodes() - 2 * tr1->getNumberOfLeaves()
            - 2 * clusters->getNumberOfInternalNodes();
}   
double PhylotreeDist::robinsonFouldsW(const TreeTemplate<Node>& trIn1, const TreeTemplate<Node>& trIn2, bool setNodesId, bool checkNames) 
            throw (bpp::Exception)
{
    //Constraints assurance  
    checkRooted(false, trIn1, trIn2);
    if (checkNames) {
        checkLeavesNames(trIn1, trIn2); 
    }
    const TreeTemplate<Node> *tr1 = setNodesId ? tools::TreesManip::createOrderedTrees(trIn1) : &trIn1;
    const TreeTemplate<Node> *tr2 = setNodesId? tools::TreesManip::createOrderedTrees(trIn2) : &trIn2;

    //The algorithm  
    if (checkNames) {
        tr1 = tr1->clone();
        tr2 = tr2->clone();
    }
    // Pay attention - changing tree to speed the algorithm
    if (!trIn1.isRooted()) PostorderTree::rootTree(const_cast<TreeTemplate<Node>&>(*tr1));
    if (!trIn2.isRooted()) PostorderTree::rootTree(const_cast<TreeTemplate<Node>&>(*tr2));

    PostorderTree pTrA(tr1->getRootNode());
    PostorderTree pTrB(tr2->getRootNode());
    ClusterTable *clusters = new ClusterTable(tr1->getNumberOfLeaves(), pTrA);
    clusters->removeUncommonElements(pTrB);
    double weight = clusters->getWeight();

    vector<const Node*> l1 = tr1->getLeaves();
    vector<const Node*> l2 = tr2->getLeaves();
    vector<double> weights;
    weights.resize(l1.size(), 0);
    for(vector<const Node*>::iterator itN = l1.begin(); itN != l1.end(); itN++) {
        weights[(*itN)->getId()] = (*itN)->getDistanceToFather();
    }
    for(vector<const Node*>::iterator itN = l2.begin(); itN != l2.end(); itN++) {
        weight += abs( weights[(*itN)->getId()] - (*itN)->getDistanceToFather());
    }        

    return weight;
}   


int PhylotreeDist::perfectMatching_splits(const TreeTemplate<Node>& trIn1, const TreeTemplate<Node>& trIn2, bool setLeavesId, bool checkNames)
    throw (Exception)
{
    checkRooted(false, trIn1, trIn2);
    if (checkNames) {      
        checkLeavesNames(trIn1, trIn2);
    }
    const TreeTemplate<Node> *tr1 = setLeavesId ? tools::TreesManip::createOrderedTrees(trIn1) : &trIn1;
    const TreeTemplate<Node> *tr2 = setLeavesId ? tools::TreesManip::createOrderedTrees(trIn2) : &trIn2;

    Splitting splitting(*tr1, *tr2);
    return getPMDistance(splitting);
}

int PhylotreeDist::perfectMatching_clusters(const TreeTemplate<Node>& trIn1, const TreeTemplate<Node>& trIn2, bool setLeavesId, bool checkNames)
    throw (Exception)
{
    checkRooted(true, trIn1, trIn2);
    if (checkNames) {
        checkLeavesNames(trIn1, trIn2);
    }        
    const TreeTemplate<Node> *tr1 = setLeavesId ? tools::TreesManip::createOrderedTrees(trIn1) : &trIn1;
    const TreeTemplate<Node> *tr2 = setLeavesId ? tools::TreesManip::createOrderedTrees(trIn2) : &trIn2;

    Clustering clustering(*tr1, *tr2);
    return getPMDistance(clustering);
}

int PhylotreeDist::perfectMatching_pairs(const TreeTemplate<Node>& trIn1, const TreeTemplate<Node>& trIn2, bool setLeavesId, bool checkNames)
    throw (Exception)
{	            
    checkRooted(true, trIn1, trIn2);
    if(checkNames) {
        if (trIn1.isMultifurcating()) throw Exception("Multifurcating tree trIn1. Trees must be bifurcating.");
        if (trIn2.isMultifurcating()) throw Exception("Multifurcating tree trIn2. Trees must be bifurcating.");

        checkLeavesNames(trIn1, trIn2);
    }
    const TreeTemplate<Node> *tr1 = setLeavesId ? tools::TreesManip::createOrderedTrees(trIn1) : &trIn1;
    const TreeTemplate<Node> *tr2 = setLeavesId ? tools::TreesManip::createOrderedTrees(trIn2) : &trIn2;

    PairLeavesSets pairs(*tr1, *tr2);           
    return getPMDistance(pairs);
}

int PhylotreeDist::getPMDistance(ITwoTreesDescriptionElements& descriptionElements)
{	
int distance;

#ifdef HUNGARIAN
    vector<vector<int> > costMatrix;
    descriptionElements.getCostReversedMatrix(costMatrix);
    distance = Hungarian::MaxWPerfectMatchingCost(costMatrix);
    return descriptionElements.fixReversedWeightsSum(distance);  
#else              
    int size = descriptionElements.getSize(); 
    int **assigncost = new int*[size];
    for (int a = 0; a < size; a++) {
        assigncost[a] = new int[size];
    }
    int *colsol = new int[size], *rowsol = new int[size];
    int *u = new int[size], *v = new int[size];

    // The exact functionality
    descriptionElements.getCostMatrix(assigncost);  
    distance = lap(size, assigncost, rowsol, colsol, u, v ); 
    // The end of exact functionality

    for (int i = 0; i < size; i++) 
        delete[] assigncost[i];
    delete[] assigncost;
    delete[] colsol;
    delete[] rowsol;
    delete[] u;
    delete[] v;
    return distance;
#endif
}



double PhylotreeDist::getNodalDistance(INodesDist* d, const TreeTemplate<Node>& trIn1, const TreeTemplate<Node>& trIn2, bool setLeavesId, bool checkNames)
    throw (Exception)
{
    checkRooted(false, trIn1, trIn2); 
    if(checkNames) {

        checkLeavesNames(trIn1, trIn2);
    }
    const TreeTemplate<Node> *tr1 = setLeavesId ? tools::TreesManip::createOrderedTrees(trIn1) : &trIn1;
    const TreeTemplate<Node> *tr2 = setLeavesId ? tools::TreesManip::createOrderedTrees(trIn2) : &trIn2;
    d->init(trIn1, trIn2);
    return d->getDistance();
}

int PhylotreeDist::nodalDistance(const TreeTemplate<Node>& trIn1, const TreeTemplate<Node>& trIn2, bool setLeavesId, bool checkNames)
    throw (Exception)
{
    UnWeigthedNodesDist d(1);
    return getNodalDistance(&d, trIn1, trIn2, setLeavesId, checkNames);
}
double PhylotreeDist::nodalDistance_pythagorean(const TreeTemplate<Node>& trIn1, const TreeTemplate<Node>& trIn2, bool setLeavesId, bool checkNames)
    throw (Exception)
{    
    UnWeigthedNodesDist d(2);
    return getNodalDistance(&d, trIn1, trIn2, setLeavesId, checkNames);
}
double PhylotreeDist::nodalDistanceW(const TreeTemplate<Node>& trIn1, const TreeTemplate<Node>& trIn2, bool setLeavesId, bool checkNames)
    throw (Exception)
{    
    WeigthedNodesDist d(1);
    return getNodalDistance(&d, trIn1, trIn2, setLeavesId, checkNames);
}
double PhylotreeDist::nodalDistanceW_pythagorean(const TreeTemplate<Node>& trIn1, const TreeTemplate<Node>& trIn2, bool setLeavesId, bool checkNames)
    throw (Exception)
{    
    WeigthedNodesDist d(2);
    return getNodalDistance(&d, trIn1, trIn2, setLeavesId, checkNames);
}


int PhylotreeDist::quartetDistance(const TreeTemplate<Node>& trIn1, const TreeTemplate<Node>& trIn2, bool setLeavesId, bool checkNames)
    throw (Exception)
{
    checkRooted(false, trIn1, trIn2); 
    if(checkNames) {
        checkLeavesNames(trIn1, trIn2);
    }
    const TreeTemplate<Node> *tr1 = setLeavesId ? tools::TreesManip::createOrderedTrees(trIn1) : &trIn1;
    const TreeTemplate<Node> *tr2 = setLeavesId ? tools::TreesManip::createOrderedTrees(trIn2) : &trIn2;
    QuartetDistance q(*tr1, *tr2);
    return q.getDistance();
}

int PhylotreeDist::tripletsDistance(const TreeTemplate<Node>& trIn1, const TreeTemplate<Node>& trIn2, bool setLeavesId, bool checkNames)
            throw (Exception)
{            
    checkRooted(true, trIn1, trIn2);                
    if(checkNames) {
        //if (trIn1.isMultifurcating()) throw Exception("Multifurcating tree trIn1. Trees must be bifurcating.");
        //if (trIn2.isMultifurcating()) throw Exception("Multifurcating tree trIn2. Trees must be bifurcating.");   
        checkLeavesNames(trIn1, trIn2);
    }
    const TreeTemplate<Node> *tr1 = setLeavesId ? tools::TreesManip::createOrderedTrees(trIn1) : &trIn1;
    const TreeTemplate<Node> *tr2 = setLeavesId ? tools::TreesManip::createOrderedTrees(trIn2) : &trIn2;

    Triplets t(*tr1, *tr2);

#ifdef TRIPL_NAIVE
    bool tr1IsBifurc = !trIn1.isMultifurcating() && !trIn2.isMultifurcating();
    GenerMatrixTree mTr1(*tr1);
    GenerMatrixTree mTr2(*tr2);
    if (tr1IsBifurc) {
        return mTr1.getTripletsDistance_binTrees(mTr2);
    } else {
        return mTr1.getTripletsDistance(mTr2);
    }
#endif
    return t.getDistance();
    
}


} // end of namespace
