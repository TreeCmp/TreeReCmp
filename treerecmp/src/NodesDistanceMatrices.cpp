//
// File: NodesDistanceMatrices.cpp
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

#include "NodesDistanceMatrices.h"

namespace tools { 
    
double INodesDist::getDistance()
{        
    int size = tr1NDists.size();
    double distance = 0;
    for (int v1 = 0; v1 < size; v1++) {
        for (int v2 = v1+1; v2 < size; v2++) {
            double diff = std::abs(tr1NDists[v1][v2] - tr2NDists[v1][v2]);
            distance += std::pow(diff, k);
        }
    }
    distance = std::pow(distance, 1/(double)k);
    return distance;        
}



void WeigthedNodesDist::getTreeNodesDists(TreeTemplate<Node>& tr, vector<vector <double> >& trNDists)
{
    vector<Node*> leaves = tr.getLeaves();

    vector <double> tmp;
    tmp.resize(leaves.size(), 0);
    trNDists.resize(leaves.size(), tmp);

    vector<int> leavesIdLevel;
    leavesIdLevel.resize(leaves.size());
    TreesManip::getLeavesLevels(tr.getRootNode(), 0, leavesIdLevel);

    for (int i = 0; i < leaves.size(); i++) {
        for (int j = i + 1; j < leaves.size(); j++) {
            Node* ni = leaves[i];
            int idi = ni->getId();
            int li = leavesIdLevel[idi];
            Node* nj = leaves[j];
            int idj = nj->getId();
            int lj = leavesIdLevel[idj];

            if (li < lj) {
                while (li != lj) {
                    double val = nj->getDistanceToFather();
                    nj = nj->getFather();
                    lj--;
                    trNDists[idi][idj] += val;
                    trNDists[idj][idi] += val;
                }
            } else {
                while (li != lj) {
                    double val = nj->getDistanceToFather();
                    ni = ni->getFather();
                    li--;
                    trNDists[idi][idj] += val;
                    trNDists[idj][idi] += val;         
                }
            }

            while (ni != nj) {
                double val = nj->getDistanceToFather();
                val += ni->getDistanceToFather();
                ni = ni->getFather();
                nj = nj->getFather();
                trNDists[idi][idj] += val;
                trNDists[idj][idi] += val;
                li--;
            }
        }
    }
}
WeigthedNodesDist::WeigthedNodesDist(int kIn)
{
    k = kIn;
}
void WeigthedNodesDist::init(const TreeTemplate<Node> &tr1, const TreeTemplate<Node> &tr2)
{
    getTreeNodesDists(const_cast<TreeTemplate<Node>& > (tr1), tr1NDists);
    getTreeNodesDists(const_cast<TreeTemplate<Node>& > (tr2), tr2NDists);
}



void UnWeigthedNodesDist::getTreeNodesDists(TreeTemplate<Node>& tr, vector<vector <double> >& trNDists)
{
    vector<Node*> leaves = tr.getLeaves();

    vector <double> tmp;
    tmp.resize(leaves.size(), 0);
    trNDists.resize(leaves.size(), tmp);    

    vector<int> leavesIdLevel;
    leavesIdLevel.resize(leaves.size());
    TreesManip::getLeavesLevels(tr.getRootNode(), 0, leavesIdLevel);

    for (int i = 0; i < leaves.size(); i++) {
        for (int j = i + 1; j < leaves.size(); j++) {
            Node* ni = leaves[i];
            int idi = ni->getId();
            int li = leavesIdLevel[idi];
            Node* nj = leaves[j];
            int idj = nj->getId();
            int lj = leavesIdLevel[idj];

            if (li < lj) {
                while (li != lj) {                    
                    nj = nj->getFather();
                    lj--;
                    trNDists[idi][idj]++; 
                    trNDists[idj][idi]++;                    
                } 
            } else {
                while (li != lj) {
                    ni = ni->getFather();
                    li--;
                    trNDists[idi][idj]++; 
                    trNDists[idj][idi]++;                        
                } 
            }     

            while (ni != nj) {
                ni = ni->getFather();
                nj = nj->getFather();
                trNDists[idi][idj] += 2;
                trNDists[idj][idi] += 2;     
                li--;
            }            
        }
    }
}
UnWeigthedNodesDist::UnWeigthedNodesDist(int kIn)
{
    k = kIn;
}
void UnWeigthedNodesDist::init(const TreeTemplate<Node> &tr1, const TreeTemplate<Node> &tr2)
{
    getTreeNodesDists(const_cast<TreeTemplate<Node>& > (tr1), tr1NDists);
    getTreeNodesDists(const_cast<TreeTemplate<Node>& > (tr2), tr2NDists);
}
} // end of namespace

