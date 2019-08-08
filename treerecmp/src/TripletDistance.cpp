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

#include "TripletDistance.h"
using namespace bpp;
using namespace std;
namespace tools {
TreeParams::TreeParams(Node* root, int n, int l)
{  
    inSize = n;
    lSize = l;
    rootId = root->getId() - lSize;

    subTr = new int[inSize * 2];
    countSubTrSize(root); 

    inSons = new list<int>[inSize];
    setInSons(root, rootId);
}
TreeParams::~TreeParams()
{
    delete[] subTr;
    delete[] inSons;
}
int TreeParams::getResolved()
{  
    int result = 0;
    for (list<int>::iterator sIt = inSons[rootId].begin(); sIt != inSons[rootId].end(); sIt++) {
        result += getResolved_in(*sIt);
    }
    return result;
}
int TreeParams::getUnresolved(int R)
{  
    return choose3(lSize) - R;
}

int TreeParams::getResolved_in(int id)
{  
    int res = 0 ;
    res += choose2(subTr[id]) * (lSize - subTr[id]);
    for (list<int>::iterator sIt = inSons[id].begin(); sIt != inSons[id].end(); sIt++) {
        res -= choose2(subTr[*sIt]) * (lSize - subTr[id]);
        res += getResolved_in(*sIt);
    } 
    return res;
}
int TreeParams::choose2(int a)
{
    return a * (a-1) / 2;
}
int TreeParams::choose3(int a)
{
    return a * (a - 1) * (a - 2) / (2 * 3);
}

int TreeParams::countSubTrSize(Node* r)
{
    vector<Node*> s = r->getSons();
    if (s.size() == 0) {
        return 1;           // a leaf
    }
    int id = r->getId() - lSize;
    subTr[id] = 0;
    for (vector<Node*>::iterator it = s.begin(); it != s.end(); it++) {
        subTr[id] += countSubTrSize(*it);
    }
    return subTr[id];
}
//collecting sons being internal nodes. Attantion! Internal nodes have ids >= l
void TreeParams::setInSons(Node* root, int rootNewId)
{
    vector<Node*> sons = root->getSons();        
    for (vector<Node*>::iterator sonIt = sons.begin(); sonIt != sons.end(); sonIt++) {
        int sonId = (*sonIt)->getId() - lSize;
        if (sonId > 0) {
            inSons[rootNewId].push_back(sonId);
            setInSons(*sonIt, sonId);
        }          
    }
}


Triplets::Triplets(const TreeTemplate<Node>& tr1In, const TreeTemplate<Node>& tr2In)
{ 
    lSize = tr1In.getNumberOfLeaves();
    int inSize1 = tr1In.getNumberOfNodes() - lSize;
    Node* r1 = const_cast<Node*> (tr1In.getRootNode());
    trP1 = new TreeParams(r1, inSize1, lSize);
    int inSize2 = tr2In.getNumberOfNodes() - lSize;
    Node* r2 = const_cast<Node*> (tr2In.getRootNode());
    trP2 = new TreeParams(r2, inSize2, lSize);

    intersection = new int*[trP1->inSize * 2];
    for (int i = 0; i < inSize1 * 2; i++) {
        intersection[i] = new int[inSize2 * 2];
        for (int j = 0; j < inSize2 * 2 + 1; j++) {
            intersection[i][j] = -1;
        }
    }
    countIntersection(r1, r2, intersection);
}

Triplets::~Triplets()
{
    for (int i = 0; i < trP1->inSize * 2; i++) {
        delete[] intersection[i];
    }
    delete[] intersection;
    delete trP1;
    delete trP2;
}

int Triplets::getDistance()
{
    int R1 = trP1->getResolved();
    int U1 = trP1->getUnresolved(R1);
    int U2 = trP2->getUnresolved(trP2->getResolved());
    int S = getSameResolved();
    int Rr1 = getResolvedT1();
    return R1 - S + (U1 - U2) + Rr1;
}
     
int Triplets::countIntersection(Node* root1, Node* root2, int** intersection)
{
    int id1 = root1->getId() - lSize;
    int id2 = root2->getId() - lSize;
    bool internal = (min(id1, id2) >= 0) ? true : false;
    //if the intersection was already considered
    if (internal && intersection[id1][id2] != -1) return intersection[id1][id2];

    int val = 0;
    vector<Node*> s1 = root1->getSons();
    vector<Node*> s2 = root2->getSons();

    if (s1.size() == 0 && s2.size() == 0) {
        return id1 == id2;
    }

    for (vector<Node*>::iterator it1 = s1.begin(); it1 != s1.end(); it1++) {
        val += countIntersection(*it1, root2, intersection);
        for (vector<Node*>::iterator it2 = s2.begin(); it2 != s2.end(); it2++) {
            val -= countIntersection(*it1, *it2, intersection);            
        }
    }        
    for (vector<Node*>::iterator it2 = s2.begin(); it2 != s2.end(); it2++) {
        val += countIntersection(root1, *it2, intersection);
    } 

    if (internal) {
        intersection[id1][id2] = val;
    }
    return val;
}   

int Triplets::getSameResolved()
{  
    int res = 0;
    for (int i = 0; i < trP1->inSize; i++) {
        if (i != trP1->rootId) for (int j = 0; j < trP2->inSize; j++) {
            if (j != trP2->rootId) {
                int pairs = TreeParams::choose2(intersection[i][j]);
                for(list<int>::iterator itK = trP1->inSons[i].begin(); itK != trP1->inSons[i].end(); itK++) {
                    for(list<int>::iterator itL = trP2->inSons[j].begin(); itL != trP2->inSons[j].end(); itL++) {
                        pairs += TreeParams::choose2(intersection[*itK][*itL]);
                    }                        
                }                    
                for(list<int>::iterator itK = trP1->inSons[i].begin(); itK != trP1->inSons[i].end(); itK++) {
                    pairs -= TreeParams::choose2(intersection[*itK][j]);
                }
                for(list<int>::iterator itL = trP2->inSons[j].begin(); itL != trP2->inSons[j].end(); itL++) {
                    pairs -= TreeParams::choose2(intersection[i][*itL]);
                }

                pairs *= lSize - (trP1->subTr[i] + trP2->subTr[j] -  intersection[i][j]);
                res += pairs;
            }
        }
    }
    return res;
}

int Triplets::getResolvedT1()
{  
    int res = 0;
    for (int i = 0; i < trP1->inSize; i++) {
        if (i != trP1->rootId) for (int j = 0; j < trP2->inSize; j++) {
            res += ro(i, i, j);
            for(list<int>::iterator itK = trP1->inSons[i].begin(); itK != trP1->inSons[i].end(); itK++) {
                res -= ro(i, *itK, j);
            }

        }
    }
    return res;
}

int Triplets::ro(int u2, int u, int v)
{
    int res = TreeParams::choose2(intersection[u][v]) * b_nega(v, u2);
    for(list<int>::iterator itL = trP2->inSons[v].begin(); itL != trP2->inSons[v].end(); itL++) {
        res -= TreeParams::choose2(intersection[u][*itL]) * b_nega(*itL, u2);
        res -= TreeParams::choose2(intersection[u][*itL]) * (b_nega(v, u2) - b_nega(*itL, u2));
        res -= intersection[u][*itL] * b_nega(*itL, u2) * (intersection[u][v] - intersection[u][*itL]);
    }
    return res;
}
int Triplets::b_nega(int b, int a)
{
    return trP2->subTr[b] - intersection[a][b];
}


} // end of namespace

