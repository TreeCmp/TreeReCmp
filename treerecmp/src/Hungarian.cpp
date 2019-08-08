//
// File: Hungarian.cpp
// Created by: Anna Pawelczyk
// Created on: 23 May 2011, 10:53
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

#include "Hungarian.h"
namespace tools {
    
int matchSize = 0;
Hungarian::Node* Hungarian::Node::UNMATCHED = NULL;

int Hungarian::MaxWPerfectMatchingCost(const vector<vector<int> >& costsIn)
{
    //STEP 1
    BipartGraphToMatch bipartGraphToMatch(costsIn);
    bool machingFound;

    while (! bipartGraphToMatch.isPerfectMatchingFound()) {
        machingFound = false;
//STEP 2
        bipartGraphToMatch.reset();
        while (!machingFound) {
//STEP 3 & 4
            bool labelsUpdated = false;
            Node *nY = bipartGraphToMatch.getAugumentingPathYNode();
            if (nY == NULL) {
                bipartGraphToMatch.updateLabeling();
                labelsUpdated= true;
            } else {
                if (nY->match == Node::UNMATCHED) {
                    bipartGraphToMatch.augumentMatching(nY);
                        machingFound = true;
                } else {
                    bipartGraphToMatch.extendAlternatingTree(nY);
                }
            }
        }
    }
    return bipartGraphToMatch.getMatchingWeight();
}


Hungarian::BipartGraphToMatch::BipartGraphToMatch(const vector<vector<int> >& costsIn) : costs(costsIn)
{
    matchingSize = 0;
    for (int i = 0; i < costs.size(); i++) {
        Y.push_back(new Node(i));
    }
    setInitialMatching();
}


int Hungarian::BipartGraphToMatch::getMatchingWeight()
{
    int weight = 0;
    for (int i = 0; i < costs.size(); i++) weight += X[i]->label + X[i]->match->label;
    return weight;
}

Hungarian::Node* Hungarian::BipartGraphToMatch::getAugumentingPathYNode()
{
    Node *y = NULL;

    for (vector<Node*>::iterator xIt = XalternatingTreeMembers.begin(); xIt != XalternatingTreeMembers.end(); xIt++) {
        for (int yId = 0; yId < costs.size(); yId++) {
            if (costs[ (*xIt)->id ][yId] == (*xIt)->label + Y[yId]->label // check if y is a neighbour(in feasible labeling sense)of x
              && !Y[yId]->alternatingTreeMember) {	
            //an exposed vertex in Y found, so augmenting path exists
                y = Y[yId];
                break;
            }
        }
        if (y != NULL) break;
    }
    return y;
}

void Hungarian::BipartGraphToMatch::updateLabeling()
{
    int minCost = INT_MAX;
    for (int i = 0; i < Y.size(); i++) {				
        if (Y[i]->slack < minCost && Y[i]->alternatingTreeMember == false) {
            minCost = Y[i]->slack;
        }
    }
    for (int i = 0; i < Y.size(); i++) {
        if (X[i]->alternatingTreeMember) {
            X[i]->label -= minCost;
        }
        if (Y[i]->alternatingTreeMember) {
            Y[i]->label += minCost;
        } else {
            Y[i]->slack -= minCost;
        }
     }
}
void Hungarian::BipartGraphToMatch::extendAlternatingTree(Hungarian::Node* nY)
{
    nY->alternatingTreeMember = true;
    if (nY->match->alternatingTreeMember == false) {
        nY->match->alternatingTreeMember = true;
        XalternatingTreeMembers.push_back(nY->match);
    }
    nY->match->slackNode = nY;
    // Updating slacks, cause a new member was added to XaugPathMembers
    int xId = nY->match->id;
    for (int yId = 0; yId < Y.size(); yId++) {
        if (Y[yId]->slack > X[xId]->label + Y[yId]->label - costs[xId][yId]) 
            setYSlack(nY->match, Y[yId]);
    }
}

void Hungarian::BipartGraphToMatch::augumentMatching(Hungarian::Node* nY)
{
    Node* nX;
    matchingSize++;
    do {
        nX = nY->slackNode;
        nX->match = nY;
        nY->match = nX;
        nY = nX->slackNode;
    } while (nX != XalternatingTreeMembers[0]);
}

void Hungarian::BipartGraphToMatch::reset()
{
    clearAlternatingTree();
    // Picking unmatched vertex
    Node* xNode = NULL;
    for (int i = 0; xNode == NULL; i++ ) {
        if (X.at(i)->match == Node::UNMATCHED) xNode = X.at(i);
    }
    XalternatingTreeMembers.push_back(xNode);
    xNode->alternatingTreeMember = true;
    for (int yId = 0; yId < Y.size(); yId++) setYSlack(xNode, Y[yId]);
}

void Hungarian::BipartGraphToMatch::clearAlternatingTree()
{
    for (int i = 0; i < X.size(); i++) {
        X[i]->alternatingTreeMember = false;
        Y[i]->alternatingTreeMember = false;
    }
    XalternatingTreeMembers.clear();
}
void Hungarian::BipartGraphToMatch::setInitialMatching()
{
    for (int xId = 0; xId < costs.size(); xId++) {
        int maxL = 0, maxYLId;
        for (int yId = 0; yId < costs.size(); yId++) {
            if (maxL < costs[xId][yId]) {
                maxL = costs[xId][yId];
                maxYLId = yId;
            } else if (maxL == costs[xId][yId] 
              && Y[maxYLId]->match != Node::UNMATCHED 
              && Y[yId]->match == Node::UNMATCHED) {
                maxYLId = yId;
            }
        }
        X.push_back(new Node(xId, maxL));
        // Setting matching conditionally
        if (Y[maxYLId]->match == Node::UNMATCHED) {
            X[xId]->match = Y[maxYLId];
            Y[maxYLId]->match = X[xId];
            matchingSize++;
        }
    }
}
void Hungarian::BipartGraphToMatch::setYSlack(Hungarian::Node* x, Hungarian::Node* y)
{
    y->slack = x->label + y->label - costs[x->id][y->id];
    y->slackNode = x;
}

} // end of namespace
