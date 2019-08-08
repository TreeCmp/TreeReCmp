//
// File: Hungarian.h
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

#ifndef HUNGARIAN_H
#define	HUNGARIAN_H
#include <vector>
#include <climits>
#include <iostream>
using namespace std;

namespace tools {
/**
 * @brief Hungarian algorithm for minimum weight perfect matching of a bipartite graph.
 * \n Adapted from http://www.cse.ust.hk/~golin/COMP572/Notes/Matching.pdf/.
 */
class Hungarian
{
private:
    class Node
    {
    public:
        int id;
        int label;
        Node* match;
        bool alternatingTreeMember;
        int slack;	// min{l(x) + l(id) - w(x,id)} for those x in XaugPathMembers
        Node* slackNode;
        static Node* UNMATCHED;
        Node(int idIn)
        {
            nodeConstructor(idIn);
        }
        Node(int idIn, int labelIn)
        {
            nodeConstructor(idIn);
            label = labelIn;
        }
    private:
        void nodeConstructor(int idIn)
        {
            id = idIn;
            label = 0;
            match = UNMATCHED;
            alternatingTreeMember = false;
            slack = INT_MAX;
            slackNode = NULL;
        }
    };

    class BipartGraphToMatch
    {
    private:
        int matchingSize;
        vector<vector<int> > costs;
        vector<Node*> X;
        vector<Node*> Y;
        vector<Node*> XalternatingTreeMembers;	// so-called in papers the S set

    public:
        BipartGraphToMatch(const vector<vector<int> >& costsIn);
        ~BipartGraphToMatch() 
        {
            for (int i = 0; i < Y.size(); i++) delete Y.at(i);
            for (int i = 0; i < X.size(); i++) delete X.at(i);
        };

        int getMatchingWeight();

        bool isPerfectMatchingFound() {	return matchingSize == X.size(); }
        /**
         * Checking if there exists at least one node in Y that
         * - is a neighbour (in feasible labeling sense)of S
         * - does not belong to T
        */
        Node* getAugumentingPathYNode();

        void updateLabeling();

        void extendAlternatingTree(Node* nY);

        void augumentMatching(Node* nY);

        int getY(int x) {return X[x]->match->id;}

        void reset();
    private:
        void clearAlternatingTree();

        void setInitialMatching();

        void setYSlack(Node* x, Node* y);
    };

public:
    static int MaxWPerfectMatchingCost(const vector<vector<int> >& costsIn);
};

} //end of namespace
#endif	/* HUNGARIAN_H */

