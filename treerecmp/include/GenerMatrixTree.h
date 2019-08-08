//
// File: GenerMatrixTree.h
// Created by: Anna Pawelczyk
// Created on: 6 Apr 2011, 16:56
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

#ifndef GENERMATRIXTREE_H
#define	GENERMATRIXTREE_H

#include <Bpp/Phyl/TreeTemplate.h>
#include <list>
#include <cstring>
using namespace bpp;
using namespace std;
namespace tools {
/**
 * @brief Deprecated. Generational matrix for Triplets distance.
 */
class GenerMatrixTree {
private:
    //LeavesList ll;
    class LeavesList
    {
    private:
        struct Element
        {
        public:
            Element(int valIn)
            {
                val = valIn;
                next = NULL;
            }
            int val;
            Element* next;
        };

        Element* first;
        Element* last;
        int listSize;
        Element* iterator;
    public:
        LeavesList(int val)
        {
            first = new Element(val);
            last = first;
            listSize = 1;
        }
        ~LeavesList()
        {
            Element* el = first;
            while (el != NULL) {
                    Element* elNext = el->next;
                    delete el;
                    el = elNext;
            }
            delete el;
        }
        int size()
        {
                return listSize;
        }
        void concat(LeavesList* other)
        {
            last->next = other->first;
            last = other->last;
            listSize += other->size();
        }
        void resetIteration()
        {
            iterator = first;
        }
        int nextId()
        {
            if (iterator == NULL) throw Exception("LeavesList's iterator out of range");
            int val = iterator->val;
            iterator = iterator->next;
            return val;
        }
        bool iterationFinished()
        {
            return iterator == NULL;
        }
    };

    int size;
    int** matrix;
    int** helpDistMatrix; 
    int theFarestNode(int i, int j, int k);
public:
    GenerMatrixTree(const TreeTemplate<Node> &tr1);
    virtual ~GenerMatrixTree();
    int getTripletsDistance_binTrees(GenerMatrixTree& other) throw (Exception);
    int getTripletsDistance(GenerMatrixTree& other) throw (Exception); 
    void getCostReversedMatrixForPairs(GenerMatrixTree& other, vector<vector<int> >&);
    int fixReversedWeightsSum_Pairs(int dist);
    int at(int x, int y);
private:
    LeavesList* browseTree(int rootId, const Tree& tr, int level) throw (Exception);
    LeavesList* browseTree_Pair_incrInternNode(int rootId, const Tree& tr, int& nodeId);
    LeavesList* browseTree_Triplets_incrInternNode(int rootId, const Tree& tr, int& nodeId);

};
} // end of namespace

#endif	/* GENERMATRIXTREE_HPP */

