//
// File: GenerMatrixTree.cpp
// Created by: Anna Pawelczyk
// Created on: 23 May 2011, 17:54
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

#include "GenerMatrixTree.h"
namespace tools {
    
GenerMatrixTree::GenerMatrixTree(const TreeTemplate<Node> &tr)
{
    size = tr.getNumberOfLeaves();
    matrix = new int*[size];
    helpDistMatrix = new int*[size];
    for (int i = 0; i < size; i++) {
        matrix[i] = new int[size];
        helpDistMatrix[i] = new int[size];
    }

    LeavesList* l = browseTree(tr.getRootId(), tr, 1);
    delete l;
}

GenerMatrixTree::~GenerMatrixTree()
{
    for (int i = 0; i < size; i++) {
        delete[] matrix[i];
        delete[] helpDistMatrix[i];
    }
    delete[] matrix;
    delete[] helpDistMatrix;
}


GenerMatrixTree::LeavesList* GenerMatrixTree::browseTree(int rootId, const Tree& tr, int level)
        throw (Exception)
{
    vector<int> sons = tr.getSonsId(rootId);
    if (sons.size() == 0) {			// a leaf
        return new LeavesList(rootId);		
    } else {				//2 sons in bifurcating tree
        //if (sons.size() != 2) throw Exception("Trees must be bifurcating. An internal node is multifurcating");
        LeavesList* l1 = browseTree(sons.at(0), tr, level + 1);
        for (int i = 1; i <sons.size(); i++) {
            LeavesList* l2 = browseTree(sons.at(i), tr, level + 1);
            for (l1->resetIteration(); !l1->iterationFinished();) {
                int leafIdTr1 = l1->nextId();
                for (l2->resetIteration(); !l2->iterationFinished();) {
                        int leafIdTr2 = l2->nextId();
                        matrix[leafIdTr1][leafIdTr2] = level;
                        matrix[leafIdTr2][leafIdTr1] = level;
                }
            }
            l1->concat(l2);
        }
        return l1;
    }
}

int GenerMatrixTree::at(int x, int y)
{
    return matrix[x][y];
}

int GenerMatrixTree::getTripletsDistance_binTrees(GenerMatrixTree& other) 
        throw (Exception)
{
    int distance = size * (size - 1) * (size - 2) / 6;      // (size over 3) - the number of triplets in a tree
    for (int line = 0; line < size; line++) {
        for (int i = 0; i < size; i++) memset(helpDistMatrix[i], 0, size * sizeof(int));        ////????????????????????? w forze?
        list<int*> usedPatterns;

        for (int i = 0; i < size; i++) {
            if (line != i) {
                int patternPartX = this->at(line, i);
                int patternPartY = other.at(line, i);
                helpDistMatrix[patternPartX][patternPartY]++;
                if (helpDistMatrix[patternPartX][patternPartY] == 2) {
                    usedPatterns.push_back(& helpDistMatrix[patternPartX][patternPartY]);
                }
            }
        }
        for (list<int*>::iterator it = usedPatterns.begin(); it != usedPatterns.end(); it++) {
            distance -= (**it) * (**it - 1) / 2;
            **it = 0;
        }
    }
    return distance;
}
int GenerMatrixTree::getTripletsDistance(GenerMatrixTree& other) 
        throw (Exception)
{
    int distance = 0;      // (size over 3) - the number of triplets in a tree
    for (int i = 0; i < size - 2; i++) {		 
        for (int j = i+1; j < size - 1; j++) {
            for (int k = j+1; k < size; k++) {
                if (this->theFarestNode(i, j, k) != other.theFarestNode(i, j, k)) {
                    distance++;
                }
            }
        }
    }
    return distance;
}

int GenerMatrixTree::theFarestNode(int i, int j, int k)
{
    int ij = matrix[i][j];
    int jk = matrix[j][k];
    int ik = matrix[i][k];
    if (ij == jk && jk == ik) return -1;
    if (ij == jk)  return j;
    if (ij == ik)  return i;
    if (ik == jk)  return k;
}


} // end of namespace
