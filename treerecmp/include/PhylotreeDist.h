//
// File: PhylotreeDist.h
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

#ifndef PHYLOTREEDIST_H
#define	PHYLOTREEDIST_H

#include <Bpp/Phyl/TreeTemplate.h>
using namespace bpp;
#include "TreesManip.h"
#include "PostorderTree.h"
#include "GenerMatrixTree.h"
#include "ClusterTable.h"
#include "Partitioning.h"
#include "QuartetDistance.h"
#include "TripletDistance.h"
#include "NodesDistanceMatrices.h"

#include "Hungarian.h"
#include "hungarianJV/lap.h"
//#include "QPartitionList.h"
using namespace tools;

namespace dist {
/**
 * @brief Static methods to count different distances between two phylometric trees.
 */
class PhylotreeDist {    
public:
    /**
     * @brief The RobinsonFoulds distance between two unrooted or rooted, multifurcating trees with the same set of leaves.
     * \n The RobinsonFoulds metric bases on counting occurrences of:
     * \n     - for unrooted trees: bipartitions not common to both trees
     * \n     - for rooted trees: clusters not common to both trees
     * \n The constraints for the two input trees: Either the same leaves id sets numbered 0..n-1 and internal nodes ids numbered n, n+1,...  or the same leaves name sets and setLeavesId parameter true
     * \n\n Time complexity: O(n)
     * \n Algorithm adapted from William H. E. Day, "Optimal algorithms for comparing trees with labeled leaves", "Journal of Classification", 1984
     * 
     * @param[in]   tr1     First tree.
     * @param[in]   tr2     Second tree.
     * @param[in]   setNodesId (optional) TRUE if the two trees do not have the same ids for the same leaves or the ids are not numbered 0..n-1. 
     * Also true if the trees do not have the internal nodes idsas a sequence n, n+1, n+2,... Defaults to TRUE.
     * @param[in]   checkNames (optional) TRUE if check whether trees have the same leaves set. Defaults to FALSE.
     * @return      Robinson-Foulds distance
     * @throw bpp::Exception if trees have different leaves sets.
     */

    static int robinsonFoulds(const TreeTemplate<Node>& trIn1, const TreeTemplate<Node>& trIn2, bool setNodesId = true, bool checkNames = false) 
            throw (bpp::Exception);	
    /**
     * @brief The RobinsonFoulds distance between two branch-weighted unrooted multifurcating trees with the same set of leaves.
     * \n The RobinsonFoulds metric bases on counting for every tree for every bipartition on one tree that does not occure in the other tree - weights of the branches which removal leads to that bipartition. 
     * Then adding for every bopartition that occures in both trees - the difference between weights of the T1 and T2 branches which removal leads to that bipartition. 
     * 
     * \n The constraints for the two input trees: Either the same leaves id sets numbered 0..n-1 and internal nodes ids numbered n, n+1,...  or the same leaves name sets and setLeavesId parameter true
     * \n\n Time complexity: O(n)
     * \n Algorithm adapted from William H. E. Day, "Optimal algorithms for comparing trees with labeled leaves", "Journal of Classification", 1984
     * 
     * @param[in]   tr1     First unrooted weighted tree.
     * @param[in]   tr2     Second unrooted weighted tree.
     * @param[in]   setNodesId (optional) TRUE if the two trees do not have the same ids for the same leaves or the ids are not numbered 0..n-1. 
     * Also true if the trees do not have the internal nodes idsas a sequence n, n+1, n+2,... Defaults to TRUE.
     * @param[in]   checkNames (optional) TRUE if check whether trees have the same leaves set. Defaults to FALSE.
     * @return      Robinson-Foulds distance
     * @throw bpp::Exception if trees have different leaves sets.
     */
    static double  robinsonFouldsW(const TreeTemplate<Node>& trIn1, const TreeTemplate<Node>& trIn2, bool setNodesId = true, bool checkNames = false) 
            throw (bpp::Exception); 


    /**
     * @brief The Minimum Weight Perfect Matching distance between two unrooted trees with the same set of leaves
     * \n  with splits as description elements.
     * \n The algorithm performs a minimum weight perfect matching on a weighted bipartite graph, where a graph's partition corresponds to a description elements set of a tree - a set of all the tree's splits (bipartitions).
     * \n The weight on an edge that links two splits A1|B1 and A2|B2 is count as  h(A_1 | B_1, A_2 | B_2) = min(|A_1 xor A_2|, |A_1 xor B_2|) = min(|A_1| + |A_2| - 2|A_1 and A_2|, n - (|A_1|+|A_2|-2|A_1 and A_2|).

     * \n The constraints for the two input trees: Either the same leaves id sets numbered 0..n-1 and internal nodes ids numbered n, n+1,...  or the same leaves name sets and setLeavesId parameter true
     * \n\n Time complexity: O(n^3)
     * \n Algorithm adapted from D. Bogdanowicz and K. Giaro. Comparing arbitrary unrooted phylogenetic trees using generalized matching split distance. In Information Technology (ICIT), 2010 2nd International Conference on, pages 259–262. IEEE, 2010.
     *
     * \n By default the minimum weight perfect matching is computed by an external library that bases on Jonker and Volgenant algorithm. If during compilation there is a preprocessor HUNGARIAN symbol defined (-DHUNGARIAN), then there is a Knuth’s Hungarian algorithm used
     *  
     * @param[in]   tr1     First unrooted tree.
     * @param[in]   tr2     Second unrooted tree.
     * @param[in]   setLeavesId (optional) TRUE if the two trees do not have the same ids for the same leaves or the ids are not numbered 0..n-1.
     * @param[in]   checkNames (optional) TRUE if check whether trees have the same leaves set and whether is unrooted. Defaults to FALSE.
     * @return      Robinson-Foulds distance
     * @throw bpp::Exception if trees have different leaves sets.
     */
    static int perfectMatching_splits(const TreeTemplate<Node>& trIn1, const TreeTemplate<Node>& trIn2, bool setLeavesId = true, bool checkNames = false)
            throw (bpp::Exception);

    /**
     * @brief The Minimum Weight Perfect Matching distance between two rooted trees with the same set of leaves
     * \n  with clusters as description elements.
     * \n The algorithm performs a minimum weight perfect matching on a weighted bipartite graph, where a graph's partition corresponds to a description elements set of a tree - a set of all the tree's clusters.
     * \n The weight on an edge that links two clusters A1 and A2 is count as h(A1 , A2 ) = |A1 xor A2 | = |A1 | + |A2 | - 2|A1 and A2 |.
     * 
     * \n The constraints for the two input trees: Either the same leaves id sets numbered 0..n-1 and internal nodes ids numbered n, n+1,...  or the same leaves name sets and setLeavesId parameter true
     * \n\n Time complexity: O(n^3)
     * \n Algorithm main assumptions base on D. Bogdanowicz and K. Giaro. Comparing arbitrary unrooted phylogenetic trees using generalized matching split distance. In Information Technology (ICIT), 2010 2nd International Conference on, pages 259–262. IEEE, 2010.
     *
     * \n By default the minimum weight perfect matching is computed by an external library that bases on Jonker and Volgenant algorithm. If during compilation there is a preprocessor HUNGARIAN symbol defined (-DHUNGARIAN), then there is a Knuth’s Hungarian algorithm used

     * \n The minimum weight perfect matching algorithm is a Knuth's Hyngarian algorithm adapted from Mordecai J. Golin which implementation base on Bipartite matching and the hungarian method. Hong Kong University of Science and Technology Course Notes http://www.cse.ust.hk/~golin/COMP572/Notes/Matching.pdf/.
     * \n Simultaneously there is a Jonker and Volgenant algorithm that solves the linear assignment problem: there is used an external library from http://www.assignmentproblems.com/LAPJV.htm.
     * 
     * @param[in]   tr1     First rooted tree.
     * @param[in]   tr2     Second rooted tree.
     * @param[in]   setLeavesId (optional) TRUE if the two trees do not have the same ids for the same leaves or the ids are not numbered 0..n-1. Defaults to TRUE.
     * @param[in]   checkNames (optional) TRUE if check whether trees have the same leaves set and whether is rooted. Defaults to FALSE.
     * @return      Robinson-Foulds distance
     * @throw bpp::Exception if trees have different leaves sets.
     */
    static int perfectMatching_clusters(const TreeTemplate<Node>& trIn1, const TreeTemplate<Node>& trIn2, bool setLeavesId = true, bool checkNames = false)
            throw (bpp::Exception);	

    /**
     * @brief The Minimum Weight Perfect Matching distance between two rooted trees with the same set of leaves
     * \n  with set of pairs as description elements.
     * \n The algorithm performs a minimum weight perfect matching on a weighted bipartite graph, where a graph's partition corresponds to a description elements set of a tree. 
     * A description element corresponds to an internal node and is a set of all the pairs of leaves to which the node is the last common ancestor. 
     * \n The weight on an edge that links two sets of pairs A1 and A2 is the number of pairs that occur in exactly one of the two description elements’ set of unordered pairs. 
     * Formally, h(A1 , A2 ) = |(A1 or A2 ) \ (A1 and A2 |) = |A1 xor A2 |

     * \n The constraints for the two input trees: Either the same leaves id sets numbered 0..n-1 and internal nodes ids numbered n, n+1,...  or the same leaves name sets and setLeavesId parameter true
     * \n\n Time complexity: O(n^3)
     * \n Algorithm main assumptions base on D. Bogdanowicz and K. Giaro. Comparing arbitrary unrooted phylogenetic trees using generalized matching split distance. In Information Technology (ICIT), 2010 2nd International Conference on, pages 259–262. IEEE, 2010.
     * \n By default the minimum weight perfect matching is computed by an external library that bases on Jonker and Volgenant algorithm. If during compilation there is a preprocessor HUNGARIAN symbol defined (-DHUNGARIAN), then there is a Knuth’s Hungarian algorithm used	 
     *
     * \n The minimum weight perfect matching algorithm is a Knuth's Hyngarian algorithm
                             * adapted from Hong Kong University of Science and Technology Course Notes http://www.cse.ust.hk/~golin/COMP572/Notes/Matching.pdf/.
                             * Simltaneously there is used an external library from http://www.assignmentproblems.com/LAPJV.htm. It bases on Jonker and Volgenant algorithm that solves the linear assignment problem.
     * 
     * @param[in]   tr1     First rooted binary tree.
     * @param[in]   tr2     Second rooted binary tree.
     * @param[in]   setLeavesId (optional) TRUE if the two trees do not have the same ids for the same leaves or the ids are not numbered 0..n-1. Defaults to TRUE.
     * @param[in]   checkNames (optional) TRUE if check whether trees have the same leaves set and whether is binary, rooted. Defaults to FALSE.
     * @return      minimum weight perfect matching split distance
     * @throw bpp::Exception if trees have different leaves sets.
     */
    static int perfectMatching_pairs(const TreeTemplate<Node>& trIn1, const TreeTemplate<Node>& trIn2, bool setLeavesId = true, bool checkNames = true)
            throw (bpp::Exception);


     /**
     * @brief The  Quartet distance between two unrooted trees with the same set of leaves.
     * \n Counts the number of sets of four species for which the topologies differ in the two trees.
     * 
     * \n The constraints for the two input trees: Either the same leaves id sets numbered 0..n-1 and internal nodes ids numbered n, n+1,...  or the same leaves name sets and setLeavesId parameter true
     * \n\n Time complexity: O(n + |V1||V2| deg1max ^ 2 deg2max) where V is the number of internal nodes, degmax is the maximal internal node degree.
     * \n Algorithm adapted from C. Christiansen, T. Mailund, C. Pedersen, and M. Randers. Computing the quartet distance between trees of arbitrary degree. Algorithms in Bioinformatics, pages 77–88, 2005.
     *
     * @param[in] tr1 First unrooted tree.
     * @param[in] tr2 Second unrooted tree.
     * @param[in] setLeavesId (optional) TRUE if the two trees do not have the same ids for the same leaves or the ids are not numbered 0..n-1. 
     * @param[in] checkNames (optional) TRUE if check whether trees have the same leaves set and whether is unrooted. Defaults to FALSE.
     * @return The quartet distance between trees.
     * @throw bpp::Exception if trees have different leaves sets or any is rooted.
     */                 
    static int quartetDistance(const TreeTemplate<Node>& trIn1, const TreeTemplate<Node>& trIn2, bool setLeavesId = true, bool checkNames = false)
            throw (Exception);

    /**
     * @brief The Triplets distance between two rooted bifurcating trees with the same set of leaves.
     * \n Counts the number of subtrees of three taxa that are different in the 
     * trees. The distance is called in [18] a close cousin of the quartet metric for
     * unrooted trees.
     * 
     * \n The constraints for the two input trees: Either the same leaves id sets numbered 0..n-1 and internal nodes ids numbered n, n+1,...  or the same leaves name sets and setLeavesId parameter true
     * \n\n Time complexity: O(n^2)
     * \n Algorithm adapted from D.E. Critchlow, D.K. Pearl, and C. Qian. The triplets distance for rooted bifurcating phylogenetic trees. Systematic biology, 45(3):323, 1996.
     *
     * @param[in] tr1 First rooted bifurcating tree.
     * @param[in] tr2 Second rooted bifurcating tree.
     * @param[in] setLeavesId (optional) TRUE if the two trees do not have the same ids for the same leaves or the ids are not numbered 0..n-1. 
     * @param[in] checkNames (optional) TRUE if check whether trees have the same leaves set and whether is binary, rooted. Defaults to FALSE.
     * @return The triplets distance between trees.
     * @throw bpp::Exception if trees have different leaves sets or any is unrooted.
     */
    static int tripletsDistance(const TreeTemplate<Node> & trIn1, const TreeTemplate<Node> & trIn2, bool setLeavesId = true, bool checkNames = false)
            throw (Exception);


    /**
     * @brief The Nodal distance between two unrooted trees with the same set of leaves.
     * It compares the phylogenetic trees considering the branch distances between each pair of nodes. 
     * To count the distance there is a manhattan metric used.
     * \n\n Time complexity: O(n^2)
     *
     * The Nodal metric bases on comparing how much of the leafA-to-leafB distances from tree1 repeat in tree2.
     * The two trees must share a common set of leaves (checked if checkNames is true).
     * \n Algorithm adapted from J. Bluis and D.G. Shin. Nodal distance algorithm: Calculating a phylogenetic tree comparison metric. In Bioinformatics and Bioengineering, 2003. Proceedings. Third IEEE Symposium on, pages 87–94. IEEE, 2003.
     * 
     * @param[in] tr1    The first unrooted tree.
     * @param[in] tr2    The second unrooted tree.
     * @param[in] metric Type of metric (mathematical expression) to use.
     * @param[in] setLeavesId       (optional) TRUE if the two trees do not have the same ids for the same leaves or the ids are not numbered 0..n-1. 
     * @param[in] checkNames        (optional) TRUE if check whether trees have the same leaves set and whether unrooted. Defaults to FALSE.
     * @return The Nodal distance between trees.
     * @throw DifferentLeavesSetsException if trees have different leaves sets or any is rooted.
     */
    static int nodalDistance(const TreeTemplate<Node>& tr1, const TreeTemplate<Node>& tr2, bool setLeavesId = true, bool checkNames = false)
            throw (Exception);

    /**
     * @brief The Nodal distance between two unrooted trees with the same set of leaves.
     * It compares the phylogenetic trees considering the branch distances between each pair of nodes. 
     * To count the distance there is a pythagorean metric used.
     * \n\n Time complexity: O(n^2)
     *
     * The Nodal metric bases on comparing how much of the leafA-to-leafB distances from tree1 repeat in tree2.
     * The two trees must share a common set of leaves (checked if checkNames is true).
     * \n Algorithm adapted from J. Bluis and D.G. Shin. Nodal distance algorithm: Calculating a phylogenetic tree comparison metric. In Bioinformatics and Bioengineering, 2003. Proceedings. Third IEEE Symposium on, pages 87–94. IEEE, 2003.
     * 
     * @param[in] tr1    The first unrooted tree.
     * @param[in] tr2    The second unrooted tree.
     * @param[in] metric Type of metric (mathematical expression) to use.
     * @param[in] setLeavesId       (optional) TRUE if the two trees do not have the same ids for the same leaves or the ids are not numbered 0..n-1. 
     * @param[in] checkNames        (optional) TRUE if check whether trees have the same leaves set and whether unrooted. Defaults to FALSE.
     * @return The Nodal distance between trees.
     * @throw DifferentLeavesSetsException if trees have different leaves sets or any is rooted.
     */	
    static double nodalDistance_pythagorean(const TreeTemplate<Node>& tr1, const TreeTemplate<Node>& tr2, bool setLeavesId = true, bool checkNames = false)
            throw (Exception);    

    /**
     * @brief The Nodal distance between two unrooted weighted trees with the same set of leaves.
     * It compares the phylogenetic trees considering the branch distances between each pair of nodes and the weights on the branches. 
     * To count the distance there is a manhattan metric used.
     * \n\n Time complexity: O(n^2)
     *
     * For each pair of leaves in each tree there is computed a sum of weights on branches that build a path between the pair of leaves. 
     * Next for each pair of leaves their distance in trees T1 and T2 are compared with use of manhattan metric. 
     * The two trees must share a common set of leaves (checked if checkNames is true).
     * \n Algorithm adapted from J. Bluis and D.G. Shin. Nodal distance algorithm: Calculating a phylogenetic tree comparison metric. In Bioinformatics and Bioengineering, 2003. Proceedings. Third IEEE Symposium on, pages 87–94. IEEE, 2003.
     * 
     * @param[in] tr1    The first unrooted weighted tree.
     * @param[in] tr2    The second unrooted weighted tree.
     * @param[in] metric Type of metric (mathematical expression) to use.
     * @param[in] setLeavesId       (optional) TRUE if the two trees do not have the same ids for the same leaves or the ids are not numbered 0..n-1. 
     * @param[in] checkNames        (optional) TRUE if check whether trees have the same leaves set and whether unrooted. Defaults to FALSE.
     * @return The Nodal distance between trees.
     * @throw DifferentLeavesSetsException if trees have different leaves sets or any is rooted.
     */
    static double nodalDistanceW(const TreeTemplate<Node>& tr1, const TreeTemplate<Node>& tr2, bool setLeavesId = true, bool checkNames = false)
            throw (Exception);

    /**
     * @brief The Nodal distance between two unrooted weighted trees with the same set of leaves.
     * It compares the phylogenetic trees considering the branch distances between each pair of nodes and the weights on the branches. 
     * To count the distance there is a pythagorean metric used.
     * \n\n Time complexity: O(n^2)
     *
     * For each pair of leaves in each tree there is computed a sum of weights on branches that build a path between the pair of leaves. 
     * Next for each pair of leaves their distance in trees T1 and T2 are compared with use of pythagorean metric. 
     * The two trees must share a common set of leaves (checked if checkNames is true).
     * \n Algorithm adapted from J. Bluis and D.G. Shin. Nodal distance algorithm: Calculating a phylogenetic tree comparison metric. In Bioinformatics and Bioengineering, 2003. Proceedings. Third IEEE Symposium on, pages 87–94. IEEE, 2003.
     * 
     * @param[in] tr1    The first unrooted weighted tree.
     * @param[in] tr2    The second unrooted weighted tree.
     * @param[in] metric Type of metric (mathematical expression) to use.
     * @param[in] setLeavesId       (optional) TRUE if the two trees do not have the same ids for the same leaves or the ids are not numbered 0..n-1. 
     * @param[in] checkNames        (optional) TRUE if check whether trees have the same leaves set and whether unrooted. Defaults to FALSE.
     * @return The Nodal distance between trees.
     * @throw DifferentLeavesSetsException if trees have different leaves sets or any is rooted.
     */
    static double nodalDistanceW_pythagorean(const TreeTemplate<Node>& tr1, const TreeTemplate<Node>& tr2, bool setLeavesId = true, bool checkNames = false)
            throw (Exception);



private:
    static bool checkLeavesNames(const TreeTemplate<Node>& trIn1, const TreeTemplate<Node>& trIn2)
            throw (bpp::Exception);

    static bool checkRooted(bool condition, const TreeTemplate<Node>& trIn1, const TreeTemplate<Node>& trIn2)
            throw (bpp::Exception);

    static int getPMDistance(ITwoTreesDescriptionElements& descriptionElements);

    static void getTreeNodesDists(TreeTemplate<Node>& tr, vector<vector <int> >& trNDists);

    static double getNodalDistance(INodesDist *d, const TreeTemplate<Node>& tr1, const TreeTemplate<Node>& tr2, bool setLeavesId = true, bool checkNames = false)
            throw (Exception);    

};
} // end of namespace
#endif	/* PHYLOTREEDIST_H */

