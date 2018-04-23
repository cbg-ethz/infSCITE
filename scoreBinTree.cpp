/*
 * scoreBinTree.cpp
 *
 *  Created on: Mar 14, 2016
 *      Author: jahnka
 */

#include <stdbool.h>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <iostream>
#include "limits.h"
#include "scoreBinTree.h"
#include "matrices.h"
#include "trees.h"

using namespace std;

/* computes the log likelihood for a single mutation for all subtrees of the binary tree, where the expected */
/* state of the mutation can be either absent or present in the whole subtree (passed as 'state' to the function) */
double* getBinSubtreeScore(bool state, int* bft, vector<vector<int> > &childLists, int mut, int nodeCount, int m, int** obsMutProfiles, double ** logScores){
	double* score = init_doubleArray(nodeCount, 0.0);
	for(int i=nodeCount-1; i>=0; i--){
		int node = bft[i];

		if(node < m){
			score[node] = logScores[obsMutProfiles[node][mut]][state];   // for leafs the score is just P(Dij|Eij)
		}
		else{                                                          // for inner nodes the score is the sum of the scores of the children
			if(childLists.at(node).size()!=2){
				cerr << "Error: node " << node << " has " << childLists.at(node).size() << " children\n";  // tree should be binary, but isn't
			}
			score[node] = score[childLists.at(node).at(0)] + score[childLists.at(node).at(1)];
		}
	}
	return score;
}


/* Computes the best log likelihood for placing a single mutation in a given sample tree */
/* Iterates through all nodes as possible placements of the mutation to find the best one */
/* All samples below the placement of the mutation should have it, mutation can also be placed at a leaf, i.e. uniquely to the sample   */
double getBinTreeMutScore(int* bft, vector<vector<int> > &childLists, int mut, int nodeCount, int m, int** obsMutProfiles, double ** logScores){

	double bestScore = -DBL_MAX;
	double* absentScore = getBinSubtreeScore(0, bft, childLists, mut, nodeCount, m, obsMutProfiles, logScores);
	double* presentScore = getBinSubtreeScore(1, bft, childLists, mut, nodeCount, m, obsMutProfiles, logScores);

	for(int p=0; p<nodeCount; p++){
		double score = absentScore[nodeCount-1] - absentScore[p] + presentScore[p];
		bestScore = max(bestScore, score);
	}

	delete [] absentScore;
	delete [] presentScore;
	return bestScore;
}

/* Computes the maximum log likelihood of a binary tree for a given mutation matrix.  */
/* Note: No extra root necessary for binary trees */
double getBinTreeScore(int** obsMutProfiles, int n, int m, double ** logScores, int* parent){

	int nodeCount = (2*m)-1;   // number of nodes in binary tree: m leafs, m-1 inner nodes (includes already the root)
	double sumScore = 0;       // sum of maximal scores of all samples
	vector<vector<int> > childLists = getChildListFromParentVector(parent, nodeCount-1);
	int* bft = getBreadthFirstTraversal(parent, nodeCount-1);

	for(int mut=0; mut<n; mut++){                                                // sum over the optimal scores of each sample
		double score = getBinTreeMutScore(bft, childLists, mut, nodeCount, m, obsMutProfiles, logScores);
		sumScore += score;
	}

	delete [] bft;
	for(int i=0; i<childLists.size(); i++){
		childLists[i].clear();
	}
	childLists.clear();
	//cout << "score: " << sumScore << "\n";
	return sumScore;
}


