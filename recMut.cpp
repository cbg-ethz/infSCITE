/*
 * doubleMut.cpp
 *
 *  Created on: Aug 15, 2015
 *      Author: jahnka
 */



#include <stdbool.h>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <iostream>
#include <sstream>
#include <math.h>
#include <stdbool.h>
#include <vector>
#include <stdlib.h>
#include <float.h>
#include "matrices.h"
#include "trees.h"
#include "scoreTree.h"
#include "recMut.h"

using namespace std;
double epsilon_recMut = 0.000000000001;  // maximal distance to current best tree where still using the accurate score computation


/* Computes the score of a new candidate tree with one recurrent mutation. First a fast approximate score is computed, then if the new score   */
/*  is better, or slightly worse than the best score so far, a more accurate but more costly score computation is done.                    */
double scoreTree_recMut(int n, int m, double** logScores, int** dataMatrix, char type, int* parentVector, double bestTreeLogScore, int recMut){

	double approx = scoreTreeFast_recMut(n, m, logScores, dataMatrix, type, parentVector, recMut);   // approximate score

	if(approx > bestTreeLogScore-epsilon_recMut){                                                  // approximate score is close to or better
		return scoreTreeAccurate_recMut(n, m, logScores, dataMatrix, type, parentVector, recMut);   // than the current best score, use accurate
	}                                                                                // score computation

	return approx;                                                              // otherwise the approximate score is sufficient
}

/* computes an approximate score for a tree with one recurrent mutation. This is fast, but rounding errors can occur  */
double scoreTreeFast_recMut(int n, int m, double** logScores, int** dataMatrix, char type, int* parentVector, int recMut){

	double result = -DBL_MAX;
	int treeSize = n+2;
	int* bft = getBreadthFirstTraversal(parentVector, treeSize-1);        // get breadth first traversal for scoring by updating the parent score
	int backMut = getBackmutation(parentVector, recMut, n);   // which copy of the recurrent mutation is the back mutation, if none is, use -1

	if(type=='m'){
		result = maxScoreTreeFast_recMut(n, m, logScores, dataMatrix, parentVector, bft, recMut, backMut);  // score by best attachment point per sample
	}
	else if(type=='s'){
		result = sumScoreTreeFast_recMut(n, m, logScores, dataMatrix, parentVector, bft, recMut, backMut);  // score by summing over all attachment points
	}

	delete [] bft;
	return result;
}

/* computes the accurate score for a tree with one recurrent mutation. This is slower that the fast computation but fewer rounding errors occur  */
double scoreTreeAccurate_recMut(int n, int m, double** logScores, int** dataMatrix, char type, int* parentVector, int recMut){

	//print_intArray(parentVector, n+1);
	double result = -DBL_MAX;
	int* bft = getBreadthFirstTraversal(parentVector, n+1);          // tree with n+2 nodes, n mutations + one recurrent mutation + the root node
	int backMut = getBackmutation(parentVector, recMut, n);   // which copy of the recurrent mutation is the back mutation, if none is, use -1

	if(type=='m'){
		result = maxScoreTreeAccurate_recMut(n, m, logScores, dataMatrix, parentVector, bft, recMut, backMut);
	}
	else if(type=='s'){
		result = sumScoreTreeAccurate_recMut(n, m, logScores, dataMatrix, parentVector, bft, recMut, backMut);
	}

	delete [] bft;
	return result;
}

/* computes an approximate scoring for a tree using the max attachment score per sample */
double maxScoreTreeFast_recMut(int n, int m, double** logScores, int** dataMatrix, int* parent, int* bft, int recMut, int backMut){

    double treeScore = 0.0;

  	for(int sample=0; sample<m; sample++){                                                                                      // get for every sample
  		double* scores = getAttachmentScoresFast_recMut(parent,n, logScores, dataMatrix[sample], bft, recMut, backMut);   // all attachment scores
  		treeScore +=  getMaxEntry(scores, n+1);                                                                                 // pick the best score per sample and sum them up
  		delete [] scores;
  	}

  	return treeScore;    // sum over the best attachment scores of all samples is tree score
}

/* computes an approximate scoring for a tree with one recurrent mutation summing the score over all attachment points per sample */
double sumScoreTreeFast_recMut(int n, int m, double** logScores, int** dataMatrix, int* parent, int* bft, int recMut, int backMut){

	double sumTreeScore = 0.0;

	for(int sample=0; sample<m; sample++){
		double* scores = getAttachmentScoresFast_recMut(parent, n, logScores, dataMatrix[sample], bft, recMut, backMut); // attachments scores of sample to each node
		double bestMaxTreeScore = getMaxEntry(scores, n+1);                                     // best of the scores (used to compute with score differences rather than scores)

		double sumScore = 0.0;
		for(int i=0; i<=n+1; i++){                                                 // sum over all attachment scores, exp is necessary as scores are actually log scores
			sumScore += exp(scores[bft[i]]-bestMaxTreeScore);                   // subtraction of best score to calculate with score differences (smaller values)
		}
		delete [] scores;
		sumTreeScore += log(sumScore)+bestMaxTreeScore;                     // transform back to log scores and change from score differences to actual scores
	}
	return sumTreeScore;
}

/* computes the log score for the complete tree using the maxScore scheme, where the best attachment point is used for each sample */
/* This uses the accurate scoring scheme to minimize rounding errors*/
double maxScoreTreeAccurate_recMut(int n, int m, double** logScores, int** dataMatrix, int* parent, int* bft, int recMut, int backMut){

    int** treeScoreMatrix = init_intMatrix(4, 2, 0);   // to track how many counts of each score type alpha, beta, (1-alpha), etc. we have in the complete tree scoring

  	for(int sample=0; sample<m; sample++){
  		int** bestAttachmentMatrix =  getBestAttachmentScoreAccurate_recMut(parent, n, logScores, dataMatrix[sample], bft, recMut, backMut); // get matrix for attachment of single sample to tree
  		treeScoreMatrix = sumMatrices(treeScoreMatrix, bestAttachmentMatrix, 4, 2);
  		free_intMatrix(bestAttachmentMatrix);
  	}
  	double treeScore = getTrueScore(treeScoreMatrix, logScores);
  	free_intMatrix(treeScoreMatrix);
  	return treeScore;
}

/* computes the log score for the complete tree using the sumScore scheme, where likelihoods of all attachment points of a sample are added */
/* This uses the accurate scoring scheme to minimize rounding errors*/
double sumScoreTreeAccurate_recMut(int n, int m, double** logScores, int** dataMatrix, int* parent, int* bft, int recMut, int backMut){

	double sumTreeScore = 0.0;

	for(int sample=0; sample<m; sample++){
		double score = getSumAttachmentScoreAccurate_recMut(parent, n, logScores, dataMatrix[sample], bft, recMut, backMut);
		//cout << score << "\t";
		sumTreeScore += score;
	}
	//cout << "\n";
	return sumTreeScore;
}


/* computes the sum score for attaching a sample to all nodes in a tree with one recurrent mutation */
double getSumAttachmentScoreAccurate_recMut(int* parent, int n, double** logScores, int* obsMutState, int* bft, int recMut, int backMut){

	int treeSize = n+2;
	int*** attachmentScoreMatrix;                                   // matrix to keep attachment scores for each sample (not summing up components to avoid rounding errors)

	attachmentScoreMatrix = getAttachmentMatrices_recMut(parent, n, obsMutState, bft, recMut, backMut);
	double* attachmentScore = getTrueScores(attachmentScoreMatrix, treeSize-1, logScores);                   // get the true attachment scores from the attachment matrices
	double bestScore = getMaxEntry(attachmentScore, n+1);                                            // identify best attachment score
	double sumScore = 0.0;
	for(int parent = 0; parent<treeSize; parent++){                                                        // get score for attaching to the other nodes in the tree
		sumScore += exp(attachmentScore[parent]-bestScore);
	}
	delete_3D_intMatrix(attachmentScoreMatrix, n+2);
	delete [] attachmentScore;
	return log(sumScore)+bestScore;
}

/* This function gets a tree with exactly one recurrent mutation and tests if it is a back mutation. */
/* It returns the node id of the back mutation (either recMut or n) if it exists, and -1 if elsewise */
int getBackmutation(int* parent, int recMut, int n){

	int copy = n;        // id of the second copy of recurrent mutation, this is always n, as the first n nodes run from 0 to n-1
	int root = n+1;      // id of the root, this is here shifted to n+1, as n is now the second copy of the recurrent mutation
	int backMut = -1;   // the id of the back mutation if one exists; takes [values -1, recMut, or n (the id of the copy)].

	int node = recMut;                // test if the first copy of the recurrent mutation is a back mutation
	while(node!=root){
		if(parent[node]==copy){
			backMut = recMut;
			break;
		}
		node = parent[node];
	}

	node = copy;                          // test if the second copy of the recurrent mutation is a back mutation
	while(node!=root){
		if(parent[node]==recMut){
			backMut = copy;
			break;
		}
		node = parent[node];
	}
	return backMut;        // return id of back mutation if one exists, otherwise return -1 (in case of parallel mutation)
}

/* computes the attachment scores of a sample to all nodes in the tree */
double* getAttachmentScoresFast_recMut(int*parent, int n, double** logScores, int* obsMutProfile, int*bft, int recMut, int backMut){

	int root = n+1;                                                        // the id of the root, here shifted to n+1, as n is now the second copy of the recurrent mutation
	int treeSize = n+2;                                                      // number of nodes in tree: n mutations + one recurrent mutation + root
	double* attachmentScore = init_doubleArray(treeSize, -DBL_MAX);
	attachmentScore[root] = rootAttachementScore(n, logScores, obsMutProfile);
	for(int i=1; i<treeSize; i++){                                                              // try all attachment points (nodes in the mutation tree)
		int node = bft[i];
		int obsMutState;
		if(node==n){
			obsMutState = obsMutProfile[recMut];
		}
		else{
			obsMutState = obsMutProfile[node];
		}
		attachmentScore[node] = attachmentScore[parent[node]];               // start from the score for attaching to the parent of the node
		if(node==backMut){                                                   // in case the node is the back mutation
			attachmentScore[node] -= logScores[obsMutState][1];         // then the expected mutation state of the recurrent mutation changes from 1 to 0
			attachmentScore[node] += logScores[obsMutState][0];
		}
		else{
			attachmentScore[node] -= logScores[obsMutState][0];         // otherwise the expected mutation state changes from 0 to 1
			attachmentScore[node] += logScores[obsMutState][1];
		}
	}
	return attachmentScore;
}

/* computes the best attachment score for the given sample to the given tree (with one recurrent mutation) using the accurate approach: */
/* this means that the score here is a matrix with counts of the different log score types to minimize rounding errors      */
int** getBestAttachmentScoreAccurate_recMut(int* parent, int n, double** logScores, int* obsMutState, int* bft, int recMut, int backMut){

	int treeSize = n+2;
	int*** attachmentScoreMatrix;       // list of matrices to keep best attachment scores for each sample (score is not summed up but matrix counting different components)

	attachmentScoreMatrix = getAttachmentMatrices_recMut(parent, n, obsMutState, bft, recMut, backMut);

	double bestScore =  -DBL_MAX;
	int bestScoreIndex = -1;

	for(int i=0; i<treeSize; i++){                                                   // now get true scores for each attachment point and find the best of them
		double newScore = getTrueScore(attachmentScoreMatrix[i], logScores);
		if(bestScore <= newScore){
			bestScoreIndex = i;
			bestScore = newScore;
		}
	}
	int** bestScoreMatrix = NULL;
	bestScoreMatrix = deepCopy_intMatrix(attachmentScoreMatrix[bestScoreIndex], 4, 2);
	delete_3D_intMatrix(attachmentScoreMatrix, treeSize);
	return bestScoreMatrix;
}


/* This computes the scores for attaching the given sample to every node in the given tree with exactly one recurrent */
/* mutation. This is the accurate scoring where each score is a matrix with counts of the different types of log scores */
/* P("observed state"|"expected state according to tree") */
int*** getAttachmentMatrices_recMut(int* parent, int n, int* obsMutState, int* bft, int recMut, int backMut){

	int root = n+1;     // the id of the root, this is here shifted to n+1, as n is now the second copy of the recurrent mutation
	int treeSize = n+2;
	int*** attachmentScoreMatrix = new int**[treeSize];      // list of matrices to keep scores for attaching the sample to the different nodes (not summing up components to avoid rounding errors)

	attachmentScoreMatrix[root] = init_intMatrix(4, 2, 0);      // start with attaching sample to root, i.e no mutation is present
	for(int mut=0; mut<n; mut++){
		attachmentScoreMatrix[root][obsMutState[mut]][0]++;
	}

	for(int i=1; i<treeSize; i++){               // get scores for attaching to the other nodes; due to bft this happens in an order such that the parent node has already been processed
		int node = bft[i];
		attachmentScoreMatrix[node] = deepCopy_intMatrix(attachmentScoreMatrix[parent[node]], 4, 2);   // start from the matrix of the parent of the node

		if(node==backMut){                                              // in case the node is the back mutation
			attachmentScoreMatrix[node][obsMutState[recMut]][1]--;   // then the expected mutation status of the recurrent mutation changes from 1 to 0
			attachmentScoreMatrix[node][obsMutState[recMut]][0]++;
		}
		else{
			int mut = node;
			if(node==n){ mut = recMut;}               // the second copy of the recurrent mutation has the same observed mutation state as the first copy
			//cout << "node = " << node << "\n";
			//cout << "obsState[mut] = " << obsMutState[mut] << "\n";
			attachmentScoreMatrix[node][obsMutState[mut]][0]--;       // otherwise the expected mutation status changes from 0 to 1
			attachmentScoreMatrix[node][obsMutState[mut]][1]++;
		}
	}
	return attachmentScoreMatrix;
}

