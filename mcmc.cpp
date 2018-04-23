/*
 * mcmc.cpp
 *
 *  Created on: Mar 27, 2015
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
#include <random>
#include <fstream>
#include <sstream>
#include "matrices.h"
#include "treelist.h"
#include "trees.h"
#include "mcmc.h"
#include "scoreTree.h"
#include "scoreBinTree.h"
#include "rand.h"
#include "limits.h"
#include "output.h"
#include "mcmcBinTreeMove.h"
#include "mcmcTreeMove.h"
using namespace std;


double logBetaPDF(double x, double bpriora, double bpriorb){
	double logScore = lgamma(bpriora+bpriorb)+(bpriora-1)*log(x)+(bpriorb-1)*log(1-x)-lgamma(bpriora)-lgamma(bpriorb);    // f(x,a,b) = gamma(a+b)/(gamma(a)gamma(b)) * x^(a-1) * (1-x)^(b-1)
	//cout << lgamma(bpriora+bpriorb) << "   " << tgamma(bpriora+bpriorb) << "  #1\n";
	//cout << (bpriora-1)*log(x) << "  #2\n";
	//cout << (bpriorb-1)*log(1-x) << "  #3\n";
	//cout << log(tgamma(bpriora)) << "  #4\n";
	//cout << log(tgamma(bpriorb)) << "  #5\n";
	//cout << logScore << "\n";
	return logScore;
}

/* a new value for the error probability beta is sampled from a normal distribution around the current beta */
double proposeNewBeta(double currBeta, double jumpSd){
	double sampledValue = sampleNormal(0, jumpSd);
	double propBeta = currBeta+sampledValue ;                   //rnorm(1,0,jumpsd)
	if(propBeta < 0){
		propBeta = abs(propBeta);
	}
	if(propBeta > 1){
		propBeta = propBeta - 2*(propBeta-1);
	}
    return propBeta;
}

double proposeNewAlpha(double currAlpha, double jumpSd){
     double sampledValue = sampleNormal(0, jumpSd);
     double propAlpha = currAlpha;
     if(currAlpha<jumpSd){
       propAlpha = currAlpha*exp(sampledValue) ;    // for small alphas change log(alpha) instead of alpha
     } else {
       propAlpha = currAlpha+sampledValue ; //rnorm(1,0,jumpsd)
     }
     if(propAlpha < 0){
         propAlpha = abs(propAlpha);
     }
     if(propAlpha > 1){
         propAlpha = propAlpha - 2*(propAlpha-1);
     }
     return propAlpha;
}



/* samples a new value for beta from a normal distribution around the current value */
double sampleNormal(double mean, double sd) {
    double u = ((double) rand() / (RAND_MAX)) * 2 - 1;
    double v = ((double) rand() / (RAND_MAX)) * 2 - 1;
    double r = u * u + v * v;
    if (r == 0 || r > 1){
    	return sampleNormal(mean, sd);
    }
    double c = sqrt(-2 * log(r) / r);
    double value =  u * c;                       // value times sd and add the mean
    return (value * sd + mean);
}


/* prints out the current tree and beta to sample from the posterior distribution */
string sampleFromPosterior(double currTreeLogScore, int n, int* currTreeParentVec, double betaProb, double currBeta, double currScore){

	std::stringstream content;
	content << currTreeLogScore  << "\t";                 // logscore of current tree
	content << countBranches(currTreeParentVec, n);       // number of branches in current tree
	if(betaProb>0.0){
		content << "\t" << currBeta;                      // current beta
		content << "\t" << currScore;                     // current combined logscore for tree and beta
	}
	content << "\t";
	for(int i=0; i<n; i++){
		content << currTreeParentVec[i] << " ";
	}
	content << "\n";
	return content.str();
}


/* updates the minimum distance between any of the optimal trees and the true tree (if available) */
int updateMinDistToTrueTree(int* trueParentVec, int* currTreeParentVec, int length, int minDistToTrueTree, int currScore, int bestScore){

	int currDistToTrueTree = getSimpleDistance(trueParentVec, currTreeParentVec, length);

	if(currScore > bestScore){
		return currDistToTrueTree;
	}

	if(currScore == bestScore && currDistToTrueTree < minDistToTrueTree){         // the current tree is closest to the true tree among the current optimal trees
		return currDistToTrueTree;
	}

	return minDistToTrueTree;
}



/* Returns the distance between two trees, where the distance is the number of nodes having different parents in the two trees  */
int getSimpleDistance(int* trueVector, int* predVector, int length){
	int dist = 0;
	for(int i=0; i<length; i++){
		if(trueVector[i]!=predVector[i]){
			dist++;
		}
	}
	return dist;
}
