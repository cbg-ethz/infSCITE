/*
 * mcmc.h
 *
 *  Created on: Mar 27, 2015
 *      Author: jahnka
 */

#ifndef MCMC_H
#define MCMC_H


double logBetaPDF(double x, double bpriora, double bpriorb);
double proposeNewBeta(double currBeta, double jumpSd);
double proposeNewAlpha(double currAlpha, double jumpSd);
double sampleNormal(double mean, double sd);
string sampleFromPosterior(double currTreeLogScore, int n, int* currTreeParentVec, double betaProb, double currBeta, double currScore);
int updateMinDistToTrueTree(int* trueParentVec, int* currTreeParentVec, int length, int minDistToTrueTree, int currScore, int bestScore);
int getSimpleDistance(int* trueVector, int* predVector, int length);


#endif
