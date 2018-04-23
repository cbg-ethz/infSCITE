/*
 * doublet.h
 *
 *  Created on: May 25, 2016
 *      Author: jahnka
 */

#ifndef DOUBLETS_H_
#define DOUBLETS_H_

#include <vector>
#include <math.h>
#include <float.h>

using namespace std;

class Beta_Distr {

public:
    double mean, sd, alpha, beta;
    bool uniform;

public:
    Beta_Distr (double priorMean, double priorSD, bool uni){
    	mean = priorMean;
    	sd = priorSD;
    	uniform = uni;
    	if(uniform){
    		alpha = 1.0;                                    // setting alpha and beta of the beta distribution for alpha to 1 gives a uniform prior
    		beta = 1.0;
    	}
    	else{
    		alpha = ((1-mean)*mean*mean/(sd*sd)) - mean;    // <-10.13585344 turn the mean and sd into parameters of the beta distribution
    	  	beta = alpha*((1/mean)-1);                      //<-13.38666556
    	}
    };
    double logBetaPDF (double x) {return lgamma(alpha+beta)+(alpha-1)*log(x)+(beta-1)*log(1-x)-lgamma(alpha)-lgamma(beta);}
};

#define DOUBLET_UNIFORM_BETA_DISTR_beta  false  // if true, use uniform distribution for the prior of beta, sets bpriora_beta=1 and bpriorb_beta=1
#define DOUBLET_UNIFORM_BETA_DISTR_alpha false  // if true, use uniform distribution for the prior of alpha, sets bpriora_alpha=1 and bpriorb_alpha=1


std::string runMCMCnew(vector<struct treeTheta>& bestTrees, double* errorRates, int noOfReps, int noOfLoops, double gamma, std::vector<double> moveProbs, int n, int m, int** dataMatrix, char scoreType, int* trueParentVec, int step, bool sample, double chi, double priorSd_beta, double priorSd_alpha, bool useTreeList, char treeType, bool doubletModel, double doubletRate, int recMut, int z0, int z1, double alphaMoveRatio, bool fixedDoubletProb);
double scoreTreeFastWithDoublets(int*parent, int n, int m, double** logScores, int** dataMatrix, bool** ancMatrix, double& propDoubletProb, double& propRelDoubletProb, int recMut, bool fixedDoubletProb);
double getCombinedProb(double L_0, double L_1, double p);
double getDoubletProb(double L_0, double L_1, double p);
vector<vector<double> > getDoubletAttachmentScoreMatrixFast(double* singletAttachmentScore, int n, int* bft, int* parent, double** logScores, int* dataVector);
vector<vector<double> > getRecMutDoubletAttachmentScoreMatrixFast(double* singletAttachmentScore, int n, int* bft, bool** anc, int* parent, double** logScores, int* obs, int recMut);
double getSampleDoubletScoreFast(vector<vector<double> > doubletAttachmentScore, int n, bool** ancMatrix, double& relDoubletScore);
double getSampleSingletScoreFast(double* singletScores, int n, int*bft);
double sumScoreTreeFastSinglet(int n, int m, double** logScores, int** dataMatrix, int* parent);
double optimizeDoubletProb(vector<double> L_0, vector<double> L_1, int m);
std::string sampleFromPosterior_doublets(double currTreeLogScore, int n, int* currTreeParentVec, double thetaProb, double currBeta, double currScore, double currAlpha, double currRelDoubletRate, double currDoubletRate);
bool isInvalidRecMutTree(int* parent, int recMut,int copy);
bool isLinearTree(int* parent, int* bft, int n);


#endif /* DOUBLETS_H_ */
