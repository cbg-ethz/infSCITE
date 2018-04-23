/*
 * mcmcDoublet.cpp
 *
 *  Created on: May 18, 2016
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
#include "doublets.h"
#include "recMut.h"

using namespace std;


double bestAvgSingletScoreSum = -DBL_MAX;
double bestAvgRelevantDoubletScoreSum = -DBL_MAX;
double propAvgSingletScoreSum = -DBL_MAX;
double propAvgRelevantDoubletScoreSum = -DBL_MAX;
double currAvgSingletScoreSum = -DBL_MAX;
double currAvgRelevantDoubletScoreSum = -DBL_MAX;

double doubletRatePrecision = 1e-6;    // defines to what precision the doublet rate should be optimized


std::string runMCMCnew(vector<struct treeTheta>& bestTrees, double* errorRates, int noOfReps, int noOfLoops, double gamma_, vector<double> moveProbs, int n, int m, int** dataMatrix, char scoreType, int* trueParentVec, int step, bool sample, double chi, double priorSd_beta, double priorSd_alpha, bool useTreeList, char treeType, bool doubletModel, double doubletRate, int recMut, int z0, int z1, double alphaMoveRatio, bool fixedDoubletProb){

	double burnInPhase = 0.25;                   // first quarter of steps are burn in phase
	unsigned int optStatesAfterBurnIn = 0;
	int burnIn = noOfLoops*burnInPhase;
	stringstream sampleOutput;
	unsigned long long int firstOccBestTree = 0; //number of steps at which the best tree was found the first time

	Beta_Distr alpha (errorRates[0], priorSd_alpha, DOUBLET_UNIFORM_BETA_DISTR_alpha);                  // prior distribution alpha (FD)
	Beta_Distr beta (errorRates[1] + errorRates[2], priorSd_beta, DOUBLET_UNIFORM_BETA_DISTR_beta);     // prior distribution beta (AD1 + AD2)

	int parentVectorSize = n;
	if(recMut >= 0){ parentVectorSize++; }             // one more node for recurrent mutation
	if(treeType=='t'){parentVectorSize = (2*m)-2;}    // transposed case: binary tree, m leafs and m-1 inner nodes, root has no parent

	double jumpSd_beta = beta.sd/chi;             // chi: scaling of the known error rate for the MH jump; resulting jump sd
	double jumpSd_alpha = alpha.sd/chi;           // chi: scaling of the known error rate for the MH jump; resulting jump sd

	int minDistToTrueTree = INT_MAX;             // smallest distance between an optimal tree and the true (if given)
	double bestTreeLogScore = -DBL_MAX;          // log score of T in best (T,beta)
	double bestScore = -DBL_MAX;                 // log score of best combination (T, beta)
	double bestBeta = beta.mean;
	double bestAlpha = alpha.mean;
	double bestRelDoubletRate = doubletRate;
	double bestDoubletRate = doubletRate;

	for(int r=0; r<noOfReps; r++){   // repeat the MCMC, start over with random tree each time, only best score and list of best trees is kept between repetitions

		cout << "MCMC repetition " << r << "\n";
		int*   currTreeParentVec;
		if(treeType=='m'){
			if(recMut<0){
				currTreeParentVec = getRandParentVec(parentVectorSize);       // start MCMC with random tree
			}
			else{
				bool invalidDoubleMutTree = true;
				while(invalidDoubleMutTree==true){
					currTreeParentVec = getRandParentVec(parentVectorSize);
					invalidDoubleMutTree = isInvalidRecMutTree(currTreeParentVec, recMut, n);
				}
			}
		}
		else{     currTreeParentVec = getRandomBinaryTree(m);}                     // transposed case: random binary tree

		bool** currTreeAncMatrix =  parentVector2ancMatrix(currTreeParentVec,parentVectorSize);
		double** currLogScores = getLogScores(errorRates[0], errorRates[1], errorRates[2], errorRates[3]);           // compute logScores of conditional probabilities
		double currBeta = beta.mean;                                                                                  // the current AD rate
		double currAlpha = alpha.mean;
		double currDoubletProb = doubletRate; // = doubletBetaPriorMean;
		double currTreeLogScore;
		double propDoubletProb = doubletRate; // = currDoubletProb;    // this is optimized for the proposed tree in every MCMC loop
		double currRelDoubletProb = 0.0;
		double propRelDoubletProb = 0.0;

		/* score tree */
		if(treeType=='m'){
			if(doubletModel){    // doublet model with and without recurrent mutation

				currTreeLogScore = scoreTreeFastWithDoublets(currTreeParentVec, n, m, currLogScores, dataMatrix, currTreeAncMatrix, propDoubletProb, propRelDoubletProb, recMut, fixedDoubletProb);
			}
			else{
				// no doublets, no recurrent mutation
				if(recMut==-1){
					if(scoreType=='s'){
						currTreeLogScore = sumScoreTreeFastSinglet(n, m, currLogScores, dataMatrix, currTreeParentVec);
					}
					else{
						currTreeLogScore = scoreTree(n, m, currLogScores, dataMatrix, scoreType, currTreeParentVec, bestTreeLogScore);
					}
				}
				else{   // no doublets, recurrent mutation
					currTreeLogScore = scoreTree_recMut(n, m, currLogScores, dataMatrix, treeType, currTreeParentVec, bestTreeLogScore, recMut);
				}
			}
		}
		else{   // binary tree, does not work with doublets, and a recurrent mutation
			currTreeLogScore = getBinTreeScore(dataMatrix, n, m, currLogScores, currTreeParentVec);
		}
		double currBetaLogScore = (moveProbs[0]==0) ? 0.0 : beta.logBetaPDF(currBeta);            // zero if beta is fixed
		double currAlphaLogScore = (moveProbs[0]==0) ? 0.0 : alpha.logBetaPDF(currAlpha);          // zero if alpha is fixed
		double currThetaLogScore = currBetaLogScore + currAlphaLogScore + z0*log(1-currAlpha) + z1*log(currAlpha);
		double currScore = currTreeLogScore + currThetaLogScore;
		//cout << "start score: " << currScore << "\n";

		for(int it=0; it<noOfLoops; it++){                                     // run the iterations of the MCMC
			//cout << it << endl;
        	if(it % 10000 == 0){
        		cout.precision(16);
        		cout << "At mcmc repetition " << r+1 << "/" << noOfReps << ", step " << it << ": ";
        		cout << "best overall score " << bestScore << ", best tree score " << bestTreeLogScore;
        		cout << ", best beta " << bestBeta << ", best alpha " << bestAlpha;
        		cout << ", best rel. doublet rate " <<  bestRelDoubletRate << " and best doublet rate " << bestDoubletRate;
        		cout << endl;
        	}

        	bool moveAccepted = false;                               // Is the MCMC move accepted?
        	bool moveChangesTheta = changeBeta(moveProbs[0]);        // true if this move changes beta, not the tree

        	if(moveChangesTheta){                     // new theta is proposed, log scores change, tree is copy of current tree
        		double propBeta = proposeNewBeta(currBeta, jumpSd_beta);
        		double propAlpha = proposeNewAlpha(currAlpha, jumpSd_alpha);

        		if (sample_0_1() < alphaMoveRatio){
        			propBeta = currBeta;                 // alpha is changed
        		}
        		else {
        		    propAlpha = currAlpha;              // beta is changed
        		}

        		double** propLogScores = deepCopy_doubleMatrix(currLogScores, 4, 2);
        		updateLogScoresAlphaBeta(propLogScores, propBeta, propAlpha);
        		double propBetaLogScore = beta.logBetaPDF(propBeta);
				double propAlphaLogScore = alpha.logBetaPDF(propAlpha);
				double propThetaLogScore = propBetaLogScore + propAlphaLogScore + + z0*log(1-propAlpha) + z1*log(propAlpha);
        		double propTreeLogScore;

        		if(treeType=='m'){
        			if(doubletModel){
        				propTreeLogScore = scoreTreeFastWithDoublets(currTreeParentVec, n, m, propLogScores, dataMatrix, currTreeAncMatrix, propDoubletProb, propRelDoubletProb, recMut, fixedDoubletProb);   // compute the new tree score for new beta
        			}
        			else{
        				if(recMut==-1){
        					if(scoreType=='s'){
        						propTreeLogScore = sumScoreTreeFastSinglet(n, m, propLogScores, dataMatrix, currTreeParentVec);
        					}
        					else{
        						propTreeLogScore = scoreTree(n, m, propLogScores, dataMatrix, scoreType, currTreeParentVec, bestTreeLogScore);
        					}
        				}
        				else{
        					propTreeLogScore = scoreTree_recMut(n, m, propLogScores, dataMatrix, treeType, currTreeParentVec, bestTreeLogScore, recMut);

        				}
        			}
        			//cout << "proposed theta score diff: " << propThetaLogScore- currThetaLogScore << "      prop score: " << propThetaLogScore << "  curr score: " << currThetaLogScore << "\n";
        		}
        		else{
        			propTreeLogScore = getBinTreeScore(dataMatrix, n, m, propLogScores, currTreeParentVec);
        		}

        		if (sample_0_1() < exp((propTreeLogScore+propThetaLogScore-currTreeLogScore-currThetaLogScore)*gamma_)){               // the proposed move is accepted
        			moveAccepted = true;
        			free_doubleMatrix(currLogScores);
        		    currTreeLogScore  = propTreeLogScore;                                       // update score of current tree
        		    currBeta = propBeta;                                                        // the current AD rate
        		    currAlpha = propAlpha;
        		    currBetaLogScore = propBetaLogScore;
        		    currAlphaLogScore = propAlphaLogScore;
        		    currThetaLogScore = propThetaLogScore;
        		    currScore = currTreeLogScore+currThetaLogScore;                          // combined score of current tree and current beta
        		    currLogScores = propLogScores;
        		    currAvgSingletScoreSum = propAvgSingletScoreSum;
        		    currAvgRelevantDoubletScoreSum	= propAvgRelevantDoubletScoreSum;
        		    currDoubletProb = propDoubletProb;
        		}
        		else{
        			delete [] propLogScores[0];
        			delete [] propLogScores;
        		}
        	}
        	else{                                   // move changed tree
        		double nbhcorrection = 1.0;
        		int* propTreeParVec;
        		double propTreeLogScore;

        		if(treeType=='m'){

        			propTreeParVec = proposeNewTree(moveProbs, parentVectorSize, currTreeAncMatrix, currTreeParentVec, nbhcorrection);              // propose new tree and

        			if(recMut>=0){
        				while(isInvalidRecMutTree(propTreeParVec, recMut, n)){
        					delete [] propTreeParVec;
        					propTreeParVec = proposeNewTree(moveProbs, parentVectorSize, currTreeAncMatrix, currTreeParentVec, nbhcorrection);
        				}
        			}
        			bool** propTreeAncMatrix = parentVector2ancMatrix(propTreeParVec, parentVectorSize);

        			if(doubletModel){
        				propTreeLogScore = scoreTreeFastWithDoublets(propTreeParVec, n, m, currLogScores, dataMatrix, propTreeAncMatrix, propDoubletProb, propRelDoubletProb, recMut, fixedDoubletProb);

        			}
        			else{
        				if(recMut==-1){
        					if(scoreType=='s'){
        						propTreeLogScore = sumScoreTreeFastSinglet(n, m, currLogScores, dataMatrix, propTreeParVec);
        					}
        					else{
        						propTreeLogScore = scoreTree(n, m, currLogScores, dataMatrix, scoreType, propTreeParVec, bestTreeLogScore);
        					}
        				}
        				else{
        					propTreeLogScore = scoreTree_recMut(n, m, currLogScores, dataMatrix, treeType, propTreeParVec, bestTreeLogScore, recMut);
        				}
        			}
        			//cout << "proposed score: " << propTreeLogScore << " after tree move, before tree score: " << currTreeLogScore << "\n";
        			free_boolMatrix(propTreeAncMatrix);
        		}
        		else{
        			propTreeParVec = proposeNextBinTree(moveProbs, m, currTreeParentVec, currTreeAncMatrix);
        			propTreeLogScore = getBinTreeScore(dataMatrix, n, m, currLogScores, propTreeParVec);
        		}

        		//cout << propTreeLogScore << "\t" << currTreeLogScore << endl;
        		if (sample_0_1() < nbhcorrection*exp((propTreeLogScore-currTreeLogScore)*gamma_)){                    // the proposed tree is accepted
        			moveAccepted = true;
        			free_boolMatrix(currTreeAncMatrix);                                            // discard outdated tree
        			delete[] currTreeParentVec;
        			currTreeAncMatrix = parentVector2ancMatrix(propTreeParVec,parentVectorSize); // update matrix of current tree
        			currTreeParentVec = propTreeParVec;                                         // update parent vector of current tree
        			currTreeLogScore  = propTreeLogScore;                                       // update score of current tree
        			currScore = currTreeLogScore+currThetaLogScore;
        			currAvgSingletScoreSum = propAvgSingletScoreSum;
        			currAvgRelevantDoubletScoreSum	= propAvgRelevantDoubletScoreSum;
        			currDoubletProb = propDoubletProb;
        			currRelDoubletProb = propRelDoubletProb;

        			//cout  << "curr p = " << currDoubletProb << " current relevant doublet prob: " << currRelDoubletProb << "\n";
        		}
        		else{
        			delete [] propTreeParVec;            // discard proposed tree
        		}
        	}

        	/* If the true tree is given update the smallest distance between a currently best tree and the true tree */
        	if(trueParentVec){
        		minDistToTrueTree = updateMinDistToTrueTree(trueParentVec, currTreeParentVec, parentVectorSize, minDistToTrueTree, currScore, bestScore);
        	}

        	/* If the list of optimal trees is used, update it */
        	if(useTreeList){
        		updateTreeList(bestTrees, currTreeParentVec, parentVectorSize, currScore, bestScore, currBeta, currAlpha);
        	}

        	/* Sample from the posterior if required and past the burn-in phase */
        	if(sample && it>=burnIn && it % step == 0){
        		sampleOutput << sampleFromPosterior_doublets(currTreeLogScore, parentVectorSize, currTreeParentVec, moveProbs[0], currBeta, currScore, currAlpha, currRelDoubletProb, currDoubletProb);
        	}

        	/* Update best tree in case we have found a new best one */
        	if(currScore > bestScore){
        		optStatesAfterBurnIn = 0;                    // new opt state found, discard old count
        		bestTreeLogScore = currTreeLogScore;
        		bestScore = currScore;                 // log score of best combination (T, beta)
        		bestBeta = currBeta;
        		bestAlpha = currAlpha;
        		bestRelDoubletRate = currRelDoubletProb;
        		bestDoubletRate = currDoubletProb;

        		firstOccBestTree = r*noOfLoops + it;

        		//cout << firstOccBestTree << ": " << minDistToTrueTree << "  " << bestTreeLogScore <<endl;
        		//cout << firstOccBestTree << ": " << bestTreeLogScore << "\t" << currBeta << "\t" << currAlpha << endl;

        		//writeToFile(getGraphVizBinTree(bestTreeParentVec, m, label, bestTreeLogScore, leafClusterId), "/Users/jahnka/Desktop/AcetoData/CTC-Exome_vcf_180105/Br23.gv");



        		//cout << "new best score: " << bestTreeLogScore << endl;
        		//writeToFile(getGraphVizBinTree(currTreeParentVec, m, label, bestTreeLogScore, leafClusterId), "/Users/jahnka/Desktop/AcetoData/CTC-Exome_vcf_180105/Br23.gv");
        		//getMutTreeGraphVizString(currTreeParentVec, parentVectorSize, vector<string> nodeLabels, false, vector<string> sampleLabels, dataMatrix, currLogScores, false, -1);

        		//cout << "tree written\n";
        		//cout << "new best score: " << bestTreeLogScore << "   " << currDoubletProb << " current relevant doublet prob: " << currRelDoubletProb << "   all doublets rate: " << currDoubletProb << " beta: " << currBeta << "  alpha: " << currAlpha << "\n";
        	}

        	/* Update the number of MCMC steps we spent in an optimal state */
        	if(currScore == bestScore && it>=burnIn){
        		optStatesAfterBurnIn++;
        	}


        }
        delete [] currTreeParentVec;
        free_doubleMatrix(currLogScores);
        free_boolMatrix(currTreeAncMatrix);
	}                                              // last repetition of MCMC done

	unsigned int noStepsAfterBurnin = noOfReps*(noOfLoops-burnIn);
	cout.precision(17);
	cout << "best log score for tree:\t" << bestTreeLogScore <<  "\n";
	cout << "#optimal steps after burn-in:\t" << optStatesAfterBurnIn << "\n";
	cout << "total #steps after burn-in:\t" << noStepsAfterBurnin << "\n";
	cout << "%optimal steps after burn-in:\t" << (1.0*optStatesAfterBurnIn)/noStepsAfterBurnin << "\n";
	if(moveProbs[0]!=0.0){
		cout << "best value for beta:\t" << bestBeta << "\n";
		cout << "best value for alpha:\t" << bestAlpha << "\n";
		cout << "best relevant doublet rate:\t" << bestRelDoubletRate << "\n";
		cout << "best doublet rate:\t" << bestDoubletRate << "\n";
		cout << "best log score for (T, theta):\t" << bestScore << "\n";
	}

	return sampleOutput.str();
}


/**************        tree scoring        ******************/

/* This runs the MCMC for learning the tree and beta, or only the tree with a fixed beta, it samples from the posterior and/or records the optimal trees/beta */
double scoreTreeFastWithDoublets(int*parent, int n, int m, double** logScores, int** dataMatrix, bool** ancMatrix, double& propDoubletProb, double& propRelDoubletProb, int recMut, bool fixedDoubletProb){

	int treeSize = n+1;                   // for root node
	if(recMut>=0){ treeSize++; }       // in case of recurrent mutation
	int parVecSize = treeSize-1;
	int* bft = getBreadthFirstTraversal(parent, parVecSize);   // get breadth first traversal for tree
	vector<double> singletSumScores(m);                          // the sum of likelihoods for single attachments for each sample
	vector<double> doubletSumScores(m);                          // the sum of likelihood for doublet attachments for each sample
	vector<double> relDoubletScores(m);                       // the sum of likelihoods of the relevant doublets
	vector<double> L_0(m);                        // the average likelihood for singlet attachment
	vector<double> L_1(m);                        // the average likelihood for doublet attachment

	double relDoubletScore = 0.0;
	for(int sample=0; sample<m; sample++){    // get singlet and doublet probabilities for all samples

		double* singletScores;                         // attachment scores of sample as singlet to each node
		vector<vector<double> > doubletScores;   // attachment scores of sample as doublet to each node pair

		if(recMut < 0){                // no recurrent mutation
			singletScores = getAttachmentScoresFast(parent, parVecSize, logScores, dataMatrix[sample], bft);
			doubletScores = getDoubletAttachmentScoreMatrixFast(singletScores, parVecSize, bft, parent, logScores, dataMatrix[sample]);

		}
		else{
			int backMut = -1;    // parallel mutation
			if(ancMatrix[recMut][n]==1){ backMut = n ;}
			else if(ancMatrix[n][recMut]==1){ backMut = recMut ;}
			singletScores = getAttachmentScoresFast_recMut(parent, n, logScores, dataMatrix[sample], bft, recMut, backMut);
			doubletScores = getRecMutDoubletAttachmentScoreMatrixFast(singletScores, n, bft, ancMatrix, parent, logScores, dataMatrix[sample], recMut);
		}

		singletSumScores[sample] = getSampleSingletScoreFast(singletScores, parVecSize, bft);	// sum over all attachment points                                                                                                                  // there are no relevant doublets in a linear tree
		doubletSumScores[sample] = getSampleDoubletScoreFast(doubletScores, parVecSize, ancMatrix, relDoubletScore); // sum over all pairs of attachment points
		relDoubletScores[sample] = relDoubletScore;

		L_0[sample] = singletSumScores[sample] - log(treeSize);        // average over all attachment points
		L_1[sample] = doubletSumScores[sample] - 2*log(treeSize);      // average over all pairs of attachment points

		delete singletScores;
	}

	// optimize the doublet probability if the option was chosen
	if(!fixedDoubletProb){
		propDoubletProb = optimizeDoubletProb(L_0, L_1, m);    // computes the best p for this tree when averaging over sample attachments
	}

	double sumScore = 0.0;
	for(int sample=0; sample<m; sample++){
		relDoubletScore +=getDoubletProb(L_0[sample], L_1[sample], propDoubletProb) * relDoubletScores[sample];
		sumScore += getCombinedProb(L_0[sample], L_1[sample], propDoubletProb);
	}
	propRelDoubletProb = relDoubletScore/m;
	//cout << "sum score after:  " << sumScore << "\n";
	delete [] bft;
	//cout << "best p: " << best_p << "  best rel doublet rate: " << relDoubletScore << "\n";
	//getchar();
	return sumScore;
}


/* combined probability for singlet and doublet for given p: (1-p)*L_0 + p * L_1 */
double getCombinedProb(double L_0, double L_1, double p){

	double betterLogScore = max(L_0, L_1);
	return log((1-p)*exp(L_0-betterLogScore) + p*exp(L_1-betterLogScore)) + betterLogScore;
}

/* gets the probability of sample being a doublet given p */
double getDoubletProb(double L_0, double L_1, double p){
	double betterLogScore = max(L_0, L_1);
	//cout << L_0 << "  " << L_1 << "  computing new score: " << p*exp(L_1) << " / (" << (1-p)*exp(L_0) << " + " << p*exp(L_1) << ")\n";
	return p*exp(L_1-betterLogScore)/ ((1-p)*exp(L_0-betterLogScore) + p*exp(L_1-betterLogScore));
}


/**************        doublet tree scoring        ******************/

vector<vector<double> > getDoubletAttachmentScoreMatrixFast(double* singletAttachmentScore, int n, int* bft, int* parent, double** logScores, int* dataVector){

	vector<vector<double> > doubletAttachmentScore(n+1,vector<double>(n+1));
	int root = n;
	doubletAttachmentScore[root][root] = singletAttachmentScore[root];
	for(int i=1; i<=n; i++){                                      // score all attachment points in combination with root attachment and attachment pairs to the same node
		int node = bft[i];
		doubletAttachmentScore[node][root] = singletAttachmentScore[node];       // second cell attaches to root -> doublet score = singlet root score
		doubletAttachmentScore[root][node] = singletAttachmentScore[node];       // first cell attaches to root
		doubletAttachmentScore[node][node] = singletAttachmentScore[node];     // second cell attaches to same node
	}

	for(int i=1; i<=n; i++){                       // try all attachment points for the first cell of the doublet
		int first = bft[i];
		for(int j=i+1; j<=n; j++){                  // try all attachment points for the second cell (wlogs with larger bft-index)
			int second = bft[j];
			doubletAttachmentScore[first][second] = doubletAttachmentScore[first][parent[second]]; // doublet score for attaching second cell at parent of other node is known
			doubletAttachmentScore[first][second] -= logScores[dataVector[second]][0];            //  expected mutation state of the doublet for the mutation at second node
			doubletAttachmentScore[first][second] += logScores[dataVector[second]][1];            //  changes from 0 to 1; it cannot have already been 1 due to the bft and j>i
			doubletAttachmentScore[second][first] = doubletAttachmentScore[first][second];
		}
	}
	return doubletAttachmentScore;
}

vector<vector<double> > getRecMutDoubletAttachmentScoreMatrixFast(double* singletAttachmentScore, int n, int* bft, bool** anc, int* parent, double** logScores, int* obs, int recMut){

	vector<vector<double> > doubletAttachmentScore(n+2,vector<double>(n+2));
	int root = n+1;
	int smallerBftCopy = recMut;
	int biggerBftCopy = n;

	for(int i=0; i<=n; i++){
		if(bft[i] == n){
			smallerBftCopy = n;                        // the instance of the recurrent mutation that has the higher bft-index
			biggerBftCopy = recMut;                  // the instance of the recurrent mutation that has the lower bft-index
			break;
		}
		if(bft[i] == recMut){
			break;
		}
	}

	bool backmutation = false;
	if(anc[smallerBftCopy][biggerBftCopy]){ backmutation = true; }     // the two copies of the recurrent mutation are in one linage

	doubletAttachmentScore[root][root] = singletAttachmentScore[root];

	for(int i=1; i<=n+1; i++){                                          // score all attachment points in combination with root attachment and attachment pairs to the same node
		int node = bft[i];
		doubletAttachmentScore[node][root]  = singletAttachmentScore[node];     // second cell attaches to root -> doublet score = singlet root score
		doubletAttachmentScore[root][node]  = singletAttachmentScore[node];     // first cell attaches to root
		doubletAttachmentScore[node][node]  = singletAttachmentScore[node];     // second cell attaches to same node
	}

	for(int i=1; i<=n+1; i++){                 // try all non-root attachment points for the first cell of the doublet
		int smallerBftAtt = bft[i];

		if(anc[smallerBftCopy][smallerBftAtt] != anc[biggerBftCopy][smallerBftAtt]){  // recurrent mutation expected to be present due to placement of 1st copy of doublet

			for(int j=i+1; j<=n+1; j++){                                             // loop through all attachment points for the 2nd copy (wlog with larger bft-index)
				int biggerBftAtt = bft[j];                   // attachment point for 2nd copy

				if(biggerBftAtt == biggerBftCopy){  // attaching to the second copy of the recurrent mutation, no change compared to attachment to parent, as mutation already expected due to first attachment
					doubletAttachmentScore[smallerBftAtt][biggerBftAtt] = doubletAttachmentScore[smallerBftAtt][parent[biggerBftAtt]];
					doubletAttachmentScore[biggerBftAtt][smallerBftAtt] = doubletAttachmentScore[smallerBftAtt][biggerBftAtt];
				}
				else{                          // attaching to a node other than the second copy of the recurrent mutation -> default update
					doubletAttachmentScore[smallerBftAtt][biggerBftAtt] = doubletAttachmentScore[smallerBftAtt][parent[biggerBftAtt]];

					doubletAttachmentScore[smallerBftAtt][biggerBftAtt] -= logScores[obs[biggerBftAtt]][0];             //  expected mutation state of the doublet for the mutation at second node
					doubletAttachmentScore[smallerBftAtt][biggerBftAtt] += logScores[obs[biggerBftAtt]][1];             // changes from 0 to 1; it cannot have already been 1 due to the bft and j>i
					doubletAttachmentScore[biggerBftAtt][smallerBftAtt] = doubletAttachmentScore[smallerBftAtt][biggerBftAtt];
				}
			}
		}
		else{                                             // placement of first copy of doublet not enough to determine whether doublet mutation is expected to be present
			for(int j=i+1; j<=n+1; j++){                                 // loop through all attachment points for the 2nd copy (wlog with larger bft-index)
				int biggerBftAtt = bft[j];                            // attachment point for 2nd copy of doublet

				if( biggerBftAtt == biggerBftCopy && anc[smallerBftCopy][biggerBftAtt]){                                       // backmutation
					doubletAttachmentScore[smallerBftAtt][biggerBftAtt] = doubletAttachmentScore[smallerBftAtt][parent[biggerBftAtt]];
					doubletAttachmentScore[smallerBftAtt][biggerBftAtt] -= logScores[obs[recMut]][1];                 //  expected mutation state of the doublet for the mutation at second node
					doubletAttachmentScore[smallerBftAtt][biggerBftAtt] += logScores[obs[recMut]][0];                 // changes from 1 to 0 because backmutation
					doubletAttachmentScore[biggerBftAtt][smallerBftAtt] = doubletAttachmentScore[smallerBftAtt][biggerBftAtt];
				}
				else if(biggerBftAtt == smallerBftCopy || biggerBftAtt == biggerBftCopy){           // attaching to a parallel recurrent mutation
					doubletAttachmentScore[smallerBftAtt][biggerBftAtt] = doubletAttachmentScore[smallerBftAtt][parent[biggerBftAtt]];
					doubletAttachmentScore[smallerBftAtt][biggerBftAtt] -= logScores[obs[recMut]][0];             //  expected mutation state of the doublet for the mutation at second node
					doubletAttachmentScore[smallerBftAtt][biggerBftAtt] += logScores[obs[recMut]][1];             // changes from 0 to 1; it cannot have already been 1 due to the bft and j>i
					doubletAttachmentScore[biggerBftAtt][smallerBftAtt] = doubletAttachmentScore[smallerBftAtt][biggerBftAtt];
				}
				else{                // attaching to a node other than a copy of the recurrent mutation -> default update
					doubletAttachmentScore[smallerBftAtt][biggerBftAtt] = doubletAttachmentScore[smallerBftAtt][parent[biggerBftAtt]];
					doubletAttachmentScore[smallerBftAtt][biggerBftAtt] -= logScores[obs[biggerBftAtt]][0];             //  expected mutation state of the doublet for the mutation at second node
					doubletAttachmentScore[smallerBftAtt][biggerBftAtt] += logScores[obs[biggerBftAtt]][1];             // changes from 0 to 1; it cannot have already been 1 due to the bft and j>i
					doubletAttachmentScore[biggerBftAtt][smallerBftAtt] = doubletAttachmentScore[smallerBftAtt][biggerBftAtt];
				}
			}

		}
	}
	return doubletAttachmentScore;
}

/* computes L_1, the average likelihood over all doublet attachments of a sample */
/* also computes L_2, by updating relDoubletScore, */
double getSampleDoubletScoreFast(vector<vector<double> > doubletAttachmentScore, int n, bool** ancMatrix, double& relDoubletScore){

	vector<double> doubletScores;            // list of the scores of all attachment pairs
	vector<double> relevantDoubletScores;    // list of the scores only for the relevant pairs

	for(int i=0; i<=n; i++){                                             // try all attachment points for the first cell of the doublet
		doubletScores.push_back(doubletAttachmentScore[i][i]);
		for(int j=i+1; j<=n; j++){                                       // try all attachment points for the second cell
			doubletScores.push_back(doubletAttachmentScore[i][j]);        // add score to list
			doubletScores.push_back(doubletAttachmentScore[j][i]);
			if(i<n && j<n && ancMatrix[i][j]==0 && ancMatrix[j][i]==0){                        // if the doublet is relevant, add to list
				relevantDoubletScores.push_back(doubletAttachmentScore[i][j]);
				relevantDoubletScores.push_back(doubletAttachmentScore[j][i]);
			}
		}
	}

	double doubletScoreSum  = 0.0;
	double relevantDoubletScoreSum = 0.0;

	doubletScoreSum = getSumOfVectorElementsExpLog(doubletScores);                    // sum over all doublet scores
	relevantDoubletScoreSum = getSumOfVectorElementsExpLog(relevantDoubletScores);   // sum over all relevant doublet scores
	relDoubletScore = exp(relevantDoubletScoreSum-doubletScoreSum);                   // L_1,rel / L_1
	return doubletScoreSum;    // return average doublet score
}


/**************        singlet tree scoring        ******************/

/* computes the sum of likelihoods over all singlet attachments of a sample */
double getSampleSingletScoreFast(double* singletScores, int n, int*bft){
	double bestSingletScore = getMaxEntry(singletScores, n+1);                 // best of the scores (used to compute with score differences rather than scores)
	double sumSingletScore = 0.0;

	for(int i=0; i<=n; i++){                                               // sum over all attachment scores, exp is necessary as scores are actually log scores
		sumSingletScore += exp(singletScores[bft[i]]-bestSingletScore);    // subtraction of best score to calculate with score differences
		//cout << "for node " << bft[i] << ": " << singletScores[bft[i]] << "    " << exp(singletScores[bft[i]]) << "\n";
	}

	return log(sumSingletScore)+bestSingletScore;          // transform back to log scores and change from score differences to actual scores

}

/* computes an approximate scoring for a tree summing the score over all attachment points per sample */
/* this is basically the old scoring, and just used for comparison with the new scoring */
double sumScoreTreeFastSinglet(int n, int m, double** logScores, int** dataMatrix, int* parent){

	double sumTreeScore = 0.0;
	int* bft = getBreadthFirstTraversal(parent, n);   // get breadth first traversal for tree

	for(int sample=0; sample<m; sample++){
		double* scores = getAttachmentScoresFast(parent, n, logScores, dataMatrix[sample], bft); // attachments scores of sample to each node
		double bestMaxTreeScore = getMaxEntry(scores, n+1);                                     // best of the scores (used to compute with score differences rather than scores)

		double sumScore = 0.0;
		for(int i=0; i<=n; i++){                                 // sum over all attachment scores, exp is necessary as scores are actually log scores
			sumScore += exp(scores[bft[i]]-bestMaxTreeScore);   // subtraction of best score to calculate with score differences (smaller values)
		}
		delete [] scores;
		sumTreeScore += log(sumScore)+bestMaxTreeScore;      // transform back to log scores and change from score differences to actual scores

	}
	delete [] bft;
	return sumTreeScore;
}


/**************      doublet rate optimization      ******************/

/* optimize doublet probability p given the tree and sample attachment probabilities as singlets or doublets (EM type optimization)*/
double optimizeDoubletProb(vector<double> L_0, vector<double> L_1, int m){

	double p = 0.5;             // initial value
	double done = false;        // optimal p not yet found

	vector<double> sampleDoubletProb(m);
	int z=0;
	double currProbSum = 0.0;
	while(!done){
		double changed = false;
		double newProbSum = 0.0;


		z++; //cout << "round " << z++ << "\n";
		for(int sample=0; sample<m; sample++){
			double newProb = getDoubletProb(L_0[sample], L_1[sample], p);  // get the doublet probability for this sample
			newProbSum += newProb;


			sampleDoubletProb[sample] = newProb;


		}
		//cout.precision(17);
		//cout << "old/new p: " << p << "     " << newProbSum/m << "\n";
		if(abs(newProbSum-currProbSum) > doubletRatePrecision ){
		//if(newProbSum != currProbSum){    // p has not yet converged

			changed = true;
			currProbSum = newProbSum;
		}
		if(!changed){       // p has converged
			break;
		}
		p = newProbSum/m;     // take the average probability over all samples
		//getchar();
	}
	//cout  << "new p: " << p << " after " << z  << " rounds\n";
	//getchar();
	return p;    // this is the best p for this tree when averaging over sample attachments
}


/**************    sampling from posterior distribution     ******************/

/* prints out the current tree and beta to sample from the posterior distribution */
string sampleFromPosterior_doublets(double currTreeLogScore, int n, int* currTreeParentVec, double thetaProb, double currBeta, double currScore, double currAlpha, double currRelDoubletRate, double currDoubletRate){

	std::stringstream content;
	content << currTreeLogScore  << "\t";                 // logscore of current tree
	content << countBranches(currTreeParentVec, n);       // number of branches in current tree
	if(thetaProb>0.0){
		content << "\t" << currBeta;                      // current beta
		content << "\t" << currAlpha;                      // current alpha
		content << "\t" << currScore;                     // current combined logscore for tree and theta
		content << "\t" << currRelDoubletRate;                // current relevant doublet rate
		content << "\t" << currDoubletRate;                // current doublet rate
	}
	content << "\t";
	for(int i=0; i<n; i++){
		content << currTreeParentVec[i] << " ";
	}
	content << "\n";
	return content.str();
}


/**************    tree topology tests     ******************/

bool isInvalidRecMutTree(int* parent, int recMut,int copy){
	if(parent[recMut] == parent[copy]){
		return true;
	}
	if(parent[recMut]==copy){
		return true;
	}
	if(parent[copy]==recMut){
		return true;
	}
	return false;
}

/* returns true if the tree is a single chain, and false elsewise */
bool isLinearTree(int* parent, int* bft, int n){
	for(int i=0; i<n; i++){
		if(parent[bft[i]]==parent[bft[i+1]]){
			return false;
		}
	}
	return true;
}

