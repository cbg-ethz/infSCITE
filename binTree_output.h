/*
 * binTree_output.h
 *
 *  Created on: Jan 18, 2018
 *      Author: jahnka
 */

#ifndef BINTREE_OUTPUT_H_
#define BINTREE_OUTPUT_H_

using namespace std;

string getBinTreeGraphVizString(int* parentVector, int parentVectorSize, vector<string> nodeLabels, vector<string> sampleNames);
double binTreeRootScore(int** obsMutProfiles, int mut, int m, double ** logScores);
int getHighestOptPlacement(int** obsMutProfiles, int mut, int m, double ** logScores, bool** ancMatrix);
int* getHighestOptPlacementVector(int** obsMutProfiles, int n, int m, double ** logScores, bool** ancMatrix);
vector<string> getBinTreeNodeLabels(int nodeCount, int* optPlacements, int n, vector<string> geneNames, vector<string> sampleNames);
int getLcaWithLabel(int node, int* parent, vector<string> label, int nodeCount);
std::string getGraphVizBinTree(int* parents, int nodeCount, int m, vector<string> label);
string imputeGenotypesBinTree(int** obsMutProfiles, int n, int m, double ** logScores, bool** ancMatrix);

#endif /* BINTREE_OUTPUT_H_ */
