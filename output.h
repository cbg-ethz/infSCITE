/*
 * output.h
 *
 *  Created on: Oct 12, 2015
 *      Author: jahnka
 */

#ifndef OUTPUT_H_
#define OUTPUT_H_

using namespace std;

void writeToFile(string content, string fileName);
string getMutTreeGraphVizString(int* parentVector, int parentVectorSize, vector<string> nodeLabels, bool attachSamples, vector<string> sampleLabels, int** data, double** logScores, bool recMutAllowed, int recMut);
string mutTreeNodes(vector<string> nodeLabels);
string mutTreeSamples(vector<string> sampleLabels);
string mutTreeEdges(int* parentVector, int parentVectorSize);
string sampleAttachment(int sampleCount, int* parentVector, int parentVectorSize, int** data, double** logScores, bool recMutAllowed, int recMut);
vector<int> getSampleAttachmentPoints(int* parentVector, int parentVectorSize, int* obsMutStates, double** logScores, bool recMutAllowed, int recMut);
string sampleAttachmentSummary(int sampleCount, int* parentVector, int parentVectorSize, int** data, double** logScores, bool recMutAllowed, int recMut);
string imputeGenotypes(int sampleCount, int* parentVector, int parentVectorSize, int** data, double** logScores, bool recMutAllowed, int recMut);
double binTreeRootScore(int** obsMutProfiles, int mut, int m, double ** logScores);
int getHighestOptPlacement(int** obsMutProfiles, int mut, int m, double ** logScores, bool** ancMatrix);
int* getHighestOptPlacementVector(int** obsMutProfiles, int n, int m, double ** logScores, bool** ancMatrix);
std::vector<std::string> getBinTreeNodeLabels(int nodeCount, int* optPlacements, int n, std::vector<std::string> geneNames);
int getLcaWithLabel(int node, int* parent, std::vector<std::string> label, int nodeCount);
std::string getGraphVizBinTree(int* parents, int nodeCount, int m, std::vector<std::string> label);
std::string getMutTreeGraphViz(std::vector<std::string> label, int nodeCount, int m, int* parent);
void writeToFile(std::string content, std::string fileName);
std::string getGraphVizFileContentNumbers(int* parents, int n);
std::string getGraphVizFileContentNames(int* parents, int n, std::vector<std::string> geneNames, bool attachSamples, bool** ancMatrix, int m, double** logScores, int** dataMatrix);
std::string getBestAttachmentString(bool ** ancMatrix, int n, int m, double** logScores, int** dataMatrix, std::vector<std::string> geneNames);
bool** attachmentPoints(bool ** ancMatrix, int n, int m, double** logScores, int** dataMatrix);
void printParentVectors(std::vector<bool**> optimalTrees, int n, int m, double** logScores, int** dataMatrix);
void printGraphVizFile(int* parents, int n);
void printSampleTrees(std::vector<int*> list, int n, std::string fileName);
void printScoreKimSimonTree(int n, int m, double** logScores, int** dataMatrix, char scoreType);

#endif /* OUTPUT_H_ */
