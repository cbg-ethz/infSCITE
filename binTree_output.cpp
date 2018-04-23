/*
 * binTree_output.cpp
 *
 *  Created on: Jan 15, 2018
 *      Author: jahnka
 */

#include <stdio.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <float.h>
#include "output.h"
#include "scoreTree.h"
#include "matrices.h"

using namespace std;


/* Score contribution by a specific mutation when placed at the root, that means all samples should have it */
/* This is the same for all trees and can be precomputed */
double binTreeRootScore(int** obsMutProfiles, int mut, int m, double ** logScores){
	double score = 0.0;
	for(int sample=0; sample<m; sample++){
		score += logScores[obsMutProfiles[sample][mut]][1];
	}
	return score;
}

/* computes the best placement of a mutation, the highest one if multiple co-optimal placements exist*/
int getHighestOptPlacement(int** obsMutProfiles, int mut, int m, double ** logScores, bool** ancMatrix){

	int nodeCount = (2*m)-1;
	int bestPlacement = (2*m)-2;   // root
	double bestPlacementScore = binTreeRootScore(obsMutProfiles, mut, m, logScores);
	//cout << bestPlacementScore << " (root)\n";
	//print_boolMatrix(bool** array, int n, int m);
	for(int p=0; p<nodeCount-1; p++){                           // try all possible placements (nodes in the mutation tree)

		double score = 0.0;                   // score for placing mutation at a specific node
		for(int sample=0; sample<m; sample++){
			//cout << p << " " << sample << "\n";
			if(ancMatrix[p][sample] == 1){
				score += logScores[obsMutProfiles[sample][mut]][1]; // sample should have the mutation
			}
			else{
				score += logScores[obsMutProfiles[sample][mut]][0]; // sample should not have the mutation
			}
		}
		if(score > bestPlacementScore){
			bestPlacement = p;
			bestPlacementScore = score;
			//cout << bestPlacementScore << " (non-root)\n";
		}
		else if (score == bestPlacementScore && ancMatrix[p][bestPlacement] == true){
			bestPlacement = p;
		}
	}

	//if(bestPlacement == (2*m)-2){
	//	cout<< "best placed at root\n";
	//	getchar();
	//}
	return bestPlacement;
}

string imputeGenotypesBinTree(int** obsMutProfiles, int n, int m, double ** logScores, bool** ancMatrix){

	stringstream content;
	content << m << endl;
	content << n << endl;

	bool** genotypes = init_boolMatrix(m, n, 0);

	for(int node=0; node< 2*m-2; node++){
		for(int otherNode=0; otherNode<2*m-2; otherNode++){
			cout << ancMatrix[node][otherNode];
		}
		cout << endl;
	}
	int* mutPlacement = getHighestOptPlacementVector(obsMutProfiles, n, m, logScores, ancMatrix);
	for(int mut=0; mut< n; mut++){
		for(int leaf=0; leaf<m; leaf++){
			//cout << "ancMatrix[" << mutPlacement[mut] << "][" << leaf << "]" << endl;
			if(mutPlacement[mut] == 2*m-2){
				genotypes[leaf][mut] = 1;       // root
			}
			else if(ancMatrix[mutPlacement[mut]][leaf]==1){
				genotypes[leaf][mut] = 1;
			}
		}
	}

	cout << "test" << endl;
	for(int sample=0; sample<m; sample++){
		content << sample+1 << "\t";
		for(int k=0; k<n; k++){
			content << genotypes[sample][k];
		}
		content << endl;
	}

//		for(int i=0; i<parentVectorSize; i++){
//			for(int j=0; j<parentVectorSize; j++){
//				content << ancMatrix[i][j];
//			}
//			content << endl;
//		}
		return content.str();
}

/* computes the best placement of a mutation, the highest one if multiple co-optimal placements exist*/
int* getHighestOptPlacementVector(int** obsMutProfiles, int n, int m, double ** logScores, bool** ancMatrix){
	int* bestPlacements = init_intArray(n, -1);
	for(int mut=0; mut<n; mut++){
		bestPlacements[mut] = getHighestOptPlacement(obsMutProfiles, mut, m, logScores, ancMatrix);
	 }
	//print_intArray(bestPlacements, n);
	return bestPlacements;
}

vector<string> getBinTreeNodeLabels(int nodeCount, int* optPlacements, int n, vector<string> geneNames, vector<string> sampleNames){
	vector<string> v;
	int count = 0;
	int sampleCount = (nodeCount+1)/2;
	for(int i = 0; i < nodeCount; i++){
		v.push_back("");
	}

	for(int mut=0; mut<n; mut++){
		string toAppend;
		int lastBreak = v.at(optPlacements[mut]).find_last_of("\n");
		string lastLine = v.at(optPlacements[mut]).substr(lastBreak+1);
		if(v.at(optPlacements[mut]) == ""){
			toAppend = geneNames.at(mut);
			count++;
		}
		else if(lastLine.length()>40){
			toAppend = "\n " + geneNames.at(mut);
			count++;
		}
		else{
			toAppend = " " + geneNames.at(mut);
			count++;
		}
		//cout << "        " << j << "\n";
		//cout << "                     "<< optPlacements[j] << "\n";

		v.at(optPlacements[mut]) += toAppend;
	}
	if(v.at(nodeCount-1) == ""){
		v.at(nodeCount-1) = "root";
	}

//	for(int i=0; i< sampleCount; i++){
//		if(v.at(i)==""){
//			v.at(i) += sampleNames.at(i);
//		}
//	}
	for(int i = 0; i < nodeCount; i++){
		if(v.at(i).find(" ") != string::npos){
			v.at(i) = "\"" + v.at(i) + "\"";
		}
	}

//	for(int i = 0; i < nodeCount; i++){
//		int mutCount = 1;
//		mutCount += (int)std::count(v.at(i).begin(),v.at(i).end(),',');
//		if(mutCount >=10){
//			v.at(i) = "\"+" + std::to_string(mutCount) + " " + std::to_string(i) + "\"";
//		}
//	}

	//cout << "added mutations " << count << "\n";

	return v;
}

/* returns the lca of a node that has a non-empty label, the root is assumed to always have a label */
int getLcaWithLabel(int node, int* parent, vector<string> label, int nodeCount){
	int root = nodeCount -1;
	int p = parent[node];;
	while(p != root && label[p]==""){
		p = parent[p];
	}
	return p;
}


/* create the content of the GraphViz output file for a mutation tree */
string getBinTreeGraphVizString(int* parentVector, int parentVectorSize, vector<string> nodeLabels, vector<string> sampleNames){

	cout << "..............................................\n";
	for(int i=0; i<parentVectorSize; i++){
		cout << nodeLabels[i] << endl;
	}

	for(int i=0; i<sampleNames.size(); i++){
		cout << sampleNames.at(i) << endl;
	}

	stringstream str;
	int leafCount = (parentVectorSize+2)/2;
	str << "digraph g{\n";
	str << "node [fontname = \"helvetica\"]";
	str << "node [color=lightgrey, style=filled, fontcolor=black, shape=circle];  \n";   // sample nodes
	for(int i=0; i<leafCount; i++){
		cout << nodeLabels[i] << " (nodeLabel) vs. " << sampleNames[i] << " (sampleLabel)" << endl;
		if(nodeLabels[i]!=sampleNames[i]){
			str << "node [color=deeppink4, style=filled, fontcolor=white, shape=box];	\n";
			//str << i << "[label=" << nodeLabels[i] << "];\n";
			str << i << "[label=" << nodeLabels[i] << "];\n";
			str << "node [color=lightgrey, style=filled, fontcolor=black, shape=circle];  \n";   // sample nodes
			str << "L" << i << "[label=" << sampleNames[i] << "];\n";
			str << i << " -> " << "L" << i << ";\n";

		}
		else{
			str << i << "[label=" << nodeLabels[i] << "];\n";
		}
	}

	str << "node [color=deeppink4, style=filled, fontcolor=white, shape=box];	\n";   // style of tree nodes
	for(int i=leafCount; i<nodeLabels.size(); i++){
		str << i << "[label=" << nodeLabels[i] << "];\n";
		//str << i << "[label=" << i << "];\n";
	}

	for(int i=0; i<parentVectorSize; i++){
		str << parentVector[i] << " -> " << i << ";\n";
	}
	str << "}\n";
	return str.str();
}

std::string getGraphVizBinTree(int* parents, int nodeCount, int m, vector<string> label){
	std::stringstream content;
	content << "digraph G {\n";
	content << "node [color=deeppink4, style=filled, fontcolor=white, fontsize=20, fontname=Verdana];\n";
	for(int i=m; i<nodeCount-1; i++){
		if(label[i] != ""){
			int labelledLCA = getLcaWithLabel(i, parents, label, nodeCount);
			content << label[labelledLCA] << " -> " << label[i] << ";\n";
//		if(label[parents[i]] == ""){
//			content  << parents[i] << " -> ";
//		}
//		else{
//			content << label[parents[i]] << " -> ";
//		}
//		if(label[i] == ""){
//		  content  << i << ";\n";
//		}
//		else{
//			content << label[i] << ";\n";

		}
	}
	content << "node [color=lightgrey, style=filled, fontcolor=black];\n";
	for(int i=0; i<m; i++){
		int labelledLCA = getLcaWithLabel(i, parents, label, nodeCount);
		content << label[labelledLCA] << " -> " << "s" << i << ";\n";



//		if(label[parents[i]] == ""){
//			content << parents[i] << " -> ";
//		}
//		else{
//			content << label[parents[i]] << " -> ";
//		}
//
//		content << "s" << i << ";\n";


	}
	content <<  "}\n";
	return content.str();
}
