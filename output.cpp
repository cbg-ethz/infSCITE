/*
 * output.cpp
 *
 *  Created on: Oct 12, 2015
 *      Author: jahnka
 */

#include <stdio.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <float.h>
#include <math.h>
#include "output.h"
#include "scoreTree.h"
#include "trees.h"
#include "recMut.h"
#include "matrices.h"

using namespace std;




/* writes the given string to file */
void writeToFile(string content, string fileName){
	ofstream outfile;
	outfile.open (fileName.c_str());
	outfile << content;
	outfile.close();
}


/* create the content of the GraphViz output file for a mutation tree */
string getMutTreeGraphVizString(int* parentVector, int parentVectorSize, vector<string> nodeLabels, bool attachSamples, vector<string> sampleLabels, int** data, double** logScores, bool recMutAllowed, int recMut){

	//cout << "writing tree file" << endl;
	stringstream str;
	str << "digraph g{\n";
	str << "node [color=deeppink4, style=filled, fontcolor=white];	\n";   // style of tree nodes

	str << mutTreeNodes(nodeLabels);                                  // tree nodes
	str << mutTreeEdges(parentVector, parentVectorSize);              // their edges

	if(attachSamples){
		str << "node [color=lightgrey, style=filled, fontcolor=black, shape=circle, fixedsize=true];  \n";   // sample nodes
		if(sampleLabels.size() < 200){
			str << mutTreeSamples(sampleLabels);                                  // sample attachment edges
			str << sampleAttachment(sampleLabels.size(), parentVector, parentVectorSize, data, logScores, recMutAllowed, recMut);
		}
		else{
			cout << "attachment summary" << endl;
			str << sampleAttachmentSummary(sampleLabels.size(), parentVector, parentVectorSize, data, logScores, recMutAllowed, recMut);
		}
	}
	str << "}\n";
	//cout << "tree file written" << endl;
	return str.str();
}

/* Generates the output lines listing nodes of the mutation tree */
string mutTreeNodes(vector<string> nodeLabels){
	stringstream nodes;
	for(int i=0; i<nodeLabels.size(); i++){
		nodes << i << "[label=\"" << nodeLabels[i] << "\"];\n";
	}
	return nodes.str();
}


/* Generates the output lines listing sample nodes of the mutation tree */
string mutTreeSamples(vector<string> sampleLabels){
	stringstream nodes;
	for(int i=0; i<sampleLabels.size(); i++){
		nodes << "s_" << i << "[label=\"" << sampleLabels[i] << "\"];\n";
	}
	return nodes.str();
}

/* Generates the output lines listing samples attached to the mutation tree */
string mutTreeEdges(int* parentVector, int parentVectorSize){
	stringstream edges;
	for(int i=0; i<parentVectorSize; i++){
		edges << parentVector[i] << " -> " << i << ";\n";
	}
	return edges.str();
}

/* Generates the output lines listing the edges attaching the samples to the mutation tree */
string sampleAttachment(int sampleCount, int* parentVector, int parentVectorSize, int** data, double** logScores, bool recMutAllowed, int recMut){

	stringstream attachment;
	for(int i=0; i< sampleCount; i++){
		vector<int> attachmentPoints = getSampleAttachmentPoints(parentVector, parentVectorSize, data[i], logScores, recMutAllowed, recMut);
		for(int j=0; j<attachmentPoints.size(); j++){
			attachment << attachmentPoints[j] << " -> s_" << i << ";\n";
		}
	}
	return attachment.str();
}

string imputeGenotypes(int sampleCount, int* parentVector, int parentVectorSize, int** data, double** logScores, bool recMutAllowed, int recMut){
	stringstream content;
	content << sampleCount << endl;
	content << parentVectorSize << endl;

	bool** ancMatrix = parentVector2ancMatrix(parentVector, parentVectorSize);
	for(int i=0; i< sampleCount; i++){
		vector<int> attachmentPoints = getSampleAttachmentPoints(parentVector, parentVectorSize, data[i], logScores, recMutAllowed, recMut);
		for(int j=0; j<attachmentPoints.size(); j++){
			content << i+1 << "\t";
			for(int k=0; k<parentVectorSize; k++){
				content << ancMatrix[k][attachmentPoints[j]];
			}
			content << endl;
		}
	}
	for(int i=0; i<parentVectorSize; i++){
		for(int j=0; j<parentVectorSize; j++){
			content << ancMatrix[i][j];
		}
		content << endl;
	}
	return content.str();
}

/* Generates the output lines listing the edges attaching the samples to the mutation tree */
string sampleAttachmentSummary(int sampleCount, int* parentVector, int parentVectorSize, int** data, double** logScores, bool recMutAllowed, int recMut){

	stringstream attachment;
	int nodeCount = parentVectorSize+1;
	//double* attachmentCounter = init_doubleArray(nodeCount, 0.0);
	int* attachmentCounter = init_intArray(nodeCount, 0);

	// attach sample proportions to tree
//	for(int i=0; i< sampleCount; i++){
//		vector<int> attachmentPoints = getSampleAttachmentPoints(parentVector, parentVectorSize, data[i], logScores, recMutAllowed, recMut);
//		double fraction = 1.0/attachmentPoints.size();
//		for(int j=0; j<attachmentPoints.size(); j++){
//			//cout << attachmentPoints[j] << " ";
//			attachmentCounter[attachmentPoints[j]] += fraction;
//		}
//		//cout << endl;
//	}
//	for(int i=0; i<nodeCount; i++){
//		if(attachmentCounter[i]>0.0){
//			attachment << "s_" << i << "[label=" << attachmentCounter[i] << ", width=" << 50.0*attachmentCounter[i]/sampleCount << "];\n";
//			attachment << i << " -> " << "s_" << i << ";\n";
//		}
//	}

	bool** ancMatrix = parentVector2ancMatrix(parentVector, parentVectorSize);
	for(int i=0; i< sampleCount; i++){
		vector<int> attachmentPoints = getSampleAttachmentPoints(parentVector, parentVectorSize, data[i], logScores, recMutAllowed, recMut);
		//cout << "sample attachment point vector written" << endl;
		int attPoint = attachmentPoints.at(0);
		for(int j=1; j<attachmentPoints.size(); j++){
			if(attachmentPoints.at(j)==parentVectorSize){
				attPoint = attachmentPoints.at(j);
			}
			else if(ancMatrix[attachmentPoints.at(j)][attPoint]){
				attPoint = attachmentPoints.at(j);
			}
		}
		attachmentCounter[attPoint]++;
	}
	//cout << "attaching samples" << endl;
	for(int i=0; i<nodeCount; i++){
		if(attachmentCounter[i]>0){
			attachment << "s_" << i << "[label=" << attachmentCounter[i] << ", width=" << sqrt(500.0*attachmentCounter[i]/sampleCount) << "];\n";
			attachment << i << " -> " << "s_" << i << ";\n";
		}
	}

	return attachment.str();
}

/* computes the score for each attachment point for given sample */
vector<int> getSampleAttachmentPoints(int* parentVector, int parentVectorSize, int* obsMutStates, double** logScores, bool recMutAllowed, int recMut){

	double bestAttachmentScore = -DBL_MAX;
	vector<int> attachmentPoints;
	double* scores;

	if(recMutAllowed){
		scores = getAttachmentScoresFast_recMut(parentVector, parentVectorSize-1, logScores, obsMutStates, getBreadthFirstTraversal(parentVector, parentVectorSize), recMut, getBackmutation(parentVector, recMut, parentVectorSize-1));
	}
	else{
		scores = getAttachmentScoresFast(parentVector, parentVectorSize,logScores, obsMutStates, getBreadthFirstTraversal(parentVector, parentVectorSize));
	}

//	for(int i=0; i<parentVectorSize+1;i++){
//		cout << scores[i] << "\t";
//	}
//	cout << endl;
	for(int i=0; i<parentVectorSize+1;i++){
		if(scores[i] >= bestAttachmentScore){
			if(scores[i] > bestAttachmentScore){
			  	attachmentPoints.clear();                // new attachment point better than current best, empty list
			  	bestAttachmentScore = scores[i];
			}
			attachmentPoints.push_back(i);               // new attachment point is current best, add to list
		}
	}
	return attachmentPoints;
}

//	vector<int> sampleAttachmentPoints(bool ** ancMatrix, int nonRootNodeCount, int mutCount, int sampleId, double** logScores, int** dataMatrix){
//
//	    double treeScore = 0.0;
//	  	double bestAttachmentScore = -DBL_MAX;                            // best score for attaching sample
//	  	vector<int> attachmentPoints;
//
//	  	double attachmentScore = 0.0;
//	  	for(int mut=0; mut<mutCount; mut++){                                // start with attaching node to root (no genes mutated)
//	  		attachmentScore += logScores[dataMatrix[sample][mut]][0];
//	  	}
//	  	if(attachmentScore >= bestAttachmentScore){                       // new attachment point as good as old one
//	  		if(attachmentScore > bestAttachmentScore){
//	  			attachmentPoints.clear();
//	  		}
//	  		attachmentPoints.push_back(nonRootNodeCount);                 // add to list
//	  	}
//
//	  	for(int parent=0; parent<nonRootNodeCount; parent++){      // try all attachment points (genes)
//	  		attachmentScore=0.0;
//	  		for(int mut=0; mut<n; mut++){     // sum up scores for each gene, score for zero if gene is not an ancestor of parent, score for one else wise
//	  		    attachmentScore += logScores[dataMatrix[sample][mut]][ancMatrix[mut][parent]];
//	  		}
//	  		if(attachmentScore > bestAttachmentScore){
//	  			bestAttachmentScore = attachmentScore;
//	  		}
//	  	}
//	  		for(int parent=0; parent<n; parent++){      // try all attachment points (genes)
//	  		 	double attachmentScore=0.0;
//	  		 	for(int gene=0; gene<n; gene++){     // sum up scores for each gene, score for zero if gene is not an ancestor of parent, score for one else wise
//	  		 		attachmentScore += logScores[dataMatrix[sample][gene]][ancMatrix[gene][parent]];
//	  		 	}
//	  		  	if(attachmentScore == bestAttachmentScore){
//	  		  		attachment[parent][sample] = true;
//	  		  	}
//	  		}
//	  		bool rootAttachment = true;
//	  		for(int parent=0; parent<n; parent++){
//	  			if(attachment[parent][sample] == true){
//	  				rootAttachment = false;
//	  				break;
//	  			}
//	  		}
//	  		if(rootAttachment == true){
//	  			attachment[n][sample] = true;
//	  		}
//	  		treeScore += bestAttachmentScore;
//
//	}

/***********       end new code      ***********/


//string getMutTreeGraphViz(vector<string> label, int nodeCount, int m, int* parent){
//	stringstream nodes;
//	stringstream leafedges;
//	stringstream edges;
//	for(int i=0; i<m; i++){
//		if(label.at(i) != ""){
//			nodes << "s" << i << "[label=\"s" << i << "\"];\n";                 // si [label="si"];
//			nodes        << i << "[label=\"" << label.at(i) << "\"];\n";                 //   i [label="i"];
//			leafedges << "s" << i << " -> " << i << ";\n";
//			edges <<        i << " -> " << getLcaWithLabel(i, parent, label, nodeCount) << ";\n";
//		}
//		else{
//			nodes << i << "[label=\"s" << i << "\"];\n";
//			leafedges << i << " -> " << getLcaWithLabel(i, parent, label, nodeCount) << ";\n";
//		}
//	}
//
//	stringstream str;
//
//	str << "digraph g{\n";
//	str << nodes;
//	str << "node [color=deeppink4, style=filled, fontcolor=white];	\n";
//	str << edges;
//	str << "node [color=lightgrey, style=filled, fontcolor=black];  \n";
//	str << leafedges;
//	str << "}\n";
//	return str.str();
//}




/*********** */


/***************/

///* creates the content for the GraphViz file from a parent vector, using numbers as node labels (from 1 to n+1) */
//std::string getGraphVizFileContentNumbers(int* parents, int n){
//	std::stringstream content;
//	content << "digraph G {\n";
//	content << "node [color=deeppink4, style=filled, fontcolor=white];\n";
//	for(int i=0; i<n; i++){
//		content << parents[i]+1  << " -> "  << i+1 << ";\n";      // plus 1 to start gene labeling at 1 (instead of 0)
//	}
//	content <<  "}\n";
//	return content.str();
//}
//
//
///* creates the content for the GraphViz file from a parent vector, using the gene names as node labels */
//std::string getGraphVizFileContentNames(int* parents, int n, vector<string> geneNames, bool attachSamples, bool** ancMatrix, int m, double** logScores, int** dataMatrix){
//	std::stringstream content;
//	content << "digraph G {\n";
//	content << "node [color=deeppink4, style=filled, fontcolor=white];\n";
//
//	for(int i=0; i<geneNames.size()-1; i++){
//		content << "\"" << geneNames[parents[i]] << "\"" << " -> "  << "\"" << geneNames[i]  << "\""  << ";\n";
//	}
//
//	if(attachSamples==true){
//
//		content << "node [color=lightgrey, style=filled, fontcolor=black];\n";
//		std::string attachment = getBestAttachmentString(ancMatrix, n, m, logScores, dataMatrix, geneNames);
//		content << attachment;
//	}
//	content <<  "}\n";
//	return content.str();
//}
//
///* creates the content for the GraphViz file from a parent vector, using the gene names as node labels */
//std::string getGraphVizFileContentNew(int* parents, int n, vector<string> geneNames, vector<string> sampleNames, bool attachSamples, bool** ancMatrix, int m, double** logScores, int** dataMatrix){
//	std::stringstream content;
//	content << "digraph G {\n";
//	content << "node [color=deeppink4, style=filled, fontcolor=white];\n";   // mutation node style
//
//	for(int i=0; i<geneNames.size(); i++){                      // define mutation nodes
//		content << i << " [label=" << geneNames[i] << "]";
//	}
//
//	for(int i=0; i<geneNames.size()-1; i++){                 // define edges between mutations
//		content << parents[i] << " -> "  << i  << ";\n";
//	}
//
//	if(attachSamples==true){
//
//		content << "node [color=lightgrey, style=filled, fontcolor=black];\n";   // sample node style
//
//		for(int i=0; i<sampleNames.size(); i++){                            // define sample nodes
//			content << "s_" << i << " [label=" << sampleNames[i] << "]";
//		}
//
//		vector<int> attachmentPoints = getBestAttachmentPoints(ancMatrix, n, m, logScores, dataMatrix, geneNames);
//		for(int i; i< attachmentPoints.size(); i++){
//		content << attachmentPoints.at(i) << " -> " << "s_" < ;
//		}
//	}
//	content <<  "}\n";
//	return content.str();
//}
//
///* creates the attachment string for the samples, the optimal attachment points are recomputed from scratch based on error log Scores */
//vector<vector<int> > getBestAttachmentPoints(bool ** ancMatrix, int n, int m, double** logScores, int** dataMatrix, vector<string> geneNames){
//	bool** matrix = attachmentPoints(ancMatrix, n, m, logScores, dataMatrix);
//	vector<vector<int> > points;
//	for(int i=0; i<=n; i++){
//		vector<int> list = new vector<int>;
//		points.push_back(list);
//	}
//
//	for(int i=0; i<=n; i++){
//		for(int j=0; j<m; j++){
//			if(matrix[i][j]==true){
//				points.at(j).push_back(i);
//			}
//		}
//	}
//	return points;
//}

///* creates the attachment string for the samples, the optimal attachment points are recomputed from scratch based on error log Scores */
//std::string getBestAttachmentString(bool ** ancMatrix, int n, int m, double** logScores, int** dataMatrix, vector<string> geneNames){
//	bool** matrix = attachmentPoints(ancMatrix, n, m, logScores, dataMatrix);
//	std::stringstream a;
//	for(int i=0; i<=n; i++){
//		for(int j=0; j<m; j++){
//			if(matrix[i][j]==true){
//				a << geneNames[i] << " -> s" << j << ";\n";
//			}
//		}
//	}
//	return a.str();
//}

///* This is a re-computation of the best attachment points of the samples to a tree for printing the tree with sample attachment. */
///* It uses the ancestor matrix of the tree and returns a bit matrix indicating the best attachment points of each sample based on the error log scores */
//vector<int> sampleAttachmentPoints(bool ** ancMatrix, int nonRootNodeCount, int mutCount, int sampleId, double** logScores, int** dataMatrix){
//
//    double treeScore = 0.0;
//  	double bestAttachmentScore = -DBL_MAX;                            // best score for attaching sample
//  	vector<int> attachmentPoints;
//
//  	double attachmentScore = 0.0;
//  	for(int mut=0; mut<mutCount; mut++){                                // start with attaching node to root (no genes mutated)
//  		attachmentScore += logScores[dataMatrix[sample][mut]][0];
//  	}
//  	if(attachmentScore >= bestAttachmentScore){                       // new attachment point as good as old one
//  		if(attachmentScore > bestAttachmentScore){
//  			attachmentPoints.clear();
//  		}
//  		attachmentPoints.push_back(nonRootNodeCount);                 // add to list
//  	}
//
//  	for(int parent=0; parent<nonRootNodeCount; parent++){      // try all attachment points (genes)
//  		attachmentScore=0.0;
//  		for(int mut=0; mut<n; mut++){     // sum up scores for each gene, score for zero if gene is not an ancestor of parent, score for one else wise
//  		    attachmentScore += logScores[dataMatrix[sample][mut]][ancMatrix[mut][parent]];
//  		}
//  		if(attachmentScore > bestAttachmentScore){
//  			bestAttachmentScore = attachmentScore;
//  		}
//  	}
//  		for(int parent=0; parent<n; parent++){      // try all attachment points (genes)
//  		 	double attachmentScore=0.0;
//  		 	for(int gene=0; gene<n; gene++){     // sum up scores for each gene, score for zero if gene is not an ancestor of parent, score for one else wise
//  		 		attachmentScore += logScores[dataMatrix[sample][gene]][ancMatrix[gene][parent]];
//  		 	}
//  		  	if(attachmentScore == bestAttachmentScore){
//  		  		attachment[parent][sample] = true;
//  		  	}
//  		}
//  		bool rootAttachment = true;
//  		for(int parent=0; parent<n; parent++){
//  			if(attachment[parent][sample] == true){
//  				rootAttachment = false;
//  				break;
//  			}
//  		}
//  		if(rootAttachment == true){
//  			attachment[n][sample] = true;
//  		}
//  		treeScore += bestAttachmentScore;
//
//}

///* This is a re-computation of the best attachment points of the samples to a tree for printing the tree with sample attachment. */
///* It uses the ancestor matrix of the tree and returns a bit matrix indicating the best attachment points of each sample based on the error log scores */
//bool** attachmentPoints(bool ** ancMatrix, int n, int m, double** logScores, int** dataMatrix){
//
//    double treeScore = 0.0;
//    bool ** attachment = init_boolMatrix(n+1, m, false);
//
//  	for(int sample=0; sample<m; sample++){       // foreach sample
//  		double bestAttachmentScore = 0.0;     // currently best score for attaching sample
//  		for(int gene=0; gene<n; gene++){   // start with attaching node to root (no genes mutated)
//  			bestAttachmentScore += logScores[dataMatrix[sample][gene]][0];
//  		}
//  		for(int parent=0; parent<n; parent++){      // try all attachment points (genes)
//  		    double attachmentScore=0.0;
//  		    for(int gene=0; gene<n; gene++){     // sum up scores for each gene, score for zero if gene is not an ancestor of parent, score for one else wise
//  		    	attachmentScore += logScores[dataMatrix[sample][gene]][ancMatrix[gene][parent]];
//  		    }
//  		    if(attachmentScore > bestAttachmentScore){
//  		        bestAttachmentScore = attachmentScore;
//  		    }
//  		}
//  		for(int parent=0; parent<n; parent++){      // try all attachment points (genes)
//  		 	double attachmentScore=0.0;
//  		 	for(int gene=0; gene<n; gene++){     // sum up scores for each gene, score for zero if gene is not an ancestor of parent, score for one else wise
//  		 		attachmentScore += logScores[dataMatrix[sample][gene]][ancMatrix[gene][parent]];
//  		 	}
//  		  	if(attachmentScore == bestAttachmentScore){
//  		  		attachment[parent][sample] = true;
//  		  	}
//  		}
//  		bool rootAttachment = true;
//  		for(int parent=0; parent<n; parent++){
//  			if(attachment[parent][sample] == true){
//  				rootAttachment = false;
//  				break;
//  			}
//  		}
//  		if(rootAttachment == true){
//  			attachment[n][sample] = true;
//  		}
//  		treeScore += bestAttachmentScore;
//  	}
//  	return attachment;
//}


///* prints all trees in list of optimal trees to the console, first as parent vector, then as GraphViz file */
//void printParentVectors(vector<bool**> optimalTrees, int n, int m, double** logScores, int** dataMatrix){
//	for(int i=0; i<optimalTrees.size(); i++){
//		int* parents = ancMatrixToParVector(optimalTrees[i], n);
//		print_intArray(parents,n);
//		//print_boolMatrix(attachmentPoints(optimalTrees[i], n, m, logScores, dataMatrix), n, m);
//		printGraphVizFile(parents, n);
//	}
//}


///* prints the GraphViz file for a tree to the console */
//void printGraphVizFile(int* parents, int n){
//	cout << "digraph G {\n";
//	cout << "node [color=deeppink4, style=filled, fontcolor=white];\n";
//	for(int i=0; i<n; i++){
//		cout << parents[i] << " -> " << i << "\n";
//	}
//	cout << "}\n";
//}

void printSampleTrees(vector<int*> list, int n, string fileName){
	if(list.size()==0){ return;}
	std::stringstream a;
	for(int i=0; i<list.size(); i++){
		for(int j=0; j<n; j++){
			a << list[i][j];
			if(j<n-1){
				a  << " ";
			}
		}
		a << "\n";
	}
	writeToFile(a.str(), fileName);
	cout << "Trees written to: " << fileName;
}

/* prints the score of the tree predicted by the Kim&Simon approach for the given error log scores */
void printScoreKimSimonTree(int n, int m, double** logScores, int** dataMatrix, char scoreType){
	int parent[] = {2, 4, 17, 2, 9, 9, 2, 2, 4, 18, 2, 1, 2, 2, 9, 2, 2, 11};
	double KimSimonScore = scoreTree(n, m, logScores, dataMatrix, scoreType, parent, -DBL_MAX);
	cout.precision(20);
	cout << "KimSimonScore: " << KimSimonScore << "\n";
}


