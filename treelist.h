/*
 * treelist.h
 *
 *  Created on: Mar 27, 2015
 *      Author: jahnka
 */

#include <string>
//#include <iostream>
//#include <sstream>

#ifndef TREELIST_H
#define TREELIST_H

struct treeTheta
{
	int* tree;
	double alpha;
    double beta;
};

void updateTreeList(std::vector<struct treeTheta>& bestTrees, int* currTreeParentVec, int n, double currScore, double bestScore, double beta, double alpha);
void resetTreeList(std::vector<struct treeTheta>& bestTrees, int* newBestTree, int n, double beta, double alpha);
void emptyVectorFast(std::vector<struct treeTheta>& optimalTrees, int n);
void emptyTreeList(std::vector<int*>& optimalTrees, int n);
struct treeTheta createNewTreeListElement(int* tree, int n, double beta, double alpha);
bool isDuplicateTreeFast(std::vector<struct treeTheta> &optimalTrees, int* newTree, int n);

#endif
