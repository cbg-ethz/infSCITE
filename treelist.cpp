/*
 * treelist.cpp
 *
 *  Created on: Mar 27, 2015
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
#include <queue>
#include "matrices.h"
#include "treelist.h"
#include "rand.h"

using namespace std;


void updateTreeList(vector<struct treeTheta>& bestTrees, int* currTreeParentVec, int n, double currScore, double bestScore, double beta, double alpha){

	if(currScore > bestScore){
		//cout << "tree list of size " << bestTrees.size() << " emptied\n";
		resetTreeList(bestTrees, currTreeParentVec, n, beta, alpha);                              // empty the list of best trees and insert current tree

	}
	else if (currScore == bestScore){
		if(!isDuplicateTreeFast(bestTrees, currTreeParentVec, n)){               // if the same tree was not previously found
			treeTheta newElem = createNewTreeListElement(currTreeParentVec, n, beta, alpha);
			bestTrees.push_back(newElem);        // add it to list
		}
	}
}


/* removes all elements from the vector and inserts the new best tree */
void resetTreeList(vector<struct treeTheta>& bestTrees, int* newBestTree, int n, double beta, double alpha){
	emptyVectorFast(bestTrees, n);                                         // empty the list of best trees
	treeTheta newElem = createNewTreeListElement(newBestTree, n, beta, alpha);
	bestTrees.push_back(newElem);                // current tree is now the only best tree
}


/* removes all elements from the vector */
void emptyVectorFast(std::vector<struct treeTheta>& optimalTrees, int n){
    for(int i=0; i<optimalTrees.size(); i++){
    	delete [] optimalTrees[i].tree;
	}
    optimalTrees.clear();
}

/* removes all elements from the vector */
void emptyTreeList(std::vector<int*>& optimalTrees, int n){
    for(int i=0; i<optimalTrees.size(); i++){
    	delete [] optimalTrees[i];
	}
    optimalTrees.clear();
}

/* creates a new tree/beta combination */
struct treeTheta createNewTreeListElement(int* tree, int n, double beta, double alpha){
	treeTheta newElem;
	newElem.tree = deepCopy_intArray(tree, n);
	newElem.beta = beta;
	newElem.alpha = alpha;
	return newElem;
}

/* returns true if the same tree was found before */
bool isDuplicateTreeFast(std::vector<struct treeTheta> &optimalTrees, int* newTree, int n){
    for(int k=0; k<optimalTrees.size(); k++){
      bool same = true;
      for(int i=0; i<n; i++){
    	  if(newTree[i] != optimalTrees[k].tree[i]){
              same = false;
              break;
          }
      }
      if(same == true){
        return true;
      }
    }
    return false;
}
