/*
 * enum.cpp
 *
 *  Created on: Mar 27, 2015
 *      Author: jahnka
 */

#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdbool.h>
#include <vector>
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
//#include <cmath>
#include "treelist.h"
#include "trees.h"
#include "matrices.h"
#include "enum.h"

#include "rand.h"
#include "scoreTree.h"
// generate random tree with n nodes
// attach m nodes randomly to the tree
// get data matrix (m x n)
// change positions according to error rates (try different rates)
// run program with data matrix
// get score of best trees
// compare tree structures
// repeat x times
using namespace std;

/* enumerates all trees and finds the best tree */
double findBestTreeByExhaustiveSearch(int n, int m, std::vector<int*> &bestTrees, double** logScores, int** dataMatrix, char scoreType){

    int noOfTrees = pow(n+1, n-1);    // number of tree topologies
    double bestScore = -DBL_MAX;          // best score found so far
    printf("number of trees: %d\n", noOfTrees);
    printf("initial best score: %e\n", bestScore);
  	for(int treeno=1; treeno<= noOfTrees; treeno++){

  		int* code = treeNumber2prueferCode(treeno, n);
  		int* parent = prueferCode2parentVector(code, n-1);
  		double newTreeScore = scoreTree(n, m, logScores, dataMatrix, scoreType, parent, bestScore);

        if(newTreeScore > bestScore){                 // case 1: the current tree has better score than previous best tree
            bestScore = newTreeScore;                 // update best score
            printf("new best score: %e\n", newTreeScore);
            emptyTreeList(bestTrees, n);                 // empty the list of best trees (ancestor matrices)
            bestTrees.push_back(parent);           // current tree is now only best tree
        }
        else if(newTreeScore == bestScore){            // case 2: the current tree is equally good as best tree so far
            bestTrees.push_back(parent);   // add it to list (duplicates not possible due to redundancy free enumeration)
        }
        else{                                         // case 3: better tree was found already
            delete [] parent;        // discard tree
        }
        delete [] code;
 	  }
     printf("list size: %lu\n", bestTrees.size());
 	  return bestScore;
}

/* given a tree number (starting from 0), return the corresponding Pruefer code */
int* treeNumber2prueferCode(int treeno, int size){
	int codeLength = size-1;                         // the length of the pruefer code is the treesize-2, since we add the root as node, we have an offset of 1
	int base = size+1;                               // the tree has n+1 nodes, this is the base for constructing the code
	int* code = new int[codeLength];
	int currentno = treeno-1;                       // trees are numbered starting from 1
	for(int i=codeLength-1; i>=0; i--){             // transforms the treeno (base 10) to the pruefer code (base n+1)
	  	code[i] = currentno % base;
	  	currentno = currentno/base;
	 }
	return code;
}









/*bool** codeToAncMatrix_new(int treeno, int size){
	int* code = treeNumber2prueferCode(treeno, size);
	int* parent = prueferCode2parentVector(code, size-1);
	return parentVector2ancMatrix(parent, size);
}*/


/*bool** codeToAncMatrix(int treeno, int size)
{

    bool** ancMatrix = init_boolMatrix(size, size, 0);
	int codeLength = size-1;
  	int code[codeLength];
  	int base = size+1;
  	int startno = treeno-1;
  	for(int i=codeLength-1; i>=0; i--){
  		code[i] = startno % base;
  		startno = startno/base;
  	}

	// get auxiliary arrays
	int lastOcc[size+1];    // node id -> index of last occ in code
	int novelLeaf[size+1];  //   code index of last occ of leaf x -> x, -1 if not a last occ
	bool isLeaf[size+1];

	for(int i=0; i<=size;i++){ lastOcc[i] = -1;	}
	for(int i=0; i<=size;i++){ novelLeaf[i] = -1;	}
	for(int i=0; i<=size;i++){ isLeaf[i] = true;	}

	for(int i=0; i<codeLength; i++){
		//printf("code[%d] = %d\n", i, code[i]);
		lastOcc[code[i]] = i;
		//printf("pos %d, lastOcc[%d] = %d\n", i+1, code[i]+1, lastOcc[code[i]]+1);
	}
	for(int i=0; i<size; i++){
		if(lastOcc[i]>=0){
			novelLeaf[lastOcc[i]] = i;
			//printf("novel leaf[%d] = %d\n", lastOcc[i], novelLeaf[lastOcc[i]]+1);
		}

	}

	for(int i=0; i<=size; i++){
		if(lastOcc[i]!= -1){
			isLeaf[i] =false;
		}
	}

	// fill connectivity matrix
	int curr = 1;                 // currently the smallest old leaf (apart from the novel leaf)
	int novel = size+10;         // node that just became a leaf after previous step
	for(int i=0; i<= size; i++){
		if(isLeaf[i] == true){
			curr = i;             // the smallest leaf in the beginning
			break;
		}
	}
	//printf("smallest leaf: %d\n", curr);
	for(int i=0; i<codeLength; i++){
		if(novel < curr){        // new child node is novel leaf

			if(novel<size && code[i]<size){
				ancMatrix[code[i]][novel] = true;              // set new edge in matrix
			}
			//printf("a: %d -> %d\n", novel+1, code[i]+1);
			for(int j=0; j<size; j++){                  // set new connections
				if(j<size && code[i]<size){
					ancMatrix[code[i]][j] = ancMatrix[novel][j] || ancMatrix[code[i]][j];
				}
				//if(ancMatrix[code[i]][j]==true){
				//	printf("a: [%d,%d] set to 1\n", j+1, code[i]+1);
				//}
			}
			novel = size+10;                           // novel leaf is no longer novel
			if(novelLeaf[i] >= 0){
				isLeaf[novelLeaf[i]] = true;     // update leaf list if there is a new leaf
				novel = novelLeaf[i];            // update novel
				//printf("novel leaf %d", novel+1);
			}

		}
		else{             // next is smallest leaf of the old list
			if(curr<size && code[i]<size){
				if(curr<size && code[i]<size){
					ancMatrix[code[i]][curr] = true;               // set edge in matrix
					//printf("%d -> %d set\n", curr+1, code[i]+1);
				}
			}
			//printf("b: %d -> %d\n", curr+1, code[i]+1);
			for(int j=0; j<size; j++){                  // set new connections
				if(j<size && code[i]<size){
					ancMatrix[code[i]][j] = ancMatrix[curr][j] || ancMatrix[code[i]][j];
				}
				if(ancMatrix[code[i]][j]==true){
					//printf("b: [%d,%d] set to 1        left: %d,%d      right: %d,%d\n", j+1, code[i]+1, j, novel+1, j+1, code[i]+1);
				}
			}
			novel = size+10;
			if(novelLeaf[i] >= 0){
				isLeaf[novelLeaf[i]] = true;     // there is a new leaf
				novel = novelLeaf[i];            // update novel
			}

			do {                                   // update curr by finding
    				curr++;                        // the next smallest leaf
				} while(isLeaf[curr] == false);
		}
	}

		for (int a=0; a<size; a++){
		printf("%i  ", a+1);
		for (int b=0; b<size; b++){
			 ancMatrix[a][b] == false? printf("0 ") : printf("1 ");
		}
		printf("\n");
	}


	// set self-connections
	for (int a=0; a<size; a++){
		ancMatrix[a][a] = true;
		//printf("%d -> %d\n", a+1, a+1);
	}


	//for (int a=0; a<size; a++){
	//	printf("%i  ", a+1);
	//	for (int b=0; b<size; b++){
	//		 ancMatrix[a][b] == false? printf("0 ") : printf("1 ");
	//	}
	//	printf("\n");
	//}

	return ancMatrix;
}*/

//bool** initAncMatrix(int size){
 // bool** ancMatrix = calloc(size, sizeof(bool*));

//    for (int i=0; i<size; ++i)
//    {
//        ancMatrix[i] = calloc(size, sizeof(bool));
 //   }
//     for (int i=0; i<size; ++i)
//    {
 //        for (int j=0; j<size; ++j)
 //   	{
 //       	ancMatrix[i][j] = 0;
 //   	}
 //   }
 //   return ancMatrix;
//}


