/*
 * enum.h
 *
 *  Created on: Mar 27, 2015
 *      Author: jahnka
 */
#ifndef ENUM_H
#define ENUM_H

double findBestTreeByExhaustiveSearch(int n, int m, std::vector<bool**> &bestTrees, double** logScores, int** dataMatrix, char scoreType);
//bool** codeToAncMatrix(int treeno, int size);

int* treeNumber2prueferCode(int treeno, int size);
//int* prueferCode2parentVector(int* code, int codeLength);
//bool** parentVector2adjMatrix(int* parent, int length);
//bool** codeToAncMatrix_new(int treeno, int size);

//bool* getInitialQueue(int* code, int codeLength);
//int* getLastOcc(int* code, int codeLength);
//int getNextInQueue(bool* queue, int pos, int length);
//void updateQueue(int node, bool* queue, int next);
//int updateQueueCutter(int node, bool* queue, int next);

#endif
