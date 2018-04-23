/*
 * doubleMut.h
 *
 *  Created on: Aug 15, 2015
 *      Author: jahnka
 */

#ifndef RECMUT_H_
#define RECMUT_H_



double scoreTree_recMut(int n, int m, double** logScores, int** dataMatrix, char type, int* parentVector, double bestTreeLogScore, int recMut);
double scoreTreeFast_recMut(int n, int m, double** logScores, int** dataMatrix, char type, int* parentVector, int recMut);
double scoreTreeAccurate_recMut(int n, int m, double** logScores, int** dataMatrix, char type, int* parentVector, int recMut);
double maxScoreTreeFast_recMut(int n, int m, double** logScores, int** dataMatrix, int* parent, int* bft, int recMut, int backMut);
double sumScoreTreeFast_recMut(int n, int m, double** logScores, int** dataMatrix, int* parent, int* bft, int recMut, int backMut);
double maxScoreTreeAccurate_recMut(int n, int m, double** logScores, int** dataMatrix, int* parent, int* bft, int recMut, int backMut);
double sumScoreTreeAccurate_recMut(int n, int m, double** logScores, int** dataMatrix, int* parent, int* bft, int recMut, int backMut);
double getSumAttachmentScoreAccurate_recMut(int* parent, int n, double** logScores, int* obsMutState, int* bft, int recMut, int backMut);
int getBackmutation(int* parent, int recMut, int n);
double* getAttachmentScoresFast_recMut(int*parent, int n, double** logScores, int* obsMutProfile, int*bft, int recMut, int backMut);
int** getBestAttachmentScoreAccurate_recMut(int* parent, int n, double** logScores, int* obsMutState, int* bft, int recMut, int backMut);
int*** getAttachmentMatrices_recMut(int* parent, int n, int* obsMutState, int* bft, int recMut, int backMut);


#endif /* RECMUT_H_ */
