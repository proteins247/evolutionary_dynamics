#ifndef __STRUCTURELIB_H__
#define __STRUCTURELIB_H__

#include<stdio.h>
#include<math.h>
#include<stdlib.h>

#define NUMCONF         10000
#define GENES			2

extern int CoordMatrix[NUMCONF][81];
extern int AllFaces[NUMCONF*6*4][9];
extern char ContactMatrixA[NUMCONF][32]; // contact matrix
extern char ContactMatrixB[NUMCONF][32]; 

void ReadCoordMatrix(char *filename);
void GetSquare(int structid, int face, int Square[3][3]);
void RotateSquare(int code, int dest[3][3], int src[3][3]);
void MakeAllFaces(void);
void MakeContactMatrix(void);
void ReadCommondata(char * location, char * energyName);
void PrintFoldedConformation(int structid, char * conformation);
#endif
