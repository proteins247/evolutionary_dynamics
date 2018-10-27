#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include <string.h>

void InitNegOne(int *list, int len);
int IsMemberInt(int candidate, int *list, int *index);
double Mean(double *data, int col, int index1, int index2);
int Stdev(double *data, int col, int index1, int index2);
double CoeffOfVariation(double *data, int col, int index1, int index2);
//void WriteArrayToFile(double *data, char *fname, int nrows, int ncols);
int CountContacts(int conf, int residue);
int RandomBit();
int AcceptOrRejectAttempt( double diff , double Tsel );
int AcceptOrRejectAttempt2( double diff1 , double diff2 , double Tsel );
