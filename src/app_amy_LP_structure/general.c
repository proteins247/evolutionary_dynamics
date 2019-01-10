#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include <string.h>
#include<time.h>
#include<assert.h>

#include "structurelib.h"
#include "../rng.h"

void InitNegOne(int *list, int len)
{
   int i;
   for (i=0;i<len;i++) {list[i]=-1;}
   return;
}

int IsMemberInt(int candidate, int *list, int *index)
{
   int i=0, ismem;
   while (list[i] != -1)
   {
      if (candidate==list[i]) { ismem = 1; *index = -1; break; }
      else {ismem = 0;}
      i++;
   }
   *index = i;
   return ismem;
}

/* Find the mean of the values in data[index1:index2][col]
 * note that the range [index1,index2) is covered as in python
 */
double Mean(double *data, int col, int index1, int index2)
{
   double mean = 0.0;
   int ii;
   for(ii = index1; ii < index2; ii++) {mean += data[ii];}
   mean /= ( index2 - index1 );
   return mean;
}

int Stdev(double *data, int col, int index1, int index2)
{
   double mean, stdev = 0.0;
   int ii;
   mean = Mean(data, col, index1, index2);
   for(ii = index1; ii < index2; ii++) { stdev += pow( data[ii] - mean ,2 ); }
   stdev /= ( index2 - index1 );
   stdev = pow( stdev, 0.5 );
   return stdev;
}

double CoeffOfVariation(double *data, int col, int index1, int index2)
{
   double mean, stdev = 0.0, coeff;
   int ii;

   /*
   // clock
   clock_t start, stop;
   double t=0.0;
   assert((start = clock())!=-1);
   */

   mean = Mean(data, col, index1, index2);
   for(ii = index1; ii < index2; ii++) { stdev += pow( data[ii] - mean ,2 );}
   stdev /= ( index2 - index1 );
   stdev = pow( stdev, 0.5 );
   coeff = stdev / mean;
   
   /*
   // clock
   stop = clock();
   t = (double) (stop-start)/CLOCKS_PER_SEC;
   printf("%0.3G\t", t);
   */

   return coeff;
}

/*
void WriteArrayToFile(double *data, char *fname, int nrows, int ncols)
{
   FILE *myfile;
   int ii, jj;

   myfile = fopen(fname, "w");

   for( ii = 0; ii < nrows; ii++)
   {
      for( jj = 0; jj < ncols; jj++)
      {
         fprintf(myfile, "%0.4g\t", data[ii][jj]);
      } // end writting cols
      fprintf(myfile, "\n");
   } // end writing rows
   
   fclose(myfile);
   return;
}
*/

// Determines the number of contacts made by a given residue in a particular conformation
int CountContacts(int conf, int residue)
{
int ncontacts=0, ii;
for (ii=0; ii<28; ii++)
{
if (ContactMatrixA[conf][ii]==residue) ncontacts++;
else if (ContactMatrixB[conf][ii]==residue) ncontacts++;
}
return ncontacts;
}

/*
 * return a random bit (0 or 1)
 */
int RandomBit()
{
   int c;
   double t;
   t = threefryrand();
   if (t<0.5) {c=0;}
   else {c=1;}
   //printf("%f\t", t);
   return c;
}

/*
 * Decides whether to accept or reject a deleterious mutation based on MC criterion
 */
int AcceptOrRejectAttempt(double diff , double Tsel )
{
   double r, crit;
   int accept;
   
   r = threefryrand();
   crit = exp(-( diff )/Tsel);
   //printf("%f\t%f\t%f\t", diff, r, crit);

   if (r < crit ) { accept = 1; }
   else { accept = 0; }

   //printf("%f %f ", crit, r);
   return accept;
}

int AcceptOrRejectAttempt2(double diff1, double diff2, double Tsel )
{
   double r, crit1, crit2;
   int accept;
   
   r = threefryrand();
   crit1 = exp(-( diff1)/Tsel);
   crit2 = exp(-( diff2)/Tsel);
   //printf("%f\t%f\t%f\t", diff, r, crit);

   if (r < crit1 && r < crit2) { accept = 1; }
   else { accept = 0; }

   //printf("%f %f ", crit, r);
   return accept;
}
