/*************************
 Modified by Orit Peleg
 opeleg@fas.harvard.edu
 on 2/26/13
 ************************/


/* These functions rank all of the species in the OrgDB in order of size so
 * one can later select the most dominant species in the system.  This prevents
 * performing calculations on organisms with low populations that do not
 * contribute much to the system.
*/

#include <stdio.h>
#include <stdlib.h>
#include "cell.h"
#include "compare.h"

int RankSpeciesSizeDB(int rankings[], int spNum)
{
  int i;
  int arraySize[spNum];

  for(i=0;i<spNum;i++) arraySize[i] = myOrgDB[i].count;

//fprintf(stderr, "Ranking now... \n");

  RankSpeciesSize(arraySize, rankings, spNum);

  return 0;
}

int RankSpeciesSize(int arraySize[], int rankings[], int spNum)
{
  int i, j, place;

  rankings[0] = 0;
//fprintf(stderr, "Ready to Sort...\n");
  for (i=1; i<spNum; i++) {
//  fprintf(stderr, "Sorting...\n");
    place = SortSpeciesSize(arraySize, rankings, arraySize[i], 0, i-1);
//  fprintf(stderr, "place before ranking sort: %d \n", place);
    for (j=i; j>place; j--) rankings[j]=rankings[j-1];
    rankings[place] = i;
    
//  for(j=0;j<i+1;j++) fprintf(stderr, "%d ", rankings[j]);
//  fprintf(stderr, "\n");
  }
  return 0;
}

int SortSpeciesSize(int arraySize[], int rankings[], int spSize, int low, int high)
{
  if (low > high) {
    fprintf(stderr, "SortSpeciesSize error: low (%d) > high (%d)\n", low, high);
    return low;
  }

  else if (low == high) {
//  fprintf(stderr, "low: %d, high: %d \n", low, high);
    if (spSize < arraySize[rankings[low]]) return low+1;
    else return low;
  }

  else {
    int mid = (low + high)/2;
//  fprintf(stderr, "low: %d, mid: %d, high: %d \n", low, mid, high);

    if (spSize == arraySize[rankings[mid]]) return mid;
    else if (spSize > arraySize[rankings[mid]]) return SortSpeciesSize(arraySize, rankings, spSize, low, mid);
    else return SortSpeciesSize(arraySize, rankings, spSize, mid+1, high);
  }
}

// This wraps these functions as a standalone program
/*
void main(int argc, char *argv[])
{

  int i;
  int num = argc-1;
  int arraySize[num];
  int rankings[num];

  for(i=0;i<num;i++) sscanf(argv[i+1],"%d",&arraySize[i]);

//  for(i=0;i<num;i++) fprintf(stderr,"%d ", arraySize[i]);
//  fprintf(stderr,"\n");

  RankSpeciesSize(arraySize, rankings, num);

  for(i=0;i<num;i++) fprintf(stdout,"%d ", arraySize[rankings[i]]);
  fprintf(stdout,"\n");

}
*/
