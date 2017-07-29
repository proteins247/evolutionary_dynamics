#ifndef __LATTICELIB_H__
#define __LATTICELIB_H__
/*
latticelib.h 
basic routines for handling 3x3x3 compact 27-mers
(c) 2005 K.Zeldovich kzeldov@fas.harvard.edu
*/
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>

#include "define.h"

// global variables
extern int ContactMatrixLen[NUMCONF]; // actually, useless
extern char ContactMatrixA[NUMCONF][32]; // contact matrix
extern char ContactMatrixB[NUMCONF][32]; 
extern double EnergyMatrix[ALPHABET][ALPHABET];
extern double espectrum[NUMCONF]; // used in Pnat calculation

/*
void ReadContactMatrix(void)
reads the lists of contacts for all 103346 structures
Input:
none, but check the path to contact103346.dat
Output: 
ContactMatrixLen (normally unused)
ContactMatrixA, ContactMatrixB: contacts
*/
void ReadContactMatrix(char *filename);


/*
void ReadEnergyMatrix(void)
Reads amino acid energy interaction matrix. MJ96 = Miyazawa and Jernigan 1996
input: none, but check the file path
output:
EnergyMatrix[][]
*/
void ReadEnergyMatrix(char *filename);


/*
double SequenceEnergy(int *Seq, int i)
returns the energy of a sequence Seq in conformation i,
0<= i < NUMCONF
*/
double SequenceEnergy(int *Seq, int conformnumber);

/*
double GetSequenceMinEnergy(int *Seq, int *atwhatconform)
For a sequence Seq, returns the minimum energy across all conformations.
The conformation at which the minimum is reached is
returned via *atwhatconform
*/
double GetSequenceMinEnergy(int *Seq, int *atwhatconform);

/*
double GetTwoStructureQScore(int StrID1, int StrID2)
calculates the normalized structure similarity (Q-score, 0...1)
between structures StrID1, StrID2,
0<= StrID1, StrID2 < NUMCOUNF,
0<= output <= 1
*/
double GetTwoStructureQScore(int StrID1, int StrID2);

/*
Non-normalized integer structure similarity, 0..28
*/
int GetTwoStructureIntQScore(int StrID1, int StrID2);

/*
void CreateRandomSequence(int *Seq, int Len)
Creates a random (amino acid) sequence Seq of length Len,
each element of Seq is between 0 and ALPHABET-1 (inclusive)
*/
void CreateRandomSequence(int *Seq, int Len);

/* copies src to dest */
void CopySeq(int *dest, int *src, int Len);

/*
void LetterToMJCodeSeq(char *letterseq, int *Seq, int Len)
Converts an amino letter sequence letterseq of length Len 
to the integer representation according to MJ96 energy matrix ordering
Used to read the sequences.
CMFIL... -> 0 1 2 3 4....
*/
void LetterToMJCodeSeq(char *letterseq, int *Seq, int Len);

/*
void PrintAACodeSequence(char *buf, int *Seq, int Len)
buf receives a letter representation of sequence Seq of length Len
Normally used for output
*/
void PrintAACodeSequence(char *buf, int *Seq, int Len);


/*
double GetSequencePnat(int *Seq, double T, int *conform)

Returns the probability of being in the ground state
for the sequence Seq of length AASEQLEN, at temperature T.
The index of conformation of the ground state is returned 
via *conform

Alters the espectrum[] global array.
*/
double GetStructurePnat(int *Seq, double T, int structid);
double GetSequencePnat(int *Seq, double T, int *conform);

/*
void PointMutateAASequence(int *Seq, int Len)
introduces a point mutation in the amino acid sequence Seq of length Len
*/
void PointMutateAASequence(int *Seq, int Len);

/*
void SwapTwoResidues(int *Seq, int Len)
swaps two residues in the amino acid sequence Seq of length Len
*/
void SwapTwoResidues(int *Seq, int Len);


/* two sequence design procedures are below */

/*
double MutateAASeqForMaxPnat(int *aaseq, double Tenv, double Tsel)
Monte Carlo optimization of sequence aaseq to maximize Pnat 
at temperature Tenv. Tsel is the Metropolis MC selection temperature.
output : aaseq = optimized sequence
returns: Pnat of the new sequence
*/
double MutateAASeqForMaxPnat(int *aaseq, double Tenv, double Tsel);

/*
double SwapResiduesForMinE(int *aaseq, double Tsel)
Monte Carlo optimization of sequence aaseq to miminize native state energy
by swapping amino acids (=at fixed concentration)
Tsel is the Metropolis MC selection temperature.
output : aaseq = optimized sequence
returns: E0 of the new sequence
*/
double SwapResiduesForMinE(int *aaseq,  double Tsel);

double GetHydrophobicity(int *aaseq, int currLen);
double GetFracHydrophobicity(int *aaseq, int currLen);
double GetFracHydrophobicity_avil(int *aaseq, int currLen);
double GetNetCharge(int *aaseq, int currLen);
double GetPosCharge(int *aaseq, int currLen);
double GetNegCharge(int *aaseq, int currLen);
double GetFracCharge(int *aaseq, int currLen);

void ConvertAAtoCharge(int *aaseq, int *new_aaseq, int currLen);

void ConvertAAtoAVIL(int *aaseq, int *new_aaseq, int currLen);

void ConvertAAtoHydro(int *aaseq, int *new_aaseq, int currLen);

#endif
