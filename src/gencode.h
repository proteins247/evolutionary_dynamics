/*
 * gencode.h
 *
 *  Created on: Dec 25, 2009
 *      Author: mheo
 */

#ifndef GENCODE_H_
#define GENCODE_H_
#define A_STOP -1

// code as in Miyazawa,Jernigan JMB 1996
#define A_Cys 0
#define A_Met 1
#define A_Phe 2
#define A_Ile 3
#define A_Leu 4
#define A_Val 5
#define A_Trp 6
#define A_Tyr 7
#define A_Ala 8
#define A_Gly 9
#define A_Thr 10
#define A_Ser 11
#define A_Asn 12
#define A_Gln 13
#define A_Asp 14
#define A_Glu 15
#define A_His 16
#define A_Arg 17
#define A_Lys 18
#define A_Pro 19

#define N_UUU 0x0000
#define N_UUC 0x0001
#define N_UUA 0x0002
#define N_UUG 0x0003

#define N_UCU 0x0010
#define N_UCC 0x0011
#define N_UCA 0x0012
#define N_UCG 0x0013

#define N_UAU 0x0020
#define N_UAC 0x0021
#define N_UAA 0x0022
#define N_UAG 0x0023

#define N_UGU 0x0030
#define N_UGC 0x0031
#define N_UGA 0x0032
#define N_UGG 0x0033

#define N_CUU 0x0100
#define N_CUC 0x0101
#define N_CUA 0x0102
#define N_CUG 0x0103

#define N_CCU 0x0110
#define N_CCC 0x0111
#define N_CCA 0x0112
#define N_CCG 0x0113

#define N_CAU 0x0120
#define N_CAC 0x0121
#define N_CAA 0x0122
#define N_CAG 0x0123

#define N_CGU 0x0130
#define N_CGC 0x0131
#define N_CGA 0x0132
#define N_CGG 0x0133


#define N_AUU 0x0200
#define N_AUC 0x0201
#define N_AUA 0x0202
#define N_AUG 0x0203

#define N_ACU 0x0210
#define N_ACC 0x0211
#define N_ACA 0x0212
#define N_ACG 0x0213

#define N_AAU 0x0220
#define N_AAC 0x0221
#define N_AAA 0x0222
#define N_AAG 0x0223

#define N_AGU 0x0230
#define N_AGC 0x0231
#define N_AGA 0x0232
#define N_AGG 0x0233

#define N_GUU 0x0300
#define N_GUC 0x0301
#define N_GUA 0x0302
#define N_GUG 0x0303

#define N_GCU 0x0310
#define N_GCC 0x0311
#define N_GCA 0x0312
#define N_GCG 0x0313

#define N_GAU 0x0320
#define N_GAC 0x0321
#define N_GAA 0x0322
#define N_GAG 0x0323

#define N_GGU 0x0330
#define N_GGC 0x0331
#define N_GGA 0x0332
#define N_GGG 0x0333



/*
Routines for manipulating nucleotide sequences
*/

void CreateRandomNucSequence(int *Seq, int Len);
void CreateRandomNucSequence2(int *Seq, int Len);
/* creates a random sequence of int{0,1,2,3} of length Len */

int PointMutateCharNucSequence(char *Seq, int Len);
int PointMutateNucSequence(int *Seq, int Len);
/*
int PointMutateCharNucSequence(char *Seq, int Len);
int PointMutateNucSequence(int *Seq, int Len);
char/int -> always 0,1,2,3 for nucleotides, not UCAG for char!
return values:
0:  synonymous mutation (no changes in aminoacid sequence)
1:  nonsynonymous mutation
-1: mutation to STOP codon
*/

int NucSeqToAASeq(int *NucSeq, int N, int *AASeq);
int CharNucSeqToAASeq(char *NucSeq, int N, int *AASeq);
/*
translates nucleotide sequence to amino acid sequence
nucleotides: 0123 -> UCAG
amino acids: 01234... -> CMFIL... MJ96 order, see latticelib.c
N is the length of nucleotide sequence
return values:
0: conversion successful
1: stop codon encoutered. STOP codon is -1 in AASeq.
*/

void PrintAASequence(char *buf, int *Seq, int Len);
void PrintCharNucCodeSequence(char *buf, char *Seq, int Len);
void PrintNucCodeSequence(char *buf, int *Seq, int Len);
/*input: int* or char* Seq {0,1,2,3}, output char *buf UCAG */

void LetterToNucCodeSeq(char *buf, int *Seq, int Len);
/*input: char *buf UCAG, output int *Seq 0123 */

void CopyIntToCharSeq(char *dest, int *src, int Len);
/*copies Len values from src to dest*/

#endif /* GENCODE_H_ */
