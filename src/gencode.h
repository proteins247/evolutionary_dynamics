/*
 * gencode.h
 *
 *  Created on: Dec 25, 2009
 *      Author: mheo
 */

#ifndef GENCODE_H_
#define GENCODE_H_

/* code as in Miyazawa,Jernigan JMB 1996 */
typedef enum AminoAcid
{
    A_STOP = -1,
    A_Cys = 0,
    A_Met = 1,
    A_Phe = 2,
    A_Ile = 3,
    A_Leu = 4,
    A_Val = 5,
    A_Trp = 6,
    A_Tyr = 7,
    A_Ala = 8,
    A_Gly = 9,
    A_Thr = 10,
    A_Ser = 11,
    A_Asn = 12,
    A_Gln = 13,
    A_Asp = 14,
    A_Glu = 15,
    A_His = 16,
    A_Arg = 17,
    A_Lys = 18,
    A_Pro = 19
} AminoAcid;

typedef enum Codon
{
    N_UUU = 0x0000,
    N_UUC = 0x0001,
    N_UUA = 0x0002,
    N_UUG = 0x0003,

    N_UCU = 0x0010,
    N_UCC = 0x0011,
    N_UCA = 0x0012,
    N_UCG = 0x0013,

    N_UAU = 0x0020,
    N_UAC = 0x0021,
    N_UAA = 0x0022,
    N_UAG = 0x0023,

    N_UGU = 0x0030,
    N_UGC = 0x0031,
    N_UGA = 0x0032,
    N_UGG = 0x0033,

    N_CUU = 0x0100,
    N_CUC = 0x0101,
    N_CUA = 0x0102,
    N_CUG = 0x0103,

    N_CCU = 0x0110,
    N_CCC = 0x0111,
    N_CCA = 0x0112,
    N_CCG = 0x0113,

    N_CAU = 0x0120,
    N_CAC = 0x0121,
    N_CAA = 0x0122,
    N_CAG = 0x0123,

    N_CGU = 0x0130,
    N_CGC = 0x0131,
    N_CGA = 0x0132,
    N_CGG = 0x0133,


    N_AUU = 0x0200,
    N_AUC = 0x0201,
    N_AUA = 0x0202,
    N_AUG = 0x0203,

    N_ACU = 0x0210,
    N_ACC = 0x0211,
    N_ACA = 0x0212,
    N_ACG = 0x0213,

    N_AAU = 0x0220,
    N_AAC = 0x0221,
    N_AAA = 0x0222,
    N_AAG = 0x0223,

    N_AGU = 0x0230,
    N_AGC = 0x0231,
    N_AGA = 0x0232,
    N_AGG = 0x0233,

    N_GUU = 0x0300,
    N_GUC = 0x0301,
    N_GUA = 0x0302,
    N_GUG = 0x0303,

    N_GCU = 0x0310,
    N_GCC = 0x0311,
    N_GCA = 0x0312,
    N_GCG = 0x0313,

    N_GAU = 0x0320,
    N_GAC = 0x0321,
    N_GAA = 0x0322,
    N_GAG = 0x0323,

    N_GGU = 0x0330,
    N_GGC = 0x0331,
    N_GGA = 0x0332,
    N_GGG = 0x0333
} Codon;


/*
 * Routines for manipulating nucleotide sequences
 */

/* 
 * Create a random sequence of int{0,1,2,3} of length Len 
 * CreateRandomNucSequence2 avoids stop codons
 */
void CreateRandomNucSequence(int *NucSeq, int Len);
void CreateRandomNucSequence2(int *NucSeq, int Len);

/*
 * translate nucleotide sequence to amino acid sequence
 * nucleotides: 0123 -> UCAG
 * amino acids: 01234... -> CMFIL... MJ96 order, see latticelib.c
 * N is the length of nucleotide sequence
 * return values:
 * 0: conversion successful
 * -1: stop codon encoutered. STOP codon is -1 in AASeq.
 */
int NucSeqToAASeq(int *NucSeq, int N, AminoAcid *AASeq);
int CharNucSeqToAASeq(char *NucSeq, int N, AminoAcid *AASeq);

/* 
 * Conversion from nucleotide sequence to codon sequences
 */
void NucSeqToCodonSeq(int *NucSeq, int NucLen, Codon *CodonSeq);

/*
 * AASeqToNucSeq
 *
 * Reverse translate aa sequence to nucleotide sequence.  The genetic
 * code is degenerate; AASeqToNucSeq selects the most frequent codon
 * in E. coli if param `random` is 0 (false), else AASeqToNucSeq picks
 * randomly among cognate codons.
 *
 * Note that AASeq should be an array of AminoAcid (i.e. ints). Use
 * LetterToAASeq to convert from single-letter sequences to int
 * sequence.
 */
void AASeqToNucSeq(AminoAcid * AASeq, int * NucSeq, int AALen, int random);

/*
 * int PointMutateCharNucSequence(char *Seq, int Len);
 * int PointMutateNucSequence(int *Seq, int Len);
 * char/int -> always 0,1,2,3 for nucleotides, not UCAG for char!
 * return values:
 * 0:  synonymous mutation (no changes in aminoacid sequence)
 * 1:  nonsynonymous mutation
 * -1: mutation to STOP codon
 */
int PointMutateNucSequence(int *NucSeq, int Len);
int PointMutateCharNucSequence(char *NucSeq, int Len);

/* Printing and data conversion routines: */
/* input: int* or char* Seq {0,1,2,3}, output char *buf UCAG */
/* returns 0 if successful, -1 upon encountering invalid letter */
void PrintAASequence(char *buf, AminoAcid *Seq, int Len);
void PrintCharNucCodeSequence(char *buf, char *NucSeq, int Len);
void PrintNucCodeSequence(char *buf, int *NucSeq, int Len);

/*input: char *buf UCAG, output int *Seq 0123 */
int LetterToNucCodeSeq(const char *buf, int *NucSeq, int Len);

/* input: amino acid char *buf ACDEF..., output AminoAcid *Seq 0123 
 * returns 0 if successful, -1 upon encountering invalid letter
 * new function added -- VZ
 */
int LetterToAASeq(const char *buf, AminoAcid *Seq, int Len);

/*copies Len values from src to dest*/
void CopyIntToCharSeq(char *dest, int *src, int Len);


/* 
 * Additional functions to handle the relationship between codons and
 * translation speed.
 */

/*
 * ReadTranslationTimes opens file at filename and reads in
 * translation times to an int array.
 *
 * times should hold at least 0x333 + 1 ints s.t. random access using
 * enum Codon values is possible. Memory is cheap and 820 ints is only
 * 3280 bytes
 *
 * The unit for these times depends on the user's application
 */
void ReadTranslationTimes(const char * filename, unsigned int * times);

#endif /* GENCODE_H_ */
