/*
 * gencode.c
 *
 *  Genetic code and nucleotide sequences routines
 *  (c) 2005 K.Zeldovich
 *  Created on: Dec 25, 2009
 *      Author: mheo
 */

#include <stdio.h>
#include <stdlib.h>
#include "gencode.h"
#include "rng.h"
// U C A G -> 0 1 2 3


/* 
 * Convert from enum Codon to enum AminoAcid.  
 */
static AminoAcid TripletToAA(Codon triplet)
{
    int out = -2;

    switch(triplet)
    {

    case N_UUU: { out =  A_Phe; break; }
    case N_UUC: { out =  A_Phe; break; }
    case N_UUA: { out =  A_Leu; break; }
    case N_UUG: { out =  A_Leu; break; }

    case N_UCU: { out =  A_Ser; break; }
    case N_UCC: { out =  A_Ser; break; }
    case N_UCA: { out =  A_Ser; break; }
    case N_UCG: { out =  A_Ser; break; }

    case N_UAU: { out =  A_Tyr; break; }
    case N_UAC: { out =  A_Tyr; break; }
    case N_UAA: { out =  A_STOP; break; }
    case N_UAG: { out =  A_STOP; break; }

    case N_UGU: { out =  A_Cys; break; }
    case N_UGC: { out =  A_Cys; break; }
    case N_UGA: { out =  A_STOP; break; }
    case N_UGG: { out =  A_Trp; break; }

//-----------

    case N_CUU: { out =  A_Leu; break; }
    case N_CUC: { out =  A_Leu; break; }
    case N_CUA: { out =  A_Leu; break; }
    case N_CUG: { out =  A_Leu; break; }

    case N_CCU: { out =  A_Pro; break; }
    case N_CCC: { out =  A_Pro; break; }
    case N_CCA: { out =  A_Pro; break; }
    case N_CCG: { out =  A_Pro; break; }

    case N_CAU: { out =  A_His; break; }
    case N_CAC: { out =  A_His; break; }
    case N_CAA: { out =  A_Gln; break; }
    case N_CAG: { out =  A_Gln; break; }

    case N_CGU: { out =  A_Arg; break; }
    case N_CGC: { out =  A_Arg; break; }
    case N_CGA: { out =  A_Arg; break; }
    case N_CGG: { out =  A_Arg; break; }

//-----------

    case N_AUU: { out =  A_Ile; break; }
    case N_AUC: { out =  A_Ile; break; }
    case N_AUA: { out =  A_Ile; break; }
    case N_AUG: { out =  A_Met; break; }

    case N_ACU: { out =  A_Thr; break; }
    case N_ACC: { out =  A_Thr; break; }
    case N_ACA: { out =  A_Thr; break; }
    case N_ACG: { out =  A_Thr; break; }

    case N_AAU: { out =  A_Asn; break; }
    case N_AAC: { out =  A_Asn; break; }
    case N_AAA: { out =  A_Lys; break; }
    case N_AAG: { out =  A_Lys; break; }

    case N_AGU: { out =  A_Ser; break; }
    case N_AGC: { out =  A_Ser; break; }
    case N_AGA: { out =  A_Arg; break; }
    case N_AGG: { out =  A_Arg; break; }

//-----------

    case N_GUU: { out =  A_Val; break; }
    case N_GUC: { out =  A_Val; break; }
    case N_GUA: { out =  A_Val; break; }
    case N_GUG: { out =  A_Val; break; }

    case N_GCU: { out =  A_Ala; break; }
    case N_GCC: { out =  A_Ala; break; }
    case N_GCA: { out =  A_Ala; break; }
    case N_GCG: { out =  A_Ala; break; }

    case N_GAU: { out =  A_Asp; break; }
    case N_GAC: { out =  A_Asp; break; }
    case N_GAA: { out =  A_Glu; break; }
    case N_GAG: { out =  A_Glu; break; }

    case N_GGU: { out =  A_Gly; break; }
    case N_GGC: { out =  A_Gly; break; }
    case N_GGA: { out =  A_Gly; break; }
    case N_GGG: { out =  A_Gly; break; }

    }

//printf("dec %.3x to %d\n",triplet,out);

    if (out == -2)
    {
        fprintf(stderr,"genetic code error\n");
        fprintf(stderr,"triplet: %.3x\n",triplet);
        exit(-1);
    }

    return (AminoAcid)out;
}


/*
 * Convert from enum AminoAcid to enum Codon
 * if (random)
 *   select randomly among corresponding codons
 * else
 *   select most common codon in E. coli for that amino acid
 *
 */
static Codon AAtoCodon(AminoAcid aa, int random)
{
    int end = -1;
    int codons[7] = {end, end, end, end, end, end, end};      /* There are max six codons for aa */
    int num_codons = 0;

    /* codons[0] will be the most common codon */
    switch (aa)
    {
    case A_Cys: {
        codons[1] = N_UGU;
        codons[0] = N_UGC;
        break; }
    case A_Met: {
        codons[0] = N_AUG;
        break; }
    case A_Phe: {
        codons[0] = N_UUU;
        codons[1] = N_UUC;
        break; }
    case A_Ile: {
        codons[0] = N_AUU;
        codons[1] = N_AUC;
        codons[2] = N_AUA;
        break; }
    case A_Leu: {
        codons[1] = N_UUA;
        codons[2] = N_UUG;
        codons[3] = N_CUU;
        codons[4] = N_CUC;
        codons[5] = N_CUA;
        codons[0] = N_CUG;
        break; }
    case A_Val: {
        codons[1] = N_GUU;
        codons[2] = N_GUC;
        codons[3] = N_GUA;
        codons[0] = N_GUG;
        break; }
    case A_Trp: {
        codons[0] = N_UGG;
        break; }
    case A_Tyr: {
        codons[0] = N_UAU;
        codons[1] = N_UAC;
        break; }
    case A_Ala: {
        codons[0] = N_GCU;
        codons[1] = N_GCC;
        codons[2] = N_GCA;
        codons[3] = N_GCG;
        break; }
    case A_Gly: {
        codons[0] = N_GGU;
        codons[1] = N_GGC;
        codons[2] = N_GGA;
        codons[3] = N_GGG;
        break; }
    case A_Thr: {
        codons[0] = N_ACU;
        codons[1] = N_ACC;
        codons[2] = N_ACA;
        codons[3] = N_ACG;
        break; }
    case A_Ser: {
        codons[0] = N_UCU;
        codons[1] = N_UCC;
        codons[2] = N_UCA;
        codons[3] = N_UCG;
        codons[4] = N_AGU;
        codons[5] = N_AGC;
        break; }
    case A_Asn: {
        codons[0] = N_AAU;
        codons[1] = N_AAC;
        break; }
    case A_Gln: {
        codons[0] = N_CAA;
        codons[1] = N_CAG;
        break; }
    case A_Asp: {
        codons[0] = N_GAU;
        codons[1] = N_GAC;
        break; }
    case A_Glu: {
        codons[0] = N_GAA;
        codons[1] = N_GAG;
        break; }
    case A_His: {
        codons[0] = N_CAU;
        codons[1] = N_CAC;
        break; }
    case A_Arg: {
        codons[0] = N_CGU;
        codons[1] = N_CGC;
        codons[2] = N_CGA;
        codons[3] = N_CGG;
        codons[4] = N_AGA;
        codons[5] = N_AGG;
        break; }
    case A_Lys: {
        codons[0] = N_AAA;
        codons[1] = N_AAG;
        break; }
    case A_Pro: {
        codons[0] = N_CCU;
        codons[1] = N_CCC;
        codons[2] = N_CCA;
        codons[3] = N_CCG;
        break; }
    case A_STOP: {
        codons[0] = N_UAA;
        codons[1] = N_UAG;
        codons[2] = N_UGA;
        break; }
    default: {
        fprintf(stderr, "AAtoCodon error\n");
        fprintf(stderr, "aa value: %d\n", aa);
        exit(-1);    }
    }

    /* count number of codons for our aa */
    while (codons[num_codons] != end)
    {
        num_codons++;
    }

    if (random)
    {
        return (Codon)codons[(int)(threefryrand() * num_codons)];
    }
    else
    {
        return (Codon)codons[0];
    }

}


/* functions in gencode.h: */

void CreateRandomNucSequence(int *NucSeq, int Len)
{
    int i, nuc;
    for (i=0; i<Len; i++)
    {
        nuc = (int)( threefryrand() * 4 );
        NucSeq[i] = nuc;
    }
    return;
}



/* 
 * CreateRandomNucSequence2 version avoids stop codons. Len should be
 * the size of Seq and divisible by 3. (this is not checked)
 * 
 */
void CreateRandomNucSequence2(int *NucSeq, int Len)
{
    Codon t;
    AminoAcid aa;
    unsigned char n1, n2, n3;
    int i, j, nuc;
    for (i=0; i<Len/3; i++)
    {
        do
        {
            for (j=0; j<3 ;j++)
            {
                nuc = (int)(threefryrand() * 4 );
                NucSeq[3*i+j] = nuc;
            }
            n1 = NucSeq[3*i];
            n2 = NucSeq[3*i+1];
            n3 = NucSeq[3*i+2];
            t = 0;
            t = (n1 << 8 ) | (n2 << 4) | n3;
            aa = TripletToAA(t);
            //printf(" %x %x %x -> 0x%.4x\n",a1,a2,a3, a);
        } while (aa == A_STOP);
    }
    return;
}


int NucSeqToAASeq(int *NucSeq, int N, AminoAcid *AASeq)
{
    Codon t;
    unsigned char n1, n2, n3;
    int k=0, i, j=0;

    for (i=0; i<N; i+=3)
    {
        n1 = NucSeq[i];
        n2 = NucSeq[i+1];
        n3 = NucSeq[i+2];

        t = 0;
        t = (n1 << 8 ) | (n2 << 4) | n3;

//printf(" %x %x %x -> 0x%.4x\n",a1,a2,a3, a);

        AASeq[j] = TripletToAA(t);
        if (AASeq[j] == A_STOP)
        {
//fprintf(stderr,"%.3x to stop\n",a);*/
            k = -1;
        }
        if (AASeq[j] == -2)
        {
            fprintf(stderr,"NucSeqToAASeq() : Genetic Code Error!!!\n");
            exit(1);
        }
        j++;
    }

    return k;
}


/* This function only differs from above in that NucSeqToAASeq is char * */
int CharNucSeqToAASeq(char *NucSeq, int N, AminoAcid *AASeq)
{
    Codon t;
    unsigned char n1, n2, n3;
    int k=0, i, j=0;

    for (i=0; i<N; i+=3)
    {
        n1 = NucSeq[i];
        n2 = NucSeq[i+1];
        n3 = NucSeq[i+2];

        t = 0;
        t = (n1 << 8 ) | (n2 << 4) | n3;

//printf(" %x %x %x -> 0x%.4x\n",a1,a2,a3, a);

        AASeq[j] = TripletToAA(t);
        if (AASeq[j] == A_STOP)
        {
//fprintf(stderr,"%.3x to stop\n",a);*/
            k = -1;
        }
        if (AASeq[j] == -2)
        {
            fprintf(stderr,"CharNucSeqToAASeq() : Genetic Code Error!!!\n");
            exit(1);
        }
        j++;
    }

    return k;
}


void AASeqToNucSeq(AminoAcid * AASeq, int * NucSeq, int AALen, int random)
{
    Codon c;
    unsigned char n1, n2, n3;
    int i, j;

    for (i=0; i<AALen; i++)
    {
        j = i * 3;
        /* enum AminoAcid -> enum Codon */
        c = AAtoCodon(AASeq[i], random);
        /* Mask and bitshift to separate the three values */
        n1 = (c & 0x0f00) >> 8;
        n2 = (c & 0x00f0) >> 4;
        n3 = (c & 0x000f);
        NucSeq[j] = n1;
        NucSeq[j+1] = n2;
        NucSeq[j+2] = n3;
    }
}


int PointMutateNucSequence(int *NucSeq, int Len)
{
    int i, j, k;
    AminoAcid AAS1[10000], AAS2[10000];

    NucSeqToAASeq(NucSeq, Len, AAS1);

    /* Select a nucleotide */
    j = (int)( threefryrand() * Len );

    /* Keep on going until nucleotide is changed */
    do
    {
        i = (int)( threefryrand() * 4 );
    } while (NucSeq[j] == i);

    NucSeq[j] = i;

    /* Get the new AA sequence */
    k = NucSeqToAASeq(NucSeq, Len, AAS2);

    if (k == -1)                 /* mutation to stop codon */
        return -1;

    j = 0;
    for (i=0; i<Len/3; i++)
        if (AAS1[i] != AAS2[i])
        {
            j = 1;
            break;
        }

    return j;
}

/* This function only differs from PointMutateNucSequence in that Seq is char * */
int PointMutateCharNucSequence(char *NucSeq, int Len)
{
    int i, j, k;
    AminoAcid AAS1[10000], AAS2[10000];

    CharNucSeqToAASeq(NucSeq, Len, AAS1);

    /* Select a nucleotide */
    j = (int)( threefryrand() * Len );

    do
    {
        i = (int)( threefryrand() * 4 );
    } while (NucSeq[j] == i);

    NucSeq[j] = i;

    /* Get the new AA sequence */
    k = CharNucSeqToAASeq(NucSeq, Len, AAS2);

    if (k == -1)
        return -1;

    j = 0;
    for (i=0; i<Len/3; i++)
        if (AAS1[i] != AAS2[i])
        {
            j = 1;
            break;
        }

    return j;
}


void PrintAASequence(char *buf, AminoAcid *Seq, int Len)
{
    int i;
    char c;
    for(i=0;i<Len;i++)
    {
        c = 'Z';
        switch (Seq[i])
        {
        case A_Cys: { c = 'C'; break; }
        case A_Met: { c = 'M'; break; }
        case A_Phe: { c = 'F'; break; }
        case A_Ile: { c = 'I'; break; }
        case A_Leu: { c = 'L'; break; }
        case A_Val: { c = 'V'; break; }
        case A_Trp: { c = 'W'; break; }
        case A_Tyr: { c = 'Y'; break; }
        case A_Ala: { c = 'A'; break; }
        case A_Gly: { c = 'G'; break; }
        case A_Thr: { c = 'T'; break; }
        case A_Ser: { c = 'S'; break; }
        case A_Asn: { c = 'N'; break; }
        case A_Gln: { c = 'Q'; break; }
        case A_Asp: { c = 'D'; break; }
        case A_Glu: { c = 'E'; break; }
        case A_His: { c = 'H'; break; }
        case A_Arg: { c = 'R'; break; }
        case A_Lys: { c = 'K'; break; }
        case A_Pro: { c = 'P'; break; }
        case A_STOP: { c = '\0'; break; }
        default: {
            fprintf(stderr, "PrintAAsequence error: sequence val: %d", Seq[i]);
            exit(1);
        }

        }
        buf[i] = c;
    }
    buf[i] = 0;
    return;
}


/* recall a nucleotide sequence is composed of ints 1,2,3,4 */
void PrintNucCodeSequence(char *buf, int *NucSeq, int Len)
{
    int i;
    char c;
    for (i=0; i<Len; i++)
    {
        c = 'Z';
        switch (NucSeq[i])
        {
        case  0: { c = 'U'; break; }
        case  1: { c = 'C'; break; }
        case  2: { c = 'A'; break; }
        case  3: { c = 'G'; break; }
        default: {
            fprintf(stderr, "PrintNucCodeSequence error: untranslatable: %d", NucSeq[i]);
            exit(1);
        }
        }
        buf[i] = c;
    }
    buf[i] = 0;
    return;
}

void PrintCharNucCodeSequence(char *buf, char *NucSeq, int Len)
{
    int i;
    char c;
    for (i=0; i<Len; i++)
    {
        c = 'Z';
        switch (NucSeq[i])
        {
        case  0: { c = 'U'; break; }
        case  1: { c = 'C'; break; }
        case  2: { c = 'A'; break; }
        case  3: { c = 'G'; break; }
        default: {
            fprintf(stderr, "PrintNucCodeSequence error: untranslatable: %d", NucSeq[i]);
            exit(1);        }
            
        }
        buf[i] = c;
    }
    buf[i] = 0;
    return;
}

void LetterToNucCodeSeq(char *buf, int *NucSeq, int Len)
{
    int i;
    int c;
    for (i=0; i<Len; i++)
    {
        c = -1;
        switch (buf[i])
        {
        case  'U': { c = 0; break; }
        case  'C': { c = 1; break; }
        case  'A': { c = 2; break; }
        case  'G': { c = 3; break; }
        }
        NucSeq[i] = c;

    }
    return;
}


void LetterToAASeq(char *buf, AminoAcid *Seq, int Len)
{
    int i;
    AminoAcid a;
    for (i=0; i<Len; i++)
    {
        switch (buf[i])
        {
        case 'C': { a = A_Cys; break; }
        case 'M': { a = A_Met; break; }
        case 'F': { a = A_Phe; break; }
        case 'I': { a = A_Ile; break; }
        case 'L': { a = A_Leu; break; }
        case 'V': { a = A_Val; break; }
        case 'W': { a = A_Trp; break; }
        case 'Y': { a = A_Tyr; break; }
        case 'A': { a = A_Ala; break; }
        case 'G': { a = A_Gly; break; }
        case 'T': { a = A_Thr; break; }
        case 'S': { a = A_Ser; break; }
        case 'N': { a = A_Asn; break; }
        case 'Q': { a = A_Gln; break; }
        case 'D': { a = A_Asp; break; }
        case 'E': { a = A_Glu; break; }
        case 'H': { a = A_His; break; }
        case 'R': { a = A_Arg; break; }
        case 'K': { a = A_Lys; break; }
        case 'P': { a = A_Pro; break; }
        default: {
            fprintf(stderr, "LetterToAASeq error: sequence char: %c", buf[i]);
            exit(1); }
        }
        Seq[i] = a;
    }
}


void CopyIntToCharSeq(char *dest, int *src, int Len)
{
    int i;
    for (i=0; i<Len; i++)
        dest[i] = src[i];
    return;
}


void ReadTranslationTimes(char * filename, unsigned int * times)
{
    FILE *infile;
    char reading_buffer[100];
    int line_num = 0;
    Codon c;
    unsigned int time;

    infile = fopen(filename, "r");
    if (!infile)
    {
        fprintf(stderr,
                "failed to open translation speed file: %s\n", filename);
        exit(1);
    }

    /* fscanf, getline? */
    while (fgets(reading_buffer, 100, infile) != NULL)
    {
        /* we ignore all lines that don't fit the format */
        if (sscanf(reading_buffer, "%x %u", &c, &time) == 2)
        {
            times[c] = time;
            line_num++;
        }
    }
    /* there are no error checks, but ideally line_num == 61, for the
       61 non-stop codons */
}
