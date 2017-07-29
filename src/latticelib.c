/*
 latticelib.c
 basic routines for handling 3x3x3 compact 27-mers
 (c) 2005 K.Zeldovich kzeldov@fas.harvard.edu
 */
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>

#include"latticelib.h"
#include"../RNG/generator.h"

// global variables
int ContactMatrixLen[NUMCONF]; // actually, useless
char ContactMatrixA[NUMCONF][32]; // contact matrix
char ContactMatrixB[NUMCONF][32];
double EnergyMatrix[ALPHABET][ALPHABET];
double Hydrophobicity[ALPHABET];
double HydrophobicityYesNo[ALPHABET];
double HydrophobicityYesNo_avil[ALPHABET];
double Charge[ALPHABET];
double espectrum[NUMCONF]; // used in Pnat calcultion

//extern double state_3rd;

//extern FILE *globaloutfile;
//extern double targetpnat;

/*
 void ReadContactMatrix(void)
 reads the lists of contacts for all 103346 structures
 Input:
 none, but check the path to contact103346.dat
 Output:
 ContactMatrixLen (normally unused)
 ContactMatrixA, ContactMatrixB: contacts
 */
void ReadContactMatrix(char *filename)
{
    FILE *infile;
    int n, a, b, nprev, i, j, maxconform;
    
    infile=fopen(filename,"r");
    if (!infile) { fprintf(stderr, "no contact matrix\n"); exit(1); }
    i=0; j=0; nprev=0;
    while(!feof(infile))
    {
        fscanf(infile,"%d %d %d\n",&n, &a, &b);
        if (n!=nprev) {nprev=n; i++; j=0; }
        
        if (a>b) // storing only the upper part
        {
            ContactMatrixLen[i]=ContactMatrixLen[i]+1;
            ContactMatrixA[i][j] = a;
            ContactMatrixB[i][j] = b;
            j++;
        }
        
    }
    maxconform = i+1;
    
    fclose(infile);
    
    fprintf(stderr, "%d conformations read\n",maxconform);
    //fprintf(stderr,"contvectorlen %d\n", ContactMatrixLen[100]);
    
    return;
} //ReadContactMatrix


/*
 void ReadEnergyMatrix(void)
 Reads amino acid energy interaction matrix. MJ96 = Miyazawa and Jernigan 1996
 input: none, but check the file path
 output:
 EnergyMatrix[][]
 */
void ReadEnergyMatrix(char *filename)
{
    FILE *infile;
    int i,j, tmp1, tmp2;
    double tmp3;
    
    infile = fopen(filename,"r");
    if (!infile) { fprintf(stderr, "no energy matrix\n"); exit(1); }
    for(i=0;i<ALPHABET;i++)
    {
        for(j=0;j<ALPHABET;j++)
        {
            fscanf(infile,"%d %d %lf\n",&tmp1, &tmp2, &tmp3);
            EnergyMatrix[tmp1][tmp2]=tmp3+3.16;	//mean of MJ96 is -3.16
        }
    }
    fclose(infile);
    
    /*
     // debug dump
     for(i=0;i<ALPHABET;i++)for(j=0;j<ALPHABET;j++) printf("%d\t%d\t%.2f\n",i,j, EnergyMatrix[i][j]);
     */
    
    fprintf(stderr,"EnergyMatrix reads!\n");
    
    return;
} //ReadEnergyMatrix


/*
 double SequenceEnergy(int *Seq, int i)
 returns the energy of a sequence Seq in conformation i,
 0<= i < NUMCONF
 */
double SequenceEnergy(int *Seq, int conformnumber)
{
    int i,j,k;
    double e=0.;
    
    // here we imply that it is a 27-mer!!
    // in general, go nested 1..AASEQLEN and lookup ContMatrix
    
    //for(k=0; k<ContactMatrixLen[conformnumber]; k++)
    for(k=0; k<HALFCONTMATRLEN; k++)
    {
        i =  ContactMatrixA[conformnumber][k];
        j =  ContactMatrixB[conformnumber][k];
        e += EnergyMatrix[Seq[i]][Seq[j]];
    }
    return e;
}

/*
 double GetSequenceMinEnergy(int *Seq, int *atwhatconform)
 For a sequence Seq, returns the minimum energy across all conformations.
 The conformation at which the minimum is reached is
 returned via *atwhatconform
 */
double GetSequenceMinEnergy(int *Seq, int *atwhatconform)
{
    int i,j;
    double e,emin=0.;
    emin = SequenceEnergy(Seq, 0);
    j=0;
    
    for(i=1;i<NUMCONF; i++)
    {
        e = SequenceEnergy(Seq, i);
        if (e<emin) {emin=e; j = i; }
    }
    
    *atwhatconform = j;
    return emin;
}


/*
 double GetTwoStructureQScore(int StrID1, int StrID2)
 calculates the normalized structure similarity (Q-score)
 between structures StrID1, StrID2,
 0<= StrID1, StrID2 < NUMCOUNF,
 0<= output <= 1
 */
double GetTwoStructureQScore(int StrID1, int StrID2)
{
    double qscore=0;
    int i,j,k, t1,t2;
    for(i=0;i<HALFCONTMATRLEN;i++)
    {
        for(j=0;j<HALFCONTMATRLEN;j++)
        {
            t1=0; t2=0;
            for(k=0;k<HALFCONTMATRLEN;k++)
            {
                if ( (ContactMatrixA[StrID1][k]==i)&&
                    (ContactMatrixB[StrID1][k]==j)   ) t1=1;
                if ( (ContactMatrixA[StrID2][k]==i)&&
                    (ContactMatrixB[StrID2][k]==j)   ) t2=1;
            }
            qscore = qscore+t1*t2;
        }
    }
    qscore = qscore / (double) HALFCONTMATRLEN;
    return qscore;
}

/*
 Non-normalized integer structure similarity, 0..28
 */
int GetTwoStructureIntQScore(int StrID1, int StrID2)
{
    int qscore=0;
    int i,j,k, t1,t2;
    for(i=0;i<HALFCONTMATRLEN;i++)
    {
        for(j=0;j<HALFCONTMATRLEN;j++)
        {
            t1=0; t2=0;
            for(k=0;k<HALFCONTMATRLEN;k++)
            {
                if ( (ContactMatrixA[StrID1][k]==i)&&
                    (ContactMatrixB[StrID1][k]==j)   ) t1=1;
                if ( (ContactMatrixA[StrID2][k]==i)&&
                    (ContactMatrixB[StrID2][k]==j)   ) t2=1;
            }
            qscore = qscore+t1*t2;
        }
    }
    return qscore;
}

/*
 void CreateRandomSequence(int *Seq, int Len)
 Creates a random (amino acid) sequence Seq of length Len,
 each element of Seq is between 0 and ALPHABET-1 (inclusive)
 */
void CreateRandomSequence(int *Seq, int Len)
{
    int i,j;
    for(i=0;i<Len;i++)
    {
        do{
            j = (int)( ((double)JKISS()/THE_MAX) * ALPHABET );
        } while(j==ALPHABET);
        Seq[i] = j;
    }
    return;
}


void CopySeq(int *dest, int *src, int Len)
{
    int i;
    for(i=0; i<Len; i++) dest[i]=src[i];
    return;
}

/*
 void LetterToMJCodeSeq(char *letterseq, int *Seq, int Len)
 Converts an amino letter sequence letterseq of length Len
 to the integer representation according to MJ96 energy matrix ordering
 Used to read the sequences.
 CMFIL -> 0 1 2 3 4
 */
void LetterToMJCodeSeq(char *letterseq, int *Seq, int Len)
{
    int i,k;
    
    for(i=0;i<Len;i++)
    {
        k=-1;
        switch (letterseq[i])
        {
            case 'C': { k=0;  break; }
            case 'M': { k=1;  break; }
            case 'F': { k=2;  break; }
            case 'I': { k=3;  break; }
            case 'L': { k=4;  break; }
            case 'V': { k=5;  break; }
            case 'W': { k=6;  break; }
            case 'Y': { k=7;  break; }
            case 'A': { k=8;  break; }
            case 'G': { k=9;  break; }
            case 'T': { k=10; break; }
            case 'S': { k=11; break; }
            case 'N': { k=12; break; }
            case 'Q': { k=13; break; }
            case 'D': { k=14; break; }
            case 'E': { k=15; break; }
            case 'H': { k=16; break; }
            case 'R': { k=17; break; }
            case 'K': { k=18; break; }
            case 'P': { k=19; break; }
        }
        Seq[i]=k;
    }
    
    return;
}

/*
 void PrintAACodeSequence(char *buf, int *Seq, int Len)
 buf receives a letter representation of sequence Seq of length Len
 Normally used for output
 */
void PrintAACodeSequence(char *buf, int *Seq, int Len)
{
    int i;
    char c;
    for(i=0;i<Len;i++)
    {
        c = 'Z';
        switch (Seq[i])
        {
            case  0: { c = 'C'; break; }
            case  1: { c = 'M'; break; }
            case  2: { c = 'F'; break; }
            case  3: { c = 'I'; break; }
            case  4: { c = 'L'; break; }
            case  5: { c = 'V'; break; }
            case  6: { c = 'W'; break; }
            case  7: { c = 'Y'; break; }
            case  8: { c = 'A'; break; }
            case  9: { c = 'G'; break; }
            case 10: { c = 'T'; break; }
            case 11: { c = 'S'; break; }
            case 12: { c = 'N'; break; }
            case 13: { c = 'Q'; break; }
            case 14: { c = 'D'; break; }
            case 15: { c = 'E'; break; }
            case 16: { c = 'H'; break; }
            case 17: { c = 'R'; break; }
            case 18: { c = 'K'; break; }
            case 19: { c = 'P'; break; }
        }
        buf[i]=c;
    }
    buf[i]=0;
    return;
}

double GetStructurePnat(int *Seq, double T, int structid)
{
    int i;
    double z=0.0e0, emin;
    emin=SequenceEnergy(Seq,structid);
    for(i=0;i<NUMCONF;i++) {
        z+=exp(-(SequenceEnergy(Seq,i)-emin)/T);
    }
    return 1.0e0/z;
}

/*
 double GetSequencePnat(int *Seq, double T, int *conform)
 
 Returns the probability of being in the ground state
 for the sequence Seq of length AASEQLEN, at temperature T.
 The index of conformation of the ground state is returned
 via *conform
 
 Alters the espectrum[] global array.
 */
double GetSequencePnat(int *Seq, double T, int *conform)
{
    int i,j;
    double emin,z,p;
    
    emin = 1e20; j=-1;
    
    for(i=0; i<NUMCONF; i++)
    {
        espectrum[i] =  SequenceEnergy(Seq,i);
        if (espectrum[i]<emin)  { emin=espectrum[i]; j=i; }
    }
    *conform = j;
    
    // subtract emin to avoid exp(-E) overflows in partition function
    for(i=0; i<NUMCONF; i++) espectrum[i] = espectrum[i]-emin;
    emin = 0;
    
    z = 0.; for(i=0; i<NUMCONF; i++) z = z + exp(-espectrum[i]/T);
    
//    if (state_3rd>0.0){
  //      z=z+state_3rd; //2.6e27
//    }
    
    p= exp(-emin/T)/z;
    //printf("P=%f\n", p);
    return p;
}

/*
 void PointMutateAASequence(int *Seq, int Len)
 introduces a point mutation in the amino acid sequence Seq of length Len
 */
void PointMutateAASequence(int *Seq, int Len)
{
    int i,j;
    //int k;
    
    do{
        j = (int)( ((double)JKISS()/THE_MAX) * Len );
    } while(j==Len);
    
    do
    {
        do{
            i = (int)( ((double)JKISS()/THE_MAX) * ALPHABET );
        } while(i==ALPHABET);
    }
    while(Seq[j]==i);
    Seq[j]=i;
    
    return;
}


/*
 void SwapTwoResidues(int *Seq, int Len)
 swaps two residues in the amino acid sequence Seq of length Len
 */
void SwapTwoResidues(int *Seq, int Len)
{
    int i,j,k;
    
    do{
        i = (int)( ((double)JKISS()/THE_MAX) * Len );
    } while(i==Len);
    
    do{
        j = (int)( ((double)JKISS()/THE_MAX) * Len );
    } while((j==Len)||(j==i));
    
    k = Seq[i]; Seq[i]=Seq[j]; Seq[j]=k;
    
    return;
}


/*
 double MutateAASeqForMaxPnat(int *aaseq, double Tenv, double Tsel)
 Monte Carlo optimization of sequence aaseq to maximize Pnat
 at temperature Tenv. Tsel is the Metropolis MC selection temperature.
 output : aaseq = optimized sequence
 returns: Pnat of the new sequence
 */
double MutateAASeqForMaxPnat(int *aaseq, double Tenv, double Tsel)
{
    int targetseq[100];
    int aaseq2[100], aaseq3[100];
    int attempt, conform;
    //int i,j,k,ii;
    //double tt
    double  e0, e1, emax, r ;
    //FILE *outfile, *outfile2;
    //double tmean=0;
//    int struct1, struct2;
//    int changetime;
    
    emax = GetSequencePnat(aaseq, Tenv, &conform);
    e0 = emax;
//    struct1 = conform;
//    changetime = 0;
    
    CopySeq(targetseq, aaseq, AASEQLEN); // targetseq is the eventual output
    CopySeq(aaseq2, aaseq, AASEQLEN);
    
    for(attempt=0; attempt<MCATTEMPT; attempt++)
    {
        CopySeq(aaseq3, aaseq2, AASEQLEN);
        PointMutateAASequence(aaseq3, AASEQLEN);
        
        e1 = GetSequencePnat(aaseq3, Tenv, &conform );
//        struct2 = conform;
        
        //if (e1>targetpnat+0.05) break;
        
        if (e1>emax) { emax=e1;  CopySeq(targetseq, aaseq3, AASEQLEN); }
        
        if (e1>e0) { // Pnat increased after mutation, always accept
            e0 = e1; CopySeq(aaseq2, aaseq3, AASEQLEN);
            //        if (struct1!=struct2) { struct1=struct2; fprintf(globaloutfile,"%d %d\n",attempt, attempt-changetime);  changetime = attempt; }
            //        if (fabs(e1-targetpnat)<0.01) { attempt=MCATTEMPT+1;}
        }
        else { // Pnat decreased - run MC game
            r = (double)JKISS()/THE_MAX;
            if ( r < exp((e1-e0)/Tsel) ) //accept
            { e0 = e1;  CopySeq(aaseq2, aaseq3, AASEQLEN);
                //        if (struct1!=struct2) { struct1=struct2; fprintf(globaloutfile,"%d %d\n",attempt, attempt-changetime);  changetime = attempt; }
            }
        }
        
    } // end for attempt
    
    CopySeq(aaseq, targetseq, AASEQLEN);
    e1 = GetSequencePnat( aaseq, Tenv, &conform );
    
    return e1;
}


/*
 double SwapResiduesForMinE(int *aaseq, double Tsel)
 Monte Carlo optimization of sequence aaseq to miminize native state energy
 by swapping amino acids (=at fixed concentration)
 Tsel is the Metropolis MC selection temperature.
 output : aaseq = optimized sequence
 returns: E0 of the new sequence
 */
double SwapResiduesForMinE(int *aaseq,  double Tsel)
{
    int targetseq[100];
    int aaseq2[100], aaseq3[100];
    int attempt, conform;
    //int i,j,k,ii;
    //double tt
    double  e0, e1, emax, r ;
    //FILE *outfile, *outfile2;
    //double tmean=0;
    
    emax = GetSequenceMinEnergy(aaseq, &conform);
    e0 = emax;
    
    CopySeq(targetseq, aaseq, AASEQLEN); // targetseq is the eventual output
    CopySeq(aaseq2, aaseq, AASEQLEN);
    
    for(attempt=0; attempt<MCATTEMPT; attempt++)
    {
        CopySeq(aaseq3, aaseq2, AASEQLEN);
        SwapTwoResidues(aaseq3, AASEQLEN);
        
        e1 = GetSequenceMinEnergy(aaseq3,  &conform );
        
        if (e1<emax) { emax=e1;  CopySeq(targetseq, aaseq3, AASEQLEN); }
        
        if (e1<e0) { // Energy decreased after mutation, always accept
            e0 = e1; CopySeq(aaseq2, aaseq3, AASEQLEN);
        }
        else { // energy inreased - run MC game
            r = (double)JKISS()/THE_MAX;
            if ( r < exp(-(e1-e0)/Tsel) ) //accept
            { e0 = e1;  CopySeq(aaseq2, aaseq3, AASEQLEN);  }
        }
        
    } // end for attempt
    
    CopySeq(aaseq, targetseq, AASEQLEN);
    e1 = GetSequenceMinEnergy( aaseq,  &conform );
    
    return e1;
}

double GetHydrophobicity(int *aaseq, int currLen){
    int i;
    double h=0.0e0;
    for(i=0;i<currLen;i++) {
        h+=Hydrophobicity[aaseq[i]];
        //fprintf(stdout,"%6d %6d %lf\n", i, aaseq[i], Hydrophobicity[aaseq[i]]);
    }
    //fprintf(stdout, "Hydrophobicity : %lf\n", h);
    //getchar();
    return h;
}

double GetNetCharge(int *aaseq, int currLen){
    int i;
    double h=0.0e0;
    for(i=0;i<currLen;i++) {
        h+=Charge[aaseq[i]];
        //fprintf(stdout,"%6d %6d %lf\n", i, aaseq[i], Hydrophobicity[aaseq[i]]);
    }
    //fprintf(stdout, "Hydrophobicity : %lf\n", h);
    //getchar();
    return h;
}

double GetPosCharge(int *aaseq, int currLen){
    int i;
    double h=0.0e0;
    for(i=0;i<currLen;i++) {
        if (Charge[aaseq[i]] > 0.0)
            h+=Charge[aaseq[i]];
        //fprintf(stdout,"%6d %6d %lf\n", i, aaseq[i], Hydrophobicity[aaseq[i]]);
    }
    //fprintf(stdout, "Hydrophobicity : %lf\n", h);
    //getchar();
    return h;
}

double GetNegCharge(int *aaseq, int currLen){
    int i;
    double h=0.0e0;
    for(i=0;i<currLen;i++) {
        if (Charge[aaseq[i]] < 0.0)
            h+=Charge[aaseq[i]];
        //fprintf(stdout,"%6d %6d %lf\n", i, aaseq[i], Hydrophobicity[aaseq[i]]);
    }
    //fprintf(stdout, "Hydrophobicity : %lf\n", h);
    //getchar();
    return h;
}

double GetFracHydrophobicity(int *aaseq, int currLen){
    int i;
    double h=0.0e0;
    for(i=0;i<currLen;i++) {
        if (HydrophobicityYesNo[aaseq[i]] == 1.0)
            h+=1.0;
        //fprintf(stdout,"%6d %6d %lf\n", i, aaseq[i], Hydrophobicity[aaseq[i]]);
    }
    //fprintf(stdout, "Hydrophobicity : %lf\n", h);
    //getchar();
    return h/((double)currLen);
}

double GetFracHydrophobicity_avil(int *aaseq, int currLen){
    int i;
    double h=0.0e0;
    for(i=0;i<currLen;i++) {
        if (HydrophobicityYesNo_avil[aaseq[i]] == 1.0)
            h+=1.0;
        //fprintf(stdout,"%6d %6d %lf\n", i, aaseq[i], Hydrophobicity[aaseq[i]]);
    }
    //fprintf(stdout, "Hydrophobicity : %lf\n", h);
    //getchar();
    return h/((double)currLen);
}


double GetFracCharge(int *aaseq, int currLen){
    int i;
    double h=0.0e0;
    for(i=0;i<currLen;i++) {
        if ((Charge[aaseq[i]] < 0.0) || (Charge[aaseq[i]] > 0.0))
            h+=1.0;
        //fprintf(stdout,"%6d %6d %lf\n", i, aaseq[i], Hydrophobicity[aaseq[i]]);
    }
    //fprintf(stdout, "Hydrophobicity : %lf\n", h);
    //getchar();
    return h/((double)currLen);
}

void ConvertAAtoCharge(int *aaseq, int *new_aaseq, int currLen){
    int i;
    //double h=0.0e0;
    for(i=0;i<currLen;i++) {
        new_aaseq[i]=Charge[aaseq[i]];
    }
}

void ConvertAAtoAVIL(int *aaseq, int *new_aaseq, int currLen){
    int i;
    //double h=0.0e0;
    for(i=0;i<currLen;i++) {
        new_aaseq[i]=HydrophobicityYesNo_avil[aaseq[i]];
    }
}

void ConvertAAtoHydro(int *aaseq, int *new_aaseq, int currLen){
    int i;
    //double h=0.0e0;
    for(i=0;i<currLen;i++) {
        new_aaseq[i]=HydrophobicityYesNo[aaseq[i]];
    }
}

















