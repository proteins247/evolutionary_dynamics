/*
 * gencode.c
 *
 *  Genetic code and nucleotide sequences routines
 *  (c) 2005 K.Zeldovich
 *  Created on: Dec 25, 2009
 *      Author: mheo
 */

#include<stdio.h>
#include<stdlib.h>
#include"gencode.h"
#include"../RNG/generator.h"
// U C A G -> 0 1 2 3


int TripletToAA(unsigned short int triplet)
{
int out;

out = -2;

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

if (out==-2){
fprintf(stderr,"genetic code error\n");
fprintf(stderr,"triplet: %.3x\n",triplet);
}

return out;
}
int NucSeqToAASeq(int *NucSeq, int N, int *AASeq)
{
unsigned short int a;
unsigned char a1, a2, a3;
int k=0,i,j=0;

for(i=0;i<N;i+=3)
{
a1 = NucSeq[i];
a2 = NucSeq[i+1];
a3 = NucSeq[i+2];

a = 0;
a = (a1 << 8 ) | (a2 << 4) | a3;

//printf(" %x %x %x -> 0x%.4x\n",a1,a2,a3, a);

AASeq[j] = TripletToAA(a);
if (AASeq[j]==-1) { //fprintf(stderr,"%.3x to stop\n",a);*/
 k=1; }
if (AASeq[j]==-2) {
  fprintf(stderr,"NucSeqToAASeq() : Genetic Code Error!!!\n");
  exit(0);
}
j++;
}

return k;
}



int CharNucSeqToAASeq(char *NucSeq, int N, int *AASeq)
{
unsigned short int a;
unsigned char a1, a2, a3;
int k=0,i,j=0;

for(i=0;i<N;i+=3)
{
a1 = NucSeq[i];
a2 = NucSeq[i+1];
a3 = NucSeq[i+2];

a = 0;
a = (a1 << 8 ) | (a2 << 4) | a3;

//printf(" %x %x %x -> 0x%.4x\n",a1,a2,a3, a);

AASeq[j] = TripletToAA(a);
if (AASeq[j]==-1) { //fprintf(stderr,"%.3x to stop\n",a);*/
 k=1; }
if (AASeq[j]==-2) {
  fprintf(stderr,"CharNucSeqToAASeq() : Genetic Code Error!!!\n");
   //exit(0);
}
 j++;
}

return k;
}

void CreateRandomNucSequence2(int *Seq, int Len)
{
  unsigned short int a;
  unsigned char a1, a2, a3;
  int i,j, aa;
  for(i=0;i<Len/3;i++) {
    do {
      for(j=0;j<3;j++){
        do{aa = (int)(((double)JKISS()/THE_MAX) * 4 );} while(aa==4);
        Seq[3*i+j]=aa;
      }
      a1 = Seq[3*i];
      a2 = Seq[3*i+1];
      a3 = Seq[3*i+2];
      a = 0;
      a = (a1 << 8 ) | (a2 << 4) | a3;
      aa=TripletToAA(a);
      //printf(" %x %x %x -> 0x%.4x\n",a1,a2,a3, a);
    } while(aa==A_STOP);
  }
  return;
}

void CreateRandomNucSequence(int *Seq, int Len)
{
int i,j;
for(i=0;i<Len;i++) {
do{
j = (int)( ((double)JKISS()/THE_MAX) * 4 );
} while(j==4);
Seq[i]=j;
}
return;
}

int PointMutateNucSequence(int *Seq, int Len)
{
int i,j,k;
int AAS1[10000], AAS2[10000];

NucSeqToAASeq(Seq, Len, AAS1);

do{
j = (int)( ((double)JKISS()/THE_MAX) * Len );
} while(j==Len);

do
  {
  do{
  i = (int)( ((double)JKISS()/THE_MAX) * 4 );//why like this? multiple muts?? oh keep on going until mrna is changed
  } while(i==4);
  }
while(Seq[j]==i);
Seq[j]=i;


k = NucSeqToAASeq(Seq, Len, AAS2);

if (k==1) return -1;

j=0;
for(i=0; i<Len/3; i++)
if (AAS1[i]!=AAS2[i]) {j=1; break;}

return j;
}



int PointMutateCharNucSequence(char *Seq, int Len)
{
int i,j,k;
int AAS1[10000], AAS2[10000];

CharNucSeqToAASeq(Seq, Len, AAS1);

do{
j = (int)( ((double)JKISS()/THE_MAX) * Len );
} while(j==Len);

do
  {
    do{
      i = (int)( ((double)JKISS()/THE_MAX) * 4 );
      } while(i==4);
  }
while(Seq[j]==i);
Seq[j]=i;

k = CharNucSeqToAASeq(Seq, Len, AAS2);

if (k==1) return -1;

j=0;
for(i=0; i<Len/3; i++)
if (AAS1[i]!=AAS2[i]) {j=1; break;}
return j;
}



void LetterToNucCodeSeq(char *buf, int *Seq, int Len)
{
int i;
int c;
for(i=0;i<Len;i++)
{
c = -1;
switch (buf[i])
{
case  'U': { c = 0; break; }
case  'C': { c = 1; break; }
case  'A': { c = 2; break; }
case  'G': { c = 3; break; }
}
Seq[i]=c;

}
return;
}

void PrintAASequence(char *buf, int *Seq, int Len)
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
case  6: { c = 'T'; break; }
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


void PrintNucCodeSequence(char *buf, int *Seq, int Len)
{
int i;
char c;
for(i=0;i<Len;i++)
{
c = 'Z';
switch (Seq[i])
{
case  0: { c = 'U'; break; }
case  1: { c = 'C'; break; }
case  2: { c = 'A'; break; }
case  3: { c = 'G'; break; }
}
buf[i]=c;
}
buf[i]=0;
return;
}

void PrintCharNucCodeSequence(char *buf, char *Seq, int Len)
{
int i;
char c;
for(i=0;i<Len;i++)
{
c = 'Z';
switch (Seq[i])
{
case  0: { c = 'U'; break; }
case  1: { c = 'C'; break; }
case  2: { c = 'A'; break; }
case  3: { c = 'G'; break; }
}
buf[i]=c;
}
buf[i]=0;
return;
}


void CopyIntToCharSeq(char *dest, int *src, int Len)
{
int i;
for(i=0; i<Len; i++) dest[i]=src[i];
return;
}
