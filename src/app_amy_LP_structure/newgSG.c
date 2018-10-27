#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "zlib.h"
#include "../LP/gencode.h"
#include "../LP/latticelib.h"
#include "general.h"
#include "structurelib.h"
#include "../PPI/bindinglib.h"
#include "../RNG/generator.h"

void Printinfo();
void Printout();

FILE *out1;

int structid[GENES];// = {6682,6682};
int nucseq[GENES][NUCSEQLEN]; 
int aaseq[GENES][AASEQLEN];


int ii, jj, kk;

double pnat[GENES];
double pint;
double oldfitness;
double fixation;
double selection;

int pop_size, mono;
int homo; 
int bmode_label;

double Tenv = 1.0;
int sim_time;

char tempaaseq[GENES][AASEQLEN];

int bmode;
int best_bmode;

int label;
char fopbuf[100];
char rootdir[100] = "/n/regal/shakhnovich_lab/rrazban/";

char homo_list[3][6]={"homo1", "homo2", "hetero"};

void Printinfo(){
	FILE *out0;
	int surface[GENES][AASURFACELEN];
	int s_i;
	char fopbuf[100];
	char rootdir[100];

	sprintf(rootdir,"/n/regal/shakhnovich_lab/rrazban/");
 
	sprintf(fopbuf, "%s/bmode%d/%s/info.txt", rootdir, bmode_label, homo_list[homo]);
	out0 = fopen(fopbuf, "w");


//print out aa
	for (jj=0; jj<GENES; jj++){
		PrintAASequence(tempaaseq[jj], aaseq[jj], AASEQLEN);
		fprintf(out0, "%s\n", tempaaseq[jj]);
	}

	fprintf(out0, "Temp: %.3f (kT)\n", Tenv);
	fprintf(out0, "N: %d \n", pop_size);
	fprintf(out0, "mono: %d \n", mono);
	fprintf(out0, "homo: %d \n", homo);
	for (jj=0; jj<GENES; jj++) fprintf(out0, "Conform%d: %d, %.3E\n", jj, structid[jj], pnat[jj]);
	fprintf(out0, "Bmode: %d, %.3E\n", bmode, pint);
	GetSurfaceAAPositions(structid[0], structid[1], bmode, surface[0], surface[1]);

	for (s_i=0;s_i<AASURFACELEN; s_i++) fprintf(out0, " %02d", surface[0][s_i]);
	fprintf(out0, "\n");
	for (s_i=0;s_i<AASURFACELEN; s_i++) fprintf(out0, " %02d", surface[1][s_i]);

	fflush(out0);	
	fclose(out0);
}

void Openfiles(){
	sprintf(fopbuf, "%s/bmode%d/%s/dats/v%d.dat", rootdir, bmode_label, homo_list[homo], label);
    out1 = fopen(fopbuf, "w");
}

void Printout(){
	fprintf(out1, "%d", ii);
	for (jj=0; jj<GENES; jj++){
		fprintf(out1, " %.3E", pnat[jj]);
	}	
	fprintf(out1, " %.3E", pint);	
	fprintf(out1, " %.3E", oldfitness);	
	fprintf(out1, " %.3E", fixation);	
	for (jj=0; jj<GENES; jj++){
		fprintf(out1, " ");
		for (kk=0; kk<NUCSEQLEN; kk++){
		fprintf(out1, "%d", nucseq[jj][kk]);
		}
	}	
	fprintf(out1, "\n");
	fflush(out1);
}

//Generates stable AA sequence of conf specified in argv[1]
int main(int argc, char *argv[]){
	FILE *in0; char buf[100];
	int oldnucseq[GENES][NUCSEQLEN];
	double fitness; 
	int status;

	if (argc != 7){
		printf("time | pop_size | mono/dimer | homo1/homo2/hetero | label | bmode_label\n");
		printf("something is missing, bro\n");
		exit(-1);
	}

	init_KISS();

	ReadCommondata();
	sim_time = atoi(argv[1]);	//order of mag
	pop_size = atof(argv[2]);	
	mono = atoi(argv[3]);	//0 for single monomer (not both)
	homo = atoi(argv[4]);		//0 or 1 for homo (selection for which one is homo)
	label = atoi(argv[5]);
	bmode_label = atoi(argv[6]);

	printf("Reading initial sequence...\n");

	sprintf(fopbuf, "%s/bmode%d/seq.txt", rootdir, bmode_label);

	if ((in0=fopen(fopbuf,"r")) == NULL) { printf("Can't open %s file\n", fopbuf); exit(-1);}
	jj = 0;
    
	while ((status = fscanf(in0,"%s\n", buf)) != EOF) {
		if (jj == homo){
			LetterToNucCodeSeq(buf, nucseq[jj], NUCSEQLEN);
			LetterToNucCodeSeq(buf, nucseq[1-jj], NUCSEQLEN);
			break;
		}
		else{
		LetterToNucCodeSeq(buf, nucseq[jj], NUCSEQLEN);
		}
		jj++;
	}
	fclose(in0);
	
	for (jj=0; jj<GENES; jj++){
//		LetterToNucCodeSeq(tempnucseq[jj], nucseq[jj], NUCSEQLEN);
		NucSeqToAASeq(nucseq[jj], NUCSEQLEN, aaseq[jj]);
		pnat[jj] = GetSequencePnat(aaseq[jj], Tenv, structid+jj);	//predefined nat
//		pnat[jj] = GetStructurePnat(aaseq[jj], Tenv, structid[jj]);
		CopySeq(oldnucseq[jj], nucseq[jj], NUCSEQLEN);
	}

	pint = GetBindingP(aaseq[0], structid[0], aaseq[1], structid[1], &(best_bmode), Tenv);

	printf("best initial bmode: %d, %.3E\n", best_bmode, pint);
	bmode=best_bmode;
	pint = GetBindingP2(aaseq[0], structid[0], aaseq[1], structid[1], bmode, Tenv);
	if (mono==0){oldfitness = pnat[0];}
	else {oldfitness = pnat[0]*pnat[1]*pint;}
	if (label==0) Printinfo();
	Openfiles();
	Printout();

	// Find stabilizing mutations 
	for(ii=1; ii<pow(10, sim_time); ii++){
		jj = RandomBit(); 
		PointMutateNucSequence(nucseq[jj], NUCSEQLEN);
		status = NucSeqToAASeq(nucseq[jj], NUCSEQLEN, aaseq[jj]);
		if (status==1) {CopySeq(nucseq[jj], oldnucseq[jj], NUCSEQLEN);}// Reject mutations that introduce a stop codon
		else{
			pnat[jj] = GetStructurePnat(aaseq[jj], Tenv, structid[jj]);
			if (homo<2){
				CopySeq(nucseq[1-jj], nucseq[jj], NUCSEQLEN);
				NucSeqToAASeq(nucseq[1-jj], NUCSEQLEN, aaseq[1-jj]);
				pnat[1-jj] = GetStructurePnat(aaseq[1-jj], Tenv, structid[1-jj]);
			}

			pint = GetBindingP2(aaseq[0], structid[0], aaseq[1], structid[1], bmode, Tenv);
			
			if (mono==0){fitness = pnat[0];}
			else {fitness = pnat[0]*pnat[1]*pint;}
			selection = (fitness - oldfitness)/oldfitness;
			fixation = (1-exp(-2*selection))/(1-exp(-2*pop_size*selection));
			if (fixation > (double)JKISS()/THE_MAX){
				CopySeq(oldnucseq[jj], nucseq[jj], NUCSEQLEN); oldfitness=fitness; Printout();
			}
			else{CopySeq(nucseq[jj], oldnucseq[jj], NUCSEQLEN);}
		}
	}
	fclose(out1);
	return 0;
}
