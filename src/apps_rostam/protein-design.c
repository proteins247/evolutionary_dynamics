/*
 Protein design on 3x3x3 lattice proteins
 monomers and dimers
 monoclonal evolution 
 employs Kimura theory as the Metropolis criterion
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>


#include "zlib.h"
#include "../src/gencode.h"
#include "../src/latticelib.h"
#include "../src/bindinglib.h"
#include "../RNG/generator.h"

#define GENES 		2


FILE *out1;

int structid[GENES];
int nucseq[GENES][NUCSEQLEN]; 
int oldnucseq[GENES][NUCSEQLEN];
int aaseq[GENES][AASEQLEN];

int ii, jj, kk;
int pi, pj;

double pnat[GENES];
double pint;
double oldfitness;
double fixation;
double selection;

double Tenv;
int pop_size, sim_time, mon_or_dim, protein; 

char tempaaseq[GENES][AASEQLEN];

int bmode;
int best_bmode;

int label;
char fopbuf[100];
char rootdir[100] = "./";

char *dir_md_list[2]={"mon", "dim"};

void Printinfo(){
	FILE *out0;
	int surface[GENES][AASURFACELEN];
	int s_i;

	sprintf(fopbuf, "%s/%s%d/info.txt", rootdir, dir_md_list[mon_or_dim], protein);
	out0 = fopen(fopbuf, "w");
	fprintf(out0, "#General information on simulation. Only printed out when label=0\n");
	fprintf(out0, "Temp: %.3f (kT)\n", Tenv);
	fprintf(out0, "N: %d \n", pop_size);
	
	if (protein==2){
		for (jj=0; jj<GENES; jj++){
			PrintAASequence(tempaaseq[jj], aaseq[jj], AASEQLEN);
			fprintf(out0, "%s\n", tempaaseq[jj]);
			fprintf(out0, "Conform%d: %d, %.3E\n", jj, structid[jj], pnat[jj]);
		}
	}
	else{
		PrintAASequence(tempaaseq[protein], aaseq[protein], AASEQLEN);
		fprintf(out0, "%s\n", tempaaseq[protein]);
		fprintf(out0, "Conform%d: %d, %.3E\n", protein, structid[protein], pnat[protein]);
	}
	if (mon_or_dim==1){	
		fprintf(out0, "Bmode: %d, %.3E\n", bmode, pint);
		GetSurfaceAAPositions(structid[pi], structid[pj], bmode, surface[0], surface[1]);
		for (s_i=0;s_i<AASURFACELEN; s_i++) fprintf(out0, " %02d", surface[0][s_i]);
		fprintf(out0, "\n");
		for (s_i=0;s_i<AASURFACELEN; s_i++) fprintf(out0, " %02d", surface[1][s_i]);
	}
	fflush(out0);	
	fclose(out0);
}

void Openfiles(){
	struct stat st = {0};

	sprintf(fopbuf, "%s/%s%d", rootdir, dir_md_list[mon_or_dim], protein);
	if (stat(fopbuf, &st) == -1) {		// check if directory exists, otherwise make it
		mkdir(fopbuf, 0700);
	}
	sprintf(fopbuf, "%s/%s%d/v%d.dat", rootdir, dir_md_list[mon_or_dim], protein, label);
    out1 = fopen(fopbuf, "w");

	fprintf(out1, "#Output only contains accepted mutations. Thus, no output between attempted mutation 0 and 11 means that no mutation has been accepted during that time and populations from attempted mutation 0-10 have sequences identical to that outputted in attempted mutation 0\n");
	if (protein==2){
		fprintf(out1, "#Attempted mutation | Energy of protein 0 functional fold  | Energy of protein 1 functional fold | Energy of functional binding mode | Probability to be in functional fold of protein 0 | Probability to be in functional fold of protein 1 | Protability to be in functional binding mode | Fitness | Probability mutation is accepted | Amino acid sequence of protein 0 | Amino acid sequence of protein 1\n");
	}
	else{
		fprintf(out1, "#Attempted mutation | Probability to be in functional fold ");
		if (mon_or_dim==1){
			fprintf(out1, "| Probability to be in functional binding mode ");
		}
		fprintf(out1, "| Fitness | Probability mutation is accepted | DNA sequence\n");	
	}
}


void Printout(){
	fprintf(out1, "%d", ii);
	if (protein==2){
		for (jj=0; jj<GENES; jj++){
			fprintf(out1, " %.3E", SequenceEnergy(aaseq[jj], structid[jj]));
		}	
		fprintf(out1, " %.3E", GetBindingEnergy2(aaseq[0], structid[0], aaseq[1], structid[1], bmode));	
		for (jj=0; jj<GENES; jj++){
			fprintf(out1, " %.3E", pnat[jj]);
		}	
		fprintf(out1, " %.3E", pint);	
	}
	else{
		fprintf(out1, " %.3E", pnat[protein]);
		if (mon_or_dim==1){		
			fprintf(out1, " %.3E", pint);	
		}
	}
	fprintf(out1, " %.3E", oldfitness);	
	fprintf(out1, " %.3E", fixation);	
	if (protein==2){	
		for (jj=0; jj<GENES; jj++){
			fprintf(out1, " ");
			PrintAASequence(tempaaseq[jj], aaseq[jj], AASEQLEN);
			fprintf(out1, "%s ", tempaaseq[jj]);
//			for (kk=0; kk<NUCSEQLEN; kk++){
//				fprintf(out1, "%d", nucseq[jj][kk]);
//			}
		}
	}	
	else{	
		fprintf(out1, " ");
		for (kk=0; kk<NUCSEQLEN; kk++){
			fprintf(out1, "%d", nucseq[protein][kk]);
		}
	}
	fprintf(out1, "\n");
	fflush(out1);
}

void ReadCommondata(){
	char rootdir[100], file[100];
	sprintf(rootdir,"%s","../commondata");
	sprintf(file, "%s/LPforms/contact10000.dat", rootdir);
	ReadContactMatrix(file);
	sprintf(file, "%s/MJ96/energy.dat", rootdir);
	ReadEnergyMatrix(file);

	if (mon_or_dim==1){
		sprintf(file, "%s/LPforms/allfaces10000.dat", rootdir);
		ReadAllSurfaces(file);	
	}
}

double calc_fitness(){
	double fitness;
	double alpha=10.0;

	if (mon_or_dim==0){
		fitness = pnat[protein];
	}
	else if (mon_or_dim==1){
		if (protein==2){
//			fitness = pnat[0]*pnat[1]*pint;
			fitness = pnat[0]*pnat[1]*pow(pint, alpha);

		}
		else{
			fitness = pow(pnat[protein], 2)*pint;
		}
	}
	else{
		printf("Fitness is undefined\n");
		exit(0);
	}
	return fitness;
}


int main(int argc, char *argv[]){
	FILE *in0; char buf[100];
	double fitness; 
	int status;

	if (argc != 7){
		printf("\nprotein-design takes 6 arguments\n");
		printf("Temperature | Population size | Attempted mutations | Monomer/dimer | Protein | Label\n\n");
		if (strcmp(argv[1], "--help")==0 || strcmp(argv[1], "-h")==0){
			printf("Description || Possible values || Example\n");
			printf("1) Temperature in units of room temperature || any positive float || 1.0, meaning the the protein(s) fold (and bind) at root temperature\n");
			printf("2) Population size || any positive integer || 10, meaning the sequence undergoes evolution with a population size of 10 organisms\n");
			printf("3) Number of mutations attempted in orders of magnitude || any positive integer || 2, meaning 100 (10^2) mutations will be attempted\n");
			printf("4) Monomer or dimer || 0 or 1 || 0, meaning protein-design will be performed for the monomer\n");
			printf("5) Protein which undergoes monomeric/dimeric evolution || 0, 1, or 2 || 0 or 1, meaning the sequence on the first or second line undergo monomeric/homodimeric evolution; 2, meaning that heterodimeric evolution is performed\n"); 
			printf("6) Label of output file, useful for running multiple jobs at once || any positive integer || 0, the output file can be found in the appropriate directory as v0.dat\n");
		}
		exit(-1);

	}

	printf("Reading initial arguments...\n");

	Tenv = atof(argv[1]);
	pop_size = atoi(argv[2]);	
	sim_time = atoi(argv[3]);
	mon_or_dim = atoi(argv[4]);
	protein = atoi(argv[5]);	
	label = atoi(argv[6]);

	if (protein==2 && mon_or_dim==0){
		printf("Protein = 2 corresponds to heterodimer evolution and must have Monomer/dimer = 1\n");
		exit(-1);
	}

	sprintf(fopbuf, "%s/seq.txt", rootdir);
	if ((in0=fopen(fopbuf,"r")) == NULL) { printf("Can't open %s file\n", fopbuf); exit(-1);}
	printf("Successfully read sequence file: %s file\n", fopbuf);

	jj = 0;
	while ((status = fscanf(in0,"%s\n", buf)) != EOF) {
		LetterToNucCodeSeq(buf, nucseq[jj], NUCSEQLEN);
		jj++;
	}
	fclose(in0);

	ReadCommondata();
	printf("Done reading initial arguments.\n\n");

	printf("Calculating initial values...\n");
	if (protein==2){
		for (jj=0; jj<GENES; jj++){
			NucSeqToAASeq(nucseq[jj], NUCSEQLEN, aaseq[jj]);
			pnat[jj] = GetSequencePnat(aaseq[jj], Tenv, structid+jj);	
			CopySeq(oldnucseq[jj], nucseq[jj], NUCSEQLEN);
		}
		pi=0;
		pj=1;
	}
	else{
		pi=pj=protein;
		NucSeqToAASeq(nucseq[pi], NUCSEQLEN, aaseq[pi]);
		pnat[pi] = GetSequencePnat(aaseq[pi], Tenv, structid+pi);	
		//pnat[pi] = GetStructurePnat(aaseq[pi], Tenv, structid[pi]);	// set structid
		CopySeq(oldnucseq[pi], nucseq[pi], NUCSEQLEN);
	}

	if (mon_or_dim==1){
		pint = GetBindingP(aaseq[pi], structid[pi], aaseq[pj], structid[pj], &(best_bmode), Tenv);
		printf("best initial bmode: %d, %.3E\n", best_bmode, pint);
//		bmode=best_bmode;
		bmode=0;
		pint = GetBindingP2(aaseq[pi], structid[pi], aaseq[pj], structid[pj], bmode, Tenv);		//necessary if set bmode!=best_bmode
	}

	oldfitness = calc_fitness();

	Openfiles();
	Printout();
	if (label==0) Printinfo(); //generic info file
	printf("Done calculating initial values.\n\n");

	printf("Finding stabilizing mutations...\n");
	init_KISS();	//initialize random number generator

	for(ii=1; ii<pow(10, sim_time); ii++){
		if (protein==2){
			jj = JKISS()%2; //randomly choose which protein to mutate for heterodimer
		}
		else{
			jj = protein;
		}

		status = PointMutateNucSequence(nucseq[jj], NUCSEQLEN);
//		status = NucSeqToAASeq(nucseq[jj], NUCSEQLEN, aaseq[jj]);

		if (status==-1) {CopySeq(nucseq[jj], oldnucseq[jj], NUCSEQLEN);}		// Reject mutations that introduce a stop codon
		else if (status==1){
			NucSeqToAASeq(nucseq[jj], NUCSEQLEN, aaseq[jj]);
			pnat[jj] = GetStructurePnat(aaseq[jj], Tenv, structid[jj]);
			
			if (mon_or_dim==1){
			pint = GetBindingP2(aaseq[pi], structid[pi], aaseq[pj], structid[pj], bmode, Tenv);
			}			

			fitness = calc_fitness();
			selection = (fitness - oldfitness)/oldfitness;
			fixation = (1-exp(-2*selection))/(1-exp(-2*pop_size*selection));
			if (fixation > (double)JKISS()/THE_MAX){
				CopySeq(oldnucseq[jj], nucseq[jj], NUCSEQLEN); oldfitness=fitness; 
				NucSeqToAASeq(nucseq[1-jj], NUCSEQLEN, aaseq[1-jj]);		//necessary for printing out aaseq 
				Printout();
			}
			else{CopySeq(nucseq[jj], oldnucseq[jj], NUCSEQLEN);}
		}
		else{} //status=0, synonymous mutation
	}
	fclose(out1);
	printf("%.0f mutations attempted.\n", pow(10,sim_time));
	return 0;
}
