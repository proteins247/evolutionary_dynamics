/*************************
 Modified by Orit Peleg
 opeleg@fas.harvard.edu
 on 2/26/13
 ************************/

#ifndef __CELL_H__
#define __CELL_H__

#include "define.h"

#define S_DEAD -1
#define S_ALIVE 1
#define S_NEWBORN 2
#define S_TODIE 3
#define S_TODIE2 4

typedef struct ORGANISM{
    char genome[MAXGENES*NUCSEQLEN];
    int structid[MAXGENES+1];
    int genecount;
    //int action_count[4];
    int SeqID[MAXGENES*(MAXGENES-1)/2];
    
    double C[MAXGENES+1], K[MAXSTATES+1][MAXSTATES+1], F[MAXSTATES+1];
	double nF[MAXSTATES+1][MAXSTATES+1];
    double Kc[MAXGENES+1], Nch[MAXGENES+1];
    float pnat[MAXGENES], hydro[MAXGENES];
    float netcharge[MAXGENES], poscharge[MAXGENES],
        negcharge[MAXGENES], fraccharge[MAXGENES], frachydro[MAXGENES], frachydro_avil[MAXGENES];
    
    float hydro_s_partner[MAXPPIS], netcharge_s_partner[MAXPPIS], poscharge_s_partner[MAXPPIS],
    negcharge_s_partner[MAXPPIS], fraccharge_s_partner[MAXPPIS], frachydro_s_partner[MAXPPIS], frachydro_avil_s_partner[MAXPPIS];
    
    float hydro_s_hub[MAXPPIS], netcharge_s_hub[MAXPPIS], poscharge_s_hub[MAXPPIS],
    negcharge_s_hub[MAXPPIS], fraccharge_s_hub[MAXPPIS], frachydro_s_hub[MAXPPIS], frachydro_avil_s_hub[MAXPPIS];
    
    float nsi[MAXGENES+1], si[MAXGENES+1];
    float nsi2[MAXGENES+1], si2[MAXGENES+1];
    
    double minpnat, meanpnat, G0, birthrate;
    
    int ppicount, ppi_pair[MAXPPIS][2], bmode[MAXPPIS];
    float pint[MAXPPIS+1], Gij[MAXPPIS+1];
    
    int dob; // date of birth
    int numkids;
    int generation;
    float mutcount;
    float mutrate;
    int mutlevel, mutorigin;
    float synmut[MAXGENES], nonsynmut[MAXGENES], tot_mut[MAXGENES];
}
organism;

typedef struct PARAMETER{
    int seed, startcode, orgcount;
    int maxdivcycle;
    int decimthresh, decimto, initpop;
    int dumpcycle, printoutcycle, plotoutcycle, screenoutcycle, seqlogcycle;
    double Tenv, tol;
    double TLow, THigh;
    int timeLow, timeHigh;
    double birthrate, deathrate, expressrate, alpha, pnatcutoff;
    double speciessizefactor;
    double mutrate[2], mutrate0, x0, b0, fixed_mutrate;
    double mutthreshfac;
    double mutthresh;
    double bindindex;
    char input[300], targetname[100];
}
parameter;

//typedef struct GENEDB
//{
//    char genome[NUCSEQLEN];
//    int structid, count, reporg;
//    float pnat, hydro;
//}
//genedb;

typedef struct ORGDB
{
    //int id;
    char genome[MAXGENES*NUCSEQLEN];
    int structid[MAXGENES], binding[4][9], bmode[MAXPPIS];
	float pint[MAXPPIS+1]; 
	float birthrate;
	int count, reporg, genecount;
    float pnat[MAXGENES], hydro[MAXGENES];
	float C[MAXGENES];
	float F[MAXGENES];
	double nF[MAXSTATES+1][MAXSTATES+1];
	float Gij[MAXPPIS+1];
    //float meanmutrate;
}
orgdb;


void GetSequenceID(int who);
int IterativeSolver(int who);
int IterativeSolver_cha(int who);
int IterativeSolver_original(int who);
float GaussianNum();
int UpdateBirthrateStoch(parameter *myParam, int who, int func);
int UpdateMutLevel(parameter *myParam, int who, int func);
int OrgDeath(parameter *myParam, int who, int compensate);
int GeneExpress(parameter *myParam, int who, int compensate, double *dC);
int OrgGeneMutate(parameter *myParam, int who, int compensate);
int OrgChildBirth(parameter *myParam, int who, int divisioncycle, int compensate);
void KillOrganism(int who);
void PrintInitialCondition(FILE *out, parameter *myParam);
void SetupParameter(int argc, char *argv[], parameter *myParam, int *orgcount);
int UpdateMonomerConcentration(parameter *myParam, int who);
void UpdateNsi(parameter *myParam, int who);
int UpdateEquilibriumConstant(parameter *myParam, int who, int func);
int UpdateEquilibriumConstant_op(parameter *myParam, int who, int func);


organism myOrg[MAXORGANISMS];
//genedb myDB[MAXORGANISMS][MAXGENES];
orgdb myOrgDB[MAXORGANISMS];

int myOrgstatus[MAXORGANISMS];
//int nDB[MAXGENES]
int nOrgDB;

void UpdatePPISeqProps(int who);

#endif /* __CELL_H_ */















