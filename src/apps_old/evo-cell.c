/*************************
 Modified by Orit Peleg
 opeleg@fas.harvard.edu
 on 2/26/13
 ************************/


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <sys/times.h>
#include <time.h>
#include <errno.h>
#include <stddef.h>

#include "zlib.h"
#include "define.h"
#include "../LP/gencode.h"
#include "../LP/latticelib.h"
#include "../PPI/bindinglib.h"
#include "cell.h"
#include "compare.h"


char aacode[] = {
    'C', 'M', 'F', 'I', 'L',
    'V', 'W', 'Y', 'A', 'G',
    'T', 'S', 'N', 'Q', 'D',
    'E', 'H', 'R', 'K', 'P'
};

int aamat[AASEQLEN][20];
int aaseq[AASEQLEN], aaseq2[AASEQLEN];
int domi_species;

int curr_MAXGENES;
int curr_MAXPPIS;
int curr_MAXSTATES;
int allow_fold_change;
int allow_gene_exp;
int allow_unfolded_states;
int allow_chaps;
float chaps_x;
int hub_ID;
int singlish;
int POST_Proc;
int sequenceversion;
int selection;
int homo;

int ResetOrgDB(parameter *myParam, int divisioncycle);
void SetMeanOrg(organism *mean);
void AddMeanOrg(organism *mean, int who);
void GetMeanOrg(organism *mean, int orgcount);
int seq_entropy(char **seq, int N, double *entropy);
void SetVarOrg(organism *var);
void AddVarOrg(organism *var, int who);
void GetVarOrg(organism *var, int orgcount);
void GetMeanOrg(organism *var, int orgcount);


clock_t HDP_time_start, HDP_time_end;
double HDP_cpu_time_used;
clock_t HDP_time_start2, HDP_time_end2;
double HDP_cpu_time_used2;

void   WriteConfig();
void   ReadConfig();


//output variables:

void PrintOutput();
void Flushfiles();
void Closefiles();
void Openfiles();
void PrepareOutput();
void PostProcessing();
void PrintHeaders(FILE *out);
void PrintHeaders1(FILE *out);

//FILE *it_solver;
FILE *out1;
FILE *out22, *out23, *out24;
FILE *out9;

parameter myParam;
organism mean;
organism var;

/* parameters for sequence entropy calculation */
double entropy[MAXGENES][AASEQLEN], entropy_sum[MAXGENES];
char *seq[MAXGENES][MAXORGANISMS];
double entropy_sumtot;


int divisioncycle, lastdecimtime;
int mutatororigin[7];
int orgcount;
int mutatorcount, count;


char oldgenome[MAXGENES*NUCSEQLEN];

char seqbuf[800];
int sizeRank[MAXORGANISMS] = {0};
int speciesSizeTotal;

int generation;
double TIME;

int RUN_BASED_ON_CONFIG=0;

int start_divisioncycle;

/***********************************************************************************************************
 **********************************************************************************************************/
int main(int argc, char *argv[]){
   	int j; 
    int kk;
    int ii, jj, kia;
    int who, status;
    int count1, count2, newborncount;
    int fate, overflowflag;
    int deathcount, death2count, maxorgind;
    double dC;
    double dt,cum_br[MAXORGANISMS],max_br,r1,r2;
    
    
    
    
    /*initialize time evaluation:*/
    ReadConfig();
    
    printf("OPEN...\n"); 
    
    SetupParameter(argc, argv, &myParam, &orgcount);
    
    Openfiles();
    if (RUN_BASED_ON_CONFIG==0){start_divisioncycle = 0;} else{start_divisioncycle=divisioncycle;}
    
    
    divisioncycle=start_divisioncycle;
    ResetOrgDB(&myParam, divisioncycle);
    RankSpeciesSizeDB(sizeRank, nOrgDB);
   
	 
    for(j=0;j<myOrgDB[sizeRank[0]].genecount*NUCSEQLEN;j++) {
	oldgenome[j]=myOrgDB[sizeRank[0]].genome[j];
   	} 
    PrepareOutput(); PrintOutput();
	
	fprintf(stderr,"Interaction surface\n");
	fprintf(stderr, "Bmode: %d\n", myOrgDB[sizeRank[0]].bmode[0]);
    int aaseq_surface_hub[AASURFACELEN], aaseq_surface_partner[AASURFACELEN], s_i, hub_i, par_i;
	for (ii=0;ii<curr_MAXPPIS; ii++){
		hub_i=0; par_i=ii+1;
        GetSurfaceAAPositions(myOrgDB[sizeRank[0]].structid[hub_i], myOrgDB[sizeRank[0]].structid[par_i], myOrgDB[sizeRank[0]].bmode[ii],
                                  aaseq_surface_hub, aaseq_surface_partner);
        for (s_i=0;s_i<AASURFACELEN; s_i++) fprintf(stderr, " %02d", aaseq_surface_hub[s_i]);
    	fprintf(stderr, "\n");
	    for (s_i=0;s_i<AASURFACELEN; s_i++) fprintf(stderr, " %02d", aaseq_surface_partner[s_i]);
            
     }
    fprintf(stderr,"\nStart simulations...\n");

    
    // MAIN LOOP
    /***********************************************************************************************************/

    for (divisioncycle=start_divisioncycle+1; divisioncycle <= myParam.maxdivcycle; divisioncycle++) {
            switch (ALGORITHM) {
                case 0: //Muyoung
                    overflowflag = 0; count1=count2=0;
                    for(who=0; who<MAXORGANISMS; who++) {
                        if (myOrgstatus[who]==S_DEAD) continue;
                        if (myOrgstatus[who]==S_NEWBORN) continue;
                        if (myOrgstatus[who]==S_TODIE) continue;
                        if (myOrgstatus[who]==S_TODIE2) continue;
                        if (myOrgstatus[who]!=S_ALIVE) {fprintf(stderr,"bad org status\n"); exit(1);}
                        // only living ones are processed
                        
                        if (allow_gene_exp==1){
                            do { fate = (int) ((double) rand() / RAND_MAX * 3); }
                            while(fate == 3);
                            // randomly assigns a fate: 0, 1, or 2
                        }
                        else{
                            do { fate = (int) ((double) rand() / RAND_MAX * 2); }
                            while(fate == 2);
                            // randomly assigns a fate: 0, or 1
                        }
                        
                        
                        if (fate == 0) {
                            status = OrgChildBirth(&myParam, who, divisioncycle, 3);
                            if (status > 0) { count1 += status;
                            } else if(status == -9) {
                                overflowflag = 1;}
                        } else if (fate == 1) {status=OrgDeath(&myParam, who, 3);
                        } else if (fate == 2) {status=GeneExpress(&myParam, who, 3, &dC);
                            //myOrg[who].action_count[3] += sqrt(dC*dC);
                        }
                    }
                    
                    deathcount=death2count=maxorgind=newborncount=0;
                    count=0;
                    for(who=0;who<MAXORGANISMS;who++) {
                        if(myOrgstatus[who] == S_NEWBORN){ myOrgstatus[who] = S_ALIVE;newborncount++;orgcount++;}
                        if(myOrgstatus[who] == S_TODIE) { deathcount++; orgcount--; KillOrganism(who);}
                        if(myOrgstatus[who] == S_TODIE2) {death2count++;KillOrganism(who);}
                        if(myOrgstatus[who] == S_ALIVE) {count++;maxorgind=who;}
                    }
                    if(count!=orgcount) {
                        fprintf(stderr, "t=%d inconsistent orgcount %d count %d\n", divisioncycle, orgcount, count);
                        fprintf(stderr, "# of new born cells : %d %d\n", count1, newborncount);
                        fprintf(stderr, "# of killed cells : %d, %d\n", deathcount, death2count);
                        exit(1);
                    }
                    if(orgcount<=0) {fprintf(stdout, "died out at gen %d\n", divisioncycle); fprintf(out1, "died out at gen %d\n", divisioncycle);
                        //fclose(zout4);
                        exit(1);}
                    if(overflowflag) {fprintf(stdout, "organism overflow %d\n", orgcount); break;}
                    
                    /* decimation */
                    kia=0;
                    if(orgcount>myParam.decimthresh) {
                        count1=0;
                        do{
                            ii = (int)(((double)rand()/RAND_MAX)*(maxorgind+1));
                            if(myOrgstatus[ii]==S_ALIVE){KillOrganism(ii);myOrgstatus[ii]=S_DEAD;kia++;count--;}
                        } while(count>myParam.decimto);
                        orgcount=count;
                        lastdecimtime=divisioncycle;
                    }
                    
                    TIME = divisioncycle;
                    break;
                    
                    
                case 1: //Gillespie
                    generation = divisioncycle;
                    for (orgcount=0; orgcount < POPSIZE; orgcount++) {
                        cum_br[0]=myOrg[0].birthrate;
                        for(who=1; who<MAXORGANISMS; who++) {
                            cum_br[who] = cum_br[who-1] + myOrg[who].birthrate;
                        }
                        max_br=cum_br[MAXORGANISMS-1];
                        r1=( (double) rand()/RAND_MAX );
                        r2=( (double) rand()/RAND_MAX );
                        if (r1 < 1.0e-50){ r1 = 1.0e-50;}
                       // dt=-(1.0e0/max_br)*log(r1);
                        dt=-(1.0e0/max_br)*log(r1/max_br);
                        TIME=TIME+dt;
                        r2=r2*max_br;
                        for(who=0; who<MAXORGANISMS; who++) {
                            if ( r2 < cum_br[who]) { break ; }
                        }
                        if ( r2 >= max_br ){ who = MAXORGANISMS-1 ; }
                        
                        //find org to kill:
                        do{ ii = (int)( ((double) rand()/RAND_MAX) * MAXORGANISMS );
                        } while( ii==MAXORGANISMS || ii==who );
                        
                        //kill it and replace it
                        memmove(&myOrg[ii], &myOrg[who], sizeof(organism));
                        
                        if (allow_gene_exp==1){			//RMR sims
                            do{ kk = (int)( ((double) rand()/RAND_MAX) * 2 );
                            } while( kk==2 );
                            // randomly assigns a fate: 0, 1
                        }
                        else{kk = 0;} // randomly assigns a fate: 0 ...
                        
                        if(kk==0){
                            //printf("before OrgGeneMutate ... \n");
                            OrgGeneMutate(&myParam,who,2);
                        }
                        if(kk==1){
                            GeneExpress(&myParam,who,2, &dC);
                        }
                        
                        if (allow_gene_exp==1){
                            do{ jj = (int)( ((double) rand()/RAND_MAX) * 2 );
                            } while( jj==2 );
                            // randomly assigns a fate: 0, 1
                        }
                        else{jj = 0;} // randomly assigns a fate: 0 ...
                        
                        if(jj==0){
                            OrgGeneMutate(&myParam,ii,2);
                        }
                        if(jj==1){
                            GeneExpress(&myParam,ii,2, &dC);
                        }
                    }
                    
                    break;
                default:
                    break;
            }
            
            //output ADJUST
            PrepareOutput();             
            Flushfiles();
        }//end if time condition
    
    if (selection==0) WriteConfig();
    Closefiles();
    return 0;
}


/***********************************************************************************************************
 **********************************************************************************************************/
void PrepareOutput(){
    int ii, who, m, jj, j;
    
    ResetOrgDB(&myParam, divisioncycle);
    RankSpeciesSizeDB(sizeRank, nOrgDB);
	int genediff=0;
    for(j=0;j<myOrgDB[sizeRank[0]].genecount*NUCSEQLEN;j++) {
        // comparing if genome of who matches any genome in the database
        genediff+=(myOrgDB[sizeRank[0]].genome[j]-oldgenome[j]);
    }
    if(!genediff) return;

    for(j=0;j<myOrgDB[sizeRank[0]].genecount*NUCSEQLEN;j++) {
 		oldgenome[j] = myOrgDB[sizeRank[0]].genome[j];
	}
	
    /* statistics for organisms */
    m=count=mutatorcount=0;
    for(ii=0;ii<curr_MAXGENES;ii++) mutatororigin[ii]=0;
    SetMeanOrg(&mean);
	SetVarOrg(&var);
    for(who=0;who<POPSIZE;who++){
        if(myOrgstatus[who]!=S_ALIVE) continue;
        AddMeanOrg(&mean, who);
        if(myOrg[who].mutlevel == 1) mutatorcount++;
        mutatororigin[myOrg[who].mutorigin]++;
        for (ii=0; ii<curr_MAXGENES; ii++) seq[ii][m] = myOrg[who].genome+ii*NUCSEQLEN;
        
        m++;
    }
    GetMeanOrg(&mean,m);
    for(who=0;who<POPSIZE;who++){
        if(myOrgstatus[who]!=S_ALIVE) continue;
        AddVarOrg(&var, who);
    }
	GetVarOrg(&var, m);
	
    entropy_sumtot = 0.0;
    for(ii=0;ii<curr_MAXGENES;ii++) {
        seq_entropy(seq[ii], m, entropy[ii]);
        
        entropy_sum[ii]=0.0;
        for(jj=0;jj<AASEQLEN;jj++) entropy_sum[ii] += entropy[ii][jj];
        entropy_sum[ii] /= (double)AASEQLEN;
        entropy_sumtot += entropy_sum[ii];
    }
    entropy_sumtot /= (double)curr_MAXGENES;

	PrintOutput(); 
}

/***********************************************************************************************************
 **********************************************************************************************************/
void WriteConfig(){
    int i, j, who;
	char rootdir[100];
	char fopbuf[100]; 
   
    sprintf(rootdir,"/n/regal/shakhnovich_lab/rrazban/");
    sprintf(fopbuf,"%s/_config.txt", rootdir);
    FILE *fp=fopen(fopbuf, "w");
    
    // variables from evo-cell
    fprintf(fp,"%d %d %d %d %d %d %d %d\n", divisioncycle, lastdecimtime, orgcount, mutatorcount, generation, speciesSizeTotal, domi_species, nOrgDB);
    fprintf(fp,"%E\n", TIME);
    
    
    // variables from cell.h: myParam, myOrg
    fprintf(fp,"%d %d %d %d %d %d %d %d %d %d %d %d\n",myParam.seed, myParam.startcode, myParam.orgcount, myParam.maxdivcycle, myParam.decimthresh, myParam.decimto, myParam.initpop, myParam.dumpcycle, myParam.printoutcycle, myParam.plotoutcycle, myParam.screenoutcycle, myParam.seqlogcycle);
    fprintf(fp,"%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", myParam.Tenv, myParam.tol, myParam.birthrate, myParam.deathrate, myParam.expressrate, myParam.alpha, myParam.pnatcutoff, myParam.speciessizefactor, myParam.mutrate[0], myParam.mutrate[1], myParam.mutrate0, myParam.x0, myParam.b0, myParam.fixed_mutrate, myParam.mutthreshfac, myParam.mutthresh, myParam.bindindex); //19
    //fprintf(fp,"%s\n", myParam.input);
    fprintf(fp,"%s\n", myParam.targetname);
    
    for (i=0; i<7; i++){
        fprintf(fp,"%d ", mutatororigin[i]);
    }
    fprintf(fp,"\n");
    
    for (who=0; who<MAXORGANISMS; who++){
        fprintf(fp,"%d %d\n", sizeRank[who], myOrgstatus[who]);
    }
    
    // variables from cell.h: myOrg
    for (who=0; who<MAXORGANISMS; who++){
        fprintf(fp,"%d %d %d %d %d %d %d\n", myOrg[who].genecount, myOrg[who].mutlevel, myOrg[who].mutorigin, myOrg[who].ppicount, myOrg[who].dob, myOrg[who].numkids, myOrg[who].generation);
        fprintf(fp,"%lf %lf %lf %lf\n", myOrg[who].minpnat, myOrg[who].meanpnat, myOrg[who].G0, myOrg[who].birthrate);
        fprintf(fp,"%f %f\n",myOrg[who].mutcount, myOrg[who].mutrate);
        
        for (i=0; i<(MAXGENES*NUCSEQLEN); i++){
            fprintf(fp,"%d ",myOrg[who].genome[i]); //MAXGENES*NUCSEQLEN;
        }
        fprintf(fp,"\n");
        
        
        for (i=0; i<(MAXGENES+1); i++){
            fprintf(fp,"%d ",myOrg[who].structid[i]);
        }
        fprintf(fp,"\n");
        
        for (i=0; i<(MAXGENES*(MAXGENES-1)/2); i++){
            fprintf(fp,"%d ",myOrg[who].SeqID[i]);
        }
        fprintf(fp,"\n");
        
        for (i=0; i<MAXPPIS; i++){
            fprintf(fp,"%d ",myOrg[who].ppi_pair[i][0]);
        }
        fprintf(fp,"\n");
        for (i=0; i<MAXPPIS; i++){
            fprintf(fp,"%d ",myOrg[who].ppi_pair[i][1]);
        }
        fprintf(fp,"\n");
        
        for (i=0; i<MAXPPIS; i++){
            fprintf(fp,"%d ",myOrg[who].bmode[i]);
        }
        fprintf(fp,"\n");
        
        for (i=0; i<(MAXSTATES+1); i++){
            fprintf(fp,"%lf ",myOrg[who].F[i]);
        }
        fprintf(fp,"\n");
        
        for (i=0; i<(MAXSTATES+1); i++){
            for (j=0; j<(MAXSTATES+1); j++){
                fprintf(fp,"%lf ",myOrg[who].K[i][j]);
            }
            fprintf(fp,"\n");
        }
        
        for (i=0; i<(MAXGENES+1); i++){
            fprintf(fp,"%lf ",myOrg[who].C[i]);
        }
        fprintf(fp,"\n");
        
        for (i=0; i<(MAXGENES+1); i++){
            fprintf(fp,"%lf ",myOrg[who].Kc[i]);
        }
        fprintf(fp,"\n");
        
        for (i=0; i<(MAXGENES+1); i++){
            fprintf(fp,"%lf ",myOrg[who].Nch[i]);
        }
        fprintf(fp,"\n");
        
        for (i=0; i<(MAXGENES); i++){
            fprintf(fp,"%f %f %f %f %f %f %f %f %f %f %f\n",
                    myOrg[who].pnat[i],
                    myOrg[who].hydro[i],
                    myOrg[who].netcharge[i],
                    myOrg[who].poscharge[i],
                    myOrg[who].negcharge[i],
                    myOrg[who].fraccharge[i],
                    myOrg[who].frachydro[i],
                    myOrg[who].frachydro_avil[i],
                    myOrg[who].synmut[i],
                    myOrg[who].nonsynmut[i],
                    myOrg[who].tot_mut[i]);
        }

        for (i=0; i<(MAXPPIS); i++){
            fprintf(fp,"%f %f %f %f %f %f %f %f %f %f %f %f %f %f\n",
                    myOrg[who].hydro_s_hub[i],
                    myOrg[who].netcharge_s_hub[i],
                    myOrg[who].poscharge_s_hub[i],
                    myOrg[who].negcharge_s_hub[i],
                    myOrg[who].fraccharge_s_hub[i],
                    myOrg[who].frachydro_s_hub[i],
                    myOrg[who].frachydro_avil_s_hub[i],
                    
                    myOrg[who].hydro_s_partner[i],
                    myOrg[who].netcharge_s_partner[i],
                    myOrg[who].poscharge_s_partner[i],
                    myOrg[who].negcharge_s_partner[i],
                    myOrg[who].fraccharge_s_partner[i],
                    myOrg[who].frachydro_s_partner[i],
                    myOrg[who].frachydro_avil_s_partner[i]);
        }

        
        
        for (i=0; i<(MAXGENES+1); i++){
            fprintf(fp,"%f %f %f %f\n",
                    myOrg[who].nsi[i],
                    myOrg[who].si[i],
                    myOrg[who].nsi2[i],
                    myOrg[who].si2[i]);
        }
        
        
        for (i=0; i<(MAXPPIS+1); i++){
            fprintf(fp,"%f ",myOrg[who].pint[i]);
        }
        fprintf(fp,"\n");
        
        for (i=0; i<(MAXPPIS+1); i++){
            fprintf(fp,"%f ",myOrg[who].Gij[i]);
        }
        fprintf(fp,"\n");
        fprintf(fp,"\n");
    }
    
    
    fclose(fp);
}

/***********************************************************************************************************
 fscanf(fp,"%s %s %s %s %s %s\n",t_string, step_count_string, stretching_string,
 output_count_string,initial_L_string,strain_stoped_once_string );
 t = atof(t_string);
 output_count=atoi(output_count_string);
 step_count=atoi(step_count_string);
 stretching=atoi(stretching_string);
 initial_L = atof(initial_L_string);
 strain_stoped_once = atoi(strain_stoped_once_string);
 **********************************************************************************************************/
void ReadConfig(){
    int i, j, who;
    FILE *fp;
    char  divisioncycle_s[100], lastdecimtime_s[100], orgcount_s[100], mutatorcount_s[100], generation_s[100],
    speciesSizeTotal_s[100], domi_species_s[100], nOrgDB_s[100], TIME_s[100];
    char  seed_s[100], startcode_s[100], mp_orgcount_s[100], maxdivcycle_s[100], decimthresh_s[100], decimto_s[100], initpop_s[100], dumpcycle_s[100], printoutcycle_s[100], plotoutcycle_s[100], screenoutcycle_s[100], seqlogcycle_s[100];
    char Tenv_s[100], tol_s[100], birthrate_s[100], deathrate_s[100], expressrate_s[100], alpha_s[100], pnatcutoff_s[100], speciessizefactor_s[100], mutrate0_s[100], mutrate01_s[100], mutrate00_s[100], x0_s[100], b0_s[100], fixed_mutrate_s[100], mutthreshfac_s[100], mutthresh_s[100], bindindex_s[100];
    //char targetname_s[500];
    char mutatororigin_i_s[100];
    char sizeRank_who_s[100], myOrgstatus_who_s[100];
    char mo_genecount_s[100], mo_mutlevel[100], mo_mutorigin[100], mo_ppicount[100], mo_dob[100], mo_numkids[100], mo_generation[100], mo_minpnat[100], mo_meanpnat[100], mo_G0[100], mo_birthrate[100], mo_mutcount[100], mo_mutrate[100], mo_genome_i[100];
    char mo_structid_i_s[100], mo_SeqID_i_s[100], mo_ppi_pair0_i_s[100], mo_ppi_pair1_i_s[100], mo_bmode_i_s[100], mo_F_i_s[100], mo_K_ij_s[100], mo_C_i_s[100], mo_Kc_i_s[100], mo_Nch_i_s[100];
    char mo_pnat_i_s[100], mo_hydro_i_s[100], mo_netcharge_i_s[100], mo_poscharge_i_s[100], mo_negcharge_i_s[100], mo_fraccharge_i_s[100], mo_frachydro_i_s[100], mo_frachydro_avil_i_s[100], mo_synmut_i_s[100], mo_nonsynmut_i_s[100], mo_tot_mut_i_s[100], mo_nsi_i_s[100], mo_si_i_s[100], mo_nsi2_i_s[100], mo_si2_i_s[100], mo_pint_i_s[100], mo_Gij_ij_s[100];
    char mo_hydro_i_s2[100], mo_netcharge_i_s2[100], mo_poscharge_i_s2[100], mo_negcharge_i_s2[100], mo_fraccharge_i_s2[100], mo_frachydro_i_s2[100], mo_frachydro_avil_i_s2[100]; 
    
    
    
    
    fp=fopen(FILE_CONFIG, "rb+");
    if(fp == NULL) { //if file does not exist
        RUN_BASED_ON_CONFIG = 0;
    }else {
        RUN_BASED_ON_CONFIG = 1;
    }
    printf("\nRUN_BASED_ON_CONFIG = %d\n",RUN_BASED_ON_CONFIG);
    
    if (RUN_BASED_ON_CONFIG){
        printf("Reading Config file ... ");
        
        // variables from evo-cell
        fscanf(fp,"%s %s %s %s %s %s %s %s\n",divisioncycle_s, lastdecimtime_s, orgcount_s, mutatorcount_s, generation_s,
               speciesSizeTotal_s, domi_species_s, nOrgDB_s);
        divisioncycle = atoi(divisioncycle_s);
        lastdecimtime = atoi(lastdecimtime_s);
        orgcount = atoi(orgcount_s);
        mutatorcount = atoi(mutatorcount_s);
        generation = atoi(generation_s);
        speciesSizeTotal = atoi(speciesSizeTotal_s);
        domi_species = atoi(domi_species_s);
        nOrgDB = atoi(nOrgDB_s);
        
        fscanf(fp,"%s\n",TIME_s);
        TIME = atof(TIME_s);
        
        
        printf("divisioncycle = %d\n",divisioncycle);
        printf("TIME = %E\n",TIME);
        //exit(0);
        
        // variables from cell.h: myParam, myOrg
        
        
        fscanf(fp,"%s %s %s %s %s %s %s %s %s %s %s %s\n",seed_s,  startcode_s,  mp_orgcount_s,  maxdivcycle_s,  decimthresh_s,  decimto_s,  initpop_s,  dumpcycle_s,  printoutcycle_s,  plotoutcycle_s,  screenoutcycle_s,  seqlogcycle_s);
        
        myParam.seed = atoi(seed_s);
        myParam.startcode = atoi(startcode_s);
        myParam.orgcount = atoi(orgcount_s);
        myParam.maxdivcycle = atoi(maxdivcycle_s);	//gets overridden
        myParam.decimthresh = atoi(decimthresh_s);
        myParam.decimto = atoi(decimto_s);
        myParam.initpop = atoi(initpop_s);
        myParam.dumpcycle = atoi(dumpcycle_s);
        myParam.printoutcycle = atoi(printoutcycle_s);
        myParam.plotoutcycle = atoi(plotoutcycle_s);
        myParam.screenoutcycle = atoi(screenoutcycle_s);
        myParam.seqlogcycle = atoi(seqlogcycle_s);
        
        //exit(0);
        
        fscanf(fp,"%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n",Tenv_s,  tol_s, birthrate_s,  deathrate_s,  expressrate_s,  alpha_s,  pnatcutoff_s,  speciessizefactor_s,  mutrate0_s,  mutrate01_s,  mutrate00_s,  x0_s,  b0_s,  fixed_mutrate_s,  mutthreshfac_s,  mutthresh_s,  bindindex_s);
        
        myParam.Tenv = atof(Tenv_s);
        myParam.tol = atof(tol_s);
        myParam.birthrate = atof(birthrate_s);
        myParam.deathrate = atof(deathrate_s);
        myParam.expressrate = atof(expressrate_s);
        myParam.alpha = atof(alpha_s);
        myParam.pnatcutoff = atof(pnatcutoff_s);
        myParam.speciessizefactor = atof(speciessizefactor_s);
        myParam.mutrate[0] = atof(mutrate00_s);
        myParam.mutrate[1] = atof(mutrate01_s);
        myParam.mutrate0 = atof(mutrate0_s);
        myParam.x0 = atof(x0_s);
        myParam.b0 = atof(b0_s);
        myParam.fixed_mutrate = atof(fixed_mutrate_s);
        myParam.mutthreshfac = atof(mutthreshfac_s);
        myParam.mutthresh = atof(mutthresh_s);
        myParam.bindindex  = atof(bindindex_s);
        
        
        //printf("myParam.b0  = %lf\n",myParam.b0 );
        //exit(0);
        
        fscanf(fp,"%s\n", myParam.targetname);
        //printf("myParam.targetname  = %s\n",myParam.targetname);
        //exit(0);
        
        
        for (i=0; i<7; i++){
            fscanf(fp,"%s ", mutatororigin_i_s);
            mutatororigin[i] = atoi(mutatororigin_i_s);
            //printf("mutatororigin[i]   = %d\n",mutatororigin[i] );
            
        }
        //exit(0);
        fscanf(fp,"\n");
        
        for (who=0; who<MAXORGANISMS; who++){
            fscanf(fp,"%s %s\n", sizeRank_who_s, myOrgstatus_who_s);
            sizeRank[who] = atoi(sizeRank_who_s);
            myOrgstatus[who] = atoi(myOrgstatus_who_s);
            //printf("sizeRank  = %d myOrgstatus=%d\n",sizeRank[who] , myOrgstatus[who] );
        }
        
        //exit(0);
        
        // variables from cell.h: myOrg
        for (who=0; who<MAXORGANISMS; who++){
            fscanf(fp,"%s %s %s %s %s %s %s\n", mo_genecount_s, mo_mutlevel, mo_mutorigin, mo_ppicount, mo_dob, mo_numkids, mo_generation);
            myOrg[who].genecount = atoi(mo_genecount_s);
            myOrg[who].mutlevel = atoi(mo_mutlevel);
            myOrg[who].mutorigin = atoi(mo_mutorigin);
            myOrg[who].ppicount = atoi(mo_ppicount);
            myOrg[who].dob = atoi(mo_dob);
            myOrg[who].numkids = atoi(mo_numkids);
            myOrg[who].generation = atoi(mo_generation);
            
            //printf("generation  = %d\n",generation);
            //exit(0);
            
            fscanf(fp,"%s %s %s %s\n", mo_minpnat, mo_meanpnat, mo_G0, mo_birthrate);
            myOrg[who].minpnat = atof(mo_minpnat);
            myOrg[who].meanpnat = atof(mo_meanpnat);
            myOrg[who].G0 = atof(mo_G0);
            myOrg[who].birthrate = atof(mo_birthrate);
            //printf("myOrg[who].meanpnat  = %lf\n",myOrg[who].meanpnat);
            //printf("myOrg[who].birthrate  = %lf\n",myOrg[who].birthrate);
            //exit(0);
            
            
            fscanf(fp,"%s %s\n",mo_mutcount, mo_mutrate);
            myOrg[who].mutcount = atof(mo_mutcount);
            myOrg[who].mutrate = atof(mo_mutrate);
            
            //printf("myOrg[who].mutrate  = %f\n",myOrg[who].mutrate);
            //exit(0);
            
            for (i=0; i<(MAXGENES*NUCSEQLEN); i++){
                fscanf(fp,"%s ",mo_genome_i);
                myOrg[who].genome[i] = atoi(mo_genome_i);
                //printf("%s", mo_genome_i); //empty.....
                //printf("%d ", myOrg[who].genome[i]);
            }
            //printf("\n");
            fscanf(fp,"\n");
            //exit(0);
            
            for (i=0; i<(MAXGENES+1); i++){
                fscanf(fp,"%s ",mo_structid_i_s);
                myOrg[who].structid[i] = atoi(mo_structid_i_s);
            }
            fscanf(fp,"\n");
            
            for (i=0; i<(MAXGENES*(MAXGENES-1)/2); i++){
                fscanf(fp,"%s ",mo_SeqID_i_s);
                myOrg[who].SeqID[i] = atoi(mo_SeqID_i_s);
            }
            fscanf(fp,"\n");
            
            for (i=0; i<MAXPPIS; i++){
                fscanf(fp,"%s ",mo_ppi_pair0_i_s);
                myOrg[who].ppi_pair[i][0] = atoi(mo_ppi_pair0_i_s);
            }
            fscanf(fp,"\n");
            for (i=0; i<MAXPPIS; i++){
                fscanf(fp,"%s ",mo_ppi_pair1_i_s);
                myOrg[who].ppi_pair[i][1] = atoi(mo_ppi_pair1_i_s);
            }
            fscanf(fp,"\n");
            
            for (i=0; i<MAXPPIS; i++){
                fscanf(fp,"%s ",mo_bmode_i_s);
                myOrg[who].bmode[i] = atoi(mo_bmode_i_s);
            }
            fscanf(fp,"\n");
            
            for (i=0; i<(MAXSTATES+1); i++){
                fscanf(fp,"%s ",mo_F_i_s);
                myOrg[who].F[i] = atof(mo_F_i_s);
            }
            fscanf(fp,"\n");
            
            for (i=0; i<(MAXSTATES+1); i++){
                for (j=0; j<(MAXSTATES+1); j++){
                    fscanf(fp,"%s ",mo_K_ij_s);
                    myOrg[who].K[i][j] = atof(mo_K_ij_s);
                }
                fscanf(fp,"\n");
            }
            
            for (i=0; i<(MAXGENES+1); i++){
                fscanf(fp,"%s ",mo_C_i_s);
                myOrg[who].C[i] = atof(mo_C_i_s);
            }
            fscanf(fp,"\n");
            
            for (i=0; i<(MAXGENES+1); i++){
                fscanf(fp,"%s ",mo_Kc_i_s);
                myOrg[who].Kc[i] = atof(mo_Kc_i_s);
            }
            fscanf(fp,"\n");
            
            for (i=0; i<(MAXGENES+1); i++){
                fscanf(fp,"%s ",mo_Nch_i_s);
                myOrg[who].Nch[i] = atof(mo_Nch_i_s);
            }
            fscanf(fp,"\n");
            
            for (i=0; i<(MAXGENES); i++){
                fscanf(fp,"%s %s %s %s %s %s %s %s %s %s %s\n",mo_pnat_i_s, mo_hydro_i_s, mo_netcharge_i_s, mo_poscharge_i_s, mo_negcharge_i_s, mo_fraccharge_i_s, mo_frachydro_i_s, mo_frachydro_avil_i_s, mo_synmut_i_s, mo_nonsynmut_i_s, mo_tot_mut_i_s);
                
                myOrg[who].pnat[i] = atof(mo_pnat_i_s);
                myOrg[who].hydro[i]= atof(mo_hydro_i_s);
                myOrg[who].netcharge[i]= atof(mo_netcharge_i_s);
                myOrg[who].poscharge[i]= atof(mo_poscharge_i_s);
                myOrg[who].negcharge[i]= atof(mo_negcharge_i_s);
                myOrg[who].fraccharge[i]= atof(mo_fraccharge_i_s);
                myOrg[who].frachydro[i]= atof(mo_frachydro_i_s);
                myOrg[who].frachydro_avil[i]= atof(mo_frachydro_avil_i_s);
                myOrg[who].synmut[i]= atof(mo_synmut_i_s);
                myOrg[who].nonsynmut[i]= atof(mo_nonsynmut_i_s);
                myOrg[who].tot_mut[i]= atof(mo_tot_mut_i_s);
            }
            
            for (i=0; i<(MAXPPIS); i++){
                fscanf(fp,"%s %s %s %s %s %s %s %s %s %s %s %s %s %s\n", mo_hydro_i_s, mo_netcharge_i_s, mo_poscharge_i_s, mo_negcharge_i_s, mo_fraccharge_i_s, mo_frachydro_i_s, mo_frachydro_avil_i_s, mo_hydro_i_s2, mo_netcharge_i_s2, mo_poscharge_i_s2, mo_negcharge_i_s2, mo_fraccharge_i_s2, mo_frachydro_i_s2, mo_frachydro_avil_i_s2);
                
                myOrg[who].hydro_s_hub[i]= atof(mo_hydro_i_s);
                myOrg[who].netcharge_s_hub[i]= atof(mo_netcharge_i_s);
                myOrg[who].poscharge_s_hub[i]= atof(mo_poscharge_i_s);
                myOrg[who].negcharge_s_hub[i]= atof(mo_negcharge_i_s);
                myOrg[who].fraccharge_s_hub[i]= atof(mo_fraccharge_i_s);
                myOrg[who].frachydro_s_hub[i]= atof(mo_frachydro_i_s);
                myOrg[who].frachydro_avil_s_hub[i]= atof(mo_frachydro_avil_i_s);
                
                myOrg[who].hydro_s_partner[i]= atof(mo_hydro_i_s2);
                myOrg[who].netcharge_s_partner[i]= atof(mo_netcharge_i_s2);
                myOrg[who].poscharge_s_partner[i]= atof(mo_poscharge_i_s2);
                myOrg[who].negcharge_s_partner[i]= atof(mo_negcharge_i_s2);
                myOrg[who].fraccharge_s_partner[i]= atof(mo_fraccharge_i_s2);
                myOrg[who].frachydro_s_partner[i]= atof(mo_frachydro_i_s2);
                myOrg[who].frachydro_avil_s_partner[i]= atof(mo_frachydro_avil_i_s2);
            }

            
            
            for (i=0; i<(MAXGENES+1); i++){
                fscanf(fp,"%s %s %s %s\n",mo_nsi_i_s, mo_si_i_s, mo_nsi2_i_s, mo_si2_i_s);
                myOrg[who].nsi[i] = atof(mo_nsi_i_s);
                myOrg[who].si[i] = atof(mo_si_i_s);
                myOrg[who].nsi2[i] = atof(mo_nsi2_i_s);
                myOrg[who].si2[i] = atof(mo_si2_i_s);
            }
            
            
            for (i=0; i<(MAXPPIS+1); i++){
                fscanf(fp,"%s ",mo_pint_i_s);
                myOrg[who].pint[i] = atof(mo_pint_i_s);
            }
            fscanf(fp,"\n");
            
            for (i=0; i<(MAXPPIS+1); i++){
                fscanf(fp,"%s ",mo_Gij_ij_s);
                myOrg[who].Gij[i] = atof(mo_Gij_ij_s);
            }
            fscanf(fp,"\n");
            fscanf(fp,"\n");
        }
        
        
        printf("Done\n\n");
    }
    //exit(0);
}


/***********************************************************************************************************
 **********************************************************************************************************/
void PrintOutput(){
	int ii;
//	int jj;

        
/*	fprintf(out1,"%08d %d",divisioncycle, orgcount);
	for (ii=0; ii<curr_MAXGENES; ii++) fprintf(out1, " %E", mean.C[ii]);
	for (ii=0; ii<curr_MAXSTATES; ii++) fprintf(out1, " %E", mean.F[ii]);
	for (ii=0; ii<curr_MAXGENES; ii++) fprintf(out1, " %.10E", mean.pnat[ii]);
	for (ii=0; ii<curr_MAXPPIS; ii++) fprintf(out1, " %E", mean.pint[ii]);
	for (ii=0; ii<curr_MAXPPIS; ii++) fprintf(out1, " %E", mean.Gij[ii]);
        
	for (ii=0; ii<curr_MAXGENES; ii++) fprintf(out1, " %06d", myOrgDB[sizeRank[0]].structid[ii]);
	for (ii=0; ii<curr_MAXPPIS; ii++) fprintf(out1, " %03d", myOrgDB[sizeRank[0]].bmode[ii]);
        
	int aaseq_surface_hub[AASURFACELEN], aaseq_surface_partner[AASURFACELEN], s_i, hub_i, par_i;
        
	for (ii=0;ii<curr_MAXPPIS; ii++){
		hub_i=0; par_i=ii+1;
		GetSurfaceAAPositions(myOrgDB[sizeRank[0]].structid[hub_i], myOrgDB[sizeRank[0]].structid[par_i], myOrgDB[sizeRank[0]].bmode[ii],
                                  aaseq_surface_hub, aaseq_surface_partner);
            
		for (s_i=0;s_i<AASURFACELEN; s_i++) fprintf(out1, " %02d", aaseq_surface_hub[s_i]);	
		for (s_i=0;s_i<AASURFACELEN; s_i++) fprintf(out1, " %02d", aaseq_surface_partner[s_i]);
	}
        
	PrintCharNucCodeSequence(seqbuf, myOrgDB[sizeRank[0]].genome, myOrgDB[sizeRank[0]].genecount*NUCSEQLEN);
	fprintf(out1, " %s ", seqbuf);
	fprintf(out1,"\n");

        
	fprintf(out22," %d", divisioncycle);
	fprintf(out22,"%.3E",TIME);
    for (ii=0; ii<curr_MAXGENES; ii++) fprintf(out22, " %.3E", mean.C[ii]);
	for (ii=0; ii<curr_MAXSTATES; ii++) fprintf(out22, " %.3E", mean.F[ii]);
	for (ii=0; ii<curr_MAXSTATES; ii++) for (jj=ii; jj<curr_MAXSTATES; jj++) fprintf(out22, " %.3E", mean.nF[ii][jj]);
	for (ii=0; ii<curr_MAXGENES; ii++) fprintf(out22, " %.3E", mean.pnat[ii]);
	for (ii=0; ii<curr_MAXPPIS; ii++) fprintf(out22, " %.3E", mean.pint[ii]);
	for (ii=0; ii<curr_MAXPPIS; ii++) fprintf(out22, " %.3E", mean.Gij[ii]);
	fprintf(out22," %.3E",(double) domi_species/orgcount);
	fprintf(out22,"\n");
        
	fprintf(out23," %d", divisioncycle);
	fprintf(out23,"%.3E",TIME);
	for (ii=0; ii<curr_MAXGENES; ii++) fprintf(out23, " %.3E", var.C[ii]);
	for (ii=0; ii<curr_MAXSTATES; ii++) fprintf(out23, " %.3E", var.F[ii]);
	for (ii=0; ii<curr_MAXSTATES; ii++) for (jj=ii; jj<curr_MAXSTATES; jj++) fprintf(out23, " %.3E", var.nF[ii][jj]);
	for (ii=0; ii<curr_MAXGENES; ii++) fprintf(out23, " %.3E", var.pnat[ii]);
	for (ii=0; ii<curr_MAXPPIS; ii++) fprintf(out23, " %.3E", var.pint[ii]);
	for (ii=0; ii<curr_MAXPPIS; ii++) fprintf(out23, " %.3E", var.Gij[ii]);
	fprintf(out23,"\n");
        
       
	fprintf(out24," %d", divisioncycle);
	PrintCharNucCodeSequence(seqbuf, myOrgDB[sizeRank[0]].genome, myOrgDB[sizeRank[0]].genecount*NUCSEQLEN);
	fprintf(out24, " %s", seqbuf);
	fprintf(out24,"\n");
*/
	int jj;
	int aaseq[AASEQLEN], aaseq2[AASEQLEN];
    CharNucSeqToAASeq(myOrgDB[sizeRank[0]].genome+0*NUCSEQLEN,NUCSEQLEN,aaseq);
    CharNucSeqToAASeq(myOrgDB[sizeRank[0]].genome+1*NUCSEQLEN,NUCSEQLEN,aaseq2);


	fprintf(out9,"%5d", divisioncycle);
	fprintf(out9, " %.3E", SequenceEnergy(aaseq, myOrgDB[sizeRank[0]].structid[0]));
	fprintf(out9, " %.3E", SequenceEnergy(aaseq2, myOrgDB[sizeRank[0]].structid[1]));
	fprintf(out9, " %.3E", GetBindingEnergy2(aaseq, myOrgDB[sizeRank[0]].structid[0], aaseq2, myOrgDB[sizeRank[0]].structid[1], myOrgDB[sizeRank[0]].bmode[0]));	
	for (ii=0; ii<curr_MAXGENES; ii++) fprintf(out9, " %.3E", myOrgDB[sizeRank[0]].pnat[ii]);
	for (ii=0; ii<curr_MAXPPIS; ii++) fprintf(out9, " %.3E", myOrgDB[sizeRank[0]].pint[ii]);
    for (ii=0; ii<curr_MAXGENES; ii++) fprintf(out9, " %.3E", myOrgDB[sizeRank[0]].C[ii]);
	for (ii=0; ii<curr_MAXSTATES; ii++) fprintf(out9, " %.3E", myOrgDB[sizeRank[0]].F[ii]);
	for (ii=0; ii<curr_MAXSTATES; ii++) for (jj=ii; jj<curr_MAXSTATES; jj++) fprintf(out9, " %.3E", myOrgDB[sizeRank[0]].nF[ii][jj]);
	for (ii=0; ii<curr_MAXPPIS; ii++) fprintf(out9, " %.3E", myOrgDB[sizeRank[0]].Gij[ii]);
	fprintf(out9," %.3E", myOrgDB[sizeRank[0]].birthrate);
	fprintf(out9," %.2f",(double) domi_species/orgcount);
//	PrintCharNucCodeSequence(seqbuf, myOrgDB[sizeRank[0]].genome, myOrgDB[sizeRank[0]].genecount*NUCSEQLEN);
	PrintAASequence(seqbuf, aaseq, AASEQLEN);
	fprintf(out9, " %s", seqbuf);
	PrintAASequence(seqbuf, aaseq2, AASEQLEN);
	fprintf(out9, " %s", seqbuf);
	fprintf(out9,"\n");
}

void PrintHeaders(FILE *out){
	int width=9;   	

	fprintf(out, "%s ","#Time");
	fprintf(out, "%*s ",width,"E1nat");
	fprintf(out, "%*s ",width,"E2nat");
	fprintf(out, "%*s ",width,"Eint");
	fprintf(out, "%*s ",width,"P1nat");
	fprintf(out, "%*s ",width,"P2nat");
	fprintf(out, "%*s ",width,"Pint");
	fprintf(out, "%*s ",width,"C1");
	fprintf(out, "%*s ",width,"C2");
	fprintf(out, "%*s ",width,"F1");
	fprintf(out, "%*s ",width,"F2");
	fprintf(out, "%*s ",width,"F11");
	fprintf(out, "%*s ",width,"F12");
	fprintf(out, "%*s ",width,"F22");
	fprintf(out, "%*s",width,"Gf");
	fprintf(out, "%*s",width,"b");
	fprintf(out, "%*s ",width,"fraction");
	fprintf(out, "%*s",width*2,"seq1");
	fprintf(out, "%*s\n",width*2,"seq2");

}


void PrintHeaders1(FILE *out){
	int width=9;   	
	
	fprintf(out, "%*s ",width,"Dcycle");
	fprintf(out, "%*s\n",width,"RNA");
}
/***********************************************************************************************************
 **********************************************************************************************************/
void Openfiles(){
    char fopbuf[100];
    char filetype_buf[100];
    char rootdir[100];

    if (RUN_BASED_ON_CONFIG == 0){
        sprintf(filetype_buf,"w");
    }
    else {
        sprintf(filetype_buf,"a");
    }
   
    sprintf(rootdir,"/n/regal/shakhnovich_lab/rrazban/");
/*    sprintf(fopbuf,"%s/initial.dat", rootdir);
    out1=fopen(fopbuf,filetype_buf);
    PrintInitialCondition(out1,&myParam);
    fclose(out1);
    
    sprintf(fopbuf,"%s/basiclog.dat", rootdir);
    out1=fopen(fopbuf,filetype_buf);
    
    sprintf(fopbuf,"%s/phenotype.avg.dat", rootdir);
    out22=fopen(fopbuf,filetype_buf);
 	fprintf(out22,"Mean\n");
	PrintHeaders(out22);

	sprintf(fopbuf,"%s/phenotype.var.dat", rootdir);
    out23=fopen(fopbuf,filetype_buf);
 	fprintf(out23,"Variance\n");
	PrintHeaders(out23);

	sprintf(fopbuf,"%s/genotype.dat", rootdir);
	out24=fopen(fopbuf,filetype_buf);
 	fprintf(out24,"Absolute\n");
	PrintHeaders1(out24);
*/
	sprintf(fopbuf,"%s/%s", rootdir, myParam.targetname);
	out9=fopen(fopbuf,filetype_buf);
	fprintf(out9, "#Output only contains timepoints in which the most prevalent sequence among organisms in the population changes. Thus, no output between attempted mutation 0 and 11 means that the most prevalent sequence remains the same and populations from attempted mutation 0-10 have most prevalent sequences identical to that outputted in attempted mutation 0\n");
	fprintf(out9, "#Outputted values are only for one of the organisms which has the most prevalent sequence in the population\n");
	fprintf(out9, "#Labels follow Orit's 2014 paper except that F stands for free concentration and not folded concentration\n");
//	fprintf(out9, "#Time | Enat0 \t| Enat1 \t| Eint \t| Pnat0 \t| Pnat1 \t| Pint \t| C0 \t| C1 \t| F0 \t| F1 \t| F0F0 \t| F1F1 \t| Gnf \t| Gf \t| b \t| Fraction of population with most likeley sequence | Amino acid sequence of protein 0 | Amino acid sequence of protein 1\n");
	PrintHeaders(out9);
}

/***********************************************************************************************************
 **********************************************************************************************************/
void Flushfiles(){
//    fflush(out1);
    fflush(stdout);
//	fflush(it_solver);
//	fflush(out22); fflush(out23); fflush(out24);
	fflush(out9);
}

/***********************************************************************************************************
 **********************************************************************************************************/
void Closefiles(){
  //  fclose(out1);
//	fclose(it_solver);
//	fclose(out22); fclose(out23); fclose(out24);
	fclose(out9);
}


/***********************************************************************************************************
 // calculates entropy at each position of an amino acid
 // based on distribution among organisms
 **********************************************************************************************************/
int seq_entropy(char **seq, int N, double *entropy){
    int i, j;
    double f;
    
    for(j=0;j<AASEQLEN;j++)
        for(i=0;i<20;i++)
            aamat[j][i]=0;
    
    for(i=0;i<N;i++){   // N: number of copies of amino acid
        CharNucSeqToAASeq(seq[i], NUCSEQLEN, aaseq);
        for(j=0;j<AASEQLEN;j++){
            aamat[j][aaseq[j]]++;
            // counts number of each amino acid at each position
        }
    }
    
    for(j=0;j<AASEQLEN;j++){
        entropy[j]=0.0;
        //fprintf(stdout, "%4d :", j);
        for(i=0;i<20;i++){
            if(aamat[j][i] > 0) {
                f = (double) aamat[j][i]/N;
                entropy[j]-= f*logf(f);
            } else {
                f =0.0;
            }
            //fprintf(stdout, " %8.3f", f);
        }
        //fprintf(stdout,"\n");
    }
    return 0;
}


/***********************************************************************************************************
 **********************************************************************************************************/

// Updates database of all organism types (i.e. distinct genomes)
// in which each distinct genome gets its own space
int ResetOrgDB(parameter *myParam, int divisioncycle){
    int i, j, jj,genediff, who;
    //int face1, face2, rotate;
    //int surfacetmp[9],surface1[9], surface2[9];
    //double e, z, emin=1e10;
    //double T = myParam->Tenv;
    
    nOrgDB=0;
    domi_species = 0;
    
    for(who=0;who<MAXORGANISMS;who++) {
        if (myOrgstatus[who]!=S_ALIVE) continue;
        
        for(i=0;i<nOrgDB;i++) {
            if(myOrgDB[i].genecount != myOrg[who].genecount) continue;
            genediff=0;
            for(j=0;j<myOrg[who].genecount*NUCSEQLEN;j++) {
                // comparing if genome of who matches any genome in the database
                genediff+=(myOrgDB[i].genome[j]-myOrg[who].genome[j])*(myOrg[i].genome[j]-myOrg[who].genome[j]);
            }
            if(!genediff) break;
        }
        if(i == nOrgDB) {
            // if the genome does not match any previous entry, a new space is
            // made for it in OrgDB
            myOrgDB[i].genecount=myOrg[who].genecount;
            for(j=0;j<myOrg[who].genecount*NUCSEQLEN;j++) myOrgDB[i].genome[j]=myOrg[who].genome[j];
            for(j=0;j<myOrg[who].genecount;j++) {
                myOrgDB[i].structid[j]=myOrg[who].structid[j];
                myOrgDB[i].pnat[j]=myOrg[who].pnat[j];
       			myOrgDB[i].C[j]=myOrg[who].C[j];
       			myOrgDB[i].F[j]=myOrg[who].F[j];
				for (jj=j; jj<curr_MAXSTATES; jj++) myOrgDB[i].nF[j][jj]=myOrg[who].nF[j][jj];
		        myOrgDB[i].hydro[j]=myOrg[who].hydro[j];
			
            }
            
            for(j=0;j<myOrg[who].genecount-1  ;j++) {
				myOrgDB[i].pint[j]=myOrg[who].pint[j];
				myOrgDB[i].bmode[j]=myOrg[who].bmode[j];
				myOrgDB[i].Gij[j] = myOrg[who].Gij[j];
            }
			myOrgDB[i].birthrate=myOrg[who].birthrate;	
            
            myOrgDB[i].reporg = who;
            myOrgDB[i].count = 1;
            //      myOrgDB[i].meanmutrate = myOrg[who].mutrate; //mod myOrgDBMut
            if(myOrgDB[i].count > domi_species) domi_species = myOrgDB[i].count;
            nOrgDB++;
        } else {
            myOrgDB[i].count++;
            //      myOrgDB[i].meanmutrate += myOrg[who].mutrate; //mod myOrgDBMut
            if(myOrgDB[i].count > domi_species) domi_species = myOrgDB[i].count;
        }
    }
    return 0;
}


void SetMeanOrg(organism *mean){
    int ii, jj;
    //mean->genecount = mean->ppicount = 0;
    mean->dob = mean->generation = mean->numkids = 0;
    mean->G0 = mean->birthrate = mean->meanpnat = mean->minpnat = mean->mutrate = mean->mutcount = 0.0e0;
    for(ii=0;ii<curr_MAXGENES;ii++) {
        mean->si2[ii] = mean->nsi2[ii] = mean->si[ii] = mean->nsi[ii] = mean->C[ii] = mean->pnat[ii] = mean-> hydro[ii] = mean->tot_mut[ii] = mean->nonsynmut[ii] = mean->synmut[ii] = 0.0e0;
        
        mean->netcharge[ii] = 0.0e0;
        mean->poscharge[ii] = 0.0e0;
        mean->negcharge[ii] = 0.0e0;
        mean->fraccharge[ii] = 0.0e0;
        mean->frachydro[ii] = 0.0e0;
        mean->frachydro_avil[ii] = 0.0e0;
        
        mean->Nch[ii] = 0.0e0;
        mean->Kc[ii] = 0.0e0;
        
    }
    
    for(ii=0;ii<curr_MAXSTATES;ii++) {
        mean->F[ii] =  0.0e0;
        for(jj=ii;jj<curr_MAXSTATES;jj++){
			mean->K[ii][jj] = 0.0e0;
			mean->nF[ii][jj] = 0.0e0;
		}
    }
    
    if (allow_chaps==1){
        mean->F[2*curr_MAXGENES] =  0.0e0;
        mean->C[curr_MAXGENES] =  0.0e0;
    }
    
    
    for(ii=0;ii<curr_MAXPPIS;ii++) {
        mean->Gij[ii] = mean->pint[ii] = 0.0e0;
        mean->netcharge_s_partner[ii] = 0.0e0;
        mean->poscharge_s_partner[ii] = 0.0e0;
        mean->negcharge_s_partner[ii] = 0.0e0;
        mean->fraccharge_s_partner[ii] = 0.0e0;
        mean->frachydro_s_partner[ii] = 0.0e0;
        mean->frachydro_avil_s_partner[ii] = 0.0e0;
        
        mean->netcharge_s_hub[ii] = 0.0e0;
        mean->poscharge_s_hub[ii] = 0.0e0;
        mean->negcharge_s_hub[ii] = 0.0e0;
        mean->fraccharge_s_hub[ii] = 0.0e0;
        mean->frachydro_s_hub[ii] = 0.0e0;
        mean->frachydro_avil_s_hub[ii] = 0.0e0;
    }
    return;
}
 
void AddMeanOrg(organism *mean, int who){
    int ii, jj;
    mean->birthrate += myOrg[who].birthrate;
    mean->mutrate += myOrg[who].mutrate;
    mean->mutcount += myOrg[who].mutcount;
    for(ii=0;ii<curr_MAXGENES;ii++) {
        mean->C[ii] += myOrg[who].C[ii];
        mean->nsi[ii] += myOrg[who].nsi[ii];
        mean->si[ii] += myOrg[who].si[ii];
        mean->nsi2[ii] += myOrg[who].nsi2[ii];
        mean->si2[ii] += myOrg[who].si2[ii];
        mean->pnat[ii] += myOrg[who].pnat[ii];
        
        mean->hydro[ii] += myOrg[who].hydro[ii];
        mean->netcharge[ii] += myOrg[who].netcharge[ii];
        mean->poscharge[ii] += myOrg[who].poscharge[ii];
        mean->negcharge[ii] += myOrg[who].negcharge[ii];
        mean->fraccharge[ii] += myOrg[who].fraccharge[ii];
        mean->frachydro[ii] += myOrg[who].frachydro[ii];
        mean->frachydro_avil[ii] += myOrg[who].frachydro_avil[ii];
        
        mean->nonsynmut[ii] += myOrg[who].nonsynmut[ii];
        mean->tot_mut[ii] += myOrg[who].tot_mut[ii];
        mean->synmut[ii] += myOrg[who].synmut[ii];
        
        mean->Nch[ii] += myOrg[who].Nch[ii];
        mean->Kc[ii] += myOrg[who].Kc[ii];
        
        //for(jj=0;jj<curr_MAXGENES;jj++) mean->K[ii][jj] += myOrg[who].K[ii][jj];
    }
    
    for(ii=0;ii<curr_MAXSTATES;ii++) {
        mean->F[ii] += myOrg[who].F[ii];
        for(jj=ii;jj<curr_MAXSTATES;jj++){
			mean->K[ii][jj] += myOrg[who].K[ii][jj];
			mean->nF[ii][jj] += myOrg[who].nF[ii][jj];	
		} 
   }
    
    if (allow_chaps==1){
        mean->F[2*curr_MAXGENES] += myOrg[who].F[2*curr_MAXGENES];
    }
    
    if (allow_chaps==1){
        mean->C[curr_MAXGENES] += myOrg[who].C[curr_MAXGENES];
    }
    
    for(ii=0;ii<curr_MAXGENES*(curr_MAXGENES-1)/2;ii++) mean->SeqID[ii] +=myOrg[who].SeqID[ii];
    mean->G0 += myOrg[who].G0;
    for(ii=0;ii<curr_MAXPPIS;ii++) {
        mean->Gij[ii] += myOrg[who].Gij[ii];
        mean->pint[ii] += myOrg[who].pint[ii];
        
        mean->hydro_s_partner[ii] += myOrg[who].hydro_s_partner[ii];
        mean->netcharge_s_partner[ii] += myOrg[who].netcharge_s_partner[ii];
        mean->poscharge_s_partner[ii] += myOrg[who].poscharge_s_partner[ii];
        mean->negcharge_s_partner[ii] += myOrg[who].negcharge_s_partner[ii];
        mean->fraccharge_s_partner[ii] += myOrg[who].fraccharge_s_partner[ii];
        mean->frachydro_s_partner[ii] += myOrg[who].frachydro_s_partner[ii];
        mean->frachydro_avil_s_partner[ii] += myOrg[who].frachydro_avil_s_partner[ii];

        
        mean->hydro_s_hub[ii] += myOrg[who].hydro_s_hub[ii];
        mean->netcharge_s_hub[ii] += myOrg[who].netcharge_s_hub[ii];
        mean->poscharge_s_hub[ii] += myOrg[who].poscharge_s_hub[ii];
        mean->negcharge_s_hub[ii] += myOrg[who].negcharge_s_hub[ii];
        mean->fraccharge_s_hub[ii] += myOrg[who].fraccharge_s_hub[ii];
        mean->frachydro_s_hub[ii] += myOrg[who].frachydro_s_hub[ii];
        mean->frachydro_avil_s_hub[ii] += myOrg[who].frachydro_avil_s_hub[ii];

    }
    return;
}

void GetMeanOrg(organism *mean, int orgcount){
    int ii, jj;
    mean->birthrate /= (double) orgcount;
    mean->mutrate /= (double)orgcount;
    mean->mutcount /= (double)orgcount;
    for(ii=0;ii<curr_MAXGENES;ii++) {
        mean->C[ii] /= (double)orgcount;
        mean->nsi[ii] /= (double)orgcount;
        mean->si[ii] /= (double)orgcount;
        mean->nsi2[ii] /= (double)orgcount;
        mean->si2[ii] /= (double)orgcount;
        
        mean->Kc[ii] /= (double)orgcount;
        mean->Nch[ii] /= (double)orgcount;
        
        mean->pnat[ii] /= (double)orgcount;
        mean->hydro[ii] /= (double)orgcount;
        mean->netcharge[ii] /= (double)orgcount;
        mean->poscharge[ii] /= (double)orgcount;
        mean->negcharge[ii] /= (double)orgcount;
        mean->fraccharge[ii] /= (double)orgcount;
        mean->frachydro[ii] /= (double)orgcount;
        mean->frachydro_avil[ii] /= (double)orgcount;
        
        
        mean->nonsynmut[ii] /= (double)orgcount;
        mean->tot_mut[ii] /= (double)orgcount;
        mean->synmut[ii] /= (double)orgcount;
        
    }
    
    for(ii=0;ii<curr_MAXSTATES;ii++) {
        mean->F[ii] /= (double)orgcount;
        for(jj=ii;jj<curr_MAXSTATES;jj++){
			mean->K[ii][jj] /= (double)orgcount;
			mean->nF[ii][jj] /= (double)orgcount;
		}
    }
    
    if (allow_chaps==1){
        mean->F[2*curr_MAXGENES] /= (double)orgcount;
        mean->C[curr_MAXGENES] /= (double)orgcount;
        
    }
    
    
    for(ii=0;ii<curr_MAXGENES*(curr_MAXGENES-1)/2;ii++) mean->SeqID[ii] /= (double)orgcount;
    mean->G0 /= (double)orgcount;
    for(ii=0;ii<curr_MAXPPIS;ii++) {
        mean->Gij[ii] /= (double)orgcount;
        mean->pint[ii] /= (double)orgcount;
        
        mean->hydro_s_partner[ii] /= (double)orgcount;
        mean->netcharge_s_partner[ii] /= (double)orgcount;
        mean->poscharge_s_partner[ii] /= (double)orgcount;
        mean->negcharge_s_partner[ii] /= (double)orgcount;
        mean->fraccharge_s_partner[ii] /= (double)orgcount;
        mean->frachydro_s_partner[ii] /= (double)orgcount;
        mean->frachydro_avil_s_partner[ii] /= (double)orgcount;

        mean->hydro_s_hub[ii] /= (double)orgcount;
        mean->netcharge_s_hub[ii] /= (double)orgcount;
        mean->poscharge_s_hub[ii] /= (double)orgcount;
        mean->negcharge_s_hub[ii] /= (double)orgcount;
        mean->fraccharge_s_hub[ii] /= (double)orgcount;
        mean->frachydro_s_hub[ii] /= (double)orgcount;
        mean->frachydro_avil_s_hub[ii] /= (double)orgcount;

    }
    return;
}

void SetVarOrg(organism *var){	
    int ii, jj;

    for(ii=0;ii<curr_MAXGENES;ii++) {
        var->C[ii] = var->pnat[ii] = 0.0e0;
	}
    
    for(ii=0;ii<curr_MAXSTATES;ii++) {
        var->F[ii] =  0.0e0;
        for(jj=ii;jj<curr_MAXSTATES;jj++){
			var->K[ii][jj] = 0.0e0;
			var->nF[ii][jj] = 0.0e0;
		}
	}	
	for(ii=0;ii<curr_MAXPPIS;ii++) {
        var->Gij[ii] = var->pint[ii] = 0.0e0;
	}
}


void AddVarOrg(organism *var, int who){	
    int ii, jj;

    for(ii=0;ii<curr_MAXGENES;ii++) {
        var->C[ii] += pow(myOrg[who].C[ii]-mean.C[ii], 2.0);
		var->pnat[ii] += pow(myOrg[who].pnat[ii]-mean.pnat[ii], 2.0);
	}
    
    for(ii=0;ii<curr_MAXSTATES;ii++) {
        var->F[ii] +=  pow(myOrg[who].F[ii]-mean.F[ii], 2.0);
        for(jj=ii;jj<curr_MAXSTATES;jj++){
			var->K[ii][jj] += pow(myOrg[who].K[ii][jj]-mean.K[ii][jj], 2.0);
			var->nF[ii][jj] += pow(myOrg[who].nF[ii][jj]-mean.nF[ii][jj], 2.0);
		}
    }
	for(ii=0;ii<curr_MAXPPIS;ii++) {
        var->Gij[ii] += pow(myOrg[who].Gij[ii]-mean.Gij[ii], 2.0); 
		var->pint[ii] += pow(myOrg[who].pint[ii]-mean.pint[ii], 2.0);
	}

}

void GetVarOrg(organism *var, int orgcount){	
    int ii, jj;

    for(ii=0;ii<curr_MAXGENES;ii++) {
        var->C[ii] /= (double)orgcount; 
		var->pnat[ii] /= (double)orgcount;
	}
    
    for(ii=0;ii<curr_MAXSTATES;ii++) {
        var->F[ii] /=  (double)orgcount;
        for(jj=ii;jj<curr_MAXSTATES;jj++){
			var->K[ii][jj] /= (double)orgcount;
			var->nF[ii][jj] /= (double)orgcount;
		}
    }
	for(ii=0;ii<curr_MAXPPIS;ii++) {
        var->Gij[ii] /= (double)orgcount;
		var->pint[ii] /= (double)orgcount;
	}
}
