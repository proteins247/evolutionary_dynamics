/*************************
 Modified by Orit Peleg
 opeleg@fas.harvard.edu
 on 2/26/13
 comment change
 ************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "cell.h"
#include "../LP/gencode.h"
#include "../PPI/bindinglib.h"

int info;
int aaseq[AASEQLEN], aaseq2[AASEQLEN], aaseq3[AASEQLEN];

int HALF_RAND_MAX=RAND_MAX/2;
int cycle, ncycle;

double r;
double x[MAXSTATES+1], xold[MAXSTATES+1], fvec[MAXSTATES+1], Kfold[MAXGENES+1], KCh[MAXGENES+1];;
int Ustructid[MAXGENES+1][NUFOLD];

char buf[100], in_aaseq[100], temp[100];

extern double Hydrophobicity[ALPHABET];
extern double HydrophobicityYesNo[ALPHABET];
extern double HydrophobicityYesNo_avil[ALPHABET];
extern double Charge[ALPHABET];
extern int curr_MAXGENES;
extern int curr_MAXPPIS;
extern int allow_fold_change;

extern int allow_gene_exp;
extern int allow_unfolded_states;
extern int allow_chaps;
extern int hub_ID;
extern int curr_MAXSTATES;
extern float chaps_x;
extern int singlish;
extern int sequenceversion;
extern int selection;
extern int homo;

extern double TIME;

extern int RUN_BASED_ON_CONFIG;
extern int POST_Proc;


//extern FILE *it_solver;

int eff_hub_ID;
char seqbuf[800];


/***********************************************************************************************************
 **********************************************************************************************************/
void SetupParameter(int argc, char *argv[], parameter *myParam, int *orgcount){
    int who, ii, jj;
    int status; //, in_tmp2;
	//int in_tmp1;
  //  float tmp1; //, tmp2, tmp3;
    FILE *fp1;
    int bmode_temp;
    double temp_pint_max;
    int temp_bmode_max;
    
    int iii;
    int curr_initial_pop_size;
    //int in_tmp3, in_tmp4, in_tmp5;
    
    int structid[MAXGENES];
    int nucseq[MAXGENES][NUCSEQLEN];
    int bmode[MAXPPIS+1]; //, bmode_max[MAXPPIS] , surf11, surf12, surf21, surf22, surf31, surf32;
    
    double pnat[MAXGENES], hydro[MAXGENES], meanpnat, minpnat, pint[MAXPPIS]; //, pmax, pint_max[MAXPPIS];
    double netcharge[MAXGENES], poscharge[MAXGENES],
    negcharge[MAXGENES], fraccharge[MAXGENES], frachydro[MAXGENES], frachydro_avil[MAXGENES];
    
    int nucseq_hub_temp[NUCSEQLEN];
    int ii_op;
  //  char fopbuf[100];
    
    double temp_pint_sum;
    int hub_face_i;
    double min_pint_sum;
    int best_face_for_pint;
    
    
    char rootdir[200], file[200];
    
    if (argc != 5) {
        fprintf(stderr, "Usage ./EVO seed sequenceversion maxdivcycle selection\n");
        exit(-1);
    }
   
    sscanf(argv[1], "%d", &(myParam->seed));
    sscanf(argv[2], "%d", &(sequenceversion));
    sscanf(argv[3], "%d", &(myParam->maxdivcycle));
    sscanf(argv[4], "%d", &(selection));
	sprintf(myParam->targetname, "seqv%d.s%d.v%d", sequenceversion, selection, myParam->seed); 
	 
	homo=0;//hard coded, when =1 homodimer on
	allow_fold_change = 1;
	allow_gene_exp = 1;
	allow_unfolded_states = 0;
	allow_chaps = 0;
	chaps_x = 0.1;
	singlish = 1;
	curr_MAXGENES = 2;
	hub_ID = 0;	
 
    if (hub_ID==10){
        eff_hub_ID=0;
    }
    
    curr_MAXPPIS = curr_MAXGENES -1;
    
    printf("myParam.targetname=%s\n",myParam->targetname);
    printf("curr_MAXGENES=%d\n", curr_MAXGENES);
    printf("hub_ID=%d\n", hub_ID);
    printf("curr_MAXPPIS=%d\n", curr_MAXPPIS);
    printf("allow_fold_change=%d\n", allow_fold_change);
    printf("allow_gene_exp=%d\n", allow_gene_exp);
    printf("allow_unfolded_states=%d\n", allow_unfolded_states);
    printf("allow_chaps=%d\n", allow_chaps);
    printf("chaps_x=%f\n", chaps_x);
    printf("singlish=%d\n", singlish);
    printf("selection=%d\n", selection);
    printf("sequenceversion=%d\n", sequenceversion);
    printf("\n");
    
    
    if (myParam->maxdivcycle == -1) {
        POST_Proc = 1;
        printf("myParam->maxdivcycle == -1, POST_Proc = %d;\n",POST_Proc );
        
    }
    else {POST_Proc = 0;}
    
    
    sprintf(rootdir,"/n/regal/shakhnovich_lab/rrazban/%s", myParam->targetname);
    
//    sprintf(fopbuf,"%s/it_solver.dat", rootdir);
  //  it_solver=fopen(fopbuf,"w");
    
    
    if (allow_chaps==1){curr_MAXSTATES = curr_MAXGENES*2;}
    else if (allow_unfolded_states==1){curr_MAXSTATES = curr_MAXGENES*2;}
    else{curr_MAXSTATES = curr_MAXGENES;}
    
    //myOrg[who].C[MAXGENES] = 0.1;  // chaperone
    
    
    /* Parameters for Organisms */
    if (RUN_BASED_ON_CONFIG==0){
        myParam->expressrate = 0.01;
        myParam->alpha = 500.0;
        
        /* Parameters for Simulation */
        if (myParam->seed < 0) {
            myParam->seed = -myParam->seed;
        } else {
            myParam->seed *= 1024;
            myParam->seed += time(NULL);
        }
        srand(myParam->seed);
        myParam->decimthresh = 5000;
        myParam->decimto = 5000;
        myParam->dumpcycle = 100; //myParam->maxdivcycle/5; //mod myOrgDBMut
        myParam->printoutcycle = 50;
        myParam->seqlogcycle = 50; //myParam->maxdivcycle / 2000;
        
        if (myParam->seqlogcycle == 0) myParam->seqlogcycle = 1;
        //    myParam->plotoutcycle = myParam->maxdivcycle / 2000;
        //    if (myParam->plotoutcycle == 0) myParam->plotoutcycle = 1;
        myParam->screenoutcycle = myParam->maxdivcycle / 200;
        if (myParam->screenoutcycle == 0) myParam->screenoutcycle = 1;
        //myParam->maxdivcycle = 4000;
        myParam->speciessizefactor = 0.8;
        myParam->Tenv = 1.0; 
        myParam->pnatcutoff = 0.0;
        //myParam->tol = sqrt(dpmpar(1));
        myParam->mutrate[0] = 0.0001;
        myParam->mutrate[1] = 0.05;
        myParam->mutthresh = 0.01;
//        myParam->deathrate = 0.005e0;
        myParam->deathrate = 0.01; // 1/100
        myParam->fixed_mutrate = 0.001; //0.0001;
    }
    else{
        srand(myParam->seed);
    }
    
  //  sprintf(rootdir,"%s",".");
    sprintf(rootdir,"%s","/n/home12/rrazban/code/commondata");
    sprintf(file, "%s/LPforms/contact10000.dat", rootdir);
    ReadContactMatrix(file);
    sprintf(file, "%s/LPforms/allfaces10000.dat", rootdir);
    ReadAllSurfaces(file);
    sprintf(file, "%s/MJ96/energy.dat", rootdir);
    ReadEnergyMatrix(file);
   /* 
    sprintf(file, "%s/MJ96/hydrophobicity.dat", rootdir);
    fp1 = fopen(file,"r");
    while (!feof(fp1)){
        fscanf(fp1,"%d%f\n",&in_tmp1, &tmp1);
        //fprintf(stdout,"%d\t%f\n", in_tmp1, tmp1);
        Hydrophobicity[in_tmp1] = (double) tmp1;
    }
	fclose(fp1);

    sprintf(file, "%s/MJ96/hydrophobicityYesNo.dat", rootdir);
    fp1 = fopen(file,"r");
    while (!feof(fp1)){
        fscanf(fp1,"%d%f\n",&in_tmp1, &tmp1);
        //fprintf(stdout,"%d\t%f\n", in_tmp1, tmp1);
        HydrophobicityYesNo[in_tmp1] = (double) tmp1;
    }
	fclose(fp1);
    
    sprintf(file, "%s/MJ96/hydrophobicityYesNo_avil.dat", rootdir);
    fp1 = fopen(file,"r");
    while (!feof(fp1)){
        fscanf(fp1,"%d%f\n",&in_tmp1, &tmp1);
        //fprintf(stdout,"%d\t%f\n", in_tmp1, tmp1);
        HydrophobicityYesNo_avil[in_tmp1] = (double) tmp1;
    }
	fclose(fp1);
    
    sprintf(file, "%s/MJ96/charge.dat", rootdir);
    fp1 = fopen(file,"r");
    while (!feof(fp1)){
        fscanf(fp1,"%d%f\n",&in_tmp1, &tmp1);
        //fprintf(stdout,"%d\t%f\n", in_tmp1, tmp1);
        Charge[in_tmp1] = (double) tmp1;
    }
	fclose(fp1);
    */
    
    /* read initial gene sequence, Need to correct*/
    if (RUN_BASED_ON_CONFIG==0){
        fprintf(stderr, "Reading initial sequence...\n");
        sprintf(file, "%s/seqs/v%d.seq", rootdir,sequenceversion);
        if ((fp1=fopen(file,"r")) == NULL) { fprintf(stderr, "Can't open %s file\n", file); exit(-1);}
        ii = 0;
        
        while ((status = fscanf(fp1,"%s\n", buf)) != EOF) {
            if (ii<curr_MAXGENES){
                
                // op: sept 2015, read aa seq instead of nuc seq ... ?
                
                
                
                LetterToNucCodeSeq(buf, nucseq[ii], NUCSEQLEN);
                NucSeqToAASeq(nucseq[ii], NUCSEQLEN, aaseq);
                pnat[ii] = GetSequencePnat(aaseq, myParam->Tenv, structid+ii);
                fprintf(stderr,"SeqID: %6d, StrID: %6d, pnat: %6.3lf", ii, structid[ii], pnat[ii]);
                fprintf(stderr, " : %s\n", buf);
                ii++;
            }
        }
        fclose(fp1);
        
        //    //swap to make sure that the hub is in ii=0:
        //    //hub_ID
        for (ii_op=0; ii_op<NUCSEQLEN; ii_op++){nucseq_hub_temp[ii_op] = nucseq[0][ii_op];}
        for (ii_op=0; ii_op<NUCSEQLEN; ii_op++){nucseq[0][ii_op] = nucseq[eff_hub_ID][ii_op];}
        for (ii_op=0; ii_op<NUCSEQLEN; ii_op++){nucseq[eff_hub_ID][ii_op] = nucseq_hub_temp[ii_op];}
        
        //LetterToNucCodeSeq(buf, nucseq[0], NUCSEQLEN);
        NucSeqToAASeq(nucseq[0], NUCSEQLEN, aaseq);
        pnat[0] = GetSequencePnat(aaseq, myParam->Tenv, structid+0);
        
        //LetterToNucCodeSeq(buf, nucseq[hub_ID], NUCSEQLEN);
        NucSeqToAASeq(nucseq[eff_hub_ID], NUCSEQLEN, aaseq);
        pnat[eff_hub_ID] = GetSequencePnat(aaseq, myParam->Tenv, structid+eff_hub_ID);
        minpnat = 8.;
        meanpnat = 0.;
        for (ii=0; ii<curr_MAXGENES; ii++) {
            NucSeqToAASeq(nucseq[ii], NUCSEQLEN, aaseq);
            hydro[ii] = GetHydrophobicity(aaseq, AASEQLEN);
            
            frachydro[ii] = GetFracHydrophobicity(aaseq, AASEQLEN);
            frachydro_avil[ii] = GetFracHydrophobicity_avil(aaseq, AASEQLEN);
            netcharge[ii] = GetNetCharge(aaseq, AASEQLEN);
            poscharge[ii] = GetPosCharge(aaseq, AASEQLEN);
            negcharge[ii] =  GetNegCharge(aaseq, AASEQLEN);
            fraccharge[ii] = GetFracCharge(aaseq, AASEQLEN);
            
            if (pnat[ii] < minpnat) minpnat = pnat[ii];
            meanpnat += pnat[ii];
            //fprintf(stderr, "SeqID: %6d, StrID: %6d, pnat: %6.3lf, hydrophobicity: %6.3lf", ii, structid[ii], pnat[ii], hydro[ii]);
            PrintNucCodeSequence(buf, nucseq[ii], NUCSEQLEN);
            //fprintf(stderr," : %s\n", buf);
            fflush(stdout);
        }
        meanpnat /= (double)curr_MAXGENES;
        
    }
    
    
    
    
    //printf("before for(ii=0; ii < MAXPPIS; ii++) loop\n"); fflush(stdout);
    
    /* Define functional interactions, define pint and bmode: */
    if (RUN_BASED_ON_CONFIG==0){
        
        if (singlish==1) {
            //bmode[ii]= 30; //doesnt work so great..
            
            if (BEST_FACE < 0) {
                
                temp_pint_sum = 0.0;
                min_pint_sum = -1.0;
                best_face_for_pint = 0;
                
                for(hub_face_i=0; hub_face_i < 6; hub_face_i++){
                    
                    temp_pint_sum = 0.0;
                    
                    for(ii=0; ii < curr_MAXPPIS; ii++){
                        temp_pint_max = -1.0;
                        for (iii=0; iii < 24; iii++){
                            bmode_temp = iii;
                            bmode[ii]= bmode_temp + (int)( hub_face_i*24 );
                            jj=ii+1;
                            NucSeqToAASeq(nucseq[0],NUCSEQLEN,aaseq);
                            NucSeqToAASeq(nucseq[jj],NUCSEQLEN,aaseq2);
                            //use a pre-defined bmode and return pnat:
                            pint[ii] = GetBindingP2(aaseq, structid[0], aaseq2, structid[jj], bmode[ii] , myParam->Tenv);
                            if (pint[ii] > temp_pint_max){
                                temp_pint_max = pint[ii];
                                temp_bmode_max = bmode[ii];
                            }
                        }
                        bmode[ii]= temp_bmode_max;
                        jj=ii+1;
                        NucSeqToAASeq(nucseq[0],NUCSEQLEN,aaseq);
                        NucSeqToAASeq(nucseq[jj],NUCSEQLEN,aaseq2);
                        //use a pre-defined bmode and return pnat:
                        pint[ii] = GetBindingP2(aaseq, structid[0], aaseq2, structid[jj], bmode[ii] , myParam->Tenv);
                        fprintf(stderr, "PPI (%1d), ppi: %.10lE, bmode: %d\n", ii, pint[ii], bmode[ii]); fflush(stdout);
                        
                        temp_pint_sum+=pint[ii];
                    }
                    if (temp_pint_sum>min_pint_sum){
                        min_pint_sum = temp_pint_sum;
                        best_face_for_pint = hub_face_i;
                    }
                    
                }
            }
            
            else{
                best_face_for_pint = BEST_FACE;
                printf("best face (pre-defined): %d\n", BEST_FACE);
            }
            
            //update for the best face:
            for(ii=0; ii < curr_MAXPPIS; ii++){
                //printf("in for(ii=0; ii < curr_MAXPPIS; ii++) loop, ii=%d\n",ii); fflush(stdout);
                
                temp_pint_max = -1.0;
                for (iii=0; iii < 24; iii++){
                    // choose a random face and rotation.
                    //do{ bmode_temp = (int)( ( (double) rand()/RAND_MAX ) * 24 ); } while(  bmode_temp == 24 );
                    bmode_temp = iii;
                    bmode[ii]= bmode_temp + (int)( best_face_for_pint*24 );
                    jj=ii+1;
                    NucSeqToAASeq(nucseq[0],NUCSEQLEN,aaseq);
                    NucSeqToAASeq(nucseq[jj],NUCSEQLEN,aaseq2);
                    //use a pre-defined bmode and return pnat:
                    pint[ii] = GetBindingP2(aaseq, structid[0], aaseq2, structid[jj], bmode[ii] , myParam->Tenv);
                    if (pint[ii] > temp_pint_max){
                        temp_pint_max = pint[ii];
                        temp_bmode_max = bmode[ii];
                    }
                }
                if (RANDOME_MODE==0){
                    bmode[ii]= temp_bmode_max;
           			bmode[ii]=0;			//force it
			     }
                else{
                    do{bmode_temp = (int)( ( (double) rand()/RAND_MAX ) * 24 ); } while(  bmode_temp == 24 );
                    bmode[ii]= bmode_temp + (int)( best_face_for_pint*24 );
                }
                jj=ii+1;
                NucSeqToAASeq(nucseq[0],NUCSEQLEN,aaseq);
                NucSeqToAASeq(nucseq[jj],NUCSEQLEN,aaseq2);
                //use a pre-defined bmode and return pnat:
                pint[ii] = GetBindingP2(aaseq, structid[0], aaseq2, structid[jj], bmode[ii] , myParam->Tenv);
                
                //        jj=ii+1;
                //        pint[ii] = (float) GetBindingP(aaseq, structid[0], aaseq, structid[jj], &(bmode_max[ii]), myParam->Tenv);
                //        bmode[ii]=bmode_max[ii];
                
                fprintf(stderr, "PPI (%1d), ppi: %.10lE, bmode: %d\n", ii, pint[ii], bmode[ii]); fflush(stdout);
            }
            
            
        }
        
        else{ //singlish==0
            
            for(ii=0; ii < curr_MAXPPIS; ii++){
                //printf("in for(ii=0; ii < curr_MAXPPIS; ii++) loop, ii=%d\n",ii); fflush(stdout);
                
                temp_pint_max = -1.0;
                for (iii=0; iii < 24; iii++){
                    bmode_temp = iii;
                    bmode[ii]= bmode_temp + (int)( ii*24 );
                    jj=ii+1;
                    NucSeqToAASeq(nucseq[0],NUCSEQLEN,aaseq);
                    NucSeqToAASeq(nucseq[jj],NUCSEQLEN,aaseq2);
                    //use a pre-defined bmode and return pnat:
                    pint[ii] = GetBindingP2(aaseq, structid[0], aaseq2, structid[jj], bmode[ii] , myParam->Tenv);
                    if (pint[ii] > temp_pint_max){
                        temp_pint_max = pint[ii];
                        temp_bmode_max = bmode[ii];
                    }
                }
                
                if (RANDOME_MODE==0){
                    bmode[ii]= temp_bmode_max;
                }
                else{
                    do{ bmode_temp = (int)( ( (double) rand()/RAND_MAX ) * 24 ); } while(  bmode_temp == 24 );
                    bmode[ii]= bmode_temp + (int)( ii*24 );
                }
                    
                jj=ii+1;
                NucSeqToAASeq(nucseq[0],NUCSEQLEN,aaseq);
                NucSeqToAASeq(nucseq[jj],NUCSEQLEN,aaseq2);
                //use a pre-defined bmode and return pnat:
                pint[ii] = GetBindingP2(aaseq, structid[0], aaseq2, structid[jj], bmode[ii] , myParam->Tenv);
                
                //        jj=ii+1;
                //        pint[ii] = (float) GetBindingP(aaseq, structid[0], aaseq, structid[jj], &(bmode_max[ii]), myParam->Tenv);
                //        bmode[ii]=bmode_max[ii];
                
                fprintf(stderr, "PPI (%1d), ppi: %.10lE, bmode: %d\n", ii, pint[ii], bmode[ii]); fflush(stdout);
 			}
            
        }
        
    }
    
    
    
    switch (ALGORITHM) {
        case 0:
            curr_initial_pop_size = INITIAL_POPSIZE;
            break;
        case 1:
            curr_initial_pop_size = POPSIZE;
            break;
        default:
            break;
    }
    
    
    if (RUN_BASED_ON_CONFIG==0){
        TIME=0.0;
        *orgcount = 0;
        for (who=0; who<MAXORGANISMS; who++) KillOrganism(who);
        
        for (who=0; who<curr_initial_pop_size; who++) {
            //printf("who %d\n",who);
            //fprintf(stderr, "Process Org[%d]...\n", who);
            myOrgstatus[who] = S_ALIVE;
            for (ii=0;ii<curr_MAXGENES;ii++) {
                CopyIntToCharSeq(myOrg[who].genome+ii*NUCSEQLEN, nucseq[ii], NUCSEQLEN);
                myOrg[who].pnat[ii] = pnat[ii];
                myOrg[who].hydro[ii] = hydro[ii];
                myOrg[who].frachydro[ii] = frachydro[ii];
                myOrg[who].frachydro_avil[ii] = frachydro_avil[ii];
                myOrg[who].fraccharge[ii] = fraccharge[ii];
                myOrg[who].netcharge[ii] = netcharge[ii];
                myOrg[who].negcharge[ii] = negcharge[ii];
                myOrg[who].poscharge[ii] = poscharge[ii];
                
                
                //fprintf(stderr,"structid[%d]=%d\n", ii, structid[ii]);
            
                
                myOrg[who].structid[ii] = structid[ii];
            }
            myOrg[who].minpnat = minpnat;
            myOrg[who].meanpnat = meanpnat;
            myOrg[who].genecount = curr_MAXGENES;
            myOrg[who].dob = 0;
            myOrg[who].numkids = 0;
            myOrg[who].generation = 0;
            myOrg[who].mutrate = myParam->fixed_mutrate; // myParam->mutrate[0];
            myOrg[who].mutcount = 0.0;
            for (ii=0; ii<curr_MAXGENES; ii++) myOrg[who].synmut[ii] = myOrg[who].nonsynmut[ii] = myOrg[who].tot_mut[ii] = 0;
            for (ii=0; ii<curr_MAXGENES; ii++) myOrg[who].C[ii] = 0.1;
            if (allow_chaps==1){myOrg[who].C[curr_MAXGENES] = 0.1;}  // chaperone
            myOrg[who].ppicount = curr_MAXPPIS;
            myOrg[who].mutlevel = 0;
            myOrg[who].mutorigin = 0;
            
            
//            printf("ok 1 who %d\n",who);
            
            
            for (ii=0; ii<curr_MAXPPIS; ii++) {
                myOrg[who].ppi_pair[ii][0]=0, myOrg[who].ppi_pair[ii][1]=ii+1;
                myOrg[who].bmode[ii]=bmode[ii];
            }
            
            
            
            // calc birthrate
            if(who==0){
                myParam->birthrate = 1.0e0;
                myParam->b0 = 1.0;
  //              printf("ok 2 who %d\n",who);
                UpdateEquilibriumConstant(myParam, who, 0);
    //            printf("ok 3 who %d\n",who);
//				myParam->b0 = 1.0 * (myParam->deathrate / myOrg[who].birthrate);
            }
            
            UpdateEquilibriumConstant(myParam, who, 0);
            GetSequenceID(who);
            (*orgcount)++;
            UpdatePPISeqProps(who);
        }
        
    }
    
    
    fprintf(stderr, "F_i");
    for (ii=0; ii<curr_MAXSTATES; ii++) fprintf(stderr, ", %lf", myOrg[0].F[ii]);
    fprintf(stderr, "\n");
    return;
}

/***********************************************************************************************************
 **********************************************************************************************************/
// solving the law of mass action, updates x to hold the current C.
int IterativeSolver(int who){
    int i,j,k;
    double tmp, eps;
    double eps_array[MAXSTATES+1];
    int eps_flag;
    int max_iterations = 300000000;
    
    double max_error  = 10e-8;//15;
    double max_error2 = 10e-3;
    
    
    double sum,new_U,new_Ch;
    
    if (allow_unfolded_states==0 || allow_unfolded_states==2){		//RMR sims
        for(i=0;i<curr_MAXGENES;i++){
            x[i]=myOrg[who].C[i];
        }
        k=0;
        
        do {
            for(i=0;i<curr_MAXGENES;i++) xold[i]=x[i];
            
            for(i=0;i<curr_MAXGENES;i++){
                tmp=0.0;
                for(j=0;j<curr_MAXGENES;j++) tmp+=xold[j]*myOrg[who].K[i][j];
                x[i]=myOrg[who].C[i]/(1.0+tmp);
                eps_array[i] = (xold[i]-x[i])*(xold[i]-x[i])/x[i];
            }
            k++;
            eps_flag = 0;
            for(i=0;i<curr_MAXGENES;i++){
                if(eps_array[i] > max_error) eps_flag=1;
            }
            
        } while(k<=max_iterations && eps_flag == 1);
        
        for(i=0;i<curr_MAXGENES;i++){
            tmp=0.0e0;
            for(j=0;j<myOrg[who].genecount;j++) tmp+=myOrg[who].K[i][j]*x[j];
            eps_array[i] = x[i]*(1.0+tmp)-myOrg[who].C[i];
        }
        eps_flag = 0;
        for(i=0;i<curr_MAXGENES;i++){
            if(eps_array[i] > max_error2) eps_flag=1;		//what's the point of mas_error2? RMR
        }
    }
    
    else { //allow_unfolded_states==1
        for(i=0;i<curr_MAXGENES;i++){
            x[i]=myOrg[who].C[i]*myOrg[who].pnat[i]; // folded
            x[i+curr_MAXGENES]=myOrg[who].C[i]*(1.- myOrg[who].pnat[i]); // unfolded
            //printf("x[%d+curr_MAXGENES]=%f\n",i,x[i+curr_MAXGENES]);
            
            Kfold[i]=myOrg[who].pnat[i]/(1.-myOrg[who].pnat[i]);
            if (allow_chaps==1) {KCh[i]=chaps_x*Kfold[i];}
        }
        if (allow_chaps==1) {
            x[2*curr_MAXGENES]=myOrg[who].C[curr_MAXGENES];
            //printf("x[2*curr_MAXGENES]=%f\n", x[2*curr_MAXGENES]);
        } //  chaperone
        
        k=0; do {
            //for(i=0;i<curr_MAXGENES;i++){printf("CC 0 x[%d]=%E x[i+curr_MAXGENES]=%E x[2*curr_MAXGENES]=%E\n",i,x[i],x[i+curr_MAXGENES],x[2*curr_MAXGENES]);}
            
            
            for(i=0;i<curr_MAXSTATES;i++){xold[i]=x[i];}
            if (allow_chaps==1) {
                xold[2*curr_MAXGENES]=x[2*curr_MAXGENES];
            } //  chaperone
            
            eps=0.0e0;
            for(i=0;i<curr_MAXGENES;i++){
                
                if (allow_chaps==1) {
                    
                    //printf("B x[%d]=%E x[i+curr_MAXGENES]=%E x[2*curr_MAXGENES]=%E\n",i,x[i],x[i+curr_MAXGENES],x[2*curr_MAXGENES]);
                    //printf("B Kfold[i]=%E KCh[i]=%E myOrg[who].K[i+curr_MAXGENES][2*curr_MAXGENES]=%E\n", Kfold[i],KCh[i],myOrg[who].K[i+curr_MAXGENES][2*curr_MAXGENES]);
                    
                    x[i]=Kfold[i]*x[i+curr_MAXGENES]+KCh[i]*myOrg[who].K[i+curr_MAXGENES][2*curr_MAXGENES]*x[i+curr_MAXGENES]*x[2*curr_MAXGENES];
                    
                    //printf("A x[%d]=%E x[i+curr_MAXGENES]=%E x[2*curr_MAXGENES]=%E\n",i,x[i],x[i+curr_MAXGENES],x[2*curr_MAXGENES]);
                    
                    sum = x[i] + x[i+curr_MAXGENES]
                    - myOrg[who].K[i][i+curr_MAXGENES]*x[i]*x[i+curr_MAXGENES]
                    + myOrg[who].K[i+curr_MAXGENES][2*curr_MAXGENES]*
                    x[i+curr_MAXGENES]*x[2*curr_MAXGENES]; // Fi + Ui - FiUi + UiCh
                    
                    for(j=0;j<curr_MAXGENES;j++){
                        sum += myOrg[who].K[i][j]*x[i]*x[j] // FiFi
                        + myOrg[who].K[i][j+curr_MAXGENES]*x[i]*x[j+curr_MAXGENES] // FiUi
                        + myOrg[who].K[i+curr_MAXGENES][j]*x[i+curr_MAXGENES]*x[j] // UiFi
                        + myOrg[who].K[i+curr_MAXGENES][j+curr_MAXGENES]*
                        x[i+curr_MAXGENES]*x[j+curr_MAXGENES]; //UiUi
                    }
                    new_U=(myOrg[who].C[i]*x[i+curr_MAXGENES])/sum;
                    
                    //printf("new_U=%E myOrg[who].C[i]=%E x[i+curr_MAXGENES]=%E sum=%E\n",new_U,myOrg[who].C[i],x[i+curr_MAXGENES],sum);
                    x[i+curr_MAXGENES]=new_U;
                    //printf("x[%d+curr_MAXGENES]=%f\n",i,x[i+curr_MAXGENES]);
                    
                }
                else{
                    x[i]=Kfold[i]*x[i+curr_MAXGENES];
                    sum = x[i] + x[i+curr_MAXGENES] - myOrg[who].K[i][i+curr_MAXGENES]*x[i]*x[i+curr_MAXGENES];//Fi + Ui - Fi*Ui
                    for(j=0;j<curr_MAXGENES;j++){
                        sum += myOrg[who].K[i][j]*x[i]*x[j] // FiFi
                        + myOrg[who].K[i][j+curr_MAXGENES]*x[i]*x[j+curr_MAXGENES] // FiUi
                        + myOrg[who].K[i+curr_MAXGENES][j]*x[i+curr_MAXGENES]*x[j] // UiFi
                        + myOrg[who].K[i+curr_MAXGENES][j+curr_MAXGENES]*x[i+curr_MAXGENES]*x[j+curr_MAXGENES]; //UiUi
                    }
                    new_U=(myOrg[who].C[i]*x[i+curr_MAXGENES])/sum;
                    x[i+curr_MAXGENES]=new_U;
                }
                //exit(0);
            }
            
            
            if (allow_chaps==1) {
                // with chaperone
                sum=x[2*curr_MAXGENES]; //Ch
                for(i=0;i<curr_MAXGENES;i++){
                    sum += myOrg[who].K[i+curr_MAXGENES][2*curr_MAXGENES]*x[i+curr_MAXGENES]*x[2*curr_MAXGENES]; // (UiCh)
                }
                new_Ch=(myOrg[who].C[curr_MAXGENES]*x[2*curr_MAXGENES])/sum;
                x[2*curr_MAXGENES]=new_Ch;
            }
            
            
            eps_flag = 0;
            for(i=0;i<curr_MAXSTATES;i++){
                eps_array[i] = (xold[i]-x[i])*(xold[i]-x[i])/x[i];
                if(eps_array[i] > max_error) eps_flag=1;
                //printf("%E ", eps_array[i]);
            }
            //printf("\n");
            if (allow_chaps==1) {
                eps_array[2*curr_MAXGENES] = (xold[2*curr_MAXGENES]-x[2*curr_MAXGENES])*(xold[2*curr_MAXGENES]-x[2*curr_MAXGENES])/x[2*curr_MAXGENES];
                if(eps_array[2*curr_MAXGENES] > max_error) eps_flag=1;
                //printf("CHAP xold:%f xnew:%f error: %E\n",xold[2*curr_MAXGENES], x[2*curr_MAXGENES],eps_array[2*curr_MAXGENES] );
            }
            
            k++;
            
        } while(k<=max_iterations && eps_flag == 1);
        
        //printf("ok 7 who %d\n",who);exit(0);
        
        
        for (i=0;i<curr_MAXGENES;i++){
            if (allow_chaps==1) {
                sum = x[i] + x[i+curr_MAXGENES] - myOrg[who].K[i][i+curr_MAXGENES]*x[i]*x[i+curr_MAXGENES]
                + myOrg[who].K[i+curr_MAXGENES][2*curr_MAXGENES]*x[i+curr_MAXGENES]*x[2*curr_MAXGENES]; // Fi + Ui - FiUi + UiCh
                
                for(j=0;j<curr_MAXGENES;j++){
                    sum += myOrg[who].K[i][j]*x[i]*x[j] // FiFi
                    + myOrg[who].K[i][j+curr_MAXGENES]*x[i]*x[j+curr_MAXGENES] // FiUi
                    + myOrg[who].K[i+curr_MAXGENES][j]*x[i+curr_MAXGENES]*x[j] // UiFi
                    + myOrg[who].K[i+curr_MAXGENES][j+curr_MAXGENES]*x[i+curr_MAXGENES]*x[j+curr_MAXGENES]; //UiUi
                }
                
                tmp=myOrg[who].C[i]-sum;
                eps += tmp*tmp;
            }
            else{
                sum = x[i] + x[i+curr_MAXGENES] - myOrg[who].K[i][i+curr_MAXGENES]*x[i]*x[i+curr_MAXGENES];
                for(j=0;j<curr_MAXGENES;j++){
                    sum += myOrg[who].K[i][j]*x[i]*x[j] // FiFi
                    + myOrg[who].K[i][j+curr_MAXGENES]*x[i]*x[j+curr_MAXGENES] // FiUi
                    + myOrg[who].K[i+curr_MAXGENES][j]*x[i+curr_MAXGENES]*x[j] // UiFi
                    + myOrg[who].K[i+curr_MAXGENES][j+curr_MAXGENES]*x[i+curr_MAXGENES]*x[j+curr_MAXGENES]; //UiUi
                }
                tmp=myOrg[who].C[i]-sum;
                eps_array[i] += tmp*tmp;
            }
        }
        
        //printf("ok 8 who %d\n",who);
        
        for(i=0;i<curr_MAXSTATES;i++){
            if(eps_array[i] > max_error2) eps_flag=1;
        }
        
        //printf("ok 9 who %d\n",who);
        
        if (allow_chaps==1) {
            sum=x[2*curr_MAXGENES]; //Ch
            for(i=0;i<curr_MAXGENES;i++){
                sum += myOrg[who].K[i+curr_MAXGENES][2*curr_MAXGENES]*x[i+curr_MAXGENES]*x[2*curr_MAXGENES]; // (UiCh)
            }
            tmp=myOrg[who].C[curr_MAXGENES]-sum;
            eps_array[2*curr_MAXGENES] = tmp*tmp;
            if(eps_array[2*curr_MAXGENES] > max_error2) eps_flag=1;
        }
        
        //printf("ok 10 who %d\n",who);
    }
    
    //printf("ok 11 who %d\n",who);
    
    if(k>=max_iterations || eps_flag == 1) {
        printf("k=%d eps_flag=%d\n",k,eps_flag);
        for(i=0;i<curr_MAXSTATES;i++){
            //printf("ok 11.2 i %d\n",i);
            printf("%E ", eps_array[i]);
//            fprintf(it_solver, "%E ", eps_array[i]);
        }
        printf("\n");
        //printf("ok 11.1 who %d\n",who);
        
        fprintf(stderr,"IterativeSolver : Soultion not converged in %d iteration\n", k);
        fprintf(stderr,"xold : %lf %lf %lf, x : %lf %lf %lf; eps : %E %E %E\n" ,xold[0], xold[1], xold[2], x[0], x[1], x[2], eps_array[0], eps_array[1], eps_array[2]);
        
        //printf("ok 11.2 who %d\n",who);
        
        for(i=0;i<myOrg[who].genecount;i++) {
            fprintf(stderr, "x[%1d] : %10.4lf, C[%1d] : %10.4lf",i, x[i], i,myOrg[who].C[i]);
            for(j=0;j<myOrg[who].genecount;j++) {
                fprintf(stderr, ", K[%1d][%1d] : %10.4lf", i,j, myOrg[who].K[i][j]);
            }
            fprintf(stderr,"\n");
        }
        
        //printf("ok 11.3 who %d\n",who);
        
        exit(-1);
    }
    
    //printf("ok 12 who %d\n",who);
    
    
    //for(i=0;i<curr_MAXGENES;i++){printf("CC 2 x[%d]=%E x[i+curr_MAXGENES]=%E x[2*curr_MAXGENES]=%E\n",i,x[i],x[i+curr_MAXGENES],x[2*curr_MAXGENES]);}
    //exit(0);
    
    return 0;
}

/***********************************************************************************************************
 **********************************************************************************************************/
// Generates a normally distributed random number with expectation value 0 and deviation 0.1
float GaussianNum() {
    
    static int iset=0;
    static float gset;
    
    float fac,rsq,v1,v2;
    
    if (iset==0) {
        do {
            
            v1 = 2.0*drand48()-1.0;
            v2 = 2.0*drand48()-1.0;
            
            rsq = v1*v1+v2*v2;
            
        } while (rsq >= 1.0 || rsq == 0.0);
        
        //fac = sqrt(-2.0*log(rsq)/rsq);
        fac = sqrt(-0.02*log(rsq)/rsq);
        gset = v1*fac;
        iset=1;
        
        return v2*fac;
        
    } else {
        iset=0;
        return gset;
    }
}

/***********************************************************************************************************
 **********************************************************************************************************/
int UpdateBirthrateStoch(parameter *myParam, int who, int func){
    int i,j,k;
    int ii;
    float C;
    UpdateMonomerConcentration(myParam, who);
    UpdateNsi(myParam, who);
    
    
    for(k=0;k<myOrg[who].ppicount;k++) {
        i=myOrg[who].ppi_pair[k][0], j=myOrg[who].ppi_pair[k][1];
	    myOrg[who].Gij[k]=myOrg[who].F[i]*myOrg[who].F[j]*myOrg[who].K[i][j]*myOrg[who].pint[k]*myOrg[who].pnat[i]*myOrg[who].pnat[j];
    }
    
    
    myOrg[who].birthrate = myParam->b0; //10E37;
    for (ii=0; ii<myOrg[who].ppicount; ii++) {
        if (allow_unfolded_states==2){
            i=myOrg[who].ppi_pair[ii][0], j=myOrg[who].ppi_pair[ii][1];
            myOrg[who].birthrate *=myOrg[who].Gij[ii]/(myOrg[who].pnat[i]*myOrg[who].pnat[j]);
        }
        else if (hub_ID==10){
            //take into account only first pair
            if (ii==0){
                myOrg[who].birthrate *=myOrg[who].Gij[ii];
            }
            else if (ii==1){
            }
            else{
                myOrg[who].birthrate *=myOrg[who].F[ii]*myOrg[who].pnat[ii];
            }
        }
		else if (selection==0){
			myOrg[who].birthrate *= myOrg[who].Gij[ii]/(myOrg[who].K[i][j]*myOrg[who].pint[ii]);
		}
        else{
            myOrg[who].birthrate *=myOrg[who].Gij[ii];
        }
    }

    if (hub_ID==10){
        myOrg[who].birthrate *=myOrg[who].F[ii]*myOrg[who].pnat[ii];
    }
    
    
    C=0.0e0;
    for(i=0;i<myOrg[who].genecount;i++) C+=myOrg[who].C[i];
    myOrg[who].birthrate/=(1.0e0+myParam->alpha*(C-((curr_MAXGENES+1)*0.1))*(C-((curr_MAXGENES+1)*0.1)));
    
    
    if(allow_chaps==1){
        for(i=0;i<curr_MAXGENES;i++){
            myOrg[who].Kc[i] = KCh[i]*myOrg[who].K[i+curr_MAXGENES][2*curr_MAXGENES]*myOrg[who].F[2*curr_MAXGENES];
        }
    }
    return 0;
    
}


/*
 * mutorigin 0 (non-mutator),
 *           1 (n->m via mutation), 2 (n->m via gene expression), 3 (n->m via T-jump)
 *           4 (m->n via mutation), 5 (m->n via gene expression), 6 (m->n via T-jump)
 * func : 0 (mutation), 1 (gene expression), 2 (temperature jump)
 */
/***********************************************************************************************************
 **********************************************************************************************************/
int OrgDeath(parameter *myParam, int who, int compensate){
    ncycle=1;
    r = compensate*myParam->deathrate;
    
    while(r>1.0){
        r/=2;
        ncycle*=2;
    }
    
    for(cycle=0;cycle<ncycle;cycle++){
        if(((double) rand()/RAND_MAX)<r) {
            myOrgstatus[who] = S_TODIE;
            return 1;
        }
    }
    return 0;
}

/***********************************************************************************************************
 **********************************************************************************************************/
int GeneExpress(parameter *myParam, int who, int compensate, double *dC){
    int ii;
    
    r = compensate*myParam->expressrate;
    switch (ALGORITHM) {
        case 0:
            ncycle=1;
            while(r>1.0){
                r/=2;
                ncycle*=2;
            }
            break;
        case 1:
            ncycle = 1;
            break;
            
        default:
            break;
    }
    
    
    for(cycle=0;cycle<ncycle;cycle++){
        /* Select 1 gene in a cell and change the expression level of the gene */
        if(((double) rand()/RAND_MAX) < r) {
            do { ii = (int) (((double)rand()/RAND_MAX)*myOrg[who].genecount); } while(ii == myOrg[who].genecount);
            if (homo==1) {
            	do { ii = (int) (((double)rand()/RAND_MAX)*(myOrg[who].genecount-1)); } while(ii == (myOrg[who].genecount-1));
			}
            *dC=GaussianNum();
            myOrg[who].C[ii]*=(1.0e0+*dC);
			if (homo==1) {
            	myOrg[who].C[ii+1]*=(1.0e0+*dC);
       		}
		 }
    }
    
    UpdateBirthrateStoch(myParam, who, 1);
    return 0;
}

/***********************************************************************************************************
 **********************************************************************************************************/
int OrgGeneMutate(parameter *myParam, int who, int compensate){
    int ii, jj, kk;
    int mutcount=0, nonsyn_mutcount=0;
    int StrID;
   	int genecount; 
    
    r=compensate*myOrg[who].mutrate;
    switch (ALGORITHM) {
        case 0:
            ncycle=1;
            while(r>1.0){
                r/=2;
                ncycle*=2;
            }
            break;
        case 1:
            ncycle = 1;
            break;
            
        default:
            break;
    }
    
    for(cycle=0;cycle<ncycle;cycle++){
        if(((double) rand()/RAND_MAX) < r) {
			if (homo==1) {
				genecount=myOrg[who].genecount-1;
			}
			else {
				genecount=myOrg[who].genecount;
			}

            for(ii=0;ii<genecount;ii++){
                //if(((double) rand()/RAND_MAX)<0.3) continue;
                if(myOrg[who].genome[ii*NUCSEQLEN]==-1) {fprintf(stderr,"genome corruption at dup\n"); exit(1);}
                do {
                    memcpy(temp, myOrg[who].genome+ii*NUCSEQLEN,NUCSEQLEN);
                    jj = PointMutateCharNucSequence(temp, NUCSEQLEN);
                    //printf("PointMutateCharNucSequence who=%d gene=%d\n", who,ii);
                } while(jj==-1); //stop codon
                memcpy(myOrg[who].genome+ii*NUCSEQLEN,temp,NUCSEQLEN);
                if (homo==1) {
                memcpy(myOrg[who].genome+(ii+1)*NUCSEQLEN,temp,NUCSEQLEN);
				}
                if(jj==1) { //nonsynonymous
                    //printf("PointMutateCharNucSequence who=%d gene=%d NON SYNON\n", who,ii);
                    CharNucSeqToAASeq(myOrg[who].genome+ii*NUCSEQLEN, NUCSEQLEN, aaseq);
                    
                    if (allow_fold_change==1){ //muyoung
                        // op sept 2015, allow fold change but don't change structure!
                        myOrg[who].pnat[ii] = (float) GetStructurePnat(aaseq,myParam->Tenv,myOrg[who].structid[ii]);

                    }
                    else{ //murat
                        myOrg[who].pnat[ii] = (float) GetSequencePnat(aaseq,myParam->Tenv,&StrID);
                        if(StrID != myOrg[who].structid[ii]){
                            do{ kk = (int)( ( (double) rand()/RAND_MAX ) * MAXORGANISMS );
                            } while ( kk==MAXORGANISMS || kk==who );
                            memmove(&myOrg[who], &myOrg[kk], sizeof(organism));
                            nonsyn_mutcount=0;
                            //printf("PointMutateCharNucSequence who=%d gene=%d NON SYNON, killed\n", who,ii);
                            break;
                        }
                    }
                    
//disregard that homo case fails
                    myOrg[who].hydro[ii]=(float) GetHydrophobicity(aaseq, AASEQLEN);
                    myOrg[who].frachydro[ii] = (float)GetFracHydrophobicity(aaseq, AASEQLEN);


                    myOrg[who].frachydro_avil[ii] = (float)GetFracHydrophobicity_avil(aaseq, AASEQLEN);
                    myOrg[who].fraccharge[ii] = (float)GetFracCharge(aaseq, AASEQLEN);
                    myOrg[who].netcharge[ii] = (float)GetNetCharge(aaseq, AASEQLEN);
                    myOrg[who].negcharge[ii] = (float)GetNegCharge(aaseq, AASEQLEN);
                    myOrg[who].poscharge[ii] = (float)GetPosCharge(aaseq, AASEQLEN);
                    
                    nonsyn_mutcount++;
                    myOrg[who].nonsynmut[ii]++;
                    myOrg[who].tot_mut[ii]++;
                } else {
                    myOrg[who].synmut[ii]++;
                    myOrg[who].tot_mut[ii]++;
                }
                mutcount++;
            } // for-loop myOrg[who].genecount
            UpdatePPISeqProps(who);
        }
    } // for cycle
    
    
    if(nonsyn_mutcount>0) {
        UpdateEquilibriumConstant(myParam, who, 0);
    }
    
    switch (ALGORITHM) {
        case 0:
            if(nonsyn_mutcount>0) {
                myOrg[who].meanpnat = 0.;
                myOrg[who].minpnat = 8;
                
                for(ii=0;ii<myOrg[who].genecount;ii++){
                    myOrg[who].meanpnat += myOrg[who].pnat[ii];
                    if(myOrg[who].minpnat > myOrg[who].pnat[ii]) myOrg[who].minpnat = myOrg[who].pnat[ii];
                    if(myOrg[who].pnat[ii] < myParam->pnatcutoff) {
                        if(myOrgstatus[who] == S_ALIVE) {
                            myOrgstatus[who]=S_TODIE;
                        } else if(myOrgstatus[who] == S_NEWBORN) {
                            myOrgstatus[who]=S_TODIE2;
                        }
                    }
                }
                myOrg[who].meanpnat /= (myOrg[who].genecount);
                if(myOrg[who].minpnat > 1.) {fprintf(stderr,"minpnat exceeds 1.\n"); exit(1);}
                //UpdateEquilibriumConstant(myParam, who, 0);
                GetSequenceID(who);
            }
            if(mutcount>0) {
                myOrg[who].mutcount += (((float) mutcount) / (myOrg[who].genecount));
            }
            
            if(myOrgstatus[who]==S_TODIE || myOrgstatus[who]==S_TODIE2) {
                return -1;
            } else {
                return 0;
            }
            
            break;
        case 1:
            return 0;
            break;
            
        default:
            break;
    }
}


/***********************************************************************************************************
 **********************************************************************************************************/
void UpdatePPISeqProps(int who){
    int aaseq_hub[NUCSEQLEN], aaseq_partner[NUCSEQLEN], aaseq_surface_hub[AASURFACELEN],
    aaseq_surface_partner[AASURFACELEN];
    int ii, i,j;
    
    
    for (ii=0;ii<curr_MAXPPIS; ii++){
        i=myOrg[who].ppi_pair[ii][0], j=myOrg[who].ppi_pair[ii][1];
        //printf("i=%d j=%d\n", i,j);
        CharNucSeqToAASeq(myOrg[who].genome+i*NUCSEQLEN,NUCSEQLEN,aaseq_hub);
        CharNucSeqToAASeq(myOrg[who].genome+j*NUCSEQLEN,NUCSEQLEN,aaseq_partner);
        
        GetSurfaceAA(aaseq_hub, aaseq_partner, aaseq_surface_hub, aaseq_surface_partner, myOrg[who].structid[i], myOrg[who].structid[j], myOrg[who].bmode[ii]);
        
        myOrg[who].hydro_s_partner[ii]=(float) GetHydrophobicity(aaseq_surface_partner, AASURFACELEN);
        myOrg[who].frachydro_s_partner[ii] = (float)GetFracHydrophobicity(aaseq_surface_partner, AASURFACELEN);
        myOrg[who].frachydro_avil_s_partner[ii] = (float)GetFracHydrophobicity_avil(aaseq_surface_partner, AASURFACELEN);
        myOrg[who].fraccharge_s_partner[ii] = (float)GetFracCharge(aaseq_surface_partner, AASURFACELEN);
        myOrg[who].netcharge_s_partner[ii] = (float)GetNetCharge(aaseq_surface_partner, AASURFACELEN);
        myOrg[who].negcharge_s_partner[ii] = (float)GetNegCharge(aaseq_surface_partner, AASURFACELEN);
        myOrg[who].poscharge_s_partner[ii] = (float)GetPosCharge(aaseq_surface_partner, AASURFACELEN);
        
        myOrg[who].hydro_s_hub[ii]=(float) GetHydrophobicity(aaseq_surface_hub, AASURFACELEN);
        myOrg[who].frachydro_s_hub[ii] = (float)GetFracHydrophobicity(aaseq_surface_hub, AASURFACELEN);
        myOrg[who].frachydro_avil_s_hub[ii] = (float)GetFracHydrophobicity_avil(aaseq_surface_hub, AASURFACELEN);
        myOrg[who].fraccharge_s_hub[ii] = (float)GetFracCharge(aaseq_surface_hub, AASURFACELEN);
        myOrg[who].netcharge_s_hub[ii] = (float)GetNetCharge(aaseq_surface_hub, AASURFACELEN);
        myOrg[who].negcharge_s_hub[ii] = (float)GetNegCharge(aaseq_surface_hub, AASURFACELEN);
        myOrg[who].poscharge_s_hub[ii] = (float)GetPosCharge(aaseq_surface_hub, AASURFACELEN);
    }
}


/***********************************************************************************************************
 **********************************************************************************************************/
int OrgChildBirth(parameter *myParam, int who, int divisioncycle, int compensate){
    int i, output;
    
    ncycle=1;
    r=compensate*myOrg[who].birthrate;
    
    //fprintf(stdout,"BirthRate : %lf\n", r);
    
    while(r>1.0){
        r/=2;
        ncycle*=2;
    }
    
    output=0;
    for(cycle=0;cycle<ncycle;cycle++){
        if(((double) rand()/RAND_MAX)<r) {
            for(i=0;i<MAXORGANISMS;i++) {
                if(myOrgstatus[i] == S_DEAD) break;
            }
            if(i != MAXORGANISMS) {
                myOrg[who].numkids++;
                memmove(&myOrg[i], &myOrg[who], sizeof(organism));
                myOrgstatus[i] = S_NEWBORN;
                myOrg[i].dob = divisioncycle;
                myOrg[i].numkids = 0;
                myOrg[i].generation++;
                
                /* mutations occurs only when cell division occurs */
                output++;
                output+=OrgGeneMutate(myParam,who,1);
                output+=OrgGeneMutate(myParam,i,1);
            } else {
                return -9; // no free org.
            }
        }
    }
    
    return output;
}


/***********************************************************************************************************
 **********************************************************************************************************/
void KillOrganism(int who){
    int i,j;
    myOrgstatus[who] = S_DEAD;
    
    for(i=0;i<curr_MAXGENES*NUCSEQLEN;i++) myOrg[who].genome[i]=-1;
    for(i=0;i<curr_MAXGENES;i++)
    {
        myOrg[who].structid[i] = -1;
        myOrg[who].pnat[i] = -1;
        myOrg[who].hydro[i] = -1;
        myOrg[who].C[i]=0.0;
        for(j=0;j<curr_MAXGENES;j++) myOrg[who].K[i][j]=0.0;
    }
    
    //for(i=0;i<4;i++) myOrg[who].action_count[i]=0;
    myOrg[who].minpnat = -1;
    myOrg[who].meanpnat = -1;
    myOrg[who].genecount = -1;
    myOrg[who].dob = -1;
    myOrg[who].numkids = -1;
    myOrg[who].generation = -1;
    myOrg[who].mutcount = -1;
    for(i=0;i<curr_MAXGENES;i++) myOrg[who].synmut[i] = myOrg[who].nonsynmut[i] = myOrg[who].tot_mut[i] = 0;
    return;
}

/***********************************************************************************************************
 **********************************************************************************************************/
void PrintInitialCondition(FILE *out, parameter *myParam){
    fprintf(out, "------- Organism Information -------\n");
    fprintf(out, "Target Name : \t\t%s\n", myParam->targetname);
    fprintf(out, "Organism Count : \t%d\n", myParam->orgcount);
    fprintf(out, "Temperature : \t\t%.5f\n", myParam->Tenv);
    fprintf(out, "BirthRate : \t\t%.5f\n", myParam->birthrate);
    fprintf(out, "DeathRate : \t\t%.5f\n", myParam->deathrate);
    fprintf(out, "ExpressRate : \t\t%.5f\n", myParam->expressrate);
    fprintf(out, "Alpha : \t\t%.5f\n", myParam->alpha);
    fprintf(out, "------- Simulation Information -------\n");
    fprintf(out, "random # seed : \t%d\n", myParam->seed);
    fprintf(out, "MaxDivCycle : \t\t%d\n", myParam->maxdivcycle);
    fprintf(out, "ProteomeDumpCycle : \t%d\n", myParam->dumpcycle);
    fprintf(out, "PrintOutCycle : \t%d\n", myParam->printoutcycle);
    fprintf(out, "PlotOutCycle : \t\t%d\n", myParam->plotoutcycle);
    fprintf(out, "SeqLogCycle : \t\t%d\n", myParam->seqlogcycle);
    fprintf(out, "ScreenOutCycle : \t%d\n", myParam->screenoutcycle);
    fprintf(out, "Sp Size Factor : \t%lf\n", myParam->speciessizefactor);
    fprintf(out, "Wildtype Mut Rate : \t%.5f\n", myParam->mutrate[0]);
    fprintf(out, "Mutator Mut Rate : \t%.5f\n", myParam->mutrate[1]);
    fprintf(out, "Mutation Threshold : \t%lf\n", myParam->mutthresh);
    fprintf(out, "Mut Thresh Factor : \t%6.4lf\n", myParam->mutthreshfac);
    fprintf(out, "Decim. Threshold : \t%d\n", myParam->decimthresh);
    fprintf(out, "Decim. To : \t\t%d\n", myParam->decimto);
    fprintf(out, "------- Output Files -------\n");
    return;
}


/***********************************************************************************************************
 **********************************************************************************************************/
int UpdateMonomerConcentration(parameter *myParam, int who){
    int i;
	int j;
    IterativeSolver(who);
    
    for(i=0;i<curr_MAXSTATES;i++){
        myOrg[who].F[i]=(float) x[i];
    }
   
    for(i=0;i<curr_MAXSTATES;i++){
    	for(j=i;j<curr_MAXSTATES;j++){
 			myOrg[who].nF[i][j] = myOrg[who].F[i]*myOrg[who].F[j]*myOrg[who].K[i][j];
 		}
	}
    if (allow_chaps){
        myOrg[who].F[2*curr_MAXGENES]=(float) x[2*curr_MAXGENES];
    }
    
    return 0;
}


//    // non-functional interactions
//
//    for(i=0;i<MAXGENES;i++){
//        myOrg[who].NF[i]=myOrg[who].C[i]-myOrg[who].F[i]-myOrg[who].F[i+MAXGENES];
//    }
//    for(i=0;i<MAXPPIS;i++){						RMR- I dont think this is right because G comes from folded proteins
//        myOrg[who].NF[0]-=myOrg[who].Gij[i];
//    }
//    for(i=1;i<MAXGENES;i++){
//        myOrg[who].NF[i]-=myOrg[who].Gij[i];
//    }
//

/***********************************************************************************************************
 **********************************************************************************************************/
void UpdateNsi(parameter *myParam, int who){
    int i,j;
    float curr_G, curr_temp_sum;
    
    for(i=0;i<myOrg[who].genecount;i++){
        curr_G = myOrg[who].F[i] * myOrg[who].pnat[i];
        if (i==0) {
            curr_temp_sum = 0.0;
            for(j=0;j<curr_MAXPPIS;j++){
                curr_temp_sum += myOrg[who].Gij[j];
            }
            myOrg[who].nsi[i] =
            1.0 - ((1.0/(myOrg[who].C[i]*myOrg[who].pnat[i]))*(curr_G+ curr_temp_sum));
            myOrg[who].si[i] =
            ((1.0/(myOrg[who].C[i]*myOrg[who].pnat[i]))*(curr_temp_sum));
        }
        else {
            myOrg[who].nsi[i] =
            1.0 - ((1.0/(myOrg[who].C[i]*myOrg[who].pnat[i]))*(curr_G+ myOrg[who].Gij[i-1]));
            myOrg[who].si[i] =
            ((1.0/(myOrg[who].C[i]*myOrg[who].pnat[i]))*(myOrg[who].Gij[i-1]));
        }
        
    }
    
    
    for(i=0;i<myOrg[who].genecount;i++){
        if (i==0) {
            curr_temp_sum = 0.0;
            for(j=0;j<curr_MAXPPIS;j++){
                curr_temp_sum += myOrg[who].Gij[j];
            }
            myOrg[who].nsi2[i] = myOrg[who].C[i]-myOrg[who].F[i]-myOrg[who].F[i+curr_MAXGENES] - curr_temp_sum;
        }
        else {
            myOrg[who].nsi2[i] = myOrg[who].C[i]-myOrg[who].F[i]-myOrg[who].F[i+curr_MAXGENES] - myOrg[who].Gij[i-1];
        }
        
    }
    
    for(i=0;i<curr_MAXGENES;i++){
        myOrg[who].Nch[i] = 0.0;
    }
    
    if (allow_chaps==1){
        // every unfolded protein is also in a complex with chaperone
        for(i=0;i<curr_MAXGENES;i++){
            myOrg[who].nsi2[i] -= myOrg[who].F[i+curr_MAXGENES]*myOrg[who].F[2*curr_MAXGENES]*myOrg[who].K[i+curr_MAXGENES][2*curr_MAXGENES];
            myOrg[who].Nch[i] = myOrg[who].F[i+curr_MAXGENES]*myOrg[who].F[2*curr_MAXGENES]*myOrg[who].K[i+curr_MAXGENES][2*curr_MAXGENES];
        }
    }
    
    for(i=0;i<myOrg[who].genecount;i++){
        myOrg[who].si2[i] = 1.0 - myOrg[who].nsi2[i];
    }
}


/***********************************************************************************************************
 **********************************************************************************************************/
int UpdateEquilibriumConstant(parameter *myParam, int who, int func){
    int i,j,k,kk,ii,UID,UID2;
    double sum,norm,C;
    
    if (allow_unfolded_states==1){
        // Generate random conformations for unfolded ensemble
        for(i=0;i<curr_MAXGENES;i++){
            for(j=0;j<NUFOLD;j++){
                do{
                    do{ UID = (int)( ( (double) rand()/RAND_MAX ) * 1000 ); } while( UID==1000 );
                    UID2 = UID + (int)( j*1000 );
                }while( UID2==myOrg[who].structid[i] );
                Ustructid[i][j] = UID2;
            }
        }
        for(i=0;i<curr_MAXGENES;i++){
            CharNucSeqToAASeq(myOrg[who].genome+i*NUCSEQLEN,NUCSEQLEN,aaseq);
            
            // K_ii for Ui-Ui
            sum=0.0e0;
            for(k=0;k<NUFOLD;k++){
                UID=Ustructid[i][k];
                for(kk=k;kk<NUFOLD;kk++){
                    UID2=Ustructid[i][kk];
                    sum += GetBindingK(aaseq,UID,aaseq,UID2,myParam->Tenv);
                }
            }
            norm=2.0/(NUFOLD*(NUFOLD+1.0));
            myOrg[who].K[i+curr_MAXGENES][i+curr_MAXGENES] = sum*norm;
            //printf("in updateEq, myOrg[who].K[%d+MAXGENES][%d+MAXGENES]=%f\n",i,i,myOrg[who].K[i+curr_MAXGENES][i+curr_MAXGENES] );
            
            // K_ij for Ui-Ch
            sum=0.0e0;
            for(k=0;k<NUFOLD;k++){
                UID=Ustructid[i][k];
                sum += GetBindingCh(aaseq,UID,myParam->Tenv);
            }
            norm=1.0e0/NUFOLD;
            myOrg[who].K[i+curr_MAXGENES][2*curr_MAXGENES]=myOrg[who].K[2*curr_MAXGENES][i+curr_MAXGENES]= sum*norm;
            //printf("in updateEq, myOrg[who].K[%d+MAXGENES][2*MAXGENES]=%f\n",i,myOrg[who].K[i+MAXGENES][2*MAXGENES]);
            
            //printf("ok 4 who %d\n",who);
            
            for(j=i;j<curr_MAXGENES;j++){
                CharNucSeqToAASeq(myOrg[who].genome+j*NUCSEQLEN,NUCSEQLEN,aaseq2);
                
                // K_ij for Fi-Fj
                myOrg[who].K[i][j]=myOrg[who].K[j][i] = GetBindingK(aaseq,myOrg[who].structid[i],aaseq2,myOrg[who].structid[j],myParam->Tenv);
                
                // K_ij for Fi-Uj
                sum=0.0e0;
                for(k=0;k<NUFOLD;k++){
                    UID=Ustructid[j][k];
                    sum += GetBindingK(aaseq,myOrg[who].structid[i],aaseq2,UID,myParam->Tenv);
                }
                norm=1.0e0/((double)NUFOLD);
                myOrg[who].K[i][j+curr_MAXGENES]=myOrg[who].K[j+curr_MAXGENES][i] = sum*norm;
                
                // K_ij for Ui-Uj and i not equal to j
                
                if (i == j){ continue; }
                
                sum=0.0e0;
                for(k=0;k<NUFOLD;k++){
                    UID=Ustructid[i][k];
                    for(kk=0;kk<NUFOLD;kk++){
                        UID2=Ustructid[j][kk];
                        sum += GetBindingK(aaseq,UID,aaseq2,UID2,myParam->Tenv);
                    }
                }
                norm=1.0e0/((double)(NUFOLD*NUFOLD));
                myOrg[who].K[i+curr_MAXGENES][j+curr_MAXGENES]=myOrg[who].K[j+curr_MAXGENES][i+curr_MAXGENES] = sum*norm;
                
            } // j
        } // i
        
    }
    
    else{		//RMR sims
        
        for(i=0;i<curr_MAXGENES;i++){
            CharNucSeqToAASeq(myOrg[who].genome+i*NUCSEQLEN,NUCSEQLEN,aaseq);
            for(j=i;j<curr_MAXGENES;j++){
                CharNucSeqToAASeq(myOrg[who].genome+j*NUCSEQLEN,NUCSEQLEN,aaseq2);
	
                myOrg[who].K[i][j]=myOrg[who].K[j][i]= (float) GetBindingK(aaseq,myOrg[who].structid[i],aaseq2,myOrg[who].structid[j],myParam->Tenv);
                
    //            myOrg[who].K[i][j]=myOrg[who].K[j][i]= (float) GetBindingK2(aaseq,myOrg[who].structid[i],aaseq2,myOrg[who].structid[j], myOrg[who].bmode[k], myParam->Tenv);
    //            printf("genes: i=%d, j=%d, myOrg[who].K[i][j]=%lf\n", i,j,myOrg[who].K[i][j]);
                
            }
        }
    }
    
    
    
    UpdateMonomerConcentration(myParam, who);
    
    //printf("ok 5.5 who %d\n",who);
    UpdateNsi(myParam, who);
    
    for(i=0;i<myOrg[who].genecount;i++){
        //printf("gene: i=%d, myOrg[who].F[i]=%e\n", i,myOrg[who].F[i]);
    }
    
    //printf("ok 6 who %d\n",who);
    
    for(k=0;k<myOrg[who].ppicount;k++){
        i=myOrg[who].ppi_pair[k][0], j=myOrg[who].ppi_pair[k][1];
        CharNucSeqToAASeq(myOrg[who].genome+i*NUCSEQLEN,NUCSEQLEN,aaseq);
        CharNucSeqToAASeq(myOrg[who].genome+j*NUCSEQLEN,NUCSEQLEN,aaseq2);
        
        myOrg[who].pint[k] = (float) GetBindingP2(aaseq, myOrg[who].structid[i], aaseq2, myOrg[who].structid[j], myOrg[who].bmode[k], myParam->Tenv);
        
        
        if (allow_unfolded_states==1){
            myOrg[who].Gij[k] = myOrg[who].F[i]*myOrg[who].F[j]*myOrg[who].K[i][j]*myOrg[who].pint[k];
        }
        else{
            myOrg[who].Gij[k] = myOrg[who].F[i]*myOrg[who].F[j]*myOrg[who].K[i][j]*myOrg[who].pint[k]*myOrg[who].pnat[i]*myOrg[who].pnat[j];
        }
        
    }
    
    myOrg[who].birthrate = myParam->b0; //10E37;
    for (ii=0; ii<myOrg[who].ppicount; ii++) {
        if (allow_unfolded_states==2){
            i=myOrg[who].ppi_pair[ii][0], j=myOrg[who].ppi_pair[ii][1];
            myOrg[who].birthrate *=myOrg[who].Gij[ii]/(myOrg[who].pnat[i]*myOrg[who].pnat[j]);
        }
        else{
            
            if (hub_ID==10){
                //take into account only first pair
                if (ii==0){
                    myOrg[who].birthrate *=myOrg[who].Gij[ii];
                }
                else if (ii==1){
                }
                else{
                    myOrg[who].birthrate *=myOrg[who].F[ii]*myOrg[who].pnat[ii];
                }
            }
           else if (selection==0){
                myOrg[who].birthrate *= myOrg[who].Gij[ii]/(myOrg[who].K[i][j]*myOrg[who].pint[ii]);
            }
	       else{
                myOrg[who].birthrate *= myOrg[who].Gij[ii];
            }
        }
    }
    if (hub_ID==10){
        myOrg[who].birthrate *=myOrg[who].F[ii]*myOrg[who].pnat[ii];
    }
    
    
    C=0.0e0;
    for(i=0;i<curr_MAXGENES;i++) C+=myOrg[who].C[i];
    myOrg[who].birthrate/=(1.0e0+myParam->alpha*(C-((curr_MAXGENES+1)*0.1))*(C-((curr_MAXGENES+1)*0.1)));
    //  printf("myOrg[who].birthrate2 = %e\n", myOrg[who].birthrate);
    
    
    return 0;
}



void GetSequenceID(int who)
{
    int i,j, k, l, m=0;
    for(i=0;i<myOrg[who].genecount-1;i++){
        CharNucSeqToAASeq(myOrg[who].genome+i*NUCSEQLEN,NUCSEQLEN,aaseq);
        for(j=i+1;j<myOrg[who].genecount;j++) {
            CharNucSeqToAASeq(myOrg[who].genome+j*NUCSEQLEN,NUCSEQLEN,aaseq2);
            l=0;
            for(k=0;k<AASEQLEN;k++) if(aaseq[k]==aaseq2[k]) l++;
            myOrg[who].SeqID[m]=l;
            m++;
        }
    }
}
