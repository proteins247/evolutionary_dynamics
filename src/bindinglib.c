/*************************
 Modified by Orit Peleg
 opeleg@fas.harvard.edu
 on 2/26/13

 Modified by Rostam Razban
 rrazban@g.harvard.edu
 on 8/1/16
 ************************/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include"latticelib.h"
#include"bindinglib.h"

int AllFaces[NUMCONF*6*4][9];

void ReadAllSurfaces(char *filename)
{
  FILE *infile;
  int i,j,k,m,t1,t2,t3,t4;
  infile=fopen(filename,"r");
  if (!infile) {fprintf(stderr,"cannot read face list %s\n",filename); exit(1);}

  for(i=0;i<NUMCONF;i++)
    for(j=0;j<6;j++)
      for(k=0;k<4;k++)
      {
        fscanf(infile,"%d %d %d",&t1,&t2,&t3); // structid, face (0-5), rotation (0,1,2,3)
        for(m=0;m<9;m++)
        {
          fscanf(infile,"%d", &t4);
          AllFaces[24*i+4*j+k][m]=t4;
        }
      }

  return;
}


void MirrorWall(int *dest, int *src)
{
	dest[0]=src[6];
	dest[1]=src[7];
	dest[2]=src[8];

	dest[3]=src[3];
	dest[4]=src[4];
 	dest[5]=src[5];

 	dest[6]=src[0];
 	dest[7]=src[1];
	dest[8]=src[2];

	return;
}


double GetBindingEnergy(int *seq1, int struct1, int face1, int *seq2, int struct2, int face2, int rotate){
	int k, surfacetmp[9],surface1[9], surface2[9];
	double e=0;

	for(k=0; k<9; k++) surface1[k] = seq1[(int) AllFaces[24*struct1+4*face1+0][k]];
	for(k=0; k<9; k++) surfacetmp[k] = seq2[(int) AllFaces[24*struct2+4*face2+rotate][k]];
	if ((face1-face2)%2 == 0){
		MirrorWall(surface2, surfacetmp);
	}
	else{
		memcpy(surface2, surfacetmp, (sizeof(surface2)+1) * sizeof(*surface2));	//why dont i do this consistently?
	}
	for(k=0; k<9; k++) e+=EnergyMatrix[surface1[k]][surface2[k]];

	return e;
}

// rigid docking!
double GetBindingP(int *seq1, int struct1, int *seq2, int struct2, int *bmode, double T)
{
  int face1, face2, rotate;
  int k, surfacetmp[9],surface1[9], surface2[9];
  double e, z=0, emin=1e10;

  for(face1=0;face1<6;face1++)
  {
    for(k=0; k<9; k++) surface1[k] = seq1[(int) AllFaces[24*struct1+4*face1+0][k]];
    for(face2=0;face2<6;face2++)
      for(rotate=0;rotate<4;rotate++)
      {
        e=0;
        for(k=0; k<9; k++) surfacetmp[k] = seq2[(int) AllFaces[24*struct2+4*face2+rotate][k]];
    	if ((face1-face2)%2 == 0){		//this looks highly suspect... should test out with small example so that i can confirm
			MirrorWall(surface2, surfacetmp);
		}
		else{
			memcpy(surface2, surfacetmp, sizeof(surfacetmp));
		}
  for(k=0; k<9; k++) e+=EnergyMatrix[surface1[k]][surface2[k]];
        z+=exp(-e/T);
        if (e<emin) {
          emin=e;
          *bmode=face1*24+face2*4+rotate;
        }
      }
  }

  //if(z<0.00001) { fprintf(stderr,"Error!!! Partition Function z = %8.3lf\n", z); exit(1); }

  return exp(-emin/T)/z;
}

// rigid docking!
double GetBindingP2(int *seq1, int struct1, int *seq2, int struct2, int bmode, double T)
{
  int face1, face2, rotate;
  int k, surfacetmp[9],surface1[9], surface2[9];
  double e, emin, z=0.0e0;

  rotate=bmode%4;
  face1=bmode/24;
  face2=(bmode%24)/4;
  for(k=0; k<9; k++) surface1[k] = seq1[(int) AllFaces[24*struct1+4*face1+0][k]];
  for(k=0; k<9; k++) surfacetmp[k] = seq2[(int) AllFaces[24*struct2+4*face2+rotate][k]];
	if ((face1-face2)%2 == 0){
		MirrorWall(surface2, surfacetmp);
    }
	else{
		memcpy(surface2, surfacetmp, sizeof(surfacetmp));
	}

  emin=0.0e0;
  for(k=0; k<9; k++) emin+=EnergyMatrix[surface1[k]][surface2[k]];

  for(face1=0;face1<6;face1++) {
    for(k=0; k<9; k++) surface1[k] = seq1[(int) AllFaces[24*struct1+4*face1+0][k]];
    for(face2=0;face2<6;face2++)
      for(rotate=0;rotate<4;rotate++)
      {
        e=0.0e0;
        for(k=0; k<9; k++) surfacetmp[k] = seq2[(int) AllFaces[24*struct2+4*face2+rotate][k]];
		if ((face1-face2)%2 == 0){
            MirrorWall(surface2, surfacetmp);
        }
        else{
			memcpy(surface2, surfacetmp, sizeof(surfacetmp));
        }

	    for(k=0; k<9; k++) e+=EnergyMatrix[surface1[k]][surface2[k]];
        z+=exp(-e/T);
      }
  }

//  if(z<0.00001) { fprintf(stderr,"Error!!! Partition Function z = %8.3lf\n", z); exit(1); }

  return exp(-emin/T)/z;
}

/////////////////////////////////////////////////////////////////////////////////////////

void GetDoubleSurfaceAA(int *seq1, int struct1, int *seq2, int struct2, int bmode, int *surface1, int *surface2){
    int face1, face2, rotate;
    int k, surfacetmp[9];
    
    rotate=bmode%4;
    face1=bmode/24;
    face2=(bmode%24)/4;
    for(k=0; k<9; k++) surface1[k] = seq1[(int) AllFaces[24*struct1+4*face1+0][k]];
    for(k=0; k<9; k++) surfacetmp[k] = seq2[(int) AllFaces[24*struct2+4*face2+rotate][k]];
	if ((face1-face2)%2 == 0){
		MirrorWall(surface2, surfacetmp);
	}
	else{
		memcpy(surface2, surfacetmp, (sizeof(surface2)+1) * sizeof(*surface2));
	}
}

/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////

// rigid docking!
void GetSurfaceAAPositions(int struct1, int struct2, int bmode, int *surface1, int *surface2)
{
    int face1, face2, rotate;
    int k, surfacetmp[9];

    
    rotate=bmode%4;
    face1=bmode/24;
    face2=(bmode%24)/4;
    for(k=0; k<9; k++) surface1[k] = (int) AllFaces[24*struct1+4*face1+0][k];
    for(k=0; k<9; k++) surfacetmp[k] = (int) AllFaces[24*struct2+4*face2+rotate][k];
	if ((face1-face2)%2 == 0){
		MirrorWall(surface2, surfacetmp);
	}
	else{
		memcpy(surface2, surfacetmp, sizeof(surfacetmp));
	}
}


/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////

// rigid docking!
double GetBindingEnergy2(int *seq1, int struct1, int *seq2, int struct2, int bmode){
    int face1, face2, rotate;
    int k, surfacetmp[9],surface1[9], surface2[9];
    double emin;
    
    rotate=bmode%4;
    face1=bmode/24;
    face2=(bmode%24)/4;
    for(k=0; k<9; k++) surface1[k] = seq1[(int) AllFaces[24*struct1+4*face1+0][k]];
    for(k=0; k<9; k++) surfacetmp[k] = seq2[(int) AllFaces[24*struct2+4*face2+rotate][k]];
	if ((face1-face2)%2 == 0){
		MirrorWall(surface2, surfacetmp);
	}
	else{
	//	memcpy(surface2, surfacetmp, (sizeof(surface2)+1) * sizeof(*surface2));
		memcpy(surface2, surfacetmp, sizeof(surface2));
	}
    emin=0.0e0;
    for(k=0; k<9; k++){
		 emin+=EnergyMatrix[surface1[k]][surface2[k]];
    }
    return emin;
}

/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////

double GetBindingEnergyHydro(int *seq1, int struct1, int *seq2, int struct2, int bmode){
    int face1, face2, rotate;
    int k, surfacetmp[9],surface1[9], surface2[9];
    double emin;
    
    int surface1_t[9], surface2_t[9];
    
    
    rotate=bmode%4;
    face1=bmode/24;
    face2=(bmode%24)/4;
    for(k=0; k<9; k++) surface1[k] = seq1[(int) AllFaces[24*struct1+4*face1+0][k]];
    for(k=0; k<9; k++) surfacetmp[k] = seq2[(int) AllFaces[24*struct2+4*face2+rotate][k]];
	if ((face1-face2)%2 == 0){
		MirrorWall(surface2, surfacetmp);
	}
	else{
		memcpy(surface2, surfacetmp, (sizeof(surface2)+1) * sizeof(*surface2));
	}

    ConvertAAtoHydro(seq1, surface1_t, AASEQLEN);
    ConvertAAtoHydro(seq2, surface2_t, AASEQLEN);
    
    emin=0.0e0;
    for(k=0; k<9; k++) {
        if ((surface1_t[k]==1) | (surface2_t[k]==1)){
            emin+=EnergyMatrix[surface1[k]][surface2[k]]-3.16;
        }
    }
    
    return emin;
}

/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////

double GetBindingEnergyCharge(int *seq1, int struct1, int *seq2, int struct2, int bmode){
    int face1, face2, rotate;
    int k, surfacetmp[9],surface1[9], surface2[9];
    double emin;
    
    int surface1_t[9], surface2_t[9];
    
    
    rotate=bmode%4;
    face1=bmode/24;
    face2=(bmode%24)/4;
    for(k=0; k<9; k++) surface1[k] = seq1[(int) AllFaces[24*struct1+4*face1+0][k]];
    for(k=0; k<9; k++) surfacetmp[k] = seq2[(int) AllFaces[24*struct2+4*face2+rotate][k]];
	if ((face1-face2)%2 == 0){
		MirrorWall(surface2, surfacetmp);
	}
	else{
		memcpy(surface2, surfacetmp, (sizeof(surface2)+1) * sizeof(*surface2));
	}

    ConvertAAtoCharge(seq1, surface1_t, AASEQLEN);
    ConvertAAtoCharge(seq2, surface2_t, AASEQLEN);
    
    emin=0.0e0;
    for(k=0; k<9; k++) {
        if (((surface1_t[k]==1) | (surface2_t[k]==-1)) | ((surface1_t[k]==-1) | (surface2_t[k]==1))){
            emin+=EnergyMatrix[surface1[k]][surface2[k]]-3.16;
        }
    }
    
    return emin;
}


/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////



void GetSurfaceAA(int *aaseq_hub, int *aaseq_partner, int *aaseq_surface_hub, int *aaseq_surface_partner,
                  int struct_hub, int struct_partner, int bmode){
    
    int face1, face2, rotate;
    int k;
    
    rotate=bmode%4;
    face1=bmode/24;
    face2=(bmode%24)/4;
    for(k=0; k<9; k++) aaseq_surface_hub[k] = aaseq_hub[(int) AllFaces[24*struct_hub+4*face1+0][k]];
    for(k=0; k<9; k++) aaseq_surface_partner[k] = aaseq_partner[(int) AllFaces[24*struct_partner+4*face2+rotate][k]];
}


void GetSingleSurfaceAA(int *aaseq_surface, int *aaseq_full, int struct_id, int face_id, int bmode, int mirror, int rotate){
    int k;
    int surface1[9];
    int roate_value;
    
    if (rotate==1){
        roate_value = bmode%4;
    }
    else{
        roate_value = 0;
    }
    
    for(k=0; k<9; k++) aaseq_surface[k] = aaseq_full[(int) AllFaces[24*struct_id+4*face_id+roate_value][k]];
    
    if (mirror == 1){
  //      MirrorWall(surface1, aaseq_surface);
        for(k=0; k<9; k++) aaseq_surface[k] = surface1[k];
    }
}




double GetBindingP3(int *seq1, int struct1, int surf1, int *seq2, int struct2, int surf2, int *bmode, double T)
{
  int face1, face2, rotate;
  int k, surfacetmp[9],surface1[9], surface2[9];
  double e, z=0, emin=1e10;

  *bmode=-1;

  for(face1=0;face1<6;face1++)
  {
    for(k=0; k<9; k++) surface1[k] = seq1[(int) AllFaces[24*struct1+4*face1+0][k]];
    for(face2=0;face2<6;face2++)
      for(rotate=0;rotate<4;rotate++)
      {
        e=0;
        for(k=0; k<9; k++) surfacetmp[k] = seq2[(int) AllFaces[24*struct2+4*face2+rotate][k]];
		if ((face1-face2)%2 == 0){
			MirrorWall(surface2, surfacetmp);
		}
		else{
			memcpy(surface2, surfacetmp, (sizeof(surface2)+1) * sizeof(*surface2));
		}

	    for(k=0; k<9; k++) e+=EnergyMatrix[surface1[k]][surface2[k]];
        //e*=1.2;
        z+=exp(-e/T);
        if (e<emin && face1 == surf1 && face2 == surf2) {
          emin=e;
          *bmode=face1*24+face2*4+rotate;
        }
      }
  }

  if(z<0.00001) { fprintf(stderr,"Error!!! Partition Function z = %8.3lf\n", z); exit(1); }

  return exp(-emin/T)/z;
}


double GetBindingK(int *seq1, int struct1, int *seq2, int struct2, double T){
  int face1, face2, rotate;
  int k, surfacetmp[9],surface1[9], surface2[9];
  double e, z=0, emin=1e10;

  for(face1=0;face1<6;face1++)
  {
    for(k=0; k<9; k++) surface1[k] = seq1[(int) AllFaces[24*struct1+4*face1+0][k]];
    for(face2=0;face2<6;face2++)
      for(rotate=0;rotate<4;rotate++)
      {
        e=0;
        for(k=0; k<9; k++) surfacetmp[k] = seq2[(int) AllFaces[24*struct2+4*face2+rotate][k]];
		if ((face1-face2)%2 == 0){
			MirrorWall(surface2, surfacetmp);
		}
		else{
			memcpy(surface2, surfacetmp, (sizeof(surface2)+1) * sizeof(*surface2));
		}

	    for(k=0; k<9; k++) e+=EnergyMatrix[surface1[k]][surface2[k]];
        z+=exp(-e/T);
        if (e<emin) emin=e;
      }
  }

  //if(z<1.0e-12) { fprintf(stderr,"Error!!! Partition Function z = %8.3lf\n", z); exit(1); }

  return z;
}


double GetBindingCh(int *seq1, int struct1, double T)
{
    int face1, rotate;
//	int surfacetmp[9];
    int k, surface1[9], Ch_surface[9];
    double e, z=0, emin=1e10;
    
    Ch_surface[0]=7; Ch_surface[1]=8; Ch_surface[2]=5;
    Ch_surface[3]=4; Ch_surface[4]=2; Ch_surface[5]=7;
    Ch_surface[6]=5; Ch_surface[7]=8; Ch_surface[8]=2;
    
    for(face1=0;face1<6;face1++){
        for(rotate=0;rotate<4;rotate++){
            
            for(k=0; k<9; k++){ surface1[k] = seq1[(int) AllFaces[24*struct1+4*face1+rotate][k]]; }
            
            e=0;
            for(k=0; k<9; k++){ e+=EnergyMatrix[surface1[k]][Ch_surface[k]]; }
            z+=exp(-e/T);
            if (e<emin) emin=e;
            
        }
    }
    //     fprintf(stderr,"Energy is %f \n", z);
    //if(z<1.0e-12) { fprintf(stderr,"Error!!! Partition Function z = %8.3lf\n", z); exit(1); }
    
    return z;
}



double GetBindingCh_OP(int *seq1, int struct1,double T){
    int face1, face2, rotate;
//	int surfacetmp[9];
    int k,surface1[9], surface2[9];
    double e, z=0, emin=1e10;
    
    
    //Tyr, Leu, Val, Leu, Phe, Tyr, Leu, Ala, Phe
    //7, 4, 5,
    //4, 2, 7,
    //4, 8, 2
    
    surface2[0] =7; surface2[1] =4; surface2[2] =5;
    surface2[3] =4; surface2[4] =2; surface2[5] =7;
    surface2[6] =4; surface2[7] =8; surface2[8] =2;
    
    for(face1=0;face1<6;face1++){
        for(k=0; k<9; k++) surface1[k] = seq1[(int) AllFaces[24*struct1+4*face1+0][k]];
        for(face2=0;face2<6;face2++)
            for(rotate=0;rotate<4;rotate++){
                e=0;
                //for(k=0; k<9; k++) surface2[k] = seq2[(int) AllFaces[24*struct2+4*face2+rotate][k]];
                for(k=0; k<9; k++) e+=EnergyMatrix[surface1[k]][surface2[k]];
                z+=exp(-e/T);
                if (e<emin) emin=e;
            }
    }
    
    //if(z<1.0e-12) { fprintf(stderr,"Error!!! Partition Function z = %8.3lf\n", z); exit(1); }
    
    return z;
}


// still need to do this: 
double GetBindingK2(int *seq1, int struct1, int *seq2, int struct2, int bmode, double T){
    int face1, face2, rotate;
    int k, surfacetmp[9],surface1[9], surface2[9];
    double z=0, emin=1e10;    
    
    rotate=bmode%4;
    face1=bmode/24;
    face2=(bmode%24)/4;
    for(k=0; k<9; k++) surface1[k] = seq1[(int) AllFaces[24*struct1+4*face1+0][k]];
    for(k=0; k<9; k++) surfacetmp[k] = seq2[(int) AllFaces[24*struct2+4*face2+rotate][k]];
	if ((face1-face2)%2 == 0){
		MirrorWall(surface2, surfacetmp);
	}
	else{
		memcpy(surface2, surfacetmp, (sizeof(surface2)+1) * sizeof(*surface2));
	}

    emin=0.0e0;
    for(k=0; k<9; k++) emin+=EnergyMatrix[surface1[k]][surface2[k]];
    z =  exp(-emin/T);
//    for(face1=0;face1<6;face1++)
//    {
//        for(k=0; k<9; k++) surfacetmp[k] = seq1[(int) AllFaces[24*struct1+4*face1+0][k]];
//        MirrorWall(surface1, surfacetmp);
//        for(face2=0;face2<6;face2++)
//            for(rotate=0;rotate<4;rotate++)
//            {
//                e=0;
//                for(k=0; k<9; k++) surface2[k] = seq2[(int) AllFaces[24*struct2+4*face2+rotate][k]];
//                for(k=0; k<9; k++) e+=EnergyMatrix[surface1[k]][surface2[k]];
//                z+=exp(-e/T);
//                if (e<emin) emin=e;
//            }
//    }
//    
//    //if(z<1.0e-12) { fprintf(stderr,"Error!!! Partition Function z = %8.3lf\n", z); exit(1); }
//    
    return z;
}
