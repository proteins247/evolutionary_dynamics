#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>

#include"structurelib.h"
#include"../latticelib.h"

int CoordMatrix[NUMCONF][81];
int AllFaces[NUMCONF*6*4][9];
char ContactMatrixA[NUMCONF][32]; // contact matrix
char ContactMatrixB[NUMCONF][32]; 

void ReadCoordMatrix(char *filename)
{
FILE *infile;
int i,j;
int xmin,ymin,zmin;

infile=fopen(filename,"r");
if (!infile) {fprintf(stderr,"cannot read coord matrix\n"); exit(1);}

for(i=0;i<NUMCONF;i++)
for(j=0;j<81;j++)
fscanf(infile,"%d",&CoordMatrix[i][j]);

fclose(infile);

// adjust to bottom left corner @ 0,0,0
for(i=0;i<NUMCONF;i++)
{
xmin=ymin=zmin=100000;

for(j=0;j<27;j++)
{
if (CoordMatrix[i][3*j+0]<xmin) xmin=CoordMatrix[i][3*j+0];
if (CoordMatrix[i][3*j+1]<ymin) ymin=CoordMatrix[i][3*j+1];
if (CoordMatrix[i][3*j+2]<zmin) zmin=CoordMatrix[i][3*j+2];
}

for(j=0;j<27;j++)
{
CoordMatrix[i][3*j+0]-=xmin;
CoordMatrix[i][3*j+1]-=ymin;
CoordMatrix[i][3*j+2]-=zmin;
}

}
//printf("coord10000.dat read\n");
return;
}


void GetSquare(int structid, int face, int Square[3][3])
{
int i,j,k,ix,iy,flag=0;
int list[9];

if (face==0) // x=0
{
flag = 1;
k=0;
for(i=0;i<27;i++)
if (CoordMatrix[structid][3*i+0]==0) { list[k]=i; k++; }
if (k!=9) { fprintf(stderr,"getsquare error %d\n",k);exit(1); }
for(ix=0;ix<3;ix++)
for(iy=0;iy<3;iy++)
for(j=0;j<9;j++)
if ((CoordMatrix[structid][3*list[j]+1]==ix) && (CoordMatrix[structid][3*list[j]+2]==iy))
{ Square[ix][iy]=list[j]; }
}

if (face==3) // x=2
{
flag = 1;
k=0;
for(i=0;i<27;i++)
if (CoordMatrix[structid][3*i+0]==2) { list[k]=i; k++; }
if (k!=9) { fprintf(stderr,"getsquare error %d\n",k);exit(1); }
for(ix=0;ix<3;ix++)
for(iy=0;iy<3;iy++)
for(j=0;j<9;j++)
if ((CoordMatrix[structid][3*list[j]+1]==ix) && (CoordMatrix[structid][3*list[j]+2]==iy))
{ Square[ix][iy]=list[j]; }
}

if (face==1) // y=0
{
flag = 1;
k=0;
for(i=0;i<27;i++)
if (CoordMatrix[structid][3*i+1]==0) { list[k]=i; k++; }
if (k!=9) { fprintf(stderr,"getsquare error %d\n",k);exit(1); }
for(ix=0;ix<3;ix++)
for(iy=0;iy<3;iy++)
for(j=0;j<9;j++)
if ((CoordMatrix[structid][3*list[j]+0]==ix) && (CoordMatrix[structid][3*list[j]+2]==iy))
{ Square[ix][iy]=list[j]; }
}

if (face==4) // y=2
{
flag = 1;
k=0;
for(i=0;i<27;i++)
if (CoordMatrix[structid][3*i+1]==2) { list[k]=i; k++; }
if (k!=9) { fprintf(stderr,"getsquare error %d\n",k);exit(1); }
for(ix=0;ix<3;ix++)
for(iy=0;iy<3;iy++)
for(j=0;j<9;j++)
if ((CoordMatrix[structid][3*list[j]+0]==ix) && (CoordMatrix[structid][3*list[j]+2]==iy))
{ Square[ix][iy]=list[j]; }
}

if (face==2) // z=0
{
flag = 1;
k=0;
for(i=0;i<27;i++)
if (CoordMatrix[structid][3*i+2]==0) { list[k]=i; k++; }
if (k!=9) { fprintf(stderr,"getsquare error %d\n",k);exit(1); }
for(ix=0;ix<3;ix++)
for(iy=0;iy<3;iy++)
for(j=0;j<9;j++)
if ((CoordMatrix[structid][3*list[j]+0]==ix) && (CoordMatrix[structid][3*list[j]+1]==iy))
{ Square[ix][iy]=list[j]; }
}

if (face==5) // z=0
{
flag = 1;
k=0;
for(i=0;i<27;i++)
if (CoordMatrix[structid][3*i+2]==2) { list[k]=i; k++; }
if (k!=9) { fprintf(stderr,"getsquare error %d\n",k);exit(1); }
for(ix=0;ix<3;ix++)
for(iy=0;iy<3;iy++)
for(j=0;j<9;j++)
if ((CoordMatrix[structid][3*list[j]+0]==ix) && (CoordMatrix[structid][3*list[j]+1]==iy))
{ Square[ix][iy]=list[j]; }
}

if (!flag) { fprintf(stderr,"wrong face\n"); exit(1); }

return;
}


void RotateSquare(int code, int dest[3][3], int src[3][3])
{
//int ix,iy;
int flag;

if (code==0)
{
flag =1;

dest[0][0] = src[0][0];
dest[0][1] = src[0][1];
dest[0][2] = src[0][2];

dest[1][0] = src[1][0];
dest[1][1] = src[1][1];
dest[1][2] = src[1][2];

dest[2][0] = src[2][0];
dest[2][1] = src[2][1];
dest[2][2] = src[2][2];
}


if (code==1)
{
flag =1;

dest[0][0] = src[2][0]; 
dest[0][1] = src[1][0];
dest[0][2] = src[0][0];

dest[1][0] = src[2][1];
dest[1][1] = src[1][1];
dest[1][2] = src[0][1];

dest[2][0] = src[2][2];
dest[2][1] = src[1][2];
dest[2][2] = src[0][2];
}

if (code==2) // upside down
{
flag =1;

dest[0][0] = src[2][2];
dest[0][1] = src[2][1];
dest[0][2] = src[2][0];

dest[1][0] = src[1][2];
dest[1][1] = src[1][1];
dest[1][2] = src[1][0];

dest[2][0] = src[0][2];
dest[2][1] = src[0][1];
dest[2][2] = src[0][0];
}

if (code==3)
{
flag =1;

dest[0][0] = src[0][2];
dest[0][1] = src[1][2];
dest[0][2] = src[2][2];

dest[1][0] = src[0][1];
dest[1][1] = src[1][1];
dest[1][2] = src[2][1];

dest[2][0] = src[0][0];
dest[2][1] = src[1][0];
dest[2][2] = src[2][0];
}


if (!flag) { fprintf(stderr,"wrong rotation code\n"); exit(1); }
return;
}

void MakeAllFaces()
{
int i,j,k,ii,jj;
int face1[3][3];
int face2[3][3];

for(i=0;i<NUMCONF;i++)
for(j=0;j<6;j++)
for(k=0;k<4;k++)
{

GetSquare(i,j, face2);
RotateSquare(k,face1,face2);

for(ii=0;ii<3;ii++) 
{ 
   for(jj=0;jj<3;jj++) 
   { 
      AllFaces[24*i+4*j+k][3*ii+jj] = face1[ii][jj]; 
   }
}
}
return;
}

void MakeContactMatrix()		//unecessary, can have it read in. in general it seems that this file was created to not have to read in contact matrix
{
  int kk, ii, jj, mm, c=0, cont=0;
  int tmp1[3], tmp2[3];

  for(kk=0; kk<NUMCONF; kk++)
  {
   c = 0;
   for(ii=0; ii<81; ii+=3)
   {
       for(mm=0; mm<3; mm++)     { tmp1[mm] = CoordMatrix[kk][ii+mm]; } //printf("%d ", tmp1[mm]); }
       for(jj=ii+9; jj<81; jj+=3)
       {
          cont = 0;
          for(mm=0; mm<3; mm++)  { tmp2[mm] = CoordMatrix[kk][jj+mm]; } //printf("%d ", tmp2[mm]);     
          for(mm=0; mm<3; mm++)    cont += abs(tmp2[mm] - tmp1[mm]);
          if(cont==1)
          {
             ContactMatrixA[kk][c]=jj/3;
             ContactMatrixB[kk][c]=ii/3;
             c++;
          }
       }
    }
   }
   return;
}

void ReadCommondata()
{
char rootdir[100], file[100];
sprintf(rootdir,"%s","/n/home12/rrazban/code");
sprintf(file, "%s/commondata/LPforms/coord10000.dat", rootdir);
ReadCoordMatrix(file);
MakeAllFaces();
MakeContactMatrix();
sprintf(file, "%s/commondata/MJ96/energy.dat", rootdir);
ReadEnergyMatrix(file);
return;
}
