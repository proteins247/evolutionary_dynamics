#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>

#include"structurelib.h"
#include"../latticelib.h"

int CoordMatrix[NUMCONF][81];
int AllFaces[NUMCONF*6*4][9];
/* char ContactMatrixA[NUMCONF][32]; // contact matrix */
/* char ContactMatrixB[NUMCONF][32];  */

void ReadCoordMatrix(char *filename)
{
    FILE *infile;
    int i, j;
    int xmin, ymin, zmin;

    infile = fopen(filename,"r");
    if (!infile)
    {
	fprintf(stderr,"cannot read coord matrix\n");
	exit(1);
    }

    for (i = 0; i < NUMCONF; i++)
	for (j = 0; j < 81; j++)
	    fscanf(infile, "%d", &CoordMatrix[i][j]);

    fclose(infile);

    // adjust to bottom left corner @ 0,0,0
    /* 
     * Note by Victor: as far as I can tell, every lattice protein 
     * in the coordinate file coord10000.dat has a vertex at (0, 0, 0),
     * so this seems like a useless operation. 
     */
    for (i = 0; i < NUMCONF; i++)
    {
	xmin = ymin = zmin = 100000;

	for (j = 0; j < 27; j++)
	{
	    if (CoordMatrix[i][3 * j + 0] < xmin)
		xmin = CoordMatrix[i][3 * j + 0];
	    if (CoordMatrix[i][3 * j + 1] < ymin)
		ymin = CoordMatrix[i][3 * j + 1];
	    if (CoordMatrix[i][3 * j + 2] < zmin)
		zmin = CoordMatrix[i][3 * j + 2];
	}

	for (j = 0; j < 27; j++)
	{
	    CoordMatrix[i][3 * j + 0] -= xmin;
	    CoordMatrix[i][3 * j + 1] -= ymin;
	    CoordMatrix[i][3 * j + 2] -= zmin;
	}
    }
    // printf("coord10000.dat read\n");
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
/* 
 * Note by Victor: By "unecessary," I think it's meant that the file 
 * contact10000.dat exists, and there is a function to read the file
 * in latticelib. Here, ContactMatrixA and ContactMatrixB are filled
 * in using the coordinates of the 10000 lattice proteins.
 */
{
    int kk, ii, jj, mm;
    int coord1[3], coord2[3];

    for (kk = 0; kk < NUMCONF; kk++) /* For each lattice protein */
    {
	int c = 0;
	for (ii = 0; ii < 81; ii += 3) /* For each coordinate */
	{
	    /* Read x, y, z into coord1 */
	    for (mm = 0; mm < 3; mm++)
	    {
		coord1[mm] = CoordMatrix[kk][ii+mm];
	    } //printf("%d ", coord1[mm]); }

	    /* Beginning with i+3rd coordinate (i+1 and i+2 can't make contact) */
	    for (jj = ii + 9; jj < 81; jj += 3)
	    {
		/* Read x, y, z into coord2 */
		for (mm = 0; mm < 3; mm++)
		{
		    coord2[mm] = CoordMatrix[kk][jj+mm];
		} //printf("%d ", coord2[mm]);

		int contact = 0;
		for (mm = 0; mm < 3; mm++)
		{
		    contact += abs(coord2[mm] - coord1[mm]);
		}
		if (contact == 1)
		{
		    ContactMatrixA[kk][c] = jj / 3;
		    ContactMatrixB[kk][c] = ii / 3;
		    c++;
		}
	    }
	}
    }

    return;
}

void ReadCommondata(char * location, char * energyName)
{
    char path[500];

    sprintf(path, "%s/LPforms/coord10000.dat", location);
    ReadCoordMatrix(path);
    MakeAllFaces();
    MakeContactMatrix();

    if (energyName)		/* if not void */
    {
	sprintf(path, "%s/MJ96/%s.dat", location, energyName);
    }
    else
    {
	sprintf(path, "%s/MJ96/energy.dat", location);
    }
    ReadEnergyMatrix(path);
    return;
}

void PrintFoldedConformation(int structid, char * conformation)
{
    int i;
    for (i = 0; i < AASEQLEN - 1; ++i)
    {
	int x_difference = CoordMatrix[structid][3 * i + 3]
	    - CoordMatrix[structid][3 * i];
	int y_difference = CoordMatrix[structid][3 * i + 4]
	    - CoordMatrix[structid][3 * i + 1];
	int z_difference = CoordMatrix[structid][3 * i + 5]
	    - CoordMatrix[structid][3 * i + 2];
	if (x_difference == 1)
	    conformation[i] = 'F';
	else if (x_difference == -1)
	    conformation[i] = 'B';
	else if (y_difference == 1)
	    conformation[i] = 'R';
	else if (y_difference == -1)
	    conformation[i] = 'L';
	else if (z_difference == 1)
	    conformation[i] = 'U';
	else if (z_difference == -1)
	    conformation[i] = 'D';
    }
    conformation[AASEQLEN - 1] = '\0';
    return;
}
