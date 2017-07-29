/*************************
 Modified by Orit Peleg
 opeleg@fas.harvard.edu
 on 2/26/13
 ************************/


/*
lattice protein binding energy routines
*/
#ifndef __BINDINGLIB_H__
#define __BINDINGLIB_H__

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include "define.h"

extern int AllFaces[NUMCONF*6*4][9];

/*
for all structures, read the residue numbers
at each position on all 6 faces in 4 orientations per face
*/
void ReadAllSurfaces(char *filename);

/*
mirror image of the interacting surface
surface is an integer array[9]
*/
void MirrorWall(int *dest, int *src);

/*
returns binding energy of seq1 and seq2.
The conformations are struct1 and struct2 [0..NUMCONF-1],
the interacting surfaces are surface1 and surface2 [0..5],
and the orientation is rotate [0..3]
ReadEnergyMatrix from latticelib must be called before 
using this function
*/
double GetBindingEnergy(int *seq1, int struct1, int surface1, 
                        int *seq2, int struct2, int surface2, int rotate);

/*
returns the probability of seq1/struct1 and seq2/struct2 to be in a specific
bound state at temperature T. Assumes rigid docking (no conformational change
upon binding).
Purely two-protein model, no competition with other proteins
*/
double GetBindingP(int *seq1, int struct1, int *seq2, int struct2, int *bmode, double T);
double GetBindingP2(int *seq1, int struct1, int *seq2, int struct2, int bmode, double T);
double GetBindingP3(int *seq1, int struct1, int surf1, int *seq2, int struct2, int surf2, int *bmode, double T);
double GetBindingK(int *seq1, int struct1, int *seq2, int struct2, double T);
double GetBindingCh(int *seq1, int struct1,double T);
double GetBindingCh_OP(int *seq1, int struct1,double T);
double GetBindingK2(int *seq1, int struct1, int *seq2, int struct2, int bmode, double T);

double GetBindingEnergy2(int *seq1, int struct1, int *seq2, int struct2, int bmode);

double GetBindingEnergyHydro(int *seq1, int struct1, int *seq2, int struct2, int bmode);

double GetBindingEnergyCharge(int *seq1, int struct1, int *seq2, int struct2, int bmode);

double GetBindingEnergy2(int *seq1, int struct1, int *seq2, int struct2, int bmode);

void GetSurfaceAAPositions(int struct1, int struct2, int bmode, int *surface1, int *surface2);

void UpdateSurfaceAAProps(int *seq, int structid, int bmode);
void GetSurfaceAA(int *aaseq_hub, int *aaseq_partner, int *aaseq_surface_hub, int *aaseq_surface_partner, int struct_hub, int struct_partner, int bmode);
void GetSingleSurfaceAA(int *aaseq_surface, int *aaseq_full, int struct_id, int face_id,int bmode, int mirror, int rotate);

void GetDoubleSurfaceAA(int *seq1, int struct1, int *seq2, int struct2, int bmode,int *surface1, int *surface2);




#endif
