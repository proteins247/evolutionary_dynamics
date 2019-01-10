/* generateStableGene.c
 * Inherited from Amy by way of Rostam.
 * Modified a bit (Victor):
 * - RNG
 * - 
 */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "zlib.h"
#include "../gencode.h"
#include "../latticelib.h"
#include "../rng.h"
#include "general.h"
#include "structurelib.h"

//Generates stable AA sequence of conf specified in argv[1]
int main(int argc, char *argv[])
{
    if (argc != 3)
    {
	printf("Usage: generateStableGene CONF SEED \n");
	exit(-1);
    }

    ReadCommondata();
    char aaseqL[100], nucseqL[100];
    int aaseqN[100], nucseqNA[100], nucseqNC[100], status;
    double Tenv = 1.0, targetPnat = 0.99;

    FILE *out;

    int confC = atoi(argv[1]);
    set_threefry_array(atoi(argv[2]), 0, 1, 2);

    // Generate a random gene with no stop codons
    status = -1;
    while (status == -1)
    {
	CreateRandomNucSequence(nucseqNA,NUCSEQLEN);
	CopySeq(nucseqNC, nucseqNA, NUCSEQLEN);
	status = NucSeqToAASeq(nucseqNA, NUCSEQLEN, aaseqN);
	//printf("stop codon encountered\n"); 
    }

    // Get initial Pnat
    double pnatC = GetStructurePnat(aaseqN, Tenv, confC);

    // Find stabilizing mutations 
    int ii;
    for (ii = 0; ii < 3000; ii++)
    {
	PointMutateNucSequence(nucseqNA, NUCSEQLEN);
	status = NucSeqToAASeq(nucseqNA, NUCSEQLEN, aaseqN);
	if (status == -1)
	{
	    CopySeq(nucseqNA, nucseqNC, NUCSEQLEN);
	} // Reject mutations that introduce a stop codon
	else
	{
	    double pnatA = GetStructurePnat(aaseqN, Tenv, confC);
	    if (pnatA >=  pnatC)
	    {
		CopySeq(nucseqNC, nucseqNA, NUCSEQLEN);
		pnatC = pnatA;
	    } // Accept stabilizing mutation
	    else
	    {
		CopySeq(nucseqNA, nucseqNC, NUCSEQLEN);
	    } // Reject destabilizing mutation
	    if (pnatC > targetPnat)
	    {
		break;
	    } // Break if target Pnat is reached
	}
    } // end attempts

    status = NucSeqToAASeq(nucseqNC, NUCSEQLEN, aaseqN);
    pnatC = GetStructurePnat(aaseqN, Tenv, confC);
    PrintAACodeSequence(aaseqL, aaseqN, AASEQLEN); // Convert numeric AA sequence to letter AA sequence
    PrintNucCodeSequence(nucseqL, nucseqNC, NUCSEQLEN); // Convert numeric nuc sequence to letter nuc sequence

    char fname[500];
    if (pnatC > targetPnat) // Save sequence to file
    {
	printf("%d\t%f\t\n%s\n%s\n", confC, pnatC, nucseqL, aaseqL);
	sprintf(fname, "%d.seq", confC); 
	out = fopen(fname, "w");
	fprintf(out, "%f\n", pnatC); 
	fprintf(out, "%s\n", nucseqL); 
    }
    else fprintf(stderr, "%d\t%f\t%s err \n", confC, pnatC, nucseqL); // Print sequence that didn't reach target Pnat to screen
    return 0;
}
