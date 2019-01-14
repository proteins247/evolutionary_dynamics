/* generateStableGene.c : Generates stable sequence for specified conf
 * Inherited from Amy by way of Rostam.
 * Modified a bit (Victor):
 * - RNG
 * - Takes in commandline arguments
 */
#include <unistd.h>
#include <getopt.h>
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

static const int PARSE_ERROR = 1;
static const int TARGET_FAILED = 2;


int main(int argc, char *argv[])
{
    double Tenv = 1.0;
    double targetPnat = 0.99;
    char * dataDir = getenv("FOLDEVO_SHARE");
    char energyName[100] = "energy";
    static int subtractMean = 0;

    // Begin option handling (https://linux.die.net/man/3/getopt)
    int c, option_index = 0;
    static struct option long_options[] = {
	{"energy", required_argument, NULL, 'e'},
	{"target-pnat", required_argument, NULL, 'p'},
	{"subtract-mean", no_argument, &subtractMean, 1},
	{"temperature", required_argument, NULL, 't'},
	{NULL, 0, NULL, 0}
    };
    
    while ((c = getopt_long(argc, argv, "e:t:p:",
			    long_options, &option_index)) != -1)
    {
	switch (c)
	{
	case 'e':
	    strcpy(energyName, optarg);
	    break;
	case 'p':
	    targetPnat = atof(optarg);
	    break;
	case 't':
	    Tenv = atof(optarg);
	    break;
	case '?':
	    // getopt_long will print an error message
	    fprintf(stderr, "Invalid argument %s\n", argv[option_index]);
	    exit(PARSE_ERROR);
	default:
	    // we shouldn't end up here i think
	    ;
	}
    }

    int confC;
    uint64_t seed;
    if (optind + 2 == argc)
    {
	confC = atoi(argv[optind++]);
	seed = atol(argv[optind]);
    }
    else
    {
	fprintf(stderr, "Expected two positional arguments: CONF SEED\n");
	exit(PARSE_ERROR);
    }

    ReadCommondata(dataDir, energyName);

    if (subtractMean)
    {
	subtractMeanFromEnergyMatrix();
    }
    /* End argument parsing */

    char aaseqL[100], nucseqL[100];
    int aaseqN[100], nucseqNA[100], nucseqNC[100], status;

    FILE *out;

    /* Set seed */
    set_threefry_array(seed, 1337, 42, 69);

    // Generate a random gene with no stop codons
    status = -1;
    while (status == -1)
    {
	CreateRandomNucSequence(nucseqNA, NUCSEQLEN);
	CopySeq(nucseqNC, nucseqNA, NUCSEQLEN);
	status = NucSeqToAASeq(nucseqNA, NUCSEQLEN, aaseqN);
    }

    // Get initial Pnat
    double pnatC = GetStructurePnat(aaseqN, Tenv, confC);

    // Find stabilizing mutations 
    int ii;
    for (ii = 0; ii < 10000; ii++)
    {
	PointMutateNucSequence(nucseqNA, NUCSEQLEN);
	status = NucSeqToAASeq(nucseqNA, NUCSEQLEN, aaseqN);
	if (status == -1)
	{
            /* Reject mutations that introduce a stop codon */
	    CopySeq(nucseqNA, nucseqNC, NUCSEQLEN);
	}
	else
	{
	    double pnatA = GetStructurePnat(aaseqN, Tenv, confC);
	    if (pnatA >=  pnatC)
	    {
		// Accept stabilizing mutation
		CopySeq(nucseqNC, nucseqNA, NUCSEQLEN);
		pnatC = pnatA;
	    }
	    else
	    {
		// Reject destabilizing mutation
		CopySeq(nucseqNA, nucseqNC, NUCSEQLEN);
	    }
	    if (pnatC > targetPnat)
	    {
		// Break if target Pnat is reached
		break;
	    } 
	}
    } // end attempts

    status = NucSeqToAASeq(nucseqNC, NUCSEQLEN, aaseqN);
    pnatC = GetStructurePnat(aaseqN, Tenv, confC);
    // Convert numeric AA sequence to letter AA sequence
    PrintAACodeSequence(aaseqL, aaseqN, AASEQLEN);
    // Convert numeric nuc sequence to letter nuc sequence
    PrintNucCodeSequence(nucseqL, nucseqNC, NUCSEQLEN);

    double energy = SequenceEnergy(aaseqN, confC);

    // Save sequence to file if successful
    char fname[500];
    char conformation[30];
    PrintFoldedConformation(confC, conformation);
    if (pnatC > targetPnat)
    {
	/* First print information to screen. */
	printf("%d\t%f\t\n", confC, pnatC);
	printf("Folded conf energy: %f\n", energy);
	printf("%s\n", nucseqL);
	printf("%s\n", aaseqL);
	printf("Folded conformation: %s\n", conformation);

	sprintf(fname, "%d.sourceme", confC); 
	out = fopen(fname, "w");
	fprintf(out, "SEED=%lu\n", seed);
	fprintf(out, "PNAT=%f\n", pnatC);
	fprintf(out, "TDESIGN=%f\n", Tenv);
	fprintf(out, "NATIVEENERGY=%f\n", energy);
	fprintf(out, "NUCSEQ=%s\n", nucseqL);
	fprintf(out, "AASEQ=%s\n", aaseqL);
	fprintf(out, "CONFORMATION=%s\n", conformation);
	fclose(out);
    }
    else
    {
	fprintf(stderr, "Failed to reach target Pnat.\n");
	fprintf(stderr, "  Maybe try another random number seed.\n");
	fprintf(stderr, "%d\t%f\t%s\n", confC, pnatC, nucseqL); // Print sequence that didn't reach target Pnat to screen
	exit(TARGET_FAILED);
    }
    return 0;
}
