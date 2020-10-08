/* 
 * translation_schedule.c : 
 *   Give translation schedule given aa sequence.
 */
#include <unistd.h>
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../gencode.h"

static const int PARSE_ERROR = 1;
Codon STOP_CODON = N_UAA;


int main(int argc, char *argv[])
{
    char energyFile[100] = "times_scale100x.dat";
    char energyPath[1000];

    // Begin option handling (https://linux.die.net/man/3/getopt)
    int c, option_index = 0;
    static struct option long_options[] = {
	{"energy", required_argument, NULL, 'e'},
	{NULL, 0, NULL, 0}
    };
    
    while ((c = getopt_long(argc, argv, "e:h",
			    long_options, &option_index)) != -1)
    {
	switch (c)
	{
	case 'e':
	    strcpy(energyFile, optarg);
	    break;
	case 'h':
	    fprintf(stderr, "help message (not implemented)\n");
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

    char * aaLetterSequence;
    if (optind + 1 == argc)
    {
	aaLetterSequence = argv[optind];
    }
    else
    {
	fprintf(stderr, "Expected positional argument: AASEQUENCE\n");
	exit(PARSE_ERROR);
    }

    /* Load energy file */
    sprintf(energyPath, "%s/%s", getenv("FOLDEVO_SHARE"), energyFile);
    unsigned int times[821];
    ReadTranslationTimes(energyPath, times);
    
    /* Process aa sequence */
    int protein_length = strlen(aaLetterSequence);
    if (protein_length > 100)
    {
	fprintf(stderr, "Input aa sequence too long (>100)");
	exit(PARSE_ERROR);
    }
    AminoAcid aaSequence[100];
    int nucSequence[300];
    Codon codonSequence[100];
    LetterToAASeq(aaLetterSequence, aaSequence, protein_length);
    AASeqToNucSeq(aaSequence, nucSequence, protein_length, 0);
    NucSeqToCodonSeq(nucSequence, protein_length * 3, codonSequence);
    
    int index = 5;
    for (; index < protein_length; ++index)
    {
	printf("%d ", times[codonSequence[index]]);
    }
    printf("%d\n", times[STOP_CODON]);

    return 0;
}
