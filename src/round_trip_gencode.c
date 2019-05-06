#include <stdio.h>
#include "gencode.h"

/* 
 * So if go from each amino acid to a codon and then codon to aa, we
 * should get the codon back.
 *
 * This function required making TripletToAA and AAtoCodon non-static
 */

int main()
{
    AminoAcid aa;
    int i;
    int error_count = 0;
    for (i = 0; i < 100; ++i)
    {
	for (aa = -1; aa < 20; ++aa)
	{
	    Codon codon = AAtoCodon(aa, 1);
	    AminoAcid aa2 = TripletToAA(codon);
	    if (aa != aa2)
	    {
		error_count++;
	    }
	}
    }

    printf("Error count %d\n", error_count);
   
    return 0;
}
