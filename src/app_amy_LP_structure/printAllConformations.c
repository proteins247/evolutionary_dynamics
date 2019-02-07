/* 
 * printAllConformations : take the 10000 conformations and print
 * their conformations as latPack strings
 */
#include <stdio.h>

#include "../define.h"
#include "../latticelib.h"
#include "general.h"
#include "structurelib.h"


int main()
{
    char * dataDir = getenv("FOLDEVO_SHARE");
    char energyName[100] = "energy";
    char buffer[100];
    ReadCommondata(dataDir, energyName);

    int i;
    for (i = 0; i < NUMCONF; ++i)
    {
	PrintFoldedConformation(i, buffer);
	printf("%s\n", buffer);
    }

    return 0;
}
