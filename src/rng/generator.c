/*	generates a better pseudo-random number 
 	than C's built in function rand()
	from D Jones 2010 JKISS generator webpage 
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>


#include "generator.h"

static unsigned int x, y, z, c; /* Seed variables */

unsigned int devrand(void) 
{
	unsigned int r;
	FILE *fn;

	fn = fopen("/dev/urandom", "r"); 
	fread(&r, 1, 4, fn);
	fclose(fn);

	return r;
}

/* Initialise KISS generator using /dev/urandom */ 
void init_KISS() 
{
	x = devrand(); 
	while (!(y = devrand())); /* y must not be zero! */ 
	z = devrand();

/* We don’t really need to set c as well but let's anyway… */ 
/* NOTE: offset c by 1 to avoid z=c=0 */ 
	c = devrand() % 698769068 + 1; /* Should be less than 698769069 */
	fprintf(stderr, "Initalized seeds: x=%d y=%d z=%d c=%d\n", x, y, z, c);
}

unsigned int JKISS()
{
	unsigned long long t;
	x = 314527869 * x + 1234567; 
	y ^= y << 5; y ^= y >> 7; y ^= y << 22; 
	t = 4294584393ULL * z + c; c = t >> 32; z = t;

	return x + y + z;
}
