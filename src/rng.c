#include "rng.h"
#include <stdio.h>

static threefry4x64_ctr_t ctr = {{}};
static threefry4x64_key_t key = {{}};

static unsigned char randomNumberIndex = 4;
static threefry4x64_ctr_t result = {{}};


void increment_counter()
{
    static unsigned char index = 0;
    if (index == 4) {
        index = 0;
    }
    ctr.v[index++]++;
}


/* Functions implementations for rng.h */

void set_threefry_array(uint64_t uk0, uint64_t uk1,
                        uint64_t uk2, uint64_t uk3)
{
    key.v[0] = uk0;
    key.v[1] = uk1;
    key.v[2] = uk2;
    key.v[3] = uk3;
}

void set_threefry_counter(uint64_t uc0, uint64_t uc1,
                          uint64_t uc2, uint64_t uc3)
{
    ctr.v[0] = uc0;
    ctr.v[1] = uc1;
    ctr.v[2] = uc2;
    ctr.v[3] = uc3;
}

void set_threefry_result(uint64_t result0, uint64_t result1,
			 uint64_t result2, uint64_t result3,
			 unsigned char index)
{
    result.v[0] = result0;
    result.v[1] = result1;
    result.v[2] = result2;
    result.v[3] = result3;
    randomNumberIndex = index;
}


void get_rng_state(uint64_t * current_key,
                   uint64_t * current_counter,
                   uint64_t * current_result,
                   unsigned char * current_index)
{
    current_key[0] = key.v[0];
    current_key[1] = key.v[1];
    current_key[2] = key.v[2];
    current_key[3] = key.v[3];
    current_counter[0] = ctr.v[0];
    current_counter[1] = ctr.v[1];
    current_counter[2] = ctr.v[2];
    current_counter[3] = ctr.v[3];
    current_result[0] = result.v[0];
    current_result[1] = result.v[1];
    current_result[2] = result.v[2];
    current_result[3] = result.v[3];
    *current_index = randomNumberIndex;    
}

double threefryrand()
{
    if (randomNumberIndex == 4)
    {
        result = threefry4x64(ctr, key);
        randomNumberIndex = 0;
        increment_counter();
    }
    return u01fixedpt_closed_open_64_53(result.v[randomNumberIndex++]);
}

uint64_t threefryrand_int()
{
   if (randomNumberIndex == 4)
   {
        result = threefry4x64(ctr, key);
        randomNumberIndex = 0;
        increment_counter();
    }
    return result.v[randomNumberIndex++];
}

void printf_counter(const char * prefix)
{
    if (prefix)
    {
        printf("%s: %lu %lu %lu %lu %u\n", prefix,
               ctr.v[0], ctr.v[1], ctr.v[2], ctr.v[3], randomNumberIndex);
    }
    else
    {
        printf("%lu %lu %lu %lu %u\n",
               ctr.v[0], ctr.v[1], ctr.v[2], ctr.v[3], randomNumberIndex);
    }
}
