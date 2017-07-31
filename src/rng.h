#include "Random123/threefry.h"
#include "Random123/u01fixedpt.h"

/* 
 * threefryrand() generates random doubles in the range [0.0, 1.0]
 * using an algorithm from the Random123 library.
 *
 * Initialize by calling set_threefry_array with a user-supplied
 * seed (2^64 - 1 possible seeds). 
 *
 * For multi-processing situations, the user-supplied seed should
 * be different for each process if one seeks to avoid having
 * the same stream of random numbers for each thread/process.
 *
 * I'm unsure about the behavior of these functions in a
 * multi-threaded situation where there is only one instance of these
 * variables. Need to learn about thread_local (C11)
 *
 */
static threefry4x64_ctr_t ctr = {{}};
static threefry4x64_key_t key = {{}};

void increment_counter() {
    static unsigned char index = 0;
    if (index == 4) {
        index = 0;
    }
    ctr.v[index++]++;
}

void set_threefry_array(unsigned long int user_key) {
    key.v[0] = user_key;
    key.v[1] = user_key;
    key.v[2] = user_key;
    key.v[3] = user_key;
}

double threefryrand() {
    static unsigned char randomNumberIndex = 4;
    static threefry4x64_ctr_t result = {{}};
    if (randomNumberIndex == 4) {
        result = threefry4x64(ctr, key);
        randomNumberIndex = 0;
        increment_counter();
    }
    return u01fixedpt_closed_open_64_53(result.v[randomNumberIndex++]);
}

