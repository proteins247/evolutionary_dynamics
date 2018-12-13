#include "Random123/threefry.h"
#include "Random123/u01fixedpt.h"

/* 
 * threefryrand() generates random doubles in the range [0.0, 1.0)
 * using an algorithm from the Random123 library. threefryrand_int
 * returns a uint64.
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

#ifndef RNG_H_
#define RNG_H_

/* Sets RNG "seed" (state) */
void set_threefry_array(uint64_t uk0, uint64_t uk1,
                        uint64_t uk2, uint64_t uk3);

/* Returns RNG "seed" */
/* Modifies current_key; it had better be >= 4 long  */
void get_rng_state(uint64_t * current_key);

/* Return double in range [0., 1.) */
double threefryrand();

/* Return random unsigned integer of 64 bits */
uint64_t threefryrand_int();

/* diagnostic */
void printf_counter();

#endif /* RNG_H_ */
