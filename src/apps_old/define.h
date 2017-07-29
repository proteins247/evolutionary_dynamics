/*************************
 Modified by Orit Peleg
 opeleg@fas.harvard.edu
 on 2/26/13
 ************************/


/*
 * define.h
 *
 *  Created on: Dec 25, 2009
 *      Author: mheo
 */

#ifndef DEFINE_H_
#define DEFINE_H_

#define AASEQLEN        27
#define AASURFACELEN    9
#define NUCSEQLEN       81

#define NUMCONF         10000
// redefine to 103346 for the full set

#define HALFCONTMATRLEN 28
#define ALPHABET        20

#define MCATTEMPT       1000

#define INITIAL_POPSIZE 100


// !!! should be as big as POPSIZE:
#define MAXORGANISMS    100

#define MAXGENES        7
// maximum genes per organism
#define MAXPPIS			6
// maximum interactions among genes
#define MAXSTATES       14
#define NUFOLD          10
// number of unfolded conformations


//#define MAXGENES        3
//// maximum genes per organism
//#define MAXPPIS			2
//// maximum interactions among genes
//#define MAXSTATES       6
//#define NUFOLD          10
//// number of unfolded conformations

//#define MAXPRINTSPECIES  100  // maximum species for seqlog
#define POPSIZE         100

#define ALGORITHM       1
// 0 muyoung's. synchronized devisions, flex pop size
// 1 gillepsie (murat). fixed pop size

#define BEST_FACE 0
// for singlish interaction. set to a negative number to let code find best face

#define RANDOME_MODE 0
// 0 - non random (best)
// 1 - random




#define FILE_CONFIG     "_config.txt"

#endif /* DEFINE_H_ */
