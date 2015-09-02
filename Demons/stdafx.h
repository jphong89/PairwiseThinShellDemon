// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once

#include "targetver.h"

#include <stdio.h>
#include <tchar.h>

#define VAL(SA,ia,ja,na) SA[(ia)*(na)+(ja)]

#define PI 3.141592654
#define EPS (1e-6)
#define max_zenyo(a,b)            (((a) > (b)) ? (a) : (b))
#define min_zenyo(a,b)            (((a) < (b)) ? (a) : (b))

#define EUCLWEIGHT 1
#define DIR_NUM 8

#define REPELTHRESHOLD 0.5 // used for avoiding collision
#define FEATDIM (59 + 27) //feature dimension
#define GEODECISDIS 4.0 //length of the "spider leg"
#define DISTHRESHOLD 6.0 // we don't want to match things too far away
#define SIGMOIDSIGMA 1
#define SIGMOIDALPHA 50
#define SIMILARITY_SCALE 4

#define STRETCH_MIU (YOUNG/2/(1+POISSON)) 
#define STRETCH_LAMBDA (YOUNG*POISSON/(1+POISSON)/(1-2*POISSON))
#define YOUNG 2
#define POISSON 0.05
#define REGWEIGHT 20 //stretching weight
#define BENDWEIGHT 10  //bending weight
#define LINKWEIGHT 200  //structural links weight

#define GPU_ON 0

// TODO: reference additional headers your program requires here
