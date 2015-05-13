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

#define REPELTHRESHOLD 0.5
#define FEATDIM (59 + 27)
#define GEODECISDIS 4.0
#define DISTHRESHOLD 8.0
#define SIGMOIDSIGMA 1
#define SIGMOIDALPHA 50

#define STRETCH_MIU (YOUNG/2/(1+POISSON))
#define STRETCH_LAMBDA (YOUNG*POISSON/(1+POISSON)/(1-2*POISSON))
#define YOUNG 2
#define POISSON 0.05
#define REGWEIGHT 40
#define BENDWEIGHT 20
#define LINKWEIGHT 200
#define BESTMATCH 500

#define GPU_ON 0

// TODO: reference additional headers your program requires here
