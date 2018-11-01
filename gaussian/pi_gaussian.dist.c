#include <math.h>
#include <stdio.h>
#include <assert.h>
#include "pi_defs.h"
#include "polyrt.h"
#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define max(x,y)    ((x) > (y)? (x) : (y))
#define min(x,y)    ((x) < (y)? (x) : (y))

#ifndef __BLOCK_CYCLIC_BLOCK_SIZE
#define __BLOCK_CYCLIC_BLOCK_SIZE 1
#endif

int pi_0( int t1, int t2, int t4, int nprocs)
{
                int __lb0, __ub0, __p0;
                __lb0 = _LB_REPLACE_ME_DISTLOOG0t4;
                __ub0 = _UB_REPLACE_ME_DISTLOOG0t4;
#ifdef __AUTO_COMPUTE_PI
                    __p0 = polyrt_one_dim_pi(t4, __lb0, __ub0, nprocs);
                    return pi_mappings_0[__p0];
#endif
#ifdef __USE_BLOCK_CYCLIC
                    if (__is_block_cyclic[0])
                     return polyrt_one_dim_pi_block_cyclic(t4, __lb0, __ub0, nprocs, __BLOCK_CYCLIC_BLOCK_SIZE);
                    else
 #endif
                    return polyrt_one_dim_pi(t4, __lb0, __ub0, nprocs);
}
