#include <omp.h>
#include <math.h>
#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define max(x,y)    ((x) > (y)? (x) : (y))
#define min(x,y)    ((x) < (y)? (x) : (y))

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <stdint.h>
#include <math.h>

#ifdef PERFCTR
#include "papiStdEventDefs.h"
#include <papi.h>
#include "papi_defs.h"
#endif

#ifdef TIME
#define IF_TIME(foo) foo;
#else
#define IF_TIME(foo)
#endif

double rtclock()
{
    struct timezone Tzp;
    struct timeval Tp;
    int stat;
    stat = gettimeofday (&Tp, &Tzp);
    if (stat != 0) printf("Error return from gettimeofday: %d",stat);
    return(Tp.tv_sec + Tp.tv_usec*1.0e-6);
}

#ifdef HAS_DECLS
#include "decls.h"
#else
#define ROWS 5000
#define COLS 5000
#define BASE 255
#endif
/* Define our arrays */
#pragma declarations
uint8_t src[ROWS][COLS][3];
uint8_t conv[ROWS][COLS][3];
uint8_t temp[ROWS][COLS][3];
float kernelX[5];
float kernelY[5];
float prod1;
float prod2;
#pragma enddeclarations

#define __PLACE_TO_INSERT_FORWARD_DECLARATIONS

void print_matrix(uint8_t a[ROWS][COLS][3]){
  int i,j,k;
  for(i=0;i<ROWS;i++)
   for(j=0;j<COLS;j++)
     {
       for(k=0;k<3;k++)
         printf("%d,", a[i][j][k]);
        printf("\t");
     }printf("\n");
}

int main(int argc, char * argv[]) {
  long int q,w,cc,r,e;
  double t_start = 0.0, t_end = 0.0;
  double Total = 0.0;

#ifdef DEBUG
   printf("Cols = %ld Rows= %ld \n", COLS,ROWS);
#endif

srand(0);
for (q = 0; q < ROWS - 5; q++) {
    for (w = 0; w < COLS - 5; w++) {
        for (cc = 0; cc < 3; cc++) {
          src[q][w][cc]= rand()%BASE;
        }
      }
    }
//declare kernels
kernelX[0] = 1.0f; kernelX[1] = 4.0f; kernelX[2] = 6.0f; kernelX[3] = 4.0f; kernelX[4] = 1.0f;
kernelY[0] = 1.0f/256; kernelY[1] = 4.0f/256; kernelY[2] = 6.0f/256; kernelY[3] = 4.0f/256; kernelY[4] = 1.0f/256;

  IF_TIME(t_start = rtclock());

/* Copyright (C) 1991-2014 Free Software Foundation, Inc.
   This file is part of the GNU C Library.

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the GNU C Library; if not, see
   <http://www.gnu.org/licenses/>.  */
/* This header is separate from features.h so that the compiler can
   include it implicitly at the start of every compilation.  It must
   not itself include <features.h> or any other header that includes
   <features.h> because the implicit include comes before any feature
   test macros that may be defined in a source file before it first
   explicitly includes a system header.  GCC knows the name of this
   header in order to preinclude it.  */
/* glibc's intent is to support the IEC 559 math functionality, real
   and complex.  If the GCC (4.9 and later) predefined macros
   specifying compiler intent are available, use them to determine
   whether the overall intent is to support these features; otherwise,
   presume an older compiler has intent to support these features and
   define these macros by default.  */
/* wchar_t uses ISO/IEC 10646 (2nd ed., published 2011-03-15) /
   Unicode 6.0.  */
/* We do not support C11 <threads.h>.  */
  int t1, t2, t3, t4, t5, t6, t7, t8;
  int lb, ub, lbd, ubd, lb2, ub2;
  register int lbv, ubv;
/* Start of CLooG code */
for (t1=0;t1<=4994;t1++) {
  for (t2=0;t2<=5037;t2++) {
    t3 = floord(t2+1,126);
    for (t4=ceild(6*t2-6*t3-124,125);t4<=min(floord(2*t2+498*t3+500,125),floord(6*t2-6*t3+6,125));t4++) {
      if ((t2 == 21*t4-1) && (t2 >= 126*t3+41)) {
        if ((t2+1)%42 == 0) {
          prod1 += (src[ t1 + 4][ (t2-t3)][ ((-t2+126*t3+83)/42)] * kernelX[ 4]);;
          /* unknown - failure constructing stmt body */;
          prod2 += (temp[ t1][ (t2-t3-4) + 4][ ((-t2+126*t3+83)/42)] * kernelY[ 4]);;
        }
        if ((t2+22)%42 == 0) {
          for (t7=ceild(250*t2+250,21);t7<=floord(250*t2+271,21);t7++) {
            prod1 += (src[ t1 + ((-250*t2+21*t7-208)/21)][ (t2-t3)][ ((-t2+126*t3+104)/42)] * kernelX[ ((-250*t2+21*t7-208)/21)]);;
            prod2 += (temp[ t1][ (t2-t3-4) + ((-250*t2+21*t7-208)/21)][ ((-t2+126*t3+104)/42)] * kernelY[ ((-250*t2+21*t7-208)/21)]);;
          }
          prod1 += (src[ t1 + 4][ (t2-t3)][ ((-t2+126*t3+104)/42)] * kernelX[ 4]);;
          /* unknown - failure constructing stmt body */;
          prod2 += (temp[ t1][ (t2-t3-4) + 4][ ((-t2+126*t3+104)/42)] * kernelY[ 4]);;
        }
        if ((41*t2+20)%42 <= 21) {
          t6 = floord(83*t2+42*t3+104,42);
          if ((t2+1)%21 == 0) {
            /* unknown - failure constructing stmt body */;
          }
        }
        for (t6=ceild(83*t2+42*t3+125,42);t6<=2*t2-2*t3+2;t6++) {
          if ((t2+1)%21 == 0) {
            /* unknown - failure constructing stmt body */;
          }
          if ((t2+1)%21 == 0) {
            /* unknown - failure constructing stmt body */;
          }
          for (t7=4*t2-4*t3+4*t6;t7<=4*t2-4*t3+4*t6+3;t7++) {
            if ((t2+1)%21 == 0) {
              prod1 += (src[ t1 + (-4*t2+4*t3-4*t6+t7)][ (t2-t3)][ (-2*t2+2*t3+t6)] * kernelX[ (-4*t2+4*t3-4*t6+t7)]);;
            }
            if ((t2+1)%21 == 0) {
              prod2 += (temp[ t1][ (t2-t3-4) + (-4*t2+4*t3-4*t6+t7)][ (-2*t2+2*t3+t6)] * kernelY[ (-4*t2+4*t3-4*t6+t7)]);;
            }
          }
          if ((t2+1)%21 == 0) {
            prod1 += (src[ t1 + 4][ (t2-t3)][ (-2*t2+2*t3+t6)] * kernelX[ 4]);;
          }
          if ((t2+1)%21 == 0) {
            /* unknown - failure constructing stmt body */;
          }
          if ((t2+1)%21 == 0) {
            prod2 += (temp[ t1][ (t2-t3-4) + 4][ (-2*t2+2*t3+t6)] * kernelY[ 4]);;
          }
          if ((t2+1)%21 == 0) {
            /* unknown - failure constructing stmt body */;
          }
        }
      }
      if ((t2 == 21*t4+20) && (t2 >= 126*t3+20) && (t2 <= 126*t3+62)) {
        for (t6=2*t2-2*t3;t6<=floord(83*t2+42*t3+62,42);t6++) {
          if ((t2+1)%21 == 0) {
            /* unknown - failure constructing stmt body */;
          }
          if ((t2+1)%21 == 0) {
            /* unknown - failure constructing stmt body */;
          }
          for (t7=4*t2-4*t3+4*t6;t7<=4*t2-4*t3+4*t6+3;t7++) {
            if ((t2+1)%21 == 0) {
              prod1 += (src[ t1 + (-4*t2+4*t3-4*t6+t7)][ (t2-t3)][ (-2*t2+2*t3+t6)] * kernelX[ (-4*t2+4*t3-4*t6+t7)]);;
            }
            if ((t2+1)%21 == 0) {
              prod2 += (temp[ t1][ (t2-t3-4) + (-4*t2+4*t3-4*t6+t7)][ (-2*t2+2*t3+t6)] * kernelY[ (-4*t2+4*t3-4*t6+t7)]);;
            }
          }
          if ((t2+1)%21 == 0) {
            prod1 += (src[ t1 + 4][ (t2-t3)][ (-2*t2+2*t3+t6)] * kernelX[ 4]);;
          }
          if ((t2+1)%21 == 0) {
            /* unknown - failure constructing stmt body */;
          }
          if ((t2+1)%21 == 0) {
            prod2 += (temp[ t1][ (t2-t3-4) + 4][ (-2*t2+2*t3+t6)] * kernelY[ 4]);;
          }
          if ((t2+1)%21 == 0) {
            /* unknown - failure constructing stmt body */;
          }
        }
        if ((41*t2+20)%42 <= 21) {
          t6 = floord(83*t2+42*t3+104,42);
          if ((t2+1)%21 == 0) {
            /* unknown - failure constructing stmt body */;
          }
          if ((t2+1)%21 == 0) {
            /* unknown - failure constructing stmt body */;
          }
          for (t7=4*t2-4*t3+4*t6;t7<=floord(250*t2+229,21);t7++) {
            if ((t2+1)%21 == 0) {
              prod1 += (src[ t1 + (-4*t2+4*t3-4*t6+t7)][ (t2-t3)][ (-2*t2+2*t3+t6)] * kernelX[ (-4*t2+4*t3-4*t6+t7)]);;
            }
            if ((t2+1)%21 == 0) {
              prod2 += (temp[ t1][ (t2-t3-4) + (-4*t2+4*t3-4*t6+t7)][ (-2*t2+2*t3+t6)] * kernelY[ (-4*t2+4*t3-4*t6+t7)]);;
            }
          }
        }
      }
      if ((t2 <= min(t3+4994,21*t4+19)) && (t2 >= max(21*t4,t3+4))) {
        for (t6=2*t2-2*t3;t6<=min(250*t3+249,2*t2-2*t3+2);t6++) {
          /* unknown - failure constructing stmt body */;
          /* unknown - failure constructing stmt body */;
          for (t7=4*t2-4*t3+4*t6;t7<=4*t2-4*t3+4*t6+3;t7++) {
            prod1 += (src[ t1 + (-4*t2+4*t3-4*t6+t7)][ (t2-t3)][ (-2*t2+2*t3+t6)] * kernelX[ (-4*t2+4*t3-4*t6+t7)]);;
            prod2 += (temp[ t1][ (t2-t3-4) + (-4*t2+4*t3-4*t6+t7)][ (-2*t2+2*t3+t6)] * kernelY[ (-4*t2+4*t3-4*t6+t7)]);;
          }
          prod1 += (src[ t1 + 4][ (t2-t3)][ (-2*t2+2*t3+t6)] * kernelX[ 4]);;
          /* unknown - failure constructing stmt body */;
          prod2 += (temp[ t1][ (t2-t3-4) + 4][ (-2*t2+2*t3+t6)] * kernelY[ 4]);;
          /* unknown - failure constructing stmt body */;
        }
      }
      if ((t2 <= 3) && (t3 == 0) && (t4 == 0)) {
        for (t6=2*t2;t6<=2*t2+2;t6++) {
          /* unknown - failure constructing stmt body */;
          for (t7=4*t2+4*t6;t7<=4*t2+4*t6+4;t7++) {
            prod1 += (src[ t1 + (-4*t2-4*t6+t7)][ t2][ (-2*t2+t6)] * kernelX[ (-4*t2-4*t6+t7)]);;
          }
          /* unknown - failure constructing stmt body */;
        }
      }
      if ((t2 == 126*t3-1) && (t2 == 21*t4+20)) {
        if ((125*t2+125)%126 == 0) {
          /* unknown - failure constructing stmt body */;
          /* unknown - failure constructing stmt body */;
          for (t7=ceild(250*t2+166,21);t7<=floord(250*t2+229,21);t7++) {
            prod1 += (src[ t1 + ((-250*t2+21*t7-166)/21)][ ((125*t2-1)/126)][ 2] * kernelX[ ((-250*t2+21*t7-166)/21)]);;
            prod2 += (temp[ t1][ ((125*t2-505)/126) + ((-250*t2+21*t7-166)/21)][ 2] * kernelY[ ((-250*t2+21*t7-166)/21)]);;
          }
        }
      }
      if ((t2 == 21*t4+20) && (t2 >= 126*t3+83)) {
        if ((t2+1)%21 == 0) {
          /* unknown - failure constructing stmt body */;
        }
        if ((t2+1)%21 == 0) {
          /* unknown - failure constructing stmt body */;
        }
        for (t7=12*t2-12*t3;t7<=floord(250*t2+229,21);t7++) {
          if ((t2+1)%21 == 0) {
            prod1 += (src[ t1 + (-12*t2+12*t3+t7)][ (t2-t3)][ 0] * kernelX[ (-12*t2+12*t3+t7)]);;
          }
          if ((t2+1)%21 == 0) {
            prod2 += (temp[ t1][ (t2-t3-4) + (-12*t2+12*t3+t7)][ 0] * kernelY[ (-12*t2+12*t3+t7)]);;
          }
        }
      }
      if ((t2 == 126*t3+20) && (t2 == 21*t4-1)) {
        if ((125*t2+20)%126 == 0) {
          for (t7=ceild(250*t2+250,21);t7<=floord(250*t2+271,21);t7++) {
            prod1 += (src[ t1 + ((-250*t2+21*t7-208)/21)][ ((125*t2+20)/126)][ 2] * kernelX[ ((-250*t2+21*t7-208)/21)]);;
            prod2 += (temp[ t1][ ((125*t2-484)/126) + ((-250*t2+21*t7-208)/21)][ 2] * kernelY[ ((-250*t2+21*t7-208)/21)]);;
          }
          prod1 += (src[ t1 + 4][ ((125*t2+20)/126)][ 2] * kernelX[ 4]);;
          /* unknown - failure constructing stmt body */;
          prod2 += (temp[ t1][ ((125*t2-484)/126) + 4][ 2] * kernelY[ 4]);;
        }
      }
      if ((t2 == 126*t3-1) && (t2 == 21*t4-1)) {
        if ((125*t2+125)%126 == 0) {
          prod1 += (src[ t1 + 4][ ((125*t2-1)/126)][ 2] * kernelX[ 4]);;
          /* unknown - failure constructing stmt body */;
          prod2 += (temp[ t1][ ((125*t2-505)/126) + 4][ 2] * kernelY[ 4]);;
        }
      }
      if ((t2 == 21*t4-1) && (t2 <= 126*t3+20)) {
        if ((t2+1)%21 == 0) {
          /* unknown - failure constructing stmt body */;
        }
      }
      if ((t2 >= 5034) && (t3 == 39) && (t4 == 239)) {
        for (t6=2*t2-78;t6<=2*t2-76;t6++) {
          /* unknown - failure constructing stmt body */;
          for (t7=4*t2+4*t6-156;t7<=4*t2+4*t6-152;t7++) {
            prod2 += (temp[ t1][ (t2-43) + (-4*t2-4*t6+t7+156)][ (-2*t2+t6+78)] * kernelY[ (-4*t2-4*t6+t7+156)]);;
          }
          /* unknown - failure constructing stmt body */;
        }
      }
    }
  }
}
/* End of CLooG code */

  IF_TIME(t_end = rtclock());
  IF_TIME(fprintf(stdout, "time = %0.6lfs\n", t_end - t_start));

#ifdef __MPI
  if (my_rank == 0) {
#endif

#ifdef __MPI
  }
#endif

  if (fopen(".test", "r")) {
#ifdef __MPI
    if (my_rank == 0) {
      // print_matrix(out);
    }
#else
    // print_matrix(out);
#endif
  }

  return 0;
}
