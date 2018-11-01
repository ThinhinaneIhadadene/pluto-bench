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

#pragma scop
for (q = 0; q < ROWS - 5; q++) {
    for (w = 0; w < COLS - 5; w++) {
        for (cc = 0; cc < 3; cc++) {
            prod1 = 0.0f;
            for (r = 0; r < 5; r++) {
                prod1 += src[q + r][w][cc] * kernelX[r];
            }
            temp[q][w][cc] = prod1;
        }
    }
}

for (q = 0; q < ROWS - 5; q++) {
    for (w = 0; w < COLS - 5; w++) {
        for (cc = 0; cc < 3; cc++) {
            prod2 = 0.;
            for (e = 0; e < 5; e++) {
                prod2 += temp[q][w + e][cc] * kernelY[e];
            }
            conv[q][w][cc] = prod2;
        }
    }
}
#pragma endscop

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
