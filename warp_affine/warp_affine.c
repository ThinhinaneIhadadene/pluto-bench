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
#define mixf(x, y, a) ((x) + ((x) - (y)) * (a))
#define clamp(x, minval, maxval) ((x) < (minval) ? (minval) : (x) > (maxval) ? (maxval) : (x))
#define SRC_ROWS 5000
#define SRC_COLS 5000
#define DST_ROWS 5000
#define DST_COLS 5000
#define BASE 255
#endif
/* Define our arrays */
#pragma declarations
uint8_t src[SRC_ROWS][SRC_COLS];
float dst[DST_ROWS][DST_COLS];
const float a00=0.1f;
const float a01=0.1f;
const float a10=0.1f;
const float a11=0.1f;
const float b00=0.1f;
const float b10 =0.1f;
#pragma enddeclarations

#define __PLACE_TO_INSERT_FORWARD_DECLARATIONS

int main(int argc, char * argv[]) {
  long int i,j;
  double t_start = 0.0, t_end = 0.0;
  double Total = 0.0;

#ifdef DEBUG
   printf("Cols = %ld Rows= %ld \n", COLS,ROWS);
#endif

  /* Initialization */
  srand(0); // seed with a constant value to verify results
  for(i=0;i<SRC_ROWS;i++)
   for(j=0;j<SRC_COLS;j++)
       src[i][j]= rand()%BASE;

  IF_TIME(t_start = rtclock());
  //all declarations must go here, else u get errors which don't indicate this
  float o_c,o_r,r,c;
  int coord_00_c,coord_00_r,coord_01_c,coord_01_r,coord_10_c,coord_10_r,coord_11_c,coord_11_r;
  float A00,A01,A10,A11;
#pragma scop
for (int n_r = 0; n_r < DST_ROWS; n_r++) {
    for (int n_c = 0; n_c < DST_COLS; n_c++) {
        o_r = a11 * n_r + a10 * n_c + b00;
        o_c = a01 * n_r + a00 * n_c + b10;

        r = o_r - floorf(o_r);
        c = o_c - floorf(o_c);

        coord_00_r = floorf(o_r);
        coord_00_c = floorf(o_c);
        coord_01_r = coord_00_r;
        coord_01_c = coord_00_c + 1;
        coord_10_r = coord_00_r + 1;
        coord_10_c = coord_00_c;
        coord_11_r = coord_00_r + 1;
        coord_11_c = coord_00_c + 1;

        coord_00_r = clamp(coord_00_r, 0, SRC_ROWS - 1);
        coord_00_c = clamp(coord_00_c, 0, SRC_COLS - 1);
        coord_01_r = clamp(coord_01_r, 0, SRC_ROWS - 1);
        coord_01_c = clamp(coord_01_c, 0, SRC_COLS - 1);
        coord_10_r = clamp(coord_10_r, 0, SRC_ROWS - 1);
        coord_10_c = clamp(coord_10_c, 0, SRC_COLS - 1);
        coord_11_r = clamp(coord_11_r, 0, SRC_ROWS - 1);
        coord_11_c = clamp(coord_11_c, 0, SRC_COLS - 1);

        A00 = src[coord_00_r][coord_00_c];
        A10 = src[coord_10_r][coord_10_c];
        A01 = src[coord_01_r][coord_01_c];
        A11 = src[coord_11_r][coord_11_c];

        dst[n_r][n_c] = mixf(mixf(A00, A10, r), mixf(A01, A11, r), c);
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

// icc -O3 -fp-model precise heat_1d_np.c -o op-heat-1d-np -lm
