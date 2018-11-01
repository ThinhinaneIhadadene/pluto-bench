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

#define CV_DESCALE(x, n) (((x) + (1 << ((n)-1))) >> (n))
enum {
    yuv_shift  = 14,
    R2Y        = 4899,
    G2Y        = 9617,
    B2Y        = 1868,
};

/* Define our arrays */
#pragma declarations
uint8_t src[ROWS][COLS][3];
uint8_t dst[ROWS][COLS];
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
  long int i,j,c;
  double t_start = 0.0, t_end = 0.0;
  double Total = 0.0;

#ifdef DEBUG
   printf("Cols = %ld Rows= %ld \n", COLS,ROWS);
#endif

  /* Initialization */
  srand(0); // seed with a constant value to verify results
  for(i=0;i<ROWS;i++)
   for(j=0;j<COLS;j++)
     for(c=0;c<3;c++)
       src[i][j][c]= rand()%BASE;

  IF_TIME(t_start = rtclock());

#pragma scop

for (int q = 0; q < ROWS; q++) {
    for (int w = 0; w < COLS; w++) {
        dst[q][w] = CV_DESCALE((src[q][w][2] * B2Y + src[q][w][1] * G2Y + src[q][w][0] * R2Y), yuv_shift);
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
