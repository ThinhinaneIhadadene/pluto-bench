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
float kernel[3][3];
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
  //kernel
  kernel[0][0]=0;kernel[1][0]=1.0f/5;kernel[2][0]=0;
  kernel[0][1]=1.0f/5;kernel[1][1]=1.0f/5;kernel[2][1]=1.0f/5;
  kernel[0][2]=0;kernel[1][2]=1;kernel[2][2]=0;

  IF_TIME(t_start = rtclock());
  float prod;

#pragma scop
for (int q = 0; q < ROWS - 2; q++) {
    for (int w = 0; w < COLS - 2; w++) {
        for (int cc = 0; cc < 3; cc++) {
            prod = 0.0f;
            for (int kq = 0; kq < 3; kq++) {
                for (int kw = 0; kw < 3; kw++) {
                    prod += src[q + kq][w + kw][cc] * kernel[kq][kw];
                }
            }
            conv[q][w][cc] = prod;
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

// icc -O3 -fp-model precise heat_1d_np.c -o op-heat-1d-np -lm
