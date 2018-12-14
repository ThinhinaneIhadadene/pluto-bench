#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <sys/time.h>
#include <stdint.h>
#include <omp.h>

//constants
#define ROWS 2000
#define COLS 2000
#define BASE 255

//global vars
#pragma declarations
uint8_t src[ROWS][COLS][3];
float kernel[3][3];
uint8_t conv[ROWS][COLS][3];
#pragma enddeclarations

#include <unistd.h>
#include <sys/time.h>
#include <math.h>

#ifdef TIME
#define IF_TIME(foo) foo;
#else
#define IF_TIME(foo)
#endif

void init_array()
{
    int i, j, k;
    srand(0);
    for (i=0; i<ROWS; i++) {
        for (j=0; j<COLS; j++) {
            for (k=0; k<3; k++) {
                src[i][j][k] = rand()%BASE;
            }
        }
    }
}

void init_kernel(){
  kernel[0][0]=0;kernel[1][0]=1.0f/5;kernel[2][0]=0;
  kernel[0][1]=1.0f/5;kernel[1][1]=1.0f/5;kernel[2][1]=0;
  kernel[0][2]=0;kernel[1][2]=1;kernel[2][2]=0;
}


void print_array()
{
  int i,j,k;
  for (i=0; i<ROWS; i++) {
      for (j=0; j<COLS; j++) {
          for (k=0; k<3; k++) {
              printf("%d,", conv[i][j][k]);
          }
          printf("\t");
      }
      printf("\n");
  }
}

double rtclock()
{
    struct timezone Tzp;
    struct timeval Tp;
    int stat;
    stat = gettimeofday (&Tp, &Tzp);
    if (stat != 0) printf("Error return from gettimeofday: %d",stat);
    return(Tp.tv_sec + Tp.tv_usec*1.0e-6);
}

int main()
{
	 int q,w,cc,kq,kw;
   double t_start, t_end;
   float prod;

	init_array() ;
  init_kernel();
	IF_TIME(t_start = rtclock());

#pragma scop
for (q = 0; q < ROWS- 2; q++) {
    for (w = 0; w < COLS - 2; w++) {
        for (cc = 0; cc < 3; cc++) {
          conv[q][w][cc]=0;
          for ( kq = 0; kq < 3; kq++) {
              for ( kw = 0; kw < 3; kw++) {
                  conv[q][w][cc]+= src[q + kq][w + kw][cc] * kernel[kq][kw];
              }
          }
        }
    }
}
#pragma endscop

	IF_TIME(t_end = rtclock());
	IF_TIME(fprintf(stdout, "%0.6lfs\n", t_end - t_start));

    if (fopen(".test", "r")) {
#ifdef __MPI
        if (my_rank == 0) {
          //  print_array();
        }
#else
        //print_array();
#endif
    }
    return 0;
}
