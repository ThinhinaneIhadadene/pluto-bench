#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <sys/time.h>
#include <stdint.h>

#define ROWS 2000
#define COLS 2000
#define BASE 255
#pragma declarations
uint8_t src[ROWS][COLS][3];
float kernelX[5];
float kernelY[5];
uint8_t temp[ROWS][COLS][3];
uint8_t conv[ROWS][COLS][3];
#pragma enddeclarations

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
                src[i][j][k]=rand()%BASE;
            }
        }
    }
}
void init_kernels(){
      kernelX[0] = 1.0f; kernelX[1] = 4.0f; kernelX[2] = 6.0f; kernelX[3] = 4.0f; kernelX[4] = 1.0f;
      kernelY[0] = 1.0f/256; kernelY[1] = 4.0f/256; kernelY[2] = 6.0f/256; kernelY[3] = 4.0f/256; kernelY[4] = 1.0f/256;
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
	int q,w,cc,e,r;
  double t_start, t_end;
	init_array() ;

	IF_TIME(t_start = rtclock());

#pragma scop
for (q = 0; q < ROWS - 5; q++) {
    for ( w = 0; w < COLS - 5; w++) {
        for (cc = 0; cc < 3; cc++) {
            temp[q][w][cc] = src[q][w][cc] * kernelX[0]+src[q + 1][w][cc] * kernelX[1]+
            src[q + 2][w][cc] * kernelX[2]+src[q + 3][w][cc] * kernelX[3]+
            src[q + 4][w][cc] * kernelX[4];
        }
    }
}
for (q = 0; q < ROWS - 5; q++) {
    for ( w = 0; w < COLS - 5; w++) {
        for (cc = 0; cc < 3; cc++) {
            conv[q][w][cc] = temp[q][w][cc] * kernelY[0]+temp[q][w + 1][cc] * kernelY[1]+
            temp[q][w + 2][cc] * kernelY[2]+temp[q][w + 3][cc] * kernelY[3]+
            temp[q][w + 4][cc] * kernelY[4];
        }
    }
}
#pragma endscop

	IF_TIME(t_end = rtclock());
	IF_TIME(fprintf(stdout, "%0.6lfs\n", t_end - t_start));

    if (fopen(".test", "r")) {
#ifdef __MPI
        if (my_rank == 0) {
            print_array();
        }
#else
        print_array();
#endif
    }
    return 0;
}
