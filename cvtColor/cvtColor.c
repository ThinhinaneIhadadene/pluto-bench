#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <stdint.h>
#include <math.h>
#include "decls.h"

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

void init_matrix(){
  int i,j,c;
  /* Initialization */
  srand(0);
  for(i=0;i<ROWS;i++)
   for(j=0;j<COLS;j++)
     for(c=0;c<3;c++)
       src[i][j][c]= rand()%BASE;
}

int main(int argc, char * argv[]) {
  long int q,w;
  double t_start = 0.0, t_end = 0.0;
  double Total = 0.0;
  init_matrix();
#ifdef DEBUG
   printf("Cols = %ld Rows= %ld \n", COLS,ROWS);
#endif
IF_TIME(t_start = rtclock());

#pragma scop
for (q = 0; q < ROWS; q++) {
    for (w = 0; w < COLS; w++) {
        dst[q][w] = CV_DESCALE((src[q][w][2] * B2Y + src[q][w][1] * G2Y + src[q][w][0] * R2Y), yuv_shift);
    }
}
#pragma endscop

IF_TIME(t_end = rtclock());
IF_TIME(fprintf(stdout, "time = %0.6lfs\n", t_end - t_start));

#ifdef __MPI
  if (my_rank == 0) {
      //this is the master we can if we want print the result or save the image
#endif
      //do nothing we have part of images so not important
#ifdef __MPI
  }
#endif
  if (fopen(".test", "r")) {
#ifdef __MPI
    if (my_rank == 0) {
      //this is the master we can if we want print the result or save the image
    }
#else
  //do nothing we have part of images so not important
#endif
  }

  return 0;
}
