#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <stdint.h>
#include <math.h>
#include "decls.h"

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

void print_matrix(uint8_t a[N][M][3]){
  int i,j,k;
  for(i=0;i<N;i++)
   for(j=0;j<M;j++)
     {
       for(k=0;k<3;k++)
         printf("%d,", a[i][j][k]);
        printf("\t");
     }printf("\n");
}
void init_matrix(){
  int i,j,c;
  /* Initialization */
  srand(0); // seed with a constant value to verify results
  for(i=0;i<N;i++)
   for(j=0;j<M;j++)
     for(c=0;c<3;c++)
       src[i][j][c]= rand()%BASE;

}

int main(int argc, char * argv[]) {
  long int i,j,c;
  double t_start = 0.0, t_end = 0.0;
  double Total = 0.0;
  init_matrix();

#ifdef DEBUG
   printf("Cols = %ld Rows= %ld \n", M,N);
#endif

IF_TIME(t_start = rtclock());

#pragma scop
for (i = 0; i < N - 2; i++) {
    for (j = 0; j < M - 2; j++) {
        for (c = 0; c < 3; c++) {
            temp[i][j][c] = (src[i][j][c]   + src[i][j+1][c]   + src[i][j+2][c]+
                             src[i+1][j][c]                    + src[i+1][j+2][c]+
                             src[i+2][j][c] + src[i+2][j+1][c] + src[i+2][j+2][c])/((unsigned char) 8);
        }
    }
}

for (i = 0; i < N- 2; i++) {
    for (j = 0; j < M - 2; j++) {
        for (c = 0; c < 3; c++) {
            out[i][j][c] = (temp[i+1][j+1][c]-temp[i+2][j][c]) + (temp[i+2][j+1][c]-temp[i+1][j][c]);
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
      //
    }
#else
//
#endif
  }
  return 0;
}
