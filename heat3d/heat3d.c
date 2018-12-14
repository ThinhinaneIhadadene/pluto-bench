#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
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

#define GLOBAL_MY_RANK
int my_rank;

#ifdef HAS_DECLS
#include "decls.h"
#else


//dimensions
#define X 110
#define Y 120
#define Z 100
//time
#define T 20000
#endif
//global vars
#pragma declarations
float a[T][Z][Y][X];
float data[Z][Y][X];
#pragma enddeclarations

double rtclock()
{
    struct timezone Tzp;
    struct timeval Tp;
    int stat;
    stat = gettimeofday (&Tp, &Tzp);
    if (stat != 0) printf("Error return from gettimeofday: %d",stat);
    return(Tp.tv_sec + Tp.tv_usec*1.0e-6);
}

void init_array(float a[Z][Y][X]) {
  const int BASE = 255;
  int x,y,z;
  /* Initialization */
  srand(time(0)); // seed with a constant value to verify results
  for(z=0;z<Z;z++)
    for(y=0;y<Y;y++)
      for(x=0;x<X;x++)
        a[z][y][x]=rand()%BASE;
}
void kernel(float heat3d[T][Z][Y][X]){
  long int t, x,y,z;

#pragma scop
for (t = 0; t <T-1; t++) {
        for (z=1;z<Z-1;z++) {
            for (y = 1; y < Y-1; y++) {
                for (x = 1; x < X-1; x++) {
                    heat3d[t+1][z][y][x] =   0.125 * (heat3d[t][z+1][y][x] - 2.0 * heat3d[t][z][y][x] + heat3d[t][z-1][y][x])
                                        + 0.125 * (heat3d[t][z][y+1][x] - 2.0 * heat3d[t][z][y][x] + heat3d[t][z][y-1][x])
                                        + 0.125 * (heat3d[t][z][y][x-1] - 2.0 * heat3d[t][z][y][x] + heat3d[t][z][y][x+1])
                                        + heat3d[t][z][y][x];

                }
            }
        }
    }
#pragma endscop
}

int main(int argc, char * argv[]) {

/* Define our arrays */

  double t_start, t_end;
  // for timekeeping
  int ts_return = -1;
  struct timeval start, end, result;
  double tdiff = 0.0;

  init_array(data);

#ifdef DEBUG
  printf("Number of points = %ld\nNumber of timesteps = %ld\n", N, T);
#endif

#ifdef TIME
  IF_TIME(t_start = rtclock());
#endif

  //Initialization
  //copy into other dimensions
  int t,z,y,x;
for(t=0;t<T;t++)
  for(z=0;z<Z;z++)
    for(y=0;y<Y;y++)
      for(x=0;x<X;x++)
        a[t][z][y][x]=data[z][y][x];
  kernel(a);

  IF_TIME(t_end = rtclock());
  IF_TIME(fprintf(stdout, "time = %lfs\n", t_end - t_start));

#ifdef PERFCTR
  PERF_EXIT;
#endif

  if (fopen(".test", "r")) {
#ifdef __MPI
    if (my_rank == 0) {

    }
#else

#endif
  }
  return 0;
}
