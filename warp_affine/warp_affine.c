#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <sys/time.h>
#include <stdint.h>

#define mixf(x, y, a) ((x) + ((x) - (y)) * (a))
#define clamp(x, minval, maxval) ((x) < (minval) ? (minval) : (x) > (maxval) ? (maxval) : (x))
#define SRC_ROWS 5000
#define SRC_COLS 5000
#define DST_ROWS 5000
#define DST_COLS 5000
#define BASE 255

#pragma declarations
const float a00=0.1f;
const float a01=0.1f;
const float a10=0.1f;
const float a11=0.1f;
const float b00=0.1f;
const float b10 =0.1f;
uint8_t src[SRC_ROWS][SRC_COLS];
float dst[DST_ROWS][DST_COLS];
float A00,A01,A10,A11;

int coord_00_c,coord_00_r,coord_01_c,coord_01_r,coord_10_c,coord_10_r,coord_11_c,coord_11_r;
float o_r,o_c,r,c;
#pragma enddeclarations

#ifdef TIME
#define IF_TIME(foo) foo;
#else
#define IF_TIME(foo)
#endif

void init_array()
{
    int i, j;
    srand(0);
    for (i=0; i<SRC_ROWS; i++) {
        for (j=0; j<SRC_COLS; j++) {
          src[i][j]=rand()%BASE;
        }
    }
}


void print_array()
{
    int i, j;

    for (i=0; i<DST_ROWS; i++) {
        for (j=0; j<DST_COLS; j++) {
          printf("%f\t",dst[i][j]);
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
  double t_start, t_end;

	init_array() ;
	IF_TIME(t_start = rtclock());
  int n_r,n_c;

#pragma scop

    for (n_r = 0; n_r < DST_ROWS; n_r++) {
        for (n_c = 0; n_c < DST_COLS; n_c++) {
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
