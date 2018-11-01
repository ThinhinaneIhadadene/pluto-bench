#include "gaussian.dist.h"
#include <string.h>
#include "polyrt.h"
#include "pi_defs.h"
#include <assert.h>
#include <omp.h>
#include <mpi.h>
#include <limits.h>
#include <math.h>
#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define max(x,y)    ((x) > (y)? (x) : (y))
#define min(x,y)    ((x) < (y)? (x) : (y))

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

/* Copyright (C) 1991-2014 Free Software Foundation, Inc.
   This file is part of the GNU C Library.

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the GNU C Library; if not, see
   <http://www.gnu.org/licenses/>.  */
/* This header is separate from features.h so that the compiler can
   include it implicitly at the start of every compilation.  It must
   not itself include <features.h> or any other header that includes
   <features.h> because the implicit include comes before any feature
   test macros that may be defined in a source file before it first
   explicitly includes a system header.  GCC knows the name of this
   header in order to preinclude it.  */
/* glibc's intent is to support the IEC 559 math functionality, real
   and complex.  If the GCC (4.9 and later) predefined macros
   specifying compiler intent are available, use them to determine
   whether the overall intent is to support these features; otherwise,
   presume an older compiler has intent to support these features and
   define these macros by default.  */
/* wchar_t uses ISO/IEC 10646 (2nd ed., published 2011-03-15) /
   Unicode 6.0.  */
/* We do not support C11 <threads.h>.  */
  int t1, t2, t3, t4, t5, t6, t7, t8, t9;
  int lb, ub, lbd, ubd, lb2, ub2;
  register int lbv, ubv;
int _num_threads;
#pragma omp parallel
{_num_threads = omp_get_num_threads();}
#ifndef GLOBAL_MY_RANK
 int my_rank;
#endif
#ifndef __MAX_NUM_RECVS
 #define __MAX_NUM_RECVS (nprocs*_num_threads*16)
#endif
#ifndef __MAX_NUM_SENDS
 #define __MAX_NUM_SENDS (nprocs*_num_threads*8)
#endif
 int nprocs, my_start, my_end, _i, __p, proc, recv_proc, distinct_recv, __tid;
 long int _num_tasks_to_execute, _num_tasks_to_unpack;
 int count;
int req_count;
 int _lb_dist, _ub_dist, lbd_t1, ubd_t1, lbd_t2, ubd_t2, lbd_t3, ubd_t3, lbd_t4, ubd_t4, lbd_t5, ubd_t5, lbd_t6, ubd_t6, lbd_t7, ubd_t7, lbd_t8, ubd_t8, lbd_t9, ubd_t9;
 MPI_Init(NULL, NULL);
 MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
 MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
 int receiver_list[nprocs];
 for (__p=0; __p<nprocs; __p++) { receiver_list[__p] = 0; }
MPI_Request reqs[2*nprocs];
MPI_Status stats[2*nprocs];
 double t_comm_start, t_comm = 0.0, t_globalcomm = 0.0;
 double t_comp_mean, t_comp_sd;
 double t_comp_start, t_comp = 0.0, t_globalcomp = 0.0;
 double t_pack_start, t_pack = 0.0, t_globalpack = 0.0;
 double t_writeout = 0.0;
 double t_unpack_start, t_unpack = 0.0, t_globalunpack = 0.0;
 double t_local = 0.0, t_global = 0.0;
 double __total_count = 0;
 double __total_count_all = 0;
int lw_count_prod1 = 0;
int send_counts_prod1[nprocs];
int recv_counts_prod1[nprocs];
int lw_recv_counts_prod1[nprocs];
int displs_prod1[nprocs];
int displs_lw_prod1[nprocs];
int send_count_prod1 = 0;
int curr_displs_prod1[nprocs];
double *send_buf_prod1 = NULL;
size_t send_buf_size_prod1 = 0;
double *recv_buf_prod1 = NULL;
size_t recv_buf_size_prod1 = 0;
double *lw_buf_prod1 = NULL;
size_t lw_buf_size_prod1 = 0;
double *lw_recv_buf_prod1 = NULL;
size_t lw_recv_buf_size_prod1 = 0;
int lw_count_temp = 0;
int send_counts_temp[nprocs];
int recv_counts_temp[nprocs];
int lw_recv_counts_temp[nprocs];
int displs_temp[nprocs];
int displs_lw_temp[nprocs];
int send_count_temp = 0;
int curr_displs_temp[nprocs];
double *send_buf_temp = NULL;
size_t send_buf_size_temp = 0;
double *recv_buf_temp = NULL;
size_t recv_buf_size_temp = 0;
double *lw_buf_temp = NULL;
size_t lw_buf_size_temp = 0;
double *lw_recv_buf_temp = NULL;
size_t lw_recv_buf_size_temp = 0;
int lw_count_prod2 = 0;
int send_counts_prod2[nprocs];
int recv_counts_prod2[nprocs];
int lw_recv_counts_prod2[nprocs];
int displs_prod2[nprocs];
int displs_lw_prod2[nprocs];
int send_count_prod2 = 0;
int curr_displs_prod2[nprocs];
double *send_buf_prod2 = NULL;
size_t send_buf_size_prod2 = 0;
double *recv_buf_prod2 = NULL;
size_t recv_buf_size_prod2 = 0;
double *lw_buf_prod2 = NULL;
size_t lw_buf_size_prod2 = 0;
double *lw_recv_buf_prod2 = NULL;
size_t lw_recv_buf_size_prod2 = 0;
int lw_count_conv = 0;
int send_counts_conv[nprocs];
int recv_counts_conv[nprocs];
int lw_recv_counts_conv[nprocs];
int displs_conv[nprocs];
int displs_lw_conv[nprocs];
int send_count_conv = 0;
int curr_displs_conv[nprocs];
double *send_buf_conv = NULL;
size_t send_buf_size_conv = 0;
double *recv_buf_conv = NULL;
size_t recv_buf_size_conv = 0;
double *lw_buf_conv = NULL;
size_t lw_buf_size_conv = 0;
double *lw_recv_buf_conv = NULL;
size_t lw_recv_buf_size_conv = 0;
send_buf_prod1 = (double *) polyrt_max_alloc(send_buf_prod1, sizeof(double)*((+1)*floorf((1*(2+nprocs))/(float)nprocs)), &send_buf_size_prod1);
recv_buf_prod1 = (double *) polyrt_max_alloc(recv_buf_prod1, nprocs*send_buf_size_prod1, &recv_buf_size_prod1);
 for (__p=0; __p<nprocs; __p++) { displs_prod1[__p] = __p*send_buf_size_prod1/sizeof(double);}
if (fopen(".test", "r")) {lw_buf_prod1 = (double *) polyrt_max_alloc(lw_buf_prod1, sizeof(double)*((+1)*floorf((1*(2+nprocs))/(float)nprocs)), &lw_buf_size_prod1);
lw_recv_buf_prod1 = (double *) polyrt_max_alloc(lw_recv_buf_prod1, nprocs*lw_buf_size_prod1, &lw_recv_buf_size_prod1);
 for (__p=0; __p<nprocs; __p++) { displs_lw_prod1[__p] = __p*lw_buf_size_prod1/sizeof(double);}
}send_buf_temp = (double *) polyrt_max_alloc(send_buf_temp, sizeof(double)*((+1*(1)*(1)*(3))*floorf((1*(2+nprocs))/(float)nprocs)), &send_buf_size_temp);
recv_buf_temp = (double *) polyrt_max_alloc(recv_buf_temp, nprocs*send_buf_size_temp, &recv_buf_size_temp);
 for (__p=0; __p<nprocs; __p++) { displs_temp[__p] = __p*send_buf_size_temp/sizeof(double);}
if (fopen(".test", "r")) {lw_buf_temp = (double *) polyrt_max_alloc(lw_buf_temp, sizeof(double)*((+1*(1)*(1)*(3))*floorf((1*(2+nprocs))/(float)nprocs)), &lw_buf_size_temp);
lw_recv_buf_temp = (double *) polyrt_max_alloc(lw_recv_buf_temp, nprocs*lw_buf_size_temp, &lw_recv_buf_size_temp);
 for (__p=0; __p<nprocs; __p++) { displs_lw_temp[__p] = __p*lw_buf_size_temp/sizeof(double);}
}send_buf_prod2 = (double *) polyrt_max_alloc(send_buf_prod2, sizeof(double)*((+1)*floorf((1*(2+nprocs))/(float)nprocs)), &send_buf_size_prod2);
recv_buf_prod2 = (double *) polyrt_max_alloc(recv_buf_prod2, nprocs*send_buf_size_prod2, &recv_buf_size_prod2);
 for (__p=0; __p<nprocs; __p++) { displs_prod2[__p] = __p*send_buf_size_prod2/sizeof(double);}
if (fopen(".test", "r")) {lw_buf_prod2 = (double *) polyrt_max_alloc(lw_buf_prod2, sizeof(double)*((+1)*floorf((1*(2+nprocs))/(float)nprocs)), &lw_buf_size_prod2);
lw_recv_buf_prod2 = (double *) polyrt_max_alloc(lw_recv_buf_prod2, nprocs*lw_buf_size_prod2, &lw_recv_buf_size_prod2);
 for (__p=0; __p<nprocs; __p++) { displs_lw_prod2[__p] = __p*lw_buf_size_prod2/sizeof(double);}
}send_buf_conv = (double *) polyrt_max_alloc(send_buf_conv, sizeof(double)*((+1*(1)*(1)*(0 - max(0,0) + 1))*floorf((1*(2+nprocs))/(float)nprocs)), &send_buf_size_conv);
recv_buf_conv = (double *) polyrt_max_alloc(recv_buf_conv, nprocs*send_buf_size_conv, &recv_buf_size_conv);
 for (__p=0; __p<nprocs; __p++) { displs_conv[__p] = __p*send_buf_size_conv/sizeof(double);}
if (fopen(".test", "r")) {lw_buf_conv = (double *) polyrt_max_alloc(lw_buf_conv, sizeof(double)*((+1*(1)*(1)*(3))*floorf((1*(2+nprocs))/(float)nprocs)), &lw_buf_size_conv);
lw_recv_buf_conv = (double *) polyrt_max_alloc(lw_recv_buf_conv, nprocs*lw_buf_size_conv, &lw_recv_buf_size_conv);
 for (__p=0; __p<nprocs; __p++) { displs_lw_conv[__p] = __p*lw_buf_size_conv/sizeof(double);}
} polyrt_init_grid_size();
 IF_TIME(t_local = rtclock());
/* Start of CLooG code */
for (t1=0;t1<=4994;t1++) {
  for (t2=0;t2<=5037;t2++) {
    t4 = floord(t2+1,126);
    for (t5=ceild(6*t2-6*t4-124,125);t5<=min(floord(2*t2+498*t4+500,125),floord(6*t2-6*t4+6,125));t5++) {
      if ((t2 == 21*t5-1) && (t2 >= 126*t4+41)) {
        if ((t2+1)%42 == 0) {
          prod1 += (src[ t1 + 4][ (t2-t4)][ ((-t2+126*t4+83)/42)] * kernelX[ 4]);;
          /* unknown - failure constructing stmt body */;
          prod2 += (temp[ t1][ (t2-t4-4) + 4][ ((-t2+126*t4+83)/42)] * kernelY[ 4]);;
        }
        if ((t2+22)%42 == 0) {
          for (t8=ceild(250*t2+250,21);t8<=floord(250*t2+271,21);t8++) {
            prod1 += (src[ t1 + ((-250*t2+21*t8-208)/21)][ (t2-t4)][ ((-t2+126*t4+104)/42)] * kernelX[ ((-250*t2+21*t8-208)/21)]);;
            prod2 += (temp[ t1][ (t2-t4-4) + ((-250*t2+21*t8-208)/21)][ ((-t2+126*t4+104)/42)] * kernelY[ ((-250*t2+21*t8-208)/21)]);;
          }
          prod1 += (src[ t1 + 4][ (t2-t4)][ ((-t2+126*t4+104)/42)] * kernelX[ 4]);;
          /* unknown - failure constructing stmt body */;
          prod2 += (temp[ t1][ (t2-t4-4) + 4][ ((-t2+126*t4+104)/42)] * kernelY[ 4]);;
        }
        if ((41*t2+20)%42 <= 21) {
          t7 = floord(83*t2+42*t4+104,42);
          if ((t2+1)%21 == 0) {
            /* unknown - failure constructing stmt body */;
          }
        }
        for (t7=ceild(83*t2+42*t4+125,42);t7<=2*t2-2*t4+2;t7++) {
          if ((t2+1)%21 == 0) {
            /* unknown - failure constructing stmt body */;
          }
          if ((t2+1)%21 == 0) {
            /* unknown - failure constructing stmt body */;
          }
          for (t8=4*t2-4*t4+4*t7;t8<=4*t2-4*t4+4*t7+3;t8++) {
            if ((t2+1)%21 == 0) {
              prod1 += (src[ t1 + (-4*t2+4*t4-4*t7+t8)][ (t2-t4)][ (-2*t2+2*t4+t7)] * kernelX[ (-4*t2+4*t4-4*t7+t8)]);;
            }
            if ((t2+1)%21 == 0) {
              prod2 += (temp[ t1][ (t2-t4-4) + (-4*t2+4*t4-4*t7+t8)][ (-2*t2+2*t4+t7)] * kernelY[ (-4*t2+4*t4-4*t7+t8)]);;
            }
          }
          if ((t2+1)%21 == 0) {
            prod1 += (src[ t1 + 4][ (t2-t4)][ (-2*t2+2*t4+t7)] * kernelX[ 4]);;
          }
          if ((t2+1)%21 == 0) {
            /* unknown - failure constructing stmt body */;
          }
          if ((t2+1)%21 == 0) {
            prod2 += (temp[ t1][ (t2-t4-4) + 4][ (-2*t2+2*t4+t7)] * kernelY[ 4]);;
          }
          if ((t2+1)%21 == 0) {
            /* unknown - failure constructing stmt body */;
          }
        }
      }
      if ((t2 == 21*t5+20) && (t2 >= 126*t4+20) && (t2 <= 126*t4+62)) {
        for (t7=2*t2-2*t4;t7<=floord(83*t2+42*t4+62,42);t7++) {
          if ((t2+1)%21 == 0) {
            /* unknown - failure constructing stmt body */;
          }
          if ((t2+1)%21 == 0) {
            /* unknown - failure constructing stmt body */;
          }
          for (t8=4*t2-4*t4+4*t7;t8<=4*t2-4*t4+4*t7+3;t8++) {
            if ((t2+1)%21 == 0) {
              prod1 += (src[ t1 + (-4*t2+4*t4-4*t7+t8)][ (t2-t4)][ (-2*t2+2*t4+t7)] * kernelX[ (-4*t2+4*t4-4*t7+t8)]);;
            }
            if ((t2+1)%21 == 0) {
              prod2 += (temp[ t1][ (t2-t4-4) + (-4*t2+4*t4-4*t7+t8)][ (-2*t2+2*t4+t7)] * kernelY[ (-4*t2+4*t4-4*t7+t8)]);;
            }
          }
          if ((t2+1)%21 == 0) {
            prod1 += (src[ t1 + 4][ (t2-t4)][ (-2*t2+2*t4+t7)] * kernelX[ 4]);;
          }
          if ((t2+1)%21 == 0) {
            /* unknown - failure constructing stmt body */;
          }
          if ((t2+1)%21 == 0) {
            prod2 += (temp[ t1][ (t2-t4-4) + 4][ (-2*t2+2*t4+t7)] * kernelY[ 4]);;
          }
          if ((t2+1)%21 == 0) {
            /* unknown - failure constructing stmt body */;
          }
        }
        if ((41*t2+20)%42 <= 21) {
          t7 = floord(83*t2+42*t4+104,42);
          if ((t2+1)%21 == 0) {
            /* unknown - failure constructing stmt body */;
          }
          if ((t2+1)%21 == 0) {
            /* unknown - failure constructing stmt body */;
          }
          for (t8=4*t2-4*t4+4*t7;t8<=floord(250*t2+229,21);t8++) {
            if ((t2+1)%21 == 0) {
              prod1 += (src[ t1 + (-4*t2+4*t4-4*t7+t8)][ (t2-t4)][ (-2*t2+2*t4+t7)] * kernelX[ (-4*t2+4*t4-4*t7+t8)]);;
            }
            if ((t2+1)%21 == 0) {
              prod2 += (temp[ t1][ (t2-t4-4) + (-4*t2+4*t4-4*t7+t8)][ (-2*t2+2*t4+t7)] * kernelY[ (-4*t2+4*t4-4*t7+t8)]);;
            }
          }
        }
      }
      if ((t2 <= min(t4+4994,21*t5+19)) && (t2 >= max(21*t5,t4+4))) {
        for (t7=2*t2-2*t4;t7<=min(250*t4+249,2*t2-2*t4+2);t7++) {
          /* unknown - failure constructing stmt body */;
          /* unknown - failure constructing stmt body */;
          for (t8=4*t2-4*t4+4*t7;t8<=4*t2-4*t4+4*t7+3;t8++) {
            prod1 += (src[ t1 + (-4*t2+4*t4-4*t7+t8)][ (t2-t4)][ (-2*t2+2*t4+t7)] * kernelX[ (-4*t2+4*t4-4*t7+t8)]);;
            prod2 += (temp[ t1][ (t2-t4-4) + (-4*t2+4*t4-4*t7+t8)][ (-2*t2+2*t4+t7)] * kernelY[ (-4*t2+4*t4-4*t7+t8)]);;
          }
          prod1 += (src[ t1 + 4][ (t2-t4)][ (-2*t2+2*t4+t7)] * kernelX[ 4]);;
          /* unknown - failure constructing stmt body */;
          prod2 += (temp[ t1][ (t2-t4-4) + 4][ (-2*t2+2*t4+t7)] * kernelY[ 4]);;
          /* unknown - failure constructing stmt body */;
        }
      }
      if ((t2 <= 3) && (t4 == 0) && (t5 == 0)) {
        for (t7=2*t2;t7<=2*t2+2;t7++) {
          /* unknown - failure constructing stmt body */;
          for (t8=4*t2+4*t7;t8<=4*t2+4*t7+4;t8++) {
            prod1 += (src[ t1 + (-4*t2-4*t7+t8)][ t2][ (-2*t2+t7)] * kernelX[ (-4*t2-4*t7+t8)]);;
          }
          /* unknown - failure constructing stmt body */;
        }
      }
      if ((t2 == 126*t4-1) && (t2 == 21*t5+20)) {
        if ((125*t2+125)%126 == 0) {
          /* unknown - failure constructing stmt body */;
          /* unknown - failure constructing stmt body */;
          for (t8=ceild(250*t2+166,21);t8<=floord(250*t2+229,21);t8++) {
            prod1 += (src[ t1 + ((-250*t2+21*t8-166)/21)][ ((125*t2-1)/126)][ 2] * kernelX[ ((-250*t2+21*t8-166)/21)]);;
            prod2 += (temp[ t1][ ((125*t2-505)/126) + ((-250*t2+21*t8-166)/21)][ 2] * kernelY[ ((-250*t2+21*t8-166)/21)]);;
          }
        }
      }
      if ((t2 == 21*t5+20) && (t2 >= 126*t4+83)) {
        if ((t2+1)%21 == 0) {
          /* unknown - failure constructing stmt body */;
        }
        if ((t2+1)%21 == 0) {
          /* unknown - failure constructing stmt body */;
        }
        for (t8=12*t2-12*t4;t8<=floord(250*t2+229,21);t8++) {
          if ((t2+1)%21 == 0) {
            prod1 += (src[ t1 + (-12*t2+12*t4+t8)][ (t2-t4)][ 0] * kernelX[ (-12*t2+12*t4+t8)]);;
          }
          if ((t2+1)%21 == 0) {
            prod2 += (temp[ t1][ (t2-t4-4) + (-12*t2+12*t4+t8)][ 0] * kernelY[ (-12*t2+12*t4+t8)]);;
          }
        }
      }
      if ((t2 == 126*t4+20) && (t2 == 21*t5-1)) {
        if ((125*t2+20)%126 == 0) {
          for (t8=ceild(250*t2+250,21);t8<=floord(250*t2+271,21);t8++) {
            prod1 += (src[ t1 + ((-250*t2+21*t8-208)/21)][ ((125*t2+20)/126)][ 2] * kernelX[ ((-250*t2+21*t8-208)/21)]);;
            prod2 += (temp[ t1][ ((125*t2-484)/126) + ((-250*t2+21*t8-208)/21)][ 2] * kernelY[ ((-250*t2+21*t8-208)/21)]);;
          }
          prod1 += (src[ t1 + 4][ ((125*t2+20)/126)][ 2] * kernelX[ 4]);;
          /* unknown - failure constructing stmt body */;
          prod2 += (temp[ t1][ ((125*t2-484)/126) + 4][ 2] * kernelY[ 4]);;
        }
      }
      if ((t2 == 126*t4-1) && (t2 == 21*t5-1)) {
        if ((125*t2+125)%126 == 0) {
          prod1 += (src[ t1 + 4][ ((125*t2-1)/126)][ 2] * kernelX[ 4]);;
          /* unknown - failure constructing stmt body */;
          prod2 += (temp[ t1][ ((125*t2-505)/126) + 4][ 2] * kernelY[ 4]);;
        }
      }
      if ((t2 == 21*t5-1) && (t2 <= 126*t4+20)) {
        if ((t2+1)%21 == 0) {
          /* unknown - failure constructing stmt body */;
        }
      }
      if ((t2 >= 5034) && (t4 == 39) && (t5 == 239)) {
        for (t7=2*t2-78;t7<=2*t2-76;t7++) {
          /* unknown - failure constructing stmt body */;
          for (t8=4*t2+4*t7-156;t8<=4*t2+4*t7-152;t8++) {
            prod2 += (temp[ t1][ (t2-43) + (-4*t2-4*t7+t8+156)][ (-2*t2+t7+78)] * kernelY[ (-4*t2-4*t7+t8+156)]);;
          }
          /* unknown - failure constructing stmt body */;
        }
      }
    }
    t4 = floord(t2+1,126);
    ;
    for (__p=0; __p<nprocs; __p++) { receiver_list[__p] = 0; };
    t4 = floord(t2+1,126);
    sigma_prod1_0(t1,t2,t4, my_rank, nprocs, receiver_list);
    IF_TIME(t_comm_start = rtclock()); for (__p=0; __p<nprocs; __p++) { send_counts_prod1[__p] = receiver_list[__p]? send_count_prod1: 0; } MPI_Alltoall(send_counts_prod1, 1, MPI_INT, recv_counts_prod1, 1, MPI_INT, MPI_COMM_WORLD); req_count=0; for (__p=0; __p<nprocs; __p++) { if(send_counts_prod1[__p] >= 1) {IF_TIME(__total_count += send_count_prod1); MPI_Isend(send_buf_prod1, send_count_prod1, MPI_DOUBLE, __p, 123, MPI_COMM_WORLD, &reqs[req_count++]);}}for (__p=0; __p<nprocs; __p++) { if(recv_counts_prod1[__p] >= 1) { MPI_Irecv(recv_buf_prod1+displs_prod1[__p], recv_counts_prod1[__p], MPI_DOUBLE, __p, 123, MPI_COMM_WORLD, &reqs[req_count++]);}} MPI_Waitall(req_count, reqs, stats); send_count_prod1 = 0; for (__p=0; __p<nprocs; __p++) curr_displs_prod1[__p] = 0;IF_TIME(t_comm += rtclock() - t_comm_start);;
    t4 = floord(t2+1,126);
    ;
    t4 = floord(t2+1,126);
    for (__p=0; __p<nprocs; __p++) { receiver_list[__p] = 0; } sigma_temp_0(t1,t2,t4, my_rank, nprocs, receiver_list); for (__p=0; __p<nprocs; __p++) { if (receiver_list[__p] != 0) { send_count_temp = pack_temp_0(t1,t2,t4,send_buf_temp,send_count_temp); break;} };
    for (__p=0; __p<nprocs; __p++) { receiver_list[__p] = 0; };
    t4 = floord(t2+1,126);
    sigma_temp_0(t1,t2,t4, my_rank, nprocs, receiver_list);
    IF_TIME(t_comm_start = rtclock()); for (__p=0; __p<nprocs; __p++) { send_counts_temp[__p] = receiver_list[__p]? send_count_temp: 0; } MPI_Alltoall(send_counts_temp, 1, MPI_INT, recv_counts_temp, 1, MPI_INT, MPI_COMM_WORLD); req_count=0; for (__p=0; __p<nprocs; __p++) { if(send_counts_temp[__p] >= 1) {IF_TIME(__total_count += send_count_temp); MPI_Isend(send_buf_temp, send_count_temp, MPI_DOUBLE, __p, 123, MPI_COMM_WORLD, &reqs[req_count++]);}}for (__p=0; __p<nprocs; __p++) { if(recv_counts_temp[__p] >= 1) { MPI_Irecv(recv_buf_temp+displs_temp[__p], recv_counts_temp[__p], MPI_DOUBLE, __p, 123, MPI_COMM_WORLD, &reqs[req_count++]);}} MPI_Waitall(req_count, reqs, stats); send_count_temp = 0; for (__p=0; __p<nprocs; __p++) curr_displs_temp[__p] = 0;IF_TIME(t_comm += rtclock() - t_comm_start);;
    t4 = floord(t2+1,126);
    proc = pi_0(t1,t2,t4, nprocs); if ((proc != my_rank) && (recv_counts_temp[proc] > 0)) { for (__p=0; __p<nprocs; __p++) { receiver_list[__p] = 0; } sigma_temp_0(t1,t2,t4, proc, nprocs, receiver_list); for (__p=0; __p<nprocs; __p++) { if (receiver_list[__p] != 0) { curr_displs_temp[proc] = unpack_temp_0(t1,t2,t4,recv_buf_temp,displs_temp[proc],curr_displs_temp[proc]); break; } } };
    t4 = floord(t2+1,126);
    ;
    for (__p=0; __p<nprocs; __p++) { receiver_list[__p] = 0; };
    t4 = floord(t2+1,126);
    sigma_prod2_0(t1,t2,t4, my_rank, nprocs, receiver_list);
    IF_TIME(t_comm_start = rtclock()); for (__p=0; __p<nprocs; __p++) { send_counts_prod2[__p] = receiver_list[__p]? send_count_prod2: 0; } MPI_Alltoall(send_counts_prod2, 1, MPI_INT, recv_counts_prod2, 1, MPI_INT, MPI_COMM_WORLD); req_count=0; for (__p=0; __p<nprocs; __p++) { if(send_counts_prod2[__p] >= 1) {IF_TIME(__total_count += send_count_prod2); MPI_Isend(send_buf_prod2, send_count_prod2, MPI_DOUBLE, __p, 123, MPI_COMM_WORLD, &reqs[req_count++]);}}for (__p=0; __p<nprocs; __p++) { if(recv_counts_prod2[__p] >= 1) { MPI_Irecv(recv_buf_prod2+displs_prod2[__p], recv_counts_prod2[__p], MPI_DOUBLE, __p, 123, MPI_COMM_WORLD, &reqs[req_count++]);}} MPI_Waitall(req_count, reqs, stats); send_count_prod2 = 0; for (__p=0; __p<nprocs; __p++) curr_displs_prod2[__p] = 0;IF_TIME(t_comm += rtclock() - t_comm_start);;
    t4 = floord(t2+1,126);
    ;
    t4 = floord(t2+1,126);
    ;
    for (__p=0; __p<nprocs; __p++) { receiver_list[__p] = 0; };
    t4 = floord(t2+1,126);
    sigma_conv_0(t1,t2,t4, my_rank, nprocs, receiver_list);
    IF_TIME(t_comm_start = rtclock()); for (__p=0; __p<nprocs; __p++) { send_counts_conv[__p] = receiver_list[__p]? send_count_conv: 0; } MPI_Alltoall(send_counts_conv, 1, MPI_INT, recv_counts_conv, 1, MPI_INT, MPI_COMM_WORLD); req_count=0; for (__p=0; __p<nprocs; __p++) { if(send_counts_conv[__p] >= 1) {IF_TIME(__total_count += send_count_conv); MPI_Isend(send_buf_conv, send_count_conv, MPI_DOUBLE, __p, 123, MPI_COMM_WORLD, &reqs[req_count++]);}}for (__p=0; __p<nprocs; __p++) { if(recv_counts_conv[__p] >= 1) { MPI_Irecv(recv_buf_conv+displs_conv[__p], recv_counts_conv[__p], MPI_DOUBLE, __p, 123, MPI_COMM_WORLD, &reqs[req_count++]);}} MPI_Waitall(req_count, reqs, stats); send_count_conv = 0; for (__p=0; __p<nprocs; __p++) curr_displs_conv[__p] = 0;IF_TIME(t_comm += rtclock() - t_comm_start);;
    t4 = floord(t2+1,126);
    ;
  }
}
/* End of CLooG code */
 IF_TIME(t_local = rtclock() - t_local);
 IF_TIME(t_writeout = rtclock());
 
#ifndef USE_LOCAL_ARRAYS
 if (fopen(".test", "r")) {
 write_out(my_rank,nprocs,lw_buf_prod1,lw_recv_buf_prod1,displs_lw_prod1,lw_buf_temp,lw_recv_buf_temp,displs_lw_temp,lw_buf_prod2,lw_recv_buf_prod2,displs_lw_prod2,lw_buf_conv,lw_recv_buf_conv,displs_lw_conv);
} 
#endif
 IF_TIME(t_writeout = rtclock() - t_writeout);
#ifdef TIME
 char *buffer = (char *)malloc(4096 * sizeof(char));
 strcpy(buffer, "");
 if (fopen(".detailed_time_report", "r")) {
 sprintf(buffer+strlen(buffer), "Node %d\n", my_rank);
 sprintf(buffer+strlen(buffer), "Write-out time: %0.6lf s\n", t_writeout);
 sprintf(buffer+strlen(buffer), "Values communicated: %0.0lf\n", __total_count);
 sprintf(buffer+strlen(buffer), "Computation time: %0.6lf s\n", t_comp);
 sprintf(buffer+strlen(buffer), "Communication time: %0.6lf s\n", t_comm);
 sprintf(buffer+strlen(buffer), "Packing time: %0.6lf s\n", t_pack);
 sprintf(buffer+strlen(buffer), "Unpacking time: %0.6lf s\n", t_unpack);
 sprintf(buffer+strlen(buffer), "Sum of split times: %0.6lf s\n", t_comp + t_comm + t_pack + t_unpack);
 sprintf(buffer+strlen(buffer), "Total time minus write-out time: %0.6lf s\n", t_local);
 sprintf(buffer+strlen(buffer), "-------");
 fprintf(stdout, "%s\n", buffer);
 }
 MPI_Reduce(&__total_count, &__total_count_all, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
 MPI_Allreduce(&t_comp, &t_comp_mean, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
 t_comp_mean = t_comp_mean/nprocs;
 t_comp_sd = (t_comp-t_comp_mean)*(t_comp-t_comp_mean);
 MPI_Allreduce(MPI_IN_PLACE, &t_comp_sd, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
 t_comp_sd = sqrt(t_comp_sd/nprocs);
 MPI_Reduce(&t_comp, &t_globalcomp, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
 MPI_Reduce(&t_comm, &t_globalcomm, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
 MPI_Reduce(&t_pack, &t_globalpack, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
 MPI_Reduce(&t_unpack, &t_globalunpack, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
 MPI_Reduce(&t_local, &t_global, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
 if (my_rank==0) {
  strcpy(buffer, "SUMMARY\n");
  sprintf(buffer+strlen(buffer), "Write-out time spent in master node: %0.6lf s\n", t_writeout);
  sprintf(buffer+strlen(buffer), "Total communication volume across all nodes: %0.6lf GB\n", __total_count_all*8/((double)1024*1024*1024));
  sprintf(buffer+strlen(buffer), "Maximum computation time spent across all nodes: %0.6lf s\n", t_globalcomp);
  sprintf(buffer+strlen(buffer), "Mean computation time across all nodes (not per-thread): %0.6lf s\n", t_comp_mean);
  sprintf(buffer+strlen(buffer), "Standard deviation in computation time across all nodes (not per-thread): %0.6lf s\n", t_comp_sd);
  sprintf(buffer+strlen(buffer), "Load imbalance factor (standard-deviation/mean): %0.2lf %%\n", 100*t_comp_sd/t_comp_mean);
  sprintf(buffer+strlen(buffer), "Maximum communication time spent across all nodes: %0.6lf s\n", t_globalcomm);
  sprintf(buffer+strlen(buffer), "Maximum packing time spent across all nodes: %0.6lf s\n", t_globalpack);
  sprintf(buffer+strlen(buffer), "Maximum unpacking time spent across all nodes: %0.6lf s\n", t_globalunpack);
  sprintf(buffer+strlen(buffer), "Maximum total time spent across all nodes: %0.6lf s\n", t_global);
  sprintf(buffer+strlen(buffer), "-------");
  fprintf(stdout, "%s\n", buffer);
 }
 free(buffer);
#endif
 MPI_Barrier(MPI_COMM_WORLD);
 MPI_Finalize();

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
/* Copyright (C) 1991-2014 Free Software Foundation, Inc.
   This file is part of the GNU C Library.

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the GNU C Library; if not, see
   <http://www.gnu.org/licenses/>.  */
/* This header is separate from features.h so that the compiler can
   include it implicitly at the start of every compilation.  It must
   not itself include <features.h> or any other header that includes
   <features.h> because the implicit include comes before any feature
   test macros that may be defined in a source file before it first
   explicitly includes a system header.  GCC knows the name of this
   header in order to preinclude it.  */
/* glibc's intent is to support the IEC 559 math functionality, real
   and complex.  If the GCC (4.9 and later) predefined macros
   specifying compiler intent are available, use them to determine
   whether the overall intent is to support these features; otherwise,
   presume an older compiler has intent to support these features and
   define these macros by default.  */
/* wchar_t uses ISO/IEC 10646 (2nd ed., published 2011-03-15) /
   Unicode 6.0.  */
/* We do not support C11 <threads.h>.  */
int writeout_pack_prod1_0(int ts1,int ts2,int ts3,double *lw_buf_prod1,int lw_count_prod1){
int recv_proc;
  int t1, t2, t3;
  int lb, ub, lbd, ubd, lb2, ub2;
  register int lbv, ubv;
/* Start of CLooG code */
if ((ts1 >= 0) && (ts1 <= 4994) && (ts2 >= ts3) && (ts2 <= ts3+4994) && (ts2 >= 126*ts3-1) && (ts2 <= 126*ts3+124)) {
  lw_buf_prod1[lw_count_prod1++] = prod1;
}
/* End of CLooG code */
return lw_count_prod1;
}
int writeout_unpack_prod1_0(int ts1,int ts2,int ts3,double *lw_recv_buf_prod1,int displs_lw_prod1){
int recv_proc;
  int t1, t2, t3;
  int lb, ub, lbd, ubd, lb2, ub2;
  register int lbv, ubv;
/* Start of CLooG code */
if ((ts1 >= 0) && (ts1 <= 4994) && (ts2 >= ts3) && (ts2 <= ts3+4994) && (ts2 >= 126*ts3-1) && (ts2 <= 126*ts3+124)) {
  
#ifndef USE_LOCAL_ARRAYS
 prod1 = lw_recv_buf_prod1[displs_lw_prod1++]; 
#endif
;
}
/* End of CLooG code */
return displs_lw_prod1;
}
int pack_temp_0(int ts1,int ts2,int ts3,double *send_buf_temp,int send_count_temp){
int recv_proc;
  int t1, t2, t3;
  int lb, ub, lbd, ubd, lb2, ub2;
  register int lbv, ubv;
/* Start of CLooG code */
if ((ts1 >= 0) && (ts1 <= 4994) && (ts2 >= ts3) && (ts2 <= ts3+4994) && (ts2 >= 126*ts3-1) && (ts2 <= 126*ts3+124)) {
  for (t3=max(0,-2*ts2+252*ts3);t3<=min(2,-2*ts2+252*ts3+249);t3++) {
    send_buf_temp[send_count_temp++] = temp[ts1][(ts2-ts3)][t3];
  }
}
/* End of CLooG code */
return send_count_temp;
}
int unpack_temp_0(int ts1,int ts2,int ts3,double *recv_buf_temp,int displs_temp,int count){
int recv_proc;
  int t1, t2, t3;
  int lb, ub, lbd, ubd, lb2, ub2;
  register int lbv, ubv;
/* Start of CLooG code */
if ((ts1 >= 0) && (ts1 <= 4994) && (ts2 >= ts3) && (ts2 <= ts3+4994) && (ts2 >= 126*ts3-1) && (ts2 <= 126*ts3+124)) {
  for (t3=max(0,-2*ts2+252*ts3);t3<=min(2,-2*ts2+252*ts3+249);t3++) {
    temp[ts1][(ts2-ts3)][t3] = recv_buf_temp[displs_temp + count++];;
  }
}
/* End of CLooG code */
return count;
}
int writeout_pack_temp_0(int ts1,int ts2,int ts3,double *lw_buf_temp,int lw_count_temp){
int recv_proc;
  int t1, t2, t3;
  int lb, ub, lbd, ubd, lb2, ub2;
  register int lbv, ubv;
/* Start of CLooG code */
if ((ts1 >= 0) && (ts1 <= 4994) && (ts2 >= ts3) && (ts2 <= ts3+4994) && (ts2 >= 126*ts3-1) && (ts2 <= 126*ts3+124)) {
  for (t3=max(0,-2*ts2+252*ts3);t3<=min(2,-2*ts2+252*ts3+249);t3++) {
    lw_buf_temp[lw_count_temp++] = temp[ts1][(ts2-ts3)][t3];
  }
}
/* End of CLooG code */
return lw_count_temp;
}
int writeout_unpack_temp_0(int ts1,int ts2,int ts3,double *lw_recv_buf_temp,int displs_lw_temp){
int recv_proc;
  int t1, t2, t3;
  int lb, ub, lbd, ubd, lb2, ub2;
  register int lbv, ubv;
/* Start of CLooG code */
if ((ts1 >= 0) && (ts1 <= 4994) && (ts2 >= ts3) && (ts2 <= ts3+4994) && (ts2 >= 126*ts3-1) && (ts2 <= 126*ts3+124)) {
  for (t3=max(0,-2*ts2+252*ts3);t3<=min(2,-2*ts2+252*ts3+249);t3++) {
    
#ifndef USE_LOCAL_ARRAYS
 temp[ts1][(ts2-ts3)][t3] = lw_recv_buf_temp[displs_lw_temp++]; 
#endif
;
  }
}
/* End of CLooG code */
return displs_lw_temp;
}
int writeout_pack_prod2_0(int ts1,int ts2,int ts3,double *lw_buf_prod2,int lw_count_prod2){
int recv_proc;
  int t1, t2, t3;
  int lb, ub, lbd, ubd, lb2, ub2;
  register int lbv, ubv;
/* Start of CLooG code */
if ((ts1 >= 0) && (ts1 <= 4994) && (ts2 >= ts3+4) && (ts2 <= ts3+4998) && (ts2 >= 126*ts3-1) && (ts2 <= 126*ts3+124)) {
  lw_buf_prod2[lw_count_prod2++] = prod2;
}
/* End of CLooG code */
return lw_count_prod2;
}
int writeout_unpack_prod2_0(int ts1,int ts2,int ts3,double *lw_recv_buf_prod2,int displs_lw_prod2){
int recv_proc;
  int t1, t2, t3;
  int lb, ub, lbd, ubd, lb2, ub2;
  register int lbv, ubv;
/* Start of CLooG code */
if ((ts1 >= 0) && (ts1 <= 4994) && (ts2 >= ts3+4) && (ts2 <= ts3+4998) && (ts2 >= 126*ts3-1) && (ts2 <= 126*ts3+124)) {
  
#ifndef USE_LOCAL_ARRAYS
 prod2 = lw_recv_buf_prod2[displs_lw_prod2++]; 
#endif
;
}
/* End of CLooG code */
return displs_lw_prod2;
}
int writeout_pack_conv_0(int ts1,int ts2,int ts3,double *lw_buf_conv,int lw_count_conv){
int recv_proc;
  int t1, t2, t3;
  int lb, ub, lbd, ubd, lb2, ub2;
  register int lbv, ubv;
/* Start of CLooG code */
if ((ts1 >= 0) && (ts1 <= 4994) && (ts2 >= ts3+4) && (ts2 <= ts3+4998) && (ts2 >= 126*ts3-1) && (ts2 <= 126*ts3+124)) {
  for (t3=max(0,-2*ts2+252*ts3);t3<=min(2,-2*ts2+252*ts3+249);t3++) {
    lw_buf_conv[lw_count_conv++] = conv[ts1][(ts2-ts3-4)][t3];
  }
}
/* End of CLooG code */
return lw_count_conv;
}
int writeout_unpack_conv_0(int ts1,int ts2,int ts3,double *lw_recv_buf_conv,int displs_lw_conv){
int recv_proc;
  int t1, t2, t3;
  int lb, ub, lbd, ubd, lb2, ub2;
  register int lbv, ubv;
/* Start of CLooG code */
if ((ts1 >= 0) && (ts1 <= 4994) && (ts2 >= ts3+4) && (ts2 <= ts3+4998) && (ts2 >= 126*ts3-1) && (ts2 <= 126*ts3+124)) {
  for (t3=max(0,-2*ts2+252*ts3);t3<=min(2,-2*ts2+252*ts3+249);t3++) {
    
#ifndef USE_LOCAL_ARRAYS
 conv[ts1][(ts2-ts3-4)][t3] = lw_recv_buf_conv[displs_lw_conv++]; 
#endif
;
  }
}
/* End of CLooG code */
return displs_lw_conv;
}
void write_out(int my_rank,int nprocs,double *lw_buf_prod1,double *lw_recv_buf_prod1,int *displs_lw_prod1,double *lw_buf_temp,double *lw_recv_buf_temp,int *displs_lw_temp,double *lw_buf_prod2,double *lw_recv_buf_prod2,int *displs_lw_prod2,double *lw_buf_conv,double *lw_recv_buf_conv,int *displs_lw_conv)
{
int __p, proc;
int _lb_dist, _ub_dist, lbd_t1, ubd_t1, lbd_t2, ubd_t2, lbd_t3, ubd_t3, lbd_t4, ubd_t4, lbd_t5, ubd_t5, lbd_t6, ubd_t6, lbd_t7, ubd_t7, lbd_t8, ubd_t8, lbd_t9, ubd_t9;
int lw_recv_counts_prod1[nprocs];
int curr_displs_lw_prod1[nprocs];
int lw_count_prod1 = 0;
int lw_recv_counts_temp[nprocs];
int curr_displs_lw_temp[nprocs];
int lw_count_temp = 0;
int lw_recv_counts_prod2[nprocs];
int curr_displs_lw_prod2[nprocs];
int lw_count_prod2 = 0;
int lw_recv_counts_conv[nprocs];
int curr_displs_lw_conv[nprocs];
int lw_count_conv = 0;
  int t1, t2, t3, t4;
  int lb, ub, lbd, ubd, lb2, ub2;
  register int lbv, ubv;
/* Start of CLooG code */
for (t1=0;t1<=4994;t1++) {
  for (t2=0;t2<=5037;t2++) {
    t4 = floord(t2+1,126);
    if (pi_0(t1,t2,t4,nprocs) == my_rank)lw_count_prod1 = writeout_pack_prod1_0(t1,t2,t4,lw_buf_prod1,lw_count_prod1);;
    assert((nprocs == 1) || (lw_count_prod1 <= displs_lw_prod1[1])); MPI_Gather(&lw_count_prod1, 1, MPI_INT, lw_recv_counts_prod1, 1, MPI_INT, 0, MPI_COMM_WORLD);MPI_Gatherv(lw_buf_prod1, lw_count_prod1, MPI_DOUBLE, lw_recv_buf_prod1, lw_recv_counts_prod1, displs_lw_prod1, MPI_DOUBLE, 0, MPI_COMM_WORLD); lw_count_prod1 = 0; for (__p=0; __p<nprocs; __p++) curr_displs_lw_prod1[__p] = displs_lw_prod1[__p];;
    if (my_rank == 0) {
      t4 = floord(t2+1,126);
      proc = pi_0(t1, t2, t4, nprocs);if ((my_rank != proc) && (lw_recv_counts_prod1[proc] > 0)) {curr_displs_lw_prod1[proc] = writeout_unpack_prod1_0(t1,t2,t4,lw_recv_buf_prod1,curr_displs_lw_prod1[proc]);};
    }
    t4 = floord(t2+1,126);
    if (pi_0(t1,t2,t4,nprocs) == my_rank)lw_count_temp = writeout_pack_temp_0(t1,t2,t4,lw_buf_temp,lw_count_temp);;
    assert((nprocs == 1) || (lw_count_temp <= displs_lw_temp[1])); MPI_Gather(&lw_count_temp, 1, MPI_INT, lw_recv_counts_temp, 1, MPI_INT, 0, MPI_COMM_WORLD);MPI_Gatherv(lw_buf_temp, lw_count_temp, MPI_DOUBLE, lw_recv_buf_temp, lw_recv_counts_temp, displs_lw_temp, MPI_DOUBLE, 0, MPI_COMM_WORLD); lw_count_temp = 0; for (__p=0; __p<nprocs; __p++) curr_displs_lw_temp[__p] = displs_lw_temp[__p];;
    if (my_rank == 0) {
      t4 = floord(t2+1,126);
      proc = pi_0(t1, t2, t4, nprocs);if ((my_rank != proc) && (lw_recv_counts_temp[proc] > 0)) {curr_displs_lw_temp[proc] = writeout_unpack_temp_0(t1,t2,t4,lw_recv_buf_temp,curr_displs_lw_temp[proc]);};
    }
    t4 = floord(t2+1,126);
    if (pi_0(t1,t2,t4,nprocs) == my_rank)lw_count_prod2 = writeout_pack_prod2_0(t1,t2,t4,lw_buf_prod2,lw_count_prod2);;
    assert((nprocs == 1) || (lw_count_prod2 <= displs_lw_prod2[1])); MPI_Gather(&lw_count_prod2, 1, MPI_INT, lw_recv_counts_prod2, 1, MPI_INT, 0, MPI_COMM_WORLD);MPI_Gatherv(lw_buf_prod2, lw_count_prod2, MPI_DOUBLE, lw_recv_buf_prod2, lw_recv_counts_prod2, displs_lw_prod2, MPI_DOUBLE, 0, MPI_COMM_WORLD); lw_count_prod2 = 0; for (__p=0; __p<nprocs; __p++) curr_displs_lw_prod2[__p] = displs_lw_prod2[__p];;
    if (my_rank == 0) {
      t4 = floord(t2+1,126);
      proc = pi_0(t1, t2, t4, nprocs);if ((my_rank != proc) && (lw_recv_counts_prod2[proc] > 0)) {curr_displs_lw_prod2[proc] = writeout_unpack_prod2_0(t1,t2,t4,lw_recv_buf_prod2,curr_displs_lw_prod2[proc]);};
    }
    t4 = floord(t2+1,126);
    if (pi_0(t1,t2,t4,nprocs) == my_rank)lw_count_conv = writeout_pack_conv_0(t1,t2,t4,lw_buf_conv,lw_count_conv);;
    assert((nprocs == 1) || (lw_count_conv <= displs_lw_conv[1])); MPI_Gather(&lw_count_conv, 1, MPI_INT, lw_recv_counts_conv, 1, MPI_INT, 0, MPI_COMM_WORLD);MPI_Gatherv(lw_buf_conv, lw_count_conv, MPI_DOUBLE, lw_recv_buf_conv, lw_recv_counts_conv, displs_lw_conv, MPI_DOUBLE, 0, MPI_COMM_WORLD); lw_count_conv = 0; for (__p=0; __p<nprocs; __p++) curr_displs_lw_conv[__p] = displs_lw_conv[__p];;
    if (my_rank == 0) {
      t4 = floord(t2+1,126);
      proc = pi_0(t1, t2, t4, nprocs);if ((my_rank != proc) && (lw_recv_counts_conv[proc] > 0)) {curr_displs_lw_conv[proc] = writeout_unpack_conv_0(t1,t2,t4,lw_recv_buf_conv,curr_displs_lw_conv[proc]);};
    }
  }
}
/* End of CLooG code */
}
