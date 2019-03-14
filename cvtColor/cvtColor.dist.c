#include "cvtColor.dist.h"
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
  int t1, t2, t3, t4, t5;
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
 int _lb_dist, _ub_dist, lbd_t1, ubd_t1, lbd_t2, ubd_t2, lbd_t3, ubd_t3, lbd_t4, ubd_t4, lbd_t5, ubd_t5;
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
int lw_count_dst = 0;
int send_counts_dst[nprocs];
int recv_counts_dst[nprocs];
int lw_recv_counts_dst[nprocs];
int displs_dst[nprocs];
int displs_lw_dst[nprocs];
int send_count_dst = 0;
int curr_displs_dst[nprocs];
double *send_buf_dst = NULL;
size_t send_buf_size_dst = 0;
double *recv_buf_dst = NULL;
size_t recv_buf_size_dst = 0;
double *lw_buf_dst = NULL;
size_t lw_buf_size_dst = 0;
double *lw_recv_buf_dst = NULL;
size_t lw_recv_buf_size_dst = 0;
send_buf_dst = (double *) polyrt_max_alloc(send_buf_dst, sizeof(double)*((+1*(1)*(0 - max(0,0) + 1))*floorf((1*(21+nprocs))/(float)nprocs)), &send_buf_size_dst);
recv_buf_dst = (double *) polyrt_max_alloc(recv_buf_dst, nprocs*send_buf_size_dst, &recv_buf_size_dst);
 for (__p=0; __p<nprocs; __p++) { displs_dst[__p] = __p*send_buf_size_dst/sizeof(double);}
if (fopen(".test", "r")) {lw_buf_dst = (double *) polyrt_max_alloc(lw_buf_dst, sizeof(double)*((+1*(250)*(5000))*floorf((1*(21+nprocs))/(float)nprocs)), &lw_buf_size_dst);
lw_recv_buf_dst = (double *) polyrt_max_alloc(lw_recv_buf_dst, nprocs*lw_buf_size_dst, &lw_recv_buf_size_dst);
 for (__p=0; __p<nprocs; __p++) { displs_lw_dst[__p] = __p*lw_buf_size_dst/sizeof(double);}
} polyrt_init_grid_size();
 IF_TIME(t_local = rtclock());
/* Start of CLooG code */
IF_TIME(t_comp_start = rtclock());
_lb_dist=0;
_ub_dist=19;
polyrt_loop_dist(_lb_dist, _ub_dist, nprocs, my_rank, &lbd_t2, &ubd_t2);
for (t2=lbd_t2;t2<=ubd_t2;t2++) {
  for (t3=0;t3<=19;t3++) {
    for (t4=250*t2;t4<=250*t2+249;t4++) {
      lbv=250*t3;
      ubv=250*t3+249;
#pragma ivdep
#pragma vector always
      for (t5=lbv;t5<=ubv;t5++) {
        dst[ t4][ t5] = (((((src[ t4][ t5][2] * B2Y) + (src[ t4][ t5][1] * G2Y)) + (src[ t4][ t5][0] * R2Y)) + ((1) << (yuv_shift - (1)))) >> yuv_shift);;
      }
    }
  }
}
IF_TIME(t_comp += rtclock() - t_comp_start);
IF_TIME(t_pack_start = rtclock());
_lb_dist=0;
_ub_dist=19;
polyrt_loop_dist(_lb_dist, _ub_dist, nprocs, my_rank, &lbd_t2, &ubd_t2);
for (t2=lbd_t2;t2<=ubd_t2;t2++) {
  ;
}
IF_TIME(t_pack += rtclock() - t_pack_start);
for (__p=0; __p<nprocs; __p++) { receiver_list[__p] = 0; };
IF_TIME(t_comm_start = rtclock());
_lb_dist=0;
_ub_dist=19;
polyrt_loop_dist(_lb_dist, _ub_dist, nprocs, my_rank, &lbd_t2, &ubd_t2);
for (t2=lbd_t2;t2<=ubd_t2;t2++) {
  sigma_dst_0(t2, my_rank, nprocs, receiver_list);
}
IF_TIME(t_comm += rtclock() - t_comm_start);
IF_TIME(t_comm_start = rtclock()); for (__p=0; __p<nprocs; __p++) { send_counts_dst[__p] = receiver_list[__p]? send_count_dst: 0; } MPI_Alltoall(send_counts_dst, 1, MPI_INT, recv_counts_dst, 1, MPI_INT, MPI_COMM_WORLD); req_count=0; for (__p=0; __p<nprocs; __p++) { if(send_counts_dst[__p] >= 1) {IF_TIME(__total_count += send_count_dst); MPI_Isend(send_buf_dst, send_count_dst, MPI_DOUBLE, __p, 123, MPI_COMM_WORLD, &reqs[req_count++]);}}for (__p=0; __p<nprocs; __p++) { if(recv_counts_dst[__p] >= 1) { MPI_Irecv(recv_buf_dst+displs_dst[__p], recv_counts_dst[__p], MPI_DOUBLE, __p, 123, MPI_COMM_WORLD, &reqs[req_count++]);}} MPI_Waitall(req_count, reqs, stats); send_count_dst = 0; for (__p=0; __p<nprocs; __p++) curr_displs_dst[__p] = 0;IF_TIME(t_comm += rtclock() - t_comm_start);;
IF_TIME(t_unpack_start = rtclock());
for (t2=0;t2<=19;t2++) {
  ;
}
IF_TIME(t_unpack += rtclock() - t_unpack_start);
/* End of CLooG code */
 IF_TIME(t_local = rtclock() - t_local);
 IF_TIME(t_writeout = rtclock());
 
#ifndef USE_LOCAL_ARRAYS
 if (fopen(".test", "r")) {
 write_out(my_rank,nprocs,lw_buf_dst,lw_recv_buf_dst,displs_lw_dst);
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
int writeout_pack_dst_0(int ts1,double *lw_buf_dst,int lw_count_dst){
int recv_proc;
  int t1, t2;
  int lb, ub, lbd, ubd, lb2, ub2;
  register int lbv, ubv;
/* Start of CLooG code */
if ((ts1 >= 0) && (ts1 <= 19)) {
  for (t1=250*ts1;t1<=250*ts1+249;t1++) {
    for (t2=0;t2<=4999;t2++) {
      lw_buf_dst[lw_count_dst++] = dst[t1][t2];
    }
  }
}
/* End of CLooG code */
return lw_count_dst;
}
int writeout_unpack_dst_0(int ts1,double *lw_recv_buf_dst,int displs_lw_dst){
int recv_proc;
  int t1, t2;
  int lb, ub, lbd, ubd, lb2, ub2;
  register int lbv, ubv;
/* Start of CLooG code */
if ((ts1 >= 0) && (ts1 <= 19)) {
  for (t1=250*ts1;t1<=250*ts1+249;t1++) {
    for (t2=0;t2<=4999;t2++) {
      
#ifndef USE_LOCAL_ARRAYS
 dst[t1][t2] = lw_recv_buf_dst[displs_lw_dst++]; 
#endif
;
    }
  }
}
/* End of CLooG code */
return displs_lw_dst;
}
void write_out(int my_rank,int nprocs,double *lw_buf_dst,double *lw_recv_buf_dst,int *displs_lw_dst)
{
int __p, proc;
int _lb_dist, _ub_dist, lbd_t1, ubd_t1, lbd_t2, ubd_t2, lbd_t3, ubd_t3, lbd_t4, ubd_t4, lbd_t5, ubd_t5;
int lw_recv_counts_dst[nprocs];
int curr_displs_lw_dst[nprocs];
int lw_count_dst = 0;
  int t1, t2;
  int lb, ub, lbd, ubd, lb2, ub2;
  register int lbv, ubv;
/* Start of CLooG code */
for (t2=0;t2<=19;t2++) {
  if (pi_0(t2,nprocs) == my_rank)lw_count_dst = writeout_pack_dst_0(t2,lw_buf_dst,lw_count_dst);;
}
assert((nprocs == 1) || (lw_count_dst <= displs_lw_dst[1])); MPI_Gather(&lw_count_dst, 1, MPI_INT, lw_recv_counts_dst, 1, MPI_INT, 0, MPI_COMM_WORLD);MPI_Gatherv(lw_buf_dst, lw_count_dst, MPI_DOUBLE, lw_recv_buf_dst, lw_recv_counts_dst, displs_lw_dst, MPI_DOUBLE, 0, MPI_COMM_WORLD); lw_count_dst = 0; for (__p=0; __p<nprocs; __p++) curr_displs_lw_dst[__p] = displs_lw_dst[__p];;
if (my_rank == 0) {
  for (t2=0;t2<=19;t2++) {
    proc = pi_0(t2, nprocs);if ((my_rank != proc) && (lw_recv_counts_dst[proc] > 0)) {curr_displs_lw_dst[proc] = writeout_unpack_dst_0(t2,lw_recv_buf_dst,curr_displs_lw_dst[proc]);};
  }
}
/* End of CLooG code */
}
