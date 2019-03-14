#include "polyrt.h"
void sigma_dst_0(int ts1, int my_rank, int nprocs, int *receiver_list);
int writeout_pack_dst_0(int ts1,double *lw_buf_dst,int lw_count_dst);
int writeout_unpack_dst_0(int ts1,double *lw_recv_buf_dst,int displs_lw_dst);
int pi_0( int t2, int nprocs);
void write_out(int my_rank,int nprocs,double *lw_buf_dst,double *lw_recv_buf_dst,int *displs_lw_dst);
