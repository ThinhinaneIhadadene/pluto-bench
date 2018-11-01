#include "polyrt.h"
void sigma_prod1_0(int ts1, int ts2, int ts3, int my_rank, int nprocs, int *receiver_list);
int writeout_pack_prod1_0(int ts1,int ts2,int ts3,double *lw_buf_prod1,int lw_count_prod1);
int writeout_unpack_prod1_0(int ts1,int ts2,int ts3,double *lw_recv_buf_prod1,int displs_lw_prod1);
void sigma_temp_0(int ts1, int ts2, int ts3, int my_rank, int nprocs, int *receiver_list);
int pack_temp_0(int ts1,int ts2,int ts3,double *send_buf_temp,int send_count_temp);
int unpack_temp_0(int ts1,int ts2,int ts3,double *recv_buf_temp,int displs_temp,int count);
int writeout_pack_temp_0(int ts1,int ts2,int ts3,double *lw_buf_temp,int lw_count_temp);
int writeout_unpack_temp_0(int ts1,int ts2,int ts3,double *lw_recv_buf_temp,int displs_lw_temp);
void sigma_prod2_0(int ts1, int ts2, int ts3, int my_rank, int nprocs, int *receiver_list);
int writeout_pack_prod2_0(int ts1,int ts2,int ts3,double *lw_buf_prod2,int lw_count_prod2);
int writeout_unpack_prod2_0(int ts1,int ts2,int ts3,double *lw_recv_buf_prod2,int displs_lw_prod2);
void sigma_conv_0(int ts1, int ts2, int ts3, int my_rank, int nprocs, int *receiver_list);
int writeout_pack_conv_0(int ts1,int ts2,int ts3,double *lw_buf_conv,int lw_count_conv);
int writeout_unpack_conv_0(int ts1,int ts2,int ts3,double *lw_recv_buf_conv,int displs_lw_conv);
int pi_0( int t1, int t2, int t4, int nprocs);
void write_out(int my_rank,int nprocs,double *lw_buf_prod1,double *lw_recv_buf_prod1,int *displs_lw_prod1,double *lw_buf_temp,double *lw_recv_buf_temp,int *displs_lw_temp,double *lw_buf_prod2,double *lw_recv_buf_prod2,int *displs_lw_prod2,double *lw_buf_conv,double *lw_recv_buf_conv,int *displs_lw_conv);
