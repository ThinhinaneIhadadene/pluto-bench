#include <omp.h>

#define S1(n_r,n_c)	o_r = (((a11 * ( n_r)) + (a10 * ( n_c))) + b00);
#define S2(n_r,n_c)	o_c = (((a01 * ( n_r)) + (a00 * ( n_c))) + b10);
#define S3(n_r,n_c)	r = (o_r - floorf(o_r));
#define S4(n_r,n_c)	c = (o_c - floorf(o_c));
#define S5(n_r,n_c)	coord_00_r = floorf(o_r);
#define S6(n_r,n_c)	coord_00_c = floorf(o_c);
#define S7(n_r,n_c)	coord_01_r = coord_00_r;
#define S8(n_r,n_c)	coord_01_c = (coord_00_c + (1));
#define S9(n_r,n_c)	coord_10_r = (coord_00_r + (1));
#define S10(n_r,n_c)	coord_10_c = coord_00_c;
#define S11(n_r,n_c)	coord_11_r = (coord_00_r + (1));
#define S12(n_r,n_c)	coord_11_c = (coord_00_c + (1));
#define S13(n_r,n_c)	coord_00_r = ((coord_00_r < (0)) ? (0) : ((coord_00_r > ((5000) - (1))) ? ((5000) - (1)) : coord_00_r));
#define S14(n_r,n_c)	coord_00_c = ((coord_00_c < (0)) ? (0) : ((coord_00_c > ((5000) - (1))) ? ((5000) - (1)) : coord_00_c));
#define S15(n_r,n_c)	coord_01_r = ((coord_01_r < (0)) ? (0) : ((coord_01_r > ((5000) - (1))) ? ((5000) - (1)) : coord_01_r));
#define S16(n_r,n_c)	coord_01_c = ((coord_01_c < (0)) ? (0) : ((coord_01_c > ((5000) - (1))) ? ((5000) - (1)) : coord_01_c));
#define S17(n_r,n_c)	coord_10_r = ((coord_10_r < (0)) ? (0) : ((coord_10_r > ((5000) - (1))) ? ((5000) - (1)) : coord_10_r));
#define S18(n_r,n_c)	coord_10_c = ((coord_10_c < (0)) ? (0) : ((coord_10_c > ((5000) - (1))) ? ((5000) - (1)) : coord_10_c));
#define S19(n_r,n_c)	coord_11_r = ((coord_11_r < (0)) ? (0) : ((coord_11_r > ((5000) - (1))) ? ((5000) - (1)) : coord_11_r));
#define S20(n_r,n_c)	coord_11_c = ((coord_11_c < (0)) ? (0) : ((coord_11_c > ((5000) - (1))) ? ((5000) - (1)) : coord_11_c));
#define S21(n_r,n_c)	A00 = src[coord_00_r][coord_00_c];
#define S22(n_r,n_c)	A10 = src[coord_10_r][coord_10_c];
#define S23(n_r,n_c)	A01 = src[coord_01_r][coord_01_c];
#define S24(n_r,n_c)	A11 = src[coord_11_r][coord_11_c];
#define S25(n_r,n_c)	dst[ n_r][ n_c] = ((A00 + ((A00 - A10) * r)) + (((A00 + ((A00 - A10) * r)) - (A01 + ((A01 - A11) * r))) * c));

		int t1, t2, t3;

		int lb, ub, lbd, ubd, lb2, ub2;
		register int lbv, ubv;

/* Start of CLooG code */
for (t1=0;t1<=4999;t1++) {
  for (t2=0;t2<=4999;t2++) {
    S2(t1,t2);
    S6(t1,t2);
    S12(t1,t2);
    S20(t1,t2);
    S10(t1,t2);
    S18(t1,t2);
    S8(t1,t2);
    S16(t1,t2);
    S14(t1,t2);
    S4(t1,t2);
    S1(t1,t2);
    S5(t1,t2);
    S11(t1,t2);
    S19(t1,t2);
    S24(t1,t2);
    S9(t1,t2);
    S17(t1,t2);
    S22(t1,t2);
    S7(t1,t2);
    S15(t1,t2);
    S23(t1,t2);
    S13(t1,t2);
    S21(t1,t2);
    S3(t1,t2);
    S25(t1,t2);
  }
}
/* End of CLooG code */
