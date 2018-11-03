#define ROWS 5000
#define COLS 5000
#define BASE 255
#define CV_DESCALE(x, n) (((x) + (1 << ((n)-1))) >> (n))
uint8_t src[ROWS][COLS][3];
uint8_t dst[ROWS][COLS];
enum {
    yuv_shift  = 14,
    R2Y        = 4899,
    G2Y        = 9617,
    B2Y        = 1868,
};
