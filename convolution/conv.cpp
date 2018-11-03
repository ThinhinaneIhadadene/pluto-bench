#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <unistd.h>
#include <ctime>
#include <cstdint>
#include "Halide.h"
#include "halide_image_io.h"

//constants
#define ROWS 5000
#define COLS 5000
#define BASE 255

//global vars
#pragma declarations
uint8_t src[ROWS][COLS][3];
float kernel[3][3];
uint8_t conv[ROWS][COLS][3];
#pragma enddeclarations

void init_array()
{
    int i, j, k;
    srand(0);
    for (i=0; i<ROWS; i++) {
        for (j=0; j<COLS; j++) {
            for (k=0; k<3; k++) {
                src[i][j][k] = rand()%BASE;
            }
        }
    }
}

void init_kernel(){
  kernel[0][0]=0;kernel[1][0]=1.0f/5;kernel[2][0]=0;
  kernel[0][1]=1.0f/5;kernel[1][1]=1.0f/5;kernel[2][1]=0;
  kernel[0][2]=0;kernel[1][2]=1;kernel[2][2]=0;
}


void print_array()
{
  int i,j,k;
  for (i=0; i<ROWS; i++) {
      for (j=0; j<COLS; j++) {
          for (k=0; k<3; k++) {
              printf("%d,", conv[i][j][k]);
          }
          printf("\t");
      }
      printf("\n");
  }
}

int main()
{
	 int q,w,cc,kq,kw; int rows;int cols;
   double t_start, t_end;
   float prod;
   using namespace Halide;
	//init_array() ;
    Buffer<uint8_t> input = Halide::Tools::load_image("../rgb.png");
    rows=input.height();
    cols=input.width();
    for(int i=0;i<rows;i++)
      for(int j=0;j<cols;j++)
        for(int k=0;k<3;k++)
          src[i][j][k]=input(j,i,k);

  init_kernel();


  for (q = 0; q < rows- 2; q++) {
      for (w = 0; w < cols - 2; w++) {
          for (cc = 0; cc < 3; cc++) {

              conv[q][w][cc] = src[q][w][cc]*kernel[0][0]+src[q+1][w][cc]*kernel[1][0]+
              src[q+2][w][cc]*kernel[2][0]+src[q][w+1][cc]*kernel[0][1]+
              src[q+1][w+1][cc]*kernel[1][1]+
              src[q+1][w+2][cc]*kernel[1][2]+src[q+2][w][cc]*kernel[2][0]+
              src[q+2][w+1][cc]*kernel[2][1]+src[q+2][w+2][cc]*kernel[2][2];
          }
      }
  }

  Buffer<uint8_t> output(cols,rows,3);
  for(int i=0;i<rows;i++)
    for(int j=0;j<cols;j++)
      for(int k=0;k<3;k++)
        output(j,i,k)=conv[i][j][k];
  Tools::save_image(output,"conv.png");

  return 0;
}
