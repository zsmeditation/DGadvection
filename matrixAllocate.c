#include <stdlib.h>

void matrixAllocateInt(int*** Mptr, int row, int col)
{
  *Mptr = (int**)malloc(row*sizeof(int*));
  int i=0;
  for (i=0; i<row; i++)
  (*Mptr)[i] = malloc(col*sizeof(int));
}


void matrixAllocateDouble(double*** Mptr, int row, int col)
{
  *Mptr = (double**)malloc(row*sizeof(double*));
  int i=0;
  for (i=0; i<row; i++)
  (*Mptr)[i] = malloc(col*sizeof(double));
}