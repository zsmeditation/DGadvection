#include "MacVar.h"

/******************************************************************/
//Function definition: negmassinvmult
void negmassinvmult(double*** invM, double R[NbDoF], double F[])
{
  int m, Im, q;
  for (m=0; m<NbElement2D; m++) // loop over elements
  {
    for (q=0; q<NbQuadrPt2D; q++) // loop over basis fcns
    {
      Im = m*NbQuadrPt2D + q;
      F[Im] = - (*invM)[q][q] * R[Im]; // assume M is diagonal
    }
  }
}