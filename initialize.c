#include "MacVar.h"

/******************************************************************/
//Function definition: initialize
void initialize(int N, int nq2, double xyn[NbQuadrPt2D][2], double h, double L)
{
  int m = 0, i,j,iq, Ndof = N*N*nq2;

  for (j=0; j<N; j++)
  {
    for (i=0; i<N; i++)
    {
      for (iq=0; iq<nq2; iq++)
      {
        m = j*N + i; // current element count
        xyglob[nq2*m + iq][0] = h * (i + xyn[iq][0]);
        xyglob[nq2*m + iq][1] = h * (j + xyn[iq][1]);
        U[nq2*m + iq] = exactsol(L, xyglob[nq2*m + iq], 0.0);
//printf("@%d: (%le,%le)\n",nq2*m + iq, xyglob[nq2*m + iq][0],xyglob[nq2*m + iq][1]);
//printf("@%d: U0 = %le\n",nq2*m + iq, U[nq2*m + iq]);
      }
    }
  }
}