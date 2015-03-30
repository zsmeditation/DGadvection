#include "MacVar.h"

/******************************************************************/
//Function definition: calcresidual
void calcresidual(double U[NbDoF], struct ResData *param,
                  int*** I2E, double*** normals)
{
  int m, i, j, q, Im, Imj, Imq, ie;
  double Fx[param->nq2], Fy[param->nq2], u[NbQuadrPt2D] = {0.0}; 
  struct DbGradNNbyNN GPhi;
  // Input param: PhiMat, U, xq, xyq, p, nq2, wq, h, V, wq1, nq1
  int nq1 = param->nq1, nq2 = param->nq2, p = param->p;
  double h = param->h, V[2] = {param->V[0], param->V[1]};

  int elemL, elemR, edgeL, edgeR;
  double n[2];
  int IL, IR;
  struct EdgePhi ePhiL, ePhiR;
  double uL, uR, UL[NbDoF], UR[NbDoF], Vdotn, Fhat[nq1];

  double r;

/*  //test: print out 1D quadrature nodes
  printf("1d quadrature nodes\n");
  printf("xq =:");
   for (j=0; j<nq1; j++)
   {
     printf(" %le,", param->xq[j]);
   }
   printf(";\n");*/

/*  //test: print out 2D quadrature nodes
  printf("2d quadrature nodes\n");
  for (i=0; i<nq1; i++)
  {
    printf("i=%d:", i);
    for (j=0; j<nq1; j++)
    {
      printf(" (%le,%le)", param->xyq[i*nq1+j][0],param->xyq[i*nq1+j][1]);
    }
    printf(";\n");
  }*/

/*//test: print out spatial order
printf("spatial order p = %d\n", p);*/

  // contributions from element interiors; surface integral
  gbasis(param->xq, param->xyq, p, &GPhi); // basis fcn gradient matrix
  // Note this is gradient in reference space and needs to be converted to physical space

/*  //test: print out gbasis fcns
  printf("gbasis fcns\n");
  for (i=0; i<nq2; i++)
  {
    printf("i=%d:",i);
    for (j=0; j<nq2; j++)
    {
      printf(" %le,", GPhi.x[i][j]);
    }
    printf(";\n");
  }*/

  for (m=0; m<NbElement2D; m++) // loop over elements
  {
    for (i=0; i<nq2; i++) // loop over unknowns
    { 
      // re-initialize u	
      for (q=0; q<nq2; q++)
      {
        u[q] = 0.0;
      }

      for (q=0; q<nq2; q++)
      {
    	  Imq = m*nq2 + q;
	      for (j=0; j<nq2; j++)
	      {
	        Imj = m*nq2 + j;
	        u[q] = u[q] + param->PhiMat[q][j] * U[Imj];  
	      }
	      Fx[q] = u[q]*V[0] * param->wq[q]; 
	      Fy[q] = u[q]*V[1] * param->wq[q];
      }

      // find Fx and Fy
      Im = m*nq2 + i;
	  // contributions to residual from element interiors
	  r = 0.0;
	  for (q=0; q<nq2; q++)
	  {
	    r = r - h*h*((GPhi.x[q][i]*Fx[q] + GPhi.y[q][i]*Fy[q]) / h);
	  }
	  R[Im] = r;

// 
// printf("elem %d, basis %d, R_interior = %le\n", m,i,R[Im]);
    }
  }
  
  // contributions from interior edges
  for (ie=0; ie<2*NbElement2D; ie++)
  {
    elemL = (*I2E)[ie][0];
    elemR = (*I2E)[ie][1];
    edgeL = (*I2E)[ie][2];
    edgeR = (*I2E)[ie][3];
// printf("ie=%d, elemL=%d, elemR=%d\n", ie, elemL, elemR);
    n[0] = (*normals)[ie][0]; n[1] = (*normals)[ie][1];
    basisEdge_new(edgeL, 1, param->xq, p, &ePhiL); // counterclockwise
    basisEdge_new(edgeR, 0, param->xq, p, &ePhiR); // clockwise

//good//test
/*for (i=0; i<nq1; i++)
{
	printf("PhiL @ 1D node %d: ", i);
	for (j=0; j<nq2; j++)
	{
		printf("%le ", ePhiL.x[i][j]);
	}
	printf("\n\n");
}
for (i=0; i<nq1; i++)
{
	printf("PhiR @ 1D node %d: ", i);
	for (j=0; j<nq2; j++)
	{
		printf("%le ", ePhiR.x[i][j]);
	}
	printf("\n\n");
}*/
    
    for (i=0; i<nq2; i++)
    {
      
      for (q=0; q<nq1; q++) // loop over nodes on an edge (1D)
      {
      	uL = 0.0; uR = 0.0;
        for (j=0; j<nq2; j++) // loop over different basis fcns
        {
          IL = (elemL-1)*nq2 + j;
          IR = (elemR-1)*nq2 + j;
          uL = uL + ePhiL.x[q][j]*U[IL];
          uR = uR + ePhiR.x[q][j]*U[IR];
        }
        // calculate upwind flux
        Vdotn = V[0]*n[0] + V[1]*n[1];
        if (Vdotn > 0)
        {
          Fhat[q] = uL*Vdotn * param->wq1[q];
        }
        else
        {
          Fhat[q] = uR*Vdotn * param->wq1[q];
        }
      }

      // update residuals  
      for (q=0; q<nq1; q++)  // loop over different basis fcns
      {
        IL = (elemL-1)*nq2 + i;
        IR = (elemR-1)*nq2 + i;
        R[IL] = R[IL] + h*ePhiL.x[q][i]*Fhat[q];
        R[IR] = R[IR] - h*ePhiR.x[q][i]*Fhat[q];
// printf("@ elem %d, basis %d, R_IL = %le\n", m,i,R[IL]);
// printf("@ elem %d, basis %d, R_IR = %le\n", m,i,R[IR]);
      }
    }
  }
  return;
}