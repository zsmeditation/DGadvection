// Include default library
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
// Include own header file
#include "MacVar.h"
#include "Fcn.h"
/******************************************************************/
// Main function
int main(int argc, char *argv[])
{
  /********************/
  // Define variables
  double L = 1.0; //length of computational domain; square in x,y
  double T = 1.0; //extent of time domain
  double V[2] = {1.0, 0.0}; //V[0] and V[1] are x,y components of advection velocity
  int p = InterpolationOrder; //interpolation order of solution approximation in space
  int N = N1D;
  int Nelem = NbElement2D; //number of quad elements along 1D and of 2D
  int nq1 = NbQuadrPt1D, nq2 = NbQuadrPt2D; //number of coefficients per 1D/2D basis
  struct QuadrPt quadrature; // Legendre-Gauss quadrature points; tensor product in 2D
  
  double dt; int it; //time step size
  double Us[NbDoF], F0[NbDoF], F1[NbDoF], F2[NbDoF], F3[NbDoF];
  double /*PhiMat[nq2][nq2],*/ M[nq2][nq2], invM[nq2][nq2]; 
  struct Db1byNN Phi;
  int i, j;
  struct ResData resdata;

  double** Vert; matrixAllocateDouble(&Vert, (N1D+1)*(N1D+1), 2); //each 1x2 row represents x,y coordinates of a vertex/node
  int** E2N; matrixAllocateInt(&E2N, NbElement2D, 4);//element-to-node matrix
  int** I2E; matrixAllocateInt(&I2E, 2*NbElement2D, 4); //internalEdge-to-element matrix
  double** normals; matrixAllocateDouble(&normals, 2*NbElement2D, 2);//normal vectors
  double h = L/N; //element side size; assuming square elements
  /********************/
  // Build a quad mesh
  buildmesh(L, N, h, &Vert, &E2N, &I2E, &normals);
  /********************/ 
  // Construct basis
  quad2d(nq1, 0, 1, &quadrature);
  /********************/
  // Initialize solution
  initialize(N, nq2, quadrature.xyq, h, L);
  // Print out intial solution
/*  printf("Initial condition\n");
  for (i=0; i<Nelem; i++)
  {
      	for (j=0; j<nq2; j++)
      	{
      		printf("%2.18le\n", U[i*nq2+j]);
      	}
  }*/
  /********************/
  // calculate matrix PhiMat
  for (i=0; i<nq2; i++)
  {
    basis(quadrature.xq, p, quadrature.xyq[i], &Phi);
    for (j=0; j<nq2; j++)
    {
      resdata.PhiMat[i][j] = Phi.x[j];
/*//test: print out Phi matrix
printf("i,j = %d,%d: Phi.x[j] = %le\n", i,j,Phi.x[j]);*/
    }
  }
  // calculate mass matrix M and inverse
    // assume M is diagonal
  for (i=0; i<nq2; i++)
  {
    for (j=0; j<nq2; j++)
    {
      M[i][j] = 0; // fill matrix with zeroes
      invM[i][j] = 0; 
    }
    M[i][i] = h*h*quadrature.wq[i]; // fill matrix diagonal
    invM[i][i] = 1/M[i][i]; 
  }

/*//test: print out invM matrix
  for (i=0; i<nq2; i++)
  {
  	printf("invM @ row %d: ", i);
    for (j=0; j<nq2; j++)
    {
    	printf("%le ", invM[i][j]);
    }
    printf("\n");
  }*/

  /********************/
  // pack up data
  resdata.nq1 = nq1; resdata.nq2 = nq2;
  resdata.p = p; 
  for (i=0; i<nq1; i++)
  {
    resdata.xq[i] = quadrature.xq[i];
    resdata.wq1[i] = quadrature.wq1[i];
  }
  for (i=0; i<nq2; i++)
  {
    resdata.xyq[i][0] = quadrature.xyq[i][0];
    resdata.xyq[i][1] = quadrature.xyq[i][1];
    resdata.wq[i] = quadrature.wq[i];
	
	//test: print out quadrature points & weights
	// printf("Quadrature @%d: (%le, %le) with w = %le\n", i,resdata.xyq[i][0],resdata.xyq[i][1],resdata.wq[i]);
  }
  resdata.h = h; 
  resdata.V[0] = V[0]; resdata.V[1] = V[1];
  /********************/
  // Time integration using RK4
  dt = T/Nt;
  for (it=0; it<Nt; it++)
  {
  	  for (i=0; i<NbDoF; i++) {Us[i] = U[i];} // copy array
      calcresidual(Us, &resdata, &I2E, &normals); // update R
      negmassinvmult(invM, R, F0); // calculate F0

      for (i=0; i<NbDoF; i++) {Us[i] = U[i] + 0.5*dt*F0[i];} // copy array
      calcresidual(Us, &resdata, &I2E, &normals); // update R
      negmassinvmult(invM, R, F1); // calculate F1
      
      for (i=0; i<NbDoF; i++) {Us[i] = U[i] + 0.5*dt*F1[i];} // copy array
      calcresidual(Us, &resdata, &I2E, &normals); // update R
      negmassinvmult(invM, R, F2); // calculate F2
      
      for (i=0; i<NbDoF; i++) {Us[i] = U[i] + dt*F2[i];} // copy array
      calcresidual(Us, &resdata, &I2E, &normals); // update R
      negmassinvmult(invM, R, F3); // calculate F3

      for (i=0; i<NbDoF; i++)
      {
        U[i] = U[i] + (dt/6)*(F0[i] + 2*F1[i] + 2*F2[i] + F3[i]); // update U
      }
  }
  /********************/
/*  // Print out residual
  printf("Final residual\n");
  for (i=0; i<Nelem; i++)
  {
      	for (j=0; j<nq2; j++)
      	{
      		printf("%2.18le\n", R[i*nq2+j]);
      	}
  }*/
  // Print out solution
  // printf("Final solution\n");
  for (i=0; i<Nelem; i++)
  {
      	for (j=0; j<nq2; j++)
      	{
      		printf("%2.18le\n", U[i*nq2+j]);
      	}
  }
  /********************/
  return 0;
}

/******************************************************************/
//Function definition: quad2d
void quad2d(int nq1, double a, double b, 
	        struct QuadrPt *quadrature)
{
  struct QuadrPt1D quadrature1d;
  int i,j, k = 0, nq2 = nq1*nq1;
  double xyq[nq2][2], wq[nq2];

  lgwt(nq1, a, b, &quadrature1d);
  
  for (j=0; j<nq1; j++)
  {
    for (i=0; i<nq1; i++)
    {
      quadrature->xyq[k][0] = quadrature1d.xq[i];
      quadrature->xyq[k][1] = quadrature1d.xq[j];
      quadrature->wq[k] = quadrature1d.wq1[i]*quadrature1d.wq1[j];
//good//printf("Quadrature @ %d: (%le,%le), w = %le\n", k,quadrature->xyq[k][0],quadrature->xyq[k][1],quadrature->wq[k]);
      k = k+1;
    }
    quadrature->xq[j] = quadrature1d.xq[j];
    quadrature->wq1[j] = quadrature1d.wq1[j];
  }
  
  quadrature->nq2 = nq2;
}
/******************************************************************/
//Function definition: lgwt
void lgwt(int nq1, double a, double b, 
	      struct QuadrPt1D *quadrature1d)
{
  int N1 = nq1, N2 = nq1+1;
  double xu[N1], y[N1], x[N1], w[N1], y0[N1];
  double dx;
  int i, j, iter = 0;
  double maxres, tolerance = 2.2204E-16;
  double L[N1][N2], Lp[N1]; // Legendre-Gauss Vandermonde Matrix and its derivative
  
  if ( nq1 == 1 )
  {
    dx = 1-(-1);
  }
  else
  {
    dx = (1-(-1))/(N1-1);
  }

  // Initialize quadrature points (stored in y vector)
  for (i=0; i<N1; i++) 
  {
    xu[i] = -1 + dx*i;
    y[i] = cos((2*i+1)*M_PI/(2*N1)) + 0.27/N1*sin(M_PI*xu[i]*(N1-1)/N2); //initial guess
//good//printf("@%d: xui = %le, yi = %le\n",i,xu[i],y[i]);
  }
  
  do
  {
    maxres = 0.0;
    for (i=0; i<N1; i++)
    {
       L[i][0] = 1;
       L[i][1] = y[i];
//good//printf("@(%d,%d): L=%le\n",i,0,L[i][0]);
//good//printf("@(%d,%d): L=%le\n",i,1,L[i][1]);
       for (j=2; j<N2; j++)
       {
         L[i][j] = ( (2*j-1)*y[i]*L[i][j-1] - (j-1)*L[i][j-2] )/j;
//good//printf("@(%d,%d): L=%le\n",i,j,L[i][j]);
       }

       Lp[i] = N2 * (L[i][N1-1] - y[i]*L[i][N2-1]) / (1-y[i]*y[i]);
//good//printf("@%d, Lp=%le\n",i,Lp[i]);
       y0[i] = y[i]; // store y
       y[i] = y0[i] - L[i][N2-1]/Lp[i]; // update y
//good//printf("@%d, y0=%le, y=%le\n",i,y0[i],y[i]);      
       maxres = fmax(maxres, fabs(y[i]-y0[i]));
    }
//good//iter++;
//good//printf("iter %d, maxres = %le\n",iter,maxres);
  } while (maxres>tolerance);
 
  for (i=0; i<nq1; i++)
  {
    quadrature1d->xq[i] = 0.5*( a*(1-y[i]) + b*(1+y[i]) );
    quadrature1d->wq1[i] = (b-a) / ((1-pow(y[i],2.0))*pow(Lp[i],2.0)) * (pow(N2,2.0)/pow(N1,2.0)); 
//good//printf("@%d, x(i) = %le, w(i) = %le\n",i,x[i],w[i]);
  } 
  quadrature1d->nq1 = nq1;
}

/******************************************************************/
//Function definition: exactsol
double exactsol(double L, double xyglob[2], double t)
{
  double u;
  u = exp(-100*(pow(xyglob[0]/L - 0.5, 2) + pow(xyglob[1]/L - 0.5, 2)));
  return u;
}

/******************************************************************/
//Function definition: lagrange
double lagrange(double xq[NbQuadrPt1D], int j, int p, double x)
{
  // Output: scalar value of j_th 1D lagrange basis function at location x
  // j ranges from 1 to p+1
  double phi;
  int i;
  // zero-th order; constant
  if (p == 0) 
    {
      phi = 1.0;
      return phi;
    }
  // first order (linear) or higher
  phi = 1.0; // initialize
  for (i=0; i<=p; i++)
  {
    if (i != (j-1))
      {
        phi = phi * (x - xq[i]) / (xq[j-1] - xq[i]);
      }
  }
  return phi;
}
/******************************************************************/
//Function definition: basis1d
void basis1d(double xq[NbQuadrPt1D], int p, double x, 
	         struct Db1byN *Phi)
{
  // Output: vector of values of all 1D basis fcns of order p at location x 
  int i;
  for (i=0; i<=p; i++) // loop over different basis fcns of order p
  {
    Phi->x[i] = lagrange(xq, i+1, p, x); // lagrange basis
  }
}
/******************************************************************/
//Function definition: basis
void basis(double xq[NbQuadrPt1D], int p, double xy[2], 
	       struct Db1byNN *Phi)
{
  // Output: matrix of values of all 2D basis fcns of order p at location x 
  int i,j,k,count, nq1 = p+1, nq2 = nq1*nq1;
  struct Db1byN Phix, Phiy;
  // 1D basis fcn values
  basis1d(xq, p, xy[0], &Phix);
  basis1d(xq, p, xy[1], &Phiy);
//printf("Phix: %2.2le, %2.2le\n",Phix.x[0],Phix.x[1]);
//printf("Phiy: %2.2le, %2.2le\n\n",Phiy.x[0],Phiy.x[1]);
  // 2D basis fcn values by tensor product
  for (j=0; j<nq1; j++) // loop over different phiy/ changing row
  {
    for (k=0; k<nq1; k++) // loop over different phix/ changing column
    {
      count = j*nq1 + k;
      Phi->x[count] = Phix.x[k] * Phiy.x[j];
    }
  }
}

/******************************************************************/
//Function definition: glagrange
double glagrange(double xq[NbQuadrPt1D], int j, int p, double x)
{
  int i, k;
  double den = 1.0, num = 0.0, prod, gphi;
  
  // find denominator
  for (i=0; i<=p; i++)
  {
    if (i != (j-1))
    {
      den = den * (xq[j-1] - xq[i]);
    }
  }

  // find numerator
  if (p == 0) // constant p=0
  {
      gphi = 0.0;
      return gphi;
  }
  else if (p == 1) // linear p=1
  {
    num = 1.0;
  }
  else // quadratic or higher order p>=2
  {
    for (i=0; i<=p; i++) // loop over different 1D nodes
    {
      if (i != (j-1))
      {
        prod = 1.0;
        for (k=0; k<=p; k++)
        {
          if (( k!=(j-1) ) && ( k!=i )) 
          {
            prod = prod * (x - xq[k]);
          } 
        }
        num = num + prod;
      }
    }
  }

gphi = num/den;
return gphi;
}

/******************************************************************/
//Function definition: gbasis1d
void gbasis1d(double xq[NbQuadrPt1D], int p, double x, 
	          struct Db1byN *GPhi)
{
  int i;

  for (i=0; i<=p; i++)
  {
    GPhi->x[i] = glagrange(xq, i+1, p, x);
  }
}

/******************************************************************/
//Function definition: gbasis
void gbasis(double xq[NbQuadrPt1D], double xyq[NbQuadrPt2D][2], int p, 
	        struct DbGradNNbyNN *GPhi)
{
  int i,j,k,count, nq1 = p+1, nq2 = nq1*nq1;
  struct Db1byN Phi1x, Phi1y, GPhi1x, GPhi1y;
  for (i=0; i<nq2; i++) // loop over all nodes
  {
    // find 1D basis and gbasis fcns at i_th node
    basis1d(xq, p, xyq[i][0], &Phi1x);
    basis1d(xq, p, xyq[i][1], &Phi1y);
    gbasis1d(xq, p, xyq[i][0], &GPhi1x);
    gbasis1d(xq, p, xyq[i][1], &GPhi1y);

/*//test: print out basis1d & gbasis1d
printf("basis1d @ node %d:",i);
for (j=0; j<nq1; j++)
{
  printf(" %le", Phi1x.x[j]);
}
printf(";\n");
//
printf("gbasis1d @ node %d:",i);
for (j=0; j<nq1; j++)
{
  printf(" %le", GPhi1x.x[j]);
}
printf(";\n");*/

    for (j=0; j<nq1; j++) // loop over different phiy/ changing row
    {
      for (k=0; k<nq1; k++) // loop over different phix/ changing column
      {
         count = j*nq1 + k;
         GPhi->x[i][count] = GPhi1x.x[k] *  Phi1y.x[j];
         GPhi->y[i][count] =  Phi1x.x[k] * GPhi1y.x[j];
      }
    }
//good // test
/*    printf("@ node %d, GPhi.x: ", i);
    for (k=0; k<=count; k++)
    {
       printf("%le ", GPhi.x[i][k]);	
    }
    printf("\n"); 

    printf("@ node %d, GPhi.y: ", i);
    for (k=0; k<=count; k++)
    {
       printf("%le ", GPhi.y[i][k]);	
    }
    printf("\n");
*/
  }
}

/******************************************************************/
//Function definition: RefEdge2Elem
void RefEdge2Elem(int e, double s[NbQuadrPt1D], int nq1, 
	              struct DbNby2 *xy)
{
  int i;
  switch(e)
  {
    case 1 :
      for (i=0; i<nq1; i++)
      {
        xy->xy[i][0] = s[i];
        xy->xy[i][1] = 0.0;
      }
      break;
    case 2 :
      for (i=0; i<nq1; i++)
      {
        xy->xy[i][0] = 1.0;
        xy->xy[i][1] = s[i];
      }
      break;
    case 3 :
      for (i=0; i<nq1; i++)
      {
        xy->xy[i][0] = 1.0 - s[i];
        xy->xy[i][1] = 1.0;
      }
      break;
    case 4 :
      for (i=0; i<nq1; i++)
      {
        xy->xy[i][0] = 0;
        xy->xy[i][1] = 1.0 - s[i];
      }
      break;
  }
}

/******************************************************************/
//Function definition: basisEdge_new
void basisEdge_new(int e, int dir, double xq[NbQuadrPt1D], int p,
	               struct EdgePhi *ePhi)
{
  int i, nq1=p+1, j,k,count;
  double xy[nq1];
  struct DbNby2 xyRef;
  struct Db1byNN Phi;
  
  for (i=0; i<nq1; i++) 
  {
    if (dir == 0) // clockwise: dir=0
    {
      xy[i] = 1.0 - xq[i];  
    }
    else if (dir == 1) // clockwise: dir=1
    {
      xy[i] = xq[i];  
    }
    else // dir is invalid input
    {
      return;
    }
  }

  // coordinates of edge nodes
  RefEdge2Elem(e, xy, nq1, &xyRef); 
   
  for (i=0; i<nq1; i++) // loop over edge nodes
  {
    basis(xq, p, xyRef.xy[i], &Phi);
    count = 0;
    for (j=0; j<nq1; j++) // loop over different phiy/ changing row
    {
      for (k=0; k<nq1; k++) // loop over different phix/ changing column
      {
         count = j*nq1 + k;
         ePhi->x[i][count] = Phi.x[count];
      }
    }
  }
}

