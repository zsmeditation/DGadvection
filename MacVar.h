#ifndef _DG_H
#define _DG_H

#include "config.h"

/******************************************************************/
// Definitions using macros from config.h
#define NbElement2D N1D*N1D // number of elements in 2D
#define NbQuadrPt1D (InterpolationOrder+1) // number of quadrature points in 1D
#define NbQuadrPt2D NbQuadrPt1D*NbQuadrPt1D // number of quadrature points in 2D
#define NbDoF NbElement2D*NbQuadrPt2D // number of degree of freedom

/******************************************************************/
// Declare global variables
double U[NbDoF], R[NbDoF], xyglob[NbDoF][2];

struct QuadrPt1D {int nq1; double xq[NbQuadrPt1D], wq1[NbQuadrPt1D];};
struct QuadrPt {int nq2; double xq[NbQuadrPt1D], wq1[NbQuadrPt1D], xyq[NbQuadrPt2D][2], wq[NbQuadrPt2D];};
struct Db1byN {double x[NbQuadrPt1D];};
struct Db1byNN {double x[NbQuadrPt2D];};
struct DbNNbyNN {double x[NbQuadrPt2D][NbQuadrPt2D];};

struct DbGradNNbyNN {double x[NbQuadrPt2D][NbQuadrPt2D], y[NbQuadrPt2D][NbQuadrPt2D];};
struct DbNby2 {double xy[NbQuadrPt1D][2];};

struct EdgePhi {double x[NbQuadrPt1D][NbQuadrPt2D];};
struct ResData
{
  int nq1, nq2, p;
  double xq[NbQuadrPt1D], wq1[NbQuadrPt1D];
  double xyq[NbQuadrPt2D][2], wq[NbQuadrPt2D];
  double h, V[2];
  double PhiMat[NbQuadrPt2D][NbQuadrPt2D];
};

#endif // _DG_H
