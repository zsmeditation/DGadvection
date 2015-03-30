#ifndef _DG_H_FCN
#define _DG_H_FCN

/******************************************************************/
// Declare functions
void buildmesh(double L, int N, double h, double*** Vert, int*** E2N,
               int*** I2E, double*** normals);
void quad2d(int nq1, double a, double b, 
	        struct QuadrPt *quadrature);
void lgwt(int nq1, double a, double b, 
	      struct QuadrPt1D *quadrature1d);
void initialize(int N, int nq2, double xyn[NbQuadrPt2D][2], double h, double L);
double exactsol(double L, double xyglob[2], double t);

double lagrange(double xq[NbQuadrPt1D], int j, int p, double x);
void basis1d(double xq[NbQuadrPt1D], int p, double x, 
	         struct Db1byN *Phi);
void basis(double xq[NbQuadrPt1D], int p, double xy[2], 
	       struct Db1byNN *Phi);

void calcresidual(double U[NbDoF], struct ResData *param,
                  int*** I2E, double*** normals);

double glagrange(double xq[NbQuadrPt1D], int j, int p, double x);
void gbasis1d(double xq[NbQuadrPt1D], int p, double x, 
	          struct Db1byN *GPhi);
void gbasis(double xq[NbQuadrPt1D], double xyq[NbQuadrPt2D][2], int p, 
	        struct DbGradNNbyNN *GPhi);
void RefEdge2Elem(int e, double s[NbQuadrPt1D], int nq1, 
	              struct DbNby2 *xy);
void basisEdge_new(int e, int dir, double xq[NbQuadrPt1D], int p,
	               struct EdgePhi *ePhi);
void negmassinvmult(double invM[NbQuadrPt2D][NbQuadrPt2D],double R[NbDoF], double F[]);

int matrixAllocateInt(int*** Mptr, int row, int col);
int matrixAllocateDouble(double*** Mptr, int row, int col);

#endif // // _DG_H_FCN