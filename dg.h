/******************************************************************/
// Define macros
#define N1D 8 //NbSideElem // number of elements in each dimension
#define InterpolationOrder 0 //OrderOfApprox // order of approximating polynomial
#define NbElement2D N1D*N1D // number of elements in 2D
#define NbQuadrPt1D InterpolationOrder+1 // number of quadrature points in 1D
#define NbQuadrPt2D (InterpolationOrder+1)*(InterpolationOrder+1) // number of quadrature points in 2D
#define NbDoF NbElement2D*NbQuadrPt2D // number of degree of freedom

/******************************************************************/
// Declare global variables
double Vert[(N1D+1)*(N1D+1)][2]; //each 1x2 row represents x,y coordinates of a vertex/node
int E2N[NbElement2D][4]; //element-to-node matrix
int I2E[2*NbElement2D][4]; //internalEdge-to-element matrix
double normals[2*NbElement2D][2]; //normal vectors
double h; //element side size; assuming square elements

double U[NbDoF] = {0.0}, R[NbDoF] = {0.0}, xyglob[NbDoF][2];

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
/******************************************************************/
// Declare functions
void buildmesh(double L, int N);
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

void calcresidual(double U[NbDoF], struct ResData *param);

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