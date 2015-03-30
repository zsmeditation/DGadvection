/******************************************************************/
//Function definition: buildmesh
void buildmesh(double L, int N, double h, double*** Vert, int*** E2N,
               int*** I2E, double*** normals)
{
  int k, ie, i, j, jj, kR, kU;
  double s[N+1];
  /* Note
   * * i, j correspond to row, column number.
   * * Counts are made row by row.
   */
  // 1D coordinates
  for (i=0; i<=N; i++)
  {
    s[i] = i*h; 
//good//printf("@%d: s = %le\n",s[i]);
  }
  // Vert
  k = 0; // node count
  for (i=0; i<=N; i++)
  {
    for (j=0; j<=N; j++)
    {
      (*Vert)[k][0] = s[j];
      (*Vert)[k][1] = s[i];
//good//printf("@%d: (%f, %f)\n;",k, (*Vert)[k][0],(*Vert)[k][1]);
      k = k+1;
    }
  }
  // E2N
  k = 0; // element count
  for (i=0; i<N; i++)
  {
    for (j=0; j<N; j++)
    {
      // Note that the element number starts from 1 and cannot be 0
      jj = i*(N+1) + (j+1);
      (*E2N)[k][0] = jj;
      (*E2N)[k][1] = jj+1;
      (*E2N)[k][2] = jj+N+2;
      (*E2N)[k][3] = jj+N+1;
//good//printf("@%d: (%d, %d, %d, %d)\n;",k, (*E2N)[k][0],(*E2N)[k][1],(*E2N)[k][2],(*E2N)[k][3]);
      k = k+1;
    }
  }
  // I2E
  k = 0; // element count
  ie = 0; // edge count
  for (i=0; i<N; i++)
  {
    for (j=0; j<N; j++)
    {
      k = k+1;
      kR = k+1; if(j==(N-1)){kR = kR-N;}
      kU = k+N; if(i==(N-1)){kU = k-N*(N-1);}
      (*I2E)[ie][0] = k;
      (*I2E)[ie][1] = kR;
      (*I2E)[ie][2] = 2;
      (*I2E)[ie][3] = 4;
      (*normals)[ie][0] = 1;
      (*normals)[ie][1] = 0;
// printf("@%d: (%d, %d), (%d, %d), normal(	x) = %1.1le\n;",ie, (*I2E)[ie][0],(*I2E)[ie][1],(*I2E)[ie][2],(*I2E)[ie][3],(*normals)[ie][0]);
      ie = ie+1;
      (*I2E)[ie][0] = k;
      (*I2E)[ie][1] = kU;
      (*I2E)[ie][2] = 3;
      (*I2E)[ie][3] = 1;
      (*normals)[ie][0] = 0;
      (*normals)[ie][1] = 1;
// printf("@%d: (%d, %d), (%d, %d), normal(x) = %1.1le\n;",ie, (*I2E)[ie][0],(*I2E)[ie][1],(*I2E)[ie][2],(*I2E)[ie][3],(*normals)[ie][0]);
      ie = ie+1;
    }
  }
  return;
}