/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Double Exchange Hamiltonian - A high-performance Exact diagonalization program

   Copyright (c) 2024-, Vijay Gopal CHILKURI
   LICENCE : GNU GPL V2
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Standard Double Exchange eigenproblem corresponding to the topology given in the file.\n\n"
  "The command line options are:\n"
  "  -f <file>, where <file> = name of the graphml file.\n\n";

#include <slepceps.h>

#include "doubleexchange.h"
#include "readgraphmllib.h"

int main(int argc,char **argv)
{
  Mat            A;           /* problem matrix */
  Mat            S2;          /* S2 oper matrix */
  EPS            eps;         /* eigenproblem solver context */
  EPSType        type;
  PetscReal      error,tol,re,im;
  PetscScalar    kr,ki, dot;
  Vec            xr,xi, vs2;
  PetscLogDouble t1,t2,tt1,tt2;
  PetscReal normfin;
  PetscReal xymatfin = 0.0;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,(char*)0,help));

  char        graphmlFileName[PETSC_MAX_PATH_LEN]; /* input file name */
  PetscBool       flg;
  PetscCall(PetscOptionsGetString(NULL, NULL, "-f", graphmlFileName, sizeof(graphmlFileName), &flg));
  PetscCheck(flg, PETSC_COMM_WORLD, PETSC_ERR_USER, "Must indicate graphml file with the -f option");
  FILE* graphmlFile = fopen(graphmlFileName, "r");

  if (graphmlFile == NULL) {
    fprintf(stderr, "Error opening the file.\n");
    return 1;  // Return an error code or use another error handling method.
  }

  igraph_t graph;
  igraph_empty(&graph, 0, IGRAPH_DIRECTED);
  igraph_integer_t num_vertices;

  if (readGraphMLFile(graphmlFile, &graph)) {
    // Successfully read the graph, now you can work with 'graph'.
    num_vertices = igraph_vcount(&graph);
  }

  // Assume configAlpha and configBeta are sorted lists of all possible alpha and beta configurations
  size_t norb = num_vertices;
  size_t nalpha = norb/2;
  size_t nbeta = norb/2;
  int natomax = 100;

  size_t sizeAlpha = binomialCoeff(norb, nalpha);
  size_t sizeBeta = binomialCoeff(norb, nbeta);

  size_t* configAlpha = malloc(sizeAlpha * sizeof(size_t));
  size_t* configBeta = malloc(sizeBeta * sizeof(size_t));

  generateConfigurations(norb, nalpha, configAlpha, &sizeAlpha);
  generateConfigurations(norb, nbeta, configBeta, &sizeBeta);

  // Sort the lists for binary search
  qsort(configAlpha, sizeAlpha, sizeof(size_t), compare);
  qsort(configBeta, sizeBeta, sizeof(size_t), compare);

  // Declare a matrix of size 3 x 4
  int rows = sizeAlpha * sizeBeta;
  int cols = rows;
  //double** matrix = declare_matrix(rows, cols);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Define Hamiltonian
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  double J =  1.0;
  double t = -1.0;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Compute the operator matrix that defines the eigensystem, Ax=kx
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscInt       n=sizeAlpha*sizeBeta,i,Istart,Iend,nev,maxit,its,nconv;
  PetscInt		   ncv, mpd;

  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"====================================================="));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\nDEHamv2: Double Exchange Eigenproblem, n=%" PetscInt_FMT "\n",n));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"=====================================================\n\n"));
  //PetscCall(MatCreate(PETSC_COMM_WORLD,&A));
  //PetscCall(MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n));
  //PetscCall(MatCreateAIJ(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,n,n,3*num_vertices,NULL,3*num_vertices,NULL,&A));
  //PetscCall(MatMPIAIJSetPreallocation(A,3*num_vertices,NULL,3*num_vertices,NULL));
  ////PetscCall(MatSetFromOptions(A));
  ////PetscCall(MatSetUp(A));
  ///*
  // * Matrix for the S2 operator
  //  */
  //PetscCall(MatCreate(PETSC_COMM_WORLD,&S2));
  //PetscCall(MatSetSizes(S2,PETSC_DECIDE,PETSC_DECIDE,n,n));
  //PetscCall(MatCreateAIJ(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,n,n,4*num_vertices,NULL,4*num_vertices,NULL,&S2));
  //PetscCall(MatMPIAIJSetPreallocation(S2,4*num_vertices,NULL,4*num_vertices,NULL));

  //PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\n Number of vertices: %ld",(long)num_vertices));
  //PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\n Number of configurations alpha = %ld, beta = %ld\n",sizeAlpha,sizeBeta));
  //PetscCall(PetscTime(&tt1));

  ///*
  // * Initialize Hamiltonian
  //  */
  //PetscCall(MatGetOwnershipRange(A,&Istart,&Iend));
  //for (i=Istart;i<Iend;i++) {

  //  igraph_vector_t MElist;
  //  igraph_vector_init(&MElist, 0);
  //  igraph_vector_t Jdetlist;
  //  igraph_vector_init(&Jdetlist, 0);

  //  int diag = getHubbardDiag(i, configAlpha, sizeAlpha, configBeta, sizeBeta);
  //  PetscCall(MatSetValue(A,i,i,(double)diag,INSERT_VALUES));
  //  //matrix[i][i] = (double)diag*U;
  //  getAllHubbardMEs(i, &MElist, &Jdetlist, configAlpha, sizeAlpha, configBeta, sizeBeta, &graph);
  //  for (int j = 0; j < igraph_vector_size(&Jdetlist); ++j) {
  //    int Jid = VECTOR(Jdetlist)[j];
  //    //matrix[i][Jid] = t*VECTOR(MElist)[j];
  //    PetscCall(MatSetValue(A,i,Jid,t*(double)VECTOR(MElist)[j],INSERT_VALUES));
  //  }

  //  igraph_vector_destroy(&MElist);
  //  igraph_vector_destroy(&Jdetlist);
  //}
  //PetscCall(PetscTime(&tt2));
  //PetscCall(PetscPrintf(PETSC_COMM_WORLD," Time used to build the matrix: %f\n",tt2-tt1));


  //PetscCall(PetscTime(&tt1));
  //PetscCall(MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY));
  //PetscCall(MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY));
  //PetscCall(PetscTime(&tt2));
  //PetscCall(PetscPrintf(PETSC_COMM_WORLD," Time used to assemble the matrix: %f\n",tt2-tt1));

  ///*
  // * Initialize S2 operator
  //  */
  //PetscCall(MatGetOwnershipRange(S2,&Istart,&Iend));
  //for (i=Istart;i<Iend;i++) {

  //  igraph_vector_t MElist;
  //  igraph_vector_init(&MElist, 0);
  //  igraph_vector_t Jdetlist;
  //  igraph_vector_init(&Jdetlist, 0);

  //  getS2Operator(i, &MElist, &Jdetlist, configAlpha, sizeAlpha, configBeta, sizeBeta, &graph, num_vertices, natomax);
  //  for (int j = 0; j < igraph_vector_size(&Jdetlist); ++j) {
  //    int Jid = VECTOR(Jdetlist)[j];
  //    //matrix[i][Jid] = VECTOR(MElist)[j];
  //    //printf(" %d %10.5f \n",Jid, VECTOR(MElist)[j]);
  //    PetscCall(MatSetValue(S2,i,Jid,t*(double)VECTOR(MElist)[j],INSERT_VALUES));
  //  }

  //  igraph_vector_destroy(&MElist);
  //  igraph_vector_destroy(&Jdetlist);
  //}
  //PetscCall(PetscTime(&tt2));
  //PetscCall(PetscPrintf(PETSC_COMM_WORLD," Time used to build the S2 operator: %f\n",tt2-tt1));


  //PetscCall(PetscTime(&tt1));
  //PetscCall(MatAssemblyBegin(S2,MAT_FINAL_ASSEMBLY));
  //PetscCall(MatAssemblyEnd(S2,MAT_FINAL_ASSEMBLY));
  //PetscCall(PetscTime(&tt2));
  //PetscCall(PetscPrintf(PETSC_COMM_WORLD," Time used to assemble the matrix: %f\n",tt2-tt1));

  //PetscCall(MatCreateVecs(A,NULL,&xr));
  //PetscCall(MatCreateVecs(A,NULL,&xi));
  //PetscCall(MatCreateVecs(A,NULL,&vs2));

  //// Save file
  ////save_matrix(matrix, rows, cols, "/tmp/benzene_c.csv");

  ///* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //              Create the eigensolver and set various options
  //   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ///*
  //   Create eigensolver context
  //*/
  //PetscCall(EPSCreate(PETSC_COMM_WORLD,&eps));

  ///*
  //   Set operators. In this case, it is a standard eigenvalue problem
  //*/
  //PetscCall(EPSSetOperators(eps,A,NULL));
  //PetscCall(EPSSetProblemType(eps,EPS_HEP));
  //PetscCall(EPSSetWhichEigenpairs(eps,EPS_SMALLEST_REAL));

  ///*
  //   Set solver parameters at runtime
  //*/
  //tol = 1.e-9;
  //maxit = 10000000;
  //PetscCall(EPSSetTolerances(eps,tol,maxit));
  ////ncv  = 9;
  ////mpd  = 10;
  ////nev  = 4;
  ////PetscCall(EPSSetDimensions(eps,nev,PETSC_DECIDE,PETSC_DECIDE));
  //PetscCall(EPSSetFromOptions(eps));

  ///* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //                    Solve the eigensystem
  //   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  //PetscCall(PetscTime(&t1));
  //PetscCall(EPSSolve(eps));
  //PetscCall(PetscTime(&t2));
  //PetscCall(PetscPrintf(PETSC_COMM_WORLD," Time used to Solve EVP: %f\n",t2-t1));

  ///*
  //   Optional: Get some information from the solver and display it
  //*/
  //PetscCall(EPSGetIterationNumber(eps,&its));
  //PetscCall(PetscPrintf(PETSC_COMM_WORLD," Number of iterations of the method: %" PetscInt_FMT "\n",its));
  //PetscCall(EPSGetType(eps,&type));
  //PetscCall(PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n\n",type));
  //PetscCall(EPSGetDimensions(eps,&nev,NULL,NULL));
  //PetscCall(PetscPrintf(PETSC_COMM_WORLD," Number of requested eigenvalues: %" PetscInt_FMT "\n",nev));
  //PetscCall(EPSGetTolerances(eps,&tol,&maxit));
  //PetscCall(PetscPrintf(PETSC_COMM_WORLD," Stopping condition: tol=%.4g, maxit=%" PetscInt_FMT "\n",(double)tol,maxit));


  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  /*
     Get number of converged approximate eigenpairs
  */
  nconv = 0;
  //PetscCall(EPSGetConverged(eps,&nconv));
  //PetscCall(PetscPrintf(PETSC_COMM_WORLD," Number of converged eigenpairs: %" PetscInt_FMT "\n\n",nconv));

  if (nconv>0) {
    /*
       Display eigenvalues and relative errors
    */
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,
         "           k          ||Ax-kx||/||kx||         S2    \n"
         "   ----------------- ---------------------------------\n"));

    for (i=0;i<nconv;i++) {
      /*
        Get converged eigenpairs: i-th eigenvalue is stored in kr (real part) and
        ki (imaginary part)
      */
      PetscCall(EPSGetEigenpair(eps,i,&kr,&ki,xr,xi));
      /*
         Compute the relative error associated to each eigenpair
      */
      PetscCall(EPSComputeError(eps,i,EPS_ERROR_RELATIVE,&error));

#if defined(PETSC_USE_COMPLEX)
      re = PetscRealPart(kr);
      im = PetscImaginaryPart(kr);
#else
      re = kr;
      im = ki;
#endif

      /*
       * Get Spin S2 value
       */
      PetscCall(MatMult(S2, xr, vs2));
      PetscCall(VecDot(xr, vs2, &dot));

      if (im!=0.0) PetscCall(PetscPrintf(PETSC_COMM_WORLD," %9f%+9fi %12g\n",(double)re,(double)im,(double)error));
      else PetscCall(PetscPrintf(PETSC_COMM_WORLD,"   %12f       %12g       %12f\n",(double)re,(double)error,(double)abs(dot)));
    }
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\n"));
  }

  // Don't forget to destroy the graph when you're done.
  igraph_destroy(&graph);

  // Close the file when you're done with it.
  fclose(graphmlFile);

  free(configAlpha);
  free(configBeta);

  // Free the memory allocated for the matrix
  //for (int i = 0; i < rows; i++) {
  //  free(matrix[i]);
  //}
  //free(matrix);

  /*
     Free work space
  */
  //PetscCall(EPSDestroy(&eps));
  //PetscCall(MatDestroy(&A));
  //PetscCall(VecDestroy(&xr));
  //PetscCall(VecDestroy(&xi));
  //PetscCall(VecDestroy(&vs2));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"===========================================\n"));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"HubHam: Success !!!                      \n"));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"===========================================\n"));
  PetscCall(SlepcFinalize());
  return 0;
}
