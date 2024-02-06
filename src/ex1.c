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

  /*
    Enter the options for the number of holes and
    the total number of alpha electrons.
  */
  size_t nholes;
  size_t nalpha;
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-nh",&nholes,NULL));
  PetscCheck(nholes, PETSC_COMM_WORLD, PETSC_ERR_USER, "Must indicate the number of holes with the -nh option");
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-na",&nalpha,NULL));
  PetscCheck(nalpha, PETSC_COMM_WORLD, PETSC_ERR_USER, "Must indicate the number of alpha e- with the -na option");
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
  size_t nsites = num_vertices;
  size_t norb   = nsites*2; // Two orbitals on each site in the DE model
  size_t nelec = norb - nholes;
  size_t nbeta = nelec - nalpha;
  int natomax = 100;

  size_t sizeCFG  = binomialCoeff(nsites, (nsites-nholes));
  size_t sizeCSF  = binomialCoeff(nelec, nalpha);
  size_t sizeTotal  = sizeCFG*sizeCSF;

  size_t* configList = malloc(sizeCFG * sizeof(size_t));
  size_t* csfList    = malloc(sizeCSF * sizeof(size_t));

  generateConfigurations(nsites, (nsites-nholes), configList, &sizeCFG);

  // Sort the lists for binary search
  qsort(configList, sizeCFG, sizeof(size_t), compare);
  printf(" Configuration List # = %ld \n",sizeCFG);
  for( int i=0; i<sizeCFG; ++i ) {
    printBits(configList[i], nsites);
  }

  generateConfigurations(nelec, nalpha, csfList, &sizeCSF);

  // Sort the lists for binary search
  qsort(csfList, sizeCSF, sizeof(size_t), compare);
  printf(" CSF           List # = %ld \n",sizeCSF);
  //for( int i=0; i<sizeCSF; ++i ) {
  //  printBits(csfList[i], nelec);
  //}

  printf(" Ne=%ld Na=%ld Nb=%ld \n", nelec, nalpha, nbeta);

  // Get the global address of a CFGxCSF pair.
  size_t icfg[1];
  size_t icsf[1];
  size_t iglobalid;
  icfg[0] = 11;
  icsf[0] = 15 - 8 + 16;
  icfg[0] = configList[4];
  icsf[0] = csfList[22];
  size_t icfgid;
  size_t icsfid;

  findPositions(configList, sizeCFG, icfg, 1, &icfgid);
  findPositions(csfList, sizeCSF, icsf, 1, &icsfid);
  iglobalid = findGlobalID(icfgid, icsfid, sizeCFG, sizeCSF);
  printf(" cfgid = %ld csfid = %ld => glbid = %ld\n",icfgid, icsfid, iglobalid);
  icfgid = findCFGID(iglobalid, sizeCFG, sizeCSF);
  icsfid = findCSFID(iglobalid, sizeCFG, sizeCSF);
  printf(" cfgid = %ld csfid = %ld => glbid = %ld\n",icfgid, icsfid, iglobalid);

  igraph_vector_int_t monoCFGList;
  igraph_vector_int_init(&monoCFGList, 0);
  igraph_vector_t monoMEs;
  igraph_vector_init(&monoMEs, 0);
  generateMonoCFGs(configList, sizeCFG, csfList, sizeCSF, &graph, icfg[0], icsf[0], &monoCFGList, &monoMEs);
  printf(" Final ----- \n");
  for( int i=0; i<igraph_vector_int_size(&monoCFGList); ++i ) {
    printf(" %d >> %f \n",i, VECTOR(monoMEs)[i]);
    iglobalid = VECTOR(monoCFGList)[i];
    icfgid = findCFGID(iglobalid, sizeCFG, sizeCSF);
    icsfid = findCSFID(iglobalid, sizeCFG, sizeCSF);
    printBits(configList[icfgid], nsites);
    printBits(csfList[icsfid], nelec);
  }

  // Declare a matrix of size 3 x 4
  //int rows = sizeAlpha * sizeBeta;
  //int cols = rows;
  //double** matrix = declare_matrix(rows, cols);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Define Hamiltonian
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  double J =  0.01;
  double K =  6.0;
  double t = -1.0;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Compute the operator matrix that defines the eigensystem, Ax=kx
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscInt       n=sizeTotal,i,Istart,Iend,nev,maxit,its,nconv;
  PetscInt		   ncv, mpd;

  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"====================================================="));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\nDEHamv2: Double Exchange Eigenproblem, n=%" PetscInt_FMT "\n",n));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"=====================================================\n\n"));

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

  free(configList);

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
