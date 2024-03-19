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
  PetscScalar    kr,ki, dot, spin;
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

  /*
    Enter the options for the number of holes and
    the total number of alpha electrons.
  */
  PetscInt  nholes;
  PetscInt  nalpha;
  PetscInt  DoS2;
  PetscReal t_inp = -1.0;
  PetscReal Jme_inp =  0.030;
  PetscReal Kme_inp = -0.0;
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-nh",&nholes,NULL));
  PetscCheck(nholes, PETSC_COMM_WORLD, PETSC_ERR_USER, "Must indicate the number of holes with the -nh option");
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-na",&nalpha,NULL));
  PetscCheck(nalpha, PETSC_COMM_WORLD, PETSC_ERR_USER, "Must indicate the number of alpha e- with the -na option");
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-ht",&t_inp,&flg));
  PetscCheck(flg, PETSC_COMM_WORLD, PETSC_ERR_USER, "Must indicate t with the -ht option");
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-hj",&Jme_inp,&flg));
  PetscCheck(flg, PETSC_COMM_WORLD, PETSC_ERR_USER, "Must indicate J with the -hj option");
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-hk",&Kme_inp,&flg));
  PetscCheck(flg, PETSC_COMM_WORLD, PETSC_ERR_USER, "Must indicate K with the -hk option");
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-hs2",&DoS2,NULL));
  PetscCheck(flg, PETSC_COMM_WORLD, PETSC_ERR_USER, "Must indicate whether or not to S2 with the -hs2 option (1=true, 0=false)");

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
  size_t nelecF1 = nsites - nholes;
  int natomax = 100;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Define Hamiltonian
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  //double Jme =  0.030;
  //double Kme = -0.8;
  //double t = -1.0;
  double Jme = Jme_inp;
  double Kme = Kme_inp;
  double t = t_inp;

  size_t sizeCFG  = binomialCoeff(nsites, (nsites-nholes));
  size_t sizeCSF  = binomialCoeff(nelec, nalpha);
  size_t sizeTotal  = sizeCFG*sizeCSF;

  size_t* configList = malloc(sizeCFG * sizeof(size_t));
  size_t* csfList    = malloc(sizeCSF * sizeof(size_t));
  // List for storing JK elems
  int max_nbrs = getMaxNeighbors(&graph, nsites);
  //printf(" Max Nbrs           # = %d \n",max_nbrs);
  size_t* csfJKList  = malloc(sizeCSF * max_nbrs * sizeof(size_t));

  generateConfigurations(nsites, (nsites-nholes), configList, &sizeCFG);

  // Sort the lists for binary search
  qsort(configList, sizeCFG, sizeof(size_t), compare);
  //printf(" Configuration List # = %ld \n",sizeCFG);
  //for( int i=0; i<sizeCFG; ++i ) {
  //  printBits(configList[i], nsites);
  //}

  generateConfigurations(nelec, nalpha, csfList, &sizeCSF);
  // Sort the lists for binary search
  qsort(csfList, sizeCSF, sizeof(size_t), compare);
  //printf(" CSF           List # = %ld \n",sizeCSF);
  //int ideter[nelec];
  //getdet(1, ideter, csfList, sizeCSF, nelec);
  //ideter[2]=2;
  //ideter[3]=1;
  //long int posA;
  //adr (ideter, &posA,
  //    csfList, sizeCFG, nelec);
  //printf(" pos=%ld\n",posA);
  //getdet(2, ideter, csfList, sizeCSF, nelec);
  //for( int i=0; i<nelec; ++i ) {
  //  printf(" %d > %d\n",i,ideter[i]);
  //}

  /* Note on ordering of determinants
   * =================================
   *
   *
   * Setup of the sites
   * -------------------
   *
   * Family 1 : --  --  --  --
   * Family 2 : --  --  --  --
   *
   * Holes are allowed only on the
   * 1st family of orbitals.
   * Configurations are generated only
   * for those orbitals i.e. nsites orbitals.
   *
   * The CSFs are generated for all the electrons
   * i.e. potentially 2*nsites.
   *
   * Ordering of CSFs
   * -----------------
   *
   * The ordering of the CSFs is as follows
   *
   * -----------------
   * | Fam 1 | Fam 2 |
   * -----------------
   * |1010100|0011001|
   * -----------------
   *
   * The first part corresponds to the Family 1
   * orbital set and the second part to
   * the Family 2 orbital set.
   *
   * This makes it easier to find the phase
   * change for one electron substitution.
   *
   * Ordering of alpha and beta
   * ---------------------------
   *
   * The alpha and beta electrons are
   * ordered as follows
   *
   * a1b1 a2b2 a3b3 ...
   *
   * where a1 is alpha for orbital 1
   * and b1 is beta electron for orbital 2
   *
   * This is different from
   * a1a2a3... b1b2b3...
   * where alpha and beta are separated.
   *
   *
    */

  // Declare a matrix of size 3 x 4
  int rows = sizeTotal;
  int cols = rows;
  //double** matrix = declare_matrix(rows, cols);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Compute the operator matrix that defines the eigensystem, Ax=kx
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscInt       n=sizeTotal,i,Istart,Iend,nev,maxit,its,nconv;
  PetscInt		   ncv, mpd;

  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"====================================================="));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\nDEHamv2: Double Exchange Eigenproblem, n=%" PetscInt_FMT "\n",n));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"=====================================================\n\n"));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," [Info] Max Nbrs   \t\t %" PetscInt_FMT "\n",(size_t)max_nbrs));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," [Info] N(configurations)   \t %" PetscInt_FMT "\n",sizeCFG));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," [Info] N(CSFs)   \t\t %" PetscInt_FMT "\n",sizeCSF));
  // Symmetric Matrix
  PetscCall(MatCreate(PETSC_COMM_WORLD,&A));
  PetscCall(MatCreateSBAIJ(PETSC_COMM_WORLD,1,PETSC_DECIDE,PETSC_DECIDE,n,n,max_nbrs*nsites,NULL,max_nbrs*nsites,NULL,&A));
  //PetscCall(MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n));
  PetscCall(MatSetType ( A, MATSBAIJ));
  PetscCall(MatMPIBAIJSetPreallocation(A,1,max_nbrs*nsites,NULL,max_nbrs*nsites,NULL));
  //PetscViewer viewer; // declare a viewer object
  //PetscViewerASCIIGetStdout(PETSC_COMM_WORLD, &viewer); // get the standard output
  /*
   * Matrix for the S2 operator
    */
  if(DoS2){
    // Symmetric Matrix
    PetscCall(MatCreate(PETSC_COMM_WORLD,&S2));
    PetscCall(MatCreateSBAIJ(PETSC_COMM_WORLD,1,PETSC_DECIDE,PETSC_DECIDE,n,n,nelec*nelec,NULL,nelec*nelec,NULL,&S2));
    PetscCall(MatSetType (S2, MATSBAIJ));
    PetscCall(MatMPIBAIJSetPreallocation(S2,1,nelec*nelec,NULL,nelec*nelec,NULL));
  }


  PetscCall(PetscTime(&tt1));


  /*
   * Initialize Hamiltonian
    */
  size_t icfg[1];
  size_t icsf[1];
  size_t cfgid, csfid;
  PetscCall(MatGetOwnershipRange(A,&Istart,&Iend));
  for (i=Istart;i<Iend;i++) {

    igraph_vector_int_t monoCFGList;
    igraph_vector_int_init(&monoCFGList, 0);
    igraph_vector_t monoMEs;
    igraph_vector_init(&monoMEs, 0);
    size_t posi = i;
    cfgid = findCFGID(posi, sizeCFG, sizeCSF);
    csfid = findCSFID(posi, sizeCFG, sizeCSF);
    icfg[0] = configList[cfgid];
    icsf[0] = csfList[csfid];
    //printBits(icfg[0], nsites);
    //printBitsDE(icsf[0], nelec, nelecF1);

    generateMonoCFGs(configList, sizeCFG, csfList, sizeCSF, &graph, posi, icfg[0], icsf[0], &monoCFGList, &monoMEs, t, Jme, Kme);
    //printf(" posi=%ld \n",i);
    for (int j = 0; j < igraph_vector_int_size(&monoCFGList); ++j) {
      PetscInt Jid = VECTOR(monoCFGList)[j];
      PetscReal val = VECTOR(monoMEs)[j];
      //PetscCall(MatSetValue(A,Jid,i,(PetscReal)val,ADD_VALUES));
      if( i > Jid ) PetscCall(MatSetValue(A,Jid,i,(PetscReal)val,INSERT_VALUES));
      else          PetscCall(MatSetValue(A,i,Jid,(PetscReal)val,INSERT_VALUES));
      //matrix[posi][Jid] += val;
      //printf("(%ld, %ld) = %f\n",posi, Jid, val);
    }

    igraph_vector_int_destroy(&monoCFGList);
    igraph_vector_destroy(&monoMEs);
  }
  PetscCall(PetscTime(&tt2));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Time used to build the matrix: %f\n",tt2-tt1));


  PetscCall(PetscTime(&tt1));
  PetscCall(MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY));
  PetscCall(PetscTime(&tt2));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Time used to assemble the matrix: %f\n",tt2-tt1));
  //PetscCall(MatView(A, viewer));

  if(DoS2) {
    /*
     * Initialize S2 operator
      */
    PetscCall(MatGetOwnershipRange(S2,&Istart,&Iend));
    for (i=Istart;i<Iend;i++) {

      igraph_vector_int_t monoCFGList;
      igraph_vector_int_init(&monoCFGList, 0);
      igraph_vector_t monoMEs;
      igraph_vector_init(&monoMEs, 0);
      size_t posi = i;
      cfgid = findCFGID(posi, sizeCFG, sizeCSF);
      csfid = findCSFID(posi, sizeCFG, sizeCSF);
      icfg[0] = configList[cfgid];
      icsf[0] = csfList[csfid];

      getS2Operator(csfid, &monoMEs, &monoCFGList, csfList, sizeCSF, &graph, nelec, nelec);
      for (int j = 0; j < igraph_vector_int_size(&monoCFGList); ++j) {
        PetscInt Jid = cfgid * sizeCSF + VECTOR(monoCFGList)[j];
        PetscInt Iid = i;
        PetscReal val = VECTOR(monoMEs)[j];
        //matrix[i][Jid] = VECTOR(MElist)[j];
        //printf(" %d %10.5f \n",Jid, VECTOR(MElist)[j]);
        if( Iid > Jid) PetscCall(MatSetValue(S2,Jid,Iid,val,INSERT_VALUES));
        else           PetscCall(MatSetValue(S2,Iid,Jid,val,INSERT_VALUES));
        //PetscCall(MatSetValue(S2,Jid,i,t*(PetscReal)VECTOR(MElist)[j],INSERT_VALUES));
        //printf("(%ld, %ld) = %f\n",Iid, Jid, val);
      }

      igraph_vector_int_destroy(&monoCFGList);
      igraph_vector_destroy(&monoMEs);
    }
    PetscCall(PetscTime(&tt2));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD," Time used to build the S2 operator: %f\n",tt2-tt1));


    PetscCall(PetscTime(&tt1));
    PetscCall(MatAssemblyBegin(S2,MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(S2,MAT_FINAL_ASSEMBLY));
    PetscCall(PetscTime(&tt2));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD," Time used to assemble the matrix: %f\n",tt2-tt1));
  }

  PetscCall(MatCreateVecs(A,NULL,&xr));
  PetscCall(MatCreateVecs(A,NULL,&xi));
  PetscCall(MatCreateVecs(A,NULL,&vs2));

  // Save file
  //save_matrix(matrix, rows, cols, "/tmp/4x4_de.csv");

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  /*
     Create eigensolver context
  */
  PetscCall(EPSCreate(PETSC_COMM_WORLD,&eps));

  /*
     Set operators. In this case, it is a standard eigenvalue problem
  */
  PetscCall(EPSSetOperators(eps,A,NULL));
  PetscCall(EPSSetProblemType(eps,EPS_HEP));
  PetscCall(EPSSetWhichEigenpairs(eps,EPS_SMALLEST_REAL));

  /*
     Set solver parameters at runtime
  */
  tol = 1.e-16;
  maxit = 10000000;
  PetscCall(EPSSetTolerances(eps,tol,maxit));
  //ncv  = 9;
  //mpd  = 10;
  //nev  = 4;
  //PetscCall(EPSSetDimensions(eps,nev,PETSC_DECIDE,PETSC_DECIDE));
  PetscCall(EPSSetFromOptions(eps));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the eigensystem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(PetscTime(&t1));
  PetscCall(EPSSolve(eps));
  PetscCall(PetscTime(&t2));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Time used to Solve EVP: %f\n",t2-t1));

  /*
     Optional: Get some information from the solver and display it
  */
  PetscCall(EPSGetIterationNumber(eps,&its));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Number of iterations of the method: %" PetscInt_FMT "\n",its));
  PetscCall(EPSGetType(eps,&type));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n\n",type));
  PetscCall(EPSGetDimensions(eps,&nev,NULL,NULL));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Number of requested eigenvalues: %" PetscInt_FMT "\n",nev));
  PetscCall(EPSGetTolerances(eps,&tol,&maxit));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Stopping condition: tol=%.4g, maxit=%" PetscInt_FMT "\n",(double)tol,maxit));


  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  /*
     Get number of converged approximate eigenpairs
  */
  nconv = 0;
  PetscCall(EPSGetConverged(eps,&nconv));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Number of converged eigenpairs: %" PetscInt_FMT "\n\n",nconv));

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

      if(DoS2) {
        /*
         * Get Spin S2 value
         */
        PetscCall(MatMult(S2, xr, vs2));
        PetscCall(VecDot(xr, vs2, &dot));

        spin = solveQuad(1.0, 1.0, -1.0*dot);
      }
      else spin = 100;

      if (im!=0.0) PetscCall(PetscPrintf(PETSC_COMM_WORLD," %9f%+9fi %12g\n",(double)re,(double)im,(double)error));
      else PetscCall(PetscPrintf(PETSC_COMM_WORLD,"   %12f       %12g       %12f\n",(double)re,(double)error,(double)fabs(spin)));
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
  PetscCall(EPSDestroy(&eps));
  PetscCall(MatDestroy(&A));
  PetscCall(VecDestroy(&xr));
  PetscCall(VecDestroy(&xi));
  PetscCall(VecDestroy(&vs2));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"===========================================\n"));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"HubHam: Success !!!                      \n"));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"===========================================\n"));
  PetscCall(SlepcFinalize());
  return 0;
}
