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
#include "get_proj_general.h"
#include "utils.h"

int main(int argc,char **argv)
{
  Mat            A;           /* problem matrix */
  Mat            S2;          /* S2 oper matrix */
  Mat            tpsx;        /* The vector that stores the TPS operator */
  EPS            eps;         /* eigenproblem solver context */
  EPSType        type;
  PetscReal      error,tol,re,im;
  PetscScalar    kr,ki, dot, spin;
  Vec            xr,xi, vs2;
  PetscLogDouble t1,t2,tt1,tt2;
  PetscReal normfin;
  PetscReal xymatfin = 0.0;
  PetscMPIInt    rank, size;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,(char*)0,help));

  char        graphmlFileName[PETSC_MAX_PATH_LEN]; /* input file name */
  char        tpsblk[PETSC_MAX_PATH_LEN]; /* input file name */
  PetscBool       flg;
  PetscCall(PetscOptionsGetString(NULL, NULL, "-f", graphmlFileName, sizeof(graphmlFileName), &flg));
  PetscCheck(flg, PETSC_COMM_WORLD, PETSC_ERR_USER, "Must indicate graphml file with the -f option");
  FILE* graphmlFile = fopen(graphmlFileName, "r");
  if (graphmlFile == NULL) {
    fprintf(stderr, "Error opening the file.\n");
    return 1;  // Return an error code or use another error handling method.
  }
  MPI_Comm_size(PETSC_COMM_WORLD,&size);
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  /*
    Enter the options for the number of holes and
    the total number of alpha electrons.
  */
  PetscInt  nholes;
  PetscInt  nalpha;
  PetscInt  DoS2 = 0;
  PetscInt  DoTPS = 0;
  PetscInt  DBGPrinting = 0;
  PetscInt  DoProj = 0;
  PetscInt  DoProjPrint = 0;
  PetscInt  doRepulsion = 0;
  PetscReal t_inp = -1.0;
  PetscReal hrepVal_inp = -1.0;
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
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-pd",&DBGPrinting,NULL));
  PetscCheck(flg, PETSC_COMM_WORLD, PETSC_ERR_USER, "Must indicate debug printing with the -pd option");
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-hrep",&doRepulsion,NULL));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-hproj",&DoProj,NULL));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-hprojPrint",&DoProjPrint,NULL));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-htps",&DoTPS,&flg));
  //PetscCheck(flg, PETSC_COMM_WORLD, PETSC_ERR_USER, "Must indicate whether or not to TPS with the -htps option (1=true, 0=false)");
  if(DoTPS == 1) {
    PetscCall(PetscOptionsGetString(NULL, NULL, "-htpsblk", tpsblk, sizeof(tpsblk), &flg));
    PetscCheck(flg, PETSC_COMM_WORLD, PETSC_ERR_USER, "Indicate tps blocks with -htpsblk option which inputs pairs of numbers separated by , (e.g. 1,2,3,4)");
  }
  if(doRepulsion == 1) {
    PetscCall(PetscOptionsGetReal(NULL,NULL,"-hrepval",&hrepVal_inp,&flg));
    PetscCheck(flg, PETSC_COMM_WORLD, PETSC_ERR_USER, "Indicate hrepval blocks with -hrepval option ");
  }

  /* 
   * Read the TPS Blocks into an array
   *
   */
  size_t* TPSBlock = malloc(MAX_TPS_BLOCKS * sizeof(size_t));
  char *p = tpsblk;
  int nblk = 0;
  setTPSBlockList(p, &nblk, TPSBlock) ;
  nblk = nblk/2;

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
  size_t* csfJKList  = malloc(sizeCSF * max_nbrs * sizeof(size_t));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Compute the operator matrix that defines the eigensystem, Ax=kx
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"====================================================="));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\nDEHamv2: Double Exchange Eigenproblem, n=%" PetscInt_FMT "\n",sizeTotal));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"=====================================================\n\n"));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," [Info] Nsites     \t\t %" PetscInt_FMT "\n",(size_t)nsites));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," [Info] Nalpha     \t\t %" PetscInt_FMT "\n",(size_t)nalpha));
  if(DBGPrinting) {
    PetscCall(PetscPrintf(PETSC_COMM_WORLD," [Info] Nbeta      \t\t %" PetscInt_FMT "\n",(size_t)nbeta));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD," [Info] NTPSBlocks \t\t %" PetscInt_FMT "\n",(size_t)nblk/2));
  }
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," [Info] Max Nbrs   \t\t %" PetscInt_FMT "\n",(size_t)max_nbrs));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," [Info] N(configurations)   \t %" PetscInt_FMT "\n",sizeCFG));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," [Info] N(CSFs)   \t\t %" PetscInt_FMT "\n",sizeCSF));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," [Info] t          \t\t %10.5f |t| \n",(double)t_inp));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," [Info] J          \t\t %10.5f |t| \n",(double)Jme));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," [Info] K          \t\t %10.5f |t| \n",(double)Kme));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," [Info] Nholes     \t\t %" PetscInt_FMT "\n",(size_t)nholes));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," [Info] NCPUs      \t\t %" PetscInt_FMT "\n",(size_t)size));
  if(doRepulsion) {
    PetscCall(PetscPrintf(PETSC_COMM_WORLD," [Info] Repulsion  \t\t %10.5f |t| \n",(double)hrepVal_inp));
  }
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," [Wait] Generating CSFs... \t\t \n"));

  generateConfigurations(nsites, (nsites-nholes), configList, &sizeCFG);

  // Sort the lists for binary search
  qsort(configList, sizeCFG, sizeof(size_t), compare);

  generateConfigurations(nelec, nalpha, csfList, &sizeCSF);
  // Sort the lists for binary search
  qsort(csfList, sizeCSF, sizeof(size_t), compare);
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," [Done] Generating CSFs !  \t\t \n"));

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
   * 
   *        i
   *        | 
   *       \|/
   *        . 
   *  6  5  4  3  2  1
   *  - --------------
   *  1  0  1  1  0  1
   *  1  1  1  1  1  1
   *  ----------------
   *  7  8  9 10 11 12  
   *        . 
   *       /|\
   *        | 
   *        j  
   *
   *        i0
   *        | 
   *       \|/
   *        . 
   *  4     3  2     1
   *  - --------------
   *  1  0  1  1  0  1
   *  1  1  1  1  1  1
   *  ----------------
   *  5  6  7  8  9 10  
   *        . 
   *       /|\
   *        | 
   *        j0 
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
  //int rows = sizeTotal;
  //int cols = rows;
  //double** matrix = declare_matrix(rows, cols);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Compute the operator matrix that defines the eigensystem, Ax=kx
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscInt       n=sizeTotal,i,Istart,Iend,nev,maxit,its,nconv;
  PetscInt		   ncv, mpd;
  // To print the matrix
  PetscViewer viewer; // declare a viewer object
  PetscViewerASCIIGetStdout(PETSC_COMM_WORLD, &viewer); // get the standard output


  // Symmetric Matrix
  PetscCall(MatCreate(PETSC_COMM_WORLD,&A));
  PetscCall(MatCreateSBAIJ(PETSC_COMM_WORLD,1,PETSC_DECIDE,PETSC_DECIDE,n,n,max_nbrs*nsites,NULL,max_nbrs*nsites,NULL,&A));
  PetscCall(MatSetType ( A, MATSBAIJ));
  PetscCall(MatMPIBAIJSetPreallocation(A,1,max_nbrs*nsites,NULL,max_nbrs*nsites,NULL));

  /*
   * Matrix for the S2 operator
    */
  if(DoS2){
    // Symmetric Matrix
    PetscCall(MatCreate(PETSC_COMM_WORLD,&S2));
    PetscCall(MatCreateSBAIJ(PETSC_COMM_WORLD,1,PETSC_DECIDE,PETSC_DECIDE,n,n,nholes*nalpha*nalpha,NULL,nholes*nalpha*nalpha,NULL,&S2));
    PetscCall(MatSetType (S2, MATSBAIJ));
    PetscCall(MatMPIBAIJSetPreallocation(S2,1,nholes*nalpha*nalpha,NULL,nholes*nalpha*nalpha,NULL));
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

    generateMonoCFGs(configList, sizeCFG, csfList, sizeCSF, &graph, posi, icfg[0], icsf[0], &monoCFGList, &monoMEs, t, Jme, Kme, doRepulsion, hrepVal_inp);
    //printf(" posi=%ld \n",i);
    for (int j = 0; j < igraph_vector_int_size(&monoCFGList); ++j) {
      PetscInt Jid = VECTOR(monoCFGList)[j];
      PetscReal val = VECTOR(monoMEs)[j];
      //PetscCall(MatSetValue(A,Jid,i,(PetscReal)val,ADD_VALUES));
      //if(i==Jid) printf(" Jid=%ld val=%10.5f \n",Jid,val);
      if( i > Jid ) PetscCall(MatSetValue(A,Jid,i,(PetscReal)val,INSERT_VALUES));
      else          PetscCall(MatSetValue(A,i,Jid,(PetscReal)val,INSERT_VALUES));
      //matrix[posi][Jid] += val;
      //printf("(%ld, %ld) = %f\n",posi, Jid, val);
    }

    igraph_vector_int_destroy(&monoCFGList);
    igraph_vector_destroy(&monoMEs);
  }
  PetscCall(PetscTime(&tt2));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," [Time] used to build the matrix: %f\n",tt2-tt1));


  PetscCall(PetscTime(&tt1));
  PetscCall(MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY));
  PetscCall(PetscTime(&tt2));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," [Time] used to assemble the matrix: %f\n",tt2-tt1));
  //if(DBGPrinting) {
  //  PetscCall(MatView(A, viewer));
  //}

  if(DoS2) {
    /*
     * Initialize S2 operator
      */
    PetscCall(PetscTime(&tt1));
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
    PetscCall(PetscPrintf(PETSC_COMM_WORLD," [Time] used to build the S2 operator: %f\n",tt2-tt1));


    PetscCall(PetscTime(&tt1));
    PetscCall(MatAssemblyBegin(S2,MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(S2,MAT_FINAL_ASSEMBLY));
    PetscCall(PetscTime(&tt2));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD," [Time] used to assemble the matrix: %f\n",tt2-tt1));
  }

  PetscCall(MatCreateVecs(A,NULL,&xr));
  PetscCall(MatCreateVecs(A,NULL,&xi));
  PetscCall(MatCreateVecs(A,NULL,&vs2));

  double xdi[nsites];
  if(DoTPS) {
    /*
     * Initialize distance matrix
     */
    for(size_t j=0;j<nsites;++j) {
      if((nsites & 1)) {
        // Odd
        //xdi[j] = -((nsites-1.0)/2.0) + j;
        /*
         * New cyclic ordering
         */
        if(j < (nsites+1)/2) {
          xdi[(nsites+1)/2 - j] = (-(nsites-1.0) + 2.0*j)/2.0;
        }
        else {
          xdi[nsites-(j-(nsites/2)+0)] = (-(nsites-1.0) + 2.0*j)/2.0;
        }
      }
      else {
        // Even
        //xdi[j] = 0.5 - ((nsites-0.0)/2.0) + j;
        /*
         * New cyclic ordering
         */
        if(j < nsites/2) {
          xdi[(nsites/2-j) - 1] = (-(nsites-1.0) + 2.0*j)/2.0;
        }
        else {
          xdi[(nsites-(j-(nsites/2)+0)) - 1] = (-(nsites-1.0) + 2.0*j)/2.0;
        }
      }
    }
    //for(size_t j=0;j<nsites;++j) {
    //  printf("%12f \n",xdi[j]);
    //}
  }

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
  tol = 1.e-24;
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
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," [Time] used to Solve EVP: %f\n",t2-t1));

  PetscCall(PetscPrintf(PETSC_COMM_WORLD," [Info] EVP: \n"));
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
    if(DoTPS){
        if(DoProj){
          PetscCall(PetscPrintf(PETSC_COMM_WORLD,
               "           E          ||Ax-Ex||/||Ex||         S                  Norm^2 (Ion, Pair)                   TPS (Diag, ExDiag)  \n"
               "   ----------------- ----------------------------------------------------------------------------------------------------------\n"));
        }
        else{
          PetscCall(PetscPrintf(PETSC_COMM_WORLD,
               "           E          ||Ax-Ex||/||Ex||         S                  TPS (Diag, ExDiag)   \n"
               "   ----------------- -----------------------------------------------------\n"));
        }
    }
    else if(DoProj){
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,
         "           E          ||Ax-Ex||/||Ex||         S              Norm^2         \n"
         "   ----------------- -------------------------------------------------\n"));
    }
    else{
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,
         "           E          ||Ax-Ex||/||Ex||         S   \n"
         "   ----------------- ----------------------------------\n"));
    }

    for (i=0;i<nev;i++) {
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

        spin = solveQuad(1.0, 1.0, -1.0*fabs(dot));
      }
      else spin = 100;

      double tpstot[nblk]; 
      double tpstotdiag[nblk]; 
      double tpstotexdiag[nblk]; 
      double weightholePairtot=0.0; 
      if(DoTPS) {
        /*
         * Get TPS value
         */
        double tps[nblk]; 
        double weightholePair=0.0; 
        double tpsdiag[nblk]; 
        double tpsexdiag[nblk]; 
        for(size_t k=0;k<nblk;++k) {
          tps[k] = 0.0;
          tpsdiag[k] = 0.0;
          tpsexdiag[k] = 0.0;
        }
        PetscCall(VecGetOwnershipRange(xr,&Istart,&Iend));
        double tpsval[nblk]; 
        for (size_t j=Istart;j<Iend;j++) {
          double val;
          int isDiag;
          int holePair;
          PetscInt ix[1];
          PetscScalar y[1];
          ix[0] = (size_t)j;
          PetscCall(VecGetValues(xr, 1, ix, y));
          size_t posi = j;
          cfgid = findCFGID(posi, sizeCFG, sizeCSF);
          csfid = findCSFID(posi, sizeCFG, sizeCSF);
          icfg[0] = configList[cfgid];
          icsf[0] = csfList[csfid];
          isDiag = 0;
          getTPSOperator(icfg[0], tpsval, xdi, configList, sizeCFG, nblk, TPSBlock, &graph, nsites, nholes, &isDiag, &holePair);
          if (holePair == 1) {
            weightholePair += y[0]*y[0];
          }
          for(size_t k=0;k<nblk;++k) {
            tps[k] += y[0]*y[0]*tpsval[k];
            if(isDiag) {
              tpsdiag[k] += y[0]*y[0]*tpsval[k];
            }
            else {
              tpsexdiag[k] += y[0]*y[0]*tpsval[k];
            }
          }
          //printf(" --> %d \n",isDiag);
          if(DBGPrinting) {
            PetscCall(PetscPrintf(PETSC_COMM_WORLD,"%12f \t",y[0]));
            printf("\t CFG: ");
						for(size_t l=0; l < nsites; ++l) {
              printf(" %d ",(1<<(l)) & icfg[0]? 1 : 0);
            }
            printf("\t CSF: ");
						for(size_t l=0; l < nelec; ++l) {
              printf(" %d ",(1<<(l)) & icsf[0]? 1 : 0);
            }
            PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\n"));
          }
        }
        MPI_Reduce(&tps, &tpstot, nblk, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
        MPI_Reduce(&tpsdiag, &tpstotdiag, nblk, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
        MPI_Reduce(&tpsexdiag, &tpstotexdiag, nblk, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
        MPI_Reduce(&weightholePair, &weightholePairtot, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
      }

  int               nstates = get_nstates((int)nholes,(int)nalpha);
  double            projvecfin[nev*nstates];
  double            normhole, normholefin, normspace;
  double            projvecfin2[nev*nstates];
      if(DoProj){

  double            projvec[nev*nstates], weightproj=0.0, weightproj2=0.0;
  double            projvec2[nev*nstates];
  double projmatrix[6][6];
        int sze = 0;
        int MS2 = 0;
        double XS=0.0;
        double V=0.0;
        int colm=6;
  const int            natomax=900;
        get_proj_general(xr, 
                         &Istart, 
                         &Iend, 
                         &nsites,
                         i, 
                         projvec, 
                         projvec2, 
                         &normhole, 
                         natomax,
                         (int)nholes,
                         nalpha,
                         (int)nholes,
                         XS,
                         V,
                         colm,
                         configList, 
                         sizeCFG,
                         csfList, 
                         sizeCSF,
                         &graph);
          MPI_Reduce(&projvec, &projvecfin, nev*nstates, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
          MPI_Reduce(&projvec2, &projvecfin2, nev*nstates, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
          MPI_Reduce(&normhole, &normholefin, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
          normspace = 0.0;
          for(int ii=0;ii<nstates;++ii){
            normspace += projvecfin[(i)*nstates + ii]*projvecfin[(i)*nstates + ii]/(1.0-normholefin);
          }
      }

      //if(DBGPrinting) {
      //  PetscCall(VecView(xr, viewer));
      //}

      if (im!=0.0) PetscCall(PetscPrintf(PETSC_COMM_WORLD," %9f%+9fi %12g\n",(double)re,(double)im,(double)error));
      else PetscCall(PetscPrintf(PETSC_COMM_WORLD,"   %12f       %12g       %12f",(double)re,(double)error,(double)fabs(spin)));
      if(DoTPS) {
        if(DoProj){
          PetscCall(PetscPrintf(PETSC_COMM_WORLD,"       %8f (%8f) (%8f) ",(double)normspace,(double)normholefin,(double)weightholePairtot));
          for(size_t j=0;j<nblk;++j) {
            PetscCall(PetscPrintf(PETSC_COMM_WORLD,"       %8f ( %8f %8f ) ",(double)tpstot[j], (double)tpstotdiag[j], (double)tpstotexdiag[j]));
          }
        }
        else{
          for(size_t j=0;j<nblk;++j) {
            PetscCall(PetscPrintf(PETSC_COMM_WORLD,"       %8f ( %8f %8f ) ",(double)tpstot[j], (double)tpstotdiag[j], (double)tpstotexdiag[j]));
          }
        }
      }
      else if(DoProj) {
          PetscCall(PetscPrintf(PETSC_COMM_WORLD,"       %8f (%8f) (%8f) ",(double)normspace,(double)normholefin,(double)weightholePairtot));
      }
      PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\n"));
      if(DoProjPrint) {
        PetscPrintf(PETSC_COMM_WORLD,"----------------------------\n");
        PetscPrintf(PETSC_COMM_WORLD,"Printing projection vectors \n");
        PetscPrintf(PETSC_COMM_WORLD,"----------------------------\n");
        for(int ii=0;ii<nstates;++ii){
          PetscPrintf(PETSC_COMM_WORLD,"%d %18f \n",ii,projvecfin[(i)*nstates + ii]);
        }
        PetscPrintf(PETSC_COMM_WORLD,"------     Done      -------\n");
      }

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
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"DEHamv2: Success !!!                      \n"));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"===========================================\n"));
  PetscCall(SlepcFinalize());
  return 0;
}
