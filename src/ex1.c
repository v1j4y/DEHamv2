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
  size_t nelecF1 = nsites - nholes;
  int natomax = 100;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Define Hamiltonian
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  double Jme =  0.01;
  double Kme =  6.0;
  double t = -1.0;

  size_t sizeCFG  = binomialCoeff(nsites, (nsites-nholes));
  size_t sizeCSF  = binomialCoeff(nelec, nalpha);
  size_t sizeTotal  = sizeCFG*sizeCSF;

  size_t* configList = malloc(sizeCFG * sizeof(size_t));
  size_t* csfList    = malloc(sizeCSF * sizeof(size_t));
  // List for storing JK elems
  int max_nbrs = getMaxNeighbors(&graph, nsites);
  printf(" Max Nbrs           # = %d \n",max_nbrs);
  size_t* csfJKList  = malloc(sizeCSF * max_nbrs * sizeof(size_t));

  generateConfigurations(nsites, (nsites-nholes), configList, &sizeCFG);

  // Sort the lists for binary search
  qsort(configList, sizeCFG, sizeof(size_t), compare);
  printf(" Configuration List # = %ld \n",sizeCFG);
  for( int i=0; i<sizeCFG; ++i ) {
    printBits(configList[i], nsites);
  }

  generateConfigurations(nelec, nalpha, csfList, &sizeCSF);

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


  // Sort the lists for binary search
  qsort(csfList, sizeCSF, sizeof(size_t), compare);
  printf(" CSF           List # = %ld \n",sizeCSF);
  //for( int i=0; i<sizeCSF; ++i ) {
  //  printf(" --- %d ---\n",i);
  //  printBits(csfList[i], nelec);
  //  for( int j=i; j<sizeCSF; ++j ) {
  //    printBits(csfList[j], nelec);
  //    int excDeg = getExecDegree(csfList[i], csfList[j]);
  //    if (excDeg == 1) {
  //      size_t holes[1];
  //      size_t part[1];
  //      excDeg = getHoles_1ex(csfList[i], csfList[j], holes);
  //      excDeg = getPart_1ex(csfList[i], csfList[j], part);
  //      size_t p=part[0];
  //      size_t h=holes[0];
  //      if( (p <= nsites && h <= nsites) || (p > nsites && h > nsites) ) {
  //        igraph_bool_t res;
  //        igraph_integer_t v1, v2;
  //        v1 = ((h-1) % nsites);
  //        v2 = ((p-1) % nsites);
  //        igraph_are_connected(&graph, v1, v2, &res);
  //        if (res) printf(" %d %d = %d (h=%ld p=%ld) => %d\n", i, j, excDeg, v1, v2, res);
  //        //printf(" %d %d = %d (h=%ld p=%ld) => %d\n", i, j, excDeg, v1, v2, res);
  //      }
  //    }
  //  }
  //}

  //for( int i=0; i<sizeCFG; ++i ) {
  //  printf(" --- %d ---\n",i);
  //  printBits(configList[i], nsites);
  //  size_t Icfg = configList[i];
  //  for( int j=0; j<nsites; ++j ) {
  //    if ((Icfg >> j) & 1) {
  //      // Get the connected vertices
  //      igraph_vector_int_t orbital_id_allowed;
  //      igraph_vector_int_init(&orbital_id_allowed, 0);
  //      getConnectedVertices(&graph, (igraph_integer_t)j, &orbital_id_allowed);
  //      printf(" > %d \n",j);

  //      // Calculate J
  //      // Loop over each connected vertex
  //      for (size_t k = 0; k < igraph_vector_int_size(&orbital_id_allowed); ++k) {
  //        size_t orbital_id = VECTOR(orbital_id_allowed)[k];
  //        if((( (Icfg >> j ) & 1) & ((Icfg >> orbital_id) & 1)) & (j>orbital_id)) {
  //          // Find the real index
  //          size_t i0, j0;
  //          size_t mask = (((size_t)1 << (j+1)) - 1);
  //          i0 = j - popcnt ( mask ^ (mask & Icfg));
  //          mask = (((size_t)1 << (orbital_id+1))-1);
  //          j0 = orbital_id - popcnt ( mask ^ (mask & Icfg));
  //          //printf(" > \t %d %ld\n",j, orbital_id);
  //          if ( (i0 != j0) ) {
  //            // Loop over CSFs
  //            for( int l=0; l<sizeCSF; ++l ) {
  //              size_t Icsf = csfList[l];
  //              if( (( (Icsf >> i0 ) & 1) ^ ((Icsf >> j0) & 1))  ) {
  //                //printf(" %d | (%ld %ld)\n",l,i0,j0);
  //                //printBits(Icsf, nelec);
  //                size_t Jcsf = Icsf;
  //                Jcsf = Jcsf ^ ((size_t)1 << (i0));
  //                Jcsf = Jcsf ^ ((size_t)1 << (j0));
  //                //printBits(Jcsf, nelec);
  //              }
  //            }
  //          }

  //        }
  //      }

  //      // Calculate K
  //      // Check if there's K
  //      for( int l=0; l<sizeCSF; ++l ) {
  //        size_t Icsf = csfList[l];
  //        size_t i0, j0;
  //        size_t mask = (((size_t)1 << (j+1)) - 1);
  //        i0 = j - popcnt ( mask ^ (mask & Icfg));
  //        j0 = j + nelecF1;
  //        if( (( (Icsf >> i0 ) & 1) ^ ((Icsf >> j0) & 1))  ) {
  //          printf(" %d | (%ld %ld)\n",l,i0,j0);
  //          printBits(Icsf, nelec);
  //          size_t Jcsf = Icsf;
  //          Jcsf = Jcsf ^ ((size_t)1 << (i0));
  //          Jcsf = Jcsf ^ ((size_t)1 << (j0));
  //          printBits(Jcsf, nelec);
  //        }
  //      }
  //    }
  //  }
  //}

  //for( int i=0; i<sizeCSF; ++i ) {
  //  printf("\n --- %d ---\n",i);
  //  printBits(csfList[i], nelec);
  //  size_t Icsf = csfList[i];
  //  for( int j=0; j<nsites; ++j ) {

  //    size_t orbital_id = j + nelecF1;
  //    printf(" > \t %d %ld\n",j, orbital_id);
  //    if((( (Icsf >> j ) & 1) ^ ((Icsf >> orbital_id) & 1))) {
  //      size_t Jcsf = Icsf;
  //      Jcsf = Jcsf ^ ((size_t)1 << (j));
  //      Jcsf = Jcsf ^ ((size_t)1 << (orbital_id));
  //      printBits(Jcsf, nelec);
  //    }
  //  }
  //}

  printf(" Ne=%ld Na=%ld Nb=%ld \n", nelec, nalpha, nbeta);

  // Get the global address of a CFGxCSF pair.
  size_t icfg[1];
  size_t icsf[1];
  size_t iglobalid;
  //icfg[0] = 11;
  //icsf[0] = 15 - 8 + 16;
  //icfg[0] = configList[4];
  //icsf[0] = csfList[22];
  icfg[0] = configList[1];
  icsf[0] = csfList[1];
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
  generateMonoCFGs(configList, sizeCFG, csfList, sizeCSF, &graph, icfg[0], icsf[0], &monoCFGList, &monoMEs, Jme, Kme);
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
