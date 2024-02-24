#include "doubleexchange.h"
#include "readgraphmllib.h"

size_t get_matelem(size_t deti, size_t detj) {
  exc_number_t exij;
  determinant_t d1[1];
  determinant_t d2[1];
  d1[0] = deti;
  d2[0] = detj;
  exij = exc_degree(1, d1, d2);
  return exij;
}

// Function to calculate binomial coefficient using lgamma function
long long binomialCoeff(size_t n, size_t k) {
    return round(exp(lgamma(n+1) - lgamma(k+1) - lgamma(n-k+1)));
}

void printBits(size_t num, size_t len) {
    for (int bit = len - 1; bit >= 0; --bit) {
        printf("%ld", (num >> bit) & 1);
    }
    printf("\n");
}

void generateConfigurations(size_t norb, size_t nelec, size_t* configAll, size_t* size) {
    *size = 0;

    for (size_t i = 0; i < (1 << norb); ++i) {
        if (popcnt (i) == nelec) {
            configAll[(*size)++] = i;
        }
    }
}

// Functions to find cfg IDs from globalID
size_t findCFGID(size_t globalID, size_t nCFG, size_t nCSF) {
    size_t cfgid = (globalID / nCSF);
    return cfgid;
}

// Functions to find csf IDs from globalID
size_t findCSFID(size_t globalID, size_t nCFG, size_t nCSF) {
    size_t cfgid = (globalID / nCSF);
    size_t csfid = globalID - cfgid*nCSF;
    return csfid;
}

// Function to compare two configurations for qsort and bsearch
int compare(const void* a, const void* b) {
    return (*(size_t*)a - *(size_t*)b);
}

// Function to find the positions of a list of configurations in a sorted list
void findPositions(size_t* configList, size_t sizeList, size_t* configs, size_t sizeConfigs, size_t* positions) {
    for (size_t i = 0; i < sizeConfigs; ++i) {
        size_t* item = (size_t*) bsearch(&configs[i], configList, sizeList, sizeof(size_t), compare);
        if (item != NULL) {
            positions[i] = item - configList;
        } else {
            positions[i] = -1;
        }
    }
}

// Function to return the global ID given CFG and CSF ids
size_t findGlobalID(size_t cfgid, size_t csfid, size_t nCFG, size_t nCSF) {
  size_t globalid = cfgid * nCSF + csfid;
  return globalid;
}

void printPositions(size_t* positions, size_t size) {
    for (size_t i = 0; i < size; ++i) {
        printf("Position: %ld\n", positions[i]);
    }
}

int getPhase(size_t alphaConfig, size_t newAlphaConfig, size_t h, size_t p) {

    // Phase
    size_t nperm;

    determinant_t d1[1];
    determinant_t d2[1];
    d1[0] = alphaConfig;
    d2[0] = newAlphaConfig;
    orbital_t h1[1];
    orbital_t p2[1];
    h1[0] = h;
    p2[0] = p;
    nperm = get_nperm_single((size_t) 1, d1, d2, h1, p2);
    // size_t phase = ((size_t) 1) & nperm;
    //printf(" %llu %llu (%d, %d) nperm = %d phase=%d \n",d1[0], d2[0], i,orbital_id,nperm,phase);
    return nperm;
}

int getExecDegree(size_t detI, size_t detJ) {
    // Result
    size_t result;
    determinant_t d1[1];
    determinant_t d2[1];
    d1[0] = detI;
    d2[0] = detJ;
    result = exc_degree((size_t) 1, d1, d2);
    return result;
}

int getHoles_1ex(size_t detI, size_t detJ, size_t *holesOut) {
    // Result
    size_t nholes;

    determinant_t d1[1];
    determinant_t d2[1];
    d1[0] = detI;
    d2[0] = detJ;
    orbital_t holes[2];
    nholes = get_holes((size_t) 1, d1, d2, holes);
    holesOut[0] = holes[0];
    return nholes;
}

int getPart_1ex(size_t detI, size_t detJ, size_t *particlesOut) {
    // Result
    size_t nparticles;

    determinant_t d1[1];
    determinant_t d2[1];
    d1[0] = detI;
    d2[0] = detJ;
    orbital_t particles[2];
    nparticles = get_particles((size_t) 1, d1, d2, particles);
    particlesOut[0] = particles[0];
    return nparticles;
}

// Get maximum neighbors
int getMaxNeighbors(const igraph_t* graph, size_t nsites) {
    int max_nbrs = 0;
    for(int i=0; i<nsites; ++i) {
        int nbrs = getNumberOfConnectedVertices(graph, (igraph_integer_t)i);
        if(nbrs > max_nbrs) max_nbrs = nbrs;
    }
    return(max_nbrs);
}

// A function to declare a matrix of given size and initialize it to 0
double** declare_matrix(int rows, int cols) {
    // Allocate memory for the matrix
    double** matrix = (double**)malloc(rows * sizeof(double*));
    for (int i = 0; i < rows; i++) {
        matrix[i] = (double*)malloc(cols * sizeof(double));
    }

    // Initialize the matrix to 0
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            matrix[i][j] = 0;
        }
    }

    // Return the matrix
    return matrix;
}

// A function to fill up the non zero elements of a matrix
void fill_matrix(int** matrix, int rows, int cols) {
    // Loop through the matrix
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            // If the element is 0, assign a random value between 1 and 10
            if (matrix[i][j] == 0) {
                matrix[i][j] = rand() % 10 + 1;
            }
        }
    }
}

// A function to print a matrix
void print_matrix(int** matrix, int rows, int cols) {
    // Loop through the matrix
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            // Print the element with a space
            printf("%d ", matrix[i][j]);
        }
        // Print a new line
        printf("\n");
    }
}

// A function to save a matrix in a file in CSV format
void save_matrix(double** matrix, int rows, int cols, char* filename) {
    // Open the file in write mode
    FILE* file = fopen(filename, "w");

    // Check if the file is opened successfully
    if (file == NULL) {
        printf("Error: could not open the file %s\n", filename);
        return;
    }

    // Loop through the matrix
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            // Write the element to the file with a comma
            fprintf(file, "%10.5f,", matrix[i][j]);
        }
        // Write a new line to the file
        fprintf(file, "\n");
    }

    // Close the file
    fclose(file);
}

// Function to generate all possible alpha determinants
void generateMonoCFGs(size_t* configList, size_t sizeConfig, size_t* csfList, size_t sizeCSF, const igraph_t* graph, size_t Icfg, size_t Icsf, igraph_vector_int_t* monoCFGList, igraph_vector_t* monoMEs, double Jme, double Kme) {
    // Get the number of orbitals
    size_t norb = igraph_vcount(graph);
    size_t nholes = norb - popcnt ( Icfg );
    size_t nelec = 2*norb - nholes;
    size_t nelecF1 = norb - nholes;
    int phase = 1;
    size_t posCFG;
    findPositions(configList, sizeConfig, &Icfg, 1, &posCFG);

    // Loop over each orbital
    for (size_t i = 0; i < norb; ++i) {
        // Check if the orbital is occupied
        if ((Icfg >> i) & 1) {
            // Get the connected vertices
            igraph_vector_int_t orbital_id_allowed;
            igraph_vector_int_init(&orbital_id_allowed, 0);
            getConnectedVertices(graph, (igraph_integer_t)i, &orbital_id_allowed);

            // Loop over each connected vertex
            for (size_t j = 0; j < igraph_vector_int_size(&orbital_id_allowed); ++j) {
                size_t orbital_id = VECTOR(orbital_id_allowed)[j];

                // Check if the connected vertex is unoccupied
                if (!((Icfg >> orbital_id) & (size_t)1)) {

                    // Create a new alpha determinant by moving the electron
                    size_t Jcfg = Icfg ^ (((size_t)1 << i) | ((size_t)1 << orbital_id));

                    // Find the position of the new alpha determinant in the list and add it to alphaDeterminants
                    size_t posCFG;
                    findPositions(configList, sizeConfig, &Jcfg, 1, &posCFG);

                    // Now the CSFs and the phase

                    // Prepare the mask
                    size_t maskI = ~((size_t)1 << (orbital_id));
                    size_t Icfgmask = Icfg & maskI;
                    size_t maskJ = ~((size_t)1 << (i));
                    size_t Jcfgmask = Jcfg & maskJ;

                    // Find the equivalent ids for the CSFs
                    size_t i0, j0;
                    size_t mask = (((size_t)1 << (i+1)) - 1);
                    i0 = i - popcnt ( mask ^ (mask & Icfg));

                    mask = (((size_t)1 << (orbital_id+1))-1);
                    j0 = orbital_id - popcnt ( mask ^ (mask & Jcfg));
                    //printf(" (%ld, %ld) | i0=%ld , j0=%ld \n",i,orbital_id,i0,j0);
                    size_t Jcsf = Icsf;
                    Jcsf = Jcsf ^ ((size_t)1 << (j0));
                    Jcsf = Jcsf ^ ((size_t)1 << (i0));
                    //printBits(Jcsf, nelec);
                    size_t posCSF;
                    findPositions(csfList, sizeCSF, &Jcsf, 1, &posCSF);

                    // Add the position of the new alpha determinant to the list
                    size_t pos;
                    pos = findGlobalID(posCFG, posCSF, sizeConfig, sizeCSF);
                    igraph_vector_int_push_back(monoCFGList, pos);

                    phase = getPhase(Icsf, Jcsf, (i0)+1, (j0)+1);
                    phase = phase & 1 == 1 ? -1 : 1;

                    // Add the position of the new alpha determinant to the list
                    igraph_vector_push_back(monoMEs, phase);
                }
            }

            // Calculate J
            // Loop over each connected vertex
            for (size_t j = 0; j < igraph_vector_int_size(&orbital_id_allowed); ++j) {
                size_t orbital_id = VECTOR(orbital_id_allowed)[j];
                if((( (Icfg >> i ) & 1) & ((Icfg >> orbital_id) & 1)) & (i>orbital_id)) {
                    // Find the real index
                    size_t i0, j0;
                    size_t mask = (((size_t)1 << (i+1)) - 1);
                    i0 = i - popcnt ( mask ^ (mask & Icfg));
                    mask = (((size_t)1 << (orbital_id+1))-1);
                    j0 = orbital_id - popcnt ( mask ^ (mask & Icfg));
                    //printf(" > \t %d %ld\n",j, orbital_id);
                    if ( (i0 != j0) ) {
                        // Loop over CSFs
                        if( (( (Icsf >> i0 ) & 1) ^ ((Icsf >> j0) & 1))  ) {
                          //printf(" %d | (%ld %ld)\n",l,i0,j0);
                          //printBits(Icsf, nelec);
                          size_t Jcsf = Icsf;
                          Jcsf = Jcsf ^ ((size_t)1 << (i0));
                          Jcsf = Jcsf ^ ((size_t)1 << (j0));
                          //printBits(Jcsf, nelec);
                          size_t posCSF;
                          findPositions(csfList, sizeCSF, &Jcsf, 1, &posCSF);

                          // Add the position of the new alpha determinant to the list
                          size_t pos;
                          pos = findGlobalID(posCFG, posCSF, sizeConfig, sizeCSF);
                          igraph_vector_int_push_back(monoCFGList, pos);

                          // Add the position of the new alpha determinant to the list
                          igraph_vector_push_back(monoMEs, Jme);
                        }
                    }
                }
            }

            // Calculate K
            // Check if there's K
            size_t i0, j0;
            size_t mask = (((size_t)1 << (i+1)) - 1);
            i0 = i - popcnt ( mask ^ (mask & Icfg));
            j0 = i + nelecF1;
            if( (( (Icsf >> i0 ) & 1) ^ ((Icsf >> j0) & 1))  ) {
              //printf(" %d | (%ld %ld)\n",l,i0,j0);
              //printBits(Icsf, nelec);
              size_t Jcsf = Icsf;
              Jcsf = Jcsf ^ ((size_t)1 << (i0));
              Jcsf = Jcsf ^ ((size_t)1 << (j0));
              //printBits(Jcsf, nelec);

              size_t posCSF;
              findPositions(csfList, sizeCSF, &Jcsf, 1, &posCSF);

              // Add the position of the new alpha determinant to the list
              size_t pos;
              pos = findGlobalID(posCFG, posCSF, sizeConfig, sizeCSF);
              igraph_vector_int_push_back(monoCFGList, pos);

              // Add the position of the new alpha determinant to the list
              igraph_vector_push_back(monoMEs, Kme);
            }

            igraph_vector_int_destroy(&orbital_id_allowed);
        }
    }

}

// Main function that calculates MEs
//void getAllDoubleExchangeMEs(size_t Idet, igraph_vector_t* MElist, igraph_vector_t* Jdetlist, size_t *configCFG, size_t sizeCFG, size_t *configCSF, size_t sizeCSF, const igraph_t* graph) {
//    int phaseAlpha;
//    int phaseBeta;
//    //Find alpha and beta ids
//    size_t cfgID = findCFGID(Idet, sizeCFG, sizeCSF);
//    size_t csfID = findCSFID(Idet, sizeCFG, sizeCSF);
//
//    // Find allowed excitations
//    igraph_vector_t alphaDeterminants;
//    igraph_vector_init(&alphaDeterminants, 0);
//    igraph_vector_t alphaMEs;
//    igraph_vector_init(&alphaMEs, 0);
//    generateDeterminants(configAlpha, sizeAlpha, graph, configAlpha[alphaID], configBeta[betaID], &alphaDeterminants, &alphaMEs, 1);
//    igraph_vector_t betaDeterminants;
//    igraph_vector_init(&betaDeterminants, 0);
//    igraph_vector_t betaMEs;
//    igraph_vector_init(&betaMEs, 0);
//    generateDeterminants(configBeta, sizeBeta, graph, configBeta[betaID], configAlpha[alphaID], &betaDeterminants, &betaMEs, 0);
//
//    for (size_t j = 0; j < igraph_vector_size(&alphaDeterminants); ++j) {
//        size_t alphaJ = VECTOR(alphaDeterminants)[j];
//        phaseAlpha = VECTOR(alphaMEs)[j];
//
//        size_t foundGlobalID = findGlobalID(alphaJ, betaID, sizeAlpha);
//
//        igraph_vector_push_back(Jdetlist, foundGlobalID);
//        igraph_vector_push_back(MElist, phaseAlpha);
//    }
//    for (size_t k = 0; k < igraph_vector_size(&betaDeterminants); ++k) {
//
//        size_t betaK = VECTOR(betaDeterminants)[k];
//        phaseBeta = VECTOR(betaMEs)[k];
//
//        size_t foundGlobalID = findGlobalID(alphaID, betaK, sizeAlpha);
//
//        igraph_vector_push_back(Jdetlist, foundGlobalID);
//        igraph_vector_push_back(MElist, phaseBeta);
//    }
//    igraph_vector_destroy(&alphaDeterminants);
//    igraph_vector_destroy(&alphaMEs);
//    igraph_vector_destroy(&betaDeterminants);
//    igraph_vector_destroy(&betaMEs);
//}
