#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <x86intrin.h>

#include <igraph.h>

#include "slater_condon.h"

#define MAX_TPS_BLOCKS 100

double solveQuad(double a, double b, double c) ;

size_t get_matelem(size_t deti, size_t detj) ;

void printBits(size_t num, size_t len) ;

void printBitsDE(size_t num, size_t len, size_t nelecF1) ;

void generateConfigurations(size_t norb, size_t nelec, size_t* configAll, size_t* size) ;

// Function to find the positions of a list of configurations in a sorted list
void findPositions(size_t* configList, size_t sizeList, size_t* configs, size_t sizeConfigs, size_t* positions) ;

// Function to return the global ID given CFG and CSF ids
size_t findGlobalID(size_t cfgid, size_t csfid, size_t nCFG, size_t nCSF) ;

// Functions to find cfg IDs from globalID
size_t findCFGID(size_t globalID, size_t nCFG, size_t nCSF) ;

// Functions to find csf IDs from globalID
size_t findCSFID(size_t globalID, size_t nCFG, size_t nCSF) ;

// Function to compare two configurations for qsort and bsearch
int compare(const void* a, const void* b) ;

// Function to calculate binomial coefficient using lgamma function
long long binomialCoeff(size_t n, size_t k) ;

void printPositions(size_t* positions, size_t size) ;

int getPhase(size_t alphaConfig, size_t newAlphaConfig, size_t h, size_t p) ;

int getExecDegree(size_t detI, size_t detJ) ;

int getHoles_1ex(size_t detI, size_t detJ, size_t *holesOut) ;
int getPart_1ex(size_t detI, size_t detJ, size_t *particlesOut) ;

// Get maximum neighbors
int getMaxNeighbors(const igraph_t* graph, size_t nsites) ;

// A function to declare a matrix of given size and initialize it to 0
double** declare_matrix(int rows, int cols) ;

// A function to save a matrix in a file in CSV format
void save_matrix(double** matrix, int rows, int cols, char* filename) ;

// Function to generate all possible alpha determinants
void generateMonoCFGs(size_t* configList, size_t sizeConfig, size_t* csfList, size_t sizeCSF, const igraph_t* graph, size_t ipos, size_t Icfg, size_t Icsf, igraph_vector_int_t* monoCFGList, igraph_vector_t* monoMEs, double t, double Jme, double Kme) ;

void getdet(size_t Icsf, int *ideter, size_t* configAlpha, size_t sizeAlpha, int norb) ;

void adr (int *ideter, size_t *iii, size_t* configAlpha, size_t sizeAlpha, int norb) ;

void getS2Operator(size_t Icsf, igraph_vector_t* MElist, igraph_vector_int_t* Jdetlist, size_t *configAlpha, size_t sizeAlpha, const igraph_t* graph, int natom, int natomax) ;

// Function to get the TPS operator
void getTPSOperator(size_t detI, double *tpsval, size_t* cfgList, size_t sizeCFG, int nblk, size_t* TPSBlock, const igraph_t* graph, size_t nsites, size_t nholes) ;
