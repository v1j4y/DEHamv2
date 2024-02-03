#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <x86intrin.h>

#include <igraph.h>

#include "slater_condon.h"

size_t get_matelem(size_t deti, size_t detj) ;

void printBits(size_t num, size_t len) ;

void printBits(size_t num, size_t len) ;

void generateConfigurations(size_t norb, size_t nelec, size_t* configAll, size_t* size) ;

// Function to find the positions of a list of configurations in a sorted list
void findPositions(size_t* configList, size_t sizeList, size_t* configs, size_t sizeConfigs, size_t* positions) ;

// Function to compare two configurations for qsort and bsearch
int compare(const void* a, const void* b) ;

// Function to calculate binomial coefficient using lgamma function
long long binomialCoeff(size_t n, size_t k) ;

void printPositions(size_t* positions, size_t size) ;

int getPhase(size_t alphaConfig, size_t newAlphaConfig, size_t h, size_t p) ;

// A function to declare a matrix of given size and initialize it to 0
double** declare_matrix(int rows, int cols) ;

// A function to save a matrix in a file in CSV format
void save_matrix(double** matrix, int rows, int cols, char* filename) ;
