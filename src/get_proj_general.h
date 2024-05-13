#include <stdio.h>
#include <petsctime.h>
#include <slepceps.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

#include "functions.h"
#include "igraph.h"

void get_proj_general(Vec valxr, 
                      PetscInt *Istart, 
                      PetscInt *Iend, 
                      size_t *natom, 
                      int iroot, 
                      double *projvec, 
                      double *projvec2, 
                      double *normproj, 
                      const int natomax,
                      int sze,
                      int MS2,
                      int nholes,
                      double XS,
                      double V,
                      int colm,
                      //double *projmatrix,
                      size_t *configList,
                      size_t sizeCFG,
                      size_t *csfList,
                      size_t sizeCSF,
                      const igraph_t* graph);
