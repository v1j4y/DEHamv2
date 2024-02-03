#include <igraph.h>
#include <stdio.h>
#include <unistd.h> /* unlink */

int readGraphMLFile(FILE* file, igraph_t* graph) ;

void getConnectedVertices(const igraph_t* graph, igraph_integer_t vertex_id, igraph_vector_int_t* result) ;
