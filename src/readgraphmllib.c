#include "readgraphmllib.h"

int readGraphMLFile(FILE* file, igraph_t* graph) {
    if (igraph_read_graph_graphml(graph, file, 0) != IGRAPH_SUCCESS) {
        fprintf(stderr, "Error reading GraphML file\n");
        return 0;  // Return an error code or use another error handling method.
    }

    return 1;  // Successful graph read.
}

void getConnectedVertices(const igraph_t* graph, igraph_integer_t vertex_id, igraph_vector_int_t* result) {
    igraph_neighbors(graph, result, vertex_id, IGRAPH_ALL);
}
