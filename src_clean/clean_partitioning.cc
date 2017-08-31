//#include <cstring>
//#include <cmath>
// #ifdef CILK
//     #include <cilk/cilk.h>
// #endif
// #include <pthread.h>
// #include <metis.h>
// #include <algorithm>
// #include <assert.h>

//#include "tools.h"
//#include "permutations.h"
//#include "tree_creation.h"
//#include "partitioning.h"
//
// #include "alt_tree_creation.h"
// #include "alt_partitioning.h"

#include <metis.h>
#include <algorithm>
#include <cmath>

#include "clean_partitioning.h"
#include "partitioning.h"

extern pthread_mutex_t metisMutex;

void clean_preprocess_partition (int *elemToNode, int nbElem, int dimElem, int nbNodes, int nbPart, int *nodePart) {

    // Configure METIS & compute the node partitioning of the mesh
    int constraint = 1, objVal;
    int *graphIndex = new int [nbNodes + 1];
    int *graphValue = new int [nbElem * 2];

    if (dimElem == 2) {
        // Case: elements are edges
        clean_graph_to_csr_sym(elemToNode, graphIndex, graphValue, nbNodes, nbElem);
    } else {
        // Case: elements are triangles, thetradrons or octoadrons
        mesh_to_nodal (graphIndex, graphValue, elemToNode, nbElem, dimElem, nbNodes);
    }

    printf("\nLet's partition\n");
    // Execution is correct without mutex although cilkscreen detects many race
    // conditions. Check if the problem is solved with future version of METIS (5.0)
    pthread_mutex_lock (&metisMutex);
    METIS_PartGraphRecursive (&nbNodes, &constraint, graphIndex, graphValue, nullptr,
                              nullptr, nullptr, &nbPart, nullptr, nullptr, nullptr,
                              &objVal, nodePart);
    pthread_mutex_unlock (&metisMutex);

    delete[] graphValue, delete[] graphIndex;
}

int clean_get_nb_parts(int nbOfItems) {
    return ceil(nbOfItems / (double) MAX_ELEM_PER_PART);
}

/* -----------------------------------------------------------------------------
    Graph manipulation
----------------------------------------------------------------------------- */
void clean_graph_to_csr_sym(int *pairToNode, int *graphIndex, int *graphValue, int nbNodes, int nbPairs) {

    std::fill(graphIndex, graphIndex + nbNodes + 1, 0);

    for (int n = 0; n < nbPairs; n++){
        graphIndex[pairToNode[2*n]]++;
        graphIndex[pairToNode[2*n+1]]++;
    }

    for(int i = 0, cumsum = 0; i < nbNodes; i++){
        int temp = graphIndex[i];
        graphIndex[i] = cumsum;
        cumsum += temp;
    }
    graphIndex[nbNodes] = 2*nbPairs;

    for(int n = 0; n < nbPairs; n++){
        int row  = pairToNode[2*n];
        int dest = graphIndex[row];
        graphValue[dest] = pairToNode[2*n+1];
        graphIndex[row]++;

        row = pairToNode[2*n+1];
        dest = graphIndex[row];
        graphValue[dest] = pairToNode[2*n];
        graphIndex[row]++;
    }

    for(int i = 0, last = 0; i <= nbNodes; i++){
        int temp = graphIndex[i];
        graphIndex[i]  = last;
        last   = temp;
    }
}
