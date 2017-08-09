#ifdef TREE_CREATION

#include <cstring>
#include <cmath>
#ifdef CILK
    #include <cilk/cilk.h>
#endif
#include <pthread.h>
#include <metis.h>
#include <algorithm>
#include <assert.h>

#include "tools.h"
//#include "permutations.h"
#include "tree_creation.h"
//#include "partitioning.h"

#include "alt_tree_creation.h"
#include "alt_partitioning.h"

extern tree_t *treeHead;
extern int *elemPerm, *nodePerm, *innerNodePerm;
extern int globalNbElem, globalNbNodes;

extern pthread_mutex_t metisMutex;

void compute_partition(int *componentToNode, int nbComponent, int dimComponent, int nbNodes, int nbPart, int *nodePart) {

    // Configure METIS & compute the node partitioning of the mesh
    int constraint = 1, objVal;
    int *graphIndex = new int [nbNodes + 1];
    int *graphValue = new int [nbComponent * 2];

    graph_to_csr_sym(componentToNode, graphIndex, graphValue, nbNodes, nbComponent);

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

void compute_innerNodePartition(int *componentToNode, int *nodePart,
        int *innerNodePart, int nbComponent, int dimComponent, int nbNodes) {

    // Initialize the global node permutation
    #ifdef OMP
        #pragma omp parallel for
        for (int i = 0; i < nbNodes; i++) {
    #elif CILK
        cilk_for (int i = 0; i < nbNodes; i++) {
    #endif
        innerNodePart[i] = nodePart[i];
    }

    // update node color
    int test = 0;
    #ifdef OMP
        #pragma omp parallel for
        for (int i = 0; i < nbComponent; i++) {
    #elif CILK
        cilk_for (int i = 0; i < nbComponent; i++) {
    #endif
        int base_color = nodePart[componentToNode[i * dimComponent]];
        bool same_color = true;
        for (int j = 1; j < dimComponent; j++) {
            if (nodePart[componentToNode[i * dimComponent + j]] != base_color) {
                same_color = false;
            }
        }
        if (!same_color) {
            // Version with 2 lines
            for (int j = 0; j < dimComponent; j++) {
                innerNodePart[componentToNode[i * dimComponent + j]] = -1;
            }

            // Version with 1 line
            // int min = nbNodes;
            // int min_color = nbNodes+1;
            // for (int j = 0; j < dimComponent; j++) {
            //     int node = componentToNode[i * dimComponent + j];
            //     if (nodePart[node] < min_color) {
            //         min_color = nodePart[node];
            //         min = node;
            //     }
            // }
            // innerNodePart[min] = -1;
        }
    }

    for (int i = 0; i < nbNodes; i++) {
        if (innerNodePart[i] == -1)
            test++;
    }

    printf("Number of nodes in sep: %d\n", test);
}


int alt_create_sepToNode2 (int *sepToNode, int *elemToNode, int firstSepElem,
                      int lastSepElem, int dimElem, int firstInnerNode, int lastInnerNode, int *nbSepNodes) {


    int nbInnerNodes = lastInnerNode - firstInnerNode + 1;

    int nbNodesIn = 0;
    int nbNodesOut = 0;
    int *tmp = new int [(lastSepElem - firstSepElem + 1) * dimElem];

    for (int i = firstSepElem*dimElem, j = 0; i < (lastSepElem+1)*dimElem; i++, j++) {

        int newNode = 0, oldNode = elemToNode[i];
        bool isMine;
        if (innerNodePerm[oldNode] >= firstInnerNode && innerNodePerm[oldNode] <= lastInnerNode) {
            isMine = true;
        } else {
            isMine = false;
        }

        bool isNewIn = true, isNewOut = true;
        // Search in IN
        for (newNode = 0; newNode < nbNodesIn; newNode++) {
            if (oldNode == tmp[newNode]) {
                isNewIn = false;
                break;
            }
        }
        // Search in OUT
        if (!isMine) {
            for (newNode = nbInnerNodes; newNode < nbInnerNodes+nbNodesOut; newNode++) {
                if (oldNode == tmp[newNode]) {
                    isNewOut = false;
                    break;
                }
            }
        }
        printf("(%d) %d // (%d) (%d)\n", isMine, innerNodePerm[oldNode], isNewIn, isNewOut);
        if (isMine && isNewIn) {
            tmp[nbNodesIn] = oldNode;
            nbNodesIn++;
        }
        else if (!isMine && isNewOut) {
            tmp[nbInnerNodes+nbNodesOut] = oldNode;
            nbNodesOut++;
        }
        sepToNode[j] = newNode;
    }

    delete[] tmp;
    assert(nbNodesIn == nbInnerNodes);
    *nbSepNodes = nbNodesIn + nbNodesOut;
}

int alt_create_sepToNode3 (int *sepToNode, int *elemToNode, int firstSepElem,
                      int lastSepElem, int dimElem, int firstInnerNode, int lastInnerNode, int *nbSepNodes) {


    int nbInnerNodes = lastInnerNode - firstInnerNode + 1;

    int nbNodesOut = 0;
    int *tmp = new int [(lastSepElem - firstSepElem + 1) * dimElem];

    for (int i = firstSepElem*dimElem, j = 0; i < (lastSepElem+1)*dimElem; i++, j++) {

        int newNode = 0, oldNode = elemToNode[i];
        bool isMine;
        if (innerNodePerm[oldNode] >= firstInnerNode && innerNodePerm[oldNode] <= lastInnerNode) {
            isMine = true;
            newNode = innerNodePerm[oldNode] - firstInnerNode;
        } else {
            isMine = false;
        }

        bool isNewOut = true;

        // Search in OUT
        if (!isMine) {
            for (newNode = 0; newNode < nbNodesOut; newNode++) {
                if (oldNode == tmp[newNode]) {
                    isNewOut = false;
                    break;
                }
            }
        }

        if (!isMine && isNewOut) {
            tmp[nbNodesOut] = oldNode;
            nbNodesOut++;
        }

        if (isMine) {
            sepToNode[j] = newNode;
        } else {
            sepToNode[j] = newNode + nbInnerNodes;
        }

    }

    delete[] tmp;
    *nbSepNodes = nbInnerNodes + nbNodesOut;
}

void alt_sep_partitioning (tree_t &tree, int *elemToNode, int *nodeIndex, int globalNbElem, int dimElem,
                       int firstSepElem, int lastSepElem, int firstNode, int lastNode,
                       int firstInnerNode, int lastInnerNode, int curNode, int depth) {


    // If there is not enough element in the separator
    int nbSepElem = lastSepElem - firstSepElem + 1;
    int nbSepPart = ceil (nbSepElem / (double)MAX_ELEM_PER_PART);
    int nbNodes = lastNode - firstNode + 1;
    int nbInnerNodes = lastInnerNode - firstInnerNode + 1;
    if (nbSepPart < 2 || nbSepElem <= MAX_ELEM_PER_PART) {
    // if (nbSepPart < 2 || nbInnerNodes <= MAX_NODE_PER_PART) {

       #ifdef MULTITHREADED_COMM
           // Set the node accessed by current D&C leaf
           fill_node_owner (elemToNode, firstSepElem, lastSepElem, dimElem, firstNode,
                            lastNode, curNode, true);
       #endif

       // Initialize the leaf
       alt_init_dc_tree (tree, firstSepElem, lastSepElem, 0, firstNode, lastNode, firstInnerNode, lastInnerNode, 0, true,
                     true, depth);

       // End of recursion
       return;
    }

    // Create temporal elemToNode containing the separator elements
    printf("Interval: %d %d and parts: %d\n", firstInnerNode, lastInnerNode, nbSepPart);
    int *sepToNode = new int [nbSepElem * dimElem];
    int *sepIndex = new int [nbInnerNodes];
    for (int j = 0; j < nbInnerNodes; j++) {
        sepIndex[j] = j;
    }
    int nbSepNodes = 0;

    alt_create_sepToNode3 (sepToNode, elemToNode, firstSepElem, lastSepElem, dimElem, firstInnerNode, lastInnerNode, &nbSepNodes);

    // Configure METIS & compute the node partitioning of the separators
    int constraint = 1, objVal;
    int *nodePart   = new int [nbSepNodes];

    compute_partition(sepToNode, nbSepElem, dimElem, nbSepNodes, nbSepPart, nodePart);

    // Create the separator D&C tree
    alt_tree_creation (tree, elemToNode, sepToNode, nodePart, nullptr, nodeIndex, sepIndex, globalNbElem,
                   dimElem, 0, nbSepPart-1, firstSepElem, lastSepElem, firstNode,
                   lastNode, firstInnerNode, lastInnerNode, 0, 0, curNode, true, depth+1);

    delete[] nodePart, delete[] sepToNode;
}

/* =============================================================================
     PRIVATE METHODS
============================================================================= */
void graph_to_csr_sym(int *pairToNode, int *graphIndex, int *graphValue, int nbNodes, int nbPairs) {

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

void coo_tocsr(int *pairToNode, int *graphIndex, int *graphValue, int nbNodes, int nbPairs) {

    //compute number of non-zero entries per row of A
    std::fill(graphIndex, graphIndex + nbNodes + 1, 0);

    for (int n = 0; n < nbPairs; n++){
        graphIndex[pairToNode[2*n]]++;
    }

    for(int i = 0, cumsum = 0; i < nbNodes; i++){
        int temp = graphIndex[i];
        graphIndex[i] = cumsum;
        cumsum += temp;
    }
    graphIndex[nbNodes] = nbPairs;

    for(int n = 0; n < nbPairs; n++){
        int row  = pairToNode[2*n];
        int dest = graphIndex[row];

        graphValue[dest] = pairToNode[2*n+1];
        graphIndex[row]++;
    }

    for(int i = 0, last = 0; i <= nbNodes; i++){
        int temp = graphIndex[i];
        graphIndex[i]  = last;
        last   = temp;
    }
}

void csr_tocoo(int *pairToNode, int *graphIndex, int *graphValue, int nbNodes, int nbPairs) {

    for(int i = 0; i < nbPairs; i++) {
        pairToNode[2*i] = 0;
        pairToNode[2*i+1] = graphValue[i];
    }

    int iGlobal = 0;
    int iEnd = graphIndex[0];
    for (int i = 0; i < nbNodes; i++) {
        iEnd = graphIndex[i+1];
        while(iGlobal < iEnd) {
            pairToNode[2*iGlobal] = i;
            iGlobal++;
        }
    }
}

#endif // TREE_CREATION
