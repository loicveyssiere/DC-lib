#ifdef CILK
    #include <cilk/cilk.h>
#endif
// #include <pthread.h>
// #include <limits.h>
// #include <cstring>
// #include <iostream>
//
// #include "tools.h"
// #include "coloring.h"
// #include "permutations.h"
// #include "partitioning.h"
// #include "tree_creation.h"

#include "DC.h"
#include "clean_tree_creation.h"
#include "clean_partitioning.h"

// Global variables used in order to persist from one call to the library to another:
// the D&C tree and permutations
extern tree_t   *treeHead;
extern int      *elemPerm;
extern int      *nodePerm;

#ifdef TREE_CREATION

// Mutex to avoid race condition in merge permutations
extern pthread_mutex_t mergeMutex;

void clean_tree_creation(tree_t &tree, createArgs_t &args) {
    int nbPart = args.lastPart - args.firstPart + 1;
    int localNbElem = tree.lastElem - tree.firstElem + 1;
    int localNbNodes = tree.lastNode - tree.firstNode + 1;

    /* -- If current node is a leaf -- */
    if (nbPart < 2 || localNbElem <= MAX_ELEM_PER_PART) {

        #ifdef MULTITHREADED_COMM
            // Set the node accessed by current D&C leaf
            fill_node_owner (elemToNode, firstElem, lastElem, dimElem, firstNode,
                             lastNode, curNode, isSep);
        #endif

        tree.isLeaf = true;
        return;
    }

    /* -- If current node is a not a leaf -- */
}

void clean_init_dc_tree (tree_t &tree, int firstElem, int lastElem,
    int firstNode, int lastNode, bool isSep, bool isLeaf, int depth) {

    tree.firstElem      = firstElem;
    tree.lastElem       = lastElem;

    tree.firstNode      = firstNode;
    tree.lastNode       = lastNode;

    tree.depth          = depth;
    tree.isSep          = isSep;
    tree.isLeaf         = isLeaf;
}

void clean_init_create_args (createArgs_t &createArgs, int *elemToNode,
    int *sepToNode, int *nodePart, int firstPart, int lastPart,
    int globalNbElem, int dimElem, int globalNbNodes) {

    createArgs.elemToNode = elemToNode;
    createArgs.sepToNode = sepToNode;
    createArgs.nodePart = nodePart;

    createArgs.firstPart = firstPart;
    createArgs.lastPart = lastPart;

    createArgs.globalNbElem = globalNbElem;
    createArgs.globalNbNodes = globalNbNodes;
    createArgs.dimElem = dimElem;
}

void DC_clean_create_tree(int *elemToNode, int nbElem, int dimElem, int nbNodes) {
    // Allocate the D&C tree & the permutation functions
    elemPerm = new int [nbElem];
    nodePerm = new int [nbNodes];
    int *nodeIndex = new int[nbNodes];
    int *nodePart = new int [nbNodes];
    int nbPart = clean_get_nb_parts(nbElem);

    treeHead = new tree_t;
    createArgs_t *createArgs = new createArgs_t;

    /* -- Initialization of the tree sctucture -- */
    clean_init_dc_tree(*treeHead, 0, nbElem-1, 0, nbNodes-1, false, false, 0);

    /* -- Initialization of the main structure -- */
    clean_init_create_args(*createArgs, elemToNode, nullptr, nodePart, 0, nbPart-1, nbElem, dimElem, nbNodes);

    #ifdef MULTITHREADED_COMM
        // Initialize the node owner array
        init_node_owner(nbNodes);
    #endif

    clean_preprocess_partition(elemToNode, nbElem, dimElem, nbNodes, nbPart, nodePart);

    // Create node permutation from node partition
    //DC_create_permutation(nodePerm, nodePart, nbNodes, nbPart);

    #ifdef OMP
        #pragma omp parallel for
        for (int i = 0; i < nbNodes; i++) {
    #elif CILK
        cilk_for (int i = 0; i < nbNodes; i++) {
    #endif
        nodeIndex[i] = i;
    }

    #ifdef OMP
        #pragma omp parallel for
        for (int i = 0; i < nbElem; i++) {
    #elif CILK
        cilk_for (int i = 0; i < nbElem; i++) {
    #endif
        elemPerm[i] = i;
    }

    // Create D&C tree
    #ifdef OMP
        #pragma omp parallel
        #pragma omp single nowait
    #endif

    clean_tree_creation(*treeHead, *createArgs);

    delete[] nodePart, delete[] nodeIndex;

    // Vectorial version with coloring of the leaves of the D&C tree
    #ifdef DC_VEC
        coloring(ElemToNode, nbElem, dimElem, nbNodes);
    #endif
}

#endif // TREE_CREATION
