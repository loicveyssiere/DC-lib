#ifndef CLEAN_TREE_CREATION_H
#define CLEAN_TREE_CREATION_H

#ifdef STATS
    #include <fstream>
#endif
#include "DC.h"

typedef struct createArgs_s {
    int *elemToNode, *sepToNode, *nodePart;
    int dimElem;
    int firstPart, lastPart;
    int globalNbElem, globalNbNodes;
} createArgs_t;

// Create the D&C tree and the element permutation, and compute the intervals of nodes
// and elements at each node of the tree
void clean_tree_creation (tree_t &tree, createArgs_t &create);

void clean_init_dc_tree (tree_t &tree, int firstElem, int lastElem,
    int firstNode, int lastNode, bool isSep, bool isLeaf, int depth);

void clean_init_create_args (createArgs_t &createArgs, int *elemToNode,
    int *sepToNode, int *nodePart, int firstPart, int lastPart,
    int globalNbElem, int dimElem, int globalNbNodes);

void DC_clean_create_tree (int *elemToNode, int nbElem, int dimElem, int nbNodes);

#endif // CLEAN_TREE_CREATION_H
