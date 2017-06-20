#ifndef ALT_PARTITIONING_H
#define ALT_PARTITIONING_H

#include "DC.h"

//void alt_partitioning(int *componentToNode, int nbComponent, int dimComponent, int nbNodes);

void compute_partition(int *componentToNode, int nbComponent, int dimComponent, int nbNodes, int nbPart, int *nodePart);

void compute_innerNodePartition(int *componentToNode, int *nodePart,
        int *innerNodePart, int nbComponent, int dimComponent, int nbNodes);

void alt_sep_partitioning (tree_t &tree, int *elemToNode, int *oldBijection, int globalNbElem, int dimElem,
                       int firstSepElem, int lastSepElem, int firstNode, int lastNode,
                       int firstInnerNode, int lastInnerNode, int curNode, long depth);

int alt_create_sepToNode (int *sepToNode, int *elemToNode, int *oldBijection, int **nodeBijection, int firstSepElem, int lastSepElem, int dimElem, int *nbSepNodes);

/* =============================================================================
     PRIVATE METHODS
============================================================================= */
void coo_tocsr(int *pairToNode, int *graphIndex, int *graphValue, int nbNodes, int nbPair);

void csr_tocoo(int *pairToNode, int *graphIndex, int *graphValue, int nbNodes, int nbPair);

void graph_to_csr_sym(int *pairToNode, int *graphIndex, int *graphValue, int nbNodes, int nbPairs);

#endif // ALT_PARTITIONING_H
