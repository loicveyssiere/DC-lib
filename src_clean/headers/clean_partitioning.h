#ifndef CLEAN_PARTITIONING_H
#define CLEAN_PARTITIONING_H

#include "DC.h"

void clean_preprocess_partition (int *elemToNode, int nbElem, int dimElem, int nbNodes, int nbPart, int *nodePart);

int clean_get_nb_parts(int nbOfItems);

void clean_graph_to_csr_sym(int *pairToNode, int *graphIndex, int *graphValue, int nbNodes, int nbPairs);

#endif // CLEAN_PARTITIONING_H
