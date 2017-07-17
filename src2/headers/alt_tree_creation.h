#ifndef ALT_TREE_CREATION_H
#define ALT_TREE_CREATION_H

#include "DC.h"

void alt_tree_creation (tree_t &tree, int *elemToNode, int *sepToNode, int *nodePart,
                    int *nodePartSize, int *nodeIndex, int *sepIndex, int globalNbElem, int dimElem, int firstPart,
                    int lastPart, int firstElem, int lastElem, int firstNode,
                    int lastNode, int firstInnerNode, int lastInnerNode,
                    int sepOffset, int sepNodeOffset, int curNode, bool isSep, int depth);

void alt_create_elem_part (int *elemPart, int *nodePart, int *elemToNode, int nbElem,
                       int dimElem, int separator, int offset, int *nbLeftElem,
                       int *nbSepElem);

void alt_create_elem_part2 (int *elemPart, int *nodePart, int *elemToNode, int nbElem,
                      int dimElem, int separator, int offset, int *nbLeftElem,
                      int *nbSepElem);

void alt_create_elem_part3 (int *elemPart, int *nodePart, int *innerNodePart,
                      int *elemToNode, int *nodeIndex, int nbElem, int nbInnerNodes,
                      int dimElem, int separator, int offsetElem, int offsetNode, int firstInnerNode,
                      int *nbLeftElem, int *nbSepElem, int *nbLeftInnerNodes,
                      int *nbSepInnerNodes);

void alt_init_dc_tree (tree_t &tree, int firstElem, int lastElem, int nbSepElem,
                 int firstNode, int lastNode, int firstInnerNode,
                 int lastInnerNode, int nbSepInnerNodes, bool isSep, bool isLeaf, int depth);

void init_intervals (tree_t &tree,
                      int firstElem, int lastElem,
                      int firstNode, int lastNode,
                      int firstInnerNode, int lastInnerNode,
                      bool isSep);

    void alt_tree_creation_new (tree_t &tree, int *elemToNode, int *sepToNode,
                      int *nodePart, int *nodePartSize, int dimElem,
                      int firstPart, int lastPart, int sepOffset, int curNode);

#endif // ALT_TREE_CREATION_H
