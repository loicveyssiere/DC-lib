/*  Copyright 2014 - UVSQ
    Authors list: Loïc Thébault, Eric Petit

    This file is part of the DC-lib.

    DC-lib is free software: you can redistribute it and/or modify it under the
    terms of the GNU Lesser General Public License as published by the Free Software
    Foundation, either version 3 of the License, or (at your option) any later version.

    DC-lib is distributed in the hope that it will be useful, but WITHOUT ANY
    WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
    PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License along with
    the DC-lib. If not, see <http://www.gnu.org/licenses/>. */

#ifndef DC_H
#define DC_H

#include <stdint.h>
#include <string>
#include <fstream>

#define MAX_ELEM_NEIGHBORS 250
#ifndef MAX_ELEM_PER_PART
#define MAX_ELEM_PER_PART 4000
#endif // MAX_ELEM_PER_PART

#define DC_SUM      0
#define DC_MIN      1
#define DC_MAX      2
#define DC_PROD     3
#define DC_LAND     4
#define DC_LOR      5
#define DC_BAND     6
#define DC_BOR      7

#define DC_INT      0
#define DC_DOUBLE   1

struct DCreturnArgs_t {
    double *val;
    int size;
    int op;
    DCreturnArgs_t(int size_, int op_) {
        val = new double[size_]();
        size = size_;
        op = op_;
    }
    ~DCreturnArgs_t() {
        delete[] val;
    }
};

// D&C tree structure
typedef struct tree_s {
    int *ownedNodes, *intfIndex, *intfNodes, *intfDst;
    int nbOwnedNodes, nbIntfNodes,
        firstElem, lastElem, lastSep,
        firstNode, lastNode,
        firstInnerNode, lastInnerNode,
        firstEdge, lastEdge,
        vecOffset;
    int depth;
    bool isSep;
    bool isLeaf;
    struct tree_s *left, *right, *sep;
} tree_t;

// D&C arguments structure
typedef struct DCargs_s {
    #ifdef MULTITHREADED_COMM
        int *ownedNodes;
        int nbOwnedNodes;
    #endif
    int firstElem, lastElem,
        firstNode, lastNode,
        firstInnerNode, lastInnerNode,
        firstEdge, lastEdge,
        isSep, isLeaf, depth;
} DCargs_t;

// D&C arguments structure for multithreaded comm
typedef struct DCcommArgs_s {
    int *intfIndex, *intfNodes, *intfDst, *intfOffset;
} DCcommArgs_t;

typedef struct {
    int *index, *value;
} index_t;
typedef struct {
    int list[MAX_ELEM_NEIGHBORS];
    int size;
} list_t;
typedef struct {
    int elem, node;
} couple_t;

// Compute the average time/cycles of each (start, stop) intervals until reset call
class DC_timer
{
    public:
        // Constructor
        DC_timer ();

        // Get time of day
        double get_avg_time ();
        void reset_time ();
        void stop_time ();
        void start_time ();

        // RDTSC
        uint64_t get_avg_cycles ();
        void reset_cycles ();
        void stop_cycles ();
        void start_cycles ();

    private:
        // Get time of day
        double startTime, avgTime;
        int timeCtr;

        // RDTSC
        uint64_t startCycles, avgCycles;
        int cyclesCtr;
};

class Tree_topo {
public:
    Tree_topo(std::string stats_file, int n_extra);
    ~Tree_topo();
    void add(DCargs_t *DCargs, int *extra);
private:
    std::ofstream stream;
    std::string stats;
    int n_extra;
};

// Get time of day
double DC_get_time ();

// RDTSC
uint64_t DC_get_cycles ();

// Wrapper used to get the root of the D&C tree before calling the real tree traversal
void DC_tree_traversal (void (*userSeqFct)  (void *, DCargs_t *),
                        void (*userVecFct)  (void *, DCargs_t *),
                        void (*userCommFct) (void *, DCcommArgs_t *),
                        void *userArgs, void *userCommArgs);

void DC_alt_tree_traversal (void (*userSeqFct)  (void *, DCargs_t *, DCreturnArgs_t *),
                        void (*userVecFct)  (void *, DCargs_t *),
                        void (*userCommFct) (void *, DCcommArgs_t *),
                        void *userArgs, void *userCommArgs, DCreturnArgs_t *returnArgs);

// Create element to element array from element to node and node to element
// Two elements are neighbors if they share a node
void DC_create_elemToElem (list_t *elemToElem, index_t &nodeToElem, int *elemToNode,
                           int firstElem, int lastElem, int dimElem);

// Create node to element structure from element to node
void DC_create_nodeToElem (index_t &nodeToElem, int *elemToNode, int nbElem,
                           int dimElem, int nbNodes);

// Permute "tab" 2D array of double using node permutation
void DC_permute_double_2d_array (double *tab, int nbItem, int dimItem);

// Permute "tab" 2D array of int using "perm"
void DC_permute_int_2d_array (int *tab, int *perm, int nbItem, int dimItem,
                              int offset);

// Permute "tab" 1D array of int using node permutation
void DC_permute_int_1d_array (int *tab, int size);

// Renumber "tab" array of int using node permutation
void DC_renumber_int_array (int *tab, int size, bool isFortran);

// Create permutation array from partition array
void DC_create_permutation (int *perm, int *part, int size, int nbPart);

// Read the D&C tree and the permutation functions
void DC_read_tree (std::string &treePath, int nbElem, int nbNodes, int nbIntf,
                   int *nbMaxComm);

// Store the D&C tree and the permutation functions to a binary file
void DC_store_tree (std::string &treePath, int nbElem, int nbNodes, int nbIntf,
                    int nbMaxComm);

// Wrapper used to get the root of the D&C tree before calling the real tree finalize
void DC_finalize_tree (int *nodeToNodeRow, int *elemToNode, int *intfIndex,
                       int *intfNodes, int *intfDstOffsets, int *nbDCcomm, int nbElem,
                       int dimElem, int nbBlocks, int nbIntf, int rank);

// Create the D&C tree and the permutations
void DC_create_tree (int *elemToNode, int nbElem, int dimElem, int nbNodes);

/* =============================================================================
    ADD-ON
============================================================================= */
#define MAX_NODE_PER_PART 100

void DC_alt_create_tree(int *componentToNode, int nbComponent, int dimComponent, int nbNodes);

void DC_alt_permute_double_2d_array (double *tab, int *perm, int nbItem, int dimItem, int offset);

int DC_get_max_elem_per_part();

/* =============================================================================
    END ADD-ON
============================================================================= */

#endif
