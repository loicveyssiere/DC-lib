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

#ifdef CILK
    #include <cilk/cilk.h>
#endif
#include <pthread.h>
#include <limits.h>
#include <cstring>
#include <iostream>

#include "tools.h"
#include "coloring.h"
#include "permutations.h"
#include "partitioning.h"
#include "tree_creation.h"

// Global variables used in order to persist from one call to the library to another:
// the D&C tree and permutations
tree_t *treeHead = nullptr;
int *elemPerm = nullptr, *nodePerm = nullptr;

#ifdef TREE_CREATION

// Mutex to avoid race condition in merge permutations
pthread_mutex_t mergeMutex = PTHREAD_MUTEX_INITIALIZER;

#ifdef MULTITHREADED_COMM

// List containing the last updater of each node, and comm counter mutex
int *nodeOwner = nullptr;
pthread_mutex_t *DCcommMutex = nullptr;

// Initialize D&C tree interface for multithreaded communication
void create_multithreaded_intf (tree_t &tree, int *intfIndex, int *intfNodes,
                                int *intfDstIndex, int *nbDCcomm, int nbIntf,
                                int nbBlocks)
{
    // If there is only one domain, do nothing
    if (nbBlocks < 2) return;

    // Run through all current D&C leaf nodes
    for (int i = 0; i < tree.nbOwnedNodes; i++) {
        int node = tree.ownedNodes[i];

        // Run through all interface nodes
        for (int j = 0; j < nbIntf; j++) {
            for (int k = intfIndex[j]; k < intfIndex[j+1]; k++) {

                // Count the number of nodes on the D&C interface
                if ((node + 1) == intfNodes[k]) {
                    tree.nbIntfNodes++;
                }
            }
        }
    }

    // If there are owned nodes on the interface
    if (tree.nbIntfNodes > 0) {
        int ctr = 0;

        // Allocate the multithreaded interface
        tree.intfIndex = new int [nbIntf+1];
        tree.intfNodes = new int [tree.nbIntfNodes];
        tree.intfDst   = new int [tree.nbIntfNodes];

        // Run through all interfaces
        for (int i = 0; i < nbIntf; i++) {
            tree.intfIndex[i] = ctr;

            // Run through all nodes owned by current D&C leaf
            for (int j = 0; j < tree.nbOwnedNodes; j++) {
                int node = tree.ownedNodes[j];

                // Run through all interface nodes
                for (int k = intfIndex[i]; k < intfIndex[i+1]; k++) {
                    int intfNode = intfNodes[k];

                    // If D&C node is on the interface, add it to the D&C interface
                    if ((node + 1) == intfNode) {
                        tree.intfNodes[ctr] = intfNode;
                        tree.intfDst[ctr]   = intfDstIndex[i] - intfIndex[i] + k;
                        ctr++;
                    }
                    if (ctr == tree.nbIntfNodes) break;
                }
                if (ctr == tree.nbIntfNodes) break;
            }
            // Increment the number of D&C nodes on the current interface
            if (ctr > tree.intfIndex[i]) {
                pthread_mutex_lock (&(DCcommMutex[i]));
                nbDCcomm[i]++;
                pthread_mutex_unlock (&(DCcommMutex[i]));
            }
        }
        tree.intfIndex[nbIntf] = ctr;
    }
}

// Compute the number of nodes owned by current leaf and fill the list
void create_owned_nodes_list (tree_t &tree, int *elemToNode, int dimElem, int curNode)
{
    // Allocate node checker
    int nbNodes = tree.lastNode - tree.firstNode + 1;
    char *nodeChecker = new char [nbNodes] ();

    // Count the number of owned nodes
    int nbOwnedNodes = 0;
    if (tree.isSep) {
        for (int i = tree.firstElem * dimElem; i < (tree.lastElem+1) * dimElem; i++) {
            int node = elemToNode[i] - 1;
            if (nodeChecker[node-tree.firstNode] > 0) continue;
            nodeChecker[node-tree.firstNode]++;
            if (nodeOwner[node] == curNode) nbOwnedNodes++;
        }
    }
    else {
        for (int i = tree.firstNode; i <= tree.lastNode; i++) {
            if (nodeOwner[i] == curNode) nbOwnedNodes++;
        }
    }
    tree.nbOwnedNodes = nbOwnedNodes;

    // Allocate and initialize the list of owned nodes
    if (nbOwnedNodes > 0) {

        // Reset node checker & allocate owned nodes list
        memset (nodeChecker, 0, nbNodes * sizeof (char));
        tree.ownedNodes = new int [nbOwnedNodes];
        int ctr = 0;

        if (tree.isSep) {
            for (int i = tree.firstElem*dimElem; i < (tree.lastElem+1)*dimElem; i++) {
                int node = elemToNode[i] - 1;
                if (nodeChecker[node-tree.firstNode] > 0) continue;
                nodeChecker[node-tree.firstNode]++;

                if (nodeOwner[node] == curNode) {
                    tree.ownedNodes[ctr] = node;
                    ctr++;
                }
                if (ctr == nbOwnedNodes) break;
            }
        }
        else {
            for (int i = tree.firstNode; i <= tree.lastNode; i++) {
                if (nodeOwner[i] == curNode) {
                    tree.ownedNodes[ctr] = i;
                    ctr++;
                }
                if (ctr == nbOwnedNodes) break;
            }
        }
    }
    delete[] nodeChecker;
}

#endif

// Compute the edge interval, the list of nodes owned by each leaf of the D&C tree,
// the interface index for multithreaded communication, and the number of D&C nodes on
// the interface
void tree_finalize (tree_t &tree, int *nodeToNodeRow, int *elemToNode, int *intfIndex,
                    int *intfNodes, int *intfDstIndex, int *nbDCcomm, int dimElem,
                    int nbBlocks, int nbIntf, int curNode, int curLevel, int LRS,
                    ofstream &dcFile)
{
    // If current node is a leaf
    if (tree.left == nullptr && tree.right == nullptr) {

        // Get the first and last edges of the leaf
        //tree.firstEdge = 0;//nodeToNodeRow[tree.firstNode];
        //tree.lastEdge  = 0;//nodeToNodeRow[tree.lastNode+1] - 1;

        #ifdef STATS
            bool hasIntfNode = false;//has_intf_node (tree, intfIndex, intfNodes, nbIntf);
            fill_dc_file_leaves (tree, dcFile, curNode, curLevel, LRS, hasIntfNode);
            count_intf_stats (hasIntfNode);
        #endif

        #ifdef MULTITHREADED_COMM
            // Create the list of nodes last updated by current D&C leaf
            create_owned_nodes_list (tree, elemToNode, dimElem, curNode);

            // Create the interface of current D&C leaf
            create_multithreaded_intf (tree, intfIndex, intfNodes, intfDstIndex,
                                       nbDCcomm, nbIntf, nbBlocks);
        #endif
    }
    else {
        #ifdef STATS
            bool hasIntfNode = false;//has_intf_node (tree, intfIndex, intfNodes, nbIntf);
            fill_dc_file_nodes (tree, dcFile, curNode, curLevel, hasIntfNode);
        #endif

        // Left & right recursion
        #ifdef OMP
            #pragma omp task default (shared)
        #elif CILK
            cilk_spawn
        #endif
        tree_finalize (*tree.right, nodeToNodeRow, elemToNode, intfIndex, intfNodes,
                       intfDstIndex, nbDCcomm, dimElem, nbBlocks, nbIntf, 3*curNode+2,
                       curLevel+1, 2, dcFile);
        #ifdef OMP
            #pragma omp task default (shared)
        #endif
        tree_finalize (*tree.left, nodeToNodeRow, elemToNode, intfIndex, intfNodes,
                       intfDstIndex, nbDCcomm, dimElem, nbBlocks, nbIntf, 3*curNode+1,
                       curLevel+1, 1, dcFile);

        // Synchronization
        #ifdef OMP
            #pragma omp taskwait
        #elif CILK
            cilk_sync;
        #endif

        // Separator recursion, if it is not empty
        if (tree.sep != nullptr) {
            tree_finalize (*tree.sep, nodeToNodeRow, elemToNode, intfIndex, intfNodes,
                           intfDstIndex, nbDCcomm, dimElem, nbBlocks, nbIntf,
                           3*curNode+3, curLevel+1, 3, dcFile);
        }
    }
}

// Wrapper used to get the root of the D&C tree before calling the real tree finalize
void DC_finalize_tree (int *nodeToNodeRow, int *elemToNode, int *intfIndex,
                       int *intfNodes, int *intfDstIndex, int *nbDCcomm, int nbElem,
                       int dimElem, int nbBlocks, int nbIntf, int rank)
{
    // Create D&C tree dot file
    #ifdef STATS
        string fileName = "dcTree_" + to_string ((long long)rank) + "_" +
                           to_string ((long long)MAX_ELEM_PER_PART) + ".dot";
        ofstream dcFile (fileName, ios::out | ios::trunc);
        init_dc_file (dcFile);
    #else
        ofstream dcFile;
    #endif

    #ifdef MULTITHREADED_COMM
        DCcommMutex = new pthread_mutex_t [nbIntf];
        for (int i = 0; i < nbIntf; i++) {
            DCcommMutex[i] = PTHREAD_MUTEX_INITIALIZER;
        }
    #endif

    // Compute the edge interval, the list of nodes owned by each leaf of the D&C tree,
    // and the interface index for multithreaded communication
    #ifdef OMP
        #pragma omp parallel
        #pragma omp single nowait
    #endif
    tree_finalize (*treeHead, nodeToNodeRow, elemToNode, intfIndex, intfNodes,
                   intfDstIndex, nbDCcomm, dimElem, nbBlocks, nbIntf, 0, 0, -1,dcFile);

    #ifdef STATS
        close_dc_file (dcFile);
        store_intf_stats (nbElem, rank);
    #endif

    #ifdef MULTITHREADED_COMM
        delete[] DCcommMutex, delete[] nodeOwner;
    #endif
}

// Initialize the content of D&C tree nodes
void init_dc_tree (tree_t &tree, int firstElem, int lastElem, int nbSepElem,
                   int firstNode, int lastNode, bool isSep, bool isLeaf)
{
    tree.ownedNodes   = nullptr;
    tree.intfIndex    = nullptr;
    tree.intfNodes    = nullptr;
    tree.intfDst      = nullptr;
    tree.nbOwnedNodes = 0;
    tree.nbIntfNodes  = 0;
    tree.firstElem    = firstElem;
    tree.lastElem     = lastElem - nbSepElem;
    tree.lastSep      = lastElem;
    tree.firstNode    = firstNode;
    tree.lastNode     = lastNode;
    tree.firstEdge    = -1;
    tree.lastEdge     = -1;
    tree.vecOffset    = 0;
    tree.isSep        = isSep;
    tree.left         = nullptr;
    tree.right        = nullptr;
    tree.sep          = nullptr;

    if (isLeaf == false) {
        tree.left  = new tree_t;
        tree.right = new tree_t;
        if (nbSepElem > 0) {
            tree.sep = new tree_t;
        }
    }
}

#ifdef MULTITHREADED_COMM

// Set the last updater of each node. The owners are the closer leaves to the root.
void fill_node_owner (int *elemToNode, int firstElem, int lastElem, int dimElem,
                      int firstNode, int lastNode, int curNode, bool isSep)
{
    // MUTEX on nodeOwner ???
    // If current D&C leaf is a separator, it has no node interval
    if (isSep) {
        for (int i = firstElem * dimElem; i < (lastElem+1) * dimElem; i++) {
            int node = nodePerm[elemToNode[i]];
            nodeOwner[node] = curNode;
        }
    }
    else {
        for (int i = firstNode; i <= lastNode; i++) {
            nodeOwner[i] = curNode;
        }
    }
}

// Initialize the owner of the nodes to MAX_INT (lower is owner)
void init_node_owner (int nbNodes)
{
    nodeOwner = new int [nbNodes];
    #ifdef OMP
        #pragma omp parallel for
        for (int i = 0; i < nbNodes; i++) {
    #elif CILK
        cilk_for (int i = 0; i < nbNodes; i++) {
    #endif
        nodeOwner[i] = INT_MAX;
    }
}

#endif

// Create element partition & count left & separator elements
void create_elem_part (int *elemPart, int *nodePart, int *elemToNode, int nbElem,
                       int dimElem, int separator, int offset, int *nbLeftElem,
                       int *nbSepElem)
{
	for (int i = 0; i < nbElem; i++) {
		int node, leftCtr = 0, rightCtr = 0;
		for (int j = 0; j < dimElem; j++) {
			node = elemToNode[(i+offset)*dimElem+j];
			if (nodePart[node] <= separator) leftCtr++;
			else                             rightCtr++;
		}
		if (leftCtr == dimElem) {
			elemPart[i] = 0;
			(*nbLeftElem)++;
		}
		else if (rightCtr == dimElem) {
			elemPart[i] = 1;
		}
		else {
			elemPart[i] = 2;
			(*nbSepElem)++;
		}
	}
}

// Create the D&C tree and the element permutation, and compute the intervals of nodes
// and elements at each node of the tree
void tree_creation (tree_t &tree, int *elemToNode, int *sepToNode, int *nodePart,
                    int *nodePartSize, int globalNbElem, int dimElem, int firstPart,
                    int lastPart, int firstElem, int lastElem, int firstNode,
                    int lastNode, int sepOffset, int curNode, bool isSep)
{
    // If current node is a leaf
    int nbPart = lastPart - firstPart + 1;
    int localNbElem = lastElem - firstElem + 1;
    if (nbPart < 2 || localNbElem <= MAX_ELEM_PER_PART) {

        #ifdef MULTITHREADED_COMM
            // Set the node accessed by current D&C leaf
            fill_node_owner (elemToNode, firstElem, lastElem, dimElem, firstNode,
                             lastNode, curNode, isSep);
        #endif

        // Initialize the leaf
        init_dc_tree (tree, firstElem, lastElem, 0, firstNode, lastNode, isSep, true);

        // End of recursion
        return;
    }

    // Else, prepare next left, right & separator recursion
    int nbLeftElem = 0, nbSepElem = 0, nbLeftNodes = 0, nbRightNodes = 0;
    int separator = firstPart + (lastPart - firstPart) / 2;

    // Count the number of left & right nodes
    if (!isSep) {
        for (int i = firstPart; i <= separator; i++) {
            nbLeftNodes += nodePartSize[i];
        }
        nbRightNodes = (lastNode - firstNode + 1) - nbLeftNodes;
    }

    // Create local element partition & count left & separator elements
    int *localElemPart = new int [localNbElem];
    if (isSep) {
        create_elem_part (localElemPart, nodePart, sepToNode, localNbElem, dimElem,
                          separator, sepOffset, &nbLeftElem, &nbSepElem);
    }
    else {
        create_elem_part (localElemPart, nodePart, elemToNode, localNbElem, dimElem,
                          separator, firstElem, &nbLeftElem, &nbSepElem);
    }

    // Create local element permutation
    int *localElemPerm = new int [localNbElem];
    DC_create_permutation (localElemPerm, localElemPart, localNbElem, 3);
    delete[] localElemPart;

    // Execution is correct without mutex although cilkscreen detects a race condition
    pthread_mutex_lock (&mergeMutex);
    // Apply local element permutation to global element permutation
    merge_permutations (localElemPerm, globalNbElem, localNbElem, firstElem, lastElem);
    pthread_mutex_unlock (&mergeMutex);

    // Permute elemToNode and sepToNode with local element permutation
    DC_permute_int_2d_array (elemToNode, localElemPerm, localNbElem, dimElem,
                             firstElem);
    if (isSep) {
        DC_permute_int_2d_array (sepToNode, localElemPerm, localNbElem, dimElem,
                                 sepOffset);
    }
    delete[] localElemPerm;

    // Initialize current node
    init_dc_tree (tree, firstElem, lastElem, nbSepElem, firstNode, lastNode, isSep,
                  false);

    // Left & right recursion
    #ifdef OMP
        #pragma omp task default(shared)
    #elif CILK
        cilk_spawn
    #endif
    tree_creation (*tree.right, elemToNode, sepToNode, nodePart, nodePartSize,
                   globalNbElem, dimElem, separator+1, lastPart, firstElem+nbLeftElem,
                   lastElem-nbSepElem, firstNode+nbLeftNodes, lastNode, sepOffset+
                   nbLeftElem, 3*curNode+2, isSep);
    #ifdef OMP
        #pragma omp task default(shared)
    #endif
    tree_creation (*tree.left, elemToNode, sepToNode, nodePart, nodePartSize,
                   globalNbElem, dimElem, firstPart, separator, firstElem, firstElem+
                   nbLeftElem-1, firstNode, lastNode-nbRightNodes, sepOffset,
                   3*curNode+1, isSep);

    // Synchronization
    #ifdef OMP
        #pragma omp taskwait
    #elif CILK
        cilk_sync;
    #endif

    // D&C partitioning of separator elements
    if (nbSepElem > 0) {
        sep_partitioning (*tree.sep, elemToNode, globalNbElem, dimElem, lastElem-
                          nbSepElem+1, lastElem, firstNode, lastNode, 3*curNode+3);
    }
}

// Create the D&C tree and the permutations
void DC_create_tree (int *elemToNode, int nbElem, int dimElem, int nbNodes)
{
    // Allocate the D&C tree & the permutation functions
    treeHead = new tree_t;
    elemPerm = new int [nbElem];
    nodePerm = new int [nbNodes];

    #ifdef MULTITHREADED_COMM
        // Initialize the node owner array
        init_node_owner (nbNodes);
    #endif

    // Create the D&C tree & the permutation functions
    partitioning (elemToNode, nbElem, dimElem, nbNodes);

    // Vectorial version with coloring of the leaves of the D&C tree
    #ifdef DC_VEC
        coloring (elemToNode, nbElem, dimElem, nbNodes);
    #endif
}

#endif
