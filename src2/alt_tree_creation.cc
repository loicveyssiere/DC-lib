#include <cstring>
#include <cmath>
#ifdef CILK
    #include <cilk/cilk.h>
    #include <cilk/cilk_api.h>
#endif
#include <pthread.h>
#include <metis.h>
#include <assert.h>
#include <algorithm>

#include "permutations.h"
#include "partitioning.h"
#include "tools.h"
#include "tree_creation.h"
// Only include tree_creation to access to global variables
// tree_creation need tools unfortunately ...

#include "alt_permutations.h"
#include "alt_partitioning.h"
#include "alt_tree_creation.h"

extern tree_t *treeHead;
extern int *elemPerm, *nodePerm;
int *innerNodePerm;
int globalNbElem, globalNbNodes;

#ifdef TREE_CREATION

// Mutex to avoid race condition in merge permutations
extern pthread_mutex_t mergeMutex;

void print_node_stat(int *elemToNode, int *nodePart, int nbElem, int nbNodes) {

    int *res = new int[nbNodes]();

    for (int i = 0; i < nbElem; i++) {
        int n1 = elemToNode[2*i], n2 = elemToNode[2*i+1];
        if (nodePart[n1] != nodePart[n2]) {
            res[n1] = -1;
            res[n2] = -1;
            // if (nodePart[n1] < nodePart[n2]) {
            //     res[n1] = -1;
            // } else {
            //     res[n2] = -1;
            // }
        }
    }

    int out = 0;
    for (int i = 0; i < nbNodes; i++) {
        if (res[i] == -1) {
            out += 1;
        }
    }
    printf("NUMBER OF NODES IN SEP: %d\n", out);
}

void DC_alt_create_tree(int *componentToNode, int nbComponent, int dimComponent, int nbNodes) {

    // Allocate the D&C tree & the permutation functions
    elemPerm = new int [nbComponent];
    nodePerm = new int [nbNodes];
    innerNodePerm = new int[nbNodes];
    int *nodeIndex = new int[nbNodes];

     printf("Number of workers: %d",  __cilkrts_get_nworkers());

    #ifdef MULTITHREADED_COMM
        // Initialize the node owner array
        init_node_owner (nbNodes);
    #endif

    // Compute node partition
    int *nodePart   = new int [nbNodes];
    //int nbPart = ceil (nbNodes / (double) MAX_NODE_PER_PART);
    int nbPart = ceil (nbComponent / (double) MAX_ELEM_PER_PART);
    compute_partition(componentToNode, nbComponent, dimComponent, nbNodes, nbPart, nodePart);

    print_node_stat(componentToNode, nodePart, nbComponent, nbNodes);

    // Create node permutation from node partition
    DC_create_permutation (nodePerm, nodePart, nbNodes, nbPart);

    // Compute the number of nodes per partition
    int *nodePartSize = new int [nbPart] ();
    for (int i = 0; i < nbNodes; i++) {
        nodePartSize[nodePart[i]]++;
    }

    // Initialize the global node permutation
    globalNbNodes = nbNodes;
    #ifdef OMP
        #pragma omp parallel for
        for (int i = 0; i < nbNodes; i++) {
    #elif CILK
        cilk_for (int i = 0; i < nbNodes; i++) {
    #endif
        innerNodePerm[i] = i;
        nodeIndex[i] = i;
    }

    // Initialize the global element permutation
    globalNbElem = nbComponent;
    #ifdef OMP
        #pragma omp parallel for
        for (int i = 0; i < nbComponent; i++) {
    #elif CILK
        cilk_for (int i = 0; i < nbComponent; i++) {
    #endif
        elemPerm[i] = i;
    }

    treeHead = new tree_t;
    //init_intervals(*treeHead, 0, nbComponent-1, 0, nbNodes-1, 0, nbNodes-1, false); // added

    // Create D&C tree
    #ifdef OMP
        #pragma omp parallel
        #pragma omp single nowait
    #endif
    // alt_tree_creation_new (*treeHead, componentToNode, nullptr, nodePart,
    //     nodePartSize, dimComponent, 0, nbPart-1, 0, 0);
    alt_tree_creation (*treeHead, componentToNode, nullptr, nodePart, nodePartSize, nodeIndex, nullptr, nbComponent,
                   dimComponent, 0, nbPart-1, 0, nbComponent-1, 0, nbNodes-1, 0, nbNodes-1, 0, 0, 0, false, 0);
    delete[] nodePartSize, delete[] nodePart, delete[] nodeIndex;

    for (int i = 0; i < nbNodes; i++) {
        nodePerm[i] = innerNodePerm[i];
    }

    // Vectorial version with coloring of the leaves of the D&C tree
    #ifdef DC_VEC
        coloring (componentToNode, nbComponent, dimComponent, nbNodes);
    #endif
}

void alt_create_elem_part (int *elemPart, int *nodePart, int *elemToNode, int nbElem,
                       int dimElem, int separator, int offset, int *nbLeftElem,
                       int *nbSepElem) {

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

void alt_create_elem_part2 (int *elemPart, int *nodePart, int *elemToNode, int nbElem,
                       int dimElem, int separator, int offset, int *nbLeftElem,
                       int *nbSepElem) {

    int *maskNodePart = new int[globalNbNodes]();

    for (int i = 0; i < nbElem; i++) {
        int node, leftCtr = 0, rightCtr = 0;
        for (int j = 0; j < dimElem; j++) {
			node = elemToNode[(i+offset)*dimElem+j];
            if (nodePart[node] <= separator) {
                leftCtr++;
            } else {
                rightCtr++;
            }
        }
        if (leftCtr != dimElem && rightCtr != dimElem) {

            // Version with 2 lines
            for (int j = 0; j < dimElem; j++) {
                node = elemToNode[(i+offset)*dimElem+j];
                maskNodePart[node] = -1;
            }

            // Version with 1 line
            // int min = globalNbNodes;
            // int min_color = globalNbNodes+1;
            // for (int j = 0; j < dimElem; j++) {
            //     node = elemToNode[(i+offset)*dimElem+j];
            //     if (nodePart[node] < min_color) {
            //         min_color = nodePart[node];
            //         min = node;
            //     }
            // }
            // maskNodePart[min] = -1;
        }
    }


    int nbRight = 0;
	for (int i = 0; i < nbElem; i++) {
		int node, leftCtr = 0, rightCtr = 0;
        int globalLeftCtr = 0, globalRightCtr = 0, targetDim = dimElem;

		for (int j = 0; j < dimElem; j++) {
			node = elemToNode[(i+offset)*dimElem+j];

            if (nodePart[node] <= separator) {
                globalLeftCtr++;
                if (maskNodePart[node] != -1) {
                    leftCtr++;
                }
            } else {
                globalRightCtr++;
                if (maskNodePart[node] != -1) {
                    rightCtr++;
                }
            }
            if (maskNodePart[node] == -1) {
                targetDim--;
            }
		}

        // Compute the real binary partition ...
        if (globalLeftCtr == targetDim) {
            elemPart[i] = 0;
			(*nbLeftElem)++;
        } else if (globalRightCtr == targetDim) {
            elemPart[i] = 1;
            nbRight++;
        } else {
            elemPart[i] = 2;
			(*nbSepElem)++;
        }

	}
    assert(nbRight + (*nbSepElem) + (*nbLeftElem) == nbElem);
    delete[] maskNodePart;
}

void alt_create_elem_part3 (int *elemPart, int *nodePart, int *innerNodePart,
                        int *elemToNode, int *nodeIndex, int nbElem, int nbInnerNodes,
                        int dimElem, int separator, int offsetElem, int offsetNode, int firstInnerNode,
                        int *nbLeftElem, int *nbSepElem, int *nbLeftInnerNodes,
                        int *nbSepInnerNodes) {

    int *maskNodePart = new int[globalNbNodes];
    for (int i = 0; i < globalNbNodes; i++) {maskNodePart[i] = globalNbNodes;}

    for (int j = 0; j < nbInnerNodes; j++) {

        int node = nodeIndex[j + offsetNode];
        if (nodePart[node] <= separator) {
            innerNodePart[j] = 0;
        } else {
            innerNodePart[j] = 1;
        }
        maskNodePart[node] = j;
    }

    for (int i = 0; i < nbElem; i++) {
        int node, innerNode, leftCtr = 0, rightCtr = 0;
        for (int j = 0; j < dimElem; j++) {
            node = elemToNode[(i+offsetElem)*dimElem+j];
            if (nodePart[node] <= separator) {
                leftCtr++;
            } else {
                rightCtr++;
            }
        }

        // Put node in the separator
        if (leftCtr != dimElem && rightCtr != dimElem) {

            // Version with 2 lines
            for (int j = 0; j < dimElem; j++) {
                node = elemToNode[(i+offsetElem)*dimElem+j];
                innerNode = maskNodePart[node];
                if (innerNode >=0 && innerNode < nbInnerNodes) {
                    innerNodePart[innerNode] = 2;
                } else {
                    // do nothing
                }
            }

            // Version with 1 line
            // int node_min = globalNbNodes, color_min = globalNbNodes;
            // for (int j = 0; j < dimElem; j++) {
            //     node = elemToNode[(i+offsetElem)*dimElem+j];
            //     if (nodePart[node] < color_min) {
            //         color_min = nodePart[node];
            //         node_min = node;
            //     }
            // }
            // innerNode = maskNodePart[node_min];
            // if (innerNode >= 0 && innerNode < nbInnerNodes) {
            //     innerNodePart[innerNode] = 2;
            // }
            // else if (innerNode < 0) {
            //     assert(false);
            // }

        }
    }

    // Compute number of inner nodes in each part
    int nbRight = 0;
    for (int j = 0; j < nbInnerNodes; j++) {
        //printf("%d\n", j);
        if (innerNodePart[j] == 0) {
            (*nbLeftInnerNodes)++;
        } else if(innerNodePart[j] == 1) {
            nbRight++;
        } else if(innerNodePart[j] == 2) {
            (*nbSepInnerNodes)++;
        } else {
            printf("%d %d\n", j, innerNodePart[j]);
            assert(false);
        }
    }
    assert(nbRight + (*nbSepInnerNodes) + (*nbLeftInnerNodes) == nbInnerNodes);

    nbRight = 0;
    for (int i = 0; i < nbElem; i++) {
        int node, innerNode;
        int leftCtr = 0, rightCtr = 0, sepCtr = 0;
        int globalLeftCtr = 0, globalRightCtr = 0;

        for (int j = 0; j < dimElem; j++) {
            node = elemToNode[(i+offsetElem)*dimElem+j];
            //innerNode = innerNodePerm[node + firstInnerNode] - offsetNode;
            innerNode = maskNodePart[node];

            if (innerNode >= 0) {
                if (innerNode >= nbInnerNodes) {
                    if (nodePart[node] <= separator) {
                        globalLeftCtr++;
                    } else {
                        globalRightCtr++;
                    }
                } else {
                    if (innerNodePart[innerNode] == 0) {
                        leftCtr++;
                    } else if (innerNodePart[innerNode] == 1) {
                        rightCtr++;
                    } else if (innerNodePart[innerNode] == 2) {
                        sepCtr++;
                    } else {
                        assert(false);
                    }
                }
            } else {
                assert(false);
            }
        }

        assert(globalLeftCtr + globalRightCtr + leftCtr + rightCtr + sepCtr == dimElem);

        if (leftCtr > 0) {
            elemPart[i] = 0;
            (*nbLeftElem)++;
        } else if (rightCtr > 0) {
            elemPart[i] = 1;
            nbRight++;
        } else if (sepCtr > 0) {
            elemPart[i] = 2;
            (*nbSepElem)++;
        } else if (globalLeftCtr == dimElem) {
            elemPart[i] = 0;
            (*nbLeftElem)++;
        } else if (globalRightCtr == dimElem) {
            elemPart[i] = 1;
            nbRight++;
        } else {
            elemPart[i] = 0;
            (*nbLeftElem)++;
            //assert(false);
        }


        //Compute the real binary partition ...
        // if (leftCtr != 0 && rightCtr != 0) {
        //     assert(false);
        // }
        // else if (leftCtr != 0 && targetDim == 0) {
        //     assert(false);
        // } else if (rightCtr != 0 && targetDim == 0) {
        //     assert(false);
        // } else if (leftCtr != 0 && rightCtr == 0) {
        //     elemPart[i] = 0;
        //     (*nbLeftElem)++;
        // } else if (leftCtr == 0 && rightCtr != 0) {
        //     elemPart[i] = 1;
        //     nbRight++;
        // } else if (targetDim == 0) {
        //     elemPart[i] = 2;
        //     (*nbSepElem)++;
        // }
        //  else {
        //     assert(false);
        // }
    }

    assert(nbRight + (*nbSepElem) + (*nbLeftElem) == nbElem);
    delete[] maskNodePart;

}

void select_partition ( int *elemPart, int *innerNodePart, int *nodePart, int *elemToNode,
                        int dimElem, int nbElem, int nbNodes, int separator,
                        int offsetElem, int offsetNode, // offset
                        int *nbLeftElem, int *nbSepElem, int *nbRightElem, // Elements
                        int *nbLeftNodes, int *nbSepNodes, int *nbRightNodes // Nodes
                    ) {

    printf("Enter in select_partition ...\n");
    // init inner node partition

	for (int i = 0; i < nbElem; i++) {
		int node, index, leftCtr = 0, rightCtr = 0;
		for (int j = 0; j < dimElem; j++) {
			node = elemToNode[(i+offsetElem)*dimElem+j];
			if (nodePart[node] <= separator) leftCtr++;
			else                             rightCtr++;

		}

        // For elements
		if (leftCtr == dimElem) {
			elemPart[i] = 0;
			(*nbLeftElem)++;
		}
		else if (rightCtr == dimElem) {
			elemPart[i] = 1;
            (*nbRightElem)++;
		}
		else {
			elemPart[i] = 2;
			(*nbSepElem)++;
		}

        for (int j = 0; j < dimElem; j++) {
            node = elemToNode[(i+offsetElem)*dimElem+j];
            index = nodePerm[node];
            assert(index >= 0);
            assert(index < nbNodes);
            if (leftCtr == dimElem) {
                innerNodePart[index] = std::max(innerNodePart[index], 0);
            } else if (rightCtr == dimElem) {
                innerNodePart[index] = std::max(innerNodePart[index], 1);
            } else {
                innerNodePart[index] = std::max(innerNodePart[index], 2);
            }
		}
	}

    // compute number of nodes of each partitions
    int add = 0;
    for (int i = 0; i < nbNodes; i++) {
        if (innerNodePart[i] == 0) {
            (*nbLeftNodes)++;
        } else if (innerNodePart[i] == 1) {
            (*nbRightNodes)++;
        } else if (innerNodePart[i] == 2) {
            (*nbSepNodes)++;
        } else {
            add++;
        }
    }
    printf("%d\n", add);
    assert(add == 0);
    assert(*nbLeftElem + *nbRightElem + *nbSepElem == nbElem);
    assert(*nbLeftNodes + *nbRightNodes + *nbSepNodes == nbNodes);
}

void alt_init_dc_tree (tree_t &tree, int firstElem, int lastElem, int nbSepElem,
                   int firstNode, int lastNode, int firstInnerNode,
                   int lastInnerNode, int nbSepInnerNodes, bool isSep, bool isLeaf) {
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
    tree.firstInnerNode = firstInnerNode;
    tree.lastInnerNode = lastInnerNode;
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
        if (nbSepElem > 0 || nbSepInnerNodes > 0) {
            tree.sep = new tree_t;
        } else {
            //assert(false);
        }
    }
}

void alt_tree_creation (tree_t &tree, int *elemToNode, int *sepToNode, int *nodePart,
                    int *nodePartSize, int *nodeIndex, int *sepIndex, int globalNbElem, int dimElem, int firstPart,
                    int lastPart, int firstElem, int lastElem, int firstNode,
                    int lastNode, int firstInnerNode, int lastInnerNode,
                    int sepOffset, int sepNodeOffset, int curNode, bool isSep, long depth) {

    printf("Enter in tree_creation ... %d - %d // %d - %d and depth: %ld\n", firstElem, lastElem, firstInnerNode, lastInnerNode, depth);
    //assert(sepOffset == firstElem);

    // If current node is a leaf
    int nbPart = lastPart - firstPart + 1;
    int localNbElem = lastElem - firstElem + 1;
    int localNbNodes = lastNode - firstNode + 1;
    int localNbInnerNodes = lastInnerNode - firstInnerNode + 1;
    assert(localNbInnerNodes >= 0);
    //if (nbPart < 2 || localNbInnerNodes <= MAX_NODE_PER_PART) {
    if (nbPart < 2 || localNbElem <= MAX_ELEM_PER_PART) {

        #ifdef MULTITHREADED_COMM
            // Set the node accessed by current D&C leaf
            fill_node_owner (elemToNode, firstElem, lastElem, dimElem, firstNode,
                             lastNode, curNode, isSep);
        #endif

        // Initialize the leaf
        alt_init_dc_tree (tree, firstElem, lastElem, 0, firstNode, lastNode, firstInnerNode, lastInnerNode, 0, isSep, true);

        // End of recursion
        return;
    }

    // Else, prepare next left, right & separator recursion
    int nbLeftElem = 0, nbSepElem = 0, nbLeftNodes = 0, nbRightNodes = 0;
    int nbLeftInnerNodes = 0, nbSepInnerNodes = 0;
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
    int *localInnerNodePart = new int [localNbInnerNodes];
    for (int i = 0; i < localNbInnerNodes; i++) {localInnerNodePart[i] = -1;}
    if (isSep) {
        // alt_create_elem_part2 (localElemPart, nodePart, sepToNode, localNbElem, dimElem,
        //                   separator, sepOffset, &nbLeftElem, &nbSepElem);
        // for (int i = 0; i < localNbInnerNodes; i++) {printf("%d ", nodeIndex[i]);}
        alt_create_elem_part3 (localElemPart, nodePart, localInnerNodePart,
                              sepToNode, sepIndex, localNbElem, localNbInnerNodes,
                              dimElem, separator, sepOffset, sepNodeOffset, firstInnerNode, &nbLeftElem, &nbSepElem,
                              &nbLeftInnerNodes, &nbSepInnerNodes);
    }
    else {
        // alt_create_elem_part2 (localElemPart, nodePart, elemToNode, localNbElem, dimElem,
        //                   separator, firstElem, &nbLeftElem, &nbSepElem);
        alt_create_elem_part3 (localElemPart,nodePart, localInnerNodePart,
                            elemToNode, nodeIndex, localNbElem, localNbInnerNodes, dimElem,
                            separator, firstElem, firstInnerNode, 0, &nbLeftElem, &nbSepElem,
                            &nbLeftInnerNodes, &nbSepInnerNodes);
    }

    // Create local element permutation
    int *localElemPerm = new int [localNbElem];
    int *localInnerNodePerm = new int[localNbInnerNodes];
    DC_create_permutation (localElemPerm, localElemPart, localNbElem, 3);
    DC_create_permutation (localInnerNodePerm, localInnerNodePart, localNbInnerNodes, 3);
    delete[] localElemPart;
    delete[] localInnerNodePart;


    // Execution is correct without mutex although cilkscreen detects a race condition
    pthread_mutex_lock (&mergeMutex);
    // Apply local element permutation to global element permutation
    alt_merge_permutations (elemPerm, localElemPerm, globalNbElem, localNbElem, firstElem, lastElem);
    alt_merge_permutations (innerNodePerm, localInnerNodePerm, globalNbNodes, localNbInnerNodes, firstInnerNode, lastInnerNode);
    pthread_mutex_unlock (&mergeMutex);

    // Permute elemToNode and sepToNode with local element permutation
    DC_permute_int_2d_array (elemToNode, localElemPerm, localNbElem, dimElem, firstElem);
    DC_permute_int_2d_array (nodeIndex, localInnerNodePerm, localNbInnerNodes, 1, firstInnerNode);
    if (isSep) {
        DC_permute_int_2d_array (sepToNode, localElemPerm, localNbElem, dimElem, sepOffset);
        DC_permute_int_2d_array (sepIndex, localInnerNodePerm, localNbInnerNodes, 1, sepNodeOffset);
    }
    delete[] localElemPerm;
    delete[] localInnerNodePerm;

    // Initialize current node
    alt_init_dc_tree (tree, firstElem, lastElem, nbSepElem, firstNode, lastNode, firstInnerNode, lastInnerNode, nbSepInnerNodes, isSep,
                  false);

    // Left & right recursion
    #ifdef OMP
        #pragma omp task default(shared)
    #elif CILK
        cilk_spawn
    #endif
    alt_tree_creation (*tree.right, elemToNode, sepToNode, nodePart, nodePartSize, nodeIndex, sepIndex,
                   globalNbElem, dimElem, separator+1, lastPart, firstElem+nbLeftElem,
                   lastElem-nbSepElem, firstNode+nbLeftNodes, lastNode,
                   firstInnerNode+nbLeftInnerNodes, lastInnerNode-nbSepInnerNodes,
                   sepOffset+nbLeftElem, sepNodeOffset+nbLeftInnerNodes, 3*curNode+2, isSep, depth*10+1);
    #ifdef OMP
        #pragma omp task default(shared)
    #endif
    alt_tree_creation (*tree.left, elemToNode, sepToNode, nodePart, nodePartSize, nodeIndex, sepIndex,
                   globalNbElem, dimElem, firstPart, separator, firstElem, firstElem+
                   nbLeftElem-1, firstNode, lastNode-nbRightNodes,
                   firstInnerNode, firstInnerNode+nbLeftInnerNodes-1, sepOffset, sepNodeOffset,
                   3*curNode+1, isSep, depth*10+2);

    // Synchronization
    #ifdef OMP
        #pragma omp taskwait
    #elif CILK
        cilk_sync;
    #endif

    // D&C partitioning of separator elements
    if (nbSepElem > 0 || nbSepInnerNodes > 0) {
        alt_sep_partitioning (*tree.sep, elemToNode, nodeIndex, globalNbElem, dimElem, lastElem-
                          nbSepElem+1, lastElem, firstNode, lastNode,
                          lastInnerNode-nbSepInnerNodes+1, lastInnerNode, 3*curNode+3, depth);
    }

}

/* -------------------------------------------------------------------------- */
void init_intervals (tree_t &tree,
                        int firstElem, int lastElem,
                        int firstNode, int lastNode,
                        int firstInnerNode, int lastInnerNode,
                        bool isSep) {

    tree.ownedNodes   = nullptr;
    tree.intfIndex    = nullptr;
    tree.intfNodes    = nullptr;
    tree.intfDst      = nullptr;
    tree.nbOwnedNodes = 0;
    tree.nbIntfNodes  = 0;
    tree.firstElem    = firstElem;
    tree.lastElem     = lastElem;
    tree.firstNode    = firstNode;
    tree.lastNode     = lastNode;
    tree.firstInnerNode = firstInnerNode;
    tree.lastInnerNode = lastInnerNode;
    tree.firstEdge    = -1;
    tree.lastEdge     = -1;
    tree.vecOffset    = 0;
    tree.isSep        = isSep;
    tree.left         = nullptr;
    tree.right        = nullptr;
    tree.sep          = nullptr;
}

void alt_tree_creation_new (tree_t &tree, int *elemToNode, int *sepToNode,
                    int *nodePart, int *nodePartSize, int dimElem,
                    int firstPart, int lastPart, int sepOffset, int curNode) {

    printf("Enter in NEW tree_creation ...\n");

    /* -------------------------------------------------------------------------
    INITIALIZATION
    ------------------------------------------------------------------------- */

    int firstElem = tree.firstElem;
    int lastElem = tree.lastElem;
    int firstNode = tree.firstNode;
    int lastNode = tree.lastNode;
    int firstInnerNode = tree.firstInnerNode;
    int lastInnerNode = tree.lastInnerNode;
    int isSep = tree.isSep;

    assert(isSep == false);

    // If current node is a leaf
    int nbPart = lastPart - firstPart + 1;
    int localNbElem = lastElem - firstElem + 1;
    int localNbNodes = lastNode - firstNode + 1;
    if (nbPart < 2 || localNbElem <= MAX_ELEM_PER_PART) {

        #ifdef MULTITHREADED_COMM
            // Set the node accessed by current D&C leaf
            fill_node_owner (elemToNode, firstElem, lastElem, dimElem, firstNode,
                             lastNode, curNode, isSep);
        #endif

        // STOP
        return;
    }

    /* -------------------------------------------------------------------------
    PARTITIONING
    ------------------------------------------------------------------------- */

    // Else, prepare next left, right & separator recursion
    int nbLeftElem = 0, nbSepElem = 0, nbRightElem = 0;
    int nbLeftNodes = 0, nbSepNodes = 0, nbRightNodes = 0;
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
    //int *localNodePart = new int [localNbNodes];
    if (isSep) {
        alt_create_elem_part (localElemPart, nodePart, sepToNode, localNbElem, dimElem,
                          separator, sepOffset, &nbLeftElem, &nbSepElem);
        // select_partition (localElemPart, localNodePart, nodePart, sepToNode,
        //                       dimElem, localNbElem, localNbNodes, separator,
        //                       sepOffset, firstInnerNode,
        //                       &nbLeftElem, &nbSepElem, &nbRightElem,
        //                       &nbLeftNodes, &nbSepNodes, &nbRightNodes);
    }
    else {
        alt_create_elem_part (localElemPart, nodePart, elemToNode, localNbElem, dimElem,
                          separator, firstElem, &nbLeftElem, &nbSepElem);
        // select_partition (localElemPart, localNodePart, nodePart, elemToNode,
        //                       dimElem, localNbElem, localNbNodes, separator,
        //                       firstElem, firstInnerNode,
        //                       &nbLeftElem, &nbSepElem, &nbRightElem,
        //                       &nbLeftNodes, &nbSepNodes, &nbRightNodes);
    }

    printf("parts: %d + %d + %d // %d + %d + %d\n",nbLeftElem, nbSepElem, nbRightElem, nbLeftNodes, nbSepNodes, nbRightNodes);

    /* -------------------------------------------------------------------------
    PERMUTATIONS
    ------------------------------------------------------------------------- */

    // Create local element permutation
    int *localElemPerm = new int [localNbElem];
    //int *localNodePerm = new int [localNbNodes];
    DC_create_permutation (localElemPerm, localElemPart, localNbElem, 3);
    //DC_create_permutation (localNodePerm, localNodePart, localNbNodes, 3);
    delete[] localElemPart;
    //delete[] localNodePart;

    // Execution is correct without mutex although cilkscreen detects a race condition
    pthread_mutex_lock (&mergeMutex);
    // Apply local element permutation to global element permutation
    alt_merge_permutations (elemPerm, localElemPerm, globalNbElem, localNbElem, firstElem, lastElem);
    //alt_merge_permutations (nodePerm, localNodePerm, globalNbNodes, localNbNodes, firstNode, lastNode);
    pthread_mutex_unlock (&mergeMutex);

    // Permute elemToNode and sepToNode with local element permutation
    DC_permute_int_2d_array (elemToNode, localElemPerm, localNbElem, dimElem, firstElem);

    if (isSep) {
        DC_permute_int_2d_array (sepToNode, localElemPerm, localNbElem, dimElem,
                                 sepOffset);
    }
    delete[] localElemPerm;
    //delete[] localNodePerm;

    /* -------------------------------------------------------------------------
    RECURSIVE CALLS
    ------------------------------------------------------------------------- */

    tree.left  = new tree_t;
    init_intervals(*tree.left, firstElem, firstElem+nbLeftElem-1, firstNode, lastNode-nbRightNodes, 0, 0, isSep);
    tree.right = new tree_t;
    init_intervals(*tree.right, firstElem+nbLeftElem, lastElem-nbSepElem, firstNode+nbLeftNodes, lastNode, 0, 0, isSep);

    // Left & right recursion
    #ifdef OMP
        #pragma omp task default(shared)
    #elif CILK
        cilk_spawn
    #endif
    alt_tree_creation_new (*tree.right, elemToNode, sepToNode, nodePart,
                        nodePartSize, dimElem, separator+1, lastPart,
                        sepOffset+nbLeftElem, 3*curNode+2);
    #ifdef OMP
        #pragma omp task default(shared)
    #endif
    alt_tree_creation_new (*tree.left, elemToNode, sepToNode, nodePart,
                        nodePartSize, dimElem, firstPart, separator,
                        sepOffset, 3*curNode+1);

    // Synchronization
    #ifdef OMP
        #pragma omp taskwait
    #elif CILK
        cilk_sync;
    #endif


}

#endif // TREE_CREATION
