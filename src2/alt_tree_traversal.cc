#ifdef CILK
    #include <cilk/cilk.h>
#endif
#include <stdlib.h>

#include "alt_tree_traversal.h"

extern tree_t *treeHead;

int DC_get_max_elem_per_part() {
    return MAX_ELEM_PER_PART;
}

// Follow the D&C tree to execute the given function in parallel
void alt_tree_traversal (void (*userSeqFct)  (void *, DCargs_t *, DCreturnArgs_t *),
                     void (*userVecFct)  (void *, DCargs_t *),
                     void (*userCommFct) (void *, DCcommArgs_t *),
                     void *userArgs, void *userCommArgs, DCreturnArgs_t *returnArgs, tree_t &tree)
{
    // If current node is a leaf, call the appropriate function
    if (tree.left == nullptr && tree.right == nullptr) {

        // Initialize the D&C arguments
        DCargs_t DCargs;
        DCargs.firstNode = tree.firstNode;
        DCargs.lastNode  = tree.lastNode;
        DCargs.firstEdge = tree.firstEdge;
        DCargs.firstInnerNode = tree.firstInnerNode;
        DCargs.lastInnerNode = tree.lastInnerNode;
        DCargs.lastEdge  = tree.lastEdge;
        DCargs.isSep     = tree.isSep;
        DCargs.isLeaf     = tree.isLeaf;
        DCargs.depth     = tree.depth;
        #ifdef MULTITHREADED_COMM
            DCargs.nbOwnedNodes = tree.nbOwnedNodes;
            DCargs.ownedNodes   = tree.ownedNodes;
        #endif

        #ifdef DC_VEC
            // Call user vectorial function on full colors
            DCargs.firstElem = tree.firstElem;
            DCargs.lastElem  = tree.vecOffset;
            userVecFct (userArgs, &DCargs);

            // Call user sequential function on other colors
            DCargs.firstElem = tree.vecOffset+1;
            DCargs.lastElem  = tree.lastElem;
            userSeqFct (userArgs, &DCargs, nullptr);
        #else
            // Call user sequential function
            DCargs.firstElem = tree.firstElem;
            DCargs.lastElem  = tree.lastElem;
            userSeqFct (userArgs, &DCargs, returnArgs);
        #endif

        // Call user communication function
        #ifdef MULTITHREADED_COMM
            if (tree.nbIntfNodes > 0) {
                DCcommArgs_t DCcommArgs;
                DCcommArgs.intfIndex = tree.intfIndex;
                DCcommArgs.intfNodes = tree.intfNodes;
                DCcommArgs.intfDst   = tree.intfDst;
                userCommFct (userCommArgs, &DCcommArgs);
            }
        #endif
    }
    else {
        DCreturnArgs_t *returnLeft = nullptr;
        DCreturnArgs_t *returnRight = nullptr;
        DCreturnArgs_t *returnSep = nullptr;
        if (returnArgs != nullptr) {
            returnLeft = new DCreturnArgs_t(returnArgs->size, returnArgs->op);
            returnRight = new DCreturnArgs_t(returnArgs->size, returnArgs->op);
            returnSep = new DCreturnArgs_t(returnArgs->size, returnArgs->op);
        }
        #ifdef OMP
            // Left & right recursion
            #pragma omp task default(shared)
            alt_tree_traversal (userSeqFct, userVecFct, userCommFct, userArgs,
                            userCommArgs, returnRight, *tree.right);
            #pragma omp task default(shared)
            alt_tree_traversal (userSeqFct, userVecFct, userCommFct, userArgs,
                            userCommArgs, returnLeft, *tree.left);
            // Synchronization
            #pragma omp taskwait
        #elif CILK
            // Left & right recursion
            cilk_spawn
            alt_tree_traversal (userSeqFct, userVecFct, userCommFct, userArgs,
                            userCommArgs, returnRight, *tree.right);
            alt_tree_traversal (userSeqFct, userVecFct, userCommFct, userArgs,
                            userCommArgs, returnLeft, *tree.left);
            // Synchronization
            cilk_sync;
        #endif

        // Separator recursion, if it is not empty
        if (tree.sep != nullptr) {
            alt_tree_traversal (userSeqFct, userVecFct, userCommFct, userArgs,
                            userCommArgs, returnSep, *tree.sep);
        }

        // Reduction
        if (returnArgs != nullptr) {
            for (int i = 0; i < returnArgs->size; i++) {
                returnArgs->val[i] = returnRight->val[i] + returnLeft->val[i];
                if (tree.sep != nullptr) {
                    returnArgs->val[i] += returnSep->val[i];
                }
            }
        }

    }
}

// Wrapper used to get the root of the D&C tree before calling the real tree traversal
void DC_alt_tree_traversal (void (*userSeqFct)  (void *, DCargs_t *, DCreturnArgs_t *),
                        void (*userVecFct)  (void *, DCargs_t *),
                        void (*userCommFct) (void *, DCcommArgs_t *),
                        void *userArgs, void *userCommArgs, DCreturnArgs_t *returnArgs)
{
    #ifdef OMP
        #pragma omp parallel
        #pragma omp single nowait
    #endif
    alt_tree_traversal (userSeqFct, userVecFct, userCommFct, userArgs, userCommArgs, returnArgs,
                    *treeHead);

}
