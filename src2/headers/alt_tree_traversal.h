#ifndef ALT_TREE_TRAVERSAL_H
#define ALT_TREE_TRAVERSAL_H

#include "DC.h"

// Follow the D&C tree to execute the given function in parallel
void alt_tree_traversal (void (*userSeqFct)  (void *, DCargs_t *, DCreturnArgs_t *),
                     void (*userVecFct)  (void *, DCargs_t *),
                     void (*userCommFct) (void *, DCcommArgs_t *),
                     void *userArgs, void *userCommArgs, DCreturnArgs_t *returnArgs, tree_t &tree);


// Wrapper used to get the root of the D&C tree before calling the real tree traversal
void DC_alt_tree_traversal (void (*userSeqFct)  (void *, DCargs_t *, DCreturnArgs_t *),
                        void (*userVecFct)  (void *, DCargs_t *),
                        void (*userCommFct) (void *, DCcommArgs_t *),
                        void *userArgs, void *userCommArgs, DCreturnArgs_t *returnArgs);

#endif // ALT_TREE_TRAVERSAL_H
