#ifndef ALT_PERMUTATIONS_H
#define ALT_PERMUTATIONS_H

#include "DC.h"

void alt_merge_permutations (int *globalPerm, int *localPerm, int globalNb, int localNb, int first, int last);

void DC_alt_permute_double_2d_array (double *tab, int *perm, int nbItem, int dimItem, int offset);

#endif // ALT_PERMUTATIONS_H
