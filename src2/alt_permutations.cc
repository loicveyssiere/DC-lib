extern int *elemPerm, *nodePerm, *innerNodePerm;

void alt_merge_permutations (int *globalPerm, int *localPerm, int globalNb, int localNb, int first, int last) {

    int ctr = 0;
    for (int i = 0; i < globalNb; i++) {
        int dst = globalPerm[i];
        if (dst >= first && dst <= last) {
            globalPerm[i] = localPerm[dst-first] + first;
            ctr++;
        }
        if (ctr == localNb)	break;
    }
}


void DC_alt_permute_double_2d_array (double *tab, int *perm, int nbItem, int dimItem, int offset)
{
    char *checkPerm = new char [nbItem] ();
    double  *tmpSrc    = new double  [dimItem];
    double  *tmpDst    = new double  [dimItem];

    // If no permutation is given, default behavior is to use D&C elemPerm
    if (perm == nullptr) perm = elemPerm;

    for (int i = 0; i < nbItem; i++) {
        if (checkPerm[i] == 1) continue;

        int init = i, src = i, dst;
        for (int j = 0; j < dimItem; j++) {
            tmpSrc[j] = tab[(i+offset)*dimItem+j];
        }
        do {
            dst = perm[src];
            for (int j = 0; j < dimItem; j++) {
                tmpDst[j] = tab[(dst+offset)*dimItem+j];
                tab[(dst+offset)*dimItem+j] = tmpSrc[j];
                tmpSrc[j] = tmpDst[j];
            }
            src = dst;
            checkPerm[src] = 1;
        }
        while (src != init);
    }
    delete[] tmpDst, delete[] tmpSrc, delete[] checkPerm;
}
