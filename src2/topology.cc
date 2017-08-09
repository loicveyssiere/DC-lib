#include "DC.h"

using namespace std;

Tree_topo::Tree_topo(string stats_file, int n_extra) {

    this->stats = stats_file;
    this->n_extra = n_extra;
    stream.open(stats);

    stream << "isSep,isLeaf,depth,nbElem,nbInner,nbNode";
    for (int i = 0; i < this->n_extra; i++) {
        stream << ",t" << (i+1);
    }
    stream << endl;

}

Tree_topo::~Tree_topo() {
    stream.close();
}

void Tree_topo::add(DCargs_t *DCargs, int *extra) {
    int firstInnerNode = DCargs->firstInnerNode;
    int lastInnerNode = DCargs->lastInnerNode;
    int diff = DCargs->lastElem - DCargs->firstElem +1;
    int innerDiff = lastInnerNode - firstInnerNode+1;
    int nodeDiff = DCargs->lastNode - DCargs->firstNode+1;
    int isSep = DCargs->isSep;
    int isLeaf = DCargs->isLeaf;
    int depth = DCargs->depth;

    stream << isSep << "," << isLeaf << "," << depth << ","
           << diff << "," << innerDiff << "," << nodeDiff;
   for (int i = 0; i < n_extra; i++) {
       stream << "," << extra[i];
   }
   stream << endl;
}
