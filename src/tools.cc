/*  Copyright 2013 - UVSQ
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

#ifdef XMPI
    #include <mpi.h>
#elif GASPI
    #include <GASPI.h>
#endif
#ifdef CILK
    #include <cilk/cilk.h>
#endif
#include <iostream>
#include <stdio.h>
#include <sys/time.h>

#include "tools.h"

extern tree_t *treeHead;

/*****************************************************************************/
/***********                        Timer                          ***********/
/*****************************************************************************/

// Get time of day
double DC_get_time ()
{
    struct timeval tv;
    gettimeofday (&tv, NULL);
    return (double)tv.tv_sec + (double) tv.tv_usec * 1.e-6;
}

// RDTSC
uint64_t DC_get_cycles ()
{
    uint64_t a, d;
    __asm__ volatile ("rdtsc" : "=a" (a), "=d" (d));
    return (d << 32) | a;
}

// Constructor
DC_timer::DC_timer () : avgTime (0), timeCtr (0), avgCycles (0), cyclesCtr (0)
{}

// Return the average time of day
double DC_timer::get_avg_time ()
{
    return avgTime;
}

// Reset the average time of day
void DC_timer::reset_time ()
{
    avgTime = 0;
    timeCtr = 0;
}

// Stop time of day timer
void DC_timer::stop_time ()
{
    double stopTime = DC_get_time ();
    avgTime = timeCtr * avgTime + (stopTime - startTime);
    timeCtr++;
    avgTime /= timeCtr;
}

// Start time of day timer
void DC_timer::start_time ()
{
    startTime = DC_get_time ();
}

// Return the average cycle counter
uint64_t DC_timer::get_avg_cycles ()
{
    return avgCycles;
}

// Reset the average cycles counter
void DC_timer::reset_cycles ()
{
    avgCycles = 0;
    cyclesCtr = 0;
}

// Stop cycle counter
void DC_timer::stop_cycles ()
{
    uint64_t stopCycles = DC_get_cycles ();
    avgCycles = cyclesCtr * avgCycles + (stopCycles - startCycles);
    cyclesCtr++;
    avgCycles /= cyclesCtr;
}

// Start cycle counter
void DC_timer::start_cycles ()
{
    startCycles = DC_get_cycles ();
}

/*****************************************************************************/
/***********                      Quick sort                       ***********/
/*****************************************************************************/

// Sort by ascending node values couple arrays using parallel quick sort
void quick_sort (couple_t *tab, int begin, int end)
{
    int left = begin - 1, right = end + 1, pivot = tab[begin].node;

    if (begin < end) {
        while (1)
        {
            do right--; while(tab[right].node > pivot);
            do left++;  while(tab[left].node  < pivot);

            if (left < right) {
                couple_t tmp = tab[left];
                tab[left]    = tab[right];
                tab[right]   = tmp;
            }
            else break;
        }
        #ifdef OMP
            #pragma omp task default (shared)
            quick_sort (tab, begin, right);
            #pragma omp task default (shared)
            quick_sort (tab, right+1, end);
            #pragma omp taskwait
        #elif CILK
            cilk_spawn
            quick_sort (tab, begin, right);
            quick_sort (tab, right+1, end);
            cilk_sync;
        #endif
    }
}

/*****************************************************************************/
/***********                  Vectorization stats                  ***********/
/*****************************************************************************/

// Compute the number of colors used per leaf and the number of elements per color
void leaf_coloring_stat (ofstream &colorPerLeaf, int *elemPerColor, int *colorPart,
                         int nbElem, int nbColors)
{
	static int curLeaf = 0;
	colorPerLeaf << curLeaf << " " << nbColors << endl;
	curLeaf++;

    int *colorSize = new int [nbColors] ();
    for (int i = 0; i < nbElem; i++) {
        colorSize[colorPart[i]]++;
		if (colorSize[colorPart[i]] > MAX_ELEM_PER_PART) {
			cerr << "Error: too much elements per color.\n";
			exit (EXIT_FAILURE);
		}
	}
	for (int i = 0; i < nbColors; i++) {
		elemPerColor[colorSize[i]-1] += colorSize[i];
	}
    delete[] colorSize;
}

// Fill the elemToElem & the elemPerColor files
void coloring_stat (int *elemPerColor, int nbElem)
{
	string fileName = "elemPerColor.csv";
	ofstream elemColor (fileName, ios::out | ios::trunc);
	elemColor << "colorSize nbElem\n";
	for (int i = 0; i < MAX_ELEM_PER_PART; i++) {
		elemColor << i+1 << " " << elemPerColor[i] << endl;
	}
	elemColor.close ();
}

// Compute the number of elements per leaf of the D&C tree
void leaf_dc_stat (tree_t &tree, ofstream &elemPerLeaf)
{
	if (tree.left == nullptr && tree.right == nullptr) {
		static int curLeaf = 0;
		elemPerLeaf << curLeaf << " " << tree.lastSep-tree.firstElem+1 << endl;
		curLeaf++;
	}
	else {
		leaf_dc_stat (*tree.left, elemPerLeaf);
		leaf_dc_stat (*tree.right, elemPerLeaf);
		if (tree.sep != nullptr) {
			leaf_dc_stat (*tree.sep, elemPerLeaf);
		}
	}
}

// Fill the elemPerLeaf file
void dc_stat ()
{
	string fileName = "elemPerLeaf.csv";
	ofstream elemPerLeaf (fileName, ios::out | ios::trunc);
	elemPerLeaf << "leafNb nbElem\n";
	leaf_dc_stat (*treeHead, elemPerLeaf);
	elemPerLeaf.close ();
}

/*****************************************************************************/
/***********                   Interfaces stats                    ***********/
/*****************************************************************************/

// Mutex to avoid race condition when incrementing leaf counters
pthread_mutex_t intfStatsMutex = PTHREAD_MUTEX_INITIALIZER;
int nbLeaves = 0, nbIntfLeaves = 0;

// Count the total number of leaves and the number of interface leaves
void count_intf_stats (bool hasIntfNode)
{
    pthread_mutex_lock (&intfStatsMutex);
    nbLeaves++;
    if (hasIntfNode) nbIntfLeaves++;
    pthread_mutex_unlock (&intfStatsMutex);
}

// Store statistics on the number of interface leaves
void store_intf_stats (int nbElem, int rank)
{
    string fileName = "intfRatio";
    ofstream intfFile (fileName, ios::out | ios::trunc);
    int localStats[2], globalStats[2];
    localStats[0] = nbElem;
    localStats[1] = nbIntfLeaves * 100 / nbLeaves;

    #ifdef XMPI
        MPI_Reduce (localStats, globalStats, 2, MPI_INT, MPI_MAX, 0,
                    MPI_COMM_WORLD);
    #elif GASPI
        gaspi_allreduce (localStats, globalStats, 2, GASPI_OP_MAX,
                         GASPI_TYPE_INT, GASPI_GROUP_ALL, GASPI_BLOCK);
    #endif

    if (rank == 0) {
        intfFile << "nbElements intfRatio\n" << globalStats[0] << " " << globalStats[1]
                 << endl;
    }
}

/*****************************************************************************/
/***********                   D&C tree dot file                   ***********/
/*****************************************************************************/

// Fill the leaves of the D&C tree dot file
void fill_dc_file_leaves (tree_t &tree, ofstream &dcFile, int curNode, int curLevel,
                          int LRS, bool hasIntfNode)
{
	dcFile << "\t" << curNode << " [label=\"" << curNode << "\\n["
		   << tree.firstElem << "," << tree.lastElem << "]\\n["
		   << tree.firstNode << "," << tree.lastNode
           << "]\\n[" << tree.firstInnerNode << "," << tree.lastInnerNode
           << "]\\n[" << tree.lastElem-tree.firstElem+1 << " - " << tree.lastInnerNode-tree.firstInnerNode+1
           << "]\"";

	if      (LRS == 1) dcFile << ", color=turquoise4];\n";
	else if (LRS == 2) dcFile << ", color=lightskyblue];\n";
	else if (LRS == 3) dcFile << ", color=grey];\n";
	else			   dcFile << ", shape=circle, color=red];\n";
}

// Fill the nodes of the D&C tree dot file
void fill_dc_file_nodes (tree_t &tree, ofstream &dcFile, int curNode, int curLevel,
                         bool hasIntfNode)
{
	dcFile << "\t" << curNode << " -> {" << 3*curNode+1 << "; " << 3*curNode+2;
	if (tree.sep != nullptr) dcFile << "; " << 3*curNode+3 << ";}\n\t";
	else			         dcFile << ";}\n\t";
	dcFile << curNode << " [label=\"" << curNode << "\\n[" << tree.firstElem
		   << "," << tree.lastElem << "," << tree.lastSep << "]\\n["
           << tree.firstNode << "," << tree.lastNode << "]\", style=rounded";
    dcFile << "];\n";
}

// Detect if given D&C node has nodes on the interface
bool has_intf_node (tree_t &tree, int *intfIndex, int *intfNodes, int nbIntf)
{
    for (int i = 0; i < nbIntf; i++) {
        for (int j = intfIndex[i]; j < intfIndex[i+1]; j++) {
            int intfNode = intfNodes[j] - 1;
            for (int k = tree.firstNode; k <= tree.lastNode; k++) {
                if (intfNode == k) {
                    return true;
                }
            }
        }
    }
    return false;
}

// Close the D&C tree dot file
void close_dc_file (ofstream &dcFile)
{
	dcFile << "}\n";
	dcFile.close ();
}

// Initialize the D&C tree dot file with default layout
void init_dc_file (ofstream &dcFile)
{
	if (!dcFile) cerr << "Error opening dcFile.\n";

	dcFile << "digraph RecursiveTree {\n\t-1 [label=\"Max elements per partition : "
           << MAX_ELEM_PER_PART << "\", shape=plaintext];\n"
		   << "\tnode [shape=box, style=\"rounded,filled\"];\n";
}

/*****************************************************************************/
/***********                       2D matrix                       ***********/
/*****************************************************************************/

// Determine for each tree node if it's a left, right or separator node
void recursive_2d_matrix (tree_t &tree, int *checkNode, int *elemToNode, int dimElem,
                          int LRS)
{
	if (tree.left == nullptr && tree.right == nullptr) {
		for (int i = tree.firstElem; i <= tree.lastSep; i++) {
			for (int j = 0; j < dimElem; j++) {
				int node = elemToNode[i*dimElem+j] - 1;
				checkNode[node] = LRS;
			}
		}
	}
	else {
		if (LRS == 3) {
			recursive_2d_matrix (*tree.left,  checkNode, elemToNode, dimElem, 3);
			recursive_2d_matrix (*tree.right, checkNode, elemToNode, dimElem, 3);
			if (tree.sep != nullptr) {
				recursive_2d_matrix (*tree.sep,   checkNode, elemToNode, dimElem, 3);
			}
		}
		else {
			recursive_2d_matrix (*tree.left,  checkNode, elemToNode, dimElem, 1);
			recursive_2d_matrix (*tree.right, checkNode, elemToNode, dimElem, 2);
			if (tree.sep != nullptr) {
				recursive_2d_matrix (*tree.sep,   checkNode, elemToNode, dimElem, 3);
			}
		}
	}
}

// Create a 2D representation of a given CSR matrix
void DC_create_2d_matrix (int *elemToNode, int *nodeToNodeRow, int *nodeToNodeColumn,
                          int nbNodes, int dimElem)
{
	int *checkNode = new int [nbNodes] ();
    string fileName = "2Dmatrix.csv";
	recursive_2d_matrix (*treeHead, checkNode, elemToNode, dimElem, -1);
	ofstream map (fileName, ios::out | ios::trunc);
	map << "X,Y,color\n";

	for (int i = 0; i < nbNodes; i++) {
		for (int j = nodeToNodeRow[i]; j < nodeToNodeColumn[i+1]; j++) {
			map << i << "," << nodeToNodeColumn[j]-1;
			if     (checkNode[i] == 1 || checkNode[nodeToNodeColumn[j]-1] == 1)
				map << ",1\n";
			else if(checkNode[i] == 2 || checkNode[nodeToNodeColumn[j]-1] == 2)
				map << ",2\n";
			else if(checkNode[i] == 3 && checkNode[nodeToNodeColumn[j]-1] == 3)
				map << ",3\n";
		}
	}

	delete[] checkNode;
}
