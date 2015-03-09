/*  Copyright 2014 - UVSQ
    Authors list: Loïc Thébault, Eric Petit

    This file is part of the D&C library.

    D&C library is free software: you can redistribute it and/or modify it under the
    terms of the GNU Lesser General Public License as published by the Free Software
    Foundation, either version 3 of the License, or (at your option) any later version.

    D&C library is distributed in the hope that it will be useful, but WITHOUT ANY
    WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
    PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License along with
    the D&C library. If not, see <http://www.gnu.org/licenses/>. */

#ifdef TREE_CREATION

#include <string.h>
#include <math.h>
#include <cilk/cilk.h>
#include <pthread.h>
#include <metis.h>

#include "tools.h"
#include "permutations.h"
#include "tree_creation.h"
#include "partitioning.h"

extern tree_t *treeHead;
extern int *elemPerm, *nodePerm;

pthread_mutex_t metisMutex = PTHREAD_MUTEX_INITIALIZER;

// Create a nodal graph from a tetrahedron mesh (created from METIS)
void mesh_to_nodal (int *graphIndex, int *graphValue, int *elemToNode,
                    int nbElem, int dimElem, int nbNodes)
{
    int nEdges, *nPtr, *nInd, *marker;

    nPtr = new int [nbNodes + 1] ();
    for (int i = 0; i < dimElem * nbElem; i++) {
        nPtr[elemToNode[i]]++;
    }
    for (int i = 1; i < nbNodes; i++) {
        nPtr[i] += nPtr[i-1];
    }
    for (int i = nbNodes; i > 0; i--) {
        nPtr[i] = nPtr[i-1];
    }
    nPtr[0] = 0;

    nInd = new int [nPtr[nbNodes]];
    for (int i = 0; i < nbElem; i++) {
        for (int j = 0; j < dimElem; j++) {
            nInd[nPtr[elemToNode[i*dimElem+j]]++] = i;
        }
    }
    for (int i = nbNodes; i > 0; i--) {
        nPtr[i] = nPtr[i-1];
    }
    nPtr[0] = 0;

    marker = new int [nbNodes];
    memset (marker, -1, nbNodes * sizeof (int));
    nEdges = graphIndex[0] = 0;
    for (int i = 0; i < nbNodes; i++) {
        marker[i] = i;
        for (int j = nPtr[i]; j < nPtr[i+1]; j++) {
            int jj = dimElem * nInd[j];
            for (int k = 0; k < dimElem; k++) {
                int kk = elemToNode[jj];
                if (marker[kk] != i) {
                    marker[kk] = i;
                    graphValue[nEdges++] = kk;
                }
                jj++;
            }
        }
        graphIndex[i+1] = nEdges;
    }
    delete[] marker, delete[] nInd, delete[] nPtr;
}

// Create temporal sepToNode array containing separator elements indexed
// contiguously from 0 to nbSepElem and return the number of separator nodes
int create_sepToNode (int *sepToNode, int *elemToNode, int firstSepElem,
                      int lastSepElem, int dimElem)
{
    int nbSepNodes = 0;
    int *sepRenum = new int [(lastSepElem - firstSepElem + 1) * dimElem];

	for (int i = firstSepElem * dimElem, j = 0; i < (lastSepElem+1)*dimElem;
             i++, j++) {
		int newNode, oldNode = elemToNode[i];
		bool isNew = true;
		for (newNode = 0; newNode < nbSepNodes; newNode++) {
			if (oldNode == sepRenum[newNode]) {
				isNew = false;
				break;
			}
		}
		if (isNew) {
            sepRenum[nbSepNodes] = oldNode;
            nbSepNodes++;
        }
		sepToNode[j] = newNode;
	}
    delete[] sepRenum;
	return nbSepNodes;
}

// D&C partitioning of separators with more than MAX_ELEM_PER_PART elements
void sep_partitioning (tree_t &tree, int *elemToNode, int globalNbElem, int dimElem,
                       int firstSepElem, int lastSepElem
#ifdef STATS
					   , ofstream &dcFile, int curNode)
#else
                       )
#endif
{
    // If there is not enough element in the separator, fill the D&C tree & exit
    int nbSepElem = lastSepElem - firstSepElem + 1;
    int nbSepPart = ceil (nbSepElem / (double)MAX_ELEM_PER_PART);
    if (nbSepPart < 2 || nbSepElem <= MAX_ELEM_PER_PART) {
        init_dc_tree (tree, firstSepElem, lastSepElem, 0, -1, -1, true);
#ifdef STATS
        fill_dc_file_leaves (dcFile, curNode, firstSepElem, lastSepElem, 3);
#endif
        return;
    }

    // Create temporal elemToNode containing the separator elements
    int *sepToNode = new int [nbSepElem * dimElem];
    int nbSepNodes = create_sepToNode (sepToNode, elemToNode, firstSepElem,
                                       lastSepElem, dimElem);

    // Configure METIS & compute the node partitioning of the separators
    int constraint = 1, objVal;
    int *graphIndex = new int [nbSepNodes + 1];
    int *graphValue = new int [nbSepNodes * 15];
    int *nodePart   = new int [nbSepNodes];
    mesh_to_nodal (graphIndex, graphValue, sepToNode, nbSepElem, dimElem, nbSepNodes);
    // Execution is correct without mutex although cilkscreen detects many race
    // conditions. Check if the problem is solved with future version of METIS (5.0)
    pthread_mutex_lock (&metisMutex);
    METIS_PartGraphRecursive (&nbSepNodes, &constraint, graphIndex, graphValue,
                              nullptr, nullptr, nullptr, &nbSepPart, nullptr, nullptr,
                              nullptr, &objVal, nodePart);
    pthread_mutex_unlock (&metisMutex);
    delete[] graphValue, delete[] graphIndex;

    // Create the separator D&C tree
    recursive_tree_creation (tree, elemToNode, sepToNode, nodePart, nullptr,
                             globalNbElem, dimElem, 0, nbSepPart-1, firstSepElem,
                             lastSepElem, -1, -1, 0
    #ifdef STATS
                             , dcFile, curNode, -1);
    #else
                             );
    #endif
    delete[] nodePart, delete[] sepToNode;
}

// Divide & Conquer partitioning of elemToNode array
void partitioning (int *elemToNode, int nbElem, int dimElem, int nbNodes)
{
	// Fortran to C elemToNode conversion
    #ifdef OMP
        #pragma omp parallel for
        for (int i = 0; i < nbElem * dimElem; i++) {
    #else    
        cilk_for (int i = 0; i < nbElem * dimElem; i++) {
    #endif    
    elemToNode[i]--;
    }

    // Configure METIS & compute the node partitioning of the mesh
    int nbPart = ceil (nbElem / (double)MAX_ELEM_PER_PART);
	int constraint = 1, objVal;
    int *graphIndex = new int [nbNodes + 1];
    int *graphValue = new int [nbNodes * 15];
	int *nodePart   = new int [nbNodes];
    mesh_to_nodal (graphIndex, graphValue, elemToNode, nbElem, dimElem, nbNodes);
    METIS_PartGraphRecursive (&nbNodes, &constraint, graphIndex, graphValue,
                              nullptr, nullptr, nullptr, &nbPart, nullptr, nullptr,
                              nullptr, &objVal, nodePart);
    delete[] graphValue, delete[] graphIndex;

	// Initialize the global element permutation
    #ifdef OMP
        #pragma omp parallel for
    	for (int i = 0; i < nbElem; i++) {
    #else
    	cilk_for (int i = 0; i < nbElem; i++) {
    #endif
		elemPerm[i] = i;
	}
    // Compute the number of nodes per partition
    int *nodePartSize = new int [nbPart] ();
    for (int i = 0; i < nbNodes; i++) {
        nodePartSize[nodePart[i]]++;
    }

	// Create D&C tree & its dot file
    #ifdef STATS
        string fileName = "dcTree_" +
                          to_string ((long long)MAX_ELEM_PER_PART) + ".dot";
    	ofstream dcFile (fileName, ios::out | ios::trunc);
    	init_dc_file (dcFile, nbPart);
    	recursive_tree_creation (*treeHead, elemToNode, nullptr, nodePart,
                                 nodePartSize, nbElem, dimElem, 0, nbPart-1, 0,
                                 nbElem-1, 0, nbNodes-1, 0, dcFile, 0, -1);
    	close_dc_file (dcFile);
    	dc_stat ();
    #else
		#ifdef OMP
		{
			#pragma omp parallel
			{
				#pragma omp single nowait
				{
					recursive_tree_creation (*treeHead, elemToNode, nullptr, nodePart,
											nodePartSize, nbElem, dimElem, 0, nbPart-1, 0,
											nbElem-1, 0, nbNodes-1, 0);
				}
			}
		}
		#else
			recursive_tree_creation (*treeHead, elemToNode, nullptr, nodePart,
										nodePartSize, nbElem, dimElem, 0, nbPart-1, 0,
										nbElem-1, 0, nbNodes-1, 0);
		#endif	
    #endif
    delete[] nodePartSize;

	// Create node permutation from node partition
	DC_create_permutation (nodePerm, nodePart, nbNodes, nbPart);
	delete[] nodePart;

	// C to Fortran elemToNode conversion
    #ifdef OMP
        #pragma omp parallel for
	    for (int i = 0; i < nbElem * dimElem; i++) {
    #else
    	cilk_for (int i = 0; i < nbElem * dimElem; i++) {
    #endif	
	elemToNode[i]++;
	}
}

#endif
