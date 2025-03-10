/**
 * @file lookuptable.cpp
 * @author Adrian Del Maestro.
 *
 * @brief LookupTable class implementation.
 */

#include "lookuptable.h"
#include "communicator.h"
#include "path.h"

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// LOOKUP TABLE CLASS --------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/**************************************************************************//**
 * Initilialize the nearest neighbor lookup table.
 *
 * We partition the simulation cell into a grid of boxes whose edge length is
 * defined by the global potential cutoff radius rc.  All date structures 
 * are initialized to be empty.  We must call an updateBeads method externally
 * in order to place the worldlines in the lookup table.
******************************************************************************/
LookupTable::LookupTable(const Container *_boxPtr, const int _numLookupTimeSlices, 
        const int _numParticles) : 
    boxPtr(_boxPtr),
    numLookupTimeSlices(_numLookupTimeSlices)
{
    /* Setup the nearest neighbor grid list */
    setupNNGrid();

    /* We begin with the assumption of evenly distributed particles, therefore
     * there could be numLabels in each cell, allowing for some breathing room */
    for (int i = 0; i < NDIM; i++) 
        hashSize[i] = numNNGrid[i];
    hashSize[NDIM] = numLookupTimeSlices;
    hashSize[NDIM+1] = static_cast<int>(_numParticles/totNumGridBoxes) + 5;

    /* Resize and initialize the main hash array */
    hash.resize(hashSize);
    hash.fill(XXX);

    /* Resize and initialize the grid and bead label and list arrays */
    resizeList(_numParticles);
    grid.fill(make_array<iVec>(XXX));
    beadLabel.fill(XXX);
    beadList.fill({XXX, XXX});
    fullBeadList.fill({XXX, XXX});
    beadSep.fill(dVec{});

    /* Initialize the cutoff^2 */
    rc2 = constants()->rc2();
}

/**************************************************************************//**
 * Free all memory.
******************************************************************************/
LookupTable::~LookupTable() {
}

/**************************************************************************//**
 *  Setup the nearest neighbor grid and lookup tables that will be used
 *  to speed the evaluation of the potential.
******************************************************************************/
void LookupTable::setupNNGrid() {

    /* We determine the number of grid boxes in each dimension and compute
     * their size */
    totNumGridBoxes = 1;
    for (int i = 0; i < NDIM; i++) {
        numNNGrid[i] = static_cast<int>(floor((boxPtr->side[i] / constants()->rc()) + EPS));

        /* Make sure we have at least one grid box */
        if (numNNGrid[i] < 1)
            numNNGrid[i] = 1;

        /* Compute the actual size of the grid */
        sizeNNGrid[i] = boxPtr->side[i] / (1.0 * numNNGrid[i]);

        /* Determine the total number of grid boxes */
        totNumGridBoxes *= numNNGrid[i];
    }

    /* Resize and initialize the numLabels array*/
    std::array <int,NDIM+1> initNumLabels;
    for (int i = 0; i < NDIM; i++) 
        initNumLabels[i] = numNNGrid[i];
    initNumLabels[NDIM] = constants()->numTimeSlices();
    numLabels.resize(initNumLabels);
    numLabels.fill(0);

    /* Now we set the array which holds the nearest neighbors of each
     * grid box which will be used in calculating the potential.  We need
     * to take into account that we may have periodic boundary conditions */
    numNN = ipow(3,NDIM);                           // The total number of NN
    numUniqueNN = (int) floor(0.5*(numNN-1) + EPS); // The unique NN
    std::array <int,NDIM+1> init;

    /* Get the init vector used to resize data structures */
    for (int i = 0; i < NDIM; i++) 
        init[i] = numNNGrid[i];
    init[NDIM] = numNN;

    /* The truncated list of NN indices with no double back references and the 
     * full list of NN indices */
    gridNNReduced.resize(init);
    gridNN.resize(init);

    /* This is somewhat complicated.  Basically we want to construct the list of 
     * nearest neighbors of a given grid box for general dimension.  This consists
     * of moving 'forward', 'zero' and 'back' in each dimension */
    DynamicArray <iVec,1> nnShift(numNN);
    nnShift.fill(iVec{}); 
    /* The shift vector */
    std::array<int,3> shift{-1,0,1};

    /* For each of the unique nearest neighbors, we construct shift vectors 
     * this includes all redundancies and even a zero shift, this is taken
     * into account later */
    int m = 0;
    for (int n = 0; n < numNN; n++) {
        for (int i = 0; i < NDIM; i++) {
            int index = ((int) floor(1.0*n/pow(3.0,1.0*(NDIM-1-i)) + EPS)) % 3;
            nnShift(m)[i] = shift[index];
        }
        m++;
    }

    /* We loop over each grid box, constructing the nearest neighbor list */
    for (int n = 0; n < totNumGridBoxes; n++) {

        /* Get the NDIM-vec coordiantes of the box n */
        gIndex = gridIndex(n);

        /* Copy over into the nn index */
        for (int i = 0; i < NDIM; i++) 
            nnIndex[i] = gIndex[i];
            
        /* Now we go through the nnShift array, and find the nearest neighbor
         * coordinates of each grid box */
        for (m = 0; m < numNN; m++) {
            nnIndex[NDIM] = m;
            gridNN(nnIndex) = gIndex + nnShift(m);

            /* Here we enforce periodic boundary conditions where applicable,
             * otherwise we set the outlier to -1 */
            for (int i = 0; i < NDIM; i++) {
                if (gridNN(nnIndex)[i] == numNNGrid[i])
                    gridNN(nnIndex)[i] = static_cast<int>(boxPtr->periodic[i]-1);
                else if (gridNN(nnIndex)[i] == -1) 
                    gridNN(nnIndex)[i] = static_cast<int>(boxPtr->periodic[i]*numNNGrid[i] - 1);
            } // end i

        } // end m
    
    } // end n

    /* For small grid boxes, we need to make sure that periodic boundary conditions
     * don't cause double counting by having nearest neighbors to one side and the
     * other corresponding to the same box.  We simply set such neighbors to -1 and
     * have logic in the potential class to skip these. */
    iVec neg,dup;
    neg.fill(-1);
    /* We go through each grid box */
    for (int n = 0; n < totNumGridBoxes; n++) {
        gIndex = gridIndex(n);
        for (int i = 0; i < NDIM; i++) 
            nnIndex[i] = gIndex[i];

        /* Now, for each grid box, we start from the end of the NN list
         * and work backwards, setting any duplicates to -1 */
        for (m = numNN-1; m >= 0; m--) {
            nnIndex[NDIM] = m;
            dup = gridNN(nnIndex); 
            for (int p = m-1; p >= 0; p--) {
                nnIndex[NDIM] = p;
                /* Check for a duplicate */
                if ((gridNN(nnIndex) == dup)) 
                    gridNN(nnIndex) = neg;
            } // end p
        } // end m
    } //end n

    /* Now we copy over the NN indices to the reduced list before we remove the
     * possibility of double counting */
    gridNNReduced = gridNN;

    /* This is more complicated.  We need to eliminate any double counting.  We do 
     * this by following all links in the gridNN table one at a time, and making 
     * sure that we never have a back link.  i.e. if (10)->(01) then we can't allow
     * an entry for (01) which points back to (10). One needs to used gridNNReduced,
     * if you wish to compute some quantity by summing over all grid boxes. */

    /* We go through each grid box one at a time */
    for (int n = 0; n < totNumGridBoxes; n++) {
        gIndex = gridIndex(n);

        /* We only consider the 'unique' nearest neighbors */
        for (m = (numUniqueNN+1); m < numNN; m++) {

            for (int i = 0; i < NDIM; i++) 
                nnIndex[i] = gIndex[i];
            nnIndex[NDIM] = m;

            /* Get the coordinate of the nearest neighbor grid box */
            dup = gridNNReduced(nnIndex); 

            /* We only follow links for real grid boxes */
            if (noneEquals(dup, XXX)) {

                for (int i = 0; i < NDIM; i++) 
                    nnIndex[i] = dup[i];

                /* Go through the unique nn of the target box and make sure
                 * we don't have any entries pointing back to the reference
                 * box */
                for (int p = (numUniqueNN+1); p < numNN; p++) {
                    nnIndex[NDIM] = p;
                    if ((gridNNReduced(nnIndex) == gIndex))
                        gridNNReduced(nnIndex) = neg;
                } // end p

            } // end any

        } // end m
    } // end n


}

/**************************************************************************//**
 *  Print the NN Lookup table.
******************************************************************************/
void LookupTable::printGrid() {

    /* This outputs the NN grid table */
    for (int n = 0; n < totNumGridBoxes; n++) {
        iVec tIndex;
        tIndex = gridIndex(n);
        for (int i = 0; i < NDIM; i++) 
            std::cout << tIndex[i];
        std::cout << "  ";
    }
    std::cout << std::endl;

    for (int n = 0; n < totNumGridBoxes; n++) 
        std::cout << "===  ";
    std::cout << std::endl;

    for (int m = 0; m < numNN; m++) {
        nnIndex[NDIM] = m;
        for (int n = 0; n < totNumGridBoxes; n++) {
            iVec tIndex;
            tIndex = gridIndex(n);
            for (int i = 0; i < NDIM; i++) 
                nnIndex[i] = tIndex[i];

            for (int i = 0; i < NDIM; i++) {
                if (gridNNReduced(nnIndex)[i] == -1)
                    std::cout << "x";
                else
                    std::cout << gridNNReduced(nnIndex)[i];
            }
            std::cout << "  ";
        }
        std::cout << std::endl;
    }
}

/**************************************************************************//**
 *  Update the full nearest neighbor for a set of fixed particles, which are
 *  never updated and have the same position at all time slices.
******************************************************************************/
void LookupTable::updateGrid(const DynamicArray <dVec,1> &fixedPos) {

    numLabels.fill(0);
    beadLocator beadIndex;
    beadIndex[0] = 0;
    for (int n = 0; n < fixedPos.extents()[0]; ++n) {

        beadIndex[1] = n;

        /* First we figure out which grid box the particle is currently in */
        grid(beadIndex) = gridIndex(fixedPos(n));

        /* Get the current number of bead labels */
        int label = numLabels(numLabelIndex(beadIndex));

        /* We need to resize the hash array if we have too many beads in a single
         * grid box */
        if (hashSize[NDIM+1]==label) {
            hashSize[NDIM+1]++;
            hash.resizeAndPreserve(hashSize);
        }

        /* Update the hash table */
        hash(hashIndex(beadIndex,label)) = beadIndex[1];

        /* Assign the reverse lookup */
        beadLabel(beadIndex) = label;

        /* Increment the number of labels */
        numLabels(numLabelIndex(beadIndex))++;

    } // n
}

/**************************************************************************//**
 *  We update the full nearest neighbor grid by filling up all data arrays at
 *  all time slices.
******************************************************************************/
void LookupTable::updateGrid(const Path &path) {

    numLabels.fill(0);
    beadLocator beadIndex;
    for (beadIndex[0] = 0; beadIndex[0] < constants()->numTimeSlices(); ++beadIndex[0]) {
        for (beadIndex[1] = 0; beadIndex[1] < path.getNumParticles(); ++beadIndex[1]) {

            if (path.worm.beadOn(beadIndex)) {

                /* First we figure out which grid box the particle is currently in */
                grid(beadIndex) = gridIndex(path(beadIndex));

                /* Get the current number of bead labels */
                int label = numLabels(numLabelIndex(beadIndex));

                /* We need to resize the hash array if we have too many beads in a single
                 * grid box */
                if (hashSize[NDIM+1]==label) {
                    hashSize[NDIM+1]++;
                    hash.resizeAndPreserve(hashSize);
                }

                /* Update the hash table */
                hash(hashIndex(beadIndex,label)) = beadIndex[1];

                /* Assign the reverse lookup */
                beadLabel(beadIndex) = label;

                /* Increment the number of labels */
                numLabels(numLabelIndex(beadIndex))++;

            } //beadOn
        } // ptcl
    } // slice

}


/**************************************************************************//**
 *  Update the NN grid with a new bead position.
******************************************************************************/
void LookupTable::updateBead(const beadLocator &beadIndex, const dVec &pos) {

    /* Find out the new grid index */
    gIndex = gridIndex(pos);

    /* If the new position, is in the same grid box, we don't have to do anything */
    if (!(gIndex == grid(beadIndex))) {

        /* Delete the current bead from the grid */
        delBead(beadIndex);
        
        /* Add the new bead to the grid, using the already computed grid index. */

        /* Update the grid index */
        grid(beadIndex) = gIndex;

        /* Get the numLabel index */
        nI = numLabelIndex(beadIndex);

        /* Get the current number of bead labels */
        int label = numLabels(nI);

        /* We need to resize the hash array if we have too many beads in a single
         * grid box */
        if (hashSize[NDIM+1]==label) {
            hashSize[NDIM+1]++;
            hash.resizeAndPreserve(hashSize);
        }

        /* Update the hash table */
        hash(hashIndex(beadIndex,label)) = beadIndex[1];

        /* Assign the reverse lookup */
        beadLabel(beadIndex) = label;

        /* Increment the number of labels */
        numLabels(nI)++;
    }
}

/**************************************************************************//**
 *  Add a single bead to the NN grid and position pos.
******************************************************************************/
void LookupTable::addBead(const beadLocator &beadIndex, const dVec &pos) {

    /* First we figure out which grid box the particle is located in */

    grid(beadIndex) = gridIndex(pos);

    /* Get the numLabel index */
    nI = numLabelIndex(beadIndex);

    /* Get the current number of bead labels */
    int label = numLabels(nI);

    /* We need to resize the hash array if we have too many beads in a single
     * grid box */
    if (hashSize[NDIM+1]==label) {
        hashSize[NDIM+1]++;
        hash.resizeAndPreserve(hashSize);
    }

    /* Update the hash table */
    hash(hashIndex(beadIndex,label)) = beadIndex[1];

    /* Assign the reverse lookup */
    beadLabel(beadIndex) = label;

    /* Increment the number of labels */
    numLabels(nI)++;
}

/**************************************************************************//**
 *  Remove a single bead from the NN grid.
******************************************************************************/
void LookupTable::delBead(const beadLocator &beadIndex) {

    swap[0] = beadIndex[0];

    /* Get the current label of beadIndex, and the label of the last bead in 
     * the grid box. */ 
    int label = beadLabel(beadIndex);
    nI = numLabelIndex(beadIndex);
    int lastLabel = numLabels(nI)-1;

    /* Swap the last bead in the box, with the bead to be removed */
    hI = hashIndex(beadIndex,lastLabel);
    if (label < lastLabel) {
        swap[1] = hash(hI);
        hash(hI) = XXX;
        hI[NDIM+1] = label;
        hash(hI) = swap[1];
        beadLabel(swap) = label;
    }
    else
        hash(hI) = XXX;

    /* Reset the label and grid */
    beadLabel(beadIndex) = XXX;
    grid(beadIndex).fill(XXX);
    
    /* Decrement the number of labels */
    numLabels(nI)--;
}

/**************************************************************************//**
 *  Update the NN lookup table and the array of beadLocators containing
 *  all beads which 'interact' with the supplied bead1.
******************************************************************************/
void LookupTable::updateInteractionList(const Path &path, const beadLocator &bead1) {

    dVec sep;

    /* Reset the number of beads */
    numBeads = 0;

    /* Get the NDIM-vec coordiantes of the box where the bead resides */
    for (int i = 0; i < NDIM; i++)
        nnIndex[i] = grid(bead1)[i];

    /* We go through all the neighbouring grid boxes of grid(bead1) and
     * assemble the interaction list */
    beadLocator bead2;
    bead2[0] = bead1[0];
    for (int n = 0; n < numNN; n++) {
        nnIndex[NDIM] = n;

        /* Get the grid index of the interacting box */
        gIndex = gridNN(nnIndex);

        /* Make sure we don't access any illegal grid boxes */
        if (noneEquals(gIndex, XXX)) {

            int maxNL = numLabels(numLabelIndex(gIndex,bead1[0]));
            hI = hashIndex(gIndex,bead1[0],0);

            for (int label = 0; label < maxNL; label++) {

                /* Get the interacting bead */
                hI[NDIM+1] = label;
                bead2[1] = hash(hI);

                /* Eliminate self-interactions */
                if (!(bead1 == bead2)) {

                    sep = path.getSeparation(bead2,bead1);

                    /* If we are within the cutoff distance, add the bead to the list, 
                     * and store their separation.*/
                    if (dot(sep,sep) < rc2) {
                        beadList(numBeads) = bead2;
                        beadSep(numBeads) = sep;
                        numBeads++;
                    }
                } // bead2 != bead1

            } // label

        }  // skip any illegal grid boxes

    } // end n
}

/**************************************************************************//**
 *  Fill up the fullBeadList array with a list of beads in the same grid box
 *  as the supplied beadIndex and its nearest neighbors at the supplied time 
 *  slice.
******************************************************************************/
void LookupTable::updateFullInteractionList(const beadLocator &beadIndex, const int slice) {

    /* Copy over the location of the central grid box */
    for (int i = 0; i < NDIM; i++)
        nnIndex[i] = grid(beadIndex)[i];

    /* Now we loop over the central box, plus all nearest neighbors, filling
     * up the beadList and incrementing the numListBeads counter */
    fullBeadList.fill({XXX, XXX});
    fullNumBeads = 0;
    for (int nn = 0; nn < numNN; nn++) {
        nnIndex[NDIM] = nn;

        /* Get the grid index of the nearest neighbor box */
        gIndex = gridNN(nnIndex);

        if (noneEquals(gIndex, XXX)) {

            /* Get the hash table index, and max number of labels*/
            int maxNL = numLabels(numLabelIndex(gIndex,slice));
            hI = hashIndex(gIndex,slice,0);

            for (int label = 0; label < maxNL; label++) {

                /* Get the interacting bead */
                hI[NDIM+1] = label;
                fullBeadList(fullNumBeads) = {slice, hash(hI)};
                fullNumBeads++;

            } // label

        }  // skip out of bounds grid boxes

    } // end nn
}

/**************************************************************************//**
 *  Fill up the fullBeadList array with a list of beads in the supplied grid box
 *  indexed by its number and all its nearest neighbors at the supplied time
 *  slice.
******************************************************************************/
void LookupTable::updateFullInteractionList(const int gNumber, const int slice) {

    gIndex = gridIndex(gNumber);

    /* Copy over the location of the central grid box */
    for (int i = 0; i < NDIM; i++)
        nnIndex[i] = gIndex[i];

    /* Now we loop over the central box, plus all nearest neighbors, filling
     * up the beadList and incrementing the numListBeads counter */
    fullBeadList.fill({XXX, XXX});
    fullNumBeads = 0;
    for (int nn = 0; nn < numNN; nn++) {
        nnIndex[NDIM] = nn;

        /* Get the grid index of the nearest neighbor box */
        gIndex = gridNN(nnIndex);

        if (noneEquals(gIndex, XXX)) {

            /* Get the hash table index, and max number of labels*/
            int maxNL = numLabels(numLabelIndex(gIndex,slice));
            hI = hashIndex(gIndex,slice,0);

            for (int label = 0; label < maxNL; label++) {

                /* Get the interacting bead */
                hI[NDIM+1] = label;
                fullBeadList(fullNumBeads) = {slice, hash(hI)};
                fullNumBeads++;

            } // label

        }  // skip out of bounds grid boxes

    } // end nn
}

/**************************************************************************//**
 *  Given two beadIndices, determine if the beads lie in neighboring grid
 *  boxes.
******************************************************************************/
bool LookupTable::gridNeighbors(const beadLocator &bead1, const beadLocator &bead2) {

    /* Copy over the bead1 grid index coordinates */
    for (int i = 0; i < NDIM; i++) 
        nnIndex[i] = grid(bead1)[i];

    /* Get the bead2 grid index */
    gIndex = grid(bead2);

    /* Run over all nearest neighbors of bead1 and see if we have found bead2's
     * grid */
    for (int nn = 0; nn < numNN; nn++) {
        nnIndex[NDIM] = nn;
        if (gIndex == gridNN(nnIndex)) {
            return true;
        }
    }

    return false;
}
