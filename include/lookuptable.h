/**
 * @file lookuptable.h 
 * @author Adrian Del Maestro
 * @date 03.16.2009
 *
 * @brief LookupTable class definition
 */

#ifndef LOOKUPTABLE_H 
#define LOOKUPTABLE_H 

#include "common.h"
#include "constants.h"
#include "container.h"

class Path;

// ========================================================================  
// LookupTable Class
// ========================================================================  
/**
 * The particle (bead) lookup table.
 *
 * Holds the lookup table and methods for finding where in the simulation
 * cell a particle is located, and which particles are close enough
 * to interact. Accomplishes this by partitioning the simulation cell
 * and storing the box at ever time slice where a bead is located.
 */
class LookupTable {

    public:
        LookupTable(const Container *, const int, const int);
        ~LookupTable();

        const Container *boxPtr;            ///< The simulation cell

        int numUniqueNN;                    ///< The number of unique nearest neighbors of each box
        int numNN;                          ///< The total number of nearest neighbors of each box

	blitz::Array <beadLocator,1> beadList;     ///< The cutoff dynamic list of interacting beads
	blitz::Array <beadLocator,1> fullBeadList; ///< The full dynamic list of interacting beads
	blitz::Array <dVec,1> beadSep;             ///< The separation between beads
        int numBeads;                       ///< The cutoff number of active beads in beadList;
        int fullNumBeads;                   ///< The full number of active beads in beadList;

        /** Return the number of NN grid boxes */
        iVec getNumNNGrid() {return numNNGrid;} 
        /** Return the total number of grid boxes */
        int getTotNumGridBoxes() {return totNumGridBoxes;}

        /* Update the NN table interaction list */
        void updateInteractionList(const Path &, const beadLocator &);
        void updateFullInteractionList(const beadLocator &, const int);
        void updateFullInteractionList(const int, const int);
        void updateGrid(const Path &);
        void updateGrid(const blitz::Array <dVec,1>&);

        void printGrid();

        /* Remove a bead from the grid */
        void delBead(const beadLocator &);
        void addBead(const beadLocator &, const dVec &);
        void updateBead(const beadLocator &, const dVec &);

        /* Returns the grid index where the suplied position resides */
        inline iVec gridIndex(const dVec &);

        /* Returns the grid number where the suplied position resides */
        inline int gridNumber(const dVec &);

        /* Converts a grid number to a grid index */
        inline iVec gridIndex(int);             

        /* Converts a grid index to a grid number */
        inline int gridNumber(const iVec &); 

        /** Resize the bead and grid lists */
        void resizeList(int _numParticles) {
            beadList.resize(_numParticles);
            fullBeadList.resize(_numParticles);
            beadSep.resize(_numParticles);
            grid.resizeAndPreserve(numLookupTimeSlices,_numParticles);
            beadLabel.resizeAndPreserve(numLookupTimeSlices,_numParticles);
        }

        /* Determine if two positions are in neighboring grid boxes */
        bool gridNeighbors(const beadLocator&, const beadLocator&);

        /** Determine if two beads are in the same grid box */
        bool gridShare(const beadLocator &bead1, const beadLocator &bead2) {
            return all(grid(bead1)==grid(bead2));
        }

    private:
        int numLookupTimeSlices;            // The number of time slices that 
        int totNumGridBoxes;                // The total number of nn grid boxes

        iVec numNNGrid;                     // The number of nn grid boxes in each direction
        iVec gIndex;                        // A commonly used grid box index vector

        double rc2;                         // A local copy of the potential cutoff squared

	std::array<int,NDIM+1> nnIndex;     // Comonly used nn index vector
	std::array<int,NDIM+1> nI;          // Used for indexing numLabels
	std::array<int,NDIM+2> hI;          // Used for indexing hash

	blitz::Array <iVec,NDIM+1> gridNN;         // The nearest neighbors of each grid box 
	blitz::Array <iVec,NDIM+1> gridNNReduced;  // The nn reduced to contain no back links

        dVec sizeNNGrid;                    // The size of the nn grid boxes in each direction

	blitz::Array <int,NDIM+2> hash;            // The main worldline lookup array
	blitz::Array <iVec,2> grid;                // The grid index of a bead
	blitz::Array <int,2> beadLabel;            // The label of a bead in a cell
	blitz::Array <int,NDIM+1> numLabels;       // The number of beads in a cell

	std::array <int,NDIM+2> hashSize;   // The current size of the hash array

        beadLocator swap;                   // Used when deleting/updating beads

        void setupNNGrid();                 // We initialize the nearest neighbor grid

        /* Return tiny vectors suitable for indexing the numLabel and hash array */
        inline std::array <int,NDIM+1> numLabelIndex(const beadLocator&);
        inline std::array <int,NDIM+1> numLabelIndex(const iVec&, const int);
        inline std::array <int,NDIM+2> hashIndex(const beadLocator&, const int);
        inline std::array <int,NDIM+2> hashIndex(const iVec&, const int, const int);
};

// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// INLINE FUNCTION DEFINITIONS
// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
/**
 * Given a particle position, return the grid index for the nearest
 * neighbor lookup table.
 */
inline iVec LookupTable::gridIndex(const dVec &pos) {
    iVec index;         // The grid index in each dimension
    for (int i = 0; i < NDIM; i++) {  
        index[i] = static_cast<int>( abs( pos[i] + 0.5 * boxPtr->side[i] - EPS ) 
                / (sizeNNGrid[i] + EPS) );
        PIMC_ASSERT(index[i]<numNNGrid[i]);
    }
    return index;
}

// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
/** 
 * Given a particle position, return the grid number for the nearest
 * neighbor lookup table. 
 */
inline int LookupTable::gridNumber(const dVec &pos) {
    int gNumber = 0;
    for (int i = 0; i < NDIM; i++) {  
        int scale = 1;
        for (int j = i+1; j < NDIM; j++) 
            scale *= numNNGrid[j];
        gNumber += scale * 
            static_cast<int>(abs( pos[i] + 0.5 * boxPtr->side[i] - EPS ) / (sizeNNGrid[i] + EPS));
    }
    PIMC_ASSERT(gNumber<totNumGridBoxes);
    return gNumber;
}

// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
/** 
 * Given the number of a grid box, returns its coordinates as a NDIM-vector.
 */
inline iVec LookupTable::gridIndex(int n) {
    iVec _gridIndex;
    for (int i = 0; i < NDIM; i++) {
        int scale = 1;
        for (int j = i+1; j < NDIM; j++) 
            scale *= numNNGrid[j];
        _gridIndex[i] = (n/scale) % numNNGrid[i];
        PIMC_ASSERT(_gridIndex[i]<numNNGrid[i]);
    }
    return _gridIndex;
}

// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
/**
 * Given the grid index of a box, return its grid number. 
 */
inline int LookupTable::gridNumber(const iVec &index) {
    int gNumber = 0;
    for (int i = 0; i < NDIM; i++) {
        int scale = 1;
        for (int j = i+1; j < NDIM; j++) 
            scale *= numNNGrid[j];
        gNumber += index[i]*scale;
    }
    PIMC_ASSERT(gNumber<totNumGridBoxes);
    return gNumber;
}

// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
/** 
 * Given a bead locator, return a NDIM+1 iVector used for indexing the numLabel
 * array. 
 */
inline std::array <int,NDIM+1> LookupTable::numLabelIndex(const beadLocator &beadIndex) {
	std::array<int,NDIM+1> index;
    for (int i = 0; i < NDIM; i++)
        index[i] = grid(beadIndex)[i];
    index[NDIM] = beadIndex[0];
    return index;
}

// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
/** 
 * Given a grid Index, and slice, return a NDIM+1 iVector used for indexing the numLabel
 * array.
 */
inline std::array <int,NDIM+1> LookupTable::numLabelIndex(const iVec &gI, const int slice) {
	std::array<int,NDIM+1> index;
    for (int i = 0; i < NDIM; i++)
        index[i] = gI[i];
    index[NDIM] = slice;
    return index;
}

// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
/** 
 * Given a bead locator and label, return a NDIM+2 iVector used for indexing the hash
 * array.
 */
inline std::array <int,NDIM+2> LookupTable::hashIndex(const beadLocator &beadIndex, 
        const int label) {
	std::array<int,NDIM+2> index;
    for (int i = 0; i < NDIM; i++)
        index[i] = grid(beadIndex)[i];
    index[NDIM]   = beadIndex[0];
    index[NDIM+1] = label;
    return index;
}

// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
/** 
 * Given a grid Index, slice and label, return a NDIM+2 iVector used for indexing the hash
 * array.
 */
inline std::array <int,NDIM+2> LookupTable::hashIndex(const iVec &gI, const int slice, 
        const int label) {
	std::array<int,NDIM+2> index;
    for (int i = 0; i < NDIM; i++)
        index[i] = gI[i];
    index[NDIM]   = slice;
    index[NDIM+1] = label;
    return index;
}

#endif
