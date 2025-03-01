/**
 * @file path.h
 * @author Adrian Del Maestro
 * @date 10.14.2008
 *
 * @brief Path class definition.
 */

#ifndef PATH_H 
#define PATH_H

#include "common.h"
#include "constants.h"
#include "container.h"
#include "worm.h"

class LookupTable;

// ========================================================================  
// Path Class
// ========================================================================  
/** 
 * The space-time trajectories.
 *
 * Holds the actual particle wordlines, consisting of a fixed number of 
 * time slices for each particle with periodic boundary conditions in 
 * imaginary time.
 */
class Path {

    public:
        Path (const Container *, LookupTable &, int, const DynamicArray<dVec,1>&, int numberBroken = 0);
        ~Path();

        /* Path* clone() const{ return new Path(*this); } */

        const int numTimeSlices;        ///< A local constant copy of the number of time slices
        int breakSlice;                 ///< The location of the break in the path (0=>no break)
        std::vector<int> brokenWorldlinesL;   ///< A list of particles with broken worldlines on left of break
        std::vector<int> brokenWorldlinesR;   ///< A list of particles with broken worldlines on right of break
        std::vector<int> closedWorldlines;   ///< A list of particles with closed worldlines on left of break

        const Container *boxPtr;        ///< A constant reference to the container class
        Worm worm;                      ///< Details on the worm

        LookupTable &lookup;            ///< A reference to the nearest neighbor lookup table.

	DynamicArray <int,1> numBeadsAtSlice;  ///< The number of active beads at a given time slice

        /** Get the size of the worldline array */
        int getNumParticles() const {return beads.extents()[1];}

        /** The number of active particles */
        int getTrueNumParticles() const {return ( worm.getNumBeadsOn() / numTimeSlices );}

        /** Operator Overloading to skip having to specifically grab .beads  */
        const dVec& operator() (int slice, int ptcl) const { 
            PIMC_ASSERT(slice>=0 && slice < numTimeSlices);
            return beads(slice,ptcl); }
        /** Operator Overloading to skip having to specifically grab .beads  */
        dVec& operator() (int slice, int ptcl) { 
            PIMC_ASSERT(slice>=0 && slice < numTimeSlices);
            return beads(slice,ptcl); }
        /** Operator Overloading to skip having to specifically grab .beads  */
        const dVec& operator() (const beadLocator &beadIndex) const { 
            PIMC_ASSERT(beadIndex[0]>=0 && beadIndex[0]<numTimeSlices);
            return beads(beadIndex); }
        /** Operator Overloading to skip having to specifically grab .beads  */
        dVec& operator() (const beadLocator &beadIndex) { 
            PIMC_ASSERT(beadIndex[0]>=0 && beadIndex[0]<numTimeSlices);
            return beads(beadIndex); }

        /** Return the velocity between two time slices of a given particle as a ndim-vector */
        dVec getVelocity(const beadLocator&) const;

        /** Return the separation std::vector between two particles in the same timeslice */
        dVec getSeparation(const beadLocator&, const beadLocator&) const;
        
        /** Return the pointer to the beads array data */
        auto get_beads_data_pointer() const;

        /** Return the extent of the beads array data */
        auto get_beads_extents() const;

        /** Output bead-link info, used for debugging */
        template<class Tstream> void printLinks(Tstream &);

        /** Output the world-line configurations in a generic format */
        void outputConfig(int) const;

        /** Move one link forward in imaginary time */
        beadLocator& next(int slice, int ptcl) {return nextLink(slice,ptcl);}
        /** Move one link forward in imaginary time */
        const beadLocator& next(int slice, int ptcl) const {return nextLink(slice,ptcl);}
        /** Move one link forward in imaginary time */
        beadLocator& next(const beadLocator &beadIndex) {return nextLink(beadIndex);}
        /** Move one link forward in imaginary time */
        const beadLocator& next(const beadLocator &beadIndex) const {return nextLink(beadIndex);}

        /** Move one link backward in imaginary time */
        beadLocator& prev(int slice, int ptcl) {return prevLink(slice,ptcl);}
        /** Move one link backward in imaginary time */
        const beadLocator& prev(int slice, int ptcl) const {return prevLink(slice,ptcl);}
        /** Move one link backward in imaginary time */
        beadLocator& prev(const beadLocator &beadIndex) {return prevLink(beadIndex);}
        /** Move one link backward in imaginary time */
        const beadLocator& prev(const beadLocator &beadIndex) const {return prevLink(beadIndex);}

        /** Move an integer number of links forward in imaginary time */
        beadLocator next(int,int,int) const;        
        beadLocator next(const beadLocator&,int) const; 

        /** Move an integer number of links backward in imaginary time */
        beadLocator prev(int,int,int) const;        
        beadLocator prev(const beadLocator&,int) const; 

        /** Add a bead to the worldline configuration at a given slice */
        beadLocator addBead(const int, const dVec &);
        /** Add a bead at the next time slice */
        beadLocator addNextBead(const beadLocator&, const dVec &);
        /** Add a bead at the previous time slice */
        beadLocator addPrevBead(const beadLocator&, const dVec &);

        /** Remove a bead from the world-line configuration */
        void delBead(const beadLocator&);
        /** Delete a bead and move forwards */
        beadLocator delBeadGetNext(const beadLocator&);
        /** Delete a bead and move backwards */
        beadLocator delBeadGetPrev(const beadLocator&);
    
        /** Break the link to right of bead **/
        void breakLink(const beadLocator&);
        /** Make a link between beads **/
        void makeLink(const beadLocator&,const beadLocator&);
        /** Break the link to right of bead t center slice AND update lists **/
        void removeCenterLink(const beadLocator&);
        /** Make a link between beads at center slice AND update lists **/
        void addCenterLink(const beadLocator&,const beadLocator&);
        /** Checks to see if worldline is broken**/
        bool isBroken(const beadLocator&) const;
        /** Returns factor for broken worldines**/
        double breakFactor(const beadLocator&,const beadLocator&) const;
        /** Checks to see if bead is in subregion A/B at break slice + 1 **/
        bool inSubregionA(const beadLocator&) const;
        bool inSubregionB(const beadLocator&) const;
        /** Check if only subregion worldlines are broken, for debugging **/
        bool checkSubregionLinks() const;

        /** Update the position of a bead in the worldine configuration */
        void updateBead(const beadLocator&, const dVec&);

        /** Used when debugging worm configurations */
        void printWormConfig(DynamicArray <beadLocator,1> &);

        /** Initialize any loaded state by left packing the array */
        void leftPack();
    
        /** Reset broken/closed worldline std::vectors **/
        void resetBrokenClosedVecs();

    private:
        friend class PathIntegralMonteCarlo;        // Friends for I/O

	DynamicArray<dVec,2> beads;                        // The wordline array
	DynamicArray<beadLocator,2> prevLink, nextLink;    // Bead connection matrices

        beadLocator lastBeadIndex;                  // Holds the index of the last bead on a slice
};

/* inline Path* new_clone(Path const& other){ */
/*     return other.clone(); */
/* } */

// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// INLINE FUNCTION DEFINITIONS
// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

/** Return the separation std::vector between two particles in the same timeslice */
inline dVec Path::getSeparation(const beadLocator &bead1, const beadLocator &bead2) const {
    dVec sep;
    sep = (*this)(bead1) - (*this)(bead2);
    boxPtr->putInBC(sep);
    return sep;
}

// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

/** Return the velocity between two time slices of a given particle as a ndim-vector */
inline dVec Path::getVelocity (const beadLocator &beadIndex) const {
    dVec vel;

    constexpr std::array<int,2> compareArr = { XXX, XXX };
    if (beadIndex==compareArr || next(beadIndex)==compareArr) {
        vel.fill(0.0);
        return (vel);
    }

    /* Calculate the 'velocity' and implement periodic boundary conditions */
    vel = beads(next(beadIndex)) - beads(beadIndex);
    boxPtr->putInBC(vel);

    return (vel);
}

// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

/** Return the pointer to the first element in the beads array */
inline auto Path::get_beads_data_pointer() const {
    return (*this).beads.data();
}

// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

/** Return the extents of the beads array */
inline auto Path::get_beads_extents() const {
    return (*this).beads.extents();
}

// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

/** Move an integer number of links forward in imaginary time */
inline beadLocator Path::next(int slice, int ptcl, int numLinks) const {
    PIMC_ASSERT(slice>=0 && slice<numTimeSlices && ptcl>=0);
    beadLocator bI;
    bI = slice,ptcl;
    for (int m = 0; m < numLinks; m++)
        bI = next(bI);
    return bI;
}

// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

/** Move an integer number of links forward in imaginary time */
inline beadLocator Path::next(const beadLocator &beadIndex, int numLinks) const {
    PIMC_ASSERT(beadIndex[0]>=0 && beadIndex[0]<numTimeSlices && beadIndex[1]>=0);
    beadLocator bI;
    bI = beadIndex;
    for (int m = 0; m < numLinks; m++)
        bI = next(bI);
    return bI;
}

// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

/** Move an integer number of links backward in imaginary time */
inline beadLocator Path::prev(int slice, int ptcl, int numLinks) const {
    PIMC_ASSERT(slice>=0 && slice<numTimeSlices && ptcl>=0);
    beadLocator bI;
    bI = slice,ptcl;
    for (int m = 0; m < numLinks; m++)
        bI = prev(bI);
    return bI;
}

// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

/** Move an integer number of links backward in imaginary time */
inline beadLocator Path::prev(const beadLocator &beadIndex, int numLinks) const {
    PIMC_ASSERT(beadIndex[0]>=0 && beadIndex[0]<numTimeSlices && beadIndex[1]>=0);
    beadLocator bI;
    bI = beadIndex;
    for (int m = 0; m < numLinks; m++)
        bI = prev(bI);
    return bI;
}

// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

/** Output bead-link info, used for debugging.*/
template<class Tstream>
void Path::printLinks(Tstream &outStream) {
    int numParticles = getNumParticles();
    for (int m = numTimeSlices-1; m >= 0; m--) {
        beadLocator beadIndex;
        for (int n = 0; n < numParticles; n++) {
            beadIndex = m,n;
            outStream << std::setw(2) << prevLink(beadIndex)[1] << " ";
        }
                        
        outStream << "\t";
        for (int n = 0; n < numParticles; n++) { 
            beadIndex = m,n;
            outStream << std::setw(2) << nextLink(beadIndex)[1] << " ";
        }
        outStream << std::endl;
    }
    outStream << std::endl;
}

#endif
