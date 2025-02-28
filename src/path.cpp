/**
 * @file path.cpp
 * @author Adrian Del Maestro
 * @date 10.14.2008
 * @brief Path class implementation.
 */

#include "path.h"
#include "lookuptable.h"
#include "communicator.h"

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// PATH CLASS ----------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/*************************************************************************//**
 * Constructor.
 *
 * We initialize the all the data structures in path given a random
 * initial configuration (classical) depending on the type of potential.
 * @param _boxPtr The pointer to the simulation cell
 * @param _lookup The nearest neighbor lookup table
 * @param _numTimeSlices The number of imaginary time slices
 * @param initialPos The initial configuration of particles
******************************************************************************/
Path::Path(const Container * _boxPtr, LookupTable &_lookup, int _numTimeSlices, 
        const DynamicArray<dVec,1> &initialPos,int numberBroken) :
    numTimeSlices(_numTimeSlices),
    boxPtr(_boxPtr), 
    worm(initialPos.size()),
    lookup(_lookup) 
{

    /* Get the initial number of particles */
    int numParticles = initialPos.size();

    /* Construct and initialize the array of beads */
    beads.resize(numTimeSlices,numParticles);
    beads = 0.0;

    /* Copy the initial condition at all time slices (a classical initial condition)*/
    for (int n = 0; n < numParticles; n++) {
        for (int i = 0; i < NDIM; i++)
	    beads.slice<1>(n)(i) = initialPos(n)[i];
    }

    /* Construct and initialize the prevLink and nextLink arrays */
    prevLink.resize(numTimeSlices,getNumParticles());
    nextLink.resize(numTimeSlices,getNumParticles());
    fill_with_function(prevLink.slice<0>(0), [](std::size_t j){ return static_cast<int>(j) - 1; });
    fill_with_function(prevLink.slice<0>(1), [](std::size_t j){ return static_cast<int>(j); });
    fill_with_function(nextLink.slice<0>(0), [](std::size_t j){ return static_cast<int>(j) + 1; });
    fill_with_function(nextLink.slice<0>(1), [](std::size_t j){ return static_cast<int>(j); });

    
    if (PIGS) {
    /* Here we implement the fixed boundary conditions in imaginary time. */
	fill_mdspan(prevLink.slice<0>(0), XXX);
	fill_mdspan(nextLink.slice<0>(numTimeSlices-1), XXX);

        /* Here we break worldlines at the center of the path if requested*/
        breakSlice = 0;
        if (numberBroken > 0){
            breakSlice = (numTimeSlices-1)/2;
            for (int n=0; n<numberBroken; n++){
                nextLink(breakSlice,n) = XXX;
                prevLink(breakSlice+1,n) = XXX;
                brokenWorldlinesL.push_back(n);
                brokenWorldlinesR.push_back(n);
            }
            for (int n=numberBroken; n<getNumParticles(); n++)
                closedWorldlines.push_back(n);
        }else if(constants()->spatialSubregionOn()){
            breakSlice = (numTimeSlices-1)/2;
            beadLocator beadIndex;
            beadIndex[0] =breakSlice+1;
            for( int n=0; n<getNumParticles(); n++){
                beadIndex[1] = n;
                if (inSubregionA(beadIndex)){
                    nextLink(breakSlice,n) = XXX;
                    prevLink(breakSlice+1,n) = XXX;
                    brokenWorldlinesL.push_back(n);
                    brokenWorldlinesR.push_back(n);
                }else{
                    closedWorldlines.push_back(n);
                }
            }
        }else{
            for (int n=0; n<getNumParticles(); n++)
                closedWorldlines.push_back(n);
        }
    }
    else {
        breakSlice = 0;
        /* Here we implement periodic boundary conditions in imaginary time */
	prevLink.slice<0>(0)(0) = numTimeSlices-1;
	nextLink.slice<0>(numTimeSlices-1)(0) = 0;
    }

    /* Initialize the number of active beads at each slice */
    numBeadsAtSlice.resize(numTimeSlices);
    numBeadsAtSlice = getNumParticles();

    /* Initialize the lookup table */
    lookup.resizeList(getNumParticles());
    lookup.updateGrid(*this);
}

/*************************************************************************//**
 * Destructor.
 *
 * Kill all blitz arrays
*****************************************************************************/
Path::~Path () {
    beads.free();
    prevLink.free();
    nextLink.free();
    numBeadsAtSlice.free();
}

/*************************************************************************//**
 * Initialize any loaded state by left packing the array.
 *
 * Go through all data arrays, and make sure the active beads are packed at
 * the beginning (LHS) of the array.
******************************************************************************/
void Path::leftPack() {

    beadLocator bead1,bead2;

    /* We go through each slice, and make sure the data arrays are left packed */
    for (bead1[0] = 0; bead1[0] < numTimeSlices; bead1[0]++) {
        for (bead1[1] = 0; bead1[1] < beads.extent(blitz::secondDim); bead1[1]++) {
            if (!worm.beadOn(bead1)) {
                bead2 = bead1;

                /* Find an active bead to the right of the inactive bead */
                bool foundBead = false;
                for (bead2[1] = bead1[1] + 1; bead2[1] < beads.extent(blitz::secondDim); bead2[1]++) {
                    if (worm.beadOn(bead2)) {
                        foundBead = true;
                        break;
                    }
                } //bead2
                
                /* If we have found one, perform the swap */
                if (foundBead) {

                    /* Copy the position back */
                    beads(bead1) = beads(bead2);

                    /* Swap the active bead indicators */
                    worm.addBead(bead1);
                    worm.delBead(bead2);

                    /* Update the links */
                    next(bead1) = next(bead2);
                    prev(bead1) = prev(bead2);
                    next(prev(bead2)) = bead1;
                    prev(next(bead2)) = bead1;

                    /* Zero out the old links */
                    next(bead2) = XXX;
                    prev(bead2) = XXX;

                } // foundBead

            } // bead1On
        } //bead11
    } //bead10
}

/**************************************************************************//**
 * Update the position of a bead in the worldine configuration.
 *
 * Given a beadIndex and spatial position, we assign the new bead location
 * and update its position in the lookup table.
******************************************************************************/
void Path::updateBead(const beadLocator &beadIndex, const dVec &pos) {
    lookup.updateBead(beadIndex,pos);
    beads(beadIndex) = pos;
}

/**************************************************************************//**
 * Add a bead at the next time slice.
 *
 * We add a single bead to our configuration at one advanced time slice from
 * the supplied bead.
 * @return next(prevIndex)
******************************************************************************/
beadLocator Path::addNextBead(const beadLocator &prevIndex, const dVec &pos) {

    beadLocator beadIndex;

    /* We make sure that the next bead doesn't already exist */
    PIMC_ASSERT(all(next(prevIndex)==XXX));

    /* The bead doesn't exist, so add it */ 
    int slice = prevIndex[0] + 1;
    if (slice >= numTimeSlices)
        slice -= numTimeSlices;

    beadIndex = addBead(slice,pos);

    next(prevIndex) = beadIndex;
    prev(beadIndex) = prevIndex;

    return beadIndex;
}

/**************************************************************************//**
 * Add a bead at the previous time slice.
 *
 * We add a single bead to our configuration at one previous time slice from
 * the supplied bead.
 * @return prev(nextIndex)
******************************************************************************/
beadLocator Path::addPrevBead(const beadLocator &nextIndex, const dVec &pos) {

    /* We make sure that the previous bead doesn't already exist */
    PIMC_ASSERT(all(prev(nextIndex)==XXX));

    /* The bead doesn't exist, so add it */ 
    int slice = nextIndex[0] - 1;
    if (slice < 0)
        slice += numTimeSlices;

    beadLocator beadIndex;
    beadIndex = addBead(slice,pos);

    prev(nextIndex) = beadIndex;
    next(beadIndex) = nextIndex;

    return beadIndex;
}

/**************************************************************************//**
 *  Add a bead to the worldline configuration at a given slice.
 * 
 *  We add a single bead to our configuration by searching for a free position
 *  at the end of the wordline array. If no such position exists, all data
 *  structures are resized to accomidate the new bead;
 *  @param slice The time slice where the new bead will be added.
 *  @param pos The new spatial position of the bead
 *  @return a beadLocator for the new bead
******************************************************************************/
beadLocator Path::addBead(const int slice, const dVec &pos) {

    int numWorldLines = getNumParticles();

    lastBeadIndex[0] = slice;
    lastBeadIndex[1] = numBeadsAtSlice(slice);

    /* Here we check and see if we have enough free-space to add a bead.  If
     * not we have to grow all our data structures */
    if (lastBeadIndex[1] == numWorldLines) {

        /* Resize and initialize the main data array which holds all worldline 
         * configurations */
        beads.resizeAndPreserve(numTimeSlices,numWorldLines + 1);
	fill_mdspan(beads.slice<1>(numWorldLines), 0.0);

        /* Resize and initialize the worm bead arrays which tells us
         * whether or not we have a bead present */
        worm.beads.resizeAndPreserve(numTimeSlices,numWorldLines + 1);
	fill_mdspan(worm.beads.slice<1>(numWorldLines), 0);

        /* Resize and initialize the previous and next link arrays */
        prevLink.resizeAndPreserve(numTimeSlices,numWorldLines + 1);
        nextLink.resizeAndPreserve(numTimeSlices,numWorldLines + 1);
	fill_mdspan(prevLink.slice<1>(numWorldLines), XXX);
	fill_mdspan(nextLink.slice<1>(numWorldLines), XXX);

        /* Resize the lookup table */
        lookup.resizeList(numWorldLines + 1);
    } 

    PIMC_ASSERT(!worm.beadOn(lastBeadIndex));

    /* Add and increment the number of beads */
    worm.addBead(lastBeadIndex);
    worm.incNumBeadsOn();

    /* Update the number of beads at the time slice */
    ++numBeadsAtSlice(slice);

    /* Initialize the connections */
    next(lastBeadIndex) = XXX;
    prev(lastBeadIndex) = XXX;

    /* Update the actual path and lookup table */
    beads(lastBeadIndex) = pos;
    lookup.addBead(lastBeadIndex,pos);

    /* Return the added bead */
    return lastBeadIndex;
}


/**************************************************************************//**
 *  Delete a bead and move backwards.
 * 
 *  Delete a single bead and return the index of the preivous bead in 
 *  imaginary time.
******************************************************************************/
beadLocator Path::delBeadGetPrev(const beadLocator &beadIndex) {

    /* Store the previous index */
    beadLocator prevBead;
    prevBead = prev(beadIndex);

    /* Remove the bead, and return the previous index */
    delBead(beadIndex);
    return prevBead;
}

/**************************************************************************//**
 *  Delete a bead and move forwards.
 * 
 *  Delete a single bead and return the index of the next bead in imaginary 
 *  time.
******************************************************************************/
beadLocator Path::delBeadGetNext(const beadLocator &beadIndex) {

    /* Store the next index */
    beadLocator nextBead;
    nextBead = next(beadIndex);

    /* Remove the bead, and return the next index */
    delBead(beadIndex);
    return nextBead;
}

/**************************************************************************//**
 *  Remove a bead from the world-line configuration.
 * 
 *  We delete a single bead from our configuration and update the emptyWorldLine
 *  array. At the same point, we must update the nearest neighbor lookup 
 *  table.
******************************************************************************/
void Path::delBead(const beadLocator &beadIndex) {

    /* Delete the bead from the lookup table */
    lookup.delBead(beadIndex);

    /* Reduce the number of beads at this time slice by one. */
    --numBeadsAtSlice(beadIndex[0]);

    /* Get the coordinates of the last bead */
    lastBeadIndex[0] = beadIndex[0];
    lastBeadIndex[1] = numBeadsAtSlice(beadIndex[0]);

    /* unlink */
    if (!all(next(beadIndex)==XXX))
        prev(next(beadIndex)) = XXX;
    if (!all(prev(beadIndex)==XXX))
        next(prev(beadIndex)) = XXX;

    /* If we are not already the largest bead label, perform the value
     * and linkage swap. */
    if (beadIndex[1] < numBeadsAtSlice(beadIndex[0])) {

        /* Copy over the position, updating the lookup table */
        lookup.delBead(lastBeadIndex);
        lookup.addBead(beadIndex,beads(lastBeadIndex));
        beads(beadIndex) = beads(lastBeadIndex);

        /* Swap the next/prev links */
        prev(beadIndex) = prev(lastBeadIndex);
        next(beadIndex) = next(lastBeadIndex);

        if (!all(next(lastBeadIndex)==XXX))
            prev(next(lastBeadIndex)) = beadIndex;
        if (!all(prev(lastBeadIndex)==XXX))
            next(prev(lastBeadIndex)) = beadIndex;

        /* We have to make sure that the worm is updated correctly if the swapped bead
         * is either a head or tail */
        if (all(worm.head==lastBeadIndex))
            worm.head = beadIndex;
        if (all(worm.tail==lastBeadIndex))
            worm.tail = beadIndex;
        if (all(worm.special1==lastBeadIndex))
            worm.special1 = beadIndex;
        if (all(worm.special2==lastBeadIndex))
            worm.special2 = beadIndex;
    }

    /* Unlink */
    next(lastBeadIndex) = XXX;
    prev(lastBeadIndex) = XXX;

    /* delete from the beads array */
    worm.delBead(lastBeadIndex);

    /* decrement the number of beads */
    worm.decNumBeadsOn();
}

/**************************************************************************//**
 *  Break a link to right of bead.
 ******************************************************************************/
void Path::breakLink(const beadLocator &beadIndexL) {
    beadLocator beadIndexR = next(beadIndexL);
    nextLink(beadIndexL[0],beadIndexL[1]) = XXX;
    prevLink(beadIndexR[0],beadIndexR[1]) = XXX;
}

/**************************************************************************//**
*  Make a link between beads.
******************************************************************************/
void Path::makeLink(const beadLocator &beadIndexL,const beadLocator &beadIndexR) {
    nextLink(beadIndexL[0],beadIndexL[1]) = beadIndexR;
    prevLink(beadIndexR[0],beadIndexR[1]) = beadIndexL;
}

/**************************************************************************//**
*  Break a link to right of bead AND update lists.
******************************************************************************/
void Path::removeCenterLink(const beadLocator &beadIndexL) {
    /* Get linked bead */
    beadLocator beadIndexR = next(beadIndexL);
    
    /* Break link */
    breakLink(beadIndexL);
    
    /* Update lists */
    std::vector<int>::iterator itr;
    itr = find(closedWorldlines.begin(), closedWorldlines.end(), beadIndexL[1]);
    closedWorldlines.erase(itr);
    brokenWorldlinesL.push_back(beadIndexL[1]);
    brokenWorldlinesR.push_back(beadIndexR[1]);
}

/**************************************************************************//**
*  Make a link between beads AND update lists
******************************************************************************/
void Path::addCenterLink(const beadLocator &beadIndexL,const beadLocator &beadIndexR) {
    
    /* Make new link */
    makeLink(beadIndexL,beadIndexR);
    
    /* Update lists */
    std::vector<int>::iterator itr;
    itr = find(brokenWorldlinesL.begin(), brokenWorldlinesL.end(), beadIndexL[1]);
    brokenWorldlinesL.erase(itr);
    itr = find(brokenWorldlinesR.begin(), brokenWorldlinesR.end(), beadIndexR[1]);
    brokenWorldlinesR.erase(itr);
    closedWorldlines.push_back(beadIndexL[1]);
}

/**************************************************************************//**
*  Reset broken bead lists.
******************************************************************************/

void Path::resetBrokenClosedVecs(){

    if ( PIGS && (breakSlice > 0) ){
        
        beadLocator beadIndex;
        
        /* Clear std::vectors */
        brokenWorldlinesL.clear();
        brokenWorldlinesR.clear();
        closedWorldlines.clear();
        
        /* Set brokenWorldlinesL and closedWorldlines by checking breakSlice */
        beadIndex[0] = breakSlice;
        for (int ptcl = 0; ptcl < numBeadsAtSlice(beadIndex[0]); ptcl++) {
            beadIndex[1] = ptcl;
            if( all(next(beadIndex) == XXX))
                brokenWorldlinesL.push_back(ptcl);
            else
                closedWorldlines.push_back(ptcl);
        }
        /* Set brokenWorldlinesR by checking breakSlice+1 */
        beadIndex[0] = breakSlice+1;
        for (int ptcl = 0; ptcl < numBeadsAtSlice(beadIndex[0]); ptcl++) {
            beadIndex[1] = ptcl;
            if( all(prev(beadIndex) == XXX))
                brokenWorldlinesR.push_back(ptcl);
        }
    }
};

/**************************************************************************//**
*  Check to see if worldline is broken at beadIndex
******************************************************************************/
bool Path::isBroken(const beadLocator &beadIndex) const{
    bool broken = false;
    if ( all(prev(beadIndex)==XXX) || all(next(beadIndex)==XXX) )
        broken = true;
    return broken;
}

/**************************************************************************//**
*  Return a constatnt factor for worldine breaks
******************************************************************************/
double Path::breakFactor(const beadLocator &beadIndex1,
                         const beadLocator &beadIndex2) const{
    double factor = 1.0;

    /*if ( PIGS ){
        if ( (breakSlice > 0) && (beadIndex1[0] == (breakSlice+1)) ){
            if ( isBroken(beadIndex1) || isBroken(beadIndex2) )
                factor = 0.5;
        }
    }*/
    return factor;
}

/**************************************************************************//**
*  Check if bead is in subregion A
******************************************************************************/
bool Path::inSubregionA(const beadLocator &beadIndex) const{
    bool inA = false;
    if( constants()->spatialSubregionOn() ){
        if ( (beadIndex[0]==breakSlice+1)&&
                //(abs(beads(beadIndex)[0]) < constants()->spatialSubregion()) )
                ( beads(beadIndex)[0] < constants()->spatialSubregion() ) )
            inA = true;
    }
    //cout << beadIndex[0] << '\t' << inA << '\t' << beads(beadIndex) << std::endl;
    return inA;
}

/**************************************************************************//**
*  Check if bead is in subregion B
******************************************************************************/
bool Path::inSubregionB(const beadLocator &beadIndex) const{
    bool inB = false;
    if( constants()->spatialSubregionOn() ){
        if ( (beadIndex[0]==breakSlice+1)&&(!inSubregionA(beadIndex)) )
            inB = true;
    }
    return inB;
}

/**************************************************************************//**
*  Check if bead is in subregion B
******************************************************************************/
bool Path::checkSubregionLinks() const{
    bool foundError = false;
    beadLocator beadIndex;
    
    beadIndex[0] = breakSlice+1;
    for(int n=0; n<getNumParticles(); n++){
        beadIndex[1] = n;
        if( inSubregionA(beadIndex)){
            foundError = (all(prev(beadIndex)!= XXX));
        }else if( inSubregionB(beadIndex) ){
            foundError = (all(prev(beadIndex)== XXX));
        }else{
            foundError=true;
        }
        if(foundError){
            std::cout << beadIndex[1] << '\t' << beads(beadIndex) << 't' <<
            inSubregionA(beadIndex) << '\t' << inSubregionB(beadIndex) << std::endl;
            break;
        }
    }
    return foundError;
}


/**************************************************************************//**
 *  Output the world-line configurations in a generic format.
 * 
 *  Output the worldline configuration to disk using a format suitable for 
 *  plotting. 
******************************************************************************/
#include "communicator.h"
void Path::outputConfig(int configNumber) const {

    int numParticles = getNumParticles();
    /* Output format for PIGS and PIMC are different */
#if PIGS
    /* We go through each particle/worldline */
    beadLocator beadIndex;

    /* Output the header */
    communicate()->file("wl")->stream() << format("# START_CONFIG %06d\n") % configNumber;

    /* Output the unit cell information.  It is always cubic.  Everything is scaled by
     * an overall factor for better visualization. */

    /* We output the bead block */
    for (int n = 0; n < numParticles;  n++) {
        for (int m = 0; m < numTimeSlices; m++) {
            beadIndex = m,n;
        
            communicate()->file("wl")->stream() << format("%8d %8d %8d") % beadIndex[0] 
            % beadIndex[1] % 1;

            /* Output the coordinates in 3D */
            int i;
            for (i = 0; i < NDIM; i++) {
                communicate()->file("wl")->stream() << format("%16.3E") % (beads(beadIndex)[i]);
            }
            while (i < 3) {
                communicate()->file("wl")->stream() << format("%16.3E") % 0.0;
                i++;
            }

            /* Output the bead indices of the connecting beads */
            communicate()->file("wl")->stream() << format("%8d %8d %8d %8d\n") % prev(beadIndex)[0] 
                % prev(beadIndex)[1] % next(beadIndex)[0] % next(beadIndex)[1];

            /* Advance the bead index */
            beadIndex = next(beadIndex);
        } 
    }
    communicate()->file("wl")->stream() << format("# END_CONFIG %06d\n") % configNumber;

    /* Flush the file stream */
    communicate()->file("wl")->stream().flush();
#else
    /* We go through all beads, and find the start and end bead for each
     * worldline, adding them to an array */
    DynamicArray <beadLocator,1> startBead,endBead;
    startBead.resize(numParticles);
    endBead.resize(numParticles);

    /* We sort the output by the number of beads in a worldline */
    DynamicArray <int,1> wlLength(numParticles);
    wlLength = 0;

    int numWorldLines = 0;

    /* Get the list of beads that are active in the simulation */
    DynamicArray <bool,2> doBead(numTimeSlices,numParticles);      
    doBead = blitz::cast<bool>(worm.getBeads());

    /* We go through each particle/worldline */
    int nwl = 0;
    beadLocator beadIndex;

    /* If we are off-diagonal, we start with the worm */
    if (!worm.isConfigDiagonal) {
        startBead(nwl) = worm.tail;
        endBead(nwl)   = XXX;

        /* Mark the beads as touched and increment the number of worldlines */
        beadIndex = startBead(nwl);
        do {
            doBead(beadIndex) = false;
            beadIndex = next(beadIndex);
        } while (!blitz::all(beadIndex==endBead(nwl)));

        /* We label a worm with a XXX */
        wlLength(nwl) = XXX;

        nwl++;
    } // off-diagonal

    /* Now go through eacheach worldline, and figure out how many particles
     * are involved in the permuation cycle */
    beadLocator testStart;
    for (int n = 0; n < numParticles; n++) {

        /* The initial bead to be moved */
        testStart = 0,n;

        /* We make sure we don't try to touch the same worldline twice */
        if (doBead(testStart)) {

            startBead(nwl) = testStart;

            /* Otherwise, we loop around until we find the initial bead */
            endBead(nwl) = startBead(nwl);

            /* Mark the beads as touched and increment the number of worldlines */
            beadIndex = startBead(nwl);
            int length = 1;
            do {
                doBead(beadIndex) = false;
                length++;
                beadIndex = next(beadIndex);
            } while (!blitz::all(beadIndex==endBead(nwl)));

            /* We label each trajectory by the number of particles it contains. */
            wlLength(nwl) = int(length/numTimeSlices)-1;

            nwl++;
        } // touchBeads

    } // n

    numWorldLines = nwl;

    /* Output the header */
    communicate()->file("wl")->stream() << format("# START_CONFIG %06d\n") % configNumber;

    /* Output the unit cell information.  It is always cubic.  Everything is scaled by
     * an overall factor for better visualization. */

    /* We output the bead block */
    for (int n = 0; n < numWorldLines;  n++) {
        beadIndex = startBead(n);
        do {
            communicate()->file("wl")->stream() << format("%8d %8d %8d") % 
                beadIndex[0] % beadIndex[1] % wlLength(n);

            /* Output the coordinates in 3D */
            int i;
            for (i = 0; i < NDIM; i++) {
                communicate()->file("wl")->stream() << format("%16.3E") % (beads(beadIndex)[i]);
            }
            while (i < 3) {
                communicate()->file("wl")->stream() << format("%16.3E") % 0.0;
                i++;
            }

            /* Output the bead indices of the connecting beads */
                communicate()->file("wl")->stream() << format("%8d %8d %8d %8d\n") % prev(beadIndex)[0] 
                    % prev(beadIndex)[1] % next(beadIndex)[0] % next(beadIndex)[1];

            /* Advance the bead index */
            beadIndex = next(beadIndex);
        } while (!all(beadIndex==endBead(n)));
    }
    communicate()->file("wl")->stream() << format("# END_CONFIG %06d\n") % configNumber;

    /* Flush the file stream */
    communicate()->file("wl")->stream().flush();

    /* Free up memory */
    startBead.free();
    endBead.free();
    wlLength.free();
    doBead.free();
#endif
}

/**************************************************************************//**
 *  Used when debugging worm configurations.
 * 
 *  We print out the current link and bead configuration to visualize a worm.
******************************************************************************/
void Path::printWormConfig(DynamicArray <beadLocator,1> &wormBeads) {

    int numWorldLines = getNumParticles();

    /* A shortform for the output file */
    std::fstream *outFilePtr;            
    outFilePtr = &(communicate()->file("debug")->stream());

    for (int m = numTimeSlices-1; m >= 0; m--) {
        /* First print out the beads */
        for (int n = 0; n < numWorldLines; n++) {
            beadLocator beadIndex;
            beadIndex = m,n;
            if (all(beadIndex==worm.head)) {
                if (worm.beadOn(beadIndex))
                    (*outFilePtr) << "^";
                else
                    (*outFilePtr) << "z";
            }
            
            else if (all(beadIndex==worm.tail)) {
                if (worm.beadOn(beadIndex))
                    (*outFilePtr) << "v";
                else
                    (*outFilePtr) << "y";
            }
            else
            {
                if (worm.beadOn(beadIndex))
                    (*outFilePtr) << "*";
                else
                    (*outFilePtr) << " ";
            }
        } // for n
                        
        (*outFilePtr) << "\t";
        /* Now print out the links */
        for (int n = 0; n < numWorldLines; n++) { 
            beadLocator beadIndex;
            beadIndex = m,n;
            if (all(next(beadIndex)==XXX))
                (*outFilePtr) << " ";
            else
                (*outFilePtr) << "|";
        } // for n

        (*outFilePtr) << "\t";
        /* Finally, if we are off-diagonal, print out the worm */
        if (!worm.isConfigDiagonal) {
            for (int n = 0; n < numWorldLines; n++) {
                beadLocator beadIndex;
                beadIndex = m,n;
                bool foundWormBead = false;
                std::string beadOut;
                for (int k = 0; k < wormBeads.extent(blitz::firstDim); k++) { 
                    
                    if (all(beadIndex==wormBeads(k))) { 
                        foundWormBead = true;
                        if ((k > 0) && (k < (wormBeads.extent(blitz::firstDim)-1))) {
                            if ( (wormBeads(k+1)[1] > wormBeads(k)[1]) || 
                                    (wormBeads(k-1)[1] < wormBeads(k)[1]) ) 
                                beadOut = "/";
                            else if ( (wormBeads(k+1)[1] < wormBeads(k)[1]) || 
                                    (wormBeads(k-1)[1] > wormBeads(k)[1]) ) 
                                beadOut = "\\";
                            else
                                beadOut = "|";
                        }
                        break;
                    }
                }
 
                if (foundWormBead) {
                    if (all(beadIndex==worm.head))
                        (*outFilePtr) << "^";
                    else if (all(beadIndex==worm.tail))
                        (*outFilePtr) << "v";
                    else
                        (*outFilePtr) << beadOut;
                }
                else {
                    if (worm.beadOn(beadIndex))
                        (*outFilePtr) << "~";
                    else
                        (*outFilePtr) << " ";
                }
            } // for n

        } // isConfigDiagonal
        (*outFilePtr) << std::endl;

    } // for int m
    (*outFilePtr) << std::endl;
}
