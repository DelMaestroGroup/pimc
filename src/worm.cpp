/**
 * @file worm.cpp
 * @author Adrian Del Maestro
 * @brief Worm class implementation.
 */

#include "worm.h"
#include "path.h"

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// WORM CLASS ----------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/**************************************************************************//**
 *  Constructor.
 * 
 *  Initialize the worm object. 
 *  @param numParticles The number of particles
******************************************************************************/
Worm::Worm(int numParticles) {

    int numTimeSlices = constants()->numTimeSlices();

    /* Setup the bead array */
    beads.resize(numTimeSlices,numParticles);
    beads.fill(1u);

    /* Count the initial number of beads */
    resetNumBeadsOn();

    /* The initial configuration is always diagonal */
    isConfigDiagonal = true;

    /* Max worm cost [PRE 74, 036701 (2006)] */
    maxWormCost = 4.0;

    /* Initialize the properties of the worm */
    reset();
}

/**************************************************************************//**
 * Destructor.
******************************************************************************/
Worm::~Worm () {
}

/**************************************************************************//**
 *  Reset the worm to a null state.
******************************************************************************/
void Worm::reset() {
         gap = XXX;
      length =   0;
         sep.fill(0.0);
        head.fill(XXX);
        tail.fill(XXX);
    special1.fill(XXX);
    special2.fill(XXX);
}

/**************************************************************************//**
 *  We update all worm properties for a new head and tail.
 *
 *  Given new head and tail positions, we assign them and compute the gap, 
 *  length and separation, while turning off all special flags.
******************************************************************************/
void Worm::update(Path &path, const beadLocator &newHead, 
        const beadLocator &newTail) {

    /* Assign the new head and tail */
    head = newHead;
    tail = newTail;

    /* Compute the new worm gap.  This is defined to be the scalar separation
     * between the head and tail modulo the number of time slices*/
    gap = tail[0] - head[0];
    if (gap < 0)
        gap += path.numTimeSlices;
    PIMC_ASSERT(gap>=0);

    /* Get the new length */
    beadLocator beadIndex;
    beadIndex = tail;
    length = 0;
    do {
        ++length;
        beadIndex = path.next(beadIndex);
        /* cout << head[0] << " " << head[1] << " " << tail[0] << " " << tail[1] << " " << beadIndex[0] << " " << beadIndex[1] << endl; */
    } while (!(beadIndex == head));

    /* Now we update the head-tail separation */
    sep = path.getSeparation(tail,head);

    /* Unlink the head and tail */
    path.next(head).fill(XXX);
    path.prev(tail).fill(XXX);

    /* Turn off the special beads */
    special1.fill(XXX);
    special2.fill(XXX);
}

/**************************************************************************//**
 * Compute the value of the potential action trajectory factor.
 * 
******************************************************************************/
double Worm::factor(const beadState state1, const beadLocator &bead2) const {

    beadState state2 = getState(bead2);

    if ((state1 == NONE) && (state2 == NONE))
        return 1.0;
    else if ((state1 == HEADTAIL) && (state2 == HEADTAIL))
        return 0.0;
    else 
        return 0.5;
}

/**************************************************************************//**
 *  Test to see if a supplied bead is located on a worm.
******************************************************************************/
bool Worm::foundBead(const Path &path, const beadLocator &beadIndex) {
    if (isConfigDiagonal)
        return false;
    else {
        beadLocator beadID;
        beadID = tail;
        while ( (!(beadID == beadIndex)) && (!(beadID == path.next(head))) ) {
            beadID = path.next(beadID);
        }
        return beadID == beadIndex;
    }
}
