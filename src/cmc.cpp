/**
 * @file cmc.cpp
 * @brief Classical Monte Carlo simulation
 * @author Adrian Del Maestro
 * @date 08.13.2012
 */

#include "common.h"
#include "constants.h"
#include "container.h"
#include "potential.h"
#include "cmc.h"
#include "communicator.h"

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// CLASSICAL MONTE CARLO CLASS -----------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/**************************************************************************//**
 * Constructor.
******************************************************************************/
ClassicalMonteCarlo::ClassicalMonteCarlo (PotentialBase *_externalPtr,
        PotentialBase *_interactionPtr, MTRand &_random, const Container *_boxPtr,
        DynamicArray <dVec,1> &initialPos) :
    externalPtr(_externalPtr), 
    interactionPtr(_interactionPtr),
    random(_random),
    boxPtr(_boxPtr),
    config(initialPos)
{
    /* The number of particles */
    numParticles = config.extents()[0];


    /* Set the fugacity z*/
    z = exp(constants()->mu()/constants()->T())/pow(constants()->dBWavelength(),NDIM);

    /* Compute the initial total potential energy */
    energy = getTotalEnergy();

    /* Initialize acceptance tracking */
    numMoveTotal = 0;
    numMoveAccept = 0;

    numDeleteTotal = 0;
    numDeleteAccept = 0;

    numInsertTotal = 0;
    numInsertAccept = 0;

    /* Initialize measurements */
    aveEnergy = 0.0;
    aveNumParticles = 0.0;
    aveEoN = 0.0;
}

/**************************************************************************//**
 * Destructor. 
******************************************************************************/
ClassicalMonteCarlo::~ClassicalMonteCarlo ()
{
    // empty destructor
}

/**************************************************************************//**
 * Compute the total energy.
******************************************************************************/
double ClassicalMonteCarlo::getTotalEnergy()
{
    double locEnergy = 0.0;
    sep = 0.0;
    for (int part1 = 0; part1 < numParticles; part1++) {
        locEnergy += externalPtr->V(config(part1));

        for (int part2 = part1+1; part2 < numParticles; part2++) {
            sep = config(part1)-config(part2);
            boxPtr->putInside(sep);
            locEnergy += interactionPtr->V(sep);
        }
    }
    return locEnergy;
}

/**************************************************************************//**
 * Perform the monte carlo simulation.
******************************************************************************/
void ClassicalMonteCarlo::run(uint32 numMCSteps, bool gce) {

    int numMeasure = 0;
    double x;
    for(uint32 n = 1; n < numMCSteps; n++) {
        int m = 0;
        do {
            if (gce)
                x = random.rand();
            else 
                x = 0.0;

            /* Perform updates */
            if (x < 1.0/3.0)
                moveParticle();
            else if (x < 2.0/3.0)
                deleteParticle();
            else
                insertParticle();

            /* Update observables */
            aveEnergy += energy;
            aveNumParticles += numParticles;
            aveEoN += energy/(1.0*numParticles);
            numMeasure++;

            m++;
        } while (m < numParticles);

        // if ((n % 50) == 0) 
        //    measure(numMeasure);
    }

    /* Update the config array */
    config.resizeAndPreserve(numParticles);
}

/**************************************************************************//**
 * Perform a simple positional update
 *
******************************************************************************/
void ClassicalMonteCarlo::moveParticle() {

    dVec oldPos;
    oldPos = 0.0;
    double oldV,newV;

    numMoveTotal++;

    int p = random.randInt(numParticles-1);

    /* Compute the old energy of particle p*/
    oldV = externalPtr->V(config(p));
    for (int p2 = 0; p2 < numParticles; p2++) {
        if (p != p2) {
            sep = config(p)-config(p2);
            boxPtr->putInside(sep);
            oldV += interactionPtr->V(sep);
        }
    }

    oldPos = config(p);

    /* The new random position */
    config(p) = boxPtr->randUpdate(random,oldPos);

    /* Compute the new energy of particle p*/
    newV = externalPtr->V(config(p));
    for (int p2 = 0; p2 < numParticles; p2++) {
        if (p != p2) {
            sep = config(p)-config(p2);
            boxPtr->putInside(sep);
            newV += interactionPtr->V(sep);
        }
    }

    deltaV = newV - oldV;

    /* Now the metropolis step */
    if (random.rand() < exp(-deltaV/constants()->T())) {
        energy += deltaV;
        numMoveAccept++;
    }
    else {
        config(p) = oldPos;
    }
}

/**************************************************************************//**
 * Perform a simple insert move.
 *
******************************************************************************/
void ClassicalMonteCarlo::insertParticle() {

    dVec newPos;
    newPos = 0.0;

    numInsertTotal++;

    newPos = boxPtr->randPosition(random);

    /* Compute the old energy of particle p*/
    deltaV = externalPtr->V(newPos);
    for (int p2 = 0; p2 < numParticles; p2++) {
        sep = newPos-config(p2);
        boxPtr->putInside(sep);
        deltaV += interactionPtr->V(sep);
    }

    double factor = z*boxPtr->volume/(numParticles+1);

    /* Now the metropolis step */
    if (random.rand() < factor*exp(-deltaV/constants()->T())) {
        energy += deltaV;
        if (config.extents()[0] < (numParticles+1))
            config.resizeAndPreserve(numParticles+1);
        config(numParticles) = newPos;
        numParticles++;
        numInsertAccept++;
    }
}

/**************************************************************************//**
 * Perform a simple delete move.
 *
******************************************************************************/
void ClassicalMonteCarlo::deleteParticle() {

    numDeleteTotal++;

    int p = random.randInt(numParticles-1);

    /* Compute the old energy of particle p*/
    deltaV = -externalPtr->V(config(p));
    for (int p2 = 0; p2 < numParticles; p2++) {
        if (p != p2) {
            sep = config(p)-config(p2);
            boxPtr->putInside(sep);
            deltaV -= interactionPtr->V(sep);
        }
    }

    double factor = numParticles/(z*boxPtr->volume);

    /* Now the metropolis step */
    if (random.rand() < factor*exp(-deltaV/constants()->T())) {
        energy += deltaV;
        config(p) = config(numParticles-1);
        numParticles--;
        numDeleteAccept++;
    }
}

/**************************************************************************//**
 * Perform measurements
 *
******************************************************************************/
void ClassicalMonteCarlo::measure(int &numMeasure) {
    std::cout << aveEnergy/(numMeasure) << "\t" << aveNumParticles/(numMeasure) 
         << "\t" << (3.0/2.0)*constants()->T() + aveEoN/(numMeasure)
         << "\t" << aveNumParticles/(numMeasure*boxPtr->volume)
         << "\t" << 1.0*numMoveAccept/(1.0*numMoveTotal) 
         << "\t" << 1.0*numInsertAccept/(1.0*numInsertTotal)
         << "\t" << 1.0*numDeleteAccept/(1.0*numDeleteTotal) << std::endl;
    aveEnergy = 0.0;
    aveEoN = 0.0;
    aveNumParticles = 0.0;
    numMeasure = 0;
}
