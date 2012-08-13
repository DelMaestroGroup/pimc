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
        Array <dVec,1> &initialPos) :
    externalPtr(_externalPtr), 
    interactionPtr(_interactionPtr),
    random(_random),
    boxPtr(_boxPtr),
    config(initialPos)
{
    /* The number of particles */
    numParticles = config.extent(firstDim);

    /* Get the initial configuration*/
//    config.resize(numParticles)
//    config = externalPtr->initialConfig;

    /* Set the fugacity z*/
    z = exp(constants()->mu()/constants()->T())/pow(constants()->dBWavelength(),NDIM);

    /* Compute the initial total potential energy */
    energy = 0.0;
    dVec sep;
    sep = 0.0;
    for (int part1 = 0; part1 < numParticles; part1++) {
        energy += externalPtr->V(config(part1));

        for (int part2 = part1+1; part2 < numParticles; part2++) {
            sep = config(part1)-config(part2);
            boxPtr->putInside(sep);
            energy += interactionPtr->V(sep);
        }
    }

    /* Initialize acceptance tracking */
    numUpdateTotal = 0;
    numUpdateAccept = 0;

    /* Initialize measurements */
    aveEnergy = 0.0;
}

/**************************************************************************//**
 * Destructor. 
******************************************************************************/
ClassicalMonteCarlo::~ClassicalMonteCarlo ()
{
    // empty destructor
}

/**************************************************************************//**
 * Perform the monte carlo simulation.
 *
******************************************************************************/
void ClassicalMonteCarlo::run() {
    int numMCSteps = 1000000;
    for(int n = 1; n < numMCSteps; n++) {
        updateMove();
        insertMove();
        deleteMove();

        aveEnergy += energy;
        aveNumParticles += numParticles;

        if ((n % 50) == 0) {
            measure(n);
        }
    }
}

/**************************************************************************//**
 * Perform a simple positional update
 *
******************************************************************************/
void ClassicalMonteCarlo::updateMove() {

    dVec oldPos,sep;
    oldPos = 0.0;
    sep = 0.0;
    double oldE,newE,deltaE;

    for (int particle = 0; particle < numParticles; particle++) {

        numUpdateTotal++;

        int p = random.randInt(numParticles-1);

        /* Compute the old energy of particle p*/
        oldE = externalPtr->V(config(p));
        for (int p2 = 0; p2 < numParticles; p2++) {
            if (p != p2) {
                sep = config(p)-config(p2);
                boxPtr->putInside(sep);
                oldE += interactionPtr->V(sep);
            }
        }

        oldPos = config(p);

        /* The new random position */
        config(p) = boxPtr->randUpdate(random,oldPos);

        /* Compute the new energy of particle p*/
        newE = externalPtr->V(config(p));
        for (int p2 = 0; p2 < numParticles; p2++) {
            if (p != p2) {
                sep = config(p)-config(p2);
                boxPtr->putInside(sep);
                newE += interactionPtr->V(sep);
            }
        }

        deltaE = newE - oldE;

        /* Now the metropolis step */
        if (random.rand() < exp(-deltaE/constants()->T())) {
            energy += deltaE;
            numUpdateAccept++;
        }
        else {
            config(p) = oldPos;
        }
    }
}

/**************************************************************************//**
 * Perform a simple insert move.
 *
******************************************************************************/
void ClassicalMonteCarlo::insertMove() {

    dVec newPos,sep;
    newPos = 0.0;
    sep = 0.0;
    double newE;

    numInsertTotal++;

    newPos = boxPtr->randPosition(random);

    /* Compute the old energy of particle p*/
    newE = externalPtr->V(newPos);
    for (int p2 = 0; p2 < numParticles; p2++) {
        sep = newPos-config(p2);
        boxPtr->putInside(sep);
        newE += interactionPtr->V(sep);
    }

    double factor = z*boxPtr->volume/(numParticles+1);

    /* Now the metropolis step */
    if (random.rand() < factor*exp(-newE/constants()->T())) {
        energy += newE;
        if (config.extent(firstDim) < (numParticles+1))
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
void ClassicalMonteCarlo::deleteMove() {

    dVec oldPos,sep;
    oldPos = 0.0;
    sep = 0.0;
    double oldE;

    numDeleteTotal++;

    int p = random.randInt(numParticles-1);
    oldPos = config(p);

    /* Compute the old energy of particle p*/
    oldE = externalPtr->V(oldPos);
    for (int p2 = 0; p2 < numParticles; p2++) {
        if (p != p2) {
            sep = oldPos-config(p2);
            boxPtr->putInside(sep);
            oldE += interactionPtr->V(sep);
        }
    }

    double factor = numParticles/(z*boxPtr->volume);

    /* Now the metropolis step */
    if (random.rand() < factor*exp(oldE/constants()->T())) {
        energy -= oldE;
        config(p) = config(numParticles-1);
        numParticles--;
        numDeleteAccept++;
    }
}

/**************************************************************************//**
 * Perform measurements
 *
******************************************************************************/
void ClassicalMonteCarlo::measure(int n) {
    cout << aveEnergy/(50.0) << "\t" << aveNumParticles/(50.0) 
         << "\t" << 1.0*numUpdateAccept/(1.0*numUpdateTotal) 
         << "\t" << 1.0*numInsertAccept/(1.0*numInsertTotal)
         << "\t" << 2.0*numDeleteAccept/(1.0*numDeleteTotal) << endl;
    aveEnergy = 0.0;
    aveNumParticles = 0.0;
}
