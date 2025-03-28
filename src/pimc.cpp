/**
 * @file pimc.cpp
 * @author Adrian Del Maestro
 * 
 * @brief PathIntegralMonteCarlo class implementation.
 */

#include "pimc.h"
#include "path.h"
#include "lookuptable.h"
#include "move.h"

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// PATH INTEGRAL MONTE CARLO CLASS -------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/**************************************************************************//**
*  Constructor.
*
*  Here we initialize all data structures, moves and estimators that will be
*  required with peforming a path integral quantum monte carlo simulation.
*  The initialization depends on whether or not we are restarting, or starting
*  from a user supplied state.
******************************************************************************/
PathIntegralMonteCarlo::PathIntegralMonteCarlo (boost::ptr_vector<Path> &_pathPtrVec,
        MTRand &_random, boost::ptr_vector<move_vector> &_movePtrVec,
        boost::ptr_vector<estimator_vector> &_estimatorPtrVec, const bool _startWithState) :
    random(_random),
    binSize(constants()->binSize()),
    Npaths(_pathPtrVec.size()),
    pathPtrVec(_pathPtrVec),
    path(pathPtrVec.front()),
    movePtrVec(_movePtrVec),
    move(movePtrVec.front()),
    estimatorPtrVec(_estimatorPtrVec),
    estimator(estimatorPtrVec.front())
{
    /* Are we starting from a saved state? */
    startWithState = _startWithState;
    
    /* Initialize stateStrings */
    stateStrings.resize(Npaths);
    
    /* We keep simulating until we have stored a certain number of measurements */
    numStoredBins = 0;
    
    /* Initialize the number of sweeps and updates per sweep*/
    numSteps = 0;

    /* This is required in case we have zero particles in the simulation */
    numUpdates = int(ceil(std::max(constants()->initialNumParticles(),1)*constants()->numTimeSlices() / 
                (1.0*constants()->Mbar())));
    if (numUpdates < 100)
        numUpdates = 100;
    
    /* Calculate the number of sweeps to make sure we touch every bead */
    numImagTimeSweeps = int(ceil((1.0*constants()->numTimeSlices()/(1.0*constants()->Mbar()))));
    
    /* These counters are used in the equilibration process */
    numConfig = 0;
    numDiagonal = 0;
    numCoMAttempted = 200;
    prevNumCoMAttempted = 200;
    numCoMAccepted = 0;
    numDisplaceAttempted = 200;
    numDisplaceAccepted = 0;
    numNAttempted = 0;

    numStepsAttempted = 2000*numUpdates;
    numMuAttempted = 2000*numUpdates;

    relaxmuMessage = false;
    relaxC0Message = false;
    equilMessage = false;
    equilODMessage = false;

    /* The target diagonal fraction (hard coded) */
    targetDiagFrac = 0.70;

    /* Do we need to shift C0 up or down? */
    sgnDiagFrac = 0;
    shiftC0 = 0.25;

    /* Have we found the desired C0 value for the diagonal fraction */
    foundC0 = false;

    /* The target/initial number of particles */
    N0 = constants()->initialNumParticles();

    /* Number probability distribution used when relaxing chemical potential.
     * We assume the maximum possible number is 2000 particles */
    PN.resize(2000);
    PN.fill(0);
    foundmu = false;
    muFactor = 1.0;
    sgnAveN = 0;

    /* Used for optimization of μ search */
    bestPN = 0;
    bestmu = constants()->mu();
    bestDiffAveN = 5000;
    
    /* Intialize the config number to zero */
    configNumber = 0;
    
    /* Determine the cumulative attempt move probabilities, and the indices of
     * various diagonal moves */
    double cumDiagProb = 0.0;
    double cumOffDiagProb = 0.0;
    std::string moveName;

    for (auto movePtr = move.begin(); movePtr != move.end(); ++movePtr) {

        /* Get the namne of the move and check if it is the generic diagonal
         * move */
        moveName = movePtr->getName();
        if ((moveName == "bisection") || (moveName == "staging"))
            moveName = "diagonal";

        /* Accumulate the diagonal moves */
        if ( (movePtr->operateOnConfig == DIAGONAL) || movePtr->operateOnConfig == ANY ) {
            cumDiagProb += constants()->attemptProb(moveName);
            attemptDiagProb.push_back(cumDiagProb);
        }
        else 
            attemptDiagProb.push_back(cumDiagProb);

        /* Accumulate the off-diagonal moves */
        if ( (movePtr->operateOnConfig == OFFDIAGONAL) || movePtr->operateOnConfig == ANY) {
            cumOffDiagProb += constants()->attemptProb(moveName);
            attemptOffDiagProb.push_back(cumOffDiagProb);
        }
        else
            attemptOffDiagProb.push_back(cumOffDiagProb);

        /* Find the indices of moves in case we need them */
        moveIndex[moveName] = std::distance(move.begin(), movePtr);
    }
    
    /* Make sure the move cumulative probability arrays add up to 1 */
    attemptDiagProb.back() = 1.0 + EPS;
    attemptOffDiagProb.back() = 1.0 + EPS;
    
    /* If we are restarting, or loading a state from disk, do so now. */
    if (startWithState || constants()->restart())
        loadState();

    /* Setup all the estimators for measurement i/o */
    for (auto &&estPtr : estimatorPtrVec) 
        for (auto &est : estPtr)
            est.prepare();

    /* Make a list of estimator names for the 0th estimator */
    for (auto estimatorPtr = estimator.begin(); estimatorPtr != estimator.end(); ++estimatorPtr) 
        estimatorIndex[estimatorPtr->getName()] = std::distance(estimator.begin(), estimatorPtr);
}

/**************************************************************************//**
 *  Destructor.
******************************************************************************/
PathIntegralMonteCarlo::~PathIntegralMonteCarlo () {
}

/**************************************************************************//**
 *  Performs one of various possible Monte Carlo updats on the path.
******************************************************************************/
std::string PathIntegralMonteCarlo::update(const double x, const int sweep, const int pathIdx=0) { 

    success = false;
    std::string moveName = "NONE";
    int index;

    /* Determine the index of the move to be performed */
    if (pathPtrVec[pathIdx].worm.isConfigDiagonal)
        index = std::lower_bound(attemptDiagProb.begin(),attemptDiagProb.end(),x)
            - attemptDiagProb.begin();
    else 
        index = std::lower_bound(attemptOffDiagProb.begin(),attemptOffDiagProb.end(),x)
            - attemptOffDiagProb.begin();

    /* Perform the move */
    moveName = movePtrVec[pathIdx].at(index).getName();
    success = movePtrVec[pathIdx].at(index).attemptMove();

    return moveName;
}


/**************************************************************************//**
 *  Diagonal Equilibration
 * 
 *  The diagonal-only portion of the equilibration.
******************************************************************************/
void PathIntegralMonteCarlo::equilStepDiagonal() {

    for (int n = 0; n < constants()->initialNumParticles(); n++) {
            
        double x = random.rand();
            
        /* Here we do the diagonal pre-equilibration moves, and allow for 
         * optimization of simulation parameters */
        if (x < constants()->attemptProb("center of mass")) {

            /* index for the center of mass */
            int index = moveIndex["center of mass"];

            /* perform the CoM move */
            for (auto & cmove : movePtrVec)
                cmove.at(index).attemptMove();

            /* We check how many CoM moves we have tried.  Every 200 moves, we see if we need
             * to adjust comDelta, provided we are in the pre-equilibration diagonal state. */
            if ( (move.at(index).getNumAttempted() > 0) 
                    && (move.at(index).getNumAttempted() > prevNumCoMAttempted)
                    && (move.at(index).getNumAttempted() % numCoMAttempted == 0) 
                    && (constants()->comDelta() < 0.5*(*std::min_element(path.boxPtr->side.begin(), path.boxPtr->side.end()))) ) {

                numCoMAccepted  = move.at(index).getNumAccepted() - numCoMAccepted;
                double CoMRatio = 1.0*numCoMAccepted / numCoMAttempted;
                if (CoMRatio < 0.2)
                    constants()->shiftCoMDelta(-0.6);
                else if (CoMRatio < 0.3)
                    constants()->shiftCoMDelta(-0.4);
                else if (CoMRatio < 0.4)
                    constants()->shiftCoMDelta(-0.2);
                else if (CoMRatio > 0.6)
                    constants()->shiftCoMDelta(0.2);
                else if (CoMRatio > 0.7)
                    constants()->shiftCoMDelta(0.4);
                else if (CoMRatio > 0.8)
                    constants()->shiftCoMDelta(0.6);

                /* Reset the counters */
                numCoMAccepted = move.at(index).getNumAccepted();
                prevNumCoMAttempted = move.at(index).getNumAttempted();
            } // CoM Delta Shift

        } // Center of mass move
        /* Now try a displace move */
        else if (x < constants()->attemptProb("center of mass") + constants()->attemptProb("displace")) {

            int index = moveIndex["displace"];
            for (auto & cmove : movePtrVec)
                cmove.at(index).attemptMove();

            /* We check how many displace moves we have tried.  Every numDisplaceAttempted moves, we see if we need
             * to adjust delta, provided we are in the pre-equilibration diagonal state. */
            if ( (move.at(index).getNumAttempted() > 0) 
                    && (move.at(index).getNumAttempted() % numDisplaceAttempted == 0) ) {

                numDisplaceAccepted = move.at(index).getNumAccepted() - numDisplaceAccepted;
                double displaceRatio = 1.0*numDisplaceAccepted / numDisplaceAttempted;
                if (displaceRatio < 0.2)
                    constants()->shiftDisplaceDelta(-0.6);
                else if (displaceRatio < 0.3)
                    constants()->shiftDisplaceDelta(-0.4);
                else if (displaceRatio < 0.4)
                    constants()->shiftDisplaceDelta(-0.2);
                else if (displaceRatio > 0.6)
                    constants()->shiftDisplaceDelta(0.2);
                else if (displaceRatio > 0.7)
                    constants()->shiftDisplaceDelta(0.4);
                else if (displaceRatio > 0.8)
                    constants()->shiftDisplaceDelta(0.6);

                //cout << "delta: " << constants()->delta() << " " << displaceRatio << std::endl;
                /* Reset the counters */
                numDisplaceAccepted = move.at(index).getNumAccepted();
            }
        } // displace move
        else {
            /* Attemp a diagonal path update*/
            for (int sweep = 0; sweep < numImagTimeSweeps; sweep++){
                for (auto &cmove : movePtrVec)
                    cmove.at(moveIndex["diagonal"]).attemptMove();
            }
        }

    } // initNumParticles
}

/**************************************************************************//**
 *  Relax the chemical potential to find a target number of particles.
 * 
 *  This portion of the equilibration attempts to find the appropriate chemical
 *  potential to converge onto a target number of particles.
******************************************************************************/
bool PathIntegralMonteCarlo::equilStepRelaxmu() {

    /* Print a message when starting the relaxation */
    if (!relaxmuMessage) {
        relaxmuMessage = true;
        std::cout << format("[PIMCID: %s] - Relax Chemical Potential.") % constants()->id() << std::endl;
    }

    double shiftmu = 1.0;
    double mu;

    for (int n = 0; n < numUpdates; n++) {

        /* Generate random number and run through all moves */
        double x = random.rand();
        std::string mName;
        for(uint32 p=0; p<Npaths; p++) {
            mName = update(x,n,p);
        }

        /* Relax the chemical potential to get a desired average number of
         * particles */
        int cN = round(1.0*path.worm.getNumBeadsOn()/constants()->numTimeSlices());
        int curSize = PN.size();
        if (cN >= curSize) {
            PN.resizeAndPreserve(cN+1);
            for (int n = curSize; n < cN; n++)
                PN(n) = 0;
        }

        PN(cN)++;
        numNAttempted += 1;

        if ( (numNAttempted == numMuAttempted) && (N0 > 0) ) {

            /* Output the distribution to the terminal*/
            std::cout << printHistogram();

            /* Compute the peak loation and average number of particles */
            int peakN = static_cast<int>(std::distance(PN.begin(), std::max_element(PN.begin(), PN.end())));
            int aveN = round(weighted_average(PN));

            /* If we have shifted the peak to the desired value, exit */
            if (peakN == N0) {
                std::cout << format("Converged on μ = %8.5f\n\n") % constants()->mu();
                return true;
            }
            else {
                std::string method;

                /* We make sure the peak is in a window 10 particles away from
                 * the target. */
                bool peakInWindow = ((peakN >= N0-5) && (peakN <= N0+5));

                if ((PN(N0) > 0) &&  (PN(N0+1) > 0) && (PN(N0-1) > 0) && peakInWindow) {
                    method = "Equal";
                    mu = constants()->mu() - 0.5*constants()->T()*log(1.0*PN(N0+1)/PN(N0-1));
                }
                else if ((PN(N0) > 0) &&  (PN(N0+1) > 0) && peakInWindow) {
                    method = "Down";
                    mu = constants()->mu() - constants()->T()*log(1.0*PN(N0+1)/PN(N0));
                }
                else if ((PN(N0) > 0) &&  (PN(N0-1) > 0) && peakInWindow) {
                    method = "Up";
                    mu = constants()->mu() - constants()->T()*log(1.0*PN(N0)/PN(N0-1));
                }
                else {
                    if (!inWindow) {
                        std::cout << "Reseting μ to previous best value." << std::endl;
                        mu = bestmu;
                    }
                    else {
                        if (peakN > N0) {
                            mu = constants()->mu() - shiftmu*muFactor;
                            method = "Down";
                        }
                        else {
                            mu = constants()->mu() + shiftmu*muFactor;
                            method = "Up";
                        }

                        /* Determine if we need to alter the muFactor. */
                        int sgn = ((N0-aveN) > 0) - ((N0-aveN) < 0);
                        if (sgnAveN != 0) {
                            if (sgn == sgnAveN)
                                muFactor *= 1.618;
                            else
                                muFactor /= 1.61803399;
                        }
                        sgnAveN = sgn;
                    }
                }

                /* Some optional debugging messages */
                /* std::cout << format("Shifting %5s:%10.2f%10.2f%10d%10d%10d") % method % */ 
                /*     constants()->mu() % mu % PN(N0-1) % PN(N0) % PN(N0+1); */

                constants()->setmu(mu);

                /* std::cout << format("%12d%12d%12d\n") % peakN % aveN % path.getTrueNumParticles(); */

                numNAttempted = 0;
                PN.fill(0);

            } // haven't moved peak yet

        } // numNAttempted == numMuAttempted
    } // numupdates

    return false;
}

/**************************************************************************//**
 *  Relax the worm constant C0 until we have found the desired diagonal
 *  fraction.
 * 
 *  This portion of the equilibration attempts to find the appropriate value of
 *  the worm constant C0 until we have arrived at a desired diagonal fraction
 *  fixed to be 75% here.
******************************************************************************/
bool PathIntegralMonteCarlo::equilStepRelaxC0() {

    /* Print a message when starting the relaxation */
    if (!relaxC0Message) {
        relaxC0Message = true;
        std::cout << format("[PIMCID: %s] - Relax Worm Constant.\n") % constants()->id() << std::endl;
    }

    for (int n = 0; n < numUpdates; n++) {

        /* Generate random number and run through all moves */
        double x = random.rand();
        std::string mName;
        for(uint32 p=0; p<Npaths; p++) {
            mName = update(x,n,p);
        }

        /* We accumulate the number of diagonal configurations */
        numConfig++;
        if (path.worm.isConfigDiagonal) 
            ++numDiagonal;

        /* Every numStepsAttempted steps, we check the diagonal fraction and update the
         * worm constant. */
        if ( numConfig == numStepsAttempted) {

            double diagFrac = 1.0*numDiagonal / (1.0*numConfig);

            if ( (diagFrac > (targetDiagFrac-0.05)) && (diagFrac <= (targetDiagFrac+0.05)) ) {
                std::cout << format("\nConverged on C0 = %8.5f\n\n") % constants()->C0();
                /* for (int i = 0; i < diagFracVals.size(); i++) */ 
                /*     std::cout << format("%12.5f%12.5f\n") % C0Vals[i] % diagFracVals[i]; */
                return true;
            }
            else {
                std::cout << format("%4.2f\t%8.5f\t") % diagFrac % constants()->C0();

                /* Store the values of C0 and the diagonal fraciton */
                C0Vals.push_back(constants()->C0());
                diagFracVals.push_back(diagFrac);

                /* Perform an iterative linaer regression and suggest a new
                 * value */
                constants()->setC0(linearRegressionC0());

                std::cout << format("%8.5f\t%5d\t%8.6f\n") % constants()->C0() 
                    % path.getTrueNumParticles() 
                    % (1.0*path.getTrueNumParticles()/path.boxPtr->volume);

                /* Reset the counters */
                numDiagonal = 0;
                numConfig = 0;

            } // Haven't converged yet

        } // numConfig == numStepsAttempted

    } // numUpdates

    return false;
}

/**************************************************************************//**
 * Linear Regression 
******************************************************************************/
double PathIntegralMonteCarlo::linearRegressionC0() {

    auto num = diagFracVals.size();
    double diagFrac = diagFracVals.back();

    if (num == 1) {

        sgnDiagFrac = ((targetDiagFrac-diagFrac) > 0) - ((targetDiagFrac-diagFrac) < 0);

        if (diagFrac < targetDiagFrac) {
            /* double C0new = C0Vals.back() - shiftC0; */
            /* return (C0new > 0.0) ? C0new : 1.0E-4; */
            return 0.4*C0Vals.back();
        }
        else
            return 1.6*C0Vals.back();
            /* return C0Vals.back() + shiftC0; */
    }
    else {

        /* We first find the closest values to targetC0 */
        /* double closestC0; */
        /* double minDiff = 1000; */
        /* for (int i = 0; i < num; i++) { */
        /*     double cDiff = abs(diagFracVals[i]-targetC0); */
        /*     if (cDiff < minDiff) { */
        /*         closestC0 = C0Vals[i]; */
        /*         minDiff = cDiff; */
        /*     } */
        /* } */

        double sumX,sumY,sumX2,sumXY,sum1;

        sum1 = sumX = sumY = sumXY = sumX2 = 0.0;

        for (decltype(num) i = 0; i < num; i++) {
            sumX += C0Vals[i];
            sumX2 += C0Vals[i]*C0Vals[i];
            sumY += diagFracVals[i];
            sumXY += C0Vals[i]*diagFracVals[i];
            sum1 += 1.0;
        }

        double denom = sum1*sumX2 - sumX*sumX;
            double a0 = (sumY*sumX2 - sumX*sumXY)/denom;
        double a1 = (sum1*sumXY - sumY*sumX)/denom;

        double C0guess = (targetDiagFrac-a0)/a1;

        if ((C0guess < 1E4) && ( C0guess > 1E-4) && (a1 < 0.0))
            return C0guess;
        else
        {

            /* Determine if we need to alter the muFactor. */
            
            /* This is commented out for now, return to see if this can be
             * improved. */
            /* int sgn = ((targetDiagFrac-diagFrac) > 0) - ((targetDiagFrac-diagFrac) < 0); */
            /* if (sgn == sgnDiagFrac) */
            /*     shiftC0 *= 1.618; */
            /* else */
            /*     shiftC0 /= 1.61803399; */
            /* sgnDiagFrac = sgn; */

            if (diagFrac < targetDiagFrac) {
                /* double C0new = C0Vals.back() - shiftC0; */
                /* return (C0new > 0.0) ? C0new : 1.0E-4; */
                return (1.0/1.618)*C0Vals.back();
            }
            else
                return 1.61803399*C0Vals.back();
                /* return C0Vals.back() + shiftC0; */
        }
    }
}

/**************************************************************************//**
 *  Equilibration.
 * 
 *  The equilibration method, where we perform fully diagonal moves 20% of the
 *  time, finding an optimal, COM step, then if it is desired, we attempt to 
 *  find an optimal value of μ and C0.
******************************************************************************/
void PathIntegralMonteCarlo::equilStep(const uint32 iStep, const bool relaxC0, const bool relaxmu) {

    /* How far are we into the equilibration? */
    double equilFrac = (1.0*iStep) / (1.0*constants()->numEqSteps());

    /* We begin by attempting to perform a purely diagonal equilibration */
    if (equilFrac < 0.1 && !startWithState)  {
        if (iStep == 0) {
            numRemainingSteps = constants()->numEqSteps()/10;
            if (numRemainingSteps == 0)
                numRemainingSteps = 1;
            barStepSize = numRemainingSteps/44;
            if (!barStepSize)
                barStepSize = 1;
            std::cout << "[" << std::flush;
            numBars = 0;
        }
        equilStepDiagonal();

        /* Output a progress bar */
        if ((iStep > 0) && (iStep % barStepSize == 0) && (numBars < 44)) {
            std::cout << "▇" << std::flush;
            numBars++;
        }
        
        if (iStep == numRemainingSteps-1) {
            for (int nb = numBars; nb < 44; nb++)
                std::cout << "▇" << std::flush;
            std::cout << "] \t. Diagonal Pre-Equilibration." << std::endl << std::flush;
        }

    }
    else if ((equilFrac < 0.2) && (relaxmu || relaxC0)) {

        if (!equilODMessage) {
            equilODMessage = true;
            std::cout << "[";
            numRemainingSteps = int(0.2*constants()->numEqSteps())-iStep;
            barStepSize = numRemainingSteps/44;
            if (!barStepSize)
                barStepSize = 1;
            numBars = 0;
        }

        /* Output a progress bar */
        if ((iStep % barStepSize == 0) && (numBars < 44)) {
            std::cout << "▇" << std::flush;
            numBars++;
        }

        if (iStep == constants()->numEqSteps()/5-1) {
            for (int nb = numBars; nb < 44; nb++)
                std::cout << "▇" << std::flush;
            std::cout << "] \t. Off-Diagonal Pre-Equilibration." << std::endl << std::flush;
        }

        /* If we are performing any parameter relaxtion, we do some additional 
         * off-diagonal pre-equilibration */
        for (int n = 0; n < numUpdates; n++) {

            /* Generate random number and run through all moves */
            double x = random.rand();
            std::string mName;
            for(uint32 p=0;p<Npaths;p++){
                mName = update(x,n,p);
            }
        }
    }
    else {
        /* Perform possible parameter relaxation schemes */
        if (relaxmu && !foundmu)
            foundmu = equilStepRelaxmu();
        else if (relaxC0 && !foundC0) {
            foundC0 = equilStepRelaxC0();
        }
        /* Otherwise (or after converged) equilibrate */
        else {
            if (!equilMessage) {
                equilMessage = true;
                std::cout << format("[PIMCID: %s] - Equilibration Stage.") % constants()->id() << std::endl << "[";
                numRemainingSteps = constants()->numEqSteps()-iStep;
                barStepSize = numRemainingSteps/44;
                if (!barStepSize)
                    barStepSize = 1;
                numBars = 0;
            }

            /* Output a progress bar */
            if ((iStep % barStepSize == 0) && (numBars < 44)) {
                std::cout << "▇" << std::flush;
                numBars++;
            }

            /* This is the regular equilibration portion after we have the desired 
             * value of C0 and μ. */
            for (int n = 0; n < numUpdates; n++) {

                /* Generate random number and run through all moves */
                double x = random.rand();
                std::string mName;
                for(uint32 p=0;p<Npaths;p++){
                    mName = update(x,n,p);
                }
            } // numUpdates
        }

    } // equilfrac > 0.2

    /* Save a state every binsize equilibrium steps provided we are diagonal*/
    if ( path.worm.isConfigDiagonal && (iStep > 0) && (iStep % binSize) == 0) 
        saveState();

    if ((iStep == constants()->numEqSteps()-1) && equilMessage){
        for (int nb = numBars; nb < 44; nb++)
            std::cout << "▇" << std::flush;
        std::cout << "]" << std::endl;
    }
}

/**************************************************************************//**
 *  PIMC step.
 * 
 *  This method performs the metropolis sampling, which is actually a 
 *  complicated multi-step operation which consists of various types of moves.  
 *  They can in general occur at different frequencies. We also measure all 
 *  estimators.
******************************************************************************/
void PathIntegralMonteCarlo::step() {

    std::string moveName;

    /* perform updates on each set of paths */
    for (uint32 pIdx=0; pIdx<Npaths; pIdx++) {

        /* We run through all moves, making sure that we could have touched each bead at least once */
        for (int n = 0; n < numUpdates ; n++)  {
            moveName = update(random.rand(),n,pIdx);
        }

        /* Perform all measurements */
        for (auto& est : estimatorPtrVec[pIdx])
            est.sample();
        
        /* Every binSize measurements, we output averages to disk and record the
         * state of the simulation on disk.  */
        if (estimatorPtrVec[pIdx].size() > 0){
            if (estimatorPtrVec[pIdx].front().getNumAccumulated() >= binSize) {

                for (auto& est : estimatorPtrVec[pIdx]) {
                    /* std::cout << est.getNumAccumulated() << std::endl; */
                    if (est.getNumAccumulated() >= binSize)
                        est.output();
                }
                if(Npaths==1) 
                    saveState();
                if (pIdx == 0)
                    ++numStoredBins;
            }
        }
    }
        
    //////Multi-path Estimators////////////
    if(estimatorPtrVec.size() > Npaths) {

        /* Multi-Path estimators are at the end of the estimator std::vectors */
        for (auto& est : estimatorPtrVec.back()) 
            est.sample();

        /* Every binSize measurements, we output averages to disk and record the
         * state of the simulation on disk.  */
        if (estimatorPtrVec.back().front().getNumAccumulated() >= binSize) {

            for (auto& est : estimatorPtrVec.back()) 
                if (est.getNumAccumulated() >= binSize) 
                    est.output();

            /* Save to disk or store a state file */
            saveState();

            if (estimator.size() == 0)
                ++numStoredBins;
        }
    }
}

/**************************************************************************//**
 *  Output simulation statistics to disk.
 *
 *  We perform this once at the very end of the simulation.  It saves the 
 *  details of accetance probabilities as well as estimators to a file.
******************************************************************************/
void PathIntegralMonteCarlo::finalOutput() {

    /* Output the acceptance data to the log file */
    communicate()->file("log")->stream() << std::endl;
    communicate()->file("log")->stream() << std::endl;

    communicate()->file("log")->stream() << "---------- Begin Acceptance Data ---------------" << std::endl;
    communicate()->file("log")->stream() << std::endl;
    communicate()->file("log")->stream() << format("%-29s\t:\t%7.5f\n") % "Total Rate" 
        % move.front().getTotAcceptanceRatio();
    communicate()->file("log")->stream() << std::endl;

    /* Ouptut all the move acceptance information to disk */
    std::string moveName;
    for (auto &cmove : move){

        moveName = cmove.getName();

        /* We only output levels for those moves which have a variable size */
        if (cmove.variableLength) {
            for (int n = 0; n <= constants()->b(); n++) {
                communicate()->file("log")->stream() << format("%-12s Level %-10d\t:\t%7.5f\t(%d/%d)\n") 
                    % moveName % n % cmove.getAcceptanceRatioLevel(n) 
                    % cmove.numAcceptedLevel(n) % cmove.numAttemptedLevel(n);
            }
        }
        communicate()->file("log")->stream() << format("%-29s\t:\t%7.5f\t(%d/%d)\n") % moveName
            % cmove.getAcceptanceRatio() % cmove.numAccepted % cmove.numAttempted;
        communicate()->file("log")->stream() << std::endl;
    
    }
    communicate()->file("log")->stream() << "---------- End Acceptance Data -----------------" << std::endl;

    communicate()->file("log")->stream() << std::endl;
    communicate()->file("log")->stream() << std::endl;

    /* Output the estimator statistics to the log file */
    communicate()->file("log")->stream() << "---------- Begin Estimator Data ----------------" << std::endl;
    communicate()->file("log")->stream() << std::endl;
    for (auto &cestimator : estimator) {
        communicate()->file("log")->stream() << format("%-33s\t:\t%16d\t%16d\n") % cestimator.getName()
            % cestimator.getNumSampled() % cestimator.getTotNumAccumulated();
    
    }
    communicate()->file("log")->stream() << std::endl;
    communicate()->file("log")->stream() << "---------- End Estimator Data ------------------" << std::endl;
}

/**************************************************************************//**
*  Save the state of the simulation to disk, including all path, worm,
*  move and estimator data.
******************************************************************************/
void PathIntegralMonteCarlo::saveState(const int finalSave) {

    std::stringstream stateStrStrm;
    
    /* We only update the std::string during the simulation */
    if (!finalSave) {

        for(uint32 pIdx=0; pIdx<Npaths; pIdx++){

            stateStrStrm.str("");

            /* leftpack and update the lookup table arrays */
            pathPtrVec[pIdx].leftPack();
            pathPtrVec[pIdx].lookup.updateGrid(path);

            /* We First write the current total number of world lines */
            stateStrStrm << pathPtrVec[pIdx].getNumParticles() << std::endl;

            /* Now write the total acceptance information for all moves */
            stateStrStrm << format("%16d\t%16d\n")
                % movePtrVec[pIdx].front().totAccepted % movePtrVec[pIdx].front().totAttempted;

            /* Now record the individual move acceptance information,
             * first for the diagonal, then off-diagonal*/
            for (const auto &cmove : movePtrVec[pIdx])
                stateStrStrm << format("%16d\t%16d\n") % cmove.numAccepted % cmove.numAttempted;

            /* Output the estimator sampling information */
            for (const auto &cestimator : estimatorPtrVec[pIdx])
                stateStrStrm << format("%16d\t%16d\n") % cestimator.getTotNumAccumulated()
                    % cestimator.getNumSampled();

            /* Now we output the actual path and worldline data */
            stateStrStrm << std::setprecision(16) << pathPtrVec[pIdx].beads << std::endl;
            stateStrStrm << pathPtrVec[pIdx].nextLink << std::endl;
            stateStrStrm << pathPtrVec[pIdx].prevLink << std::endl;

            /* Output the worm data */
            stateStrStrm << pathPtrVec[pIdx].worm.beads << std::endl;

            /* Save the state of the random number generator */
            uint32 randomState[random.SAVE];
            random.save(randomState);
            for (int i = 0; i < random.SAVE; i++)
                stateStrStrm << randomState[i] << " ";
            stateStrStrm << std::endl;

            /* store the state std::string */
            stateStrings[pIdx] = stateStrStrm.str();

        } // for pidx
    }

    /* Save the state file to disk */
    std::string stateFileName;
    if (constants()->saveStateFiles() || finalSave) {

        for(uint32 pIdx=0; pIdx<Npaths; pIdx++) {

            stateFileName = "state";
            if(pIdx > 0)
                stateFileName += str(format("%d") % (pIdx+1));

            /* Prepare the state file for writing */
            communicate()->file(stateFileName.c_str())->reset();

            /* Write the stateString to disk */
            communicate()->file(stateFileName.c_str())->stream() << stateStrings[pIdx];

            /* Rename and copy the file. */
            communicate()->file(stateFileName.c_str())->rename();
        }
    }
}

/**************************************************************************//**
 *  Load a classical ground state from file.
******************************************************************************/
void PathIntegralMonteCarlo::loadClassicalState(DynamicArray <dVec,2> &tempBeads,
        DynamicArray <unsigned int, 2> &tempWormBeads, int numWorldLines) {

    /* We go through each active worldline and create a new classical
     * configuration */
    int ptcl;
    beadLocator beadIndex;
    
    for(uint32 pIdx=0; pIdx<Npaths; pIdx++){
        ptcl = 0;
        for (int n = 0; n < numWorldLines; n++) {
            beadIndex = {0, n};

            /* Check if the bead is on */
            if (tempWormBeads(beadIndex)) {
                
                /* Assign the classical configuration */
            fill_mdspan(pathPtrVec[pIdx].beads.slice<1>(ptcl), tempBeads(beadIndex));
            fill_mdspan(pathPtrVec[pIdx].worm.beads.slice<1>(ptcl), 1);
                ptcl++;
            }
        }
    }
}

/**************************************************************************//**
 *  Load a quantum ground state from file.
******************************************************************************/
void PathIntegralMonteCarlo::loadQuantumState(DynamicArray <dVec,2> &tempBeads, 
        DynamicArray <beadLocator,2> &tempNextLink, DynamicArray <beadLocator,2> &tempPrevLink,
        int numTimeSlices, int tempNumWorldLines) {

    /* Prevent double counting of worldlines */
    DynamicArray <bool, 1> doBead(tempNumWorldLines);

    beadLocator startBead,beadIndex;
    beadLocator newStartBead;
    int ptcl;
    int slice;
    
    for(uint32 pIdx=0; pIdx<Npaths; pIdx++) {

        newStartBead = {0, 0};
        ptcl = 0;
        slice = 0;
        doBead.fill(true);
        
        /* Now we iterate through each worldline exactly once */
        for (int n = 0; n < tempNumWorldLines; n++) {

            /* The initial bead to be moved */
            startBead = {0, n};

            /* We make sure we don't try to touch the same worldline twice */
            if (doBead(n)) {

                beadIndex = startBead;

                /* The world line length, we simply advance until we have looped back on 
                 * ourselves. */
                slice = 0;
                newStartBead = {slice, ptcl};
                do {
                    /* We turn off any zero-slice beads we have touched */
                    if (beadIndex[0]==0)
                        doBead(beadIndex[1]) = false;

                    pathPtrVec[pIdx].beads(slice % numTimeSlices,ptcl) = tempBeads(beadIndex);
                    pathPtrVec[pIdx].worm.beads(slice % numTimeSlices,ptcl) = 1;

                    beadIndex = tempNextLink(beadIndex);
                    ++slice;

                    /* Do a forward reconnection, provided we are not at the
                     * last bead */
                    if ( ((slice % numTimeSlices) == 0) && !all(beadIndex, startBead)) {
                        pathPtrVec[pIdx].nextLink(numTimeSlices - 1, ptcl) = {0, ptcl + 1};
                        pathPtrVec[pIdx].prevLink(0, ptcl + 1) = {numTimeSlices - 1, ptcl};
                        ++ptcl;
                    }
                } while (!all(beadIndex, startBead));

                /* Now we have to add the remaining beads and perform the final
                 * reconnection */
                for (int tslice = (slice % numTimeSlices); tslice < numTimeSlices; tslice++) {
                    pathPtrVec[pIdx].beads(tslice,ptcl) = tempBeads(beadIndex);
                    pathPtrVec[pIdx].worm.beads(tslice,ptcl) = 1;
                }
                pathPtrVec[pIdx].nextLink(numTimeSlices - 1, ptcl) = newStartBead;
                pathPtrVec[pIdx].prevLink(newStartBead) = {numTimeSlices - 1, ptcl};
                ++ptcl;

            } // doBead
        } // n
    } // pIdx

}

/**************************************************************************//**
 *  Load a previous state of the simulation from disk.
******************************************************************************/
void PathIntegralMonteCarlo::loadState() {

    std::string tempString;
    
    std::string fileInitStr = "init";
    
    for( uint32 pIdx=0; pIdx<Npaths; pIdx++){
        
        if(pIdx>0){
            fileInitStr = str(format("init%d") % (pIdx+1));
        }

        /* Reset the total acceptance information */
        movePtrVec[pIdx].front().resetTotAccept();

        /* Reset all the individual move acceptance information */
        for (auto &cmove : movePtrVec[pIdx])
            cmove.resetAccept();

        /* Reset estimator sampling information */
        for (auto &cestimator : estimatorPtrVec[pIdx])
            cestimator.restart(0,0);

        /* We first read the former total number of world lines */
        int numWorldLines;
        int numTimeSlices = pathPtrVec[pIdx].numTimeSlices;
        communicate()->file(fileInitStr)->stream() >> numWorldLines;

        /* Now we skip through the input file until we find the beads matrix.  This
         * is signalled by the appearance of an open bracket "(" */
        while (!communicate()->file(fileInitStr)->stream().eof()) {
            if (communicate()->file(fileInitStr)->stream().peek() != '(') 
                getline (communicate()->file(fileInitStr)->stream(), tempString);
            else
                break;
        }

        /* Now we resize all path data members and read them from the init state file */
        pathPtrVec[pIdx].beads.resize(numTimeSlices,numWorldLines);
        pathPtrVec[pIdx].nextLink.resize(numTimeSlices,numWorldLines);
        pathPtrVec[pIdx].prevLink.resize(numTimeSlices,numWorldLines);
        pathPtrVec[pIdx].worm.beads.resize(numTimeSlices,numWorldLines);

        /* A temporary container for the beads array */
    DynamicArray <dVec,2> tempBeads;

        /* Get the worldline configuration */
        communicate()->file(fileInitStr)->stream() >> tempBeads;

        /* The temporary number of time slices */
        int tempNumTimeSlices = tempBeads.extents()[0];

        if (tempNumTimeSlices == numTimeSlices) {
        
            /* Copy over the beads array */
            pathPtrVec[pIdx].beads = tempBeads;

            /* Get the link arrays */
            communicate()->file(fileInitStr)->stream() >> pathPtrVec[pIdx].nextLink;
            communicate()->file(fileInitStr)->stream() >> pathPtrVec[pIdx].prevLink;

            /* Repeat for the worm file */
            communicate()->file(fileInitStr)->stream() >> pathPtrVec[pIdx].worm.beads;

        } // locBeads.rows() == numTimeSlices
        else {

            /* Initialize the links */
            // Get extents of DynamicArray
            std::array<std::size_t, 2> extents = pathPtrVec[pIdx].prevLink.extents();
            std::size_t rows = extents[0];
            std::size_t cols = extents[1];
            
            // Get pointers to the underlying contiguous storage.
            auto* pprev = pathPtrVec[pIdx].prevLink.data();
            auto* pnext = pathPtrVec[pIdx].nextLink.data();
            
            for (std::size_t i = 0; i < rows; ++i) {
                for (std::size_t j = 0; j < cols; ++j) {
                    std::size_t index = i * cols + j;
                    // For prevLink, assign {i-1, j}
                    pprev[index] = { static_cast<int>(i) - 1, static_cast<int>(j) };
                    // For nextLink, assign {i+1, j}
                    pnext[index] = { static_cast<int>(i) + 1, static_cast<int>(j) };
                }
            }
        
            /* Here we implement the initial periodic boundary conditions in 
             * imaginary time */
            auto pLrow = pathPtrVec[pIdx].prevLink.slice<0>(0);
            for (std::size_t j = 0; j < pLrow.extent(0); ++j) {
                pLrow(j) = { numTimeSlices - 1, static_cast<int>(j) };
            }

            auto nLrow = pathPtrVec[pIdx].nextLink.slice<0>(numTimeSlices - 1);
            for (std::size_t j = 0; j < nLrow.extent(0); ++j) {
                nLrow(j) = { 0, static_cast<int>(j) };
            }

            /* Reset the worm.beads array */
            pathPtrVec[pIdx].worm.beads.fill(0);

            /* Temporary containers for the links and worm beads */
            DynamicArray <beadLocator,2> tempNextLink;
            DynamicArray <beadLocator,2> tempPrevLink;
            DynamicArray <unsigned int,2> tempWormBeads;

            /* Get the link arrays and worm file */
            communicate()->file(fileInitStr)->stream() >> tempNextLink;
            communicate()->file(fileInitStr)->stream() >> tempPrevLink;
            communicate()->file(fileInitStr)->stream() >> tempWormBeads;

            /* Load a classical (all time slice positions equal) from the input
             * file */
            loadClassicalState(tempBeads,tempWormBeads, numWorldLines);

            /* Load a quantum initial state from a file */
            //if (tempNumTimeSlices < numTimeSlices) {
            //    loadQuantumState(tempBeads,tempNextLink,tempPrevLink,
            //            numTimeSlices,int(sum(tempWormBeads)/tempNumTimeSlices));
            //}

            /* Now we make sure all empty beads are unlinked */
            beadLocator beadIndex;
            for (int slice = 0; slice < numTimeSlices; slice++) {
                beadIndex[0] = slice;
                for (int ptcl = 0; ptcl < numWorldLines; ptcl++) {
                    beadIndex[1] = ptcl;
                    if (!pathPtrVec[pIdx].worm.beads(beadIndex)) {
                        pathPtrVec[pIdx].nextLink(beadIndex) = {XXX, XXX};
                        pathPtrVec[pIdx].prevLink(beadIndex) = {XXX, XXX};
                    }
                }
            }

        } // locBeads.rows() != numTimeSlices

        /* Load the state of the random number generator, only if we are restarting 
         * the simulation */
        if (constants()->restart()) {
            uint32 randomState[random.SAVE];
            for (int i = 0; i < random.SAVE; i++) 
                communicate()->file(fileInitStr)->stream() >> randomState[i];
            random.load(randomState);
        }

        /* Reset the number of on beads */
        pathPtrVec[pIdx].worm.resetNumBeadsOn();

        /* Go through all beads, and make sure they fit inside the simulation cell.
         * At the same time, determine how many active beads there are per slice */
        beadLocator beadIndex;
        pathPtrVec[pIdx].numBeadsAtSlice.fill(0);
        for (beadIndex[0] = 0; beadIndex[0] < numTimeSlices; ++beadIndex[0]) {
            for (beadIndex[1] = 0; beadIndex[1] < numWorldLines; ++beadIndex[1]) {
                pathPtrVec[pIdx].boxPtr->putInside(path(beadIndex));
                if (pathPtrVec[pIdx].worm.beadOn(beadIndex))
                    ++pathPtrVec[pIdx].numBeadsAtSlice(beadIndex[0]);
            } // particles
        } // slices

        /* Reset the worm parameters */
        pathPtrVec[pIdx].worm.isConfigDiagonal = true;
        pathPtrVec[pIdx].worm.reset();

        /* Resize and update the lookup table arrays */
        pathPtrVec[pIdx].lookup.resizeList(numWorldLines);
        pathPtrVec[pIdx].lookup.updateGrid(path);

        /* Reset the broken/closed worldline std::vectors */
        pathPtrVec[pIdx].resetBrokenClosedVecs();
        
        /* Close the file */
        communicate()->file(fileInitStr)->close();
        
    }
}

/**************************************************************************//**
 *  Output the worldline configuration to disk using PDB , suitable
 *  for plotting using vmd. 
 *  
 *  We must post-process the final pdb file and split it up due to connectivity 
 *  changes. 
 *  @see For the PDB specification: 
 *  http://www.wwpdb.org/documentation/format32/v3.2.html
******************************************************************************/
void PathIntegralMonteCarlo::outputPDB() {

    numParticles = path.getNumParticles();
    int numTimeSlices = path.numTimeSlices;

    configNumber++;

    /* We go through all beads, and find the start and end bead for each
     * worldline, adding them to an array */
    DynamicArray <beadLocator,1> startBead,endBead;
    startBead.resize(numParticles);
    endBead.resize(numParticles);

    /* We sort the output by the number of beads in a worldline */
    DynamicArray <int,1> wlLength(numParticles);
    wlLength.fill(0);

    /* This is the index-beadNumber mapping array */
    DynamicArray <int,2> beadNum(numTimeSlices,numParticles);
    beadNum.fill(0);

    int numWorldLines = 0;

    /* Get the list of beads that are active in the simulation */
    DynamicArray <bool,2> doBead(numTimeSlices,numParticles);      
    doBead = castArray<bool>(path.worm.getBeads());

    /* We go through each particle/worldline */
    int nwl = 0;
    int beadNumber = 0;
    for (int n = 0; n < numParticles; n++) {

        /* The initial bead to be moved */
        startBead(nwl) = {0,n};

        /* We make sure we don't try to touch the same worldline twice */
        if (doBead(startBead(nwl))) {

            /* Check to see if the start bead is on a worm.  If it is, we start
             * at the worm tail and end at its head. */
            if (path.worm.foundBead(path,startBead(nwl))) {
                startBead(nwl) = path.worm.tail;
                endBead(nwl)   = path.next(path.worm.head);
            }
            /* Otherwise, we loop around until we find the initial bead */
            else
                endBead(nwl) = startBead(nwl);

            /* Mark the beads as touched and increment the number of worldlines */
            beadLocator beadIndex;
            beadIndex = startBead(nwl);
            int length = 1;
            do {
                doBead(beadIndex) = false;
                beadNum(beadIndex) = beadNumber;
                beadNumber++;
                length++;
                beadIndex = path.next(beadIndex);
            } while (!all(beadIndex, endBead(nwl)));

            /* We label each trajectory by the number of particles it contains.
             * a worm is always given label 0 */
            if ((length % numTimeSlices) == 0)
                wlLength(nwl) = length/numTimeSlices;
            else 
                wlLength(nwl) = 0;

            nwl++;
        } // doBead

    } // n
    numWorldLines = nwl;

    /* Output the PDB header */
    communicate()->file("wl")->stream() << format("REMARK [CONFIG %04d]\n") % configNumber;

    /* Output the unit cell information.  It is always cubic.  Everything is scaled by
     * an overall factor for better visualization. */

    double scale = 10.0;                
    int i;
    communicate()->file("wl")->stream() << format("%-6s") % "CRYST1";
    for (i = 0; i < NDIM; i++) 
        communicate()->file("wl")->stream() << format("%9.3f") % (scale*path.boxPtr->side[i]);
    while (i < 3) {
        communicate()->file("wl")->stream() << format("%9.3f") % 1.0;
        i++;
    }
    communicate()->file("wl")->stream() << format("%7.2f%7.2f%7.2f %-11s%4d\n") % 90.0 % 90.0 % 90.0 % "P 1" % 1;

    /* We output the atom block */
    beadLocator beadIndex;
    for (int n = 0; n < numWorldLines;  n++) {
        beadIndex = startBead(n);
        do {
            /* We give the zero-time-slice bead a special name */
            if (beadIndex[0] == 0) {
                communicate()->file("wl")->stream() << format("%-6s%5d %-4s %03d %9s") % "ATOM" 
                    % beadNum(beadIndex) % "H0" % wlLength(n) % " ";
            }
            else {
                communicate()->file("wl")->stream() << format("%-6s%5d %-4s %03d %9s") % "ATOM" 
                    % beadNum(beadIndex) % "HE" % wlLength(n) % " ";
            }

            /* Output the coordinates in 3D */
            for (i = 0; i < NDIM; i++) {
                communicate()->file("wl")->stream() << format("%8.3f") % (scale*path(beadIndex)[i]);
            }
            while (i < 3) {
                communicate()->file("wl")->stream() << format("%8.3f") % 0.0;
                i++;
            }
            communicate()->file("wl")->stream() << format("%14s\n") % "HE";

            beadIndex = path.next(beadIndex);
        } while (!all(beadIndex, endBead(n)));
    }
    communicate()->file("wl")->stream() <<("TER\n");

    /* Now output the connect block */
    for (int n = 0; n < numWorldLines;  n++) {
        beadIndex = startBead(n);
        do {
            communicate()->file("wl")->stream() << format("%-6s%5d") % "CONECT" % beadNum(beadIndex);
            beadLocator prevIndex,nextIndex;
            prevIndex = path.prev(beadIndex);
            nextIndex = path.next(beadIndex);
            
            /* We make sure that we don't connect beads linked by periodic bondary 
             * conditions */

            /* Check the previous bead */
            if (path.worm.beadOn(prevIndex)) {
                dVec sep;
                sep = path(beadIndex) - path(prevIndex);
                if (dot(sep,sep) < path.boxPtr->rcut2)
                    communicate()->file("wl")->stream() << format("%5d") % beadNum(prevIndex);
            }

            /* Now the next bead */
            if (path.worm.beadOn(nextIndex)) {
                dVec sep;
                sep = path(nextIndex) - path(beadIndex);
                if (dot(sep,sep) < path.boxPtr->rcut2)
                    communicate()->file("wl")->stream() << format("%5d") % beadNum(nextIndex);
            }
            communicate()->file("wl")->stream() << std::endl;

            beadIndex = path.next(beadIndex);
        } while (!all(beadIndex, endBead(n)));
    }
    communicate()->file("wl")->stream() <<("END\n");
}

/**************************************************************************//**
 *  
 *  
******************************************************************************/
void PathIntegralMonteCarlo::printWormState() {

    /* We make a list of all the beads contained in the worm */
    DynamicArray <beadLocator,1> wormBeads;    // Used for debugging
    wormBeads.resize(path.worm.length+1);
    wormBeads.fill({XXX, XXX});

    /* Output the worldline configuration */
    communicate()->file("debug")->stream() << " (" << path.getTrueNumParticles() << ")" << std::endl;
    communicate()->file("debug")->stream() << "head " << path.worm.head[0] << " " << path.worm.head[1]
        << " tail " << path.worm.tail[0] << " " << path.worm.tail[1]
        << " length " << path.worm.length 
        << " gap " << path.worm.gap << std::endl;

    if (!path.worm.isConfigDiagonal) {
        beadLocator beadIndex;
        beadIndex = path.worm.tail;
        int n = 0;
        do {
            wormBeads(n) = beadIndex;
            beadIndex = path.next(beadIndex);
            ++n;
        } while(!all(beadIndex, path.next(path.worm.head)));
    }

    path.printWormConfig(wormBeads);
    path.printLinks<std::fstream>(communicate()->file("debug")->stream());
}

/**************************************************************************//**
 *  Output a histogram to the terminal.  Useful for visualizing distributions.
 *  
******************************************************************************/
std::string PathIntegralMonteCarlo::printHistogram() {

    std::string histogramDisplay = "\n";

    /* Setup the window of the histogram to be analyzed.*/
    decltype(PN.size()) start,end,window;
    window = 5;

    start = N0-window;
    if (start < 0) 
        start = 0;
    end = N0 + window;

    /* We make sure we haven't moved past the end of the probability
     * distribution array */
    if (end > PN.size()-1)
        end = PN.size()-1;

    double factor = 50.0/(*std::max_element(PN.begin(),PN.end()));
    int sumPN = 0;
    for (decltype(start) n = start; n <= end; n++) {
        int numStars = int(PN(n)*factor);
        histogramDisplay +=  str(format("%03d|") % n);
        for (int i = 0; i < numStars; i++)
            histogramDisplay +=  "▇";
        histogramDisplay += "\n";
        sumPN += PN(n);
    }
    histogramDisplay += "\n";

    /* if (sumPN > bestPN) { */
    /*     bestPN = sumPN; */
        /* bestmu = constants()->mu(); */
        /* std::cout << "setting mu" << std::endl; */
    /* } */

    double aveN = (weighted_average(PN));

    double diffAveN = abs(N0-aveN);

    if (diffAveN < bestDiffAveN) {
        bestDiffAveN = diffAveN;
        bestmu = constants()->mu();
    }

    inWindow = (sumPN > 0);
    if (bestPN == 0) 
        inWindow = true;

    return histogramDisplay;
}
