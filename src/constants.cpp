/** 
 * @file constants.cpp
 * @author Adrian Del Maestro
 *
 * @brief ConstantParameters class implementation.
 */

#include "constants.h"
#include <time.h>

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// CONSTANT PARAMETERS  CLASS ------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/**************************************************************************//**
 *  An empty constructor which simply sets all constants to null.
******************************************************************************/
ConstantParameters::ConstantParameters() 
{ 
    /* empty constructor */
}

/**************************************************************************//**
 *  Initialize all constants from command line, XMl and defaults. 
 *
 *  We initialize all constant parameters used in the simulation.  The value of
 *  lambda = hbar^2/2 m k_B is computed in units where lenghts are measured in 
 *  angstroms and energies in kelvin.
******************************************************************************/
void ConstantParameters::initConstants(po::variables_map &params) {

    /* We use boost to generate a UUID for the simulation */
    if (params["restart"].empty()) {

        id_ = boost::uuids::to_string(boost::uuids::random_generator()());

        /* Add a possible user specified label */
        string label_ = params["label"].as<string>();
        if (label_.length() > 12)
            label_ = label_.substr(0,12);
        
        id_.replace(id_.end()-label_.length(),id_.end(),label_);
        restart_ = false;
    }
    else {
        id_ = params["restart"].as<string>();
        restart_ = true;
    }

    /* Are we starting from a supplied state file? */
    startWithState_ = !params["start_with_state"].as<string>().empty();
    
    /* Set the wall clock state */
    if (params["wall_clock"].empty()) {
        wallClockOn_ = false;
        wallClock_ = 0;
    } 
    else {
        wallClockOn_ = true;
        /* Set wallClock_ in seconds*/
        wallClock_ = uint32( floor(params["wall_clock"].as<double>()*3600));
    }

    /* Are we working in the grand canonical ensemble? */
    canonical_ = !params["canonical"].empty();

    /* Are we saving a state file every bin? */
    saveStateFiles_ = params["no_save_state"].empty();

    /* Do we want variable length diagonal updates? */
    varUpdates_ = params["var_updates"].empty();
    
    /* Set the particle number window */
    window_ = canonical_ && !params["window"].empty();
    if (window_)
        windowWidth_ = params["window"].as<int>();
    else
        windowWidth_ = 0;
    
    /* Set the ensemble weight */
    gaussianEnsemble_ = canonical_ && !params["gaussian_window_width"].empty();
    if (gaussianEnsemble_)
        gaussianEnsembleSD_ = params["gaussian_window_width"].as<double>();
    else
        gaussianEnsembleSD_ = 0.0;

    /* The maximum winding number sampled */
    maxWind_ = params["max_winding"].as<int>();

    /* Assigned values */
    b_  = int (ceil(log(1.0*params["update_length"].as<int>()) / log(2.0)-EPS));

    /* We need to make sure b_ < numTimeSlices */
    while (ipow(2,b_) >= params["number_time_slices"].as<int>())
        b_--;

    /* Assigned values */
    Mbar_           = params["update_length"].as<int>();
    T_              = params["temperature"].as<double>();
    imagTimeLength_ = params["imaginary_time_length"].as<double>();
    mu_             = params["chemical_potential"].as<double>();
    m_              = params["mass"].as<double>();
    lambda_         = 24.24 / m_;
    rc_             = params["potential_cutoff"].as<double>();
    rc2_            = rc_*rc_;
    C0_             = params["worm_constant"].as<double>();
    numTimeSlices_  = params["number_time_slices"].as<int>();
    aCC_            = params["carbon_carbon_dist"].as<double>();
    if (PIGS)
        tau_       = 1.0/((numTimeSlices_-1)*T_);
    else
        tau_       = 1.0/(numTimeSlices_*T_);
    V_             = params["volume"].as<double>();
    L_             = params["side"].as<dVec>()[NDIM-1];
    numEqSteps_    = params["number_eq_steps"].as<uint32>();
    binSize_       = params["bin_size"].as<uint32>();

    graphenelut3d_file_prefix_ = params["graphenelut3d_file_prefix"].as<string>();
    isf_input_                 = params["isf_input"].as<string>();
    isf_input_type_            = params["isf_input_type"].as<string>();
    virialWindow_              = params["virial_window"].as<int>();
    
    initialNumParticles_ = params["number_particles"].as<int>();
    numBroken_ = params["number_broken"].as<int>();

    spatialSubregionOn_ = !params["spatial_subregion"].empty();
    if (spatialSubregionOn_)
        spatialSubregion_ = params["spatial_subregion"].as<double>();

    endFactor_ = params["end_factor"].as<double>();
    Npaths_ = params["number_paths"].as<int>();

    intPotentialType_ = params["interaction"].as<string>();
    extPotentialType_ = params["external"].as<string>();
    waveFunctionType_ = params["wavefunction"].as<string>();
    actionType_       = params["action"].as<string>();

    /* Computed values */
    dBWavelength_ = 2.0*sqrt(M_PI * lambda_ / T_);
    comDelta_ = 0.04*dBWavelength_;
    displaceDelta_ = 0.04*dBWavelength_;
    getC();

    /* Set the move probabilities */

    /* At present, the pigs code has only diagonal moves */
    if (PIGS) {
        attemptProb_["open"] = 0.0;
        attemptProb_["insert"] = 0.0;
        attemptProb_["close"] = 0.0;
        attemptProb_["advance head"] = 0.0;
        attemptProb_["recede head"] = 0.0;
        attemptProb_["advance tail"] = 0.0;
        attemptProb_["recede tail"] = 0.0;
        attemptProb_["remove"] = 0.0;
        attemptProb_["swap head"] = 0.0;
        attemptProb_["swap tail"] = 0.0;
        attemptProb_["diagonal"] = 0.6;
        attemptProb_["center of mass"] = 0.1;
        attemptProb_["displace"] = 0.0;
        attemptProb_["end staging"] = 0.3;
        attemptProb_["mid staging"] = 0.0;
        attemptProb_["swap break"] = 0.0;
    }
    else {
        attemptProb_["open"] = 0.4;
        attemptProb_["insert"] = 0.4;
        attemptProb_["close"] = 0.15;
        attemptProb_["advance head"] = 0.075;
        attemptProb_["recede head"] = 0.075;
        attemptProb_["advance tail"] = 0.075;
        attemptProb_["recede tail"] = 0.075;
        attemptProb_["remove"] = 0.15;
        attemptProb_["swap head"] = 0.10;
        attemptProb_["swap tail"] = 0.10;
        attemptProb_["diagonal"] = 0.19;
        attemptProb_["center of mass"] = 0.01;
        attemptProb_["displace"] = 0.0;
        attemptProb_["end staging"] = 0.0;
        attemptProb_["swap break"] = 0.0;
        attemptProb_["mid staging"] = 0.0;
    }

    double totProb = attemptProb_["close"] + attemptProb_["advance head"] + attemptProb_["recede head"]
        + attemptProb_["advance tail"] + attemptProb_["recede tail"] + attemptProb_["remove"]
        + attemptProb_["swap head"] + attemptProb_["swap tail"] + attemptProb_["diagonal"] 
        + attemptProb_["center of mass"] + attemptProb_["displace"] + attemptProb_["end staging"] 
        + attemptProb_["mid staging"]+attemptProb_["swap break"];

    if (abs(totProb - 1.0) > EPS) {
        cout << "Close + AdvanceHead + RecedeHead + AdvanceTail + RecedeTail + Remove + SwapHead " 
             << "+ SwapTail + Diagonal + CoM Probability != 1" << endl;
        cout << totProb << endl;
        exit(EXIT_FAILURE);
    }
    PIMC_ASSERT(totProb-1.0 < EPS);

    totProb = attemptProb_["open"] + attemptProb_["insert"] + attemptProb_["diagonal"]
       + attemptProb_["center of mass"] + attemptProb_["displace"] + attemptProb_["swap break"] 
       + attemptProb_["end staging"] + attemptProb_["mid staging"];
    
    if (abs(totProb - 1.0) > EPS) {
        cout << "Open + Insert + Diagonal + CoM Probability != 1" << endl;
        cout << totProb << endl;
        exit(EXIT_FAILURE);
    }
    PIMC_ASSERT(totProb-1.0 < EPS);
}

/**************************************************************************//**
 *  This public method returns an instance of the constant object,  Only one
 *  can ever exist at a time.
******************************************************************************/
ConstantParameters* ConstantParameters::getInstance ()
{   
    static ConstantParameters inst;
    return &inst;
}
