/** 
 * @file setup.cpp
 * @author Adrian Del Maestro
 * @date 04.06.2011
 * 
 * @brief Implements the setup stage of the simulation.
 */

#include "setup.h"
#include "container.h"
#include "constants.h"
#include "communicator.h"
#include "potential.h"
#include "wavefunction.h"
#include "action.h"
#include "move.h"
#include "estimator.h"

/**************************************************************************//**
 * Create a comma separated list from a vector of strings
 *
 * @param option the stl vector of options
 * @return a comma separated list of options
******************************************************************************/
string getList(const vector<string> &options) {

    ostringstream optionList;
    std::copy(options.begin(),options.end()-1, 
            std::ostream_iterator<string>(optionList, ", "));
    optionList << options.back();
    return optionList.str();
}

/**************************************************************************//**
 * Overload << to print a vector of strings to the terminal 
******************************************************************************/
ostream& operator<<(ostream& os, const vector<string>& vals)
{
    for (auto const& val :vals) {
        if (val != vals.back())
            os << val << ", ";
        else
            os << val;
    } 
    return os;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// PARAMETER CLASS -----------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/**************************************************************************//**
 * Split a string at a supplied delimeter.
 *
 * @param s the string to be split
 * @param delim the string delimiter
 * @return a vector of strings
******************************************************************************/
vector<string> Parameters::split(const string &s, char delim) {
    vector<string> elems;
    stringstream ss(s);
    string item;
    while (getline(ss, item, delim))
        elems.push_back(item);
    return elems;
}

/**************************************************************************//**
 * Print all values in a parameter map.
 *
 * @see http://stackoverflow.com/questions/21008893/boost-program-options-iterate-over-variables-map
******************************************************************************/
void Parameters::print() {

    for(auto const& par : params ) {
        if (!par.second.empty()) {
            string key = par.first;
            cout << key << ":\t";
            if (type.at(key) == typeid(bool))
                cout << params.count(key);
            else
                extract[par.first](par.second);
            cout << endl;
        }
    }
}

/**************************************************************************//**
 * Insert all parameters into a pointer map of options.
 *
 * At present we can accept bool, double, int, uint32, string and
 * vector<string> options.  This needs to be modified for each new option data 
 * type.
******************************************************************************/
void Parameters::setupCommandLine(boost::ptr_map<string,po::options_description> &_options)  {

    /* We iterate through the parameters map */
    for (auto & par : params) {

        /* get the key */
        string key = par.first;

        /* Get the correct option class */
        po::options_description &option = _options[pClass[key]];

        /* construct a command line label out of the long and short name */
        string label = key;
        if (shortName[key] != "")
            label += "," + shortName[key];

        /* insert the option for each possible data type */
        if (type.at(key) == typeid(bool))
            option.add_options()(label.c_str(),helpMessage[key].c_str());
            /* option.add_options()(label.c_str(),po::bool_switch()->default_value(false),helpMessage[key].c_str()); */
        else if (type.at(key) == typeid(double))
            option.add_options()(label.c_str(),initValue<double>(key),helpMessage[key].c_str());
        else if (type.at(key) == typeid(int))
            option.add_options()(label.c_str(),initValue<int>(key),helpMessage[key].c_str());
        else if (type.at(key) == typeid(uint32))
            option.add_options()(label.c_str(),initValue<uint32>(key),helpMessage[key].c_str());
        else if (type.at(key) == typeid(string))
            option.add_options()(label.c_str(),initValue<string>(key),helpMessage[key].c_str());
        else if (type.at(key) == typeid(vector<string>)) {

            if (state[key] == DEFAULTED) {
                /* a temporary copy of the vector of strings */
                vector<string> ops = par.second.as<vector<string>>();
                option.add_options()(label.c_str(), po::value<vector<string>>()->
                        default_value(ops, getList(ops))->composing(), helpMessage[key].c_str());
            }
            else
                option.add_options()(label.c_str(), 
                        po::value<vector<string>>()->composing(), helpMessage[key].c_str());
        }
        else
            cerr << "insertOption Failed to find a valid type.";
    }
}

/**************************************************************************//**
 * Update the state of all parameters 
 *
 * We check both the commmand line options as well as any parameters specified
 * in an XML file.  The command line takes precedence.
 *
 * !!NB!! XML options needed to be correctly sorted in the file or they will
 * be ignored.  This should be fixed to exit gracefully with an error message.
******************************************************************************/
void Parameters::update(int argc, char *argv[], po::options_description &cmdLineOptions) {

    /* create a temporary parameter map */
    po::variables_map cmdparams;                    

    /* Get the values from the command line */
    po::store(po::parse_command_line(argc, argv, cmdLineOptions), cmdparams);
    po::notify(cmdparams);
    
    /* Now we go through our local parameter map and adjust values, keeping in
     * mind that we store boolean values */
    for (auto & par : params) {

        /* get the key and value */
        string key = par.first;
        auto &val = par.second;

        /* If the parameter exists, copy it over */
        if (cmdparams.count(key)) {

            val = cmdparams[key];

            /* update the state */
            if (cmdparams[key].empty())
                state[key] = UNSET; 
            else if (cmdparams[key].defaulted())
                state[key] = DEFAULTED;
            else
                state[key] = SET;
            }
    }

    /* Now we load a potential xml file from disk. It is optional so we use a
     * try/catch. */
    if (!params["param_file"].empty()) {

        try {
            pt::ptree xmlParams;                        
            pt::read_xml(params["param_file"].as<string>(),xmlParams);
            xmlParams = xmlParams.get_child("pimc");

            /* Go through our local parameter map and determine what needs to be
             * set from the xml file */
            for (auto & par : params) {

                /* get the key */
                string key = par.first;

                /* the parameter class string */
                string paramClass = pClass[key] + "_parameters";

                /* Get a clean xml node */
                pt::ptree xml;
                if (getNode(xmlParams,xml,paramClass)) {

                    if (type.at(key) == typeid(bool)){
                        set<bool>(key,xml);
                    }
                    else if (type.at(key) == typeid(double)) {
                        set<double>(key,xml);
                    }
                    else if (type.at(key) == typeid(int)) {
                        set<int>(key,xml);
                    }
                    else if (type.at(key) == typeid(uint32)) {
                        set<uint32>(key,xml);
                    }
                    else if (type.at(key) == typeid(string)) {
                        set<string>(key,xml);
                    }
                    else if (type.at(key) == typeid(vector<string>)) {

                        string pluralKey = key + "s";
                        if (xml.count(pluralKey) != 0) {

                            /* a temporary copy of the vector of strings */
                            vector<string> ops;

                            /* get the options */
                            for (auto const &opName : xml.get_child(pluralKey)) 
                                ops.push_back(opName.second.data());

                            /* update the parameter */
                            set<vector<string>>(key,ops);
                            state[key] = SET;
                        }
                    }
                    else
                        cerr << "xmlInsertOption Failed to find a valid type.";
                } //if getNode
            } //for par

        }
        catch(...) {
            cerr << "Cannot read parameter file: " << params["param_file"].as<string>() << endl;
        }
    }
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// SETUP CLASS ---------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
//
/**************************************************************************//**
 * Setup the program_options variables.
 *
 * We initialize all variables and define the names of all allowed interaction 
 * and external potentials.
******************************************************************************/
Setup::Setup() :
    params(),
    cmdLineOptions("Command Line Options")
{
    /* Initialize the option class names */
    optionClassNames = {"simulation","physical","cell","algorithm","potential",
        "measurement"};

    /* Define the allowed interaction potential names */
    interactionPotentialName = {"aziz", "szalewicz", "delta", "lorentzian", "sutherland", 
        "hard_sphere", "hard_rod", "free", "delta1D", "harmonic", "dipole"};
    interactionNames = getList(interactionPotentialName);

    /* Define the allowed external  potential names */
    externalPotentialName = {"free", "harmonic", "osc_tube", "lj_tube", "plated_lj_tube",
        "hard_tube", "hg_tube", "fixed_aziz", "gasp_prim", "fixed_lj", "graphene", "graphenelut",
         "graphenelut3d", "graphenelut3dgenerate", "graphenelut3dtobinary", "graphenelut3dtotext"};
    externalNames = getList(externalPotentialName);

    /* Define the allowed action names */
    actionName = {"primitive", "li_broughton", "gsf", "pair_product"};
    actionNames = getList(actionName);

    /* Define the allowed trial wave function names */
    waveFunctionName = {"constant", "sech", "jastrow", "lieb", "sutherland"};
    waveFunctionNames = getList(waveFunctionName);

    /* Define the allowed random number generator names */
    randomGeneratorName = {"boost_mt19937","std_mt19937", "pimc_mt19937"};
    randomGeneratorNames = getList(randomGeneratorName);

    /* Get the allowed estimator names */
    estimatorName = estimatorFactory()->getNames();
    vector<string> multiEstimatorName = multiEstimatorFactory()->getNames();
    estimatorName.insert(estimatorName.end(), multiEstimatorName.begin(), 
            multiEstimatorName.end());
    estimatorNames = getList(estimatorName);

    /* Get the allowed move names */
    moveName = moveFactory()->getNames(); 
    moveNames = getList(moveName);
}

/**************************************************************************//**
 * Setup all possible command line and xml parameters.
 *
******************************************************************************/
void Setup::initParameters() {

    /* Instantiate all the option classes */
    for(string oClass : optionClassNames) {
        string optionClassName = oClass;
        optionClassName[0] = toupper(optionClassName[0]);
        optionClassName += " Options";
        optionClasses.insert(oClass,new po::options_description(optionClassName));
    }

    /* Initialize the simulation options */
    string oClass = "simulation";
    params.add<bool>("help,h","produce help message",oClass);
    params.add<bool>("version","output repo version",oClass);
    params.add<bool>("validate","validate command line or xml options",oClass);
    params.add<bool>("dimension","output currently compiled dimension",oClass);
    params.add<int>("output_config,o","number of output configurations",oClass,0);
    params.add<uint32>("process,p","process or cpu number",oClass,0);
    params.add<string>("restart,R","restart running simulation with PIMCID",oClass);
    params.add<double>("wall_clock,W","set wall clock limit in hours",oClass);
    params.add<string>("start_with_state,s", "start simulation with a supplied state file.",oClass,"");
    params.add<bool>("no_save_state","Only save a state file at the end of a simulation",oClass);
    params.add<bool>("estimator_list","Output a list of estimators in xml format.",oClass);
    params.add<bool>("update_list","Output a list of updates in xml format.",oClass);
    params.add<string>("label","a label to append to all estimator files.",oClass,"");
    params.add<string>("rng,G",str(format("random number generator type:\n%s") % randomGeneratorNames).c_str(),oClass,"pimc_mt19937");
    params.add<string>("param_file","a valid path to the parameters input xml file.",oClass);

    /* Initialize the cell options */
    oClass = "cell";
    params.add<string>("geometry,b","simulation cell type [prism,cylinder]",oClass,"prism");
    params.add<double>("size,L","linear system size [angstroms]",oClass);
    params.add<double>("Lx","system size in x-direction [angstroms]",oClass);
    params.add<double>("Ly","system size in y-direction [angstroms]",oClass);
    params.add<double>("Lz","system size in z-direction [angstroms]",oClass);
    params.add<double>("radius,r","tube radius [angstroms]",oClass);

    /* Initialize the potential options */
    oClass = "potential";
    params.add<string>("interaction,I",str(format("interaction potential type:\n%s") % interactionNames).c_str(),oClass,"aziz");
    params.add<string>("external,X",str(format("external potential type:\n%s") % externalNames).c_str(),oClass,"free");
    params.add<double>("scattering_length,a","scattering length [angstroms]",oClass,1.0);
    params.add<double>("delta_width","delta function potential width",oClass,1.0E-3);
    params.add<double>("delta_strength,c","delta function potential integrated strength",oClass,10.0);
    params.add<double>("interaction_strength,g","interaction parameter",oClass,1.0);
    params.add<double>("omega","harmonic interaction potential frequency",oClass,1.0);
    params.add<double>("lj_sigma","Lennard-Jones hard-core radius [angstroms]",oClass,2.74);
    params.add<double>("lj_epsilon","Lennard-Jones energy scale [kelvin]",oClass,16.2463);
    params.add<double>("lj_width","Radial with of LJ plated cylinder material [angstroms]",oClass);
    params.add<double>("lj_density","Density LJ plated cylinder material [angstroms^(-3)]",oClass);
    params.add<double>("hourglass_radius","differential radius for hourglass potential [angstroms]",oClass,0.0);
    params.add<double>("hourglass_width","constriction width for hourglass potential [angstroms]",oClass,0.0);
    params.add<string>("fixed,f","input file name for fixed atomic positions.",oClass,"");
    params.add<double>("potential_cutoff,l","interaction potential cutoff length [angstroms]",oClass);
    params.add<double>("empty_width_y,y","how much space (in y-) around Gasparini barrier",oClass);
    params.add<double>("empty_width_z,z","how much space (in z-) around Gasparini barrier",oClass);

    /* These are graphene potential options */
    params.add<double>("strain","strain of graphene lattice in y-direction",oClass,0.00);
    params.add<double>("poisson","Poisson's ratio for graphene",oClass,0.165);
    params.add<double>("carbon_carbon_dist,A","Carbon-Carbon distance for graphene",oClass,1.42);
    params.add<string>("graphenelut3d_file_prefix","GrapheneLUT3D file prefix <prefix>serialized.{dat|txt}",oClass,"");

    /* Initialize the physical options */
    oClass = "physical";
    params.add<bool>("canonical","perform a canonical simulation",oClass);
    params.add<double>("mass,m","particle mass [amu]",oClass,4.0030);
    params.add<double>("density,n",str(format("initial density [angstroms^(-%d)]") % NDIM).c_str(),oClass);
    params.add<int>("number_particles,N","number of particles",oClass);
    params.add<double>("temperature,T","temperature [kelvin]",oClass);
    params.add<string>("wavefunction",str(format("trial wave function type:\n%s") 
                % waveFunctionNames).c_str(),oClass,"constant");
    params.add<double>("R_LL_wfn","length scale of the lieb liniger wave function",oClass,0.0);
    params.add<double>("k_LL_wfn","wave number of the lieb liniger wave function",oClass,0.0);
    params.add<double>("end_factor","end bead potential action multiplicatave factor",oClass,1.0);
    params.add<double>("chemical_potential,u","chemical potential [kelvin]",oClass,0.0);
    params.add<int>("number_paths","number of paths",oClass,1);

    /* Initialize the algorithm options */
    oClass = "algorithm";
    params.add<bool>("relax","perform a worm constant relaxation",oClass);
    params.add<bool>("relaxmu", "perform a chemical potential relaxation to target a fixed density",oClass);
    params.add<int>("number_time_slices,P","number of time slices",oClass);
    params.add<int>("window","set particle number window",oClass);
    params.add<double>("gaussian_window_width", "set gaussian ensemble weight",oClass);
    params.add<double>("imaginary_time_step,t", "imaginary time step [kelvin^(-1)]",oClass);
    params.add<double>("imaginary_time_length","total path length in imaginary time [kelvin^(-1)]",oClass);
    params.add<double>("worm_constant,C", "worm acceptance constant",oClass,1.0);
    params.add<double>("com_delta,D", "center of mass update radius[angstroms]",oClass);
    params.add<double>("displace_delta,d", "displace update radius [angstroms]",oClass);
    params.add<int>("update_length,M", "non-local update length, (Mbar), -1 sets maximal value",oClass);
    params.add<string>("action", str(format("action type:\n%s") % actionNames).c_str(),oClass,"gsf");
    params.add<bool>("full_updates", "perform full staging updates",oClass);
    params.add<bool>("var_updates", "perform variable length diagonal updates",oClass);
    params.add<bool>("staging", "perform staging instead of bisection for the diagonal update",oClass);
    params.add<int>("max_winding", "the maximum winding number to sample",oClass,1);

    /* Initialize the measurement options */
    oClass = "measurement";
    params.add<uint32>("number_eq_steps,E", "number of equilibration steps",oClass,1);
    params.add<int>("number_bins_stored,S", "number of estimator bins stored",oClass,1);
    params.add<uint32>("bin_size", "number of updates per bin",oClass,100);
    params.add<double>("estimator_radius,w", "maximum radius for cylinder estimators",oClass,2.0); 
    params.add<int>("virial_window,V", "centroid virial energy estimator window",oClass,5);
    params.add<int>("number_broken", "number of broken world-lines",oClass,0);
    params.add<double>("spatial_subregion", "define a spatial subregion",oClass);

    /* The updates, measurement defaults, and ensemble can depend on PIGS vs PIMC */
    vector<string> estimatorsToMeasure;
    vector<string> movesToPerform;
    if (PIGS) {
        params.set<bool>("canonical",true);
        estimatorsToMeasure = {EnergyEstimator::name, TimeEstimator::name};
        movesToPerform = {CenterOfMassMove::name, StagingMove::name, EndStagingMove::name,
            DisplaceMove::name};         
    }
    else {
        estimatorsToMeasure = {EnergyEstimator::name, NumberParticlesEstimator::name,
            TimeEstimator::name, DiagonalFractionEstimator::name};

        movesToPerform = {CenterOfMassMove::name, BisectionMove::name, OpenMove::name,
            CloseMove::name, InsertMove::name, RemoveMove::name, AdvanceHeadMove::name, 
            RecedeHeadMove::name, AdvanceTailMove::name, RecedeTailMove::name, SwapHeadMove::name,
            SwapTailMove::name};

    }
    params.add<vector<string>>("update", str(format("Monte Carlo updates to be performed:\n%s") 
                % moveNames).c_str(), "algorithm",movesToPerform);
    params.add<vector<string>>("estimator,e", str(format("estimators to be measured:\n%s") 
                % estimatorNames).c_str(),"measurement",estimatorsToMeasure);

}

/**************************************************************************//**
 * Make a list of possible option strings
 *
 * @param option the stl vector of options
 * @return a comma separated list of options
******************************************************************************/
string Setup::getXMLOptionList(const vector<string> &options, const string tag) {

    string optionList = "";
    for(auto name : options)
        optionList += str(format("<%s>\n    %s\n<\\%s>\n") % tag % name % tag);
    return optionList;
}

/**************************************************************************//**
 * Define all command line options and get them from the command line.
 *
 * We use boost::program options to get simulation parameters from the 
 * command line.
 * @param argc number of command line arguments
 * @param argv command line string
******************************************************************************/
void Setup::getOptions(int argc, char *argv[])
{
    /* Initialize all possible options */
    initParameters();

    /* Insert all options to the appropriate optionClass map for reading
     * from the command line. */
    params.setupCommandLine(optionClasses);

    /* Add all classes to the complete set */
    for(string oClass : optionClassNames)
        cmdLineOptions.add(optionClasses[oClass]);

    /* Update the parameter map from the commamand line and XML */
    params.update(argc,argv,cmdLineOptions);
}

/**************************************************************************//**
 * Parse the command line options for obvious errors and return values.
 *
 * Here we go through the commmand line options and test for any problems.  
 * This probably needs more work to test all possible outcomes.
 * @return true if we exit, false if we continue
******************************************************************************/
bool Setup::parseOptions() {

    /* Do we need help? */
    if (params("help")) {
        cout << cmdLineOptions << endl;
        return true;
    }

    /* Output the dimension the code was compiled with then exit */
    if (params("dimension")) {
        cout << endl << format("Code was compiled for a %d-dimensional system.") % NDIM 
            << endl << endl;
        return true;
    }

    /* Output the git version that the code was compiled with */
    if (params("version")) {
        cout << endl << format("Code was compiled with repo version %s.") % REPO_VERSION
            << endl << endl;
        return true;
    }

    /* Print a list of allowed estimators in xml format */
    if (params("estimator_list")) {
        cout << endl << "The following estimators can be included in an xml input file." << endl;
        cout << getXMLOptionList(estimatorName, "estimator") << endl << endl;
        return true;
    }

    /* Have we defined a temperature for a PIMC simulation?*/
    if (!params("temperature") && !PIGS ) {
        cerr << endl << "ERROR: No temperature defined!" << endl << endl;
        cerr << "Action: specify temperature (T)" << endl;
        return true;
    }

    /* Have we mistakingly defined one for a PIGS simulation? */
    if (params("temperature") && PIGS) {
        cerr << endl << "ERROR: Did you mean to define a non-zero temperature for PIGS?" << endl << endl;
        cerr << "Action: remove temperature (T)" << endl;
        return true;
    }

    /* Have we physically defined a simulation cell? */
    dVec side;
    definedCell = false;
    if (params("size")) {
        definedCell = true;
        side = params["size"].as<double>();
        params.set<dVec>("side",side);

    }
    else if ((params("Lx") + params("Ly") + params("Lz")) == NDIM) {
        definedCell = true;

        /* The side labels */
        vector<string> sides = {"Lx","Ly","Lz"};
        for (int i = 0; i < NDIM; i++)
            side[i] = params[sides[i]].as<double>();
        params.set<dVec>("side",side);
    }

    /* Make sure we have defined enough options to create the simulation cell */
    if (!( (params("density") && params("number_particles")) ||
         (definedCell && params("number_particles")) ||
         (definedCell && params("density")) ) ) {
        cerr << endl << "ERROR: Cannot create the simulation cell!" << endl << endl;
        cerr << "Action: define a valid simulation cell." << endl;
        cerr << "Need: [number_particles (N) AND density (n)] OR " << endl;
        cerr << "      [number_particles (N) AND size (L) or Lx,Ly,Lz] OR" << endl;
        cerr << "      [size (L) or Lx,Ly,Lz AND density (n)]" << endl << endl;
        cerr << optionClasses["cell"] << endl;
        return true;
    }

    /* Make sure we have selected a valid cell type */
    if (!( (params["geometry"].as<string>() == "cylinder")  || 
           (params["geometry"].as<string>() == "prism") ))
    {
        cerr << endl << "ERROR: Invalid simulation cell type." << endl << endl;
        cerr << "Action: change cell (b) to one of:" << endl
            << "\t[prism,cylinder]" << endl;
        return true;
    }

    /* Can we create the worldlines? */
    if (!((params("imaginary_time_length") && params("number_time_slices")) || 
               (params("imaginary_time_length") && params("imaginary_time_step")) || 
               (params("number_time_slices") && params("imaginary_time_step")) ||
               (params("temperature") && (params("number_time_slices") || 
               (params("imaginary_time_step")))) ) ) {
        cerr << endl << "ERROR: Cannot create imaginary time paths!" << endl << endl;
        cerr << "Action: define imaginary_time_length AND imaginary_time_step (t) OR" << endl;
        cerr << "        define imaginary_time_length AND number_time_steps (P) OR" << endl; 
        cerr << "        define number_time_slices (P) AND imaginary_time_step (t)" << endl << endl;
        cerr << "        define temperature (T) AND (number_time_slices (P) OR imaginary_time_step (t))" << endl << endl;
        cerr << optionClasses["algorithm"] << endl;
        return true;
    }

    /* Make sure we have selected a valid interaction potential */
    if (std::find(interactionPotentialName.begin(), interactionPotentialName.end(), 
                params["interaction"].as<string>()) == interactionPotentialName.end()) {
        cerr << endl << "ERROR: Invalid interaction potential!" << endl << endl;
        cerr << "Action: set interaction (I) to one of:" << endl
             << "\t[" << interactionNames << "]" <<  endl;
        return true;
    }

    /* Make sure we have selected a valid external potential */
    if (std::find(externalPotentialName.begin(), externalPotentialName.end(), 
                params["external"].as<string>()) == externalPotentialName.end()) {
        cerr << endl << "ERROR: Invalid external potential!" << endl << endl;
        cerr << "Action: set external (X) must to one of:" << endl
             << "\t[" << externalNames << "]" << endl;
        return true;
    }

    /* Make sure we have selected a valid action type */
    if (std::find(actionName.begin(), actionName.end(), 
                params["action"].as<string>()) == actionName.end()) {
        cerr << endl << "ERROR: Invalid action!" << endl << endl;
        cerr << "Action: set action to one of:" << endl
             << "\t[" << actionNames << "]" <<  endl;
        return true;
    }
    
    /* Make sure we have selected a valid trial wave fucntion */
    if (std::find(waveFunctionName.begin(), waveFunctionName.end(), 
                params["wavefunction"].as<string>()) == waveFunctionName.end()) {
        cerr << endl << "ERROR: Invalid trial wave function!" << endl << endl;
        cerr << "Action: set wave function to one of :" << endl
        << "\t[" << waveFunctionNames << "]" << endl;
        return true;
    }

    /* Make sure we use the pair product approximation for discontinuous
     * potentials */
    if (((params["interaction"].as<string>() == "hard_sphere") ||
        (params["interaction"].as<string>() == "hard_rod") ||
        (params["interaction"].as<string>() == "delta1D") ) &&
        (params["action"].as<string>() != "pair_product") ) {
        cout << endl;
        cerr << "ERROR: Need to use the pair product approximation with discontinuous potentials!"; 
        cout << endl << endl;
        cerr << "Action: change the action to pair_product." << endl;
        return 1;
    }

    /* We can only use the cylinder potentials for a 3D system */
    if ((params["external"].as<string>().find("tube") != string::npos) && (NDIM != 3)) {
        cerr << endl << "ERROR: Can only use tube potentials for a 3D system!" << endl << endl;
        cerr << "Action: change the potential or recompile with ndim=3." << endl;
        return 1;
    }

    /* Need to specify a radius for the tube potentials */
    if ( (params["external"].as<string>().find("tube") != string::npos) && (!params("radius")) ) {
        cerr << endl << "ERROR: Incomplete specification for external potential!" << endl << endl;
        cerr << "Action: specfiy a radius (r) for the tube potentials." << endl;
        return 1;
    }

    /* Need to specify a y- barrier width scale factor for Gasparini potential */
    if ( (params["external"].as<string>().find("gasp_prim") != string::npos) && 
            (!params("empty_width_y")) ) {
        cerr << endl << "ERROR: Incomplete specification for external potential!" << endl << endl;
        cerr << "Action: specify a y- scale factor (y) for the Gasparini potential." << endl;
        return 1;
    }

    /* Need to specify a z- barrier width scale factor for Gasparini potential */
    if ( (params["external"].as<string>().find("gasp_prim") != string::npos) && 
            (!params("empty_width_z")) ) {
        cerr << endl << "ERROR: Incomplete specification for external potential!" << endl << endl;
        cerr << "Action: specify a z- scale factor (z) for the Gasparini potential." << endl;
        return 1;
    }

    /* We can only use the hard sphere potential in a 3D system */
    if ((params["interaction"].as<string>().find("hard_sphere") != string::npos) && (NDIM != 3)) {
        cerr << endl << "ERROR: Can only use hard sphere potentials for a 3D system!" << endl << endl;
        cerr << "Action: change the potential or recompile with ndim=3." << endl;
        return 1;
    }

    /* We can only use the hard rod potential in a 1D system */
    if ((params["interaction"].as<string>().find("hard_rod") != string::npos) && (NDIM != 1)) {
        cerr << endl << "ERROR: Can only use hard rod potentials for a 1D system!" << endl << endl;
        cerr << "Action: change the potential or recompile with ndim=1." << endl;
        return 1;
    }
    
    /* We can only use the delta1D potential in a 1D system */
    if ((params["interaction"].as<string>().find("delta1D") != string::npos) && (NDIM != 1)) {
        cerr << endl << "ERROR: Can only use delta1D potentials for a 1D system!" << endl << endl;
        cerr << "Action: change the potential or recompile with ndim=1." << endl;
        return 1;
    }

    /* If a list of estimators has been supplied, we need to verify */
    for (string name : params["estimator"].as<vector<string>>()){
        if (std::find(estimatorName.begin(), estimatorName.end(), name) == estimatorName.end()) {
            cerr << "ERROR: Tried to measure a non-existent estimator: " << name << endl;
            cerr << "Action: set estimator to one of:" << endl
                << "\t[" << estimatorNames<< "]" <<  endl;
            return true;
        }
    }

    /* Make sure we don't measure any pigs estimators if we have T > 0 */
    /* if (!PIGS) { */
    /*     for (string name : params["estimator"].as<vector<string>>()) { */
    /*         bool foundPIGSEstimator = (name.find("pigs") != string::npos); */
    /*         if (foundPIGSEstimator) { */
    /*             cerr << "ERROR: Tried to measure a PIGS estimator when T > 0: " << name << endl; */
    /*             cerr << "Action: remove " << name << " estimator." <<  endl; */
    /*             return true; */
    /*         } */
    /*     } */
    /* } */
    /* /1* Make sure all our estimators are pigs estimators *1/ */
    /* else { */
    /*     for (string name : params["estimator"].as<vector<string> >()) { */
    /*         bool foundPIGSEstimator = (name.find("pigs") != string::npos) */
    /*             || (name.find("time") != string::npos); */

    /*         if (!foundPIGSEstimator) { */
    /*             cerr << "ERROR: Tried to measure a non-PIGS estimator when T = 0: " << name << endl; */
    /*             cerr << "Action: remove " << name << " estimator." <<  endl; */
    /*             return true; */
    /*         } */
    /*     } */
    /* } */

    /* Only measure cylinder estimators when we have a cylinder cell */
    if (params["geometry"].as<string>() == "prism") {
        for (string name : params["estimator"].as<vector<string> >()) {
            bool foundCylinderEstimator = (name.find("cylinder") != string::npos);
            if (foundCylinderEstimator) {
                cerr << "ERROR: Tried to measure a cylinder estimator in a prism: " << name << endl;
                cerr << "Action: remove " << name << " estimator." <<  endl;
                return true;
            }
        }
    }

    /* Validate the move list */
    for (string name : params["update"].as<vector<string>>()){
        if (std::find(moveName.begin(), moveName.end(), name) == moveName.end()) {
            cerr << "ERROR: Tried to perform a non-existent move: " << name << endl;
            cerr << "Action: set move to one of:" << endl
                << "\t[" << moveNames << "]" <<  endl;
            return true;
        }
    }

    if (params("validate")) {
        cerr << "SUCCESS: All command line and/or xml options have been verified." << endl;
        cerr << "Action: remove --validate flag to proceed with simulation." << endl;
        return true;
    }

    return false;
}

/**************************************************************************//**
 * Return the random seed.
 *
 * We add the process number to a fixed initial random seed.
 * @param startSeed The fixed initial seed
 * @return A seed shifted by the process number
******************************************************************************/
uint32 Setup::seed (const uint32 startSeed) {
    return startSeed + params["process"].as<uint32>();
}

/**************************************************************************//**
 * Setup the simulation cell.
 *
 * We setup the simulation cell, and return a pointer to a container opject
 * with a type that depends on the specified simulation cell.
******************************************************************************/
Container * Setup::cell() {

    /* Initialize the local box pointer */
    Container *boxPtr = NULL;

    /* Setup a cylindrical simulation cell */
    if (params["geometry"].as<string>() == "cylinder") {

        double radius = params["radius"].as<double>();

        /* We need to adjust the maximal possible radius for a hour glass
         * potential that bows outwards */
        if ( (params["external"].as<string>() == "hg_tube") && 
             (params["hourglass_radius"].as<double>() > 0) )
            radius += params["hourglass_radius"].as<double>();


        if (definedCell && params("number_particles"))
            boxPtr = new Cylinder(radius,params["side"].as<dVec>()[NDIM-1]);
        else if (definedCell && params("density")) {
            boxPtr = new Cylinder(radius,params["side"].as<dVec>()[NDIM-1]);
            params.set<int>("number_particles", int(boxPtr->volume*params["density"].as<double>()));
        }
        else
            boxPtr = new Cylinder(params["density"].as<double>(), 
                    radius,params["number_particles"].as<int>());
    }
    /* Setup a hyperprism */
    else if (params["geometry"].as<string>() == "prism") {

        /* we determine if we are using a non-periodic cell for 
         * the graphene potential */
        iVec periodic;
        periodic = 1;
        if (params["external"].as<string>().find("graphene") != string::npos)
            periodic[NDIM-1] = 0;

        if (definedCell && params("number_particles")) 
            boxPtr = new Prism(params["side"].as<dVec>(),periodic);
        else if (definedCell && params("density")) {
            boxPtr = new Prism(params["side"].as<dVec>(),periodic);
            params.set<int>("number_particles", int(boxPtr->volume*params["density"].as<double>()));
        }
        else
            boxPtr = new Prism(params["density"].as<double>(),params["number_particles"].as<int>());
    }

    /* Add the volume to the parameter map */
    params.set<double>("volume",boxPtr->volume);
    params.set<dVec>("side",boxPtr->side);

    return boxPtr;
}

    /* if ((params["external"].as<string>().find("tube") != string::npos) && (NDIM != 3)) { */
/**************************************************************************//**
 * Setup the worldlines.
 *
 * Depending on whether we have defined the size of the imaginary time step
 * tau or the number of time slices we setup the imaginary time extent of the
 * worldlines.
******************************************************************************/
bool Setup::worldlines() {

    int numTimeSlices,numDeltaTau;
    double tau,imaginaryTimeLength;
    bool pathBreak = (params["number_broken"].as<int>() > 0 || 
            params("spatial_subregion"));

    /* !!NB!! We separate the worldline building into two sections, one for PIGS and
     * one for PIMC.  This can probably be cleaned up in the future. */
    if (!PIGS) {

        /* Set the imaginary time length*/
        if (!params("imaginary_time_length"))
            params.set<double>("imaginary_time_length",1.0/params["temperature"].as<double>());

        /* We determine if we have fixed P or tau.  We require that the number of time 
         * slices is always even. */
        if (!params("number_time_slices")) {
            tau = params["imaginary_time_step"].as<double>();
            numTimeSlices = static_cast<int>(1.0/(params["temperature"].as<double>() * tau) + EPS);
            if ((numTimeSlices % 2) != 0)
                numTimeSlices--;
            params.set<int>("number_time_slices",numTimeSlices);
        }
        else {
            numTimeSlices = params["number_time_slices"].as<int>();
            if ((numTimeSlices % 2) != 0)
                numTimeSlices--;
            tau = 1.0/(params["temperature"].as<double>() * numTimeSlices);
            params.set<int>("number_time_slices",numTimeSlices);
            params.set<double>("imaginary_time_step",tau);
        }
    }
    /* This is for PIGS (T=0) simulations */
    else {
        /* Fixing the total imaginary time and the number of time slices */
        if ( params("imaginary_time_length") && params("number_time_slices") ) {

            /* We make sure we have an odd number of time slices */
            numTimeSlices = params["number_time_slices"].as<int>();
            if ((numTimeSlices % 2) == 0)
                numTimeSlices++;
            /* Make sure the center of the path is an odd (even) slice for a
             closed (open) path calucation */
            numDeltaTau = (numTimeSlices-1)/2;
            if ( ((pathBreak)&&(numDeltaTau%2==0)) || ((!pathBreak)&&(numDeltaTau%2==1))    )
                numTimeSlices += 2;

            tau = params["imaginary_time_length"].as<double>() / (numTimeSlices-1);
            params.set<int>("number_time_slices",numTimeSlices);
            params.set<double>("imaginary_time_step",tau);
        }
        /* Fixing the total imaginary time and the time step.  We need to make sure
         * we can get an odd number of slices */
        else if ( params("imaginary_time_length") && params("imaginary_time_step") ) {
            tau = params["imaginary_time_step"].as<double>();
            numTimeSlices = static_cast<int>((params["imaginary_time_length"].as<double>() / tau) + EPS)+1;

            /* We make sure we have an odd number of time slices */
            if ((numTimeSlices % 2) == 0)
                numTimeSlices++;

            /* Make sure the center of the path is an odd (even) slice for a
             closed (open) path calucation */
            numDeltaTau = (numTimeSlices-1)/2;
            if ( ((pathBreak)&&(numDeltaTau%2==0)) || ((!pathBreak)&&(numDeltaTau%2==1))    )
                numTimeSlices += 2;
            
            imaginaryTimeLength = (numTimeSlices-1)*tau;
            params.set<double>("imaginary_time_length",imaginaryTimeLength);
            
            params.set<int>("number_time_slices",numTimeSlices);
        }
        /* Fixing the number of time steps and the size of the time step.  */
        else {
            /* We make sure we have an odd number of time slices */
            numTimeSlices = params["number_time_slices"].as<int>();
            if ((numTimeSlices % 2) == 0)
                numTimeSlices++;
            /* Make sure the center of the path is an odd (even) slice for a
             closed (open) path calucation */
            numDeltaTau = (numTimeSlices-1)/2;
            if ( ((pathBreak)&&(numDeltaTau%2==0)) || ((!pathBreak)&&(numDeltaTau%2==1))    )
                numTimeSlices += 2;

            /* Set the updated number of time slices and imaginary time length */
            params.set<int>("number_time_slices",numTimeSlices);
            imaginaryTimeLength = (numTimeSlices-1)*params["imaginary_time_step"].as<double>();
            params.set<double>("imaginary_time_length",imaginaryTimeLength);
        }

        /* We set the effective temperature to 1.0/imagTimeLength */
        params.set<double>("temperature",1.0/params["imaginary_time_length"].as<double>());
    }

    /* If we haven't fixed Mbar, do so now. We use the value defined in
     * PRE 74, 036701 (2006). We make sure it is always >= 8 */
    int Mbar = 0;
    if (!params("update_length")) {
        /* Compute the average inter-particle separation */
        double ell = pow((1.0*params["number_particles"].as<int>() /
                    params["volume"].as<double>()),-1.0/(1.0*NDIM));
        double lambda = 24.24 / params["mass"].as<double>();
        Mbar = int(ell*ell/(16*lambda*params["imaginary_time_step"].as<double>()));
        if (Mbar < 2)
            Mbar = 2;
        else if (Mbar > numTimeSlices)
            Mbar = numTimeSlices/2;
        params.set<int>("update_length",Mbar);
    }
    else 
    {
        /* Any negative integer sets the maximum Mbar, useful for free particles */
        if (params["update_length"].as<int>() < 0)
            Mbar = numTimeSlices - 2;
        else
            Mbar = params["update_length"].as<int>();
    }

    /* Now we make sure it is even and not too large*/
    if (Mbar % 2)
        Mbar++;

    params.set<int>("update_length",Mbar);

    if (Mbar > numTimeSlices) {
        cerr << Mbar << " " << numTimeSlices << endl;
        cerr << endl << "ERROR: Update length > number time slices!" << endl << endl;
        cerr << "Action: Increase number_time_slices (P) OR" <<  endl;
        cerr << "        Decrease update_length (M) OR"  << endl;
        cerr << "        Decrease imaginary_time_step (t) OR" << endl; 
        cerr << "        Increase imaginary_time_length" << endl; 
        return true;
    }

    return false;
}

/**************************************************************************//**
 * Setup the simulation constants.
 *
 * Fix all simulation constants.
******************************************************************************/
void Setup::setConstants() {

    /* At present, we need to make sure that if a pair_product action has been
     * selected, that we turn off the cuttoff by making it the size of the box */
    if (params["action"].as<string>() == "pair_product" || 
            !params("potential_cutoff") )
        params.set<double>("potential_cutoff",params["side"].as<dVec>()[NDIM-1]);

    /* Set the required constants */
    constants()->initConstants(params());

    /* If we have specified either the center of mass or displace shift on the
     * command line, we update their values. */
    if (params("com_delta"))
        constants()->setCoMDelta(params["com_delta"].as<double>());
    else {
        params.set<double>("com_delta",constants()->comDelta());
    }
    if (params("displace_delta"))
        constants()->setDisplaceDelta(params["displace_delta"].as<double>());
    else {
        params.set<double>("displace_delta",constants()->displaceDelta());
    }

}

/**************************************************************************//**
 * Setup the communicator.
 *
 * Initialize the communicator, we need to know if we are outputting any config
 * files to disk.  The files are labelled differently depending on whether we
 * are in the canonical or grand-canonical ensemble.  We also need to initialize
 * a possible initial state file and a fixed position file.  Also, since the
 * value of tau we might specifiy at the command line is not the actual one
 * used in the simulation (since the number of time slices must be an integer)
 * we pass it to the communicator for propper labelling of output files.
******************************************************************************/
void Setup::communicator() {
        
    communicate()->init(params["imaginary_time_step"].as<double>(),
            (params["output_config"].as<int>() > 0),params["start_with_state"].as<string>(),
            params["fixed"].as<string>());
}

/*************************************************************************//**
 * Setup the interaction potential.
 *
 * Based on the user's choice we create a new interaction potential pointer
 * which is returned to the main program.  
******************************************************************************/
PotentialBase * Setup::interactionPotential(const Container* boxPtr) {

    PotentialBase *interactionPotentialPtr = NULL;
    if (constants()->intPotentialType() == "free")
        interactionPotentialPtr = new FreePotential();
    else if (constants()->intPotentialType() == "delta")
        interactionPotentialPtr = new DeltaPotential(params["delta_width"].as<double>(),
                params["delta_strength"].as<double>());
    else if (constants()->intPotentialType() == "sutherland")
        interactionPotentialPtr = new SutherlandPotential(params["interaction_strength"].as<double>());
    else if (constants()->intPotentialType() == "hard_sphere")
        interactionPotentialPtr = new HardSpherePotential(params["scattering_length"].as<double>());
    else if (constants()->intPotentialType() == "hard_rod")
        interactionPotentialPtr = new HardRodPotential(params["scattering_length"].as<double>());
    else if (constants()->intPotentialType() == "delta1D")
        interactionPotentialPtr = new Delta1DPotential(params["delta_strength"].as<double>());
    else if (constants()->intPotentialType() == "lorentzian")
        interactionPotentialPtr = new LorentzianPotential(params["delta_width"].as<double>(),
                params["delta_strength"].as<double>());
    else if (constants()->intPotentialType() == "aziz")
        interactionPotentialPtr = new AzizPotential(boxPtr);
    else if (constants()->intPotentialType() == "szalewicz")
        interactionPotentialPtr = new SzalewiczPotential(boxPtr);
    else if (constants()->intPotentialType() == "harmonic")
        interactionPotentialPtr = new HarmonicPotential(params["omega"].as<double>());
    else if (constants()->intPotentialType() == "dipole")
        interactionPotentialPtr = new DipolePotential();

    return interactionPotentialPtr;
}

/*************************************************************************//**
 * Setup the external potential.
 *
 * Based on the user's choice we create a new external potential pointer
 * which is returned to the main program.  
******************************************************************************/
PotentialBase * Setup::externalPotential(const Container* boxPtr) {

    PotentialBase *externalPotentialPtr = NULL;

    if (constants()->extPotentialType() == "harmonic")
        externalPotentialPtr = new HarmonicPotential(params["omega"].as<double>());
    else if (constants()->extPotentialType() == "free")
        externalPotentialPtr = new FreePotential();
    else if (constants()->extPotentialType() == "osc_tube")
        externalPotentialPtr = new HarmonicCylinderPotential(params["radius"].as<double>());
    else if (constants()->extPotentialType() == "hg_tube")
        externalPotentialPtr = new LJHourGlassPotential(params["radius"].as<double>(),
                params["hourglass_radius"].as<double>(), params["hourglass_width"].as<double>());
    else if (constants()->extPotentialType() == "plated_lj_tube")
        externalPotentialPtr = new PlatedLJCylinderPotential(params["radius"].as<double>(),
                params["lj_width"].as<double>(), params["lj_sigma"].as<double>(),
                params["lj_epsilon"].as<double>(), params["lj_density"].as<double>());
    else if (constants()->extPotentialType() == "lj_tube")
        externalPotentialPtr = new LJCylinderPotential(params["radius"].as<double>());
    else if (constants()->extPotentialType() == "hard_tube") 
        externalPotentialPtr = new HardCylinderPotential(params["radius"].as<double>());
    else if (constants()->extPotentialType() == "single_well")
        externalPotentialPtr = new SingleWellPotential();
    else if (constants()->extPotentialType() == "fixed_aziz") 
        externalPotentialPtr = new FixedAzizPotential(boxPtr);
    else if (constants()->extPotentialType() == "gasp_prim")
        externalPotentialPtr = new Gasparini_1_Potential(params["empty_width_z"].as<double>(),
                params["empty_width_y"].as<double>(),boxPtr);
    else if (constants()->extPotentialType() == "graphene") 
        externalPotentialPtr = new GraphenePotential(params["strain"].as<double>(),
                params["poisson"].as<double>(),
                params["carbon_carbon_dist"].as<double>(),
                params["lj_sigma"].as<double>(),
                params["lj_epsilon"].as<double>()
                                );
    else if (constants()->extPotentialType() == "graphenelut") 
        externalPotentialPtr = new GrapheneLUTPotential(params["strain"].as<double>(),
                params["poisson"].as<double>(),
                params["carbon_carbon_dist"].as<double>(),
                params["lj_sigma"].as<double>(),
                params["lj_epsilon"].as<double>(),
                boxPtr
                                );
    else if (constants()->extPotentialType() == "fixed_lj") 
        externalPotentialPtr = new FixedPositionLJPotential(params["lj_sigma"].as<double>(), 
                params["lj_epsilon"].as<double>(),boxPtr);
    else if (constants()->extPotentialType() == "graphenelut3d") 
        externalPotentialPtr = new GrapheneLUT3DPotential(
            params["graphenelut3d_file_prefix"].as<string>(),
            boxPtr
        );
    else if (constants()->extPotentialType() == "graphenelut3dgenerate") 
        externalPotentialPtr = new GrapheneLUT3DPotentialGenerate(
            params["strain"].as<double>(),
            params["poisson"].as<double>(),
            params["carbon_carbon_dist"].as<double>(),
            params["lj_sigma"].as<double>(),
            params["lj_epsilon"].as<double>(),
            boxPtr
        );
    else if (constants()->extPotentialType() == "graphenelut3dtobinary") 
        externalPotentialPtr = new GrapheneLUT3DPotentialToBinary(
            params["graphenelut3d_file_prefix"].as<string>(),
            boxPtr
        );
    else if (constants()->extPotentialType() == "graphenelut3dtotext") 
        externalPotentialPtr = new GrapheneLUT3DPotentialToText(
            params["graphenelut3d_file_prefix"].as<string>(),
            boxPtr
        );

    return externalPotentialPtr;
}

/*************************************************************************//**
* Setup the trial wave function.
*
* Based on the user's choice we create a new trial wave function  pointer
* which is returned to the main program.
******************************************************************************/
WaveFunctionBase * Setup::waveFunction(const Path &path, LookupTable &lookup) {
    
    WaveFunctionBase *waveFunctionPtr = NULL;

    if (constants()->waveFunctionType() == "constant")
        waveFunctionPtr = new WaveFunctionBase(path,lookup);
    else if (constants()->waveFunctionType() == "sech")
        waveFunctionPtr = new SechWaveFunction(path,lookup);
    else if (constants()->waveFunctionType() == "jastrow")
        waveFunctionPtr = new JastrowWaveFunction(path,lookup);
    else if (constants()->waveFunctionType() == "lieb")
        waveFunctionPtr = new LiebLinigerWaveFunction(path,lookup);
    else if (constants()->waveFunctionType() == "sutherland")
        waveFunctionPtr = new SutherlandWaveFunction(path,lookup,
                params["interaction_strength"].as<double>());
    
    return waveFunctionPtr;
}

/*************************************************************************//**
* Setup the action.
*
* Based on the user's choices we create a new action pointer which is returned
* to the main program.
******************************************************************************/
ActionBase * Setup::action(const Path &path, LookupTable &lookup, 
        PotentialBase * externalPotentialPtr, 
        PotentialBase * interactionPotentialPtr, 
        WaveFunctionBase * waveFunctionPtr) {

    ActionBase *actionPtr = NULL;

    /* Non Local Actions */
    if (constants()->actionType() == "pair_product") {
        actionPtr = new NonLocalAction(path,lookup,externalPotentialPtr,
                interactionPotentialPtr,waveFunctionPtr,false,constants()->actionType());   
    }
    /* Local Actions */
    else {

        /* The factors needed for local actions. */
        TinyVector <double,2> VFactor;      
        TinyVector <double,2> gradVFactor;  

        VFactor = 1.0;
        gradVFactor = 0.0;
        int period = 1;

        if (constants()->actionType() == "gsf") {

            /* There are different pre-factors for even/odd slices, we use the 
             * Prokof'ev version. I have taken alpha = 0 here. */
            VFactor[0] = 2.0/3.0;
            VFactor[1] = 4.0/3.0;

            double alpha = 0.0;
            gradVFactor[0] = 2.0*alpha/9.0;
            gradVFactor[1] = 2.0*(1.0-alpha)/9.0;
            period = 2;
        }
        else if (constants()->actionType() == "li_broughton") {
            gradVFactor = 1.0 / 12.0;
            period = 2;
        }

        /* Do we want to use full staging moves? */
        bool local = !params("full_updates");

        actionPtr = new LocalAction(path,lookup,externalPotentialPtr,
                interactionPotentialPtr,waveFunctionPtr,VFactor,gradVFactor,local,
                constants()->actionType(),constants()->endFactor(),period);
    }

    return actionPtr;
}

/*************************************************************************//**
 * Define the Monte Carlo updates that will be performed
 *
 * @param path A reference to the paths
 * @param actionPtr The action in use
 * @param random The random number generator 
 * @return a list of Monte Carlo updates
******************************************************************************/
boost::ptr_vector<MoveBase>* Setup::moves(Path &path, ActionBase *actionPtr, 
        MTRand &random) {

    /* Adapt the move list based on defined parameters */
    vector<string> updates = params["update"].as<vector<string>>();
    if (PIGS) {

        /* If we have broken paths */
        if (params["number_broken"].as<int>() > 0) {

            /* Add the swap break move */
            if (std::find(updates.begin(), updates.end(), SwapBreakMove::name) 
                    == updates.end()) 
                updates.push_back(SwapBreakMove::name);

            params.set<vector<string>>("update",updates);
            constants()->setAttemptProb("diagonal",0.5);
            constants()->setAttemptProb("swap break",0.1);
        
        } else if (params("spatial_subregion")) {

            /* Add the mid staging move */
            if (std::find(updates.begin(), updates.end(), MidStagingMove::name) 
                    == updates.end()) 
                updates.push_back(MidStagingMove::name);

            params.set<vector<string>>("update",updates);
            constants()->setAttemptProb("diagonal",0.5);
            constants()->setAttemptProb("mid staging",0.1);
        }
    }
    else {

        /* We potentially replace bisection with staging */
        if ( params("full_updates") || params("staging") ||
                (params["action"].as<string>() == "pair_product") ||
                (params["max_winding"].as<int>() > 1)  ) 
           {

            /* Replace bisection with staging */
            auto it = std::find(updates.begin(), updates.end(), BisectionMove::name);
            if (it == updates.end()) 
                updates.push_back(StagingMove::name);
            else
                updates.at(std::distance(updates.begin(),it)) = StagingMove::name;

            params.set<vector<string>>("update",updates);
        }
    }

    /* Create a new list of moves that will be returned */
    boost::ptr_vector<MoveBase>* movePtr = new boost::ptr_vector<MoveBase>();

    /* Instatiate the moves */
    for (auto& name : params["update"].as<vector<string>>())
        movePtr->push_back(moveFactory()->Create(name,path,actionPtr,random));
    
    return movePtr;
}

/*************************************************************************//**
 * Create a list of estimators to be measured
 *
 * @see https://stackoverflow.com/questions/39165432/find-last-element-in-stdvector-which-satisfies-a-condition
 * @param path A reference to the paths
 * @param actionPtr The action in use
 * @param random The random number generator 
 * @return a list of estimators
******************************************************************************/
boost::ptr_vector<EstimatorBase> * Setup::estimators(Path &path, 
        ActionBase *actionPtr, MTRand &random) {

    /* Create the list of estimator pointers */
    boost::ptr_vector<EstimatorBase>* estimatorPtr = new boost::ptr_vector<EstimatorBase>();

    /* Instatiate the single path estimators */
    for (auto& name : params["estimator"].as<vector<string>>()) {

        if (name.find("multi") == string::npos)
            estimatorPtr->push_back(estimatorFactory()->Create(name,path,
                        actionPtr,random,params["estimator_radius"].as<double>()));
    }

    /* We determine where a line break is needed for all estimators writing to
     * a common estimator file */
    for (const auto & common : {"estimator","cyl_estimator"}) {
        auto ePtr = find_if(estimatorPtr->rbegin(), estimatorPtr->rend(), 
                [common](EstimatorBase &e) { return e.getLabel() == common; });
        if (ePtr != estimatorPtr->rend()) 
            ePtr->addEndLine();
    }

    return estimatorPtr;
}

/*************************************************************************//**
* Create a list of double path estimators to be measured
*
* @param path A reference to the paths
* @param actionPtr The action in use
* @return a list of double path estimators
******************************************************************************/
boost::ptr_vector<EstimatorBase> * Setup::estimators(
        boost::ptr_vector<Path> &pathPtrVec,
        boost::ptr_vector<ActionBase> &actionPtrVec, MTRand &random) {
    
    /* Create the list of estimator pointers */
    boost::ptr_vector<EstimatorBase>* multiEstimatorPtr = new boost::ptr_vector<EstimatorBase>();

    /* Instatiate the multi path estimators */
    for (auto& name : params["estimator"].as<vector<string>>()) {

        if (name.find("multi") != string::npos)
            multiEstimatorPtr->push_back(multiEstimatorFactory()->Create(name,pathPtrVec[0],
                        pathPtrVec[1],&actionPtrVec[0],&actionPtrVec[1],random,
                        params["estimator_radius"].as<double>()));
    }
    
    return multiEstimatorPtr;
}

/*************************************************************************//**
 * Compare a char array with an option name
 *
 * @param option The commmand line string
 * @param target The target option name
******************************************************************************/
void Setup::cleanCommandLineOptions(int argc, char*argv[], vector<string> &cmdArg,
        vector<string> &cmdSep, vector<string> &cmdVal) {

    /* Get the program name */
    cmdArg.push_back(argv[0]);
    cmdSep.push_back("");
    cmdVal.push_back("");

    for (int n = 1; n < argc; n++) {

        string carg(argv[n]);

        /* Determine if we have a long or short option */
        auto numDashes = std::count(carg.begin(), carg.end(), '-');

        /* Short option */
        if (numDashes == 1) {
            if (carg.length() <= 2) {
                cmdArg.push_back(carg);
                /* make sure we aren't at the end of the list */
                if ((n+1) < argc) {
                    cmdSep.push_back(" ");
                    cmdVal.push_back(string(argv[n+1]));
                }
                else {
                    cmdSep.push_back("");
                    cmdVal.push_back("");
                }
                ++n;
            }
            else {
                cmdArg.push_back(carg.substr(0,2));
                cmdSep.push_back(" ");
                cmdVal.push_back(carg.substr(2,carg.length()-2));
            }
        }
        /* Long option */
        else if (numDashes == 2) {

            /* Do we have an equal sign? */
            auto posEqual = carg.find("="); 
            if (posEqual != string::npos) {
                cmdArg.push_back(carg.substr(0,posEqual));
                cmdSep.push_back("=");
                string nextArg(carg.substr(posEqual+1));

                /* Do we have a space in our value? */
                if (nextArg.find(" ") != string::npos)  
                    cmdVal.push_back(string("\"") + nextArg + string("\""));
                else
                    cmdVal.push_back(nextArg);
            }
            else {
                cmdArg.push_back(carg);

                /* make sure we aren't at the end of the list */
                if ((n+1) < argc) {
                    string nextArg(argv[n+1]);

                    /* is this option a boolean switch? */
                    if (nextArg.find('-') == string::npos) {

                        cmdSep.push_back(" ");

                        /* Do we have a space in our value? */
                        if (nextArg.find(" ") != string::npos)  
                            cmdVal.push_back(string("\"") + nextArg + string("\""));
                        else
                            cmdVal.push_back(nextArg);
                        ++n;
                    }
                    else {
                        cmdSep.push_back("");
                        cmdVal.push_back("");
                    }
                }
                else {
                    cmdSep.push_back("");
                    cmdVal.push_back("");
                }
            }
        }
    }
}

/*************************************************************************//**
 * Output the simulation parameters to a log file.
 *
 * After we have finished equilibrating, we output all the simulation
 * parameters to disk in addition to a command that can be used to 
 * restart the simulation.
 * @param argc The number of command line arguments
 * @param argv The commmand line string
 * @param _seed The random seed
 * @param boxPtr A pointer to the container
 * @param nnGrid The lookup table nearest neighbor grid
******************************************************************************/
void Setup::outputOptions(int argc, char *argv[], const uint32 _seed, 
        const Container *boxPtr, const iVec &nnGrid) {

    /* Pre-process the command line options to make them suitable for log
     * output */
    vector<string> cmdArg, cmdSep, cmdVal;
    cleanCommandLineOptions(argc,argv,cmdArg,cmdSep,cmdVal);

    communicate()->file("log")->stream() << endl << "# ";

    /* Construct the command that would be required to restart the simulation */
    bool outputC0 = false;
    bool outputmu = false;
    bool outputD = false;
    bool outputd = false;
    bool outputEstimator = false;
    bool outputUpdate = false;

    for (uint32 n = 0; n < cmdArg.size(); ++n) {

        auto arg = cmdArg[n];

        if ( (arg == "-s") || (arg == "--start_with_state") )  {
            /* do nothing */
        }
        else if ( (arg == "-C") || (arg == "worm_constant") ) {
            communicate()->file("log")->stream() << format("-C %21.15e ") % constants()->C0();
            outputC0 = true;
        }
        else if ( (arg == "-u") || (arg == "chemical_potential") ) {
            communicate()->file("log")->stream() << format("-u %21.15e ") % constants()->mu();
            outputmu = true;
        }
        else if ((arg == "-p") || (arg == "--process")) {
            communicate()->file("log")->stream() << format("-p %03d ") % params["process"].as<uint32>();
        }
        else if ((arg == "-D") || (arg == "--com_delta")) {
            communicate()->file("log")->stream() << format("-D %21.15e ") % constants()->comDelta();
            outputD = true;
        }
        else if ((arg == "-d") || (arg == "--displace_delta")) {
            communicate()->file("log")->stream() << format("-d %21.15e ") % constants()->displaceDelta();
            outputd = true;
        }
        else if ((arg == "--estimator") || (arg == "-e")) {
            outputEstimator = true;
        }
        else if (arg == "--update") {
            outputUpdate = false;
        }
        else 
            communicate()->file("log")->stream() << cmdArg[n] << cmdSep[n] << cmdVal[n] << " ";
    }

    /* Output the restart flag */
    communicate()->file("log")->stream() << format("-R %s ") % constants()->id();

    /* If we haven't specified the worm constant, output it now */
    if (!outputC0)
        communicate()->file("log")->stream() << format("-C %21.15e ") % constants()->C0();

    /* If we haven't specified the chemical potential , output it now */
    if (!outputmu)
        communicate()->file("log")->stream() << format("-u %21.15e ") % constants()->mu();

    /* If we haven't specified the center of mass Delta, output it now */
    if (!outputD)
        communicate()->file("log")->stream() << format("-D %21.15e ") % constants()->comDelta();
    
    /* If we haven't specified the displace delta, output it now */
    if (PIGS && !outputd)
        communicate()->file("log")->stream() << format("-d %21.15e ") % constants()->displaceDelta();

    /* If we specified estimators, add them to the restart string */
    if (outputEstimator) {
        for (const auto &estName : params["estimator"].as<vector<string>>()) {
            string wrapper("");
            if (estName.find(" ") != string::npos)  
                wrapper = "\"";
            communicate()->file("log")->stream() << "--estimator=" << wrapper << estName << wrapper << " ";
        }
    }

    /* If we specified updates, add them to the restart string */
    if (outputUpdate) {
        for (const auto &mvName : params["updates"].as<vector<string>>())  {
            string wrapper("");
            if (mvName.find(" ") != string::npos)  
                wrapper = "\"";
            communicate()->file("log")->stream() << "--update=" << wrapper << mvName << wrapper << " ";
        }
    }

    communicate()->file("log")->stream() << endl << endl;
    communicate()->file("log")->stream() << "---------- Begin Simulation Parameters ----------" << endl;
    communicate()->file("log")->stream() << endl;

    /* record the full command line string */
    communicate()->file("log")->stream() << format("%-24s\t:\t") % "Command String";

    for (uint32 n = 0; n < cmdArg.size(); n++)
        communicate()->file("log")->stream() << cmdArg[n] << cmdSep[n] << cmdVal[n] << " ";
    communicate()->file("log")->stream() << endl; 

    if (constants()->canonical())
        communicate()->file("log")->stream() << format("%-24s\t:\t%s\n") % "Ensemble" % "canonical";
    else
        communicate()->file("log")->stream() << format("%-24s\t:\t%s\n") % "Ensemble" % "grand canonical";

    if (PIGS)
        communicate()->file("log")->stream() << format("%-24s\t:\t%s\n") % "Simulation Type" % "PIGS";
    else
        communicate()->file("log")->stream() << format("%-24s\t:\t%s\n") % "Simulation Type" % "PIMC";
    communicate()->file("log")->stream() << format("%-24s\t:\t%s\n") % "Action Type" % params["action"].as<string>();
    communicate()->file("log")->stream() << format("%-24s\t:\t%s\n") % "Number of paths" % params["number_paths"].as<int>();
    communicate()->file("log")->stream() << format("%-24s\t:\t%s\n") % "Interaction Potential" % 
        params["interaction"].as<string>();

    /* Ouptut a possible delta function width and strength */
    if ( (params["interaction"].as<string>().find("delta") != string::npos) ||
        (params["interaction"].as<string>().find("lorentzian") != string::npos) ) {
        communicate()->file("log")->stream() << format("%-24s\t:\t%8.3e\n") % "Delta Width"
            % params["delta_width"].as<double>();
        communicate()->file("log")->stream() << format("%-24s\t:\t%-7.2f\n") % "Delta Strength"
            % params["delta_strength"].as<double>();
    }

    /* Ouptut a possible sutherland model interaction strength*/
    if (params["interaction"].as<string>().find("sutherland") != string::npos) {
        communicate()->file("log")->stream() << format("%-24s\t:\t%8.3e\n") % "Interaction Strength"
            % params["interaction_strength"].as<double>();
    }

    /* Output harmonic interaction frequecy */
    if ( (params["interaction"].as<string>().find("harmonic") != string::npos) ) {
        communicate()->file("log")->stream() << format("%-24s\t:\t%-7.2f\n") % "Harmonic Int. Freq."
            % params["omega"].as<double>();
    }

    /* Output a possible scattering length */
    if ( (params["interaction"].as<string>().find("hard_sphere") != string::npos) ||
         (params["interaction"].as<string>().find("hard_rod") != string::npos) ) {
        communicate()->file("log")->stream() << format("%-24s\t:\t%-7.2f\n") % "Scattering Length" 
            % params["scattering_length"].as<double>();
    }
    communicate()->file("log")->stream() << format("%-24s\t:\t%s\n") 
        % "External Potential" % params["external"].as<string>();

    /* output a possible hourglass radius and width for the hourglass potential. */
    if (params["external"].as<string>().find("hg_tube") != string::npos) {
        communicate()->file("log")->stream() << format("%-24s\t:\t%-7.2f\n") 
            % "HourGlass Radius" % params["hourglass_radius"].as<double>();
        communicate()->file("log")->stream() << format("%-24s\t:\t%-7.2f\n") 
            % "HourGlass Width" % params["hourglass_width"].as<double>();
    }

    /* output a possible carbon-carbon distance, strain value and LJ
     * parameters for the graphene potential */
    if (params["external"].as<string>().find("graphene") != string::npos) {
        communicate()->file("log")->stream() << format("%-24s\t:\t%-7.2f\n") 
            % "Carbon Carbon Distance" % params["carbon_carbon_dist"].as<double>();
        communicate()->file("log")->stream() << format("%-24s\t:\t%-7.2f\n") 
            % "Graphene Strain %" % params["strain"].as<double>();
        communicate()->file("log")->stream() << format("%-24s\t:\t%-7.2f\n") 
            % "Graphene Poission Ratio %" % params["poisson"].as<double>();
        communicate()->file("log")->stream() << format("%-24s\t:\t%-7.2f\n") 
            % "Graphene-Carbon LJ Sigma" % params["lj_sigma"].as<double>();
        communicate()->file("log")->stream() << format("%-24s\t:\t%-7.2f\n") 
            % "Graphene-Carbon LJ Epsilon" % params["lj_epsilon"].as<double>();
    }

    /* output possible paramters of the plated LJ cylinder potential */
    if (params["external"].as<string>().find("plated") != string::npos) {
        communicate()->file("log")->stream() << format("%-24s\t:\t%-12.5e\n") 
            % "Plating Radial Width" % params["lj_width"].as<double>();
        communicate()->file("log")->stream() << format("%-24s\t:\t%-12.5e\n") 
            % "Plating LJ Sigma" % params["lj_sigma"].as<double>();
        communicate()->file("log")->stream() << format("%-24s\t:\t%-12.5e\n") 
            % "Plating LJ Epsilon" % params["lj_epsilon"].as<double>();
        communicate()->file("log")->stream() << format("%-24s\t:\t%-12.5e\n") 
            % "Plating LJ Density" % params["lj_density"].as<double>();
    }

    if (PIGS) {
        communicate()->file("log")->stream() << format("%-24s\t:\t%s\n") 
            % "Wavefunction Type" % params["wavefunction"].as<string>();
        communicate()->file("log")->stream() <<
            format("%-24s\t:\t%7.5f\n") % "End Factor" % params["end_factor"].as<double>();
        /* Output possible wave function parameters */
        if ( (params["wavefunction"].as<string>().find("lieb") != string::npos) ) {
            communicate()->file("log")->stream() << format("%-24s\t:\t%-7.2f\n") % "Wavefunction length scale"
                % params["R_LL_wfn"].as<double>();
            communicate()->file("log")->stream() << format("%-24s\t:\t%-7.4f\n") % "Wavefunction wave number"
                % params["k_LL_wfn"].as<double>();
        }
    }
    communicate()->file("log")->stream() << 
        format("%-24s\t:\t%7.5f\n") % "Temperature" % params["temperature"].as<double>();
    communicate()->file("log")->stream() << 
        format("%-24s\t:\t%7.5f\n") % "Chemical Potential" % constants()->mu();
    communicate()->file("log")->stream() << 
        format("%-24s\t:\t%7.5f\n") % "Particle Mass" % params["mass"].as<double>();
    communicate()->file("log")->stream() << 
        format("%-24s\t:\t%d\n") % "Number Time Slices" % constants()->numTimeSlices();
    communicate()->file("log")->stream() << 
        format("%-24s\t:\t%7.5f\n") % "Specified Imaginary Time Step" % params["imaginary_time_step"].as<double>();
    communicate()->file("log")->stream() << 
        format("%-24s\t:\t%7.5f\n") % "Imaginary Time Step" % constants()->tau();
    communicate()->file("log")->stream() << 
        format("%-24s\t:\t%7.5f\n") % "Imaginary Time Length" % constants()->imagTimeLength();
    communicate()->file("log")->stream() << 
        format("%-24s\t:\t%d\n") % "Initial Number Particles" % params["number_particles"].as<int>();
    communicate()->file("log")->stream() << 
        format("%-24s\t:\t%7.5f\n") % "Initial Density" 
        % (1.0*params["number_particles"].as<int>()/boxPtr->volume);
    communicate()->file("log")->stream() <<
        format("%-24s\t:\t%d\n") % "Num. Broken World-lines" % params["number_broken"].as<int>();
    if ( constants()->spatialSubregionOn()){
        communicate()->file("log")->stream() <<
            format("%-24s\t:\t%d\n") % "Spatial Subregion" % params["spatial_subregion"].as<double>();
    }
    communicate()->file("log")->stream() << 
        format("%-24s\t:\t%s\n") % "Container Type" % boxPtr->name;
    communicate()->file("log")->stream() << format("%-24s\t:\t") % "Container Dimensions";
    for (int i = 0; i < NDIM; i++) {
        communicate()->file("log")->stream() << format("%7.5f") % boxPtr->side[i];
        if (i < (NDIM-1))
            communicate()->file("log")->stream() << " x ";
        else
            communicate()->file("log")->stream() << endl;
    }
    communicate()->file("log")->stream() << format("%-24s\t:\t%7.5f\n") % "Container Volume" 
        % boxPtr->volume;
    communicate()->file("log")->stream() << format("%-24s\t:\t") % "Lookup Table";
    for (int i = 0; i < NDIM; i++) {
        communicate()->file("log")->stream() << format("%d") % nnGrid[i];
        if (i < (NDIM-1))
            communicate()->file("log")->stream() << " x ";
        else
            communicate()->file("log")->stream() << endl;
    }
    communicate()->file("log")->stream() <<
        format("%-24s\t:\t%d\n") % "Maximum Winding Sector" % params["max_winding"].as<int>();
    communicate()->file("log")->stream() << format("%-24s\t:\t%7.5f\n") % "Initial Worm Constant" % 
        params["worm_constant"].as<double>();
    communicate()->file("log")->stream() << format("%-24s\t:\t%7.5f\n") % "Worm Constant" % constants()->C0();
    communicate()->file("log")->stream() << format("%-24s\t:\t%7.5f\n") % "Initial CoM Delta" % params["com_delta"].as<double>();
    communicate()->file("log")->stream() << format("%-24s\t:\t%7.5f\n") % "CoM Delta" % constants()->comDelta();
    if (PIGS) {
        communicate()->file("log")->stream() << format("%-24s\t:\t%7.5f\n") % "Initial Displace Delta" % params["displace_delta"].as<double>();
        communicate()->file("log")->stream() << format("%-24s\t:\t%7.5f\n") % "Displace Delta" % constants()->displaceDelta();
    }
    communicate()->file("log")->stream() << 
        format("%-24s\t:\t%d\n") % "Bisection Parameter" % constants()->b();
    communicate()->file("log")->stream() << 
        format("%-24s\t:\t%d\n") % "Update Length" % constants()->Mbar();
    communicate()->file("log")->stream() << 
        format("%-24s\t:\t%7.5f\n") % "Potential Cutoff Length" % params["potential_cutoff"].as<double>();
    communicate()->file("log")->stream() <<
        format("%-24s\t:\t%d\n") % "Bin Size" % params["bin_size"].as<uint32>();
    communicate()->file("log")->stream() << 
        format("%-24s\t:\t%d\n") % "Number EQ Steps" % params["number_eq_steps"].as<uint32>();
    communicate()->file("log")->stream() << 
        format("%-24s\t:\t%d\n") % "Number Bins Stored" % params["number_bins_stored"].as<int>();
    communicate()->file("log")->stream() << format("%-24s\t:\t%d\n") % "Random Number Seed" % _seed;
    if (params["virial_window"].as<int>() == 0){
        communicate()->file("log")->stream() << 
            format("%-24s\t:\t%d\n") % "Virial Window" % constants()->numTimeSlices();
    }
    else{
        communicate()->file("log")->stream() << 
            format("%-24s\t:\t%d\n") % "Virial Window" % params["virial_window"].as<int>();
    }
    communicate()->file("log")->stream() << endl;
    communicate()->file("log")->stream() << "---------- End Simulation Parameters ------------" << endl;
}
