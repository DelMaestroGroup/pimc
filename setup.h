/**
 * @file setup.h 
 * @author Adrian Del Maestro
 * @date 04.06.2011
 *
 * @brief Parses command line input and sets up the details of the simulation.
 */

#ifndef SETUP_H 
#define SETUP_H

#include "common.h"
#include <boost/program_options.hpp>
namespace po = boost::program_options;

class Container; 
class PotentialBase;
// ========================================================================  
// Setup Class
// ========================================================================  
/** 
 * Setup the simulation.
 *
 * This class uses boost::program_options to parse command line input 
 * specifying all possible simulation options.  In the future I would like
 * to implement the use of configuration files for static options as well.
 * We then parse the options, checking for correctness, and define all
 * simulation constants.  We setup the container and potential pointers
 * and finally provide a method to output all the parameters to a log
 * file on disk.
 * 
 * @see http://www.boost.org/doc/libs/release/doc/html/program_options.html
 */
class Setup {
	public:
		Setup();

		/* Get the options from the command line */
		void getOptions(int, char*[]);
		/* Parse the options and check for errors (needs to be more complete) */
		bool parseOptions();

		/* Setup the worldlines */
		bool worldlines();

		/* Setup the physical simulation cell */
		Container *cell();

		/* Setup the simulation constants */
		void setConstants();
		
		/* Setup the communicator */
		void communicator();

		/* Define the random seed */
		uint32 seed(const uint32);
		
		/* Output all options to disk */
		void outputOptions(int, char*[], const uint32, const Container*, const iVec&);

		/* Setup the interaction potential */
		PotentialBase * interactionPotential();

		/* Setup the external potential */
		PotentialBase * externalPotential(const Container*);

        po::variables_map params;       			///< The command line parameter map

	private:
		vector<string> interactionPotentialName;	///< The allowed interaction potential names
		vector<string> externalPotentialName;		///< The allowed external potential names

		string interactionNames;					///< The interaction output list
		string externalNames;						///< The external output list

		dVec side;									///< The sides of the simulation cell
		double volume;								///< A local copy of the volume.

		bool definedCell;							///< The user has physically set the sim. cell

        po::options_description generalOptions; 	///< General options
        po::options_description cellOptions;		///< Simulation cell options
        po::options_description potentialOptions;	///< Potential options
        po::options_description physicalOptions;	///< Physical parameter options
        po::options_description algorithmicOptions;	///< Algorithmic options
        po::options_description cmdLineOptions;		///< All options

		/* Two template methods used to deal with the wonky nature of setting 
		 * parameters programattically */
		template<typename TType> void setOption(string, const TType);
		template<typename TType> void insertOption(string, const TType);

        /* Ensure the correct directory stucture has been created */
//        void createOutputDirectory();
};

/** Use boost::program_options funky method for setting parameters */
template<typename Ttype>
void Setup::setOption(string key, const Ttype val) {
	po::variables_map::iterator it(params.find(key));
	po::variable_value & v(it->second);
	v.value() = val;
}

/** Use boost::program_options funky method for inserting parameters */
template<typename Ttype>
void Setup::insertOption(string key, const Ttype val) {
	params.insert(std::make_pair(key, po::variable_value(val, false)));
	po::notify(params);
}

#endif
