/**
 * @file communicator.h
 * @author Adrian Del Maestro
 * @date 03.25.2009	
 *
 * @brief Communicator class definition.
 */

#ifndef COMMUNICATOR_H 
#define COMMUNICATOR_H

#include "common.h"
#include "constants.h"
#include <cstring>
#include <fstream>

// ========================================================================  
// Communicator Class
// ========================================================================  
/** 
 * All files used for input/output.
 *
 * Holds information on input and output files that will need to be accessed 
 * throughout the simulation and allows their access via the singleton 
 * design pattern.
 * @see http://en.wikipedia.org/wiki/Singleton_pattern
 */
class Communicator
{
	public:
		static Communicator* getInstance();

		/* Initialize the output files */
		void init(const double, const bool,const string, const string);

		/* Get methods */
		fstream & logFile() {return logFile_;}				///< Get log file
		fstream & estimatorFile() {return estimatorFile_;}	///< Get estimator file
		fstream & superFile() {return superFile_;}			///< Get superfluid file
		fstream & stateFile() {return stateFile_;}			///< Get state file
		fstream & initFile() {return initFile_;}			///< Get initialization file
		fstream & fixedFile() {return fixedFile_;}			///< Get fixed file
		fstream & wlFile() {return wlFile_;}				///< Get worldline file
		fstream & permCycleFile() {return permCycleFile_;}	///< Get permutation cycle file
		fstream & obdmFile() {return obdmFile_;}			///< Get one body density matrix file
		fstream & pairFile() {return pairFile_;}			///< Get pair correlation file
		fstream & radialFile() {return radialFile_;}		///< Get radial density file
		fstream & wormFile() {return wormFile_;}			///< Get worm properties file
		fstream & numberFile() {return numberFile_;}		///< Get number distribution file

		fstream & debugFile();

		/* The cylinder estimator files */
		fstream & cylEstimatorFile() {return cylEstimatorFile_;} ///< Get the cylinder estimator file
		fstream & cylSuperFile() {return cylSuperFile_;}		 ///< Get the cylinder superfluid file
		fstream & cylNumberFile() {return cylNumberFile_;}		 ///< Get the cylinder number file 
		fstream & cylObdmFile() {return cylObdmFile_;}			 ///< Get the cylinder OBDM file
		fstream & cylPairFile() {return cylPairFile_;}			 ///< Get the cylinder pair CF file
		fstream & cylPotentialFile() { return cylPotentialFile_;}///< Get the cylinder potential file

		/* Reset the state file */
		void resetStateFile(ios_base::openmode);	

		/* Reset the fixed file */
		void resetFixedFile();	
		
	protected:
		Communicator();									
		~Communicator();								
		Communicator(const Communicator&);				///< Copy constructor
		Communicator& operator= (const Communicator&);	///< Singleton equals

	private:
		fstream logFile_;			// The logging simulation parameters file
		fstream estimatorFile_;		// Single value estimators
		fstream superFile_;			// Superfluid estimators
		fstream debugFile_;			// Debugging output
		fstream stateFile_;			// Diagonal path information
		fstream initFile_;			// Initial configuration state
		fstream fixedFile_;			// Location of fixed non-updatable particles
		fstream wlFile_;			// Particle worldlines
		fstream permCycleFile_;		// Permutation cycle statistics
		fstream obdmFile_;			// One body density matrix
		fstream pairFile_;			// Pair correlation function
		fstream radialFile_;		// Radial density
		fstream wormFile_;			// Worm properties
		fstream numberFile_;		// The number distribution

		/* All the cylinder estimator files */
		fstream cylEstimatorFile_;
		fstream cylSuperFile_;			
		fstream cylNumberFile_;		
		fstream cylObdmFile_;			
		fstream cylPairFile_;			
		fstream cylPotentialFile_;			

		string ensemble;			// The type of ensemble
		string dataName;			// The labelling scheme of the output files

		string stateName;			// The name of the state file
		string fixedName;			// The name of the fixed file

		map <string,fstream*> file;	    // The file map
		map <string,fstream*> cylFile;	// The cylinder file map

		void openFile(const string , fstream *, ios_base::openmode);
};


/**************************************************************************//**
 *  Global public access to the communcator singleton.
******************************************************************************/
inline Communicator* communicate() {
	Communicator *temp = Communicator::getInstance();
	return temp;
}
#endif

