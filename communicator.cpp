/** 
 * @file communicator.cpp
 * @author Adrian Del Maestro
 *
 * @brief Communicator class implementation.
 */

#include "communicator.h"
//#include <boost/filesystem.hpp>

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// COMMUNICATOR CLASS --------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/**************************************************************************//**
 *  Constructor.
******************************************************************************/
Communicator::Communicator() : logFile_(),estimatorFile_(),superFile_(),debugFile_(),
	stateFile_(),initFile_(), fixedFile_(), wlFile_(),permCycleFile_(),obdmFile_(),
	pairFile_(),radialFile_(), wormFile_(), numberFile_(), cylEstimatorFile_(), 
	cylSuperFile_(), cylNumberFile_(), cylObdmFile_(), cylPairFile_(),cylPotentialFile_()
{ 

}

/**************************************************************************//**
 *  We go through all files, and if any are open, we close them.
******************************************************************************/
Communicator::~Communicator() { 

	for (map<string,fstream*>::iterator filePtr = file.begin(); 
			filePtr != file.end(); filePtr++) {
		if ((filePtr->second)->is_open())
			(filePtr->second)->close();
	}

	for (map<string,fstream*>::iterator filePtr = cylFile.begin(); 
			filePtr != cylFile.end(); filePtr++) {
		if ((filePtr->second)->is_open())
			(filePtr->second)->close();
	}

}

/**************************************************************************//**
 *  Test and attempt to open a supplied file.
******************************************************************************/
void Communicator::openFile(const string fileName, fstream *_file, 
		ios_base::openmode mode) {

	/* Convert the filename to a c string, and open the file */ 
	_file->open(fileName.c_str(), mode);
	if (!(*_file)) {
		cerr << "Unable to process file: " << fileName << endl;
		exit(1);
	}
}

/**************************************************************************//**
 *  Make sure that the current PIMCID doesn't already exist.
 *  
 *  We use boost::filesystem to make sure that the current ID doesn't yet
 *  exist.  If it does, we increment the ID until we have a unique one.
******************************************************************************/
//void Communicator::getUniqueID() {
//
//
//    /* Create a test name (we use the log file here) and see if there are any
//     * matches. */
//    string testName;
//    testName = str(format("OUTPUT/*log-*-%09u.dat") % constants()->id());
//    boost::filesystem::path logpath(testName);
//
//    cout << logpath << endl;
//
//    while(boost::filesystem::exists(logpath)) {
//
//        /* Increment the ID */
//        constants()->incID();
//
//        /* Try again*/
//        testName = str(format("OUTPUT/*log-*-%09u.dat") % constants()->id());
//        logpath = boost::filesystem::path(testName);
//    }
//}

/**************************************************************************//**
 *  Initialize all input/output files.
 *  
 *  If we are in the grand-canonical ensemble, then initialize all output
 *  files with the form gce-xxx-T-L-mu-tau-ID.dat, whereas if we are in
 *  the canonical ensemble we label as ce-xxx-T-N-n-tau-ID.dat.
******************************************************************************/
void Communicator::init(const double _tau, const bool outputWorldline, const
        string _initName, const string _fixedName) {

	/* Fill up the file map */
	file["log"]       = &logFile_;
	file["estimator"] = &estimatorFile_;
	file["super"]     = &superFile_;
	file["pcycle"]    = &permCycleFile_;
	file["obdm"]      = &obdmFile_;
	file["pair"]      = &pairFile_;
	file["worm"]      = &wormFile_;
	file["number"]    = &numberFile_;
	file["state"]     = &stateFile_;

	/* We only include the radial density estimator if NDIM > 1 */
	if (NDIM > 1)
		file["radial"]    = &radialFile_;

	if (outputWorldline) 
		file["wl"]  = &wlFile_;

	/* We only add cylinder files if we have a cylindrical potential */
	if (constants()->extPotentialType().find("tube") != string::npos) {
		cylFile["estimator"] = &cylEstimatorFile_;
		cylFile["super"]     = &cylSuperFile_;
		cylFile["number"]    = &cylNumberFile_;
		cylFile["obdm"]      = &cylObdmFile_;
		cylFile["pair"]      = &cylPairFile_;
		cylFile["potential"] = &cylPotentialFile_;
	}

	if (!constants()->canonical()) {
		ensemble = "gce";
		dataName = str(format("%06.3f-%07.3f-%+08.3f-%7.5f-%09u") % constants()->T() 
				% constants()->L() % constants()->mu() % _tau % constants()->id());
	}
	else {
		ensemble = "ce";
		dataName = str(format("%06.3f-%04d-%06.3f-%7.5f-%09d") % constants()->T()
				% constants()->initialNumParticles() 
				% (1.0*constants()->initialNumParticles()/constants()->V()) 
				% _tau % constants()->id());
	}

	/* Depending out wheter or not we are restarting the simulations, the open mode
	 * changes and first line changes. */
	ios_base::openmode mode;
	string header;
	if (constants()->restart()) {
		mode = ios::out|ios::app;
		header = "";
	}
	else {
		mode = ios::out;
		header =  str(format("# PIMCID: %09d\n") % constants()->id());
	}

	/* Now we go through all files and open them up */
	string fileName;
	for (map<string,fstream*>::iterator filePtr = file.begin(); 
			filePtr != file.end(); filePtr++) {
		fileName = str(format("OUTPUT/%s-%s-%s.dat") % ensemble % filePtr->first % dataName);
		openFile(fileName,filePtr->second,mode);
		*(filePtr->second) << header;
	}

	/* Now do the same for any possible cylinder files */
	for (map<string,fstream*>::iterator filePtr = cylFile.begin(); 
			filePtr != cylFile.end(); filePtr++) {
		fileName = str(format("OUTPUT/CYLINDER/%s-%s-%s.dat") % ensemble % filePtr->first % dataName);
		openFile(fileName,filePtr->second,mode);
		*(filePtr->second) << header;
	}

	/* Save the name of the state file */
	stateName = str(format("OUTPUT/%s-state-%s.dat") % ensemble % dataName);

	/* We define a possible init file */
	if (constants()->restart()) {
		file["init"]  = &initFile_;
		fileName = str(format("OUTPUT/%s-%s-%s.dat") % ensemble % "state" % dataName);
		openFile(fileName,&initFile_,ios::in);
	}
	else if (!_initName.empty()) {
		file["init"] = &initFile_;
		openFile(_initName,&initFile_,ios::in);
	}

	/* Initialize a possible fixed coordinate file */
	if (!_fixedName.empty()) {
		fixedName = _fixedName;
		file["fixed"] = &fixedFile_;
		openFile(fixedName,&fixedFile_,ios::in);
	}
}

/**************************************************************************//**
 * Get the debug file if it is needed.
******************************************************************************/
fstream & Communicator::debugFile() {
	if (!debugFile_.is_open()) {
		file["debug"] = &debugFile_;
	    string fileName;
		fileName = str(format("OUTPUT/%s-%s-%s.dat") % ensemble % "debug" % dataName);
		openFile(fileName,&debugFile_,ios::out);
		debugFile_ << str(format("# PIMCID: %09d\n") % constants()->id());
	}

	return debugFile_;
}

/**************************************************************************//**
 *  Reset the state file and prepare based on mode.
******************************************************************************/
void Communicator::resetStateFile(ios_base::openmode mode) {
	/* If the file is open, close it before re-opening */
	if (stateFile_.is_open()) 
		stateFile_.close();
	openFile(stateName,&stateFile_,mode);
}

/**************************************************************************//**
 *  Reset the fixed file and prepare based on mode.
******************************************************************************/
void Communicator::resetFixedFile() {

	/* If the file is open, close it before re-opening */
	if (fixedFile_.is_open()) 
		fixedFile_.close();
	openFile(fixedName,&fixedFile_,ios::in);
}

/**************************************************************************//**
 *  This public method gets an instance of the Communicator object,  only one
 *  can ever exist at a time.
******************************************************************************/
Communicator* Communicator::getInstance ()
{   
	static Communicator inst;
	return &inst;
}
