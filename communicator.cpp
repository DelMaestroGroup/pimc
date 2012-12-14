/** 
 * @file communicator.cpp
 * @author Adrian Del Maestro
 *
 * @brief Communicator class implementation.
 */

#include "communicator.h"
#include <boost/filesystem.hpp>

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// FILE CLASS ----------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/**************************************************************************//**
 *  Constructor.
 *
 *  Given a file type/label, ensemble type, data label and out directory,
 *  create a unique file name and backup file name.
 *
 *  @param _type The file type: state, log, etc
 *  @param _data The unique data string identifier
 *  @param ensemble ce: canonical, gce: grand canonical
 *  @param outDir The output directory
******************************************************************************/
File::File(string _type, string _data, string ensemble, string outDir) {

    /* The file name */
    name = str(format("%s/%s-%s-%s.dat") % outDir % ensemble % _type % _data);

    /* Create a backup name */
    bakname = str(format("%s/%s-%s-%s.bak") % outDir % ensemble % _type % _data);
}

/**************************************************************************//**
 *  Constructor.
 *
 *  Create a filename from a string.
 *
 *  @param _name A file name.
******************************************************************************/
File::File(string _name) {

    name = _name;
    bakname = "";

}

/**************************************************************************//**
 *  Close the file.
******************************************************************************/
void File::close() {
    if (rwfile.is_open())
        rwfile.close();
}

/**************************************************************************//**
 *  Open the file.
 *
 *  @param mode A valid file rw mode
******************************************************************************/
void File::open(ios_base::openmode mode) {

	/* Convert the filename to a c string, and open the file */ 
	rwfile.open(name.c_str(), mode);
	if (!rwfile) {
		cerr << "Unable to process file: " << name << endl;
		exit(EXIT_FAILURE);
	}
}

/**************************************************************************//**
 *  Open the file.
 *
 *  @param mode A valid file rw mode
 *  @param _name A valid file name
******************************************************************************/
void File::open(ios_base::openmode mode, string _name) {

	/* Convert the filename to a c string, and open the file */ 
	rwfile.open(_name.c_str(), mode);
	if (!rwfile) {
		cerr << "Unable to process file: " << _name << endl;
		exit(EXIT_FAILURE);
	}
}

/**************************************************************************//**
 *  Reset a file.
 *  
 *  This method is used to prepare a file for writing that is meant to be
 *  overwritten.  In order to be safe, we write to a .bak version of the file. 
******************************************************************************/
void File::reset() {

    /* Close the current file */
    close();

    /* Open a backup file and replace any content if it exists */
    open(ios::out|ios::trunc,bakname);
}

/**************************************************************************//**
 *  Rename a file.
 *  
 *  After we have performed a write to a .bak file, we rename it to the .dat
 *  file using boost::filesystem.
******************************************************************************/
void File::rename() {

    close();

    /* Create the filesystem paths */
    boost::filesystem::path datpath(name);
    boost::filesystem::path bakpath(bakname);

    /* Perform the rename */
    boost::filesystem::rename(bakpath, datpath);

}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// COMMUNICATOR CLASS --------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------


/**************************************************************************//**
 *  Make sure that the current PIMCID doesn't already exist.
 *  
 *  We use boost::filesystem to make sure that the current ID doesn't yet
 *  exist.  If it does, we increment the ID until we have a unique one.
******************************************************************************/
void CCommunicator::getUniqueID(const double _tau) {

    /* Create a test name (we use the log file here) and see if there are any
     * matches. */
    string logName;
    logName = str(format("OUTPUT/%s-log-%s.dat") % ensemble % dataName);
    boost::filesystem::path logPath(logName);

    while(boost::filesystem::exists(logPath)) {

        /* Increment the ID */
        constants()->incid();

        /* Try again */
        logName = str(format("OUTPUT/%s-log-%s.dat") % ensemble % dataName);
        logPath = boost::filesystem::path(logName);

    }

    /* Now we need to update the dataName variable to reflect the new PIMCID */
	if (!constants()->canonical()) {
		dataName = str(format("%06.3f-%07.3f-%+08.3f-%7.5f-%09u") % constants()->T() 
				% constants()->L() % constants()->mu() % _tau % constants()->id());
	}
	else {
		dataName = str(format("%06.3f-%04d-%06.3f-%7.5f-%09d") % constants()->T()
				% constants()->initialNumParticles() 
				% (1.0*constants()->initialNumParticles()/constants()->V()) 
				% _tau % constants()->id());
	}
}

/**************************************************************************//**
 *  Initialize all input/output files.
 *  
 *  If we are in the grand-canonical ensemble, then initialize all output
 *  files with the form gce-xxx-T-L-mu-tau-ID.dat, whereas if we are in
 *  the canonical ensemble we label as ce-xxx-T-N-n-tau-ID.dat.
******************************************************************************/
void CCommunicator::init(double _tau, bool outputWorldline, string _initName,
        string _fixedName)
{

    /* Set local class variables */
    baseDir = "OUTPUT"; 
    initName = _initName;
    fixedName = _fixedName;

    /* Determine the ensemble and unique parameter file string or dataname */
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

    /* Check to make sure the correct directory structure for OUTPUT files is
     * in place.  We do this using boost::filesystem */
    boost::filesystem::path outputPath(baseDir);
    boost::filesystem::create_directory(outputPath);

    /* If we have cylinder output files, add the required directory. */
	if (constants()->extPotentialType().find("tube") != string::npos) {
        boost::filesystem::path cylPath(baseDir + "/CYLINDER");
        boost::filesystem::create_directory(cylPath);
    }

	/* Depending on whether or not we are restarting the simulations, the open mode
	 * changes and header line changes. */
	if (constants()->restart()) {
		mode = ios::out|ios::app;
		header = "";
	}
	else {
		mode = ios::out;

        /* Make sure we have a unique PIMCID using the log file */
    //    getUniqueID(_tau);

        /* Update the header line */
		header =  str(format("# PIMCID: %09d\n") % constants()->id());
	}
}

/**************************************************************************//**
 * Initialze a file based on a type.
******************************************************************************/
void CCommunicator::initFile(string type) {

    /* Check a possible initialization file */
    if (type == "init") {

        /* There are only two reasons we would need an init file, either we are
         * restarting, or starting from a given initialization file */
        if (constants()->restart())
            file_["init"] = new File("state",dataName,ensemble,baseDir);
        else 
            file_["init"] = new File(initName);
        
        file_["init"]->open(ios::in);
    }
	/* Initialize a possible fixed coordinate file */
    else if (type == "fixed") {
		file_[type] = new File(fixedName);
        file_[type]->open(ios::in);
    }
    /* All other file types act normally */
    else {
        string outDir = baseDir;
        string ctype = type;

        /* Deal with possible cylinder output files */
        if (type.find("cyl_") != string::npos) {
            outDir = baseDir + "/CYLINDER";
            ctype.erase(0,4);
        }

        /* Construct the file and open it */
        file_[type] = new File(ctype,dataName,ensemble,outDir);
        file_[type]->open(mode);
        
        /* Write the header line */
        (*file_[type])() << header;
    }
}


/**************************************************************************//**
 *  This public method gets an instance of the Communicator object,  only one
 *  can ever exist at a time.
******************************************************************************/
CCommunicator* CCommunicator::getInstance ()
{   
	static CCommunicator inst;
	return &inst;
}


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
	pairFile_(),radialFile_(), wormFile_(), numberFile_(), positionFile_(),
    windDensFile_(), cylEstimatorFile_(), cylSuperFile_(), cylNumberFile_(),
    cylObdmFile_(), cylPairFile_(),cylPotentialFile_()
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
		exit(EXIT_FAILURE);
	}
}

/**************************************************************************//**
 *  Make sure that the current PIMCID doesn't already exist.
 *  
 *  We use boost::filesystem to make sure that the current ID doesn't yet
 *  exist.  If it does, we increment the ID until we have a unique one.
******************************************************************************/
void Communicator::getUniqueID(const double _tau) {

    /* Create a test name (we use the log file here) and see if there are any
     * matches. */
    string logName;
    logName = str(format("OUTPUT/%s-log-%s.dat") % ensemble % dataName);
    boost::filesystem::path logPath(logName);

    while(boost::filesystem::exists(logPath)) {

        /* Increment the ID */
        constants()->incid();

        /* Try again */
        logName = str(format("OUTPUT/%s-log-%s.dat") % ensemble % dataName);
        logPath = boost::filesystem::path(logName);
    }

    /* Now we need to update the dataName variable to reflect the new PIMCID */
	if (!constants()->canonical()) {
		dataName = str(format("%06.3f-%07.3f-%+08.3f-%7.5f-%09u") % constants()->T() 
				% constants()->L() % constants()->mu() % _tau % constants()->id());
	}
	else {
		dataName = str(format("%06.3f-%04d-%06.3f-%7.5f-%09d") % constants()->T()
				% constants()->initialNumParticles() 
				% (1.0*constants()->initialNumParticles()/constants()->V()) 
				% _tau % constants()->id());
	}
}

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
    file["position"]  = &positionFile_;
    file["windDens"]  = &windDensFile_;

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

    /* Check to make sure the correct directory structure for OUTPUT files is
     * in place.  We do this using boost::filesystem */
    boost::filesystem::path outputPath("OUTPUT");
    boost::filesystem::create_directory(outputPath);

    /* If we have cylinder output files, add the required directory:
     * OUTPUT/CYLINDER */
	if (constants()->extPotentialType().find("tube") != string::npos) {
        boost::filesystem::path cylPath("OUTPUT/CYLINDER");
        boost::filesystem::create_directory(cylPath);
    }

	/* Depending on whether or not we are restarting the simulations, the open mode
	 * changes and header line changes. */
	ios_base::openmode mode;
	string header;
	if (constants()->restart()) {
		mode = ios::out|ios::app;
		header = "";
	}
	else {
		mode = ios::out;

        /* Now we make sure we have a unique PIMCID using the log file */
        getUniqueID(_tau);

        /* Update the header line */
		header =  str(format("# PIMCID: %09d\n") % constants()->id());
	}

    /* The catch-all filename */
	string fileName; 

	/* We define a possible initialization file */
	if (constants()->restart()) {
		file["init"]  = &initFile_;
		fileName = str(format("OUTPUT/%s-%s-%s.dat") % ensemble % "state" % dataName);

        /* Check to make sure that the file exists */
        boost::filesystem::path initName(fileName);
        if (!boost::filesystem::exists(initName)) {
            cerr << "Trying to restart from an unknown PIMCID: " << constants()->id() << endl;
            cerr << fileName << " does not exist!" << endl;
            exit(EXIT_FAILURE);
        }
		openFile(fileName,&initFile_,ios::in);
	}
	else if (!_initName.empty()) {
		file["init"] = &initFile_;
		openFile(_initName,&initFile_,ios::in);
	}

	/* Now we go through all files and open them up except for a possible init
     * file */
	for (map<string,fstream*>::iterator filePtr = file.begin(); 
			filePtr != file.end(); filePtr++) {
        /* Check to make sure we are not init */
        if (filePtr->first != "init") {
            fileName = str(format("OUTPUT/%s-%s-%s.dat") % ensemble % filePtr->first % dataName);
            openFile(fileName,filePtr->second,mode);
            *(filePtr->second) << header;
        }
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
	/* Save the name of the position density file */
	positionName = str(format("OUTPUT/%s-position-%s.dat") % ensemble % dataName);
	/* Save the name of the winding number density file */
	windDensName = str(format("OUTPUT/%s-windDens-%s.dat") % ensemble % dataName);

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
 *  Reset the Particle Position Density Histogram output file.
******************************************************************************/
void Communicator::resetPositionFile(ios_base::openmode mode) {
	/* If the file is open, close it before re-opening */
	if (positionFile_.is_open()) 
		positionFile_.close();
	openFile(positionName,&positionFile_,mode);
}

/**************************************************************************//**
 *  Reset the Winding Number Density Histogram output file.
******************************************************************************/
void Communicator::resetWindDensFile(ios_base::openmode mode) {
	/* If the file is open, close it before re-opening */
	if (windDensFile_.is_open()) 
		windDensFile_.close();
	openFile(windDensName,&windDensFile_,mode);
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
