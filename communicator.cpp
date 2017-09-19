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
File::File(string _name) : name(_name), bakname() {

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

    /* Write the generic header to the file */
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
void Communicator::init(double _tau, bool outputWorldline, string _initName,
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
        getUniqueID(_tau);

        /* Update the header line */
        header =  str(format("# PIMCID: %09d\n") % constants()->id());
    }
}

/**************************************************************************//**
 * Initialze a file based on a type.
******************************************************************************/
void Communicator::initFile(string type) {

    /* Check a possible initialization file */
    if (type.find("init") != string::npos ) {

        /* We need to determine the name of the state file.  i.e. does it need
         * an integer appended after it? */
        string stateName = "state";

        /* If we have a numerical label, append it to the name of state */
        if (type != "init")
            stateName += stateName.substr(4,string::npos);

        /* There are only two reasons we would need an init file, either we are
         * restarting, or starting from a given initialization file */
        if (constants()->restart())
            file_.insert(type, new File(stateName,dataName,ensemble,baseDir));
        else 
            file_.insert(type, new File(initName));
        
        file_.at(type).open(ios::in);
    }
    /* Initialize a possible fixed coordinate file */
    else if (type == "fixed") {
        file_.insert(type,new File(fixedName));
        file_.at(type).open(ios::in);
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
        file_.insert(type, new File(ctype,dataName,ensemble,outDir));
        file_.at(type).open(mode);
        
        /* Write the header line */
        file_.at(type).stream() << header;
    }
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
