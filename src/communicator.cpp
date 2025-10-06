/** 
 * @file communicator.cpp
 * @author Adrian Del Maestro
 *
 * @brief Communicator class implementation.
 */

#include "communicator.h"

/* Filesystem is tricky as it is not yet widely supported.  We try to address
 * that here. */
#if defined(__clang__) || (defined(__GNUC__) && (__GNUC__ > 8))
    #include <filesystem>
    namespace fs = std::filesystem;
#else
#if defined(__GNUC__) && (__GNUC__ > 6)
        #include <experimental/filesystem>
        namespace fs = std::experimental::filesystem;
#endif
#endif

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// BASE FILE CLASS -----------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
BaseFile::BaseFile(std::string _name) : name(_name), bakname() {}

void BaseFile::rename(){
	close();
	fs::rename(bakname.c_str(), name.c_str());
}

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
 *  @param _data The unique data std::string identifier
 *  @param ensemble ce: canonical, gce: grand canonical
 *  @param outDir The output directory
******************************************************************************/
File::File(std::string _type, std::string _data, std::string ensemble, std::string outDir) {

    /* The file name */
	std::string newName = str(format("%s/%s-%s-%s.dat") % outDir % ensemble % _type % _data);
	setName(newName);

    /* Create a backup name */
	std::string newBakname = str(format("%s/%s-%s-%s.bak") % outDir % ensemble % _type % _data);
	setBakname(newBakname);

    /* Determine if the file already exists */
    exists_ = fs::exists(name);

    /* Has the file been prepared for writing? */
    prepared_ = false;
}

/**************************************************************************//**
 *  Constructor.
 *
 *  Create a filename from a std::string.
 *
 *  @param _name A file name.
******************************************************************************/
File::File(std::string _name) : BaseFile(_name) {
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
void File::open(OpenMode modeEnum) {
	std::ios_base::openmode mode = OpenModeTranslator(modeEnum);

    /* Convert the filename to a c std::string, and open the file */ 
    rwfile.open(name.c_str(), mode);
    if (!rwfile) {
        std::cerr << "Unable to process file: " << name << std::endl;
        exit(EXIT_FAILURE);
    }
}

// This was done to provide a common interface between plain-text and HDF5 files.
// OpenMode enum is in communicator.h
std::ios_base::openmode File::OpenModeTranslator(OpenMode modeEnum){
	std::ios_base::openmode mode = std::ios::in; // default to read
	switch (modeEnum) {
		case OpenMode::READ:
			mode = std::ios::in;
			break;
		case OpenMode::READ_WRITE:
			mode = std::ios::in | std::ios::out;
			break;
		case OpenMode::WRITE:
			mode = std::ios::out;
			break;
		case OpenMode::WRITE_TRUNC:
			mode = std::ios::out | std::ios::trunc;
			break;
		case OpenMode::WRITE_APPEND:
			mode = std::ios::out | std::ios::app;
			break;
	}
	return mode;
}

/**************************************************************************//**
 *  Open the file.
 *
 *  @param mode A valid file rw mode
 *  @param _name A valid file name
******************************************************************************/
void File::open(OpenMode modeEnum, std::string _name) {
	std::ios_base::openmode mode = OpenModeTranslator(modeEnum);

    /* Convert the filename to a c std::string, and open the file */ 
    rwfile.open(_name.c_str(), mode);
    if (!rwfile) {
        std::cerr << "Unable to process file: " << _name << std::endl;
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
    open(OpenMode::WRITE_TRUNC,bakname);

    /* Write the generic header to the file */
}

/**************************************************************************//**
 *  Rename a file.
 *  
 *  After we have performed a write to a .bak file, we rename it to .dat
******************************************************************************/
void File::rename() {

    close();

    /* Perform the rename */
    fs::rename(bakname.c_str(), name.c_str());
}

void File::write(std::string data) {
	stream() << data;	
	stream().flush();
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// HDF5 FILE CLASS --------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
/**************************************************************************//**
 *  Constructor.
 *
 *  Given a file type/label, ensemble type, data label and out directory,
 *  create a unique file name and backup file name.
 *
 *  @param _type The file type: state, log, etc
 *  @param _data The unique data std::string identifier
 *  @param ensemble ce: canonical, gce: grand canonical
 *  @param outDir The output directory
******************************************************************************/
HDF5File::HDF5File(std::string _type, std::string _data, std::string ensemble, std::string outDir) {

    /* The file name */
	std::string newName = str(format("%s/%s-%s-%s.hdf5") % outDir % ensemble % _type % _data);
	setName(newName);

    /* Create a backup name */
	std::string newBakname = str(format("%s/%s-%s-%s.bak") % outDir % ensemble % _type % _data);
	setBakname(newBakname);

    /* Determine if the file already exists */
    exists_ = fs::exists(name);

    /* Has the file been prepared for writing? */
    prepared_ = false;
	/* Initialize file handle to invalid */
	file_id = H5I_INVALID_HID;
}

HDF5File::HDF5File(std::string _name) : BaseFile(_name) {
	file_id = H5I_INVALID_HID;
}

// This was done to provide a common interface between plain-txt and HDF5 files. 
// OpenMode enum is in communicator.h
unsigned HDF5File::OpenModeTranslator(OpenMode modeEnum){
	unsigned mode = H5F_ACC_RDONLY; // default to read-only
	switch (modeEnum) {
		case OpenMode::READ:
			mode = H5F_ACC_RDONLY;
			break;
		case OpenMode::READ_WRITE:
			mode = H5F_ACC_RDWR;
			break;
		case OpenMode::WRITE:
			mode = H5F_ACC_TRUNC;
			break;
		case OpenMode::WRITE_TRUNC:
			mode = H5F_ACC_TRUNC;
			break;
		case OpenMode::WRITE_APPEND:
			mode = H5F_ACC_TRUNC;
			std::cerr << "ERROR: Attempting to append to an HDF5File.\n";
			exit(EXIT_FAILURE);
			break;
	}
	return mode;
}

void HDF5File::open(OpenMode modeEnum) {
	unsigned mode = OpenModeTranslator(modeEnum);
	if (fs::exists(name)){
		// file already exists
		file_id = H5Fopen(name.c_str(), mode, H5P_DEFAULT);
	} else {
		// file does not yet exist
		file_id = H5Fcreate(name.c_str(), mode, H5P_DEFAULT, H5P_DEFAULT);
	}

	if (file_id == H5I_INVALID_HID){
        std::cerr << "Unable to process file: " << name << "\n";
        exit(EXIT_FAILURE);
	}
}

void HDF5File::open(OpenMode modeEnum, std::string _name) {
	unsigned mode = OpenModeTranslator(modeEnum);
	file_id = H5Fcreate(_name.c_str(), mode, H5P_DEFAULT, H5P_DEFAULT);
	if (file_id == H5I_INVALID_HID){
        std::cerr << "Unable to process file: " << _name << "\n";
        exit(EXIT_FAILURE);
	}

}

// ... template implementations moved to include/communicator_impl.tpp

void HDF5File::reset(){
	close();
	open(OpenMode::WRITE_TRUNC, bakname);
}

void HDF5File::rename(){
	fs::rename(bakname.c_str(), name.c_str());
}

/*
 * This method takes a string which signifies the target key.
 */
// ... template implementations moved to include/communicator_impl.tpp

void HDF5File::writeState(const StateObject state){
	writeHDF5(state.worldOrder);
	writeHDF5(state.acceptanceInfo);
	writeHDF5(state.moveAcceptance);
	writeHDF5(state.estimatorSampling);
	writeHDF5(state.paths);
	writeHDF5(state.next);
	writeHDF5(state.prev);
	writeHDF5(state.worm);
	writeHDF5(state.randomState);
}

StateObject HDF5File::readState(){
	StateObject state;

	// Read datasets and populate the StateObject.data fields. Names in the
	// StateObject omit the leading '/'.
	try {
		// Helper to obtain dataset dimensions (will open file if needed)
		auto getDims = [&](const char *dname)->std::vector<hsize_t> {
			bool openedHere = false;
			if (file_id == H5I_INVALID_HID) { open(OpenMode::READ); openedHere = true; }
			hid_t did = H5Dopen2(file_id, dname, H5P_DEFAULT);
			if (did == H5I_INVALID_HID) {
				if (openedHere) H5Fclose(file_id);
				throw std::runtime_error(std::string("Unable to open dataset for dims: ") + dname);
			}
			hid_t sid = H5Dget_space(did);
			int ndims = H5Sget_simple_extent_ndims(sid);
			std::vector<hsize_t> dims(ndims);
			H5Sget_simple_extent_dims(sid, dims.data(), NULL);
			H5Sclose(sid);
			H5Dclose(did);
			if (openedHere) H5Fclose(file_id);
			return dims;
		};

		state.worldOrder.data = readDataset<uint32>("/worldOrder");
		state.worldOrder.name = "worldOrder";
		state.worldOrder.type = DataType::UINT_32;
		try { state.worldOrder.dims = getDims("/worldOrder"); } catch(...) {}

		state.acceptanceInfo.data = readDataset<uint32>("/acceptanceInfo");
		state.acceptanceInfo.name = "acceptanceInfo";
		state.acceptanceInfo.type = DataType::UINT_32;
		try { state.acceptanceInfo.dims = getDims("/acceptanceInfo"); } catch(...) {}

		state.moveAcceptance.data = readDataset<uint32>("/moveAcceptance");
		state.moveAcceptance.name = "moveAcceptance";
		state.moveAcceptance.type = DataType::UINT_32;
		try { state.moveAcceptance.dims = getDims("/moveAcceptance"); } catch(...) {}

		state.estimatorSampling.data = readDataset<uint32>("/estimatorSampling");
		state.estimatorSampling.name = "estimatorSampling";
		state.estimatorSampling.type = DataType::UINT_32;
		try { state.estimatorSampling.dims = getDims("/estimatorSampling"); } catch(...) {}

		state.paths.data = readDataset<double>("/paths");
		state.paths.name = "paths";
		state.paths.type = DataType::DOUBLE;
		state.paths.dims = getDims("/paths");

		state.next.data = readDataset<double>("/next");
		state.next.name = "next";
		state.next.type = DataType::DOUBLE;
		try { state.next.dims = getDims("/next"); } catch(...) {}

		state.prev.data = readDataset<double>("/prev");
		state.prev.name = "prev";
		state.prev.type = DataType::DOUBLE;
		try { state.prev.dims = getDims("/prev"); } catch(...) {}

	// worm and randomState use uint32 in this codebase
	state.worm.data = readDataset<uint32>("/worm");
	state.worm.name = "worm";
	state.worm.type = DataType::UINT_32;
	try { state.worm.dims = getDims("/worm"); } catch(...) {}

	state.randomState.data = readDataset<uint32>("/randomState");
	state.randomState.name = "randomState";
	state.randomState.type = DataType::UINT_32;
	try { state.randomState.dims = getDims("/randomState"); } catch(...) {}

		std::cout << "Finished reading state from file: " << name << std::endl;

	} catch (...) {
		std::cerr << "Exception while reading state from HDF5 file: " << name << std::endl;
		exit(EXIT_FAILURE);
	}

	return state;
}


// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// COMMUNICATOR CLASS --------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/**************************************************************************//**
 *  Initialize all input/output files.
 *  
 *  If we are in the grand-canonical ensemble, then initialize all output
 *  files with the form gce-xxx-T-L-mu-tau-ID.dat, whereas if we are in
 *  the canonical ensemble we label as ce-xxx-T-N-n-tau-ID.dat.
******************************************************************************/
void Communicator::init(double _tau, bool outputWorldline, std::string _initName,
        std::string _fixedName)
{

    /* Set local class variables */
    baseDir = "OUTPUT"; 
    initName = _initName;
    fixedName = _fixedName;
    tau = _tau;


    /* Determine the ensemble and unique parameter file std::string or dataname */
    if (!constants()->canonical()) {
        ensemble = "gce";
        dataName = str(format("%06.3f-%07.3f-%+08.3f-%7.5f-%s") % constants()->T() 
                % constants()->L() % constants()->mu() % tau % constants()->id());
    }
    else {
        ensemble = "ce";
        dataName = str(format("%06.3f-%04d-%06.3f-%7.5f-%s") % constants()->T()
                % constants()->initialNumParticles() 
                % (1.0*constants()->initialNumParticles()/constants()->V()) 
                % tau % constants()->id());
    }

    /* Check to make sure the correct directory structure for OUTPUT files is
     * in place. */ 
    fs::path outputPath(baseDir);
    fs::create_directory(outputPath);

    /* If we have cylinder output files, add the required directory. */
    if (constants()->extPotentialType().find("tube") != std::string::npos) {
        fs::path cylPath(baseDir + "/CYLINDER");
        fs::create_directory(cylPath);
    }

    /* A header line for the files */
    header =  str(format("# PIMCID: %s\n") % constants()->id());

    /* Depending on whether or not we are restarting the simulations, the open mode
     * changes. */
    if (constants()->restart()) {
        mode = OpenMode::WRITE_APPEND;  // This open mode is not allowed for HDF5 files.
    }
    else {
        mode = OpenMode::WRITE;
    }
}

/**************************************************************************//**
 * Initialze a file based on a type.
******************************************************************************/
void Communicator::initFile(std::string type) {

    /* Check a possible initialization file */
    if (type.find("init") != std::string::npos ) {
        /* We need to determine the name of the state file.  i.e. does it need
         * an integer appended after it? */
        std::string stateName = "state";

        /* If we have a numerical label, append it to the name of state */
        if (type != "init")
            stateName += stateName.substr(4, std::string::npos);

        /* There are only two reasons we would need an init file, either we are
         * restarting, or starting from a given initialization file 
		 */
        if (constants()->restart())
            file_.insert(type, new File(stateName, dataName, ensemble, baseDir));
        else 
            file_.insert(type, new File(initName));
        
        file_.at(type).open(OpenMode::READ);
    }
    /* Initialize a possible fixed coordinate file */
    else if (type == "fixed") {
        file_.insert(type,new File(fixedName));
        file_.at(type).open(OpenMode::READ);
    }
    /* All other file types act normally */
    else {
        std::string outDir = baseDir;
        std::string ctype = type;

        /* Deal with possible cylinder output files */
        if (type.find("cyl_") != std::string::npos) {
            outDir = baseDir + "/CYLINDER";
            ctype.erase(0,4);
        }

		if (type.find("state") != std::string::npos) {
			file_.insert(type, new File(ctype, dataName, ensemble, outDir));
			file_.at(type).open(mode);
		} else {
			/* Construct the file and open it */
			file_.insert(type, new File(ctype,dataName,ensemble,outDir));
			file_.at(type).open(mode);
			/* Write the header line if the file doesn't exist */
			if (!file_.at(type).exists())
				file_.at(type).write(header);

		}

    }
}

/**************************************************************************//**
 * Initialize an HDF5 file based on a type.
******************************************************************************/
void Communicator::initHDF5File(std::string type) {
    /* Check a possible initialization file */
    if (type.find("init") != std::string::npos ) {
        /* We need to determine the name of the state file.  i.e. does it need
         * an integer appended after it? */
        std::string stateName = "state";

        /* If we have a numerical label, append it to the name of state */
        if (type != "init")
            stateName += stateName.substr(4, std::string::npos);

        /* There are only two reasons we would need an init file, either we are
         * restarting, or starting from a given initialization file 
		 */
		if (constants()->restart()){
			std::cout << "initializing a new HDF5File object. dataName: " << dataName.c_str() << std::endl;
            hdf5File_.insert(type, new HDF5File(stateName, dataName, ensemble, baseDir));
		}
		else {
			std::cout << "initializing a new HDF5File object." << std::endl;
			hdf5File_.insert(type, new HDF5File(initName));
		}
	}
    else if (type == "fixed") {
    }
    else {
        std::string outDir = baseDir;
        std::string ctype = type;
		if (type.find("state") != std::string::npos) {
			hdf5File_.insert(type, new HDF5File(ctype, dataName, ensemble, outDir));
			hdf5File_.at(type).open(mode);
		} else {
			
		}
	}
}

/**************************************************************************//**
 * Update the data name and rename any existing files
******************************************************************************/
void Communicator::updateNames() {

    /* Determine the ensemble and unique parameter file std::string or dataname */
    if (!constants()->canonical()) {
        dataName = str(format("%06.3f-%07.3f-%+08.3f-%7.5f-%s") % constants()->T() 
                % constants()->L() % constants()->mu() % tau % constants()->id());
    }
    else {
        dataName = str(format("%06.3f-%04d-%06.3f-%7.5f-%s") % constants()->T()
                % constants()->initialNumParticles() 
                % (1.0*constants()->initialNumParticles()/constants()->V()) 
                % tau % constants()->id());
    }

    /* Perform the rename for each file in the map */
    for (auto const& [key, filePtr] : file_)
    {

        std::string oldName(filePtr->getName());

        /* Replace with the new data name, we need to do this for both name and
         * backup name. */

		filePtr->replaceName(dataName);
		filePtr->replaceBakname(dataName);

        /* Perform the rename */
        fs::rename(oldName.c_str(), filePtr->getName().c_str());
    }

    /* Perform the rename for each file in the map */
    for (auto const& [key, filePtr] : hdf5File_)
    {

        std::string oldName(filePtr->name);

        /* Replace with the new data name, we need to do this for both name and
         * backup name. 
		 */
		filePtr->name = dataName + ".hdf5";
		filePtr->bakname = dataName + ".bak";

        /* Perform the rename */
        fs::rename(oldName.c_str(), filePtr->name.c_str());
    }
}

StateObject Communicator::createStateObject(const char *fileInitStr){
	StateObject newStateObject;

	// Create the full filename to check extension
	std::string filename = std::string(fileInitStr);
	
	// Check if the filename already has .hdf5 extension
	std::string hdf5Filename;
	if (filename.size() >= 5 && filename.substr(filename.size() - 5) == ".hdf5") {
		hdf5Filename = filename;  // Use as-is if already has .hdf5
	} else {
		hdf5Filename = filename + ".hdf5";  // Append .hdf5 if not present
	}
	
	// Check if HDF5 file exists
	if (fs::exists(hdf5Filename)) {
		std::cout << "Loading HDF5 state file: " << hdf5Filename << std::endl;
		try {
			// Create a temporary HDF5File object directly with the specified filename
			HDF5File tempFile(hdf5Filename);
			tempFile.open(OpenMode::READ);
			newStateObject = tempFile.readState();
			std::cout << "Successfully loaded HDF5 state" << std::endl;
		} catch (const std::exception& e) {
			std::cerr << "Error: Failed to load HDF5 state file " << hdf5Filename << ": " << e.what() << std::endl;
			std::exit(EXIT_FAILURE);
		}
	} else {
		std::cerr << "Error: State file not found: " << hdf5Filename << std::endl;
		std::exit(EXIT_FAILURE);
	}
		
	return newStateObject;
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
