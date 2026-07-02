/** 
 * @file communicator.cpp
 * @author Adrian Del Maestro
 *
 * @brief Communicator class implementation.
 */

#include "communicator.h"

/* Filesystem is tricky as it is not yet widely supported.  We try to address
 * that here. */
#if defined(USE_STD_FILESYSTEM)
    #include <filesystem>
    namespace fs = std::filesystem;
#elif defined(USE_EXPERIMENTAL_FILESYSTEM)
    #include <experimental/filesystem>
    namespace fs = std::experimental::filesystem;
#else
    #error "No filesystem support found"
#endif

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
    name = str(format("%s/%s-%s-%s.dat") % outDir % ensemble % _type % _data);

    /* Create a backup name */
    bakname = str(format("%s/%s-%s-%s.bak") % outDir % ensemble % _type % _data);

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
File::File(std::string _name) : name(_name), bakname() {

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
void File::open(std::ios_base::openmode mode) {

    /* Convert the filename to a c std::string, and open the file */ 
    rwfile.open(name.c_str(), mode);
    if (!rwfile) {
        std::cerr << "Unable to process file: " << name << std::endl;
        exit(EXIT_FAILURE);
    }
}

/**************************************************************************//**
 *  Open the file.
 *
 *  @param mode A valid file rw mode
 *  @param _name A valid file name
******************************************************************************/
void File::open(std::ios_base::openmode mode, std::string _name) {

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
    open(std::ios::out|std::ios::trunc,bakname);

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
    if (constants<std::string>("external").find("tube") != std::string::npos) {
        fs::path cylPath(baseDir + "/CYLINDER");
        fs::create_directory(cylPath);
    }

    /* A header line for the files */
    header =  str(format("# PIMCID: %s\n") % constants()->id());

    /* Depending on whether or not we are restarting the simulations, the open mode
     * changes. */
    if (constants()->restart()) {
        mode = std::ios::out|std::ios::app;
    }
    else {
        mode = std::ios::out;
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
         * restarting, or starting from a given initialization file */
        if (constants()->restart())
            file_.insert(type, new File(stateName,dataName,ensemble,baseDir));
        else 
            file_.insert(type, new File(initName));
        
        file_.at(type).open(std::ios::in);
    }
    /* Initialize a possible fixed coordinate file */
    else if (type == "fixed") {
        file_.insert(type,new File(fixedName));
        file_.at(type).open(std::ios::in);
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

        /* Construct the file and open it */
        file_.insert(type, new File(ctype,dataName,ensemble,outDir));
        file_.at(type).open(mode);

        /* Write the header line if the file doesn't exist */
        if (!file_.at(type).exists())
            file_.at(type).stream() << header;
    }
}

/**************************************************************************//**
 * Update the data name and rename any existing files
******************************************************************************/
void Communicator::updateNames() {

    /* We create a new dataName based on the posibility of updated paramters. */

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

        /* We don't want to update init or fixed files */
        if ( (key.find("init") == std::string::npos)  && 
             (key.find("fixed") == std::string::npos) )  {

            std::string oldName(filePtr->name);

            /* Replace with the new data name, we need to do this for both name and
             * backup name. */
            filePtr->name.replace(filePtr->name.end()-dataName.length()-4,filePtr->name.end()-4,dataName);
            filePtr->bakname.replace(filePtr->bakname.end()-dataName.length()-4,filePtr->bakname.end()-4,dataName);

            /* Perform the rename */
            fs::rename(oldName.c_str(), filePtr->name.c_str());
        }
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

#ifdef ENABLE_HDF5
void Communicator::initH5File(std::string type) {

    if (type.find("init") != std::string::npos) {
        std::string stateName = "state";
        if (type != "init") stateName += type.substr(4);
        if (constants()->restart())
            h5file_.insert(type, new H5File(stateName, dataName, ensemble, baseDir));
        else
            h5file_.insert(type, new H5File(initName));   // user-provided .h5 file
        h5file_.at(type).open(H5F_ACC_RDONLY);
    }
    else {
        std::string outDir = baseDir;
        std::string ctype  = type;
        if (type.find("cyl_") != std::string::npos) {
            outDir = baseDir + "/CYLINDER";
            ctype.erase(0, 4);
        }
        h5file_.insert(type, new H5File(ctype, dataName, ensemble, outDir));
        /* note: don't open here — saveState() calls reset() to prep the bakfile */
    }
}


// ---------------------------------------------------------------------------
// H5File CLASS --------------------------------------------------------------
// ---------------------------------------------------------------------------

/* Type mapping */
template <> hid_t H5File::h5type<int>()           { return H5T_NATIVE_INT; }
template <> hid_t H5File::h5type<unsigned int>()  { return H5T_NATIVE_UINT; }
template <> hid_t H5File::h5type<long>()          { return H5T_NATIVE_LONG; }
template <> hid_t H5File::h5type<unsigned long>() { return H5T_NATIVE_ULONG; }   
template <> hid_t H5File::h5type<float>()         { return H5T_NATIVE_FLOAT; }
template <> hid_t H5File::h5type<double>()        { return H5T_NATIVE_DOUBLE; }

/**************************************************************************//**
 *  Constructor with type/data/ensemble/outDir, mirroring File.
******************************************************************************/
H5File::H5File(std::string _type, std::string _data,
               std::string ensemble, std::string outDir)
{
    name    = str(format("%s/%s-%s-%s.h5")     % outDir % ensemble % _type % _data);
    bakname = str(format("%s/%s-%s-%s.bak.h5") % outDir % ensemble % _type % _data);
    exists_ = fs::exists(name);
}

H5File::H5File(std::string _name) : name(_name), bakname() {
    exists_ = fs::exists(name);
}

void H5File::close() {
    if (file_id >= 0) {
        H5Fclose(file_id);
        file_id = -1;
    }
}

void H5File::open(unsigned mode) {
    close();
    if (mode == H5F_ACC_TRUNC) {
        file_id = H5Fcreate(name.c_str(), H5F_ACC_TRUNC,
                            H5P_DEFAULT, H5P_DEFAULT);
    } else {
        file_id = H5Fopen(name.c_str(), mode, H5P_DEFAULT);
    }
    if (file_id < 0) {
        std::cerr << "Unable to open HDF5 file: " << name << std::endl;
        exit(EXIT_FAILURE);
    }
}

/* Open the .bak.h5 truncated, ready for a fresh write. */
void H5File::reset() {
    close();
    file_id = H5Fcreate(bakname.c_str(), H5F_ACC_TRUNC,
                        H5P_DEFAULT, H5P_DEFAULT);
    if (file_id < 0) {
        std::cerr << "Unable to create HDF5 backup: " << bakname << std::endl;
        exit(EXIT_FAILURE);
    }
}

void H5File::rename() {
    close();
    fs::rename(bakname.c_str(), name.c_str());
}

void H5File::createGroup(const std::string& path) {
    hid_t gid = H5Gcreate(file_id, path.c_str(),
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (gid >= 0) H5Gclose(gid);
}

bool H5File::hasObject(const std::string& path) {
    return H5Lexists(file_id, path.c_str(), H5P_DEFAULT) > 0;
}

/* Scalar write/read */
template <typename T>
void H5File::writeScalar(const std::string& path, T value) {
    hid_t space = H5Screate(H5S_SCALAR);
    hid_t dset  = H5Dcreate(file_id, path.c_str(), h5type<T>(), space,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dset, h5type<T>(), H5S_ALL, H5S_ALL, H5P_DEFAULT, &value);
    H5Dclose(dset);
    H5Sclose(space);
}

template <typename T>
T H5File::readScalar(const std::string& path) {
    T value{};
    hid_t dset = H5Dopen(file_id, path.c_str(), H5P_DEFAULT);
    H5Dread(dset, h5type<T>(), H5S_ALL, H5S_ALL, H5P_DEFAULT, &value);
    H5Dclose(dset);
    return value;
}

/* Array write/read — caller supplies dims and a flat buffer */
template <typename T>
void H5File::writeArray(const std::string& path,
                        const T* data,
                        const std::vector<hsize_t>& dims)
{
    hid_t space = H5Screate_simple(static_cast<int>(dims.size()),
                                   dims.data(), nullptr);
    hid_t dset  = H5Dcreate(file_id, path.c_str(), h5type<T>(), space,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dset, h5type<T>(), H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    H5Dclose(dset);
    H5Sclose(space);
}

template <typename T>
void H5File::readArray(const std::string& path, T* data) {
    hid_t dset = H5Dopen(file_id, path.c_str(), H5P_DEFAULT);
    H5Dread(dset, h5type<T>(), H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    H5Dclose(dset);
}

std::vector<hsize_t> H5File::getDims(const std::string& path) {
    hid_t dset  = H5Dopen(file_id, path.c_str(), H5P_DEFAULT);
    hid_t space = H5Dget_space(dset);
    int rank = H5Sget_simple_extent_ndims(space);
    std::vector<hsize_t> dims(rank);
    H5Sget_simple_extent_dims(space, dims.data(), nullptr);
    H5Sclose(space);
    H5Dclose(dset);
    return dims;
}

/* Explicit instantiations for the types used in saveState/loadState */
template void   H5File::writeScalar<int>(const std::string&, int);
template void   H5File::writeScalar<unsigned int>(const std::string&, unsigned int);
template void   H5File::writeScalar<double>(const std::string&, double);
template int    H5File::readScalar<int>(const std::string&);
template double H5File::readScalar<double>(const std::string&);
template void   H5File::writeArray<int>(const std::string&, const int*, const std::vector<hsize_t>&);
template void   H5File::writeArray<unsigned int>(const std::string&, const unsigned int*, const std::vector<hsize_t>&);
template void   H5File::writeArray<double>(const std::string&, const double*, const std::vector<hsize_t>&);
template void   H5File::readArray<int>(const std::string&, int*);
template void   H5File::readArray<unsigned int>(const std::string&, unsigned int*);
template void   H5File::readArray<double>(const std::string&, double*);
template void   H5File::writeArray<unsigned long>(const std::string&, const unsigned long*, const std::vector<hsize_t>&);
template void   H5File::readArray<unsigned long>(const std::string&, unsigned long*);
#endif
