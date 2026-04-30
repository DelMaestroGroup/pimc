/**
 * @file communicator.h
 * @author Adrian Del Maestro
 * @date 03.25.2009 
 *
 * @brief Class definitions for all file input/output
 */

#ifndef COMMUNICATOR_H 
#define COMMUNICATOR_H

#include "common.h"
#include "constants.h"
#include <cstring>
#include <fstream>

#ifdef ENABLE_HDF5
#include <hdf5.h>
#include <vector>
#endif

// ========================================================================  
// File Class
// ========================================================================  
/** 
 * A basic input/output file class.
 *
 */
class File
{
    public:
    
        File(std::string, std::string, std::string, std::string);
        File(std::string);
        ~File() {close();}

        /* Return the file stream */
        std::fstream & stream() { 
            return rwfile;
        }

        /* Open the primary file */
        void open(std::ios_base::openmode);

        /* Reset and rename the primary file */
        void reset();
        void rename();
        void prepare() {prepared_ = true;}
        bool prepared() {return prepared_;}

        /* Close the file if open */
        void close();       

        bool exists() {return exists_;}    ///< did the file exist before opening?

    protected:
        friend class Communicator;    // Friends for I/O

        std::string name;        // The File name
        std::string bakname;     // The backup name

        bool exists_;       // Does the file exist? Check on creation.
        bool prepared_;      // Has the file already been prepared for writing?

        std::fstream rwfile;     // The i/o file object

        /* An alternate open which takes a filename */
        void open(std::ios_base::openmode, std::string);

};

#ifdef ENABLE_HDF5
class H5File
{
    public:

        H5File(std::string _type, std::string _data,
               std::string ensemble, std::string outDir);
        H5File(std::string _name);
        ~H5File() { close(); }

        /* Open with one of: H5F_ACC_RDONLY, H5F_ACC_RDWR, H5F_ACC_TRUNC */
        void open(unsigned mode);

        void reset();
        void rename();

        void close();

        hid_t id() const { return file_id; }
        bool  isOpen() const { return file_id >= 0; }
        bool  exists() const { return exists_; }


        template <typename T>
        void writeScalar(const std::string& path, T value);

        template <typename T>
        T readScalar(const std::string& path);

        template <typename T>
        void writeArray(const std::string& path,
                        const T* data,
                        const std::vector<hsize_t>& dims);

        template <typename T>
        void readArray(const std::string& path, T* data);

        std::vector<hsize_t> getDims(const std::string& path);

        void createGroup(const std::string& path);
        bool hasObject(const std::string& path);

    protected:
        friend class Communicator;

        std::string name;
        std::string bakname;

        bool   exists_  = false;
        hid_t  file_id  = -1;

    private:
        /* HDF5 type map — specialized in the .cpp */
        template <typename T> static hid_t h5type();
};
#endif

// ========================================================================  
// Communicator Class
// ========================================================================  
/** 
 * Performs input/output.
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

        /** Initialize the output files */
        void init(double, bool, std::string, std::string);

        /** Get method returning file object */
        File *file(std::string type) {
            if (!file_.count(type))
                initFile(type);
            return &file_.at(type);
        }

#ifdef ENABLE_HDF5
        // HDF5 accessor 
        H5File *h5file(std::string type) {
            if (!h5file_.count(type))
                initH5File(type);
            return &h5file_.at(type);
        }
#endif

        void updateNames();

    protected:
        Communicator() {}                                   
        Communicator(const Communicator&);              ///< Copy constructor
        Communicator& operator= (const Communicator&);  ///< Singleton equals

    private:
        std::ios_base::openmode mode;    // The file i/o mode

        std::string ensemble;            // The type of ensemble
        std::string dataName;            // The labelling scheme of the output files
        std::string header;              // A unique file header

        std::string initName;      // A possible initial file name
        std::string fixedName;     // A posible fixed file name
        std::string baseDir;       // The output base directory

        double tau;          // A local copy of the actual imaginary time step.

        boost::ptr_map<std::string,File> file_; // The file map

        /* Initialize a input/output file */
        void initFile(std::string);

#ifdef ENABLE_HDF5
        void initH5File(std::string);
        boost::ptr_map<std::string, H5File> h5file_;
#endif
};


/**************************************************************************//**
 *  Global public access to the communcator singleton.
******************************************************************************/
inline Communicator* communicate() {
    Communicator *temp = Communicator::getInstance();
    return temp;
}
#endif

