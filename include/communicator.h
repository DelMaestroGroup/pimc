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
};


/**************************************************************************//**
 *  Global public access to the communcator singleton.
******************************************************************************/
inline Communicator* communicate() {
    Communicator *temp = Communicator::getInstance();
    return temp;
}
#endif

