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
    
        File(string, string, string, string);
        File(string);
        ~File() {close();}

        /* Return the file stream */
		fstream & stream() { 
            return rwfile;
        }

        /* Open the primary file */
		void open(ios_base::openmode);

        /* Reset and rename the primary file */
        void reset();
        void rename();

        /* Close the file if open */
        void close();       

    protected:
        string name;        // The File name
        string bakname;     // The backup name

        fstream rwfile;     // The i/o file object

        /* An alternate open which takes a filename */
		void open(ios_base::openmode,string);

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
		void init(double,bool,string,string);

        /** Get method returning file object */
        File *file(string type) {
            if (!file_.count(type))
                initFile(type);
            return &file_.at(type);
        }

	protected:
		Communicator() {}									
		Communicator(const Communicator&);				///< Copy constructor
		Communicator& operator= (const Communicator&);	///< Singleton equals

	private:
        ios_base::openmode mode;    // The file i/o mode

		string ensemble;			// The type of ensemble
		string dataName;			// The labelling scheme of the output files
        string header;              // A unique file header

        string initName;      // A possible initial file name
        string fixedName;     // A posible fixed file name
        string baseDir;       // The output base directory

        boost::ptr_map<string,File> file_; // The file map

        /* Makes sure we have a unique PIMCID */
        void getUniqueID(const double);

        /* Initialize a input/output file */
        void initFile(string);
};


/**************************************************************************//**
 *  Global public access to the communcator singleton.
******************************************************************************/
inline Communicator* communicate() {
	Communicator *temp = Communicator::getInstance();
	return temp;
}
#endif

