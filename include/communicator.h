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
#include "hdf5.h"
#include "hdf5_data_object.h"
#include "state_object.h"

// ========================================================================  
// An enum class for abstracting the open modes.
// ========================================================================  
enum class OpenMode
{
	READ,
	READ_WRITE,
	WRITE,
	WRITE_APPEND,
	WRITE_TRUNC,
};

// ========================================================================  
// BaseFile Class
// ========================================================================  
/** 
 * An interface for all file related classes to inherit from
 *
 */
class BaseFile
{
    public:
		BaseFile() = default;
        BaseFile(std::string);
		virtual ~BaseFile() = default;
		virtual void open(OpenMode modeEnum) = 0;
		virtual void reset(){
			close();
			open(OpenMode::WRITE_TRUNC, bakname);
		}
		virtual void rename();
		virtual void prepare() = 0;
		virtual void write(std::string data) = 0;  // This needs to be overloaded in HDF5File
		virtual void flush() = 0;
		virtual bool prepared() {return prepared_;}
		virtual void close() = 0;	
		virtual bool exists() {return exists_;}

		virtual std::string getName() {return name;}
		virtual std::string getBakname() {return bakname;}
		virtual void setName(std::string newName){name = newName;}
		virtual void setBakname(std::string newName){bakname = newName;} 

		virtual void replaceName(std::string newName) = 0;
		virtual void replaceBakname(std::string newName) = 0;

        virtual std::fstream & stream() { 
            return rwfile;
        }

	protected:
        std::string name;        // The File name
        std::string bakname;     // The backup name
        std::fstream rwfile;     // The i/o file object

        bool exists_;       // Does the file exist? Check on creation.
        bool prepared_;      // Has the file already been prepared for writing?
							 //
        virtual void open(OpenMode modeEnum, std::string) = 0;

};


// ========================================================================  
// File Class
// ========================================================================  
/** 
 * A basic input/output file class.
 *
 */
class File : public BaseFile
{
    public:
    
        File(std::string, std::string, std::string, std::string);
        File(std::string);
        ~File() override {close();}

        /* Return the file stream */
        std::fstream & stream() override { 
            return rwfile;
        }

		void flush() override {stream().flush();}
        /* Open the primary file */
        void open(OpenMode modeEnum) override;

        /* Reset and rename the primary file */
        void reset() override;
        void rename() override;
        void prepare()  override {prepared_ = true;}
		void write(std::string data) override;

        /* Close the file if open */
        void close() override;       

		void replaceName(std::string newName) override {
        	name.replace(name.end()-newName.length()-4, name.end()-4, newName);
		}
		void replaceBakname(std::string newName) override{
        	bakname.replace(bakname.end()-newName.length()-4, bakname.end()-4, newName);
		}


    protected:
        friend class Communicator;    // Friends for I/O

        /* An alternate open which takes a filename */
        void open(OpenMode modeEnum, std::string) override;
	
	private:
		std::ios_base::openmode OpenModeTranslator(OpenMode modeEnum);

};
// ========================================================================  
// HDF5 File Class
// ========================================================================  
/** 
 * An input/output file class for HDF5 files.
 *
 */
class HDF5File : public BaseFile
{
	public:

		HDF5File(std::string, std::string, std::string, std::string);
		HDF5File(std::string _name);
		~HDF5File() override{
			close();
		}

		std::fstream & stream() override {
			std::cerr << "ERROR: stream() is not implemented for HDF5 files.\n";
			exit(EXIT_FAILURE);
			return rwfile;
		}
		void flush() override {
			std::cerr << "ERROR: flush() is not implemented for HDF5 files.\n";
			exit(EXIT_FAILURE);
		}

		void reset() override;
		void rename() override;

		/* Open the primary file */
		void open(OpenMode modeEnum) override;



		void write(std::string data) override {
			std::cerr << "ERROR: Writing data in string format is not implemented for HDF5 Files.\n";
			exit(EXIT_FAILURE);
		};

		template<typename T>
		void writeHDF5(
			const HDF5DataObject<T>& dataObj,
			const std::string parentPath = "/"
		);

		void writeState(StateObject state);
		StateObject readState();

		template<typename T>
		std::vector<T> readDataset(const char *datasetName);

        void prepare()  override {prepared_ = true;}
        /* Close the file if open */
        void close() override{
			H5Fclose(file_id);
		}

		void replaceName(std::string newName) override {
        	name.replace(name.end()-newName.length()-5, name.end()-5, newName);
		}
		void replaceBakname(std::string newName) override{
        	bakname.replace(bakname.end()-newName.length()-4, bakname.end()-4, newName);
		}
	protected:
        friend class Communicator;    // Friends for I/O

		hid_t openDataset(std::string datasetName){
			return H5Dopen2(file_id, datasetName.c_str(), H5P_DEFAULT);
		};

		void open(OpenMode modeEnum, std::string _name) override;

	private:
		hid_t file_id;  // HDF5 file handle
		unsigned OpenModeTranslator(OpenMode modeEnum);

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
			if (!file_.count(type)){
                initFile(type);
			}
            return &file_.at(type);
        }

        /** Get method returning hdf5 file object */
        HDF5File *hdf5File(std::string type) {
			if (!hdf5File_.count(type)){
                initHDF5File(type);
			}
            return &hdf5File_.at(type);
        }

        /** Get method returning initialization file name */
        std::string getInitName() const {
            return initName;
        }

        void updateNames();

		StateObject createStateObject(const char *fileInitStr);


    protected:
        Communicator() {}                                   
        Communicator(const Communicator&);              ///< Copy constructor
        Communicator& operator= (const Communicator&);  ///< Singleton equals

    private:
        // std::ios_base::openmode mode;    // The file i/o mode
		OpenMode mode;  	// The file i/o mode

        std::string ensemble;            // The type of ensemble
        std::string dataName;            // The labelling scheme of the output files
        std::string header;              // A unique file header

        std::string initName;      // A possible initial file name
        std::string fixedName;     // A posible fixed file name
        std::string baseDir;       // The output base directory

        double tau;          // A local copy of the actual imaginary time step.

        boost::ptr_map<std::string,File> file_; // The plain text file map
        boost::ptr_map<std::string,HDF5File> hdf5File_; // The hdf5 file map

        /* Initialize a input/output file */
        void initFile(std::string);
        void initHDF5File(std::string);
};


/**************************************************************************//**
 *  Global public access to the communcator singleton.
******************************************************************************/
inline Communicator* communicate() {
    Communicator *temp = Communicator::getInstance();
    return temp;
}
#endif

#include "communicator_impl.tpp"

