#ifndef HDF5_DATA_OBJECT 
#define HDF5_DATA_OBJECT

#include <hdf5.h>
#include <string>
#include <vector>


// Enum to hold the type of data
enum class DataType {
    UINT_32,
    FLOAT,
    DOUBLE,
    CHAR,
	UNSIGNEDLONG,
	UNSIGNEDINT
};

// Struct to hold information and data for a dataset to be written to HDF5
template<typename T>
struct HDF5DataObject {
    std::string name;             // Dataset name
    DataType type;                // Data type (e.g., int, float)
    std::vector<hsize_t> dims;    // Dimensions of the dataset
	std::vector<T> data;                   // Pointer to the data (generic)

	// I should probably use the constructor to get the datatype variable automatically via the T type
    HDF5DataObject(
		const std::string& dataset_name,
		DataType datatype,
		const std::vector<hsize_t>& dimensions,
		const std::vector<T>& dataset
	) : name(dataset_name), type(datatype), dims(dimensions), data(dataset) {}

	HDF5DataObject() : name(""), type(), dims(0, 0), data(0, T()) {}

};
#endif
