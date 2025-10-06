// Implementations for HDF5File template methods moved here so they are
// available to all translation units that include communicator.h

#ifndef COMMUNICATOR_IMPL_TPP
#define COMMUNICATOR_IMPL_TPP

// Ensure this file is only included from C++ compilation units that have
// the required declarations

template<typename T>
std::vector<T> HDF5File::readDataset(const char *datasetName){
    bool openedHere = false;
    if (file_id == H5I_INVALID_HID) {
        // open file for read if not already open
        open(OpenMode::READ);
        openedHere = true;
    }

    hid_t dataset_id = H5Dopen2(file_id, datasetName, H5P_DEFAULT);

    if (dataset_id == H5I_INVALID_HID) {
        std::cerr << "Error: Unable to open dataset " << datasetName << " in file " << name << std::endl;
        if (openedHere) H5Fclose(file_id);
        exit(EXIT_FAILURE);
    }

    // Get dataspace and dimensions
    hid_t dataspace_id = H5Dget_space(dataset_id);
    int ndims = H5Sget_simple_extent_ndims(dataspace_id);
    std::vector<hsize_t> dims(ndims);
    H5Sget_simple_extent_dims(dataspace_id, dims.data(), NULL);

    hsize_t numElements = 1;
    for (int i = 0; i < ndims; ++i) numElements *= dims[i];

    // Allocate a vector for the data
    std::vector<T> buffer;
    try {
        buffer.resize(static_cast<size_t>(numElements));
    } catch (const std::bad_alloc &e) {
        std::cerr << "Unable to allocate memory for dataset " << datasetName << ": " << e.what() << std::endl;
        H5Sclose(dataspace_id);
        H5Dclose(dataset_id);
        if (openedHere) H5Fclose(file_id);
        exit(EXIT_FAILURE);
    }

    // Determine native memory datatype from the dataset type
    hid_t dtype = H5Dget_type(dataset_id);
    hid_t mem_type = H5Tget_native_type(dtype, H5T_DIR_DEFAULT);

    herr_t status = H5Dread(dataset_id, mem_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer.data());

    if (status < 0) {
        std::cerr << "Error: Unable to read dataset " << datasetName << " in file " << name << std::endl;
        H5Tclose(mem_type);
        H5Tclose(dtype);
        H5Sclose(dataspace_id);
        H5Dclose(dataset_id);
        if (openedHere) H5Fclose(file_id);
        exit(EXIT_FAILURE);
    }

    // Clean up HDF5 types and handles
    H5Tclose(mem_type);
    H5Tclose(dtype);
    H5Sclose(dataspace_id);
    H5Dclose(dataset_id);
    if (openedHere) H5Fclose(file_id);

    return buffer;
}


template<typename T>
void HDF5File::writeHDF5(
    const HDF5DataObject<T>& dataObject,
    const std::string parentPath
){
    std::string hdf5Key = parentPath + "/" + dataObject.name;

    // Create a dataspace for the dataset
    hid_t dataspace_id = H5Screate_simple(static_cast<int>(dataObject.dims.size()), dataObject.dims.data(), NULL);

    // Determine the HDF5 file datatype (storage) and the memory datatype
    // used for the call to H5Dwrite. For UINT_32 we always store as a
    // 32-bit standard little-endian type for portability, and choose the
    // memory type from the template parameter T.
    hid_t storage_type;
    hid_t mem_type;
    switch (dataObject.type) {
        case DataType::UINT_32:
            storage_type = H5T_STD_U32LE; // stable 32-bit storage
            // pick mem type based on size of T
            if (sizeof(T) == 4) mem_type = H5T_NATIVE_UINT32;
            else if (sizeof(T) == 8) mem_type = H5T_NATIVE_UINT64;
            else mem_type = H5T_NATIVE_UINT;
            break;
        case DataType::FLOAT:
            storage_type = H5T_IEEE_F32LE;
            mem_type = H5T_NATIVE_FLOAT;
            break;
        case DataType::DOUBLE:
            storage_type = H5T_IEEE_F64LE;
            mem_type = H5T_NATIVE_DOUBLE;
            break;
        case DataType::CHAR:
            storage_type = H5T_C_S1;
            mem_type = H5T_NATIVE_CHAR;
            break;
        case DataType::UNSIGNEDLONG:
            storage_type = H5T_STD_U64LE;
            mem_type = H5T_NATIVE_ULONG;
            break;
        case DataType::UNSIGNEDINT:
            storage_type = H5T_STD_U32LE;
            mem_type = H5T_NATIVE_UINT;
            break;
        default:
            std::cerr << "Unknown data type: " << static_cast<std::underlying_type<DataType>::type>(dataObject.type) << std::endl;
            H5Sclose(dataspace_id);
            exit(EXIT_FAILURE);
    }

    hid_t dataset_id;
    htri_t exists = H5Lexists(file_id, hdf5Key.c_str(), H5P_DEFAULT);
    if (exists > 0){
        dataset_id = H5Dopen2(file_id, hdf5Key.c_str(), H5P_DEFAULT);
        if (dataset_id == H5I_INVALID_HID){
            std::cerr << "Unable to open dataset: " << hdf5Key.c_str() << std::endl;
            H5Sclose(dataspace_id);
            exit(EXIT_FAILURE);
        }

        // Check on-disk datatype for compatibility with requested storage_type.
        // If it doesn't match (for example stored as float32 but we want float64)
        // then delete and re-create the dataset to avoid lossy writes.
        hid_t existing_dtype = H5Dget_type(dataset_id);
        bool compatible = false;
        if (existing_dtype != H5I_INVALID_HID) {
            // Compare datatype class and size where appropriate
            H5T_class_t existing_class = H5Tget_class(existing_dtype);
            H5T_class_t desired_class = H5Tget_class(storage_type);
            size_t existing_size = H5Tget_size(existing_dtype);
            size_t desired_size = H5Tget_size(storage_type);
            if (existing_class == desired_class && existing_size == desired_size) compatible = true;
            H5Tclose(existing_dtype);
        }

        if (!compatible) {
            // Close and remove the existing dataset then create a new one with correct storage
            H5Dclose(dataset_id);
            // Attempt to remove the existing link/dataset; ignore failure but report it
            if (H5Lexists(file_id, hdf5Key.c_str(), H5P_DEFAULT) > 0) {
                if (H5Ldelete(file_id, hdf5Key.c_str(), H5P_DEFAULT) < 0) {
                    std::cerr << "Warning: failed to delete incompatible dataset " << hdf5Key << ". Attempting to continue." << std::endl;
                }
            }

            dataset_id = H5Dcreate2(
                file_id,
                hdf5Key.c_str(),
                storage_type,
                dataspace_id,
                H5P_DEFAULT,
                H5P_DEFAULT,
                H5P_DEFAULT
            );
            if (dataset_id == H5I_INVALID_HID){
                std::cerr << "Unable to create dataset: " << hdf5Key.c_str() << std::endl;
                H5Sclose(dataspace_id);
                exit(EXIT_FAILURE);
            }
        }

        // Write the data to the dataset
        herr_t writeStatus = H5Dwrite(dataset_id, mem_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, dataObject.data.data());
        if (writeStatus < 0){
            std::cerr << "Unable to write dataset: " << hdf5Key.c_str() << std::endl;
            H5Dclose(dataset_id);
            H5Sclose(dataspace_id);
            exit(EXIT_FAILURE);
        }

    } else if (exists == 0){
        // Create the dataset
        dataset_id = H5Dcreate2(
            file_id,
            hdf5Key.c_str(),
            storage_type,
            dataspace_id,
            H5P_DEFAULT,
            H5P_DEFAULT,
            H5P_DEFAULT
        );
        if (dataset_id == H5I_INVALID_HID){
            std::cerr << "Unable to create dataset: " << hdf5Key.c_str() << std::endl;
            H5Sclose(dataspace_id);
            exit(EXIT_FAILURE);
        }

        // Write the data to the dataset
    herr_t writeStatus = H5Dwrite(dataset_id, mem_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, dataObject.data.data());
        if (writeStatus < 0){
            std::cerr << "Unable to write dataset: " << hdf5Key.c_str() << std::endl;
            H5Dclose(dataset_id);
            H5Sclose(dataspace_id);
            exit(EXIT_FAILURE);
        }
    } else {
        std::cerr << "Error while checking dataset: " << hdf5Key.c_str() << std::endl;
        H5Sclose(dataspace_id);
        exit(EXIT_FAILURE);
    }

    // Clean up resources
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);
}

#endif // COMMUNICATOR_IMPL_TPP
