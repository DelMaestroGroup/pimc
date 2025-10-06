#ifndef STATE_OBJECT 
#define STATE_OBJECT

#include "hdf5_data_object.h"

struct StateObject {
	HDF5DataObject<uint32> worldOrder;
	HDF5DataObject<uint32> acceptanceInfo;
	HDF5DataObject<uint32> moveAcceptance;
	HDF5DataObject<uint32> estimatorSampling;
	HDF5DataObject<double> paths;
	HDF5DataObject<double> next;
	HDF5DataObject<double> prev;
	HDF5DataObject<uint32> worm;
	HDF5DataObject<uint32> randomState;

	StateObject() : worldOrder(), acceptanceInfo(), moveAcceptance(), estimatorSampling(), paths(), prev(), worm(), randomState() {}
};
#endif
