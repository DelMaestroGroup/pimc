// convert_blitz_to_dynamic.cpp
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <array>
#include <algorithm>
#include <string>
#include <functional>
#include <boost/archive/binary_iarchive.hpp>
#include "blitz/array.h"
#include "dynamic_array.h"
#include "dynamic_array_serialization.h"

// Helper function to convert a blitz::Array to a DynamicArray.
template<typename T, int Rank>
DynamicArray<T, Rank> convertBlitzToDynamic(const blitz::Array<T, Rank>& blitzArray) {
    // Get the extents from blitz::Array.
    blitz::TinyVector<int, Rank> shape = blitzArray.shape();
    std::array<std::size_t, Rank> extents;
    for (int i = 0; i < Rank; ++i)
        extents[i] = static_cast<std::size_t>(shape(i));

    // Create and resize the DynamicArray.
    DynamicArray<T, Rank> dynArray;
    std::apply([&](auto&&... args){ dynArray.resize(args...); }, extents);

    // Copy the data (assuming blitz stores data contiguously).
    std::copy(blitzArray.data(), blitzArray.data() + dynArray.size(), dynArray.data());
    return dynArray;
}

void printUsage(const std::string& progName) {
    std::cerr << "Usage: " << progName << " <file_prefix>\n"
              << "  Looks for a legacy archive file named <file_prefix>_serialized.dat\n"
              << "  and writes new DynamicArray archives as <file_prefix>_serialized_V3d.dat, etc.\n";
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        printUsage(argv[0]);
        return 1;
    }

    std::string filePrefix = argv[1];
    std::string legacyFile = filePrefix + "_serialized.dat";

    try {
        // Open the legacy archive file.
        std::ifstream ifs(legacyFile, std::ios::binary);
        if (!ifs)
            throw std::runtime_error("Failed to open legacy archive file: " + legacyFile);

        // Load legacy data using Boost.Archive.
        const auto aflags = boost::archive::no_header | boost::archive::no_tracking;
        boost::archive::binary_iarchive ia(ifs, aflags);

        blitz::Array<double, 3> V3d;
        blitz::Array<double, 3> gradV3d_x;
        blitz::Array<double, 3> gradV3d_y;
        blitz::Array<double, 3> gradV3d_z;
        blitz::Array<double, 3> grad2V3d;
        blitz::Array<double, 1> LUTinfo;
        ia >> V3d >> gradV3d_x >> gradV3d_y >> gradV3d_z >> grad2V3d >> LUTinfo;

        // Convert each blitz::Array to a DynamicArray.
        auto dynV3d       = convertBlitzToDynamic(V3d);
        auto dynGradV3d_x = convertBlitzToDynamic(gradV3d_x);
        auto dynGradV3d_y = convertBlitzToDynamic(gradV3d_y);
        auto dynGradV3d_z = convertBlitzToDynamic(gradV3d_z);
        auto dynGrad2V3d  = convertBlitzToDynamic(grad2V3d);
        auto dynLUTinfo   = convertBlitzToDynamic(LUTinfo);

        // Save the new DynamicArrays using the new serialization functions.
        saveDynamicArray(dynV3d,       filePrefix + "_serialized_V3d.dat");
        saveDynamicArray(dynGradV3d_x, filePrefix + "_serialized_gradV3d_x.dat");
        saveDynamicArray(dynGradV3d_y, filePrefix + "_serialized_gradV3d_y.dat");
        saveDynamicArray(dynGradV3d_z, filePrefix + "_serialized_gradV3d_z.dat");
        saveDynamicArray(dynGrad2V3d,  filePrefix + "_serialized_grad2V3d.dat");
        saveDynamicArray(dynLUTinfo,   filePrefix + "_serialized_LUTinfo.dat");

        std::cout << "Conversion completed successfully.\n";
    }
    catch (const std::exception& ex) {
        std::cerr << "Conversion failed: " << ex.what() << "\n";
        return 1;
    }

    return 0;
}

