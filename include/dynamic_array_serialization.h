#ifndef DYNAMIC_ARRAY_SERIALIZATION_H
#define DYNAMIC_ARRAY_SERIALIZATION_H

#include "dynamic_array.h"
#include <fstream>
#include <stdexcept>
#include <array>
#include <numeric>
#include <functional>

// Save the DynamicArray to a binary file.
template <typename T, std::size_t Rank>
void saveDynamicArray(const DynamicArray<T, Rank>& arr, const std::string& filename) {
    std::ofstream ofs(filename, std::ios::binary);
    if (!ofs)
        throw std::runtime_error("Could not open file for writing: " + filename);
    
    // Write the extents (dimensions)
    auto extents = arr.extents();
    ofs.write(reinterpret_cast<const char*>(extents.data()), extents.size() * sizeof(std::size_t));
    
    // Write the total number of elements.
    std::size_t total = arr.size();
    ofs.write(reinterpret_cast<const char*>(&total), sizeof(total));
    
    // Write the raw data.
    ofs.write(reinterpret_cast<const char*>(arr.data()), total * sizeof(T));
    
    if (!ofs)
        throw std::runtime_error("Error occurred while writing data to " + filename);
}

// Load the DynamicArray from a binary file.
template <typename T, std::size_t Rank>
void loadDynamicArray(DynamicArray<T, Rank>& arr, const std::string& filename) {
    std::ifstream ifs(filename, std::ios::binary);
    if (!ifs)
        throw std::runtime_error("Could not open file for reading: " + filename);
    
    // Read extents.
    std::array<std::size_t, Rank> extents;
    ifs.read(reinterpret_cast<char*>(extents.data()), extents.size() * sizeof(std::size_t));
    if (!ifs)
        throw std::runtime_error("Error reading extents from " + filename);
    
    // Resize the array to match the extents.
    std::apply([&](auto&&... args){ arr.resize(args...); }, extents);
    
    // Read total number of elements.
    std::size_t total;
    ifs.read(reinterpret_cast<char*>(&total), sizeof(total));
    if (!ifs)
        throw std::runtime_error("Error reading total size from " + filename);
    
    if (total != arr.size())
        throw std::runtime_error("Data size mismatch in file " + filename);
    
    // Read the raw data.
    ifs.read(reinterpret_cast<char*>(arr.data()), total * sizeof(T));
    if (!ifs)
        throw std::runtime_error("Error reading data from " + filename);
}

#endif // DYNAMIC_ARRAY_SERIALIZATION_H

