/**
 * @file xyz_io.h
 * @author Timothy R. Prisk
 * @date 2026-06-05
 *
 * @brief Read initial particle positions from a plain XYZ file.
 *
 * This header declares readXYZFile(), a utility for parsing standard
 * XYZ-format files and returning the contents as a DynamicArray of dVec
 * suitable for handing to the Path constructor as an initial condition.
 */
#ifndef XYZ_IO_H
#define XYZ_IO_H

#include "common.h"     // pulls in dVec, DynamicArray<>, NDIM, and <string>

class Container;

/**
 * @brief Peek the xyz file and get the atom count from the first line.
 */
int peekXYZAtomCount(const std::string &filename);

/**
 * @brief Read positions from a plain XYZ file and return them as a 
 * DynamicArray of dVec.
 *
 * @param filename  Path to the XYZ file.
 * @param boxPtr    Pointer to the simulation container.
 *
 * @return  A DynamicArray<dVec,1> of length N holding the positions in
 *          Angstroms. Indexed 0..N-1.
 *
 * @throws std::runtime_error if the file cannot be opened, is malformed,
 *         or contains coordinates outside the simulation cell.
 */
DynamicArray<dVec,1> readXYZFile(const std::string &filename,
                                  const Container *boxPtr);

#endif // XYZ_IO_H