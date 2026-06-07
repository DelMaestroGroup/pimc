/**
 * @file xyz_io.cpp
 * @author Timothy R. Prisk
 * @date 2026-06-05
 *
 * @brief Implementation of readXYZFile().
 */
#include "xyz_io.h"
#include "container.h"     // Needed to access boxPtr->side[d]

#include <fstream>         // std::ifstream for reading the file
#include <sstream>         // std::istringstream for per-line parsing
#include <stdexcept>       // std::runtime_error for signaling failures
#include <string>          // std::string for filename and line buffers


// Peek into the XYZ file and get the atom count from the first line.
int peekXYZAtomCount(const std::string &filename)
{
    std::ifstream fin(filename);
    if (!fin) 
    {
        std::ostringstream msg;
        msg << "peekXYZAtomCount: cannot open '" << filename << "'.";
        throw std::runtime_error(msg.str());
    }

    std::string line;
    if (!std::getline(fin, line))
        throw std::runtime_error(
            "peekXYZAtomCount: empty file (no atom count).");

    int nAtoms = 0;
    std::istringstream iss(line);
    if (!(iss >> nAtoms) || nAtoms <= 0) 
    {
        std::ostringstream msg;
        msg << "peekXYZAtomCount: first line of '" << filename
            << "' must be a positive integer atom count; got: '"
            << line << "'.";
        throw std::runtime_error(msg.str());
    }

    return nAtoms;
}


// Tolerence check. Is this atom inside the simulation box?
static constexpr double XYZ_BOUNDS_TOL = 1.0e-6;

// Read the XYZ file and return the positions as a DynamicArray<dVec,1>.
DynamicArray<dVec,1> readXYZFile(const std::string &filename,
                                  const Container *boxPtr)
{
    // Open the file.
    std::ifstream fin(filename);
    if (!fin) 
    {
        std::ostringstream msg;
        msg << "readXYZFile: cannot open '" << filename << "'.";
        throw std::runtime_error(msg.str());
    }

    // Get the number of atoms from the first line.
    int nAtoms = 0;
    {
        std::string line;
        if (!std::getline(fin, line))
            throw std::runtime_error("readXYZFile: empty file (no atom count).");
        std::istringstream iss(line);
        if (!(iss >> nAtoms) || nAtoms <= 0) 
        {
            std::ostringstream msg;
            msg << "readXYZFile: first line of '" << filename
                << "' must be a positive integer atom count; got: '"
                << line << "'.";
            throw std::runtime_error(msg.str());
        }
    }

    // Ignore the second line (comment).
    {
        std::string comment;
        if (!std::getline(fin, comment))
            throw std::runtime_error(
                "readXYZFile: file ended before comment line.");
    }

    // Read N atoms from the remaining lines. Each line has format:
    //   <element> <x> <y> <z>
    // where <element> is read into a string and thrown away, and the three
    // coordinates are parsed as doubles in Angstroms.
    DynamicArray<dVec,1> positions(nAtoms);

    for (int i = 0; i < nAtoms; ++i) 
    {
        std::string line;
        if (!std::getline(fin, line)) 
        {
            std::ostringstream msg;
            msg << "readXYZFile: file '" << filename
                << "' ended after " << i << " atoms; expected " << nAtoms << ".";
            throw std::runtime_error(msg.str());
        }

        // Parse element + NDIM coordinates from the line.
        std::istringstream iss(line);
        std::string element;
        dVec pos;

        if (!(iss >> element)) 
        {
            std::ostringstream msg;
            msg << "readXYZFile: could not parse element column for atom "
                << i << " in '" << filename << "'. Line: '" << line << "'.";
            throw std::runtime_error(msg.str());
        }

        for (int d = 0; d < NDIM; ++d) 
        {
            if (!(iss >> pos[d])) {
                std::ostringstream msg;
                msg << "readXYZFile: could not parse coordinate " << d
                    << " for atom " << i << " in '" << filename
                    << "'. Line: '" << line << "'.";
                throw std::runtime_error(msg.str());
            }
        }

        /* Bounds check: the simulation cell is centered at the origin, so
         * each coordinate must lie in [-side[d]/2, +side[d]/2], within the
         * XYZ_BOUNDS_TOL tolerance. The user is responsible for ensuring
         * the box (set via -L or computed from -n) is large enough to
         * contain all atoms. */
        for (int d = 0; d < NDIM; ++d) 
        {
            const double half = 0.5 * boxPtr->side[d];
            if (pos[d] < -half - XYZ_BOUNDS_TOL ||
                pos[d] >  half + XYZ_BOUNDS_TOL) 
                {
                std::ostringstream msg;
                msg << "readXYZFile: atom " << i
                    << " is outside the simulation cell in dimension " << d
                    << ". Coordinate = " << pos[d]
                    << " A, box half-extent = " << half << " A. "
                    << "Check that -L (or -n and -N) define a cell large enough "
                    << "to contain all atoms in '" << filename << "'.";
                throw std::runtime_error(msg.str());
            }
        }

        positions(i) = pos;
    }

    return positions;
}