/**
 * @file common.h
 * @author Adrian Del Maestro
 * @date 10.13.2008
 * @brief Global common header with shared dependencies and methods.
 * @details We use blitz++ for all arrays and used boost::format to simplify
 * input and output.
 * @see http://sourceforge.net/projects/blitz/
 * @see http://www.boost.org/doc/libs/1_34_1/libs/format/index.html 
 */

#ifndef COMMON_H 
#define COMMON_H

#include <iomanip>
#include <iostream>
#include <string>
#include <cmath>
#include <cassert>
#include <vector>
#include <map>

#include <blitz/array.h>
#include <boost/format.hpp>

/* Debugging librarys and definitions. There is only an effect if
 * the code is compiled with PIMC_DEBUG on.*/
#ifdef PIMC_DEBUG

/** Switch on the debug flag for Blitz++ */
#define BZ_DEBUG

/** A PIMC debugging message */
#define PIMC_DEBUG_MESSAGE(X)\
{ std::cerr << "[pimc] " << __FILE__ << ": " << __LINE__ << " " << X << std::endl; }

/** Rename assert method */
#define PIMC_ASSERT(X) assert(X)

#else // !PIMC_DEBUG

/** A PIMC debugging message */
#define PIMC_DEBUG_MESSAGE(X)

/** Rename assert method */
#define PIMC_ASSERT(X) 

#endif // ifdef PIMC_DEBUG 

/* We either specify the number of dimensions of the simulation at compile time,
 * or it defaults to 1D. */
#ifndef NDIM
#define NDIM 1 ///< Number of spatial dimnsions
#endif

/* We Default to turning on the NN lookup table. Comment this line out for testing
 * small systems where the lookup table might actually slow things down. */
#define NN_TABLE ///< Turn on the nearest neighbor particle lookup table

#include "MersenneTwister.h"

#define NPCFSEP 400 	///< The number of spatial separations to be used in the pair correlation function
#define NOBDMSEP 50	 	///< The number of spatial separations to be used in the one body density matrix
#define NRADSEP 200	 	///< The number of spatial separations to be used in the radial density
#define EPS 1.0E-7 		///< A small number
#define BIG 1.0E30		///< A big number
#define XXX	-999		///< Used to refer to a nonsense value

using namespace std;
using namespace blitz;
using boost::format;

/** Unsigned integer type, at least 32 bits */
typedef unsigned long uint32;  		

/** A NDIM-vector of type double */
typedef TinyVector<double,NDIM> dVec;

/** A NDIM-vector of type integer*/
typedef TinyVector<int,NDIM> iVec;

/** time-slice,bead-number world line index */
typedef TinyVector<int,2> beadLocator;

/** Integer array iterator */
typedef Array<int,1>::iterator intIter;

/** Constant integer array iterator */
typedef Array<int,1>::const_iterator cintIter;

/** beadLocator array iterator */
typedef Array<beadLocator,1>::iterator beadIter;

/** Each bead can have three possible states */
enum beadState {HEADTAIL,SPECIAL,NONE};

/** Return the integer value of a number raised to a power */
inline int ipow (int base, int power) {
	return static_cast<int>(floor(pow(1.0*base,1.0*power) + EPS));
}

/** Minimum of two inputs */
template<typename Ttype>
inline Ttype& min(Ttype& x, Ttype& y) { return (x < y ? x : y); } 

/** Maximum of two inputs */
template<typename Ttype>
inline Ttype& max(Ttype& x, Ttype& y) { return (x > y ? x : y); } 

/****************************************************************************/
/* Binary Search Functions
 * These functions provide templated binary search functionality which
 * returns an iterator to the position in an array where we would insert
 * a provided value.
 */
/*****************************************************************************/

/** The distance between two iterators */
template<class Titer>
inline void distanceIter(Titer first, Titer last, unsigned int &n)
{ while (first != last) { ++first; ++n;} }

/** Advance an iterator n steps */
template<class Titer>
void advance( Titer &i, unsigned int n) { while (n--) ++i; }

/** Return an iterator which points to the position where we 
 * would insert x into a sorted array */
template<class Titer, typename Ttype> 
Titer lowerBound(Titer first, Titer last, Ttype x) {
	unsigned int len = 0;
	unsigned int half;
	distanceIter(first, last, len);
	Titer middle;
	while (len > 0) {
		half = len >> 1;
		middle = first;
		advance(middle, half);
		if (*middle < x) {
			first = middle;
			++first;
			len = len - half - 1;
		}
		else
			len = half;
	}
	return first;
}

#endif

