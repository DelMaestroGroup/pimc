/**
 * @file common.h
 * @author Adrian Del Maestro
 * @date 10.13.2008
 * @brief Global common header with shared dependencies and methods.
 * @details We use DynamicArray and std::array for all arrays and used boost::format to simplify
 * input and output.
 * @see include/dynamic_array.h
 * @see http://www.boost.org/doc/libs/release/libs/format/index.html 
 */

#ifndef COMMON_H 
#define COMMON_H

#include <iomanip>
#include <iostream>
#include <string>
#include <cmath>
#include <cassert>
#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include <sstream>
#include <chrono>
#include <tuple>
#include <optional>
#include <functional> // has std::bind

#include "dynamic_array.h"

#define BZ_HAVE_BOOST_SERIALIZATION
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

#include <boost/format.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/ptr_container/ptr_map.hpp>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include <boost/math/tools/minima.hpp> // find minima using Brent's method

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

/* Determine if we are performing a PIGS (T=0) or PIMC (T>0) simulation */
#ifndef PIGS
#define PIGS false
#endif

/* Determine if we are studying boltzmannons */
#ifndef BOLTZMANNONS
#define BOLTZMANNONS false
#endif

/* Used for getting the repo version number into the code */
#ifndef REPO_VERSION
#define REPO_VERSION "none"
#endif

/* We Default to turning on the NN lookup table. Comment this line out for testing
 * small systems where the lookup table might actually slow things down. */
#define NN_TABLE ///< Turn on the nearest neighbor particle lookup table

#include "MersenneTwister.h"

#define NPCFSEP 50      ///< Spatial separations to be used in the pair correlation function
#define NOBDMSEP 50     ///< Spatial separations to be used in the one body density matrix
#define NRADSEP 200     ///< Spatial separations to be used in the radial density
#define NGRIDSEP 51     ///< Spatial separations to be used in each dimension of the particle position grid
#define EPS 1.0E-7      ///< A small number
#define DBL_EPS std::numeric_limits<double>::epsilon() //< Smallest double
#define BIG 1.0E30      ///< A big number
#define LBIG 69.07755279 ///< The log of a big number
#define XXX -1          ///< Used to refer to a nonsense beadIndex

using boost::format;

/** Unsigned integer type, at least 32 bits */
typedef unsigned long uint32;

/** A NDIM x NDIM matrix of type double */
typedef std::array<std::array<double, NDIM>, NDIM> dMat;

/** A NDIM-vector of type double */
typedef std::array<double,NDIM> dVec;

/** A NDIM-vector of type integer*/
typedef std::array<int,NDIM> iVec;

/** time-slice,bead-number world line index */
typedef std::array<int,2> beadLocator;

/** Integer array iterator */
typedef DynamicArray<int,1>::iterator intIter;

/** Constant integer array iterator */
typedef DynamicArray<int,1>::const_iterator cintIter;

/** beadLocator array iterator */
typedef DynamicArray<beadLocator,1>::iterator beadIter;

/** Each bead can have three possible states */
enum beadState {HEADTAIL,SPECIAL,NONE};

/** Each move can operate on only the digaonal ensemble, only the
 * off-diagonal ensemble, or both */
enum ensemble {DIAGONAL, OFFDIAGONAL, ANY};

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

/** A python-like enumerate 
 * @see https://www.reedbeta.com/blog/python-like-enumerate-in-cpp17/ 
 * */
template <typename T, 
          typename TIter = decltype(std::begin(std::declval<T>())),
          typename = decltype(std::end(std::declval<T>()))>
constexpr auto enumerate(T && iterable)
{
    struct iterator
    {
        size_t i;
        TIter iter;
        bool operator != (const iterator & other) const { return iter != other.iter; }
        void operator ++ () { ++i; ++iter; }
        auto operator * () const { return std::tie(i, *iter); }
    };
    struct iterable_wrapper
    {
        T iterable;
        auto begin() { return iterator{ 0, std::begin(iterable) }; }
        auto end() { return iterator{ 0, std::end(iterable) }; }
    };
    return iterable_wrapper{ std::forward<T>(iterable) };
}

template <typename T>
T weighted_average(const std::vector<T>& x) {

    auto [numerator, denominator] = std::accumulate(
        x.begin(), x.end(),
        std::make_pair(T{0}, T{0}),
        [i = size_t{0}](std::pair<T, T> acc, T val) mutable {
            acc.first += static_cast<T>(i) * val; // Weighted sum
            acc.second += val; // Sum of weights
            ++i;
            return acc;
        }
    );

    return denominator != T{0} ? numerator / denominator : T{0}; // Avoid division by zero
}

template <typename ResultType, typename VectorType, typename MatrixType, std::size_t Dimension>
void apply_matrix_vector_product(
    std::array<ResultType, Dimension>& resultVector,
    const std::array<VectorType, Dimension>& inputVector,
    const std::array<std::array<MatrixType, Dimension>, Dimension>& transformationMatrix)
{
    for (std::size_t row = 0; row < Dimension; ++row) {
        resultVector[row] += std::inner_product(
            transformationMatrix[row].begin(),
            transformationMatrix[row].end(),
            inputVector.begin(),
            ResultType{0}
        );
    }
}

template<typename Container>
double dot(const Container& a, const Container& b) {
    return std::inner_product(a.begin(), a.end(), b.begin(), 0.0);
}

template<typename Container, typename Predicate>
bool all(const Container& c, Predicate pred) {
    return std::all_of(c.begin(), c.end(), pred);
}

#endif
