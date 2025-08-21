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
#include <utility>

#include "array_math.h"
#include "dynamic_array_math.h"

#include <boost/format.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/ptr_container/ptr_map.hpp>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include <boost/math/tools/minima.hpp> // find minima using Brent's method

/* Debugging librarys and definitions. There is only an effect if
 * the code is compiled with PIMC_DEBUG on.*/
#ifdef PIMC_DEBUG

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

// Weighted average function for 1D containers.
template <typename Container>
auto weighted_average(const Container& container) -> double {
    // Get the extent (assumes a 1D container).
    const std::size_t n = container.extents()[0];
    
    // Use double for accumulation.
    double weightedTotal = 0.0;
    double total = 0.0;
    
    for (std::size_t i = 0; i < n; ++i) {
        auto value = container(i);
        weightedTotal += static_cast<double>(i) * static_cast<double>(value);
        total += static_cast<double>(value);
    }
    
    // FIXME not sure if I should throw here instead of returning 0.0
    return total == 0.0 ? 0.0 : weightedTotal / total;
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

template <typename T, std::size_t N>
bool all(const std::array<T, N>& arr, T value) {
    for (const auto& elem : arr) {
        if (elem != value) {
            return false;
        }
    }
    return true;
}

// Helper function to generate an array filled with a given value
template <typename T, std::size_t N, std::size_t... Is>
constexpr std::array<T, N> make_array_impl(const T& value, std::index_sequence<Is...>) {
    return {{ (static_cast<void>(Is), value)... }};
}

// Helper to make and fill array of given type and size
template <typename T, std::size_t N>
constexpr std::array<T, N> make_array(const T& value) {
    return make_array_impl<T, N>(value, std::make_index_sequence<N>{});
}

// Overload that deduces type and size from a given array type.
template <typename ArrayType>
constexpr ArrayType make_array(const typename ArrayType::value_type& value) {
    return make_array_impl<typename ArrayType::value_type,
                           std::tuple_size<ArrayType>::value>(
                           value,
                           std::make_index_sequence<std::tuple_size<ArrayType>::value>{});
}

// FIXME These stream overloads should be backwards compatible. Should probably move to communicator.h
// Overload for printing any std::array<T, N>
template<typename T, std::size_t N>
std::ostream& operator<<(std::ostream& os, const std::array<T, N>& arr) {
    os << "(";
    for (std::size_t i = 0; i < N; ++i) {
        os << arr[i];
        if (i < N - 1) {
            os << ",";
        }
    }
    os << ")";
    return os;
}

// Overload for printing 2D DynamicArray
template<typename T>
std::ostream& operator<<(std::ostream& os, const DynamicArray<T, 2>& arr) {
    auto extents = arr.extents();
    os << "(" << 0 << "," << extents[0] - 1 << ") x ("
       << 0 << "," << extents[1] - 1 << ")\n";
    os << "[ ";
    for (std::size_t i = 0; i < extents[0]; ++i) {
        for (std::size_t j = 0; j < extents[1]; ++j) {
            os << arr(i, j) << " ";
        }
        if (i < extents[0] - 1)
            os << "\n  ";
    }
    os << "]" << std::endl;
    return os;
}

// Overload for streaming into std::array<T, N>
template<typename T, std::size_t N>
std::istream& operator>>(std::istream& is, std::array<T, N>& a) {
    char ch;
    // Expect an opening '('.
    is >> ch;
    if (ch != '(') {
        is.setstate(std::ios::failbit);
        return is;
    }
    for (std::size_t i = 0; i < N; ++i) {
        is >> a[i];
        if (i < N - 1) {
            // Expect a comma.
            is >> ch;
            if (ch != ',') {
                is.setstate(std::ios::failbit);
                return is;
            }
        }
    }
    // Expect a closing ')'.
    is >> ch;
    if (ch != ')') {
        is.setstate(std::ios::failbit);
    }
    return is;
}

// Overload for streaming into 2D DynamicArray
template<typename T>
std::istream& operator>>(std::istream& is, DynamicArray<T, 2>& arr)
{
    // Should start with: (0,extents[0]-1) x (0,extents[1]-1)
    int rowStart, rowEnd, colStart, colEnd;
    char ch;

    // Read the first '('
    is >> ch;                 // should be '('
    if (!is || ch != '(') {
        is.setstate(std::ios::failbit);
        return is;
    }

    // Read rowStart, then comma
    is >> rowStart >> ch;     // rowStart, then ','
    if (!is || ch != ',') {
        is.setstate(std::ios::failbit);
        return is;
    }

    // Read rowEnd, then ')'
    is >> rowEnd >> ch;       // rowEnd, then ')'
    if (!is || ch != ')') {
        is.setstate(std::ios::failbit);
        return is;
    }

    // Read 'x'
    is >> ch;                 // should be 'x'
    if (!is || ch != 'x') {
        is.setstate(std::ios::failbit);
        return is;
    }

    // Read '(' for columns
    is >> ch;                 // should be '('
    if (!is || ch != '(') {
        is.setstate(std::ios::failbit);
        return is;
    }

    // Read colStart, then comma
    is >> colStart >> ch;     // colStart, then ','
    if (!is || ch != ',') {
        is.setstate(std::ios::failbit);
        return is;
    }

    // Read colEnd, then ')'
    is >> colEnd >> ch;       // colEnd, then ')'
    if (!is || ch != ')') {
        is.setstate(std::ios::failbit);
        return is;
    }

    // Resize arr to (rowEnd - rowStart + 1) x (colEnd - colStart + 1)
    arr.resize(std::size_t(rowEnd - rowStart + 1),
               std::size_t(colEnd - colStart + 1));

    // Next we expect a line beginning with '['
    // so skip forward until we see '[':
    while (is >> ch) {
        if (ch == '[')
            break;
    }
    if (!is) {
        is.setstate(std::ios::failbit);
        return is;
    }

    // Now read arr's contents, row by row
    auto extents = arr.extents();
    for (std::size_t i = 0; i < extents[0]; ++i) {
        for (std::size_t j = 0; j < extents[1]; ++j) {
            if (!(is >> arr(i, j))) {
                is.setstate(std::ios::failbit);
                return is;
            }
        }
    }

    // Consume characters until the closing ']' is found
    while (is >> ch) {
        if (ch == ']')
            break;
    }
    if (!is) {
        is.setstate(std::ios::failbit);
        return is;
    }

    return is;
}

// Returns true if no element in the container equals the given value.
template <typename Container, typename T>
constexpr bool noneEquals(const Container& container, const T& value) {
    return std::none_of(container.begin(), container.end(), [&value](const auto& elem) {
        return elem == value;
    });
}

// Returns true if any element in the container equals the given value.
template <typename Container, typename T>
constexpr bool anyEquals(const Container& container, const T& value) {
    return std::any_of(container.begin(), container.end(), [&value](const auto& elem) {
        return elem == value;
    });
}

// Returns true if all elements in the container equal the given value.
template <typename Container, typename T>
constexpr bool allEquals(const Container& container, const T& value) {
    return std::all_of(container.begin(), container.end(), [&value](const auto& elem) {
        return elem == value;
    });
}

// Helper function that performs elementwise comparison using a fold expression.
// The parameter pack I... represents a compile-time sequence of indices.
template <typename T, std::size_t N, std::size_t... I>
constexpr bool all_impl(const std::array<T, N>& a, const std::array<T, N>& b, std::index_sequence<I...>) {
    // The fold expression below expands to a series of element comparisons:
    // (a[0] == b[0]) && (a[1] == b[1]) && ... && (a[N-1] == b[N-1])
    return ((a[I] == b[I]) && ...);
}

// Primary function that compares two std::array objects elementwise.
// It generates a compile-time sequence of indices [0, 1, ..., N-1] using std::make_index_sequence.
template <typename T, std::size_t N>
constexpr bool all(const std::array<T, N>& a, const std::array<T, N>& b) {
    // std::make_index_sequence<N>{} creates an index sequence for all valid array indices.
    // The helper function all_impl is then called with this sequence to perform the comparisons.
    return all_impl(a, b, std::make_index_sequence<N>{});
}

#endif

