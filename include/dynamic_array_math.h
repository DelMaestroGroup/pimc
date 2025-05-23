#ifndef DYNAMIC_ARRAY_MATH_H
#define DYNAMIC_ARRAY_MATH_H

#include "dynamic_array.h"
#include <algorithm>
#include <numeric>
#include <functional>
#include <stdexcept>
#include <type_traits>
#include <array>
#include <iterator>

//------------------------------------------------------------------------------
// Helper: Create a DynamicArray from a given extents array.
//------------------------------------------------------------------------------
namespace dynamic_array_math {

template<typename T, std::size_t Rank, std::size_t... I>
DynamicArray<T, Rank> create_dynamic_array_from_extents_impl(
    const std::array<std::size_t, Rank>& extents,
    std::index_sequence<I...>)
{
    // Unpack the extents into the DynamicArray constructor.
    return DynamicArray<T, Rank>( extents[I]... );
}

template<typename T, std::size_t Rank>
DynamicArray<T, Rank> create_dynamic_array_from_extents(
    const std::array<std::size_t, Rank>& extents)
{
    return create_dynamic_array_from_extents_impl<T, Rank>(
        extents, std::make_index_sequence<Rank>{});
}

} // namespace dynamic_array_math

//------------------------------------------------------------------------------
// Element-wise Operations
//------------------------------------------------------------------------------

// Element-wise addition
template<typename T, typename U, std::size_t Rank>
auto operator+(const DynamicArray<T, Rank>& a, const DynamicArray<U, Rank>& b)
    -> DynamicArray<typename std::common_type<T, U>::type, Rank>
{
    if (a.extents() != b.extents())
        throw std::invalid_argument("DynamicArray dimensions must match for addition.");
    using result_t = typename std::common_type<T, U>::type;
    DynamicArray<result_t, Rank> result =
        dynamic_array_math::create_dynamic_array_from_extents<result_t, Rank>(a.extents());
    std::transform(a.begin(), a.end(), b.begin(), result.begin(), std::plus<result_t>());
    return result;
}

// Element-wise subtraction
template<typename T, typename U, std::size_t Rank>
auto operator-(const DynamicArray<T, Rank>& a, const DynamicArray<U, Rank>& b)
    -> DynamicArray<typename std::common_type<T, U>::type, Rank>
{
    if (a.extents() != b.extents())
        throw std::invalid_argument("DynamicArray dimensions must match for subtraction.");
    using result_t = typename std::common_type<T, U>::type;
    DynamicArray<result_t, Rank> result =
        dynamic_array_math::create_dynamic_array_from_extents<result_t, Rank>(a.extents());
    std::transform(a.begin(), a.end(), b.begin(), result.begin(), std::minus<result_t>());
    return result;
}

// Element-wise multiplication
template<typename T, typename U, std::size_t Rank>
auto operator*(const DynamicArray<T, Rank>& a, const DynamicArray<U, Rank>& b)
    -> DynamicArray<typename std::common_type<T, U>::type, Rank>
{
    if (a.extents() != b.extents())
        throw std::invalid_argument("DynamicArray dimensions must match for element-wise multiplication.");
    using result_t = typename std::common_type<T, U>::type;
    DynamicArray<result_t, Rank> result =
        dynamic_array_math::create_dynamic_array_from_extents<result_t, Rank>(a.extents());
    std::transform(a.begin(), a.end(), b.begin(), result.begin(), std::multiplies<result_t>());
    return result;
}

// Element-wise division
template<typename T, typename U, std::size_t Rank>
auto operator/(const DynamicArray<T, Rank>& a, const DynamicArray<U, Rank>& b)
    -> DynamicArray<typename std::common_type<T, U>::type, Rank>
{
    if (a.extents() != b.extents())
        throw std::invalid_argument("DynamicArray dimensions must match for element-wise division.");
    using result_t = typename std::common_type<T, U>::type;
    DynamicArray<result_t, Rank> result =
        dynamic_array_math::create_dynamic_array_from_extents<result_t, Rank>(a.extents());
    std::transform(a.begin(), a.end(), b.begin(), result.begin(), std::divides<result_t>());
    return result;
}

// Unary minus (element-wise negation)
template<typename T, std::size_t Rank>
auto operator-(const DynamicArray<T, Rank>& a)
    -> DynamicArray<T, Rank>
{
    DynamicArray<T, Rank> result =
        dynamic_array_math::create_dynamic_array_from_extents<T, Rank>(a.extents());
    std::transform(a.begin(), a.end(), result.begin(), std::negate<T>());
    return result;
}

//------------------------------------------------------------------------------
// Compound Assignment Operators
//------------------------------------------------------------------------------

// Compound assignment with another DynamicArray
template<typename T, typename U, std::size_t Rank>
DynamicArray<T, Rank>& operator+=(DynamicArray<T, Rank>& a, const DynamicArray<U, Rank>& b) {
    if (a.extents() != b.extents())
        throw std::invalid_argument("DynamicArray dimensions must match for addition assignment.");
    auto it_a = a.begin();
    auto it_b = b.begin();
    for (; it_a != a.end(); ++it_a, ++it_b)
        *it_a += *it_b;
    return a;
}

template<typename T, typename U, std::size_t Rank>
DynamicArray<T, Rank>& operator-=(DynamicArray<T, Rank>& a, const DynamicArray<U, Rank>& b) {
    if (a.extents() != b.extents())
        throw std::invalid_argument("DynamicArray dimensions must match for subtraction assignment.");
    auto it_a = a.begin();
    auto it_b = b.begin();
    for (; it_a != a.end(); ++it_a, ++it_b)
        *it_a -= *it_b;
    return a;
}

template<typename T, typename U, std::size_t Rank>
DynamicArray<T, Rank>& operator*=(DynamicArray<T, Rank>& a, const DynamicArray<U, Rank>& b) {
    if (a.extents() != b.extents())
        throw std::invalid_argument("DynamicArray dimensions must match for multiplication assignment.");
    auto it_a = a.begin();
    auto it_b = b.begin();
    for (; it_a != a.end(); ++it_a, ++it_b)
        *it_a *= *it_b;
    return a;
}

template<typename T, typename U, std::size_t Rank>
DynamicArray<T, Rank>& operator/=(DynamicArray<T, Rank>& a, const DynamicArray<U, Rank>& b) {
    if (a.extents() != b.extents())
        throw std::invalid_argument("DynamicArray dimensions must match for division assignment.");
    auto it_a = a.begin();
    auto it_b = b.begin();
    for (; it_a != a.end(); ++it_a, ++it_b)
        *it_a /= *it_b;
    return a;
}

//------------------------------------------------------------------------------
// Scalar Operations
//------------------------------------------------------------------------------

// Multiply DynamicArray by a scalar
template<typename T, typename U, std::size_t Rank>
auto operator*(const DynamicArray<T, Rank>& a, const U& scalar)
    -> DynamicArray<typename std::common_type<T, U>::type, Rank>
{
    using result_t = typename std::common_type<T, U>::type;
    DynamicArray<result_t, Rank> result =
        dynamic_array_math::create_dynamic_array_from_extents<result_t, Rank>(a.extents());
    std::transform(a.begin(), a.end(), result.begin(),
                   [scalar](const T& x) { return x * scalar; });
    return result;
}

// Multiply scalar by DynamicArray
template<typename T, typename U, std::size_t Rank>
auto operator*(const U& scalar, const DynamicArray<T, Rank>& a)
    -> DynamicArray<typename std::common_type<T, U>::type, Rank>
{
    return a * scalar;
}

// Divide DynamicArray by a scalar
template<typename T, typename U, std::size_t Rank>
auto operator/(const DynamicArray<T, Rank>& a, const U& scalar)
    -> DynamicArray<typename std::common_type<T, U>::type, Rank>
{
    using result_t = typename std::common_type<T, U>::type;
    DynamicArray<result_t, Rank> result =
        dynamic_array_math::create_dynamic_array_from_extents<result_t, Rank>(a.extents());
    std::transform(a.begin(), a.end(), result.begin(),
                   [scalar](const T& x) { return x / scalar; });
    return result;
}

// Divide scalar by DynamicArray (element-wise: scalar divided by each element)
template<typename T, typename U, std::size_t Rank>
auto operator/(const U& scalar, const DynamicArray<T, Rank>& a)
    -> DynamicArray<typename std::common_type<T, U>::type, Rank>
{
    using result_t = typename std::common_type<T, U>::type;
    DynamicArray<result_t, Rank> result =
        dynamic_array_math::create_dynamic_array_from_extents<result_t, Rank>(a.extents());
    std::transform(a.begin(), a.end(), result.begin(),
                   [scalar](const T& x) { return scalar / x; });
    return result;
}

// Add scalar to DynamicArray
template<typename T, typename U, std::size_t Rank>
auto operator+(const DynamicArray<T, Rank>& a, const U& scalar)
    -> DynamicArray<typename std::common_type<T, U>::type, Rank>
{
    using result_t = typename std::common_type<T, U>::type;
    DynamicArray<result_t, Rank> result =
        dynamic_array_math::create_dynamic_array_from_extents<result_t, Rank>(a.extents());
    std::transform(a.begin(), a.end(), result.begin(),
                   [scalar](const T& x) { return x + scalar; });
    return result;
}

// Add DynamicArray to scalar
template<typename T, typename U, std::size_t Rank>
auto operator+(const U& scalar, const DynamicArray<T, Rank>& a)
    -> DynamicArray<typename std::common_type<T, U>::type, Rank>
{
    return a + scalar;
}

// Subtract scalar from DynamicArray
template<typename T, typename U, std::size_t Rank>
auto operator-(const DynamicArray<T, Rank>& a, const U& scalar)
    -> DynamicArray<typename std::common_type<T, U>::type, Rank>
{
    using result_t = typename std::common_type<T, U>::type;
    DynamicArray<result_t, Rank> result =
        dynamic_array_math::create_dynamic_array_from_extents<result_t, Rank>(a.extents());
    std::transform(a.begin(), a.end(), result.begin(),
                   [scalar](const T& x) { return x - scalar; });
    return result;
}

// Subtract DynamicArray from scalar
template<typename T, typename U, std::size_t Rank>
auto operator-(const U& scalar, const DynamicArray<T, Rank>& a)
    -> DynamicArray<typename std::common_type<T, U>::type, Rank>
{
    using result_t = typename std::common_type<T, U>::type;
    DynamicArray<result_t, Rank> result =
        dynamic_array_math::create_dynamic_array_from_extents<result_t, Rank>(a.extents());
    std::transform(a.begin(), a.end(), result.begin(),
                   [scalar](const T& x) { return scalar - x; });
    return result;
}

// Compound assignment for scalar addition
template<typename T, typename U, std::size_t Rank>
DynamicArray<T, Rank>& operator+=(DynamicArray<T, Rank>& a, const U& scalar) {
    for (auto& x : a)
        x += scalar;
    return a;
}

// Compound assignment for scalar subtraction
template<typename T, typename U, std::size_t Rank>
DynamicArray<T, Rank>& operator-=(DynamicArray<T, Rank>& a, const U& scalar) {
    for (auto& x : a)
        x -= scalar;
    return a;
}

// Compound assignment for scalar multiplication
template<typename T, typename U, std::size_t Rank>
DynamicArray<T, Rank>& operator*=(DynamicArray<T, Rank>& a, const U& scalar) {
    for (auto& x : a)
        x *= scalar;
    return a;
}

// Compound assignment for scalar division
template<typename T, typename U, std::size_t Rank>
DynamicArray<T, Rank>& operator/=(DynamicArray<T, Rank>& a, const U& scalar) {
    for (auto& x : a)
        x /= scalar;
    return a;
}

//------------------------------------------------------------------------------
// Reduction Functions
//------------------------------------------------------------------------------

// Sum of all elements in a DynamicArray.
template<typename T, std::size_t Rank>
T sum(const DynamicArray<T, Rank>& a) {
    return std::accumulate(a.begin(), a.end(), T());
}

// Product of all elements in a DynamicArray.
template<typename T, std::size_t Rank>
T product(const DynamicArray<T, Rank>& a) {
    return std::accumulate(a.begin(), a.end(), T(1), std::multiplies<T>());
}

// Minimum element in a DynamicArray.
template<typename T, std::size_t Rank>
T min(const DynamicArray<T, Rank>& a) {
    if (a.size() == 0)
        throw std::invalid_argument("DynamicArray is empty; cannot compute minimum.");
    return *std::min_element(a.begin(), a.end());
}

// Maximum element in a DynamicArray.
template<typename T, std::size_t Rank>
T max(const DynamicArray<T, Rank>& a) {
    if (a.size() == 0)
        throw std::invalid_argument("DynamicArray is empty; cannot compute maximum.");
    return *std::max_element(a.begin(), a.end());
}

// Dot product for 1D DynamicArrays (vectors).
template<typename T>
T dot(const DynamicArray<T, 1>& a, const DynamicArray<T, 1>& b) {
    if (a.size() != b.size())
        throw std::invalid_argument("Vector sizes must match for dot product.");
    T result = T();
    auto it_a = a.begin();
    auto it_b = b.begin();
    for (; it_a != a.end(); ++it_a, ++it_b)
        result += (*it_a * *it_b);
    return result;
}

//------------------------------------------------------------------------------
// Matrix Multiplication (for 2D DynamicArrays)
//------------------------------------------------------------------------------

// Conventional matrix multiplication: if A is MxN and B is NxP,
// then the product is an MxP matrix.
// (Note: This is provided as a free function named 'matmul' to avoid
//  conflict with element-wise multiplication.)
template<typename T, typename U, typename V>
auto matmul(const DynamicArray<T, 2>& A, const DynamicArray<U, 2>& B)
    -> DynamicArray<typename std::common_type<T, U, V>::type, 2>
{
    auto A_extents = A.extents();
    auto B_extents = B.extents();
    if (A_extents[1] != B_extents[0])
        throw std::invalid_argument("Matrix multiplication dimension mismatch.");
    using result_t = typename std::common_type<T, U, V>::type;
    DynamicArray<result_t, 2> result(A_extents[0], B_extents[1]);
    for (std::size_t i = 0; i < A_extents[0]; ++i) {
        for (std::size_t j = 0; j < B_extents[1]; ++j) {
            result(i, j) = result_t();
            for (std::size_t k = 0; k < A_extents[1]; ++k)
                result(i, j) += A(i, k) * B(k, j);
        }
    }
    return result;
}

#endif // DYNAMIC_ARRAY_MATH_H

