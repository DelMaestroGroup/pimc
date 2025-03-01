#ifndef ARRAY_MATH_H
#define ARRAY_MATH_H

#include <array>
#include <cstddef>
#include <type_traits>

// ------------------------
// Vector operations
// ------------------------

// Element-wise addition
template<typename T, typename U, std::size_t N>
inline std::array<typename std::common_type<T, U>::type, N>
operator+(const std::array<T, N>& a, const std::array<U, N>& b) {
    using result_t = typename std::common_type<T, U>::type;
    std::array<result_t, N> result;
    for (std::size_t i = 0; i < N; ++i) {
        result[i] = a[i] + b[i];
    }
    return result;
}

// Element-wise subtraction
template<typename T, typename U, std::size_t N>
inline std::array<typename std::common_type<T, U>::type, N>
operator-(const std::array<T, N>& a, const std::array<U, N>& b) {
    using result_t = typename std::common_type<T, U>::type;
    std::array<result_t, N> result;
    for (std::size_t i = 0; i < N; ++i) {
        result[i] = a[i] - b[i];
    }
    return result;
}

// Element-wise multiplication
template<typename T, typename U, std::size_t N>
inline std::array<typename std::common_type<T, U>::type, N>
operator*(const std::array<T, N>& a, const std::array<U, N>& b) {
    using result_t = typename std::common_type<T, U>::type;
    std::array<result_t, N> result;
    for (std::size_t i = 0; i < N; ++i) {
        result[i] = a[i] * b[i];
    }
    return result;
}

// Element-wise division
template<typename T, typename U, std::size_t N>
inline std::array<typename std::common_type<T, U>::type, N>
operator/(const std::array<T, N>& a, const std::array<U, N>& b) {
    using result_t = typename std::common_type<T, U>::type;
    std::array<result_t, N> result;
    for (std::size_t i = 0; i < N; ++i) {
        result[i] = a[i] / b[i];
    }
    return result;
}

// Unary minus operator (negation)
template<typename T, std::size_t N>
inline std::array<T, N> operator-(const std::array<T, N>& a) {
    std::array<T, N> result;
    for (std::size_t i = 0; i < N; ++i) {
        result[i] = -a[i];
    }
    return result;
}

// Compound assignment operators for element-wise operations
// Note: The left-hand side type remains T.
template<typename T, typename U, std::size_t N>
inline std::array<T, N>& operator+=(std::array<T, N>& a, const std::array<U, N>& b) {
    for (std::size_t i = 0; i < N; ++i) {
        a[i] += b[i];
    }
    return a;
}

template<typename T, typename U, std::size_t N>
inline std::array<T, N>& operator-=(std::array<T, N>& a, const std::array<U, N>& b) {
    for (std::size_t i = 0; i < N; ++i) {
        a[i] -= b[i];
    }
    return a;
}

template<typename T, typename U, std::size_t N>
inline std::array<T, N>& operator*=(std::array<T, N>& a, const std::array<U, N>& b) {
    for (std::size_t i = 0; i < N; ++i) {
        a[i] *= b[i];
    }
    return a;
}

template<typename T, typename U, std::size_t N>
inline std::array<T, N>& operator/=(std::array<T, N>& a, const std::array<U, N>& b) {
    for (std::size_t i = 0; i < N; ++i) {
        a[i] /= b[i];
    }
    return a;
}

// Scalar operations (element-wise with a scalar)

// Multiply array by a scalar
template<typename T, typename U, std::size_t N>
inline std::array<typename std::common_type<T, U>::type, N>
operator*(const std::array<T, N>& a, const U& scalar) {
    using result_t = typename std::common_type<T, U>::type;
    std::array<result_t, N> result;
    for (std::size_t i = 0; i < N; ++i) {
        result[i] = a[i] * scalar;
    }
    return result;
}

// Multiply scalar by an array
template<typename T, typename U, std::size_t N>
inline std::array<typename std::common_type<T, U>::type, N>
operator*(const T& scalar, const std::array<U, N>& a) {
    return a * scalar; // reuse the above definition
}

// Divide array by a scalar
template<typename T, typename U, std::size_t N>
inline std::array<typename std::common_type<T, U>::type, N>
operator/(const std::array<T, N>& a, const U& scalar) {
    using result_t = typename std::common_type<T, U>::type;
    std::array<result_t, N> result;
    for (std::size_t i = 0; i < N; ++i) {
        result[i] = a[i] / scalar;
    }
    return result;
}

// Divide scalar by an array
template<typename T, typename U, std::size_t N>
inline std::array<typename std::common_type<T, U>::type, N>
operator/(const T& scalar, const std::array<U, N>& a) {
    using result_t = typename std::common_type<T, U>::type;
    std::array<result_t, N> result;
    for (std::size_t i = 0; i < N; ++i) {
        result[i] = scalar / a[i];
    }
    return result;
}

// Add array to scalar
template<typename T, typename U, std::size_t N>
inline std::array<typename std::common_type<T, U>::type, N>
operator+(const std::array<T, N>& a, const U& scalar) {
    using result_t = typename std::common_type<T, U>::type;
    std::array<result_t, N> result;
    for (std::size_t i = 0; i < N; ++i) {
        result[i] = a[i] + scalar;
    }
    return result;
}

// Add scalar to array
template<typename T, typename U, std::size_t N>
inline std::array<typename std::common_type<T, U>::type, N>
operator+(const T& scalar, const std::array<U, N>& a) {
    return a + scalar;
}

// Subtract scalar from array
template<typename T, typename U, std::size_t N>
inline std::array<typename std::common_type<T, U>::type, N>
operator-(const std::array<T, N>& a, const U& scalar) {
    using result_t = typename std::common_type<T, U>::type;
    std::array<result_t, N> result;
    for (std::size_t i = 0; i < N; ++i) {
        result[i] = a[i] - scalar;
    }
    return result;
}

// Subtract array from scalar
template<typename T, typename U, std::size_t N>
inline std::array<typename std::common_type<T, U>::type, N>
operator-(const T& scalar, const std::array<U, N>& a) {
    using result_t = typename std::common_type<T, U>::type;
    std::array<result_t, N> result;
    for (std::size_t i = 0; i < N; ++i) {
        result[i] = scalar - a[i];
    }
    return result;
}

// Compound assignment for scalar operations (left-hand side remains type T)
template<typename T, typename U, std::size_t N>
inline std::array<T, N>& operator*=(std::array<T, N>& a, const U& scalar) {
    for (std::size_t i = 0; i < N; ++i) {
        a[i] *= scalar;
    }
    return a;
}

template<typename T, typename U, std::size_t N>
inline std::array<T, N>& operator/=(std::array<T, N>& a, const U& scalar) {
    for (std::size_t i = 0; i < N; ++i) {
        a[i] /= scalar;
    }
    return a;
}

template<typename T, std::size_t N>
inline T dot(const std::array<T, N>& a, const std::array<T, N>& b) {
    T result = T();
    for (std::size_t i = 0; i < N; ++i) {
        result += a[i] * b[i];
    }
    return result;
}

template<typename T, std::size_t N>
inline T sum(const std::array<T, N>& a) {
    T result = T();
    for (std::size_t i = 0; i < N; ++i) {
        result += a[i];
    }
    return result;
}

// ------------------------
// Matrix operations
// ------------------------

// Define a type alias for rectangular matrices: a Matrix of T with ROWS rows and COLS columns.
template<typename T, std::size_t ROWS, std::size_t COLS>
using Matrix = std::array<std::array<T, COLS>, ROWS>;

// ------------------------
// Element-wise Operations for matrices
// ------------------------

// Element-wise addition of two matrices with the same dimensions.
template<typename T, typename U, std::size_t ROWS, std::size_t COLS>
inline std::array<std::array<typename std::common_type<T, U>::type, COLS>, ROWS>
operator+(const Matrix<T, ROWS, COLS>& A, const Matrix<U, ROWS, COLS>& B)
{
    using result_t = typename std::common_type<T, U>::type;
    std::array<std::array<result_t, COLS>, ROWS> result;
    for (std::size_t i = 0; i < ROWS; ++i)
        for (std::size_t j = 0; j < COLS; ++j)
            result[i][j] = A[i][j] + B[i][j];
    return result;
}

// Element-wise subtraction of two matrices with the same dimensions.
template<typename T, typename U, std::size_t ROWS, std::size_t COLS>
inline std::array<std::array<typename std::common_type<T, U>::type, COLS>, ROWS>
operator-(const Matrix<T, ROWS, COLS>& A, const Matrix<U, ROWS, COLS>& B)
{
    using result_t = typename std::common_type<T, U>::type;
    std::array<std::array<result_t, COLS>, ROWS> result;
    for (std::size_t i = 0; i < ROWS; ++i)
        for (std::size_t j = 0; j < COLS; ++j)
            result[i][j] = A[i][j] - B[i][j];
    return result;
}

// Unary minus operator (matrix negation)
template<typename T, std::size_t ROWS, std::size_t COLS>
inline Matrix<T, ROWS, COLS> operator-(const Matrix<T, ROWS, COLS>& A)
{
    Matrix<T, ROWS, COLS> result;
    for (std::size_t i = 0; i < ROWS; ++i)
        for (std::size_t j = 0; j < COLS; ++j)
            result[i][j] = -A[i][j];
    return result;
}

// Compound assignment for addition.
template<typename T, typename U, std::size_t ROWS, std::size_t COLS>
inline Matrix<T, ROWS, COLS>& operator+=(Matrix<T, ROWS, COLS>& A,
                                           const Matrix<U, ROWS, COLS>& B)
{
    for (std::size_t i = 0; i < ROWS; ++i)
        for (std::size_t j = 0; j < COLS; ++j)
            A[i][j] += B[i][j];
    return A;
}

// Compound assignment for subtraction.
template<typename T, typename U, std::size_t ROWS, std::size_t COLS>
inline Matrix<T, ROWS, COLS>& operator-=(Matrix<T, ROWS, COLS>& A,
                                           const Matrix<U, ROWS, COLS>& B)
{
    for (std::size_t i = 0; i < ROWS; ++i)
        for (std::size_t j = 0; j < COLS; ++j)
            A[i][j] -= B[i][j];
    return A;
}

// ------------------------
// Scalar Operations for matrices
// ------------------------

// Multiply a matrix by a scalar.
template<typename T, typename U, std::size_t ROWS, std::size_t COLS>
inline std::array<std::array<typename std::common_type<T, U>::type, COLS>, ROWS>
operator*(const Matrix<T, ROWS, COLS>& A, const U& scalar)
{
    using result_t = typename std::common_type<T, U>::type;
    std::array<std::array<result_t, COLS>, ROWS> result;
    for (std::size_t i = 0; i < ROWS; ++i)
        for (std::size_t j = 0; j < COLS; ++j)
            result[i][j] = A[i][j] * scalar;
    return result;
}

// Multiply a scalar by a matrix.
template<typename T, typename U, std::size_t ROWS, std::size_t COLS>
inline std::array<std::array<typename std::common_type<T, U>::type, COLS>, ROWS>
operator*(const U& scalar, const Matrix<T, ROWS, COLS>& A)
{
    return A * scalar;
}

// Divide a matrix by a scalar.
template<typename T, typename U, std::size_t ROWS, std::size_t COLS>
inline std::array<std::array<typename std::common_type<T, U>::type, COLS>, ROWS>
operator/(const Matrix<T, ROWS, COLS>& A, const U& scalar)
{
    using result_t = typename std::common_type<T, U>::type;
    std::array<std::array<result_t, COLS>, ROWS> result;
    for (std::size_t i = 0; i < ROWS; ++i)
        for (std::size_t j = 0; j < COLS; ++j)
            result[i][j] = A[i][j] / scalar;
    return result;
}

// Divide a scalar by a matrix.
template<typename T, typename U, std::size_t ROWS, std::size_t COLS>
inline std::array<std::array<typename std::common_type<T, U>::type, COLS>, ROWS>
operator/(const U& scalar, const Matrix<T, ROWS, COLS>& A)
{
    using result_t = typename std::common_type<T, U>::type;
    std::array<std::array<result_t, COLS>, ROWS> result;
    for (std::size_t i = 0; i < ROWS; ++i)
        for (std::size_t j = 0; j < COLS; ++j)
            result[i][j] = scalar / A[i][j];
    return result;
}

// Add a matrix to a scalar.
template<typename T, typename U, std::size_t ROWS, std::size_t COLS>
inline std::array<std::array<typename std::common_type<T, U>::type, COLS>, ROWS>
operator+(const Matrix<T, ROWS, COLS>& A, const U& scalar)
{
    using result_t = typename std::common_type<T, U>::type;
    std::array<std::array<result_t, COLS>, ROWS> result;
    for (std::size_t i = 0; i < ROWS; ++i)
        for (std::size_t j = 0; j < COLS; ++j)
            result[i][j] = A[i][j] + scalar;
    return result;
}

// Add a scalar to a matrix.
template<typename T, typename U, std::size_t ROWS, std::size_t COLS>
inline std::array<std::array<typename std::common_type<T, U>::type, COLS>, ROWS>
operator+(const U& scalar, const Matrix<T, ROWS, COLS>& A)
{
    return A + scalar;
}

// Subtract a scalar from a matrix.
template<typename T, typename U, std::size_t ROWS, std::size_t COLS>
inline std::array<std::array<typename std::common_type<T, U>::type, COLS>, ROWS>
operator-(const Matrix<T, ROWS, COLS>& A, const U& scalar)
{
    using result_t = typename std::common_type<T, U>::type;
    std::array<std::array<result_t, COLS>, ROWS> result;
    for (std::size_t i = 0; i < ROWS; ++i)
        for (std::size_t j = 0; j < COLS; ++j)
            result[i][j] = A[i][j] - scalar;
    return result;
}

// Subtract a matrix from a scalar.
template<typename T, typename U, std::size_t ROWS, std::size_t COLS>
inline std::array<std::array<typename std::common_type<T, U>::type, COLS>, ROWS>
operator-(const U& scalar, const Matrix<T, ROWS, COLS>& A)
{
    using result_t = typename std::common_type<T, U>::type;
    std::array<std::array<result_t, COLS>, ROWS> result;
    for (std::size_t i = 0; i < ROWS; ++i)
        for (std::size_t j = 0; j < COLS; ++j)
            result[i][j] = scalar - A[i][j];
    return result;
}

// Compound assignment for scalar multiplication.
template<typename T, typename U, std::size_t ROWS, std::size_t COLS>
inline Matrix<T, ROWS, COLS>& operator*=(Matrix<T, ROWS, COLS>& A, const U& scalar)
{
    for (std::size_t i = 0; i < ROWS; ++i)
        for (std::size_t j = 0; j < COLS; ++j)
            A[i][j] *= scalar;
    return A;
}

// Compound assignment for scalar division.
template<typename T, typename U, std::size_t ROWS, std::size_t COLS>
inline Matrix<T, ROWS, COLS>& operator/=(Matrix<T, ROWS, COLS>& A, const U& scalar)
{
    for (std::size_t i = 0; i < ROWS; ++i)
        for (std::size_t j = 0; j < COLS; ++j)
            A[i][j] /= scalar;
    return A;
}

// ------------------------
// Matrix Multiplication
// ------------------------

// Conventional matrix multiplication.
// If A is of dimensions M x N and B is of dimensions N x P,
// then the product A * B is a matrix of dimensions M x P.
// To allow mixed types in A, B, and the sum we introduce three type parameters.
template<typename T, typename U, typename V, std::size_t M, std::size_t N, std::size_t P>
inline std::array<std::array<typename std::common_type<T, U, V>::type, P>, M>
operator*(const Matrix<T, M, N>& A, const std::array<std::array<U, P>, N>& B)
{
    using result_t = typename std::common_type<T, U, V>::type;
    std::array<std::array<result_t, P>, M> result{};
    for (std::size_t i = 0; i < M; ++i)
    {
        for (std::size_t j = 0; j < P; ++j)
        {
            result[i][j] = result_t(); // value-initialize to 0
            for (std::size_t k = 0; k < N; ++k)
            {
                result[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return result;
}

// ------------------------
// Hadamard Product
// ------------------------

// Element-wise (Hadamard) multiplication of two matrices of the same size.
template<typename T, typename U, std::size_t ROWS, std::size_t COLS>
inline std::array<std::array<typename std::common_type<T, U>::type, COLS>, ROWS>
hadamard(const Matrix<T, ROWS, COLS>& A, const Matrix<U, ROWS, COLS>& B)
{
    using result_t = typename std::common_type<T, U>::type;
    std::array<std::array<result_t, COLS>, ROWS> result;
    for (std::size_t i = 0; i < ROWS; ++i)
        for (std::size_t j = 0; j < COLS; ++j)
            result[i][j] = A[i][j] * B[i][j];
    return result;
}

#endif // ARRAY_MATH_H

