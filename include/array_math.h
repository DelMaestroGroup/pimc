#ifndef ARRAY_MATH_H
#define ARRAY_MATH_H

#include <array>
#include <cstddef>

// Element-wise addition
template<typename T, std::size_t N>
inline std::array<T, N> operator+(const std::array<T, N>& a, const std::array<T, N>& b) {
    std::array<T, N> result;
    for (std::size_t i = 0; i < N; ++i) {
        result[i] = a[i] + b[i];
    }
    return result;
}

// Element-wise subtraction
template<typename T, std::size_t N>
inline std::array<T, N> operator-(const std::array<T, N>& a, const std::array<T, N>& b) {
    std::array<T, N> result;
    for (std::size_t i = 0; i < N; ++i) {
        result[i] = a[i] - b[i];
    }
    return result;
}

// Element-wise multiplication
template<typename T, std::size_t N>
inline std::array<T, N> operator*(const std::array<T, N>& a, const std::array<T, N>& b) {
    std::array<T, N> result;
    for (std::size_t i = 0; i < N; ++i) {
        result[i] = a[i] * b[i];
    }
    return result;
}

// Element-wise division
template<typename T, std::size_t N>
inline std::array<T, N> operator/(const std::array<T, N>& a, const std::array<T, N>& b) {
    std::array<T, N> result;
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

template<typename T, std::size_t N>
inline std::array<T, N>& operator+=(std::array<T, N>& a, const std::array<T, N>& b) {
    for (std::size_t i = 0; i < N; ++i) {
        a[i] += b[i];
    }
    return a;
}

template<typename T, std::size_t N>
inline std::array<T, N>& operator-=(std::array<T, N>& a, const std::array<T, N>& b) {
    for (std::size_t i = 0; i < N; ++i) {
        a[i] -= b[i];
    }
    return a;
}

template<typename T, std::size_t N>
inline std::array<T, N>& operator*=(std::array<T, N>& a, const std::array<T, N>& b) {
    for (std::size_t i = 0; i < N; ++i) {
        a[i] *= b[i];
    }
    return a;
}

template<typename T, std::size_t N>
inline std::array<T, N>& operator/=(std::array<T, N>& a, const std::array<T, N>& b) {
    for (std::size_t i = 0; i < N; ++i) {
        a[i] /= b[i];
    }
    return a;
}

// Scalar operations (element-wise with a scalar)

// Multiply array by a scalar
template<typename T, std::size_t N>
inline std::array<T, N> operator*(const std::array<T, N>& a, const T& scalar) {
    std::array<T, N> result;
    for (std::size_t i = 0; i < N; ++i) {
        result[i] = a[i] * scalar;
    }
    return result;
}

// Multiply scalar by an array
template<typename T, std::size_t N>
inline std::array<T, N> operator*(const T& scalar, const std::array<T, N>& a) {
    return a * scalar; // reuse the above definition
}

// Divide array by a scalar
template<typename T, std::size_t N>
inline std::array<T, N> operator/(const std::array<T, N>& a, const T& scalar) {
    std::array<T, N> result;
    for (std::size_t i = 0; i < N; ++i) {
        result[i] = a[i] / scalar;
    }
    return result;
}

// Compound assignment for scalar operations

template<typename T, std::size_t N>
inline std::array<T, N>& operator*=(std::array<T, N>& a, const T& scalar) {
    for (std::size_t i = 0; i < N; ++i) {
        a[i] *= scalar;
    }
    return a;
}

template<typename T, std::size_t N>
inline std::array<T, N>& operator/=(std::array<T, N>& a, const T& scalar) {
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

#endif // ARRAY_MATH_H
