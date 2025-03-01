#include <vector>
#include <array>
#include <cstddef>
//#include <mdspan>
#include "mdspan.h"
#include <numeric>
#include <stdexcept>
#include <algorithm>
#include <utility>
#include <type_traits>

// Helper: fill a 1D mdspan with a linear function f(index).
template <typename Mdspan, typename Func>
void fill_with_function(Mdspan view, Func f) {
    for (std::size_t i = 0; i < view.extent(0); ++i) {
        view(i) = f(i);
    }
}

// Helper to fill a slice
template <typename Mdspan>
void fill_mdspan(Mdspan view, const typename Mdspan::value_type& value) {
    std::size_t total = 1;
    for (std::size_t d = 0; d < Mdspan::rank(); ++d) {
        total *= view.extent(d);
    }
    std::fill(view.data(), view.data() + total, value);
}

// Helper to compute row‐major strides given extents.
template <std::size_t Rank>
std::array<std::size_t, Rank> compute_strides(const std::array<std::size_t, Rank>& extents) {
    std::array<std::size_t, Rank> strides{};
    strides[Rank - 1] = 1;
    for (std::size_t i = Rank - 1; i-- > 0; ) {
        strides[i] = strides[i + 1] * extents[i + 1];
    }
    return strides;
}

template <typename T, std::size_t Rank>
class DynamicArray {
public:
    using dextents_type  = std::dextents<std::size_t, Rank>;
    using mdspan_type    = std::mdspan<T, dextents_type>;
    using iterator       = typename std::vector<T>::iterator;
    using const_iterator = typename std::vector<T>::const_iterator;

    // Constructor: accepts exactly Rank extents (e.g., DynamicArray<double,2> arr(3,4);)
    template <typename... Extents,
              typename = std::enable_if_t<(sizeof...(Extents) == Rank) &&
                                          (std::conjunction_v<std::is_convertible<Extents, std::size_t>...>)>>
    DynamicArray(Extents... extents)
      : extents_{ static_cast<std::size_t>(extents)... },
        data_( (static_cast<std::size_t>(extents) * ...) ),
        view_(data_.data(), make_dextents(extents_))
    { }

    // Default constructor (empty array)
    DynamicArray() = default;

    // Copy constructor: reinitializes the mdspan view over the new vector storage.
    DynamicArray(const DynamicArray& other)
      : extents_(other.extents_),
        data_(other.data_),
        view_(data_.data(), make_dextents(extents_))
    { }

    // Move constructor
    DynamicArray(DynamicArray&& other) noexcept
      : extents_(std::move(other.extents_)),
        data_(std::move(other.data_)),
        view_(data_.data(), make_dextents(extents_))
    { }

    // Copy assignment operator
    DynamicArray& operator=(const DynamicArray& other) {
        if (this != &other) {
            extents_ = other.extents_;
            data_    = other.data_;
            view_    = mdspan_type(data_.data(), make_dextents(extents_));
        }
        return *this;
    }

    // Move assignment operator
    DynamicArray& operator=(DynamicArray&& other) noexcept {
        if (this != &other) {
            extents_ = std::move(other.extents_);
            data_    = std::move(other.data_);
            view_    = mdspan_type(data_.data(), make_dextents(extents_));
        }
        return *this;
    }

    // Resize the array without preserving data.
    template <typename... Extents, typename = std::enable_if_t<(sizeof...(Extents) == Rank)>>
    void resize(Extents... new_extents) {
        extents_ = { static_cast<std::size_t>(new_extents)... };
        std::size_t new_size = (static_cast<std::size_t>(new_extents) * ...);
        data_.resize(new_size);
        view_ = mdspan_type(data_.data(), make_dextents(extents_));
    }

    // Resize and preserve data in the overlapping region.
    // Data is preserved for indices in each dimension up to std::min(old, new) extent.
    template <typename... Extents, typename = std::enable_if_t<(sizeof...(Extents) == Rank)>>
    void resizeAndPreserve(Extents... new_extents) {
        std::array<std::size_t, Rank> newExtents = { static_cast<std::size_t>(new_extents)... };

        // Compute the product (total number of elements) for the new extents.
        std::size_t newTotal = std::accumulate(newExtents.begin(), newExtents.end(), size_t{1}, std::multiplies<>());

        // Prepare new storage.
        std::vector<T> newData(newTotal);

        // Compute the overlapping region extents.
        std::array<std::size_t, Rank> commonExtents;
        for (std::size_t i = 0; i < Rank; ++i) {
            commonExtents[i] = std::min(extents_[i], newExtents[i]);
        }

        // Compute total number of elements in the common region.
        std::size_t commonTotal = std::accumulate(commonExtents.begin(), commonExtents.end(), size_t{1}, std::multiplies<>());

        // For each element in the overlapping region, copy from old to new.
        for (std::size_t idx = 0; idx < commonTotal; ++idx) {
            // Convert linear index (in common region) to multi-index.
            auto multiIndex = linearToMultiIndex(idx, commonExtents);
            // Compute corresponding linear indices in the old and new arrays.
            std::size_t oldIndex = multiIndexToLinear(multiIndex, extents_);
            std::size_t newIndex = multiIndexToLinear(multiIndex, newExtents);
            newData[newIndex] = data_[oldIndex];
        }

        // Update the array.
        extents_ = newExtents;
        data_ = std::move(newData);
        view_ = mdspan_type(data_.data(), make_dextents(extents_));
    }

    // Element access operators (non-const and const) with exactly Rank indices.
    template <typename... Indices, typename = std::enable_if_t<(sizeof...(Indices) == Rank)>>
    T& operator()(Indices... indices) {
        return view_(static_cast<std::size_t>(indices)...);
    }

    template <typename... Indices, typename = std::enable_if_t<(sizeof...(Indices) == Rank)>>
    const T& operator()(Indices... indices) const {
        return view_(static_cast<std::size_t>(indices)...);
    }

    // Element access operators (non-const and const) using containter with exactly Rank size
    template <typename Container, typename = std::enable_if_t<Container::size() == Rank>>
    T& operator()(const Container& indices) {
        return (*this)(indices[0], indices[1]);  // Assuming Rank==2.
    }

    template <typename Container, typename = std::enable_if_t<Container::size() == Rank>>
    const T& operator()(const Container& indices) const {
        return (*this)(indices[0], indices[1]);  // Assuming Rank==2.
    }

    // Returns a pointer to the underlying contiguous storage.
    T* data() { return data_.data(); }
    const T* data() const { return data_.data(); }

    // Returns the current mdspan view.
    mdspan_type view() { return view_; }
    mdspan_type view() const { return view_; }

    // Returns the total number of elements.
    std::size_t size() const { return data_.size(); }

    // Returns the extents (dimensions) as a std::array.
    std::array<std::size_t, Rank> extents() const { return extents_; }

    // Fill the array with a specified value.
    void fill(const T& value) {
        std::fill(data_.begin(), data_.end(), value);
    }

    // Swap the contents with another DynamicArray.
    void swap(DynamicArray& other) noexcept {
        std::swap(extents_, other.extents_);
        data_.swap(other.data_);
        view_ = mdspan_type(data_.data(), make_dextents(extents_));
        other.view_ = mdspan_type(other.data_.data(), make_dextents(other.extents_));
    }

    // Iterator access to the underlying storage.
    auto begin() { return data_.begin(); }
    auto end() { return data_.end(); }
    auto begin() const { return data_.begin(); }
    auto end() const { return data_.end(); }

    // Fix the dimension FixedDim to a given index and return a lower-rank mdspan view.
    // Usage: For a 2D array, arr.slice<1>(n) returns a 1D mdspan corresponding to arr(:, n).
    template <std::size_t FixedDim>
    auto slice(std::size_t index) const {
        static_assert(FixedDim < Rank, "FixedDim out of range");
        // Compute strides for our current extents.
        const auto strides = compute_strides(extents_);
        // Check index validity.
        if (index >= extents_[FixedDim])
            throw std::out_of_range("slice index out of range");

        // Compute the offset in the linear data array.
        T* new_ptr = data_.data() + index * strides[FixedDim];

        // Build new extents by dropping the FixedDim.
        std::array<std::size_t, Rank - 1> newExtents{};
        for (std::size_t i = 0, j = 0; i < Rank; ++i) {
            if (i != FixedDim) {
                newExtents[j++] = extents_[i];
            }
        }
        // Create dextents for the new rank.
        auto new_dextents = make_dextents(newExtents);
        using new_mdspan_t = std::mdspan<T, std::dextents<std::size_t, Rank - 1>>;
        return new_mdspan_t(new_ptr, new_dextents);
    }

private:
    std::vector<T> data_;
    std::array<std::size_t, Rank> extents_{};
    mdspan_type view_{ nullptr, dextents_type{} };

    // Helper: convert an extents array to a std::dextents object.
    static dextents_type make_dextents(const std::array<std::size_t, Rank>& ext) {
        return make_dextents_impl(ext, std::make_index_sequence<Rank>{});
    }

    template <std::size_t... I>
    static dextents_type make_dextents_impl(const std::array<std::size_t, Rank>& ext, std::index_sequence<I...>) {
        return dextents_type{ ext[I]... };
    }

    // Overload for lower-rank extents.
    template <std::size_t R>
    static std::dextents<std::size_t, R> make_dextents(const std::array<std::size_t, R>& ext) {
        return make_dextents_impl_lower(ext, std::make_index_sequence<R>{});
    }

    template <std::size_t R, std::size_t... I>
    static std::dextents<std::size_t, R> make_dextents_impl_lower(const std::array<std::size_t, R>& ext, std::index_sequence<I...>) {
        return std::dextents<std::size_t, R>{ ext[I]... };
    }

    // Helper: Convert a multi-index (given as an array of indices) to a linear index,
    // assuming row-major order.
    static std::size_t multiIndexToLinear(const std::array<std::size_t, Rank>& indices,
                                            const std::array<std::size_t, Rank>& extents) {
        std::size_t linear = 0;
        for (std::size_t i = 0; i < Rank; ++i) {
            linear = linear * extents[i] + indices[i];
        }
        return linear;
    }

    // Helper: Convert a linear index into a multi-index given extents (row-major order).
    static std::array<std::size_t, Rank> linearToMultiIndex(std::size_t linear,
                                                            const std::array<std::size_t, Rank>& extents) {
        std::array<std::size_t, Rank> indices{};
        for (std::size_t i = Rank; i-- > 0; ) {
            indices[i] = linear % extents[i];
            linear /= extents[i];
        }
        return indices;
    }
};

// Helper: Cast all elements of a DynamicArray<T, Rank> to type To.
// This returns a new DynamicArray<To, Rank> with the same extents.
template<typename To, typename T, std::size_t Rank>
DynamicArray<To, Rank> castArray(const DynamicArray<T, Rank>& input) {
    // Get extents from the input (assumes extents() returns a std::array<size_t, Rank>)
    auto ext = input.extents();
    // Construct the output array using extents. For a 2D array:
    // (For higher ranks, you'll need to unpack ext appropriately.)
    if constexpr(Rank == 2) {
        DynamicArray<To, 2> output(ext[0], ext[1]);
        std::transform(input.begin(), input.end(), output.begin(),
                       [](const T &value) { return static_cast<To>(value); });
        return output;
    } else {
        // For higher ranks, implement similar logic.
        // You might iterate over all elements using the total size.
        DynamicArray<To, Rank> output;
        output.resize(ext[0], ext[1]); // Example for 2D; extend for other ranks.
        std::transform(input.begin(), input.end(), output.begin(),
                       [](const T &value) { return static_cast<To>(value); });
        return output;
    }
}

