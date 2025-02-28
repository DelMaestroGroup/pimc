#include <vector>
#include <array>
#include <cstddef>
#include <mdspan>
#include <numeric>
#include <stdexcept>
#include <algorithm>
#include <utility>
#include <type_traits>

template <typename T, std::size_t Rank>
class DynamicArray {
public:
    using dextents_type = std::dextents<std::size_t, Rank>;
    using mdspan_type    = std::mdspan<T, dextents_type>;

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

