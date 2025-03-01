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

    // Element access operators
    // 1. Multi-index overload: expects exactly Rank separate indices.
    template <typename... Indices,
              typename = std::enable_if_t<(sizeof...(Indices) == Rank)>>
    T& operator()(Indices... indices) {
        std::size_t linear = computeLinearIndex(indices...);
        return data_[linear];
    }

    template <typename... Indices,
              typename = std::enable_if_t<(sizeof...(Indices) == Rank)>>
    const T& operator()(Indices... indices) const {
        std::size_t linear = computeLinearIndex(indices...);
        return data_[linear];
    }

    // 2. Container overload: accepts a container (like std::array) with Rank elements.
    template <typename Container,
              typename = std::enable_if_t<
                  std::tuple_size<std::decay_t<Container>>::value == Rank>>
    T& operator()(const Container& indices) {
        return callWithContainer(indices, std::make_index_sequence<Rank>{});
    }

    template <typename Container,
              typename = std::enable_if_t<
                  std::tuple_size<std::decay_t<Container>>::value == Rank>>
    const T& operator()(const Container& indices) const {
        return callWithContainer(indices, std::make_index_sequence<Rank>{});
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
    std::array<std::size_t, Rank> extents_{};
    std::vector<T> data_;
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

    // Helper: Compute the linear index from multi-dimensional indices.
    template <typename... Indices>
    std::size_t computeLinearIndex(Indices... indices) const {
        std::array<std::size_t, Rank> idx = { static_cast<std::size_t>(indices)... };
        return multiIndexToLinear(idx, extents_);
    }

    // Helper: Unpack container indices into individual indices (non-const overload).
    template <typename Container, std::size_t... I>
    T& callWithContainer(const Container& indices, std::index_sequence<I...>) {
        return (*this)(indices[I]...);
    }

    // Helper: Unpack container indices into individual indices (const overload).
    template <typename Container, std::size_t... I>
    const T& callWithContainer(const Container& indices, std::index_sequence<I...>) const {
        return (*this)(indices[I]...);
    }
};

// Specialization for bool.
template <std::size_t Rank>
class DynamicArray<bool, Rank> {
public:
    using dextents_type = std::dextents<std::size_t, Rank>;
    using mdspan_type   = std::mdspan<bool, dextents_type>;

    //=== Proxy for bool element access =======================================
    class BoolReference {
        unsigned char* ptr_;
    public:
        explicit BoolReference(unsigned char* ptr) : ptr_(ptr) {}
        operator bool() const { return *ptr_ != 0; }
        BoolReference& operator=(bool value) {
            *ptr_ = value ? 1 : 0;
            return *this;
        }
        BoolReference& operator=(const BoolReference& other) {
            *ptr_ = static_cast<bool>(other) ? 1 : 0;
            return *this;
        }
    };

    //=== Custom iterator types ===============================================
    class iterator {
    public:
        using iterator_category = std::random_access_iterator_tag;
        using value_type        = bool;
        using difference_type   = std::ptrdiff_t;
        using reference         = BoolReference;

        explicit iterator(unsigned char* p) : ptr_(p) {}
        reference operator*() const { return BoolReference(ptr_); }
        iterator& operator++() { ++ptr_; return *this; }
        iterator operator++(int) { iterator tmp(*this); ++ptr_; return tmp; }
        iterator& operator--() { --ptr_; return *this; }
        iterator operator--(int) { iterator tmp(*this); --ptr_; return tmp; }
        iterator& operator+=(difference_type n) { ptr_ += n; return *this; }
        iterator& operator-=(difference_type n) { ptr_ -= n; return *this; }
        iterator operator+(difference_type n) const { return iterator(ptr_ + n); }
        iterator operator-(difference_type n) const { return iterator(ptr_ - n); }
        difference_type operator-(const iterator& other) const { return ptr_ - other.ptr_; }
        bool operator==(const iterator& other) const { return ptr_ == other.ptr_; }
        bool operator!=(const iterator& other) const { return ptr_ != other.ptr_; }
        bool operator<(const iterator& other) const { return ptr_ < other.ptr_; }
        bool operator>(const iterator& other) const { return ptr_ > other.ptr_; }
        bool operator<=(const iterator& other) const { return ptr_ <= other.ptr_; }
        bool operator>=(const iterator& other) const { return ptr_ >= other.ptr_; }
        reference operator[](difference_type n) const { return BoolReference(ptr_ + n); }
    private:
        unsigned char* ptr_;
    };

    class const_iterator {
    public:
        using iterator_category = std::random_access_iterator_tag;
        using value_type        = bool;
        using difference_type   = std::ptrdiff_t;
        using reference         = bool;
        using pointer           = const bool*;

        explicit const_iterator(const unsigned char* p) : ptr_(p) {}
        reference operator*() const { return (*ptr_ != 0); }
        const_iterator& operator++() { ++ptr_; return *this; }
        const_iterator operator++(int) { const_iterator tmp(*this); ++ptr_; return tmp; }
        const_iterator& operator--() { --ptr_; return *this; }
        const_iterator operator--(int) { const_iterator tmp(*this); --ptr_; return tmp; }
        const_iterator& operator+=(difference_type n) { ptr_ += n; return *this; }
        const_iterator& operator-=(difference_type n) { ptr_ -= n; return *this; }
        const_iterator operator+(difference_type n) const { return const_iterator(ptr_ + n); }
        const_iterator operator-(difference_type n) const { return const_iterator(ptr_ - n); }
        difference_type operator-(const const_iterator& other) const { return ptr_ - other.ptr_; }
        bool operator==(const const_iterator& other) const { return ptr_ == other.ptr_; }
        bool operator!=(const const_iterator& other) const { return ptr_ != other.ptr_; }
        bool operator<(const const_iterator& other) const { return ptr_ < other.ptr_; }
        bool operator>(const const_iterator& other) const { return ptr_ > other.ptr_; }
        bool operator<=(const const_iterator& other) const { return ptr_ <= other.ptr_; }
        bool operator>=(const const_iterator& other) const { return ptr_ >= other.ptr_; }
        reference operator[](difference_type n) const { return *(ptr_ + n) != 0; }
    private:
        const unsigned char* ptr_;
    };

    using iterator_type = iterator;
    using const_iterator_type = const_iterator;

    //=== Constructors =========================================================
    template <typename... Extents,
              typename = std::enable_if_t<(sizeof...(Extents) == Rank) &&
                                            (std::conjunction_v<std::is_convertible<Extents, std::size_t>...>)>>
    DynamicArray(Extents... extents)
      : extents_{ static_cast<std::size_t>(extents)... },
        data_( (static_cast<std::size_t>(extents) * ...), 0 ),
        view_( reinterpret_cast<bool*>(data_.data()), make_dextents(extents_) )
    { }

    DynamicArray() = default;

    DynamicArray(const DynamicArray& other)
      : extents_(other.extents_),
        data_(other.data_),
        view_( reinterpret_cast<bool*>(data_.data()), make_dextents(extents_) )
    { }

    DynamicArray(DynamicArray&& other) noexcept
      : extents_(std::move(other.extents_)),
        data_(std::move(other.data_)),
        view_( reinterpret_cast<bool*>(data_.data()), make_dextents(extents_) )
    { }

    DynamicArray& operator=(const DynamicArray& other) {
        if (this != &other) {
            extents_ = other.extents_;
            data_ = other.data_;
            view_ = mdspan_type(reinterpret_cast<bool*>(data_.data()), make_dextents(extents_));
        }
        return *this;
    }

    DynamicArray& operator=(DynamicArray&& other) noexcept {
        if (this != &other) {
            extents_ = std::move(other.extents_);
            data_ = std::move(other.data_);
            view_ = mdspan_type(reinterpret_cast<bool*>(data_.data()), make_dextents(extents_));
        }
        return *this;
    }

    //=== Resize functions =====================================================
    template <typename... Extents, typename = std::enable_if_t<(sizeof...(Extents) == Rank)>>
    void resize(Extents... new_extents) {
        extents_ = { static_cast<std::size_t>(new_extents)... };
        std::size_t new_size = (static_cast<std::size_t>(new_extents) * ...);
        data_.resize(new_size, 0);
        view_ = mdspan_type(reinterpret_cast<bool*>(data_.data()), make_dextents(extents_));
    }

    template <typename... Extents, typename = std::enable_if_t<(sizeof...(Extents) == Rank)>>
    void resizeAndPreserve(Extents... new_extents) {
        std::array<std::size_t, Rank> newExtents = { static_cast<std::size_t>(new_extents)... };
        std::size_t newTotal = std::accumulate(newExtents.begin(), newExtents.end(), size_t{1}, std::multiplies<>());
        std::vector<unsigned char> newData(newTotal, 0);

        std::array<std::size_t, Rank> commonExtents;
        for (std::size_t i = 0; i < Rank; ++i)
            commonExtents[i] = std::min(extents_[i], newExtents[i]);

        std::size_t commonTotal = std::accumulate(commonExtents.begin(), commonExtents.end(), size_t{1}, std::multiplies<>());
        for (std::size_t idx = 0; idx < commonTotal; ++idx) {
            auto multiIndex = linearToMultiIndex(idx, commonExtents);
            std::size_t oldIndex = multiIndexToLinear(multiIndex, extents_);
            std::size_t newIndex = multiIndexToLinear(multiIndex, newExtents);
            newData[newIndex] = data_[oldIndex];
        }
        extents_ = newExtents;
        data_ = std::move(newData);
        view_ = mdspan_type(reinterpret_cast<bool*>(data_.data()), make_dextents(extents_));
    }

    //=== Element access =======================================================
    template <typename... Indices,
              typename = std::enable_if_t<(sizeof...(Indices) == Rank)>>
    BoolReference operator()(Indices... indices) {
        std::size_t linear = computeLinearIndex(indices...);
        return BoolReference(data_.data() + linear);
    }

    template <typename... Indices,
              typename = std::enable_if_t<(sizeof...(Indices) == Rank)>>
    bool operator()(Indices... indices) const {
        std::size_t linear = computeLinearIndex(indices...);
        return data_[linear] != 0;
    }

    template <typename Container,
              typename = std::enable_if_t<std::tuple_size<std::decay_t<Container>>::value == Rank>>
    auto operator()(const Container& indices) {
        return callWithContainer(indices, std::make_index_sequence<Rank>{});
    }

    template <typename Container,
              typename = std::enable_if_t<std::tuple_size<std::decay_t<Container>>::value == Rank>>
    auto operator()(const Container& indices) const {
        return callWithContainer(indices, std::make_index_sequence<Rank>{});
    }

    //=== Data and view access ================================================
    bool* data() { return reinterpret_cast<bool*>(data_.data()); }
    const bool* data() const { return reinterpret_cast<const bool*>(data_.data()); }

    mdspan_type view() { return view_; }
    mdspan_type view() const { return view_; }

    std::size_t size() const { return data_.size(); }
    std::array<std::size_t, Rank> extents() const { return extents_; }

    void fill(bool value) {
        std::fill(data_.begin(), data_.end(), value ? 1 : 0);
    }

    void swap(DynamicArray& other) noexcept {
        std::swap(extents_, other.extents_);
        data_.swap(other.data_);
        view_ = mdspan_type(reinterpret_cast<bool*>(data_.data()), make_dextents(extents_));
        other.view_ = mdspan_type(reinterpret_cast<bool*>(other.data_.data()), make_dextents(other.extents_));
    }

    //=== Iterator support =====================================================
    iterator begin() { return iterator(data_.data()); }
    iterator end() { return iterator(data_.data() + data_.size()); }
    const_iterator begin() const { return const_iterator(data_.data()); }
    const_iterator end() const { return const_iterator(data_.data() + data_.size()); }
    const_iterator cbegin() const { return const_iterator(data_.data()); }
    const_iterator cend() const { return const_iterator(data_.data() + data_.size()); }

    //=== Slice: Return a lower-rank mdspan view ===============================
    // Non-const overload.
    template <std::size_t FixedDim>
    auto slice(std::size_t index) -> std::mdspan<bool, std::dextents<std::size_t, Rank - 1>> {
        static_assert(FixedDim < Rank, "FixedDim out of range");
        const auto strides = compute_strides(extents_);
        if (index >= extents_[FixedDim])
            throw std::out_of_range("slice index out of range");
        bool* new_ptr = data() + index * strides[FixedDim];
        std::array<std::size_t, Rank - 1> newExtents{};
        for (std::size_t i = 0, j = 0; i < Rank; ++i) {
            if (i != FixedDim)
                newExtents[j++] = extents_[i];
        }
        auto new_dextents = make_dextents(newExtents);
        return std::mdspan<bool, std::dextents<std::size_t, Rank - 1>>(new_ptr, new_dextents);
    }

    // Const overload.
    template <std::size_t FixedDim>
    auto slice(std::size_t index) const -> std::mdspan<const bool, std::dextents<std::size_t, Rank - 1>> {
        static_assert(FixedDim < Rank, "FixedDim out of range");
        const auto strides = compute_strides(extents_);
        if (index >= extents_[FixedDim])
            throw std::out_of_range("slice index out of range");
        const bool* new_ptr = data() + index * strides[FixedDim];
        std::array<std::size_t, Rank - 1> newExtents{};
        for (std::size_t i = 0, j = 0; i < Rank; ++i) {
            if (i != FixedDim)
                newExtents[j++] = extents_[i];
        }
        auto new_dextents = make_dextents(newExtents);
        return std::mdspan<const bool, std::dextents<std::size_t, Rank - 1>>(new_ptr, new_dextents);
    }

private:
    std::array<std::size_t, Rank> extents_{};
    std::vector<unsigned char> data_;
    mdspan_type view_{ nullptr, dextents_type{} };

    //=== Helpers ==============================================================
    static dextents_type make_dextents(const std::array<std::size_t, Rank>& ext) {
        return make_dextents_impl(ext, std::make_index_sequence<Rank>{});
    }
    template <std::size_t... I>
    static dextents_type make_dextents_impl(const std::array<std::size_t, Rank>& ext, std::index_sequence<I...>) {
        return dextents_type{ ext[I]... };
    }
    template <std::size_t R>
    static std::dextents<std::size_t, R> make_dextents(const std::array<std::size_t, R>& ext) {
        return make_dextents_impl_lower(ext, std::make_index_sequence<R>{});
    }
    template <std::size_t R, std::size_t... I>
    static std::dextents<std::size_t, R> make_dextents_impl_lower(const std::array<std::size_t, R>& ext, std::index_sequence<I...>) {
        return std::dextents<std::size_t, R>{ ext[I]... };
    }
    static std::size_t multiIndexToLinear(const std::array<std::size_t, Rank>& indices,
                                            const std::array<std::size_t, Rank>& extents) {
        std::size_t linear = 0;
        for (std::size_t i = 0; i < Rank; ++i)
            linear = linear * extents[i] + indices[i];
        return linear;
    }
    static std::array<std::size_t, Rank> linearToMultiIndex(std::size_t linear,
                                                            const std::array<std::size_t, Rank>& extents) {
        std::array<std::size_t, Rank> indices{};
        for (std::size_t i = Rank; i-- > 0; ) {
            indices[i] = linear % extents[i];
            linear /= extents[i];
        }
        return indices;
    }
    template <typename... Indices>
    std::size_t computeLinearIndex(Indices... indices) const {
        std::array<std::size_t, Rank> idx = { static_cast<std::size_t>(indices)... };
        return multiIndexToLinear(idx, extents_);
    }
    template <typename Container, std::size_t... I>
    auto callWithContainer(const Container& indices, std::index_sequence<I...>) {
        return (*this)(indices[I]...);
    }
    template <typename Container, std::size_t... I>
    auto callWithContainer(const Container& indices, std::index_sequence<I...>) const {
        return (*this)(indices[I]...);
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

