/**
 * @file set.h
 * @brief Set class implementation for the Expressive Set Accumulator
 */

#ifndef ACC_SET_H
#define ACC_SET_H

#include <unordered_set>
#include <cstdint>
#include <cstddef>

namespace acc {

/**
 * @brief Represents a set of elements
 */
class Set {
private:
    std::unordered_set<uint64_t> elements;

public:
    // Constructors
    Set() = default;
    Set(const Set& other) = default;
    Set& operator=(const Set& other) = default;
    explicit Set(const std::unordered_set<uint64_t>& elems);

    // Print
    void print() const;
    
    // Set operations
    Set intersection(const Set& other) const;
    Set union_with(const Set& other) const;
    Set difference(const Set& other) const;
    bool is_subset_of(const Set& other) const;
    
    // Element operations
    void insert(uint64_t element);
    void remove(uint64_t element);
    bool contains(uint64_t element) const;
    
    // Iterators and size
    std::unordered_set<uint64_t>::const_iterator begin() const;
    std::unordered_set<uint64_t>::const_iterator end() const;
    size_t size() const;
    bool empty() const;
    
    // Comparison operators
    bool operator==(const Set& other) const;
    bool operator!=(const Set& other) const;
};

} // namespace acc

#endif // ACC_SET_H