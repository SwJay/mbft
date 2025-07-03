/**
 * @file set.cpp
 * @brief Implementation of the Set class for the Expressive Set Accumulator
 */

#include "acc/set.h"
#include <iostream>

namespace acc {

Set::Set(const std::unordered_set<uint64_t>& elems) : elements(elems) {}

Set Set::intersection(const Set& other) const {
    Set result;
    
    // Use the smaller set to iterate through
    const Set& smaller = (this->size() < other.size()) ? *this : other;
    const Set& larger = (this->size() < other.size()) ? other : *this;
    
    // Add elements from the smaller set that are also in the larger set
    for (const auto& elem : smaller.elements) {
        if (larger.contains(elem)) {
            result.insert(elem);
        }
    }
    
    return result;
}

void Set::print() const {
    for (const auto& elem : elements) {
        std::cout << elem << " ";
    }
    std::cout << std::endl;
}


Set Set::union_with(const Set& other) const {
    // Start with a copy of the larger set
    const Set& larger = (this->size() > other.size()) ? *this : other;
    const Set& smaller = (this->size() > other.size()) ? other : *this;
    
    Set result(larger.elements);
    
    // Add all elements from the smaller set
    for (const auto& elem : smaller.elements) {
        result.insert(elem);
    }
    
    return result;
}

Set Set::difference(const Set& other) const {
    Set result;
    
    // Add elements from this set that are not in the other set
    for (const auto& elem : this->elements) {
        if (!other.contains(elem)) {
            result.insert(elem);
        }
    }
    
    return result;
}

bool Set::is_subset_of(const Set& other) const {
    // Check if all elements in this set are also in the other set
    for (const auto& elem : this->elements) {
        if (!other.contains(elem)) {
            return false;
        }
    }
    
    return true;
}

void Set::insert(uint64_t element) {
    elements.insert(element);
}

void Set::remove(uint64_t element) {
    elements.erase(element);
}

bool Set::contains(uint64_t element) const {
    return elements.find(element) != elements.end();
}

std::unordered_set<uint64_t>::const_iterator Set::begin() const {
    return elements.begin();
}

std::unordered_set<uint64_t>::const_iterator Set::end() const {
    return elements.end();
}

size_t Set::size() const {
    return elements.size();
}

bool Set::empty() const {
    return elements.empty();
}

bool Set::operator==(const Set& other) const {
    return elements == other.elements;
}

bool Set::operator!=(const Set& other) const {
    return !(*this == other);
}

} // namespace acc 