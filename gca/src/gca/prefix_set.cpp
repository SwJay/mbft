/**
 * @file prefix_set.cpp
 * @brief Implementation of the PrefixSet class
 */

#include <cmath>
#include <bitset>
#include <cassert>
#include <iostream>
#include "gca/prefix_set.h"

namespace gca {

PrefixSet::PrefixSet(const SortedList& sortedList, size_t prefixSize, uint32_t VMax, uint32_t MaxId, uint32_t unit):
    columnId(sortedList.getColumnId()),
    prefixSize(prefixSize),
    VMax(VMax),
    MaxId(MaxId),
    unit(unit) {
    // Get the transactions from the sorted list
    const std::vector<mht::Transaction>& allTransactions = sortedList.getTransactions();
    
    // Add the first prefixSize transactions to the prefix set
    for (size_t i = 0; i < prefixSize && i < allTransactions.size(); ++i) {
        transactions.push_back(allTransactions[i]);
    }
}

const std::vector<mht::Transaction>& PrefixSet::getTransactions() const {
    return transactions;
}

size_t PrefixSet::size() const {
    return transactions.size();
}

size_t PrefixSet::getPrefixSize() const {
    return prefixSize;
}

size_t PrefixSet::getSizeADS() const {
    size_t sizeADS = sizeof(columnId) + sizeof(prefixSize) + sizeof(VMax) + sizeof(MaxId);
    for (const mht::Transaction& transaction : transactions) {
        sizeADS += transaction.getSizeADS();
    }
    return sizeADS;
}

acc::Set PrefixSet::toAccumulatorSet() const {
    // Create an empty accumulator set
    acc::Set accSet;
    
    // Add the extended attribute values to the accumulator set
    std::vector<uint64_t> extendedValues = getExtendedAttributeValues();
    for (uint64_t value : extendedValues) {
        accSet.insert(value);
    }
    
    return accSet;
}

std::vector<uint64_t> PrefixSet::getExtendedAttributeValues() const {
    std::vector<uint64_t> extendedValues;
    
    // For each transaction in the prefix set
    for (int i = 0; i < transactions.size(); i++) {
        // Calculate the extended attribute value
        uint64_t extendedValue = calculateExtendedAttributeValue(i);
        
        extendedValues.push_back(extendedValue);
    }
    
    return extendedValues;
}

uint64_t PrefixSet::calculateExtendedAttributeValue(uint32_t transactionId) const {
    // Get the transaction
    const mht::Transaction& transaction = transactions[transactionId];

    // Get the attribute value
    uint32_t attributeValue = transaction.getAmount();
    attributeValue = attributeValue / unit;
    transactionId = transactionId % MaxId;
    
    // Numerical attribute is merely amount, cid=1
    assert(columnId == 0);

    // Calculate the extended attribute value
    // The extended value is (cid || value || tid)
    // where cid is the column ID, value is the attribute value, and tid is the transaction ID
    
    // Determine the number of bits for each component
    // const size_t cidBits = std::ceil(std::log2(d_t));
    const size_t valueBits = std::ceil(std::log2(VMax)); // Assuming 32-bit attribute values
    const size_t tidBits = std::ceil(std::log2(MaxId));
    
    // Combine the components into a single 64-bit value
    uint64_t extendedValue = 0;
    
    // Set the column ID bits
    // extendedValue |= static_cast<uint64_t>(columnId) << (valueBits + tidBits);
    // std::cout << "columnId: " << columnId << std::endl;
    
    // Set the attribute value bits
    extendedValue |= static_cast<uint64_t>(attributeValue) << tidBits;
    
    
    // Set the transaction ID bits
    extendedValue |= static_cast<uint64_t>(transactionId);
    // if (transactionId == transactions.size() - 1) {
    //     std::cout << "attributeValue: " << attributeValue << ", transactionId: "
    //     << transactionId <<  ", extendedValue: " << extendedValue << std::endl;
    // }
    
    return extendedValue;
}

} // namespace gca
