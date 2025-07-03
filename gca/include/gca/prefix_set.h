/**
 * @file prefix_set.h
 * @brief Prefix set implementation for the GCA²-tree
 */

#ifndef GCA_PREFIX_SET_H
#define GCA_PREFIX_SET_H

#include <vector>

#include "acc/set.h"
#include "mht/transaction.h"
#include "gca/sorted_list.h"


namespace gca {

/**
 * @brief Represents a prefix set of transactions
 * 
 * A prefix set contains the first j transactions of a sorted list.
 */
class PrefixSet {
public:
    /**
     * @brief Constructor for PrefixSet
     * @param sortedList The sorted list to create the prefix set from
     * @param prefixSize The size of the prefix (number of transactions to include)
     */
    PrefixSet(const SortedList& sortedList, size_t prefixSize, uint32_t VMax, uint32_t MaxId, uint32_t unit);

    /**
     * @brief Get the transactions in the prefix set
     * @return The transactions
     */
    const std::vector<mht::Transaction>& getTransactions() const;

    /**
     * @brief Get the size of the prefix set
     * @return The size of the prefix set
     */
    size_t size() const;

    /**
     * @brief Get the prefix size
     * @return The prefix size
     */
    size_t getPrefixSize() const;

    /**
     * @brief Convert the prefix set to an accumulator set
     * @return The accumulator set
     */
    acc::Set toAccumulatorSet() const;

    /**
     * @brief Get the extended attribute values for all transactions in the prefix set
     * @return The extended attribute values
     */
    std::vector<uint64_t> getExtendedAttributeValues() const;

    /**
     * @brief Get the size of the prefix set
     * @return The size of the prefix set
     */
    size_t getSizeADS() const;

private:
    /**
     * @brief Calculate the extended attribute value for a transaction
     * @param transactionId The transaction ID
     * @return The extended attribute value
     */
    uint64_t calculateExtendedAttributeValue(uint32_t transactionId) const;

    size_t prefixSize;
    
    std::vector<mht::Transaction> transactions;
    
    uint32_t columnId;

    uint32_t VMax;
    uint32_t MaxId;
    uint32_t unit;
};

} // namespace gca

#endif // GCA_PREFIX_SET_H
