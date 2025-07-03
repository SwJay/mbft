/**
 * @file sorted_list.h
 * @brief Sorted list implementation for the GCA²-tree
 */

#ifndef GCA_SORTED_LIST_H
#define GCA_SORTED_LIST_H

#include <vector>
#include <string>

#include "mht/transaction.h"

namespace gca {

/**
 * @brief Represents a sorted list of transactions for a specific attribute
 */
class SortedList {
public:
    /**
     * @brief Constructor for SortedList
     */
    SortedList(uint32_t columnId = 0);

    /**
     * @brief Add a transaction to the sorted list
     * @param transaction The transaction to add
     */
    void addTransaction(const mht::Transaction& transaction);

    /**
     * @brief Sort the transactions in the list
     */
    void sort();

    /**
     * @brief Get the transactions in the sorted list
     * @return The transactions
     */
    const std::vector<mht::Transaction>& getTransactions() const;

    /**
     * @brief Get the attribute name
     * @return The attribute name
     */
    const std::string& getAttributeName() const;

    /**
     * @brief Get the attribute type
     * @return The attribute's data type
     */
    const std::string& getAttributeType() const;

    /**
     * @brief Get the value of the attribute for a transaction at a specific index
     * @param index The index of the transaction
     * @return The value of the attribute
     */
    uint32_t getValueAt(size_t index) const;

    /**
     * @brief Get the transaction at a specific index
     * @param index The index of the transaction
     * @return The transaction
     */
    const mht::Transaction& getTransactionAt(size_t index) const;

    /**
     * @brief Get the size of the sorted list
     * @return The size of the sorted list
     */
    size_t size() const;

    /**
     * @brief Get the size of the sorted list
     * @return The size of the sorted list
     */
    size_t getSizeADS() const;

    /**
     * @brief Get the column ID for this sorted list
     * @return The column ID
     */
    uint32_t getColumnId() const;

private:
    uint32_t columnId;    
    std::vector<mht::Transaction> transactions;
};

} // namespace gca

#endif // GCA_SORTED_LIST_H
