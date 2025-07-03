/**
 * @file sorted_list.cpp
 * @brief Implementation of the SortedList class
 */

#include <algorithm>
#include <cassert>
#include <stdexcept>

#include "gca/sorted_list.h"

namespace gca {

SortedList::SortedList(uint32_t columnId): columnId(columnId) {
}

void SortedList::addTransaction(const mht::Transaction& transaction) {
    transactions.push_back(transaction);
}

void SortedList::sort() {
    std::sort(transactions.begin(), transactions.end(),
        [this](const mht::Transaction& a, const mht::Transaction& b) {
            // Get the value of the attribute for each transaction
            uint32_t valueA = 0;
            uint32_t valueB = 0;

            // Numerical attribute is merely amount, cid=0
            assert(columnId == 0);
            
            valueA = a.getAmount();
            valueB = b.getAmount();

            return valueA < valueB;
        });
}

const std::vector<mht::Transaction>& SortedList::getTransactions() const {
    return transactions;
}

uint32_t SortedList::getColumnId() const {
    return columnId;
}

uint32_t SortedList::getValueAt(size_t index) const {
    if (index >= transactions.size()) {
        return 0;
    }

    // Numerical attribute is merely amount, cid=0
    assert(columnId == 0);
    
    return transactions[index].getAmount();
}

size_t SortedList::size() const {
    return transactions.size();
}

size_t SortedList::getSizeADS() const {
    size_t sizeADS = sizeof(columnId);
    for (const mht::Transaction& transaction : transactions) {
        sizeADS += transaction.getSizeADS();
    }
    return sizeADS;
}

const mht::Transaction& SortedList::getTransactionAt(size_t index) const {
    if (index >= transactions.size()) {
        throw std::out_of_range("Index out of range");
    }
    return transactions[index];
}

} // namespace gca
