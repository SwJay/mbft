#include <iostream>
#include <utility>

#include "mht/transaction.h"

namespace mht {

void Transaction::display() const
{
    std::cout << "------ Transaction::display() ------" << std::endl;
    std::cout << "txid: " << txid << std::endl;
    std::cout << "timestamp: " << timestamp << std::endl;
    std::cout << "amount: " << amount << std::endl;
    std::cout << "addresses: ";
    for (const string& address : addresses) {
        std::cout << address << " ";
    }
    std::cout << std::endl;
    std::cout << "txnType: " << txnType << std::endl;
    std::cout << "------------------------------------" << std::endl;
}

size_t Transaction::getSizeADS() const {
    size_t sizeADS = sizeof(txid) + sizeof(timestamp) + sizeof(amount);
    for (string addr : addresses) {
        sizeADS += addr.size();
    }
    return sizeADS;
}

} // namespace mht