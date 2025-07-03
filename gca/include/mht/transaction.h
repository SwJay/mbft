#ifndef TRANSACTION_H
#define TRANSACTION_H

#include <cstdint>
#include <string>
#include <vector>

using std::string;
using std::vector;

namespace mht {

class Transaction{
public:
    Transaction() = default;
    Transaction(uint32_t txid, uint32_t timestamp, uint32_t amount, vector<string> addresses, string txnType = "default"):
        txid(txid),
        timestamp(timestamp),
        amount(amount),
        addresses(addresses),
        txnType(txnType) {};
    ~Transaction() = default;

    void display() const;

    uint32_t getTxid() const { return txid; }
    uint32_t getTimestamp() const { return timestamp; }
    uint32_t getAmount() const { return amount; }
    vector<string> getAddresses() const { return addresses; }
    string getTxnType() const { return txnType; }

    size_t getSizeADS() const;

    uint32_t txid;
    uint32_t timestamp;
    uint32_t amount;
    vector<string> addresses;
    string txnType; // Placeholder for transaction type
};

} // namespace mht

#endif //TRANSACTION_H
