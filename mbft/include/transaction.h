#ifndef TRANSACTION_H
#define TRANSACTION_H

//#include <cstdint>
//#include <set>
//#include <string>

#include "utils.h"

class Transaction{
public:
    Transaction() = default;
    Transaction(uint32_t txid, uint32_t timestamp, uint32_t amount, vector<string> addresses);
    Transaction(uint32_t txid, uint32_t timestamp, uint32_t amount, vector<string> addresses, vector<double> sketches);

    void display() const;

    uint32_t getTxid() const;
    uint32_t getTimestamp() const;
    double getValue(AggType aggType) const;
    uint32_t getAmount() const;
    vector<string> getAddresses() const;
    vector<double> getSketches() const;

    size_t getSizeADS() const;

    uint32_t txid;
    uint32_t timestamp;
    uint32_t amount;
    vector<string> addresses;
    vector<double> sketches;
};

#endif //TRANSACTION_H
