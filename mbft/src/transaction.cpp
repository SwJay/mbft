#include "transaction.h"

#include <iostream>
#include <utility>

Transaction::Transaction(uint32_t txidIn, uint32_t timestampIn, uint32_t amountIn, vector<string> addressesIn):
    txid(txidIn),
    timestamp(timestampIn),
    amount(amountIn),
    addresses(move(addressesIn)),
    sketches(0)
{
}

Transaction::Transaction(uint32_t txidIn, uint32_t timestampIn, uint32_t amountIn, vector<string> addressesIn, vector<double> sketchesIn):
    txid(txidIn),
    timestamp(timestampIn),
    amount(amountIn),
    addresses(move(addressesIn)),
    sketches(move(sketchesIn))
{
}

void Transaction::display() const
{
    cout << "------ Transaction::display() ------" << endl;
    cout << "txid: " << txid << endl;
    cout << "timestamp: " << timestamp << endl;
    cout << "amount: " << amount << endl;
    cout << "addresses: ";
    for (const string& address : addresses) {
        cout << address << " ";
    }
    cout << endl;
    cout << "sketches: ";
    for (const double& sketch : sketches) {
        cout << sketch << " ";
    }
    cout << endl;
}

uint32_t Transaction::getTxid() const
{
    return txid;
}

uint32_t Transaction::getTimestamp() const
{
    return timestamp;
}

double Transaction::getValue(AggType aggType) const
{
    if (aggType) {  // sketches
        return sketches[aggType - 1];
    }
    return amount;
}

uint32_t Transaction::getAmount() const
{
    return amount;
}

vector<string> Transaction::getAddresses() const
{
    return addresses;
}

vector<double> Transaction::getSketches() const
{
    return sketches;
}

size_t Transaction::getSizeADS() const {
    size_t sizeADS = sizeof(txid) + sizeof(timestamp) + sizeof(amount);
    for (string addr : addresses) {
        sizeADS += addr.size();
    }
    return sizeADS;
}