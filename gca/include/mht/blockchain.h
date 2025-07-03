#ifndef BLOCKCHAIN_H
#define BLOCKCHAIN_H

#include <cstdint>
#include <vector>
#include <set>
#include <string>
#include <pbc/pbc.h>

#include "acc/accumulator.h"
#include "mht/merkle.h"
#include "mht/transaction.h"
#include "gca/gca_tree.h"
#include "gca/query_defs.h"

namespace mht {

class TxnPool {
public:
    TxnPool(unsigned int nSeed, std::string inPath, uint32_t VMax, int maxLines = -1);

    std::string getAddr(unsigned int posTxn, unsigned int posAddr) const;
    int getSize() const;
    uint32_t ComputeTxnId(std::string strTxn);

    // Property
    unsigned int nSeed;
    std::vector<Transaction> pool;
    uint32_t unit;
};

class Block {
public:
    Block() = default;
    Block(uint32_t timestamp, const std::vector<Transaction>& vTxns, acc::AccPublicKey& pk, unsigned int nSeed, uint32_t VMax, uint32_t MaxId, uint32_t unit);

    ~Block();

    std::pair<acc::Set, gca::DimProof> querySingleInter(const gca::Query& query);
    std::pair<acc::Set, gca::DimProof> queryMergeInter(const gca::Query& query);

    bool verifySingleInter(const gca::Query& query, const gca::DimProof& proof);
    bool verifyMergeInter(const gca::Query& query, const gca::DimProof& proof);

    // Getters
    uint32_t getRootHash() const;
    std::vector<Transaction> getTxns() const;
    uint32_t getSizeADS() const;
    unsigned int getSeed() const;

    void display() const;

    // Properties
    unsigned int nSeed;
    uint32_t timestamp;
    std::vector<Transaction> vTransactions;
    std::unique_ptr<gca::GCATree> gcaTree; // Placeholder for GCA² tree
};

class Chain {
public:
    Chain() = default; // Keep default constructor for potential use, but main one initializes keys

    // Chain(string key_path, string pbc_param_path, unsigned int nSeed, uint64_t universe_size, uint32_t VMax, uint32_t MaxId, uint32_t mergeThreshold=0);

    Chain(acc::AccPublicKey& pk, unsigned int nSeed, uint32_t VMax, uint32_t MaxId, uint32_t mergeThreshold=0);

    ~Chain();

    void buildChain(const TxnPool& txnPool);

    std::pair<uint64_t, gca::VO> query(gca::Query q);
    bool verify(gca::Query q, uint64_t result, gca::VO vo);

    // Getters
    Block* getBlock(uint32_t timestamp) const;
    uint32_t getSize() const;
    double getSizeADS() const;
    void display() const;

    // Property
    unsigned int nSeed;
    std::vector<Block*> chain;
    acc::AccPublicKey& publicKey; // Store the accumulator public key
    uint32_t mergeThreshold;

    uint32_t VMax = 0;
    uint32_t MaxId = 0;
};

} // namespace mht

#endif //BLOCKCHAIN_H
