#include <pbc/pbc.h>                     // Include PBC library for pairing

#include <iostream>
#include <random>
#include <regex>
#include <cassert>
#include <fstream>
#include <map> // Added for std::map
#include <set> // Include for std::set
#include <stdexcept> // Include for std::runtime_error
#include <string> // Include for std::hash, std::stoi
#include <cassert> // Include for assert

#include "gca/gca_tree.h" // Include GCATree header
#include "acc/keys.h"      // Include accumulator keys
#include "acc/accumulator.h" // Include accumulator base
#include "mht/blockchain.h"

// Use gca namespace for GCATree
using namespace gca;
using namespace std;

/*
 * TxnPool Implementation
 */

namespace mht {

TxnPool::TxnPool(unsigned int nSeedIn, string inPath, uint32_t VMax, int maxLines):
    nSeed(nSeedIn)
{
    ifstream infile(inPath);
    string sTxn, tmp;
    regex reg(R"(^(.*)\s\[(.*)\]\s\{(.*)\}$)");
    regex regAddr(R"(0x[0-9a-f]+)");
    smatch m, ma;
    uint32_t txid;
    uint32_t timestamp;
    uint32_t amount;
    vector<string> addresses;

    double varRand;
    uniform_real_distribution<double> uniformDist(0,1);

    if (infile.is_open()) {
        int lineCount = 0;
        uint32_t max_amount = 0;
        uint32_t min_amount = UINT32_MAX;
        while (getline(infile, sTxn) && (maxLines < 0 || lineCount < maxLines)) {
            // Parse transaction data

            txid = ComputeTxnId(sTxn);
            // txid = lineCount;

            regex_match(sTxn, m, reg);
            timestamp = stoi(m[1].str());
            amount = stoi(m[2].str());
            max_amount = max(max_amount, amount);
            min_amount = min(min_amount, amount);

            // Parse addresses
            assert(addresses.empty());
            string addrs = m[3].str();
            auto pos = addrs.cbegin();
            auto end = addrs.cend();
            for (; regex_search(pos, end, ma, regAddr); pos = ma.suffix().first) {
                addresses.emplace_back(ma.str());
            }
            if (addresses.size() == 1) {
                addresses.emplace_back(addresses[0]);
            }

            // Create transaction
            pool.emplace_back(Transaction(txid, timestamp, amount, addresses));
            addresses.clear();
            lineCount++;
        }
        assert(max_amount < VMax * min_amount);
        unit = min_amount;
        std::cout << "VMax: " << VMax << std::endl;
        std::cout << "min_amount: " << min_amount << std::endl;
        std::cout << "max_amount: " << max_amount << std::endl;
    } else {
         cerr << "Error opening transaction pool file: " << inPath << endl;
         // Decide if throwing an error is appropriate here
         // throw runtime_error("Failed to open transaction pool file.");
    }
}

// Get address from a specific transaction and position
string TxnPool::getAddr(unsigned int posTxn, unsigned int posAddr) const {
    if (posTxn < pool.size() && posAddr < pool[posTxn].addresses.size()) {
        return pool[posTxn].addresses[posAddr];
    }
    return ""; // Return empty string if indices are out of bounds
}

// Get the total number of transactions in the pool
int TxnPool::getSize() const {
    return pool.size();
}

// Compute a simple hash-based transaction ID
// NOTE: This uses std::hash, which is not cryptographically secure.
// For a real blockchain, a proper cryptographic hash (e.g., SHA-256) should be used.
uint32_t TxnPool::ComputeTxnId(string strTxn) {
    return static_cast<uint32_t>(hash<string>{}(strTxn));
}

/*
 * Block Implementation
 */

Block::Block(uint32_t timestampIn, const vector<Transaction>& vTxns, acc::AccPublicKey& pk, unsigned int nSeedIn, uint32_t VMax, uint32_t MaxId, uint32_t unit):
    timestamp(timestampIn),
    nSeed(nSeedIn),
    vTransactions(vTxns),
    gcaTree(make_unique<GCATree>(pk, vTxns, nSeedIn, VMax, MaxId, unit)){}

Block::~Block(){}

// Display block information (basic)
void Block::display() const {
    cout << "Block Timestamp: " << timestamp << endl;
    cout << "  Number of Transactions: " << vTransactions.size() << endl;
}

// Get the root of the Merkle Hash Tree over transactions
uint32_t Block::getRootHash() const {
    return gcaTree->getRootHash();
}

// Get the vector of transactions in the block
vector<Transaction> Block::getTxns() const {
    return vTransactions;
}

// Get the size of the ADS (placeholder - GCA² tree size is calculated elsewhere)
uint32_t Block::getSizeADS() const {
    // This should ideally return the size contributed by the GCA² tree for this block.
    // Placeholder: Could estimate based on aggrDigest structure if needed.
    return gcaTree->getSizeADS();
}

// Get the seed used for this block
unsigned int Block::getSeed() const {
    return nSeed;
}


/*
 * Chain Implementation
 */

// Chain::Chain(string key_path, string pbc_param_path, unsigned int nSeedIn, uint64_t universe_size, uint32_t VMax, uint32_t MaxId, uint32_t mergeThreIn):
//     nSeed(nSeedIn),
//     publicKey(acc::Accumulator::genkey(key_path, pbc_param_path, universe_size)),
//     mergeThreshold(mergeThreIn),
//     VMax(VMax),
//     MaxId(MaxId)
// {
//     cout << "Accumulator Keys Generated and Initialized." << endl;
// }

Chain::Chain(acc::AccPublicKey& pk, unsigned int nSeed, uint32_t VMax, uint32_t MaxId, uint32_t mergeThreIn):
    nSeed(nSeed),
    publicKey(pk),
    mergeThreshold(mergeThreIn),
    VMax(VMax),
    MaxId(MaxId)
{
    cout << "Accumulator Public Key Initialized." << endl;
}

Chain::~Chain() {
    cout << "Destroying Chain..." << endl;
    for (Block* block : chain) {
        delete block; // This will call Block destructor
    }
    chain.clear();

    cout << "Chain destroyed." << endl;
}

// Build the blockchain using Merkle Hash Trees and GCA² Trees
void Chain::buildChain(const TxnPool& txnPool) {
    if (txnPool.pool.empty()) {
        cout << "Transaction pool is empty. No blocks to build." << endl;
        return;
    }

    // Group transactions by timestamp
    map<uint32_t, vector<Transaction>> txnsByTimestamp;
    for (const auto& txn : txnPool.pool) {
        txnsByTimestamp[txn.timestamp].push_back(txn);
        // VMax = max(VMax, txn.amount);
    }
    // MaxId = txnPool.pool.size();

    cout << "Building chain for " << txnsByTimestamp.size() << " blocks..." << endl;

    // Create blocks for each timestamp group
    for (const auto& pair : txnsByTimestamp) {
        uint32_t timestamp = pair.first;
        const vector<Transaction>& blockTxns = pair.second;
        if (!blockTxns.empty()) {
            // cout << "Processing block for timestamp: " << timestamp << " with " << blockTxns.size() << " transactions." << endl;

            // Create a new block (builds GCA^2 Tree)
            Block* newBlock = new Block(timestamp, blockTxns, publicKey, nSeed, VMax, MaxId, txnPool.unit);

            chain.push_back(newBlock);
            // cout << "  Block added to chain." << endl;
        }
    }

    cout << "Finished building MHT chain with " << chain.size() << " blocks." << endl;
}

// Get a block by its timestamp (simple linear search)
Block* Chain::getBlock(uint32_t timestamp) const {
    for (Block* block : chain) {
        if (block->timestamp == timestamp) {
            return block;
        }
    }
    return nullptr; // Block not found
}

// Get the number of blocks in the chain
uint32_t Chain::getSize() const {
    return chain.size();
}

// Get the total size of ADSs in the chain (placeholder)
double Chain::getSizeADS() const {
    double totalSize = 0.0;
    for (const Block* block : chain) {
        totalSize += block->getSizeADS(); // Summing placeholder sizes for now
    }
    return totalSize;
}

// Display chain information (basic)
void Chain::display() const {
    cout << "Chain Information:" << endl;
    cout << "  Number of Blocks: " << chain.size() << endl;
    cout << "  Seed: " << nSeed << endl;
    // Optionally display info for each block
    // cout << "--- Blocks ---" << endl;
    // for (const Block* block : chain) {
    //     block->display();
    // }
    // cout << "--------------" << endl;
}

pair<uint64_t, VO> Chain::query(Query query) {
    VO vo(publicKey.pairing);
    uint64_t res = 0;
    vector<acc::Set> setDimResVec;

    uint32_t currentBlk = query.beginBlock;
    uint32_t endBlk = query.endBlock;
    assert(1 <= currentBlk && currentBlk <= endBlk && endBlk <= getSize());

    while (currentBlk <= endBlk) {
        if (mergeThreshold && currentBlk % mergeThreshold == 1 && currentBlk + mergeThreshold - 1 <= endBlk) {    // Search Combined MBF-tree
            currentBlk += mergeThreshold - 1;
            auto [setDimRes, proof] = getBlock(currentBlk)->queryMergeInter(query);

            vo.dimProofs.push(std::move(proof));
            setDimResVec.push_back(std::move(setDimRes));
            vo.accDimRes.push(acc::Accumulator::setup(setDimRes, publicKey));
        } else {
            // std::cout << "Chain.QuerySingleInter start on blk: " << currentBlk << std::endl;
            auto [setDimRes, proof] = getBlock(currentBlk)->querySingleInter(query);
            // std::cout << "Chain.QuerySingleInter finished." << std::endl;

            // std::cout << "setDimRes.size(): " << setDimRes.size() << std::endl;

            vo.dimProofs.push(std::move(proof));
            setDimResVec.push_back(std::move(setDimRes));
            vo.accDimRes.push(acc::Accumulator::setup(setDimRes, publicKey));
        }
        // VO.pop() loop
        currentBlk += 1;
    }

    std::cout << "Dim finished." << std::endl;

    // Nested Intersevtion
    // No intersection, since only single cid = 1 for attribute amount

    // Nested Union

    acc::Set setUnionRes = setDimResVec[0];
    // setUnionRes.print();
    for (uint32_t i = 1; i < setDimResVec.size(); i++) {
        // auto [res, proof] = acc::Accumulator::query_nested_union(setUnionRes, setDimResVec[i], publicKey);
        // vo.nestedUnionRes.push(res);
        // vo.nestedUnionProofs.push(proof);
        setUnionRes = setUnionRes.union_with(setDimResVec[i]);
        // setUnionRes.print();
    }
    // if (vo.nestedUnionRes.empty()) {
    //     vo.nestedUnionRes.push(acc::Accumulator::setup(setUnionRes, publicKey));
    // }

    // std::cout << "setUnionRes.size(): " << setUnionRes.size() << std::endl;
    vo.unionRes = acc::Accumulator::setup(setUnionRes, publicKey);
    
    std::cout << "Union finished." << std::endl;

    // Nested Range
    const size_t valueBits = std::ceil(std::log2(VMax));
    const size_t tidBits = std::ceil(std::log2(MaxId));

    uint64_t l = static_cast<uint64_t>(query.cid) << (valueBits + tidBits);
    uint64_t r = (static_cast<uint64_t>(query.cid + 1) << (valueBits + tidBits)) - 1;
    // std::cout << "l: " << l << std::endl;
    // std::cout << "r: " << r << std::endl;

    tie(vo.nestedRangeRes, vo.nestedRangeProof) = acc::Accumulator::query_nested_range(setUnionRes, l, r, publicKey);

    acc::Set setRangeRes;
    for (const auto& element : setUnionRes) {
        if (element >= l && element <= r) {
            setRangeRes.insert(element);
        }
    }

    // std::cout << "setRangeRes.size(): " << setRangeRes.size() << std::endl;

    std::cout << "Range finished." << std::endl;

    // Aggr
    switch (query.agg) {
        case AggType::COUNT: {
            acc::CountProof count_proof(publicKey.pairing);
            tie(res, count_proof) = acc::Accumulator::query_count(setRangeRes, publicKey);
            vo.aggrProof = std::move(count_proof);
            break;
        }
        case AggType::MAX: {
            acc::MaxProof max_proof(publicKey.pairing);
            tie(res, max_proof) = acc::Accumulator::query_max(setRangeRes, publicKey);
            vo.aggrProof = std::move(max_proof);

            // res = res >> tidBits;
            // res = res & ((1 << valueBits) - 1);
            break;
        }
        case AggType::MIN: {
            acc::MinProof min_proof(publicKey.pairing);
            tie(res, min_proof) = acc::Accumulator::query_min(setRangeRes, publicKey);
            vo.aggrProof = std::move(min_proof);

            // res = res >> tidBits;
            // res = res & ((1 << valueBits) - 1);
            break;
        }
        case AggType::SUM: {
            acc::SumProof sum_proof(publicKey.pairing);
            tie(res, sum_proof) = acc::Accumulator::query_sum(setRangeRes, publicKey);
            vo.aggrProof = std::move(sum_proof);
            break;
        }
        case AggType::AVG: {
            // AVG暂时注释掉
            break;
        }
    }

    std::cout << "Aggr finished." << std::endl;
    // std::cout << "res: " << res << std::endl;

    return pair<uint64_t, VO>{res, vo};
}

pair<acc::Set, DimProof> Block::querySingleInter(const Query& query) {
    auto result = gcaTree->querySingleInter(query);
    return result;
}

pair<acc::Set, DimProof> Block::queryMergeInter(const Query& query) {
    // TODO: for each gca query single
    return gcaTree->querySingleInter(query);
}


bool Chain::verify(Query query, uint64_t result, VO vo) {
    uint32_t currentBlk = query.beginBlock;
    uint32_t endBlk = query.endBlock;

    while (currentBlk <= endBlk) {
        if (mergeThreshold && currentBlk % mergeThreshold == 1 && currentBlk + mergeThreshold - 1 <= endBlk) {
            currentBlk += mergeThreshold - 1;
            if (!getBlock(currentBlk)->verifyMergeInter(query, vo.dimProofs.front())) {
                return false;
            }
        } else {
            if (!getBlock(currentBlk)->verifySingleInter(query, vo.dimProofs.front())) {
                return false;
            }
        }
        vo.dimProofs.pop();  // 移除队列中的第一个元素
        currentBlk += 1;
    }
    
    // Nested Intersection
    // No Intersection, since only single cid = 1 for attribute amount

    // Nested Union
    // auto res = vo.accDimRes.front();
    // vo.accDimRes.pop();
    // while (!vo.nestedUnionProofs.empty()) {
    //     auto l_acc = res;

    //     auto r_acc = vo.accDimRes.front();
    //     vo.accDimRes.pop();
        
    //     auto proof = vo.nestedUnionProofs.front();
    //     vo.nestedUnionProofs.pop();

    //     res = vo.nestedUnionRes.front();
    //     vo.nestedUnionRes.pop();

    //     if (!acc::Accumulator::verify_nested(l_acc, r_acc, proof, res, publicKey)) {
    //         return false;
    //     }
    // }
    
    auto res = vo.unionRes;

    // Nested Range
    const size_t valueBits = std::ceil(std::log2(VMax));
    const size_t tidBits = std::ceil(std::log2(MaxId));

    uint64_t l = static_cast<uint64_t>(query.cid) << (valueBits + tidBits);
    uint64_t r = (static_cast<uint64_t>(query.cid + 1) << (valueBits + tidBits)) - 1;
    // std::cout << "l: " << l << std::endl;
    // std::cout << "r: " << r << std::endl;

    if (!acc::Accumulator::verify_nested_range(res, vo.nestedRangeProof, vo.nestedRangeRes, l, r, publicKey)) {
        return false;
    }

    // Aggr
    bool aggr_verified = false;
    switch (query.agg) {
        case AggType::COUNT: {
            acc::CountProof& count_proof = std::get<acc::CountProof>(vo.aggrProof);
            aggr_verified = acc::Accumulator::verify_aggr(vo.nestedRangeRes, count_proof, result, publicKey);
            break;
        }
        case AggType::MAX: {
            acc::MaxProof& max_proof = std::get<acc::MaxProof>(vo.aggrProof);
            aggr_verified = acc::Accumulator::verify_aggr(vo.nestedRangeRes, max_proof, result, publicKey);
            break;
        }
        case AggType::MIN: {
            acc::MinProof& min_proof = std::get<acc::MinProof>(vo.aggrProof);
            aggr_verified = acc::Accumulator::verify_aggr(vo.nestedRangeRes, min_proof, result, publicKey);
            break;
        }
        case AggType::SUM: {
            acc::SumProof& sum_proof = std::get<acc::SumProof>(vo.aggrProof);
            aggr_verified = acc::Accumulator::verify_aggr(vo.nestedRangeRes, sum_proof, result, publicKey);
            break;
        }
        default:
            aggr_verified = false;
    }
    
    if (!aggr_verified) {
        std::cout << "**** Aggr verification failed. ****" << std::endl;
    }
    return aggr_verified;
}

bool Block::verifySingleInter(const Query& query, const gca::DimProof& proof) {
    auto result = gcaTree->verifySingleInter(query, proof);
    return result;
}

bool Block::verifyMergeInter(const Query& query, const gca::DimProof& proof) {
    // TODO: for each gca query single
    return gcaTree->verifySingleInter(query, proof);
}

} // namespace mht