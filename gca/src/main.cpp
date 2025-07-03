/**
 * @file main.cpp
 * @brief Main file for the GCA²-tree implementation demonstrating integrated build.
 */

#include <iostream>
#include <set>
#include <string>
#include <vector>      // Include vector for queryResult pair
#include <stdexcept>   // For exception handling
#include <limits>      // For numeric_limits
#include <random>      // For default_random_engine
#include <cassert>
#include <cmath>
#include <chrono>
#include <fstream> // Include fstream for ofstream
#include <filesystem> // Include filesystem for create_directories

#include "gca/gca_tree.h"
#include "mht/blockchain.h" // Needed for Chain, Block, TxnPool
#include "cmdline.h"

#define MAX_AMOUNT 8000000

// Use namespaces
using namespace gca;
using namespace acc;
using namespace mht;
using namespace std;

struct Params{
    Params() = default;
    Params(const cmdline::parser& argParser);

    // Property:
    unsigned int nSeed;
    uint32_t universe_bits;
    uint32_t vmax_bits;
    uint32_t maxid_bits;
    
    uint32_t max_lines;
    
    string txn_data_path;
    string pbc_param_path;
    string key_dir;

    // Query construction
    AggType agg;
    uint32_t amount_lower;
    uint32_t amount_upper;
    uint32_t begin_block;
    uint32_t end_block;

    // Test case
    int test_case;
    int repeats;
};

Params::Params(const cmdline::parser& argParser)
{
    // Parameter initialization
    // nSeed
    nSeed = argParser.get<int>("seed");
    if (nSeed == 0) {
        nSeed = time(0);
    }
    universe_bits = argParser.get<int>("universe-bits");
    vmax_bits = argParser.get<int>("vmax-bits");
    maxid_bits = argParser.get<int>("maxid-bits");

    max_lines = argParser.get<int>("max-lines");

    txn_data_path = argParser.get<string>("txn-data-path");
    pbc_param_path = argParser.get<string>("pbc-param-path");
    key_dir = argParser.get<string>("key-dir");

    // Query construction
    agg = (AggType)argParser.get<int>("agg-type");
    amount_lower = argParser.get<uint32_t>("lower");
    amount_upper = argParser.get<uint32_t>("upper");
    begin_block = argParser.get<uint32_t>("begin");
    end_block = argParser.get<uint32_t>("end");

    test_case = argParser.get<int>("test-case");
    repeats = argParser.get<int>("repeats");
}

Query genQuery(const Params& params, AggType agg, int window, double v, TxnPool& txnPool, int maxBlock) {
    string txnType = "default"; // Example transaction type
    int cid = 0; // 0 corresponds to 'amount'
    uint32_t beginBlock, endBlock;
    uint32_t lowerBound, upperBound;
    string addr;

    default_random_engine randEng(params.nSeed);
    int blockRange = maxBlock - window;
    int valueRange = (int)(MAX_AMOUNT * v);
    
    // Random time window
    beginBlock = randEng() % blockRange + 1;
    endBlock = beginBlock + window - 1;

    // Random value range
    lowerBound = randEng() % (MAX_AMOUNT - valueRange);
    upperBound = lowerBound + valueRange;

    cout << "query: lowerBound: " << lowerBound << ", upperBound: " << upperBound << ", beginBlock: " << beginBlock << ", endBlock: " << endBlock << endl;

    return Query(txnType, agg, cid, lowerBound, upperBound, beginBlock, endBlock);
}

void runGenKeys(const Params& params) {
    uint64_t universe_size = 1 << params.universe_bits;
    Accumulator::genkey(params.key_dir, params.pbc_param_path, universe_size);
}

void runSingleQuery(const Params& params) {
    assert(params.universe_bits == params.vmax_bits + params.maxid_bits);
    uint64_t universe_size = 1 << params.universe_bits;
    uint64_t VMax = 1 << params.vmax_bits;
    uint64_t MaxId = 1 << params.maxid_bits;

    // 1. Create Chain (this initializes keys and pairing)
    cout << "Creating Accumulator Public Key..." << endl;
    AccPublicKey pk = Accumulator::genkey(params.key_dir, params.pbc_param_path, universe_size);
    cout << "Accumulator Public Key Initialized." << endl;

    cout << "Creating Chain (Seed: " << params.nSeed << ")..." << endl;
    Chain chain(pk, params.nSeed, VMax, MaxId);
    cout << "Chain created. Accumulator Initialized." << endl;

    // 2. Create Transaction Pool
    cout << "Creating Transaction Pool (Data: " << params.txn_data_path << ")..." << endl;
    TxnPool txnPool(params.nSeed, params.txn_data_path, VMax, params.max_lines);
    cout << "Transaction Pool created with " << txnPool.getSize() << " transactions." << endl;

    // 3. Build the Chain (including GCA² trees for each block)
    cout << "Building Integrated Chain..." << endl;
    chain.buildChain(txnPool);
    cout << "Chain built with " << chain.getSize() << " blocks." << endl;

    // 4. Query the Chain
    string txnType = "default"; // Example transaction type
    int cid = 0; // 0 corresponds to 'amount'
    Query query(txnType, params.agg, cid, params.amount_lower, params.amount_upper, params.begin_block, params.end_block);

    // 7. Process the Query
    cout << "Processing query..." << endl;
    uint64_t result;
    VO vo(chain.publicKey.pairing);
    tie(result, vo) = chain.query(query);
    cout << "VO size: " << vo.getSizeADS() << endl;
    cout << "Query processed. Result: " << result << endl;

    // 8. Verify the Query Result (using the block's aggrDigest)
    cout << "Verifying query result..." << endl;
    bool isValid = chain.verify(query, result, vo);

    if (isValid) {
        cout << "Query verification successful." << endl;
        switch (query.agg) {
            case AggType::COUNT: {
                cout << "Result: " << result << endl;
                break;
            }
            case AggType::SUM:
            case AggType::MAX: 
            case AggType::MIN: {
                const size_t valueBits = std::ceil(std::log2(chain.VMax));
                const size_t tidBits = std::ceil(std::log2(chain.MaxId));
                result = result >> tidBits;
                // result = result & ((1 << valueBits) - 1);
                cout << "Result: " << result << endl;
                break;
            }
        }
    } else {
        cout << "Query verification failed." << endl;
    }

    cout << "--- Test Complete ---" << endl;
}


void varyWindow(const Params& params, AccPublicKey& pk) {
    assert(params.universe_bits == params.vmax_bits + params.maxid_bits);
    uint64_t universe_size = 1 << params.universe_bits;
    uint64_t VMax = 1 << params.vmax_bits;
    uint64_t MaxId = 1 << params.maxid_bits;

    vector<int> windowSizes = {400, 800, 1200, 1600, 2000};
    double value = 0.5;
    int repeats = params.repeats;
    AggType agg = AggType::COUNT;

    try {
        std::filesystem::create_directories("../res");
    } catch (const std::filesystem::filesystem_error& e) {
        std::cerr << "Error when creating directory '../res': " << e.what() << std::endl;
        throw;
    }

    string resPath = "../res/window.csv";
    ofstream resfile(resPath, ios::out);
    resfile << "Block Size,\tQuery Time(us),\tVO Size(B),\tVerify Time(us)" << endl;

    vector<double> queryTime(windowSizes.size(), 0);
    vector<double> voSize(windowSizes.size(), 0);
    vector<double> verifyTime(windowSizes.size(), 0);

    // 1. Create Chain (this initializes keys and pairing)

    cout << "Creating Chain (Seed: " << params.nSeed << ")..." << endl;
    Chain chain(pk, params.nSeed, VMax, MaxId);
    cout << "Chain created. Accumulator Initialized." << endl;

    // 2. Create Transaction Pool
    cout << "Creating Transaction Pool (Data: " << params.txn_data_path << ")..." << endl;
    TxnPool txnPool(params.nSeed, params.txn_data_path, VMax, params.max_lines);
    cout << "Transaction Pool created with " << txnPool.getSize() << " transactions." << endl;

    // 3. Build the Chain (including GCA² trees for each block)
    cout << "Building Integrated Chain..." << endl;
    chain.buildChain(txnPool);
    cout << "Chain built with " << chain.getSize() << " blocks." << endl;

    // Query the Chain
    for (int i = 0; i < windowSizes.size(); i++) {
        int window = windowSizes[i];
        for (int j = 0; j < repeats; j++) {
            Query query = genQuery(params, agg, window, value, txnPool, chain.getSize());
            // Process the Query
            cout << "Processing query... window: " << window << " repeat: " << j << endl;
            uint64_t result;
            VO vo(chain.publicKey.pairing);
            auto start = std::chrono::system_clock::now();
            tie(result, vo) = chain.query(query);
            auto end = std::chrono::system_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
            queryTime[i] += duration.count();
            voSize[i] += vo.getSizeADS();

            cout << "Verifying query result..." << endl;
            start = std::chrono::system_clock::now();
            bool isValid = chain.verify(query, result, vo);
            assert(isValid);
            end = std::chrono::system_clock::now();
            duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
            verifyTime[i] += duration.count();
        }
        queryTime[i] /= repeats;
        voSize[i] /= repeats;
        verifyTime[i] /= repeats;
        resfile << window << ",\t" << queryTime[i] << ",\t" << voSize[i] << ",\t" << verifyTime[i] << endl;
    }
    resfile.close();
    cout << "--- Test Complete ---" << endl;
}


void varyValue(const Params& params, AccPublicKey& pk) {
    assert(params.universe_bits == params.vmax_bits + params.maxid_bits);
    uint64_t universe_size = 1 << params.universe_bits;
    uint64_t VMax = 1 << params.vmax_bits;
    uint64_t MaxId = 1 << params.maxid_bits;

    int window = 800;
    vector<double> values = {0.1, 0.3, 0.5, 0.7, 0.9};
    int repeats = params.repeats;
    AggType agg = AggType::COUNT;

    try {
        std::filesystem::create_directories("../res");
    } catch (const std::filesystem::filesystem_error& e) {
        std::cerr << "Error when creating directory '../res': " << e.what() << std::endl;
        throw;
    }

    string resPath = "../res/value.csv";
    ofstream resfile(resPath, ios::out);
    resfile << "Value,\tQuery Time(us),\tVO Size(B),\tVerify Time(us)" << endl;

    vector<double> queryTime(values.size(), 0);
    vector<double> voSize(values.size(), 0);
    vector<double> verifyTime(values.size(), 0);

    cout << "Creating Chain (Seed: " << params.nSeed << ")..." << endl;
    Chain chain(pk, params.nSeed, VMax, MaxId);
    cout << "Chain created. Accumulator Initialized." << endl;

    // 2. Create Transaction Pool
    cout << "Creating Transaction Pool (Data: " << params.txn_data_path << ")..." << endl;
    TxnPool txnPool(params.nSeed, params.txn_data_path, VMax, params.max_lines);
    cout << "Transaction Pool created with " << txnPool.getSize() << " transactions." << endl;

    // 3. Build the Chain (including GCA² trees for each block)
    cout << "Building Integrated Chain..." << endl;
    chain.buildChain(txnPool);
    cout << "Chain built with " << chain.getSize() << " blocks." << endl;

    // Query the Chain
    for (int i = 0; i < values.size(); i++) {
        double value = values[i];
        for (int j = 0; j < repeats; j++) {
            Query query = genQuery(params, agg, window, value, txnPool, chain.getSize());
            // Process the Query
            cout << "Processing query... value: " << value << " repeat: " << j << endl;
            uint64_t result;
            VO vo(chain.publicKey.pairing);
            auto start = std::chrono::system_clock::now();
            tie(result, vo) = chain.query(query);
            auto end = std::chrono::system_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
            queryTime[i] += duration.count();
            voSize[i] += vo.getSizeADS();

            cout << "Verifying query result..." << endl;
            start = std::chrono::system_clock::now();
            bool isValid = chain.verify(query, result, vo);
            assert(isValid);
            end = std::chrono::system_clock::now();
            duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
            verifyTime[i] += duration.count();
        }
        queryTime[i] /= repeats;
        voSize[i] /= repeats;
        verifyTime[i] /= repeats;
        resfile << value << ",\t" << queryTime[i] << ",\t" << voSize[i] << ",\t" << verifyTime[i] << endl;
    }
    resfile.close();
    cout << "--- Test Complete ---" << endl;
}


void varyBlkSize(const Params& params, AccPublicKey& pk) {
    assert(params.universe_bits == params.vmax_bits + params.maxid_bits);
    uint64_t universe_size = 1 << params.universe_bits;
    uint64_t VMax = 1 << params.vmax_bits;
    uint64_t MaxId = 1 << params.maxid_bits;

    int window = 800;
    double value = 0.5;
    vector<int> blkSizes = {32, 64, 128, 256, 512};
    int repeats = params.repeats;
    AggType agg = AggType::COUNT;

    try {
        std::filesystem::create_directories("../res");
    } catch (const std::filesystem::filesystem_error& e) {
        std::cerr << "Error when creating directory '../res': " << e.what() << std::endl;
        throw;
    }

    string consPath = "../res/construct.csv";
    string resPath = "../res/blkSize.csv";
    ofstream consfile(consPath, ios::out);
    ofstream resfile(resPath, ios::out);
    consfile << "Block Size,\tTime(us),\tSize(B)" << endl;
    resfile << "Block Size,\tQuery Time(us),\tVO Size(B),\tVerify Time(us)" << endl;

    double constTime = 0;
    double constSize = 0;
    vector<double> queryTime(blkSizes.size(), 0);
    vector<double> voSize(blkSizes.size(), 0);
    vector<double> verifyTime(blkSizes.size(), 0);

    // 1. Create Chain (this initializes keys and pairing)

    vector<Chain> chains;
    cout << "Creating Chain (Seed: " << params.nSeed << ")..." << endl;
    for (int i = 0; i < blkSizes.size(); i++) { 
        Chain chain(pk, params.nSeed, VMax, MaxId);
        chains.push_back(chain);
    }
    cout << "Chain created. Accumulator Initialized." << endl;

    // 2. Create Transaction Pool
    vector<TxnPool> txnPools;
    
    for (int i = 0; i < chains.size(); i++) {
        string txn_data_path = "../data/mht/eth-" + to_string(blkSizes[i]) + ".dat";
        cout << "Creating Transaction Pool (Data: " << txn_data_path << ")..." << endl;
        TxnPool txnPool(params.nSeed, txn_data_path, VMax, params.max_lines);
        txnPools.push_back(txnPool);
    }
    cout << "Transaction Pool created with " << txnPools[0].getSize() << " transactions." << endl;

    // 3. Build the Chain (including GCA² trees for each block)
    cout << "Building Integrated Chain..." << endl;
    
    for (int i = 0; i < chains.size(); i++) {
        auto start = std::chrono::system_clock::now();
        chains[i].buildChain(txnPools[i]);
        auto end = std::chrono::system_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        constTime = duration.count() / chains[i].getSize();
        constSize = chains[i].getSizeADS() / chains[i].getSize();
        consfile << blkSizes[i] << ",\t" << constTime << ",\t" << constSize << endl;
        cout << "Chain built with " << chains[i].getSize() << " blocks." << endl;
    }
    
    // Query the Chain
    for (int i = 0; i < blkSizes.size(); i++) {
        int blkSize = blkSizes[i];
        for (int j = 0; j < repeats; j++) {
            Query query = genQuery(params, agg, window, value, txnPools[i], chains[i].getSize());
            // Process the Query
            cout << "Processing query... blkSize: " << blkSize << " repeat: " << j << endl;
            uint64_t result;
            VO vo(chains[i].publicKey.pairing);
            auto start = std::chrono::system_clock::now();
            tie(result, vo) = chains[i].query(query);
            auto end = std::chrono::system_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
            queryTime[i] += duration.count();
            voSize[i] += vo.getSizeADS();

            cout << "Verifying query result..." << endl;
            start = std::chrono::system_clock::now();
            bool isValid = chains[i].verify(query, result, vo);
            assert(isValid);
            end = std::chrono::system_clock::now();
            duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
            verifyTime[i] += duration.count();
        }
        queryTime[i] /= repeats;
        voSize[i] /= repeats;
        verifyTime[i] /= repeats;
        resfile << blkSize << ",\t" << queryTime[i] << ",\t" << voSize[i] << ",\t" << verifyTime[i] << endl;
    }
    resfile.close();
    cout << "--- Test Complete ---" << endl;
}


void varyAgg(const Params& params, AccPublicKey& pk) {
    assert(params.universe_bits == params.vmax_bits + params.maxid_bits);
    uint64_t universe_size = 1 << params.universe_bits;
    uint64_t VMax = 1 << params.vmax_bits;
    uint64_t MaxId = 1 << params.maxid_bits;

    int window = 800;
    double value = 0.5;
    int repeats = params.repeats;
    vector<AggType> aggs = {AggType::MAX, AggType::COUNT, AggType::SUM, AggType::COUNT};
    try {
        std::filesystem::create_directories("../res");
    } catch (const std::filesystem::filesystem_error& e) {
        std::cerr << "Error when creating directory '../res': " << e.what() << std::endl;
        throw;
    }

    string resPath = "../res/agg-sum.csv";
    ofstream resfile(resPath, ios::out);
    resfile << "Agg,\tQuery Time(us),\tVO Size(B),\tVerify Time(us)" << endl;

    vector<double> queryTime(aggs.size(), 0);
    vector<double> voSize(aggs.size(), 0);
    vector<double> verifyTime(aggs.size(), 0);

    cout << "Creating Chain (Seed: " << params.nSeed << ")..." << endl;
    Chain chain(pk, params.nSeed, VMax, MaxId);
    cout << "Chain created. Accumulator Initialized." << endl;

    // 2. Create Transaction Pool
    cout << "Creating Transaction Pool (Data: " << params.txn_data_path << ")..." << endl;
    TxnPool txnPool(params.nSeed, params.txn_data_path, VMax, params.max_lines);
    cout << "Transaction Pool created with " << txnPool.getSize() << " transactions." << endl;

    // 3. Build the Chain (including GCA² trees for each block)
    cout << "Building Integrated Chain..." << endl;

    chain.buildChain(txnPool);
    cout << "Chain built with " << chain.getSize() << " blocks." << endl;

    // Query Chain
    for (int i = 0; i < aggs.size(); i++) {
        AggType agg = aggs[i];
        for (int j = 0; j < repeats; j++) {
            Query query = genQuery(params, agg, window, value, txnPool, chain.getSize());
            // Process the Query
            cout << "Processing query... agg: " << i << " repeat: " << j << endl;
            uint64_t result;
            VO vo(chain.publicKey.pairing);
            auto start = std::chrono::system_clock::now();
            tie(result, vo) = chain.query(query);
            auto end = std::chrono::system_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
            queryTime[i] += duration.count();
            voSize[i] += vo.getSizeADS();

            cout << "Verifying query result..." << endl;
            start = std::chrono::system_clock::now();
            bool isValid = chain.verify(query, result, vo);
            assert(isValid);
            end = std::chrono::system_clock::now();
            duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
            verifyTime[i] += duration.count();
        }
        queryTime[i] /= repeats;
        voSize[i] /= repeats;
        verifyTime[i] /= repeats;
        resfile << i << ",\t" << queryTime[i] << ",\t" << voSize[i] << ",\t" << verifyTime[i] << endl;
    }
    
    resfile.close();
    cout << "--- Test Complete ---" << endl;
}

int main(int argc, char *argv[]) {
    cout << "Starting test..." << endl;

    cmdline::parser argParser;

    argParser.add<int>("seed", 's', "Seed for randomness in Hash & BF", false, 0);
    // argParser.add<int>("universe-bits", 'q', "Universe bits", false, 6);
    argParser.add<int>("universe-bits", 'q', "Universe bits", false, 14);
    // argParser.add<int>("vmax-bits", 'v', "Max value", false, 4);
    argParser.add<int>("vmax-bits", 'v', "Max value", false, 9);
    // argParser.add<int>("maxid-bits", 'i', "Max id", false, 2);
    argParser.add<int>("maxid-bits", 'i', "Max id", false, 5);

    argParser.add<int>("max-lines", 'm', "Max lines to read", false, -1);

    // Experiment
    argParser.add<string>("txn-data-path", 'd', "Path to transaction data", false, "../../data/eth-128.dat");
    argParser.add<string>("pbc-param-path", 'p', "Path to PBC parameters", false, "../param/a.param");
    argParser.add<string>("key-dir", 'k', "Path to the public key", false, "../data/key");

    // Query requirements
    argParser.add<int>("agg-type", 'a', "Query aggregation type: COUNT, SUM, MAX, MIN, AVG", false, 0, cmdline::range(0, 3));
    argParser.add<uint32_t>("lower", 'l', "Lower bound of amount", false, 4);
    argParser.add<uint32_t>("upper", 'u', "Upper bound of amount", false, 6);
    argParser.add<uint32_t>("begin", 'b', "Start block", false, 2);
    argParser.add<uint32_t>("end", 'e', "End block", false, 3);

    /* Test Case
     * 0. GenKeys:     Generate keys for the chain
     * 1. Single shot:  Build a chain, do single query & verification
     * 2. Vary window: test performance under verying time window
     * 3. Vary value range: test performance under varying value range
     * 4. Vary block size: test performance under varying block size
     * 5. Vary aggregate type: test performance under varying aggregate type
     * */
    argParser.add<int>("test-case", 'c', "Test cases", false, 0, cmdline::range(0, 6));
    argParser.add<int>("repeats", 'r', "Number of repeats", false, 100);

    argParser.parse_check(argc, argv);
    Params params(argParser);

    uint64_t universe_size = 1 << params.universe_bits;
    cout << "Creating Accumulator Public Key..." << endl;
    AccPublicKey pk = Accumulator::genkey(params.key_dir, params.pbc_param_path, universe_size);
    cout << "Accumulator Public Key Initialized." << endl;

    switch(params.test_case) {
        case 0:
            cout << "Case 0. GenKeys: Generate keys for the chain" << endl;
            runGenKeys(params);
            break;
        case 1:
            cout << "Case 1. Single shot: Build a chain, do single query & verification" << endl;
            runSingleQuery(params);
            break;
        case 2:
            cout << "Case 2. Vary Window: Build a chain, do single query & verification" << endl;
            varyWindow(params, pk);
            break;
        case 3:
            cout << "Case 3. Vary Value: Build a chain, do single query & verification" << endl;
            varyValue(params, pk);
            break;
        case 4:
            cout << "Case 4. Vary Block Size: Build a chain, do single query & verification" << endl;
            varyBlkSize(params, pk);
            break;
        case 5:
            cout << "Case 5. Vary Agg: Build a chain, do single query & verification" << endl;
            varyAgg(params, pk);
            break;
        default:
            break;
    }

    return 0;
}