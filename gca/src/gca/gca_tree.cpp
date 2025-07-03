/**
 * @file gca_tree.cpp
 * @brief Implementation of the GCATree class
 */

#include <algorithm>
#include <unordered_set>
#include <iostream>
#include <string>
#include <cassert>

#include "mht/merkle.h"
#include "mht/transaction.h"
#include "gca/gca_tree.h"
#include "gca/triplet.h"
#include "gca/prefix_set.h"

namespace gca {

GCATree::GCATree(acc::AccPublicKey& pk, const std::vector<mht::Transaction>& transactions, unsigned int nSeed, uint32_t VMax, uint32_t MaxId, uint32_t unit):
    publicKey(pk),
    nSeed(nSeed)
{
    // Group transactions by type. We now use a single "default" type,
    std::unordered_map<std::string, std::vector<mht::Transaction>> transactionsByType;

    // All transactions go to default type
    for (const mht::Transaction& transaction : transactions) {
        transactionsByType[transaction.getTxnType()].push_back(transaction);
    }
    
    // Collect all hashes from all transaction types
    std::vector<uint32_t> allHashes;
    std::vector<std::string> processedTypes;
    
    // For each transaction type
    for (const auto& pair : transactionsByType) {
        const std::string& transactionType = pair.first;
        const std::vector<mht::Transaction>& typeTransactions = pair.second;
        
        // Create sorted lists for this transaction type
        std::vector<SortedList> sortedLists = createSortedLists(typeTransactions, transactionType);
        sortedListsMap[transactionType] = sortedLists;
        
        // Generate prefix sets from the sorted lists
        std::vector<std::vector<PrefixSet>> prefixSets = generatePrefixSets(sortedLists, VMax, MaxId, unit);
        prefixSetsMap[transactionType] = prefixSets;

        
        // Calculate digests for the prefix sets
        std::vector<std::vector<acc::AccValue>> digests = calculateDigests(prefixSets);
        digestsMap[transactionType] = digests;

        // Construct triplets from the digests
        std::vector<std::vector<Triplet>> triplets = constructTriplets(sortedLists, prefixSets, digests);
        tripletsMap[transactionType] = triplets;

        // Add this type to the list of processed types
        processedTypes.push_back(transactionType);

        // Collect hashes from this type's triplets
        for (uint32_t i = 0; i < triplets.size(); i++) {
            for (uint32_t j = 0; j < triplets[i].size(); j++) {
                allHashes.push_back(hash(transactionType, i, j, triplets[i][j], nSeed));
            }
        }
    }
    
    root = new mht::MerkleNode(allHashes, nSeed);
    sizeADS += root->getSizeADS();
}

GCATree::~GCATree() {
    delete root;
}

uint32_t GCATree::getSizeADS() {
    return sizeADS;
}

std::vector<SortedList> GCATree::createSortedLists(const std::vector<mht::Transaction>& transactions, const std::string& transactionType) {
    std::vector<SortedList> sortedLists;
    
    // Create sorted lists for queryable attributes defined in the blockchain
    const std::vector<uint32_t>& queryCids{0}; // single cid = 0 for attribtue amount
    
    for (const auto& cid : queryCids) {
        SortedList list(cid);
        for (const mht::Transaction& transaction : transactions) {
            list.addTransaction(transaction);
        }
        list.sort();
        sortedLists.push_back(list);
        sizeADS += list.getSizeADS();
    }
    
    return sortedLists;
}

std::vector<std::vector<PrefixSet>> GCATree::generatePrefixSets(const std::vector<SortedList>& sortedLists, uint32_t VMax, uint32_t MaxId, uint32_t unit) {
    std::vector<std::vector<PrefixSet>> prefixSets;
    
    // For each sorted list
    for (const SortedList& sortedList : sortedLists) {
        std::vector<PrefixSet> listPrefixSets;
        
        // Create prefix sets of increasing size
        for (size_t i = 1; i <= sortedList.size(); ++i) {
            PrefixSet prefixSet(sortedList, i, VMax, MaxId, unit);
            listPrefixSets.push_back(prefixSet);
            sizeADS += prefixSet.getSizeADS();
        }
        
        prefixSets.push_back(listPrefixSets);
    }
    
    return prefixSets;
}

std::vector<std::vector<acc::AccValue>> GCATree::calculateDigests(const std::vector<std::vector<PrefixSet>>& prefixSets) {
    std::vector<std::vector<acc::AccValue>> digests;
    
    // For each list of prefix sets
    for (const std::vector<PrefixSet>& listPrefixSets : prefixSets) {
        std::vector<acc::AccValue> listDigests;
        
        // Calculate a digest for each prefix set
        for (const PrefixSet& prefixSet : listPrefixSets) {
            // Convert the prefix set to an accumulator set
            acc::Set accSet = prefixSet.toAccumulatorSet();
            
            // Calculate the digest using the accumulator
            acc::AccValue digest = acc::AccValue::setup(accSet, publicKey);
            sizeADS += digest.getSizeADS();

            listDigests.push_back(digest);
        }
        
        digests.push_back(listDigests);
    }
    
    return digests;
}

std::vector<std::vector<Triplet>> GCATree::constructTriplets(
    const std::vector<SortedList>& sortedLists,
    const std::vector<std::vector<PrefixSet>>& prefixSets,
    const std::vector<std::vector<acc::AccValue>>& digests) {
    std::vector<std::vector<Triplet>> triplets;
    
    // For each sorted list, prefix sets, and digests
    for (size_t i = 0; i < sortedLists.size(); ++i) {
        const SortedList& sortedList = sortedLists[i];
        const std::vector<PrefixSet>& listPrefixSets = prefixSets[i];
        const std::vector<acc::AccValue>& listDigests = digests[i];
        
        std::vector<Triplet> listTriplets;
        
        // Construct a triplet for each prefix set and digest
        for (size_t j = 0; j < listPrefixSets.size(); ++j) {
            // Get the j-th value from the sorted list
            uint32_t value = sortedList.getValueAt(j);
            
            // Get the (j+1)-th value from the sorted list
            // If j is the last index, use a sentinel value
            uint32_t nextValue = (j + 1 < sortedList.size()) ? sortedList.getValueAt(j + 1) : UINT32_MAX;
            
            // Get the digest for the prefix set
            const acc::AccValue& digest = listDigests[j];
            
            // Create a triplet
            Triplet triplet(value, nextValue, digest);
            sizeADS += triplet.getSizeADS();
            listTriplets.push_back(triplet);
        }
        
        triplets.push_back(listTriplets);
    }
    
    return triplets;
}

uint32_t GCATree::hash(std::string t, uint32_t cid, uint32_t j, Triplet triplet, const unsigned int nSeed) const {
    string cid_str = std::to_string(cid);
	string j_str = std::to_string(j);
    string trip_str = std::to_string(triplet.hash(nSeed));
    
	string str = t + cid_str + j_str + trip_str;
	vector<unsigned char> vec(str.begin(), str.end());
	uint32_t hash;

    MurmurHash3_x86_32(vec.data(), vec.size(), nSeed, &hash);

    return hash;
}

uint32_t GCATree::getRootHash() const {
    return root->getHash();
}

std::tuple<int, int, DimBoundType> GCATree::findBounds(const SortedList& sortedList, uint32_t low, uint32_t high) const {
    assert(sortedList.size() > 0 && "Input array must not be empty");
    
    // Check if max element is < low
    if (sortedList.getValueAt(sortedList.size() - 1) < low) {
        return {sortedList.size(), sortedList.size(), DimBoundType::MAX_BELOW_LOW};
    }
    
    // Check if min element is > high
    if (sortedList.getValueAt(0) > high) {
        return {-1, -1, DimBoundType::MIN_ABOVE_UP};
    }
    
    DimBoundType boundType;

    // Check corner case: first element >= low
    int lowIndex;
    if (sortedList.getValueAt(0) >= low) {
        lowIndex = -1;
        boundType = DimBoundType::MIN_INCLUDED;
    } else {
        // Find lowIndex (last element < low)
        int left = 0;
        int right = sortedList.size() - 1;
        lowIndex = -1;
        
        while (left <= right) {
            int mid = left + (right - left) / 2;
            
            if (sortedList.getValueAt(mid) < low) {
                lowIndex = mid;
                left = mid + 1;
            } else {
                right = mid - 1;
            }
        }
        boundType = DimBoundType::IN_RANGE;
    }
    
    // Find highIndex (last element <= high)
    int left = 0;
    int right = sortedList.size() - 1;
    int highIndex = -1;
    
    while (left <= right) {
        int mid = left + (right - left) / 2;
        
        if (sortedList.getValueAt(mid) <= high) {
            highIndex = mid;
            left = mid + 1;
        } else {
            right = mid - 1;
        }
    }
    
    return {lowIndex, highIndex, boundType};
}

std::pair<acc::Set, DimProof> GCATree::querySingleInter(const Query& query) {
    acc::AccValue NestDiffRes(publicKey.pairing);
    DimProof proof(publicKey.pairing);
    acc::Set setRes;
    
    // Verify we have sorted lists for the queried attribute
    auto sortedListsIt = sortedListsMap.find(query.TxnType);
    if (sortedListsIt == sortedListsMap.end() || query.cid > sortedListsIt->second.size()) {
        throw std::runtime_error("Invalid attribute index");
    }

    // Get the relevant sorted list, prefix sets, and triplets
    const SortedList& sortedList = sortedListsIt->second[query.cid];
    const std::vector<PrefixSet>& prefixSets = prefixSetsMap[query.TxnType][query.cid];
    const std::vector<acc::AccValue>& digests = digestsMap[query.TxnType][query.cid];
    const std::vector<Triplet>& triplets = tripletsMap[query.TxnType][query.cid];

    // Find the bounds of the query
    std::tuple<int, int, DimBoundType> bounds = findBounds(sortedList, query.lowerBound, query.upperBound);
    proof.lowIndex = std::get<0>(bounds);
    proof.highIndex = std::get<1>(bounds);
    proof.boundType = std::get<2>(bounds);

    // std::cout << "proof.boundType: " << static_cast<int>(proof.boundType) << std::endl;

    switch (proof.boundType) {
        case DimBoundType::MAX_BELOW_LOW: {
            assert(proof.lowIndex == sortedList.size() && proof.highIndex == sortedList.size());
            proof.tripletUp = triplets[sortedList.size() - 1];
            proof.MHPathUp = root->getMerklePath(sortedList.size() - 1);

            break;
        }
        case DimBoundType::MIN_ABOVE_UP: {
            assert(proof.lowIndex == -1 && proof.highIndex == -1);
            proof.tripletLow = triplets[0];
            proof.MHPathLow = root->getMerklePath(0);

            break;
        }
        case DimBoundType::MIN_INCLUDED: {
            assert(proof.lowIndex == -1 && proof.lowIndex != proof.highIndex);
            proof.tripletLow = triplets[0];
            proof.tripletUp = triplets[proof.highIndex];
            proof.MHPathLow = root->getMerklePath(0);
            proof.MHPathUp = root->getMerklePath(proof.highIndex);

            setRes = prefixSets[proof.highIndex].toAccumulatorSet();

            break;
        }
        case DimBoundType::IN_RANGE: {
            if(!(proof.lowIndex != -1 && proof.lowIndex != proof.highIndex)){
                std::cout << "lowindex " << proof.lowIndex << std::endl;
                std::cout << "highindex " << proof.highIndex << std::endl;
            }
            // assert(proof.lowIndex != -1 && proof.lowIndex != proof.highIndex);
            proof.tripletLow = triplets[proof.lowIndex];
            proof.tripletUp = triplets[proof.highIndex];
            proof.MHPathLow = root->getMerklePath(proof.lowIndex);
            proof.MHPathUp = root->getMerklePath(proof.highIndex);
            std::tie(proof.NestDiffRes, proof.NestDiffProof) = acc::Accumulator::query_nested_difference(
                prefixSets[proof.lowIndex].toAccumulatorSet(), prefixSets[proof.highIndex].toAccumulatorSet(), publicKey);

            setRes = prefixSets[proof.highIndex].toAccumulatorSet().difference(prefixSets[proof.lowIndex].toAccumulatorSet());

            // std::cout << "GCA.setRes.size(): " << setRes.size() << std::endl;

            break;
        }
    }

    // 使用已初始化的proof对象
    // auto res =  std::pair<acc::Set, DimProof>(setRes, proof);
    return std::pair<acc::Set, DimProof>(setRes, std::move(proof));
}

bool GCATree::verifySingleInter(const Query& query, const DimProof& proof) {
    // Verify based on the bound type from the proof
    switch (proof.boundType) {
        case DimBoundType::MAX_BELOW_LOW: {            
            // Verify the maximum value is less than the lower bound
            if (proof.tripletUp.getValue() >= query.lowerBound) {
                return false;
            }
            
            // Calculate the hash of the triplet and verify the Merkle path
            uint32_t tripletHash = hash(query.TxnType, query.cid, proof.highIndex - 1, proof.tripletUp, nSeed);

            // Verify the Merkle path for the triplet
            return root->VerifyMerklePath(tripletHash, proof.MHPathUp, nSeed, root->getHash());
        }
        
        case DimBoundType::MIN_ABOVE_UP: {            
            // Verify the minimum value is greater than the upper bound
            if (proof.tripletLow.getValue() <= query.upperBound) {
                return false;
            }
            
            // Calculate the hash of the triplet and verify the Merkle path
            uint32_t tripletHash = hash(query.TxnType, query.cid, proof.lowIndex + 1, proof.tripletLow, nSeed);
            
            // Verify the Merkle path for the triplet
            return root->VerifyMerklePath(tripletHash, proof.MHPathLow, nSeed, root->getHash());
        }
        
        case DimBoundType::MIN_INCLUDED: {
            if (proof.tripletLow.getValue() < query.lowerBound) {
                return false;
            }

            if (proof.tripletUp.getValue() > query.upperBound || proof.tripletUp.getNextValue() <= query.upperBound) {
                return false;
            }
            
            // Calculate and verify the hashes for both triplets
            uint32_t lowTripletHash = hash(query.TxnType, query.cid, proof.lowIndex + 1, proof.tripletLow, nSeed);
            uint32_t upTripletHash = hash(query.TxnType, query.cid, proof.highIndex, proof.tripletUp, nSeed);
            
            // Verify both Merkle paths
            bool lowPathValid = root->VerifyMerklePath(lowTripletHash, proof.MHPathLow, nSeed, root->getHash());
            bool upPathValid = root->VerifyMerklePath(upTripletHash, proof.MHPathUp, nSeed, root->getHash());
            
            return lowPathValid && upPathValid;
        }
        
        case DimBoundType::IN_RANGE: {
            if (proof.tripletLow.getValue() >= query.lowerBound || proof.tripletLow.getNextValue() < query.lowerBound) {
                return false;
            }

            if (proof.tripletUp.getValue() > query.upperBound || proof.tripletUp.getNextValue() <= query.upperBound) {
                return false;
            }
            
            // Calculate and verify hashes for both triplets
            uint32_t lowTripletHash = hash(query.TxnType, query.cid, proof.lowIndex, proof.tripletLow, nSeed);
            uint32_t upTripletHash = hash(query.TxnType, query.cid, proof.highIndex, proof.tripletUp, nSeed);
            
            // Verify Merkle paths
            bool lowPathValid = root->VerifyMerklePath(lowTripletHash, proof.MHPathLow, nSeed, root->getHash());
            bool upPathValid = root->VerifyMerklePath(upTripletHash, proof.MHPathUp, nSeed, root->getHash());
            
            if (!lowPathValid || !upPathValid) {
                return false;
            }

            // Verify the nested difference proof
            return acc::Accumulator::verify_nested(
                const_cast<acc::AccValue&>(proof.tripletLow.getDigest()), 
                const_cast<acc::AccValue&>(proof.tripletUp.getDigest()), 
                const_cast<acc::NestedDifferenceProof&>(proof.NestDiffProof), 
                const_cast<acc::AccValue&>(proof.NestDiffRes), 
                publicKey
            );
        }
        
        default:
            return false;
    }
}


} // namespace gca