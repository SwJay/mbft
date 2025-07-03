/**
 * @file gca_tree.h
 * @brief Generic and Capacity-efficient Authenticated Aggregate tree (GCA²-tree) implementation
 */

#ifndef GCA_TREE_H
#define GCA_TREE_H

#include <vector>
#include <unordered_map>
#include <string>
#include <memory>
#include <set>

#include "acc/accumulator.h"
#include "mht/transaction.h"
#include "mht/merkle.h"
#include "gca/sorted_list.h"
#include "gca/prefix_set.h"
#include "gca/triplet.h"
#include "gca/query_defs.h"

namespace gca {

/**
 * @brief Main class for the GCA²-tree
 */
class GCATree {
public:
    GCATree(acc::AccPublicKey& pk, const std::vector<mht::Transaction>& transactions, unsigned int nSeed, uint32_t VMax, uint32_t MaxId, uint32_t unit);
    ~GCATree();

    std::pair<acc::Set, DimProof> querySingleInter(const Query& query);

    bool verifySingleInter(const Query& query, const DimProof& proof);

    uint32_t getRootHash() const;

    uint32_t getSizeADS();

private:
    // Helper methods
    std::vector<SortedList> createSortedLists(const std::vector<mht::Transaction>& transactions, const std::string& transactionType);
    std::vector<std::vector<PrefixSet>> generatePrefixSets(const std::vector<SortedList>& sortedLists, uint32_t VMax, uint32_t MaxId, uint32_t unit);
    std::vector<std::vector<acc::AccValue>> calculateDigests(const std::vector<std::vector<PrefixSet>>& prefixSets);
    std::vector<std::vector<Triplet>> constructTriplets(
        const std::vector<SortedList>& sortedLists,
        const std::vector<std::vector<PrefixSet>>& prefixSets,
        const std::vector<std::vector<acc::AccValue>>& digests);
    
    uint32_t hash(string t, uint32_t cid, uint32_t j, Triplet triplet, const unsigned int nSeed=0) const;

    std::tuple<int, int, DimBoundType> findBounds(const SortedList& sortedArray, uint32_t low, uint32_t high) const;

    // Member variables
    acc::AccPublicKey& publicKey;
    std::unordered_map<std::string, std::vector<SortedList>> sortedListsMap;
    std::unordered_map<std::string, std::vector<std::vector<PrefixSet>>> prefixSetsMap;
    std::unordered_map<std::string, std::vector<std::vector<acc::AccValue>>> digestsMap;
    std::unordered_map<std::string, std::vector<std::vector<Triplet>>> tripletsMap;
    mht::MerkleNode* root;
    unsigned int nSeed;

    uint32_t sizeADS = 0;
};

} // namespace gca

#endif // GCA_TREE_H
