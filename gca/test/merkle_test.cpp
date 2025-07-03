#include "gtest/gtest.h"
#include "mht/merkle.h"
#include <vector>
#include <string>
#include <cstdint>

// Helper function to create leaf hashes from strings
std::vector<uint32_t> create_leaf_hashes(const std::vector<std::string>& data, unsigned int seed) {
    std::vector<uint32_t> hashes;
    hashes.reserve(data.size());
    for (const auto& s : data) {
        std::vector<unsigned char> vec(s.begin(), s.end());
        uint32_t hash;
        MurmurHash3_x86_32(vec.data(), static_cast<int>(vec.size()), seed, &hash);
        hashes.push_back(hash);
    }
    return hashes;
}

// Helper function to get hashes from a vector of MerkleNode
std::vector<uint32_t> get_hashes_from_path(const std::vector<mht::MerklePathNode>& path) {
    std::vector<uint32_t> hashes;
    hashes.reserve(path.size());
    for (const auto& node : path) {
        hashes.push_back(node.hash);
    }
    return hashes;
}

// Define a test fixture for Merkle Tree tests
class MerkleTreeTest : public ::testing::Test {
protected:
    unsigned int seed = 42; // Consistent seed for hashing
};

// Test case for an even number of leaves
TEST_F(MerkleTreeTest, GetAndVerifyMerklePath_EvenLeaves) {
    std::vector<std::string> data = {"data1", "data2", "data3", "data4"};
    std::vector<uint32_t> leafHashes = create_leaf_hashes(data, seed);

    mht::MerkleNode merkleTree(leafHashes, seed);
    uint32_t rootHash = merkleTree.getHash();

    // --- Test Path for Leaf 0 ("data1") ---
    int leafIndex0 = 0;
    std::vector<mht::MerklePathNode> path0 = merkleTree.getMerklePath(leafIndex0);
    // Expected path for leaf 0 (hash(data1)): [hash(data2), hash(hash(data3)+hash(data4))]
    ASSERT_EQ(path0.size(), 2);
    // Calculate expected intermediate hashes
    uint32_t h1 = leafHashes[1]; // hash(data2)
    uint32_t h2 = leafHashes[2]; // hash(data3)
    uint32_t h3 = leafHashes[3]; // hash(data4)
    uint32_t h23 = merkleTree.Hash(h2, h3); // hash(hash(data3)+hash(data4))

    std::vector<uint32_t> expectedPathHashes0 = {h1, h23};
    std::vector<uint32_t> actualPathHashes0 = get_hashes_from_path(path0);
    EXPECT_EQ(actualPathHashes0, expectedPathHashes0) 
        << "Path hashes mismatch for leaf 0.\nExpected: " 
        << h1 << ", " << h23 << "\nGot: " 
        << (actualPathHashes0.empty() ? "empty" : std::to_string(actualPathHashes0[0]) + ", " + std::to_string(actualPathHashes0[1]));
    EXPECT_TRUE(merkleTree.VerifyMerklePath(leafHashes[leafIndex0], path0, seed, rootHash));
    merkleTree.getMerklePath(0);

    // Verification with wrong data
    EXPECT_FALSE(merkleTree.VerifyMerklePath(leafHashes[1], path0, seed, rootHash)); // Wrong leaf hash
    EXPECT_FALSE(merkleTree.VerifyMerklePath(leafHashes[leafIndex0], path0, seed, rootHash));      // Wrong index
    if (!path0.empty()) { // Modify path if not empty
        std::vector<mht::MerklePathNode> wrongPath0 = path0;
        wrongPath0[0] = mht::MerklePathNode(12345, true); // Tamper with path
        
        EXPECT_FALSE(merkleTree.VerifyMerklePath(leafHashes[leafIndex0], wrongPath0, seed, rootHash));
        std::cout << "----------------" << std::endl;
    }
    EXPECT_FALSE(merkleTree.VerifyMerklePath(leafHashes[leafIndex0], path0, seed, 12345)); // Wrong root hash
    std::cout << "----------------" << std::endl;


    // --- Test Path for Leaf 2 ("data3") ---
    std::cout << "========== Leaf 2: " << leafHashes[2] << std::endl;

    int leafIndex2 = 2;
    std::vector<mht::MerklePathNode> path2 = merkleTree.getMerklePath(leafIndex2);
    // Expected path for leaf 2 (hash(data3)): [hash(data4), hash(hash(data1)+hash(data2))]
    ASSERT_EQ(path2.size(), 2);
    uint32_t h0 = leafHashes[0]; // hash(data1)
    // h1 = leafHashes[1] // hash(data2) - already calculated
    uint32_t h01 = merkleTree.Hash(h0, h1); // hash(hash(data1)+hash(data2))
    // h3 = leafHashes[3] // hash(data4) - already calculated

    std::vector<uint32_t> expectedPathHashes2 = {h3, h01};
    EXPECT_EQ(get_hashes_from_path(path2), expectedPathHashes2);
    EXPECT_TRUE(merkleTree.VerifyMerklePath(leafHashes[leafIndex2], path2, seed, rootHash));
    EXPECT_FALSE(merkleTree.VerifyMerklePath(leafHashes[leafIndex2], path2, seed, 99999)); // Wrong root

}

// Test case for an odd number of leaves (last leaf duplicated)
TEST_F(MerkleTreeTest, GetAndVerifyMerklePath_OddLeaves) {
    std::vector<std::string> data = {"A", "B", "C"};
    std::vector<uint32_t> leafHashes = create_leaf_hashes(data, seed);

    mht::MerkleNode merkleTree(leafHashes, seed);
    uint32_t rootHash = merkleTree.getHash();

    // Calculate intermediate hashes (including duplication)
    uint32_t hA = leafHashes[0];
    uint32_t hB = leafHashes[1];
    uint32_t hC = leafHashes[2];
    uint32_t hAB = merkleTree.Hash(hA, hB);
    uint32_t hCC = merkleTree.Hash(hC, hC); // Due to duplication
    // Root should be Hash(hAB, hCC)

    // --- Test Path for Leaf 1 ("B") ---
    int leafIndex1 = 1;
    std::vector<mht::MerklePathNode> path1 = merkleTree.getMerklePath(leafIndex1);
    // Expected path for leaf 1 (hB): [hA, hCC]
    ASSERT_EQ(path1.size(), 2);
    std::vector<uint32_t> expectedPathHashes1 = {hA, hCC};
    std::vector<uint32_t> actualPathHashes1 = get_hashes_from_path(path1);
    EXPECT_EQ(actualPathHashes1, expectedPathHashes1)
        << "Path hashes mismatch for leaf 1.\nExpected: "
        << hA << ", " << hCC << "\nGot: "
        << (actualPathHashes1.empty() ? "empty" : std::to_string(actualPathHashes1[0]) + ", " + std::to_string(actualPathHashes1[1]));
    EXPECT_TRUE(merkleTree.VerifyMerklePath(leafHashes[leafIndex1], path1, seed, rootHash));
    EXPECT_FALSE(merkleTree.VerifyMerklePath(leafHashes[leafIndex1], path1, seed, hAB)); // Wrong root

    // --- Test Path for Leaf 2 ("C") ---
    int leafIndex2 = 2;
    std::vector<mht::MerklePathNode> path2 = merkleTree.getMerklePath(leafIndex2);
    // Expected path for leaf 2 (hC): [hC (the duplicated one), hAB]
    ASSERT_EQ(path2.size(), 2);
    std::vector<uint32_t> expectedPathHashes2 = {hC, hAB};
    EXPECT_EQ(get_hashes_from_path(path2), expectedPathHashes2);
    EXPECT_TRUE(merkleTree.VerifyMerklePath(leafHashes[leafIndex2], path2, seed, rootHash));
}

// Test case for a single leaf
TEST_F(MerkleTreeTest, GetAndVerifyMerklePath_SingleLeaf) {
    std::vector<std::string> data = {"lonely"};
    std::vector<uint32_t> leafHashes = create_leaf_hashes(data, seed);

    mht::MerkleNode merkleTree(leafHashes, seed);
    uint32_t rootHash = merkleTree.getHash();

    // Root hash should be the leaf hash itself
    EXPECT_EQ(rootHash, leafHashes[0]);

    int leafIndex0 = 0;
    std::vector<mht::MerklePathNode> path0 = merkleTree.getMerklePath(leafIndex0);

    // Path should be empty for a single leaf tree
    EXPECT_TRUE(path0.empty());

    // Verification should still work with an empty path
    EXPECT_TRUE(merkleTree.VerifyMerklePath(leafHashes[leafIndex0], path0, seed, rootHash));
    // Verification should fail if root hash is wrong
    EXPECT_FALSE(merkleTree.VerifyMerklePath(leafHashes[leafIndex0], path0, seed, rootHash + 1));
}

// Test case for empty input
TEST_F(MerkleTreeTest, GetAndVerifyMerklePath_EmptyTree) {
    std::vector<uint32_t> leafHashes; // Empty vector

    mht::MerkleNode merkleTree(leafHashes, seed);
    uint32_t rootHash = merkleTree.getHash(); // Should be 0 or default

    EXPECT_EQ(rootHash, 0); // Assuming default hash is 0 for empty tree

    // Getting path for any index should return empty path (or throw, depending on getLeaf)
    // Let's assume getMerklePath handles invalid index by returning empty path as implemented
    EXPECT_TRUE(merkleTree.getMerklePath(0).empty());
    EXPECT_TRUE(merkleTree.getMerklePath(-1).empty());

    // Verification doesn't make much sense for an empty tree, but test edge case
    // Trying to verify a leaf that doesn't exist against the default root (0)
    EXPECT_FALSE(merkleTree.VerifyMerklePath(123, {}, seed, rootHash));
}

// Test case for invalid leaf index
TEST_F(MerkleTreeTest, GetMerklePath_InvalidIndex) {
    std::vector<std::string> data = {"data1", "data2"};
    std::vector<uint32_t> leafHashes = create_leaf_hashes(data, seed);
    mht::MerkleNode merkleTree(leafHashes, seed);

    // Index out of bounds (positive)
    EXPECT_TRUE(merkleTree.getMerklePath(2).empty()); // Expect empty path based on current implementation
    // Index out of bounds (negative)
    EXPECT_TRUE(merkleTree.getMerklePath(-1).empty()); // Expect empty path based on current implementation
} 