#include <gtest/gtest.h>
#include "acc/accumulator.h"
#include <vector>
#include <random>

using namespace acc;

class UnionTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Generate public key
        std::string key_path = "../../data/key";
        std::string pbc_param_path = "";
        uint64_t universe_size = 128;
        pk = new AccPublicKey(Accumulator::genkey(key_path, pbc_param_path, universe_size));
    }

    void TearDown() override {
        delete pk;
    }

    AccPublicKey* pk;
};

TEST_F(UnionTest, BasicUnionTest) {
    // Create two sets
    Set set1, set2;
    std::vector<uint64_t> elements1 = {1, 2, 3, 4, 5};
    std::vector<uint64_t> elements2 = {3, 4, 5, 6, 7};
    
    for (auto elem : elements1) {
        set1.insert(elem);
    }
    for (auto elem : elements2) {
        set2.insert(elem);
    }

    // Create accumulators
    AccValue acc1 = Accumulator::setup(set1, *pk);
    AccValue acc2 = Accumulator::setup(set2, *pk);

    // Perform union operation and get proof
    auto [result, proof] = Accumulator::query_union(set1, set2, *pk);

    // Verify the proof
    EXPECT_TRUE(Accumulator::verify(acc1, acc2, proof, result, *pk));
}

TEST_F(UnionTest, EmptySetUnionTest) {
    // Create one empty set and one non-empty set
    Set set1, set2;
    std::vector<uint64_t> elements = {1, 2, 3, 4, 5};
    
    for (auto elem : elements) {
        set2.insert(elem);
    }

    // Create accumulators
    AccValue acc1 = Accumulator::setup(set1, *pk);
    AccValue acc2 = Accumulator::setup(set2, *pk);

    // Perform union operation and get proof
    auto [result, proof] = Accumulator::query_union(set1, set2, *pk);

    // Verify the proof
    EXPECT_TRUE(Accumulator::verify(acc1, acc2, proof, result, *pk));
}

TEST_F(UnionTest, LargeSetUnionTest) {
    // Create two large sets
    Set set1, set2;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<uint64_t> dis(1, 127);

    // Insert random elements into each set
    for (int i = 0; i < 75; ++i) {
        set1.insert(dis(gen));
        set2.insert(dis(gen));
    }

    // Create accumulators
    AccValue acc1 = Accumulator::setup(set1, *pk);
    AccValue acc2 = Accumulator::setup(set2, *pk);

    // Perform union operation and get proof
    auto [result, proof] = Accumulator::query_union(set1, set2, *pk);

    // Verify the proof
    EXPECT_TRUE(Accumulator::verify(acc1, acc2, proof, result, *pk));
}

TEST_F(UnionTest, InvalidProofTest) {
    // Create two sets
    Set set1, set2;
    std::vector<uint64_t> elements1 = {1, 2, 3};
    std::vector<uint64_t> elements2 = {4, 5, 6};
    
    for (auto elem : elements1) {
        set1.insert(elem);
    }
    for (auto elem : elements2) {
        set2.insert(elem);
    }

    // Create accumulators
    AccValue acc1 = Accumulator::setup(set1, *pk);
    AccValue acc2 = Accumulator::setup(set2, *pk);

    // Create a different result set
    Set result;
    result.insert(7);  // This element is not in either set

    // Generate union proof
    auto [_, proof] = Accumulator::query_union(set1, set2, *pk);

    // Verify the proof with incorrect result - should fail
    EXPECT_FALSE(Accumulator::verify(acc1, acc2, proof, result, *pk));
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
} 