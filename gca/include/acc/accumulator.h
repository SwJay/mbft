/**
 * @file accumulator.h
 * @brief Main class for the Expressive Set Accumulator
 */

#ifndef ACC_ACCUMULATOR_H
#define ACC_ACCUMULATOR_H

#include "acc/set.h"
#include "acc/keys.h"
#include "acc/acc_value.h"
#include "acc/proof.h"

namespace acc {

/**
 * @brief Main class for generating and verifying proofs
 */
class Accumulator {
public:
    // Key generation
    // static std::pair<AccSecretKey, AccPublicKey> genkey(uint64_t universe_size);
    static AccPublicKey genkey(std::string key_dir, std::string pbc_param_path, uint64_t universe_size);

    static std::pair<AccPublicKey, AccPublicKey> genkey_test(std::string key_dir, std::string pbc_param_path, uint64_t universe_size);
    
    // Digest setup
    static AccValue setup(const Set& set, AccPublicKey& pk);
    
    // Update accumulator
    // static AccValue update(const AccValue& acc, bool is_insert, uint64_t element, const AccPublicKey& pk);
    
    // Query operations
    static std::pair<Set, IntersectionProof> query_intersection(const Set& a, const Set& b, AccPublicKey& pk);
    static std::pair<Set, UnionProof> query_union(const Set& a, const Set& b, AccPublicKey& pk);
    static std::pair<Set, DifferenceProof> query_difference(const Set& a, const Set& b, AccPublicKey& pk);
    // Function queries
    static std::pair<uint64_t, CountProof> query_count(const Set& a, AccPublicKey& pk);
    static std::pair<uint64_t, SumProof> query_sum(const Set& a, AccPublicKey& pk);
    // Min and max queries
    static std::pair<uint64_t, MinProof> query_min(const Set& a, AccPublicKey& pk);
    static std::pair<uint64_t, MaxProof> query_max(const Set& a, AccPublicKey& pk);    
    // Range query
    static std::pair<Set, RangeProof> query_range(const Set& a, uint64_t l, uint64_t r, AccPublicKey& pk);

    // Nested query
    static std::pair<AccValue, NestedIntersectionProof> query_nested_intersection(const Set& a, const Set& b, AccPublicKey& pk);
    static std::pair<AccValue, NestedUnionProof> query_nested_union(const Set& a, const Set& b, AccPublicKey& pk);
    static std::pair<AccValue, NestedDifferenceProof> query_nested_difference(const Set& a, const Set& b, AccPublicKey& pk);
    static std::pair<AccValue, NestedRangeProof> query_nested_range(const Set& a, uint64_t l, uint64_t r, AccPublicKey& pk);


    // Verify proof
    static bool verify(AccValue& a_acc, AccValue& b_acc, Proof& proof, const Set& result, AccPublicKey& pk); 
    // Verify proof for numeric results
    static bool verify_aggr(AccValue& acc, AggrProof& proof, uint64_t result, AccPublicKey& pk);
    static bool verify_range(AccValue& acc, RangeProof& proof, const Set& result, uint64_t l, uint64_t r, AccPublicKey& pk);

    // Nested proof
    static bool verify_nested(AccValue& lhs_acc, AccValue& rhs_acc, NestedProof& proof, AccValue& res_acc, AccPublicKey& pk);
    static bool verify_nested_range(AccValue& acc, NestedRangeProof& proof, AccValue& res_acc, uint64_t l, uint64_t r, AccPublicKey& pk);
};

} // namespace acc

#endif // ACC_ACCUMULATOR_H 