#ifndef GCA_QUERY_DEFS_H
#define GCA_QUERY_DEFS_H

#include <cstdint>
#include <string>
#include <vector>
#include <set> // For address filters or similar
#include <memory>
#include <variant>
#include <queue>

#include "acc/proof.h"
#include "gca/triplet.h"

namespace gca {

// Enum for supported aggregate functions
enum class AggType {
    COUNT,
    SUM,
    MAX,
    MIN,
    AVG // Average might require SUM and COUNT internally
};

// Basic structure to represent a query
// Based on usage in gca/src/main.cpp and GCA² paper concepts
struct Query {
    string TxnType;
    AggType agg;        // The aggregate function (COUNT, SUM, etc.)
    uint32_t cid;         // Index or identifier of the attribute to aggregate on
    uint32_t lowerBound;   // Lower bound for the aggregation attribute range (inclusive)
    uint32_t upperBound;   // Upper bound for the aggregation attribute range (inclusive)

    // Time window
    uint32_t beginBlock;         // Start timestamp (inclusive)
    uint32_t endBlock;  // End timestamp (inclusive)

    // Constructor based on main.cpp example (adjust as needed)
    Query(string TxnType, AggType agg, uint32_t attrIdx, uint32_t lower, uint32_t upper,
          uint32_t tStart, uint32_t tEnd):
          TxnType(TxnType),
          agg(agg),
          cid(attrIdx),
          lowerBound(lower),
          upperBound(upper),
          beginBlock(tStart),
          endBlock(tEnd) {}

    // Default constructor
    Query() = default;
};

enum class DimBoundType {
    MAX_BELOW_LOW,      // Q.low > O.max
    MIN_ABOVE_UP,       // Q.up < O.min
    MIN_INCLUDED,       // Q.low <= O.min
    IN_RANGE            // Normal case
};

struct DimProof {
    DimBoundType boundType;

    int lowIndex;
    int highIndex;
    Triplet tripletLow;
    Triplet tripletUp;
    std::vector<mht::MerklePathNode> MHPathLow;
    std::vector<mht::MerklePathNode> MHPathUp;
    
    acc::NestedDifferenceProof NestDiffProof;
    acc::AccValue NestDiffRes;

    DimProof(pairing_t pairing):
        boundType(DimBoundType::IN_RANGE),
        lowIndex(0),
        highIndex(0),
        tripletLow(0, 0, acc::AccValue(pairing)),
        tripletUp(0, 0, acc::AccValue(pairing)),
        MHPathLow({}),
        MHPathUp({}),
        NestDiffProof(pairing),
        NestDiffRes(pairing) {};

    size_t getSizeADS() {
        // size_t size = sizeof(lowIndex) + sizeof(highIndex) + tripletLow.getSizeADS() + tripletUp.getSizeADS() + NestDiffProof.getSizeADS() + NestDiffRes.getSizeADS();
        // for (const mht::MerklePathNode& node : MHPathLow) {
        //     size += node.getSizeADS();
        // }
        // for (const mht::MerklePathNode& node : MHPathUp) {
        //     size += node.getSizeADS();
        // }
        return NestDiffProof.getSizeADS() + NestDiffRes.getSizeADS();
    }

};

struct VO {
    std::queue<DimProof> dimProofs;
    std::queue<acc::AccValue> accDimRes;

    acc::NestedRangeProof nestedRangeProof;
    acc::AccValue nestedRangeRes;

    std::queue<acc::NestedUnionProof> nestedUnionProofs;
    std::queue<acc::AccValue> nestedUnionRes;

    acc::AccValue unionRes;

    // 使用variant存储不同类型的AggrProof
    std::variant<acc::CountProof, acc::SumProof, acc::MinProof, acc::MaxProof> aggrProof;

    VO(pairing_t pairing):
        nestedRangeProof(pairing),
        nestedRangeRes(pairing),
        unionRes(pairing),
        // 初始化为CountProof，之后可以根据需要替换
        aggrProof(acc::CountProof(pairing)) {};
    
    size_t getSizeADS() {
        size_t size = 0;

        // Iterate over a copy of dimProofs
        if (!dimProofs.empty()) {
            std::queue<DimProof> dimProofsCopy = dimProofs;
            while (!dimProofsCopy.empty()) {
                size += dimProofsCopy.front().getSizeADS();
                dimProofsCopy.pop();
            }
        }

        // Iterate over a copy of accDimRes
        if (!accDimRes.empty()) {
            std::queue<acc::AccValue> accDimResCopy = accDimRes;
            while (!accDimResCopy.empty()) {
                size += accDimResCopy.front().getSizeADS();
                accDimResCopy.pop();
            }
        }

        std::cout << "dim size: " << size << std::endl;
        
        size += nestedRangeProof.getSizeADS();
        size += nestedRangeRes.getSizeADS();

        // // Iterate over a copy of nestedUnionProofs
        // if (!nestedUnionProofs.empty()) {
        //     std::queue<acc::NestedUnionProof> nestedUnionProofsCopy = nestedUnionProofs;
        //     while (!nestedUnionProofsCopy.empty()) {
        //         size += nestedUnionProofsCopy.front().getSizeADS();
        //         nestedUnionProofsCopy.pop();
        //     }
        // }

        // // Iterate over a copy of nestedUnionRes
        // if (!nestedUnionRes.empty()) {
        //     std::queue<acc::AccValue> nestedUnionResCopy = nestedUnionRes;
        //     while (!nestedUnionResCopy.empty()) {
        //         size += nestedUnionResCopy.front().getSizeADS();
        //         nestedUnionResCopy.pop();
        //     }
        // }

        size += unionRes.getSizeADS();

        if (std::holds_alternative<acc::CountProof>(aggrProof)) {
            size += std::get<acc::CountProof>(aggrProof).getSizeADS();
        } else if (std::holds_alternative<acc::SumProof>(aggrProof)) {
            size += std::get<acc::SumProof>(aggrProof).getSizeADS();
        } else if (std::holds_alternative<acc::MinProof>(aggrProof)) {
            size += std::get<acc::MinProof>(aggrProof).getSizeADS();
        } else if (std::holds_alternative<acc::MaxProof>(aggrProof)) {
            size += std::get<acc::MaxProof>(aggrProof).getSizeADS();
        }

        return size;
    }
};

} // namespace gca

#endif // GCA_QUERY_DEFS_H
