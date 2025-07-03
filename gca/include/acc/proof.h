/**
 * @file proof.h
 * @brief Concrete proof implementations for the Expressive Set Accumulator
 */

#ifndef ACC_PROOF_H
#define ACC_PROOF_H

#include "acc/acc_value.h"
#include "acc/set.h"

namespace acc {

// Forward declarations
class AccPublicKey;

/**
 * @brief Base class for all proofs
 */
class Proof {
public:
    virtual ~Proof() = default;
    virtual bool verify(AccValue& lhs_acc, AccValue& rhs_acc, const Set& result, AccPublicKey& pk) = 0;
};

/**
 * @brief Base class for proofs with numeric results
 */
class AggrProof : public Proof {
public:
    virtual ~AggrProof() = default;
    
    // Override verify method to handle numeric results
    bool verify(AccValue& lhs_acc, AccValue& rhs_acc, const Set& result, AccPublicKey& pk) override {
        // This method will never be called directly
        return false;
    }
    
    // New verify method for numeric results
    virtual bool verify_aggr(AccValue& acc, uint64_t result, AccPublicKey& pk) = 0;
};

/**
 * @brief Proof for count operation
 */
class CountProof : public AggrProof {
public:
    element_t g_a_s;          // g^a(s)
    
public:
    // Constructor and destructor
    CountProof(pairing_t pairing);
    CountProof(const CountProof& other);                // 拷贝构造函数
    CountProof& operator=(const CountProof& other);    // 赋值运算符
    ~CountProof() noexcept override;

    size_t getSizeADS() {
        return element_length_in_bytes(g_a_s) * sizeof(char);
    }
    
    // Verify the proof with numeric result
    bool verify_aggr(AccValue& acc, uint64_t result, AccPublicKey& pk) override;
};

/**
 * @brief Proof for sum operation
 */
class SumProof : public AggrProof {
public:
    element_t g_b_s;          // pi_1: g^b(s)
    uint64_t count;           // pi_2: A(1)
    
public:
    // Constructor and destructor
    SumProof(pairing_t pairing);
    SumProof(const SumProof& other);                // 拷贝构造函数
    SumProof& operator=(const SumProof& other);    // 赋值运算符
    ~SumProof() noexcept override;

    size_t getSizeADS() {
        return element_length_in_bytes(g_b_s) * sizeof(char) + sizeof(count);
    }
    
    // Verify the proof with numeric result
    bool verify_aggr(AccValue& acc, uint64_t result, AccPublicKey& pk) override;
};

/**
 * @brief Proof for minimum operation
 */
class MinProof : public AggrProof {
public:
    element_t pi_min;
    
public:
    // Constructor and destructor
    MinProof(pairing_t pairing);
    MinProof(const MinProof& other);                // 拷贝构造函数
    MinProof& operator=(const MinProof& other);    // 赋值运算符
    ~MinProof() noexcept override;

    size_t getSizeADS() {
        return element_length_in_bytes(pi_min) * sizeof(char);
    }
    
    // Verify the proof with numeric result
    bool verify_aggr(AccValue& acc, uint64_t result, AccPublicKey& pk) override;
};

/**
 * @brief Proof for minimum operation
 */
class MaxProof : public AggrProof {
public:
    element_t pi_max;
    
public:
    // Constructor and destructor
    MaxProof(pairing_t pairing);
    MaxProof(const MaxProof& other);                // 拷贝构造函数
    MaxProof& operator=(const MaxProof& other);    // 赋值运算符
    ~MaxProof() noexcept override;

    size_t getSizeADS() {
        return element_length_in_bytes(pi_max) * sizeof(char);
    }
    
    // Verify the proof with numeric result
    bool verify_aggr(AccValue& acc, uint64_t result, AccPublicKey& pk) override;
};

/**
 * @brief Proof for intersection operation
 */
class IntersectionProof : public Proof {
public:
    element_t I_r;          // g^{\sum_{i∈I} r^i}
    element_t I_r_beta;     // g^{beta·\sum_{i∈I} r^i}
    element_t Q_s_r;        // g^{q(r,s)}
    element_t Q_s_r_delta;  // g^{delta·q(r,s)} (renamed from q_x_y_gamma)
    element_t L_r;          // g^{\sum_{i∈I} r^i/r}
    
public:
    // Constructor and destructor
    IntersectionProof(pairing_t pairing);
    IntersectionProof(const IntersectionProof& other);                // 拷贝构造函数
    IntersectionProof& operator=(const IntersectionProof& other);    // 赋值运算符
    ~IntersectionProof() noexcept override;
    
    // Verify the proof
    bool verify(AccValue& lhs_acc, AccValue& rhs_acc, const Set& result, AccPublicKey& pk) override;
};

/**
 * @brief Proof for Union operation
 */
class UnionProof : public Proof {
public:
    element_t I_r;          // g^{\sum_{i∈I} r^i}
    element_t I_r_beta;     // g^{beta·\sum_{i∈I} r^i}
    element_t Q_s_r;        // g^{q(r,s)}
    element_t Q_s_r_delta;  // g^{delta·q(r,s)} (renamed from q_x_y_gamma)
    element_t L_r;          // g^{\sum_{i∈I} r^i/r}
    
    element_t U_r;          // g^{\sum_{i∈U} r^i} (accumulator for difference set)
    
public:
    // Constructor and destructor
    UnionProof(pairing_t pairing);
    UnionProof(const UnionProof& other);                // 拷贝构造函数
    UnionProof& operator=(const UnionProof& other);    // 赋值运算符
    ~UnionProof() noexcept override;
    
    // Verify the proof
    bool verify(AccValue& lhs_acc, AccValue& rhs_acc, const Set& result, AccPublicKey& pk) override;
};

/**
 * @brief Proof for Difference operation
 */
class DifferenceProof : public Proof {
public:
    element_t I_r;          // g^{\sum_{i∈I} r^i}
    element_t I_r_beta;     // g^{beta·\sum_{i∈I} r^i}
    element_t Q_s_r;        // g^{q(r,s)}
    element_t Q_s_r_delta;  // g^{delta·q(r,s)} (renamed from q_x_y_gamma)
    element_t L_r;          // g^{\sum_{i∈I} r^i/r}
    element_t D_r;          // g^{\sum_{i∈D} r^i} (accumulator for difference set)
    
public:
    // Constructor and destructor
    DifferenceProof(pairing_t pairing);
    DifferenceProof(const DifferenceProof& other);                // 拷贝构造函数
    DifferenceProof& operator=(const DifferenceProof& other);    // 赋值运算符
    ~DifferenceProof() noexcept override;
    
    // Verify the proof
    bool verify(AccValue& lhs_acc, AccValue& rhs_acc, const Set& result, AccPublicKey& pk) override;
};

/**
 * @brief Proof for range query
 */
class RangeProof : public Proof {
public:
    element_t B_s;              // g^B(s)
    element_t B_s_alpha;        // g^{alpha·B(s)}
    element_t B_s_r;            // g^B(s,r)
    element_t B_s_r_gamma;      // g^{gamma·B(s,r)}
    element_t D_s;              // g^D(s)
    element_t D_s_alpha;        // g^{alpha·D(s)}
    element_t Z_s_r;            // g^Z(s,r)

    IntersectionProof pi_BC;    // pi_BC
    IntersectionProof pi_BD;    // pi_BD
    IntersectionProof pi_CD;    // pi_CD
    MaxProof pi_1;              // pi_1
    MinProof pi_2;              // pi_2
    
    uint64_t max_B;             // max_B
    uint64_t min_D;             // min_D

    AccValue d_B;
    AccValue d_D;
    
public:
    // Constructor and destructor
    RangeProof(pairing_t pairing);
    RangeProof(const RangeProof& other);                // 拷贝构造函数
    RangeProof& operator=(const RangeProof& other);    // 赋值运算符
    ~RangeProof() noexcept override;

    // Override verify method to handle numeric results
    bool verify(AccValue& lhs_acc, AccValue& rhs_acc, const Set& result, AccPublicKey& pk) override {
        // This method will never be called directly
        return false;
    }
    
    // New verify method for numeric results
    bool verify_range(AccValue& acc, const Set& result, uint64_t l, uint64_t r, AccPublicKey& pk);
};

/**
 * @brief Proof for nested difference operation
 */
class NestedProof : public Proof {
public:
    virtual ~NestedProof() = default;
    
    // Override verify method to handle numeric results
    bool verify(AccValue& lhs_acc, AccValue& rhs_acc, const Set& result, AccPublicKey& pk) override {
        // This method will never be called directly
        return false;
    }
    
    // New verify method for numeric results
    virtual bool verify_nested(AccValue& lhs_acc, AccValue& rhs_acc, AccValue& res_acc, AccPublicKey& pk) = 0;
};


/**
 * @brief Proof for nested intersection operation
 */
class NestedIntersectionProof : public NestedProof {
public:
    element_t I_r_beta;     // g^{beta·\sum_{i∈I} r^i}
    element_t Q_s_r;        // g^{q(r,s)}
    element_t Q_s_r_delta;  // g^{delta·q(r,s)} (renamed from q_x_y_gamma)
    element_t L_r;          // g^{\sum_{i∈I} r^i/r}

    element_t Z_s_r;
    element_t I_s_r_gamma;

    element_t I_s_alpha;     // g^{alpha·\sum_{i∈I} s^i}
    element_t Q_r_s;        // g^{q(r,s)}
    element_t Q_r_s_delta;  // g^{delta·q(r,s)} (renamed from q_x_y_gamma)
    element_t L_s;          // g^{\sum_{i∈I} s^i/s}
    
    element_t Z_r_s;
    element_t I_r_s_gamma;

public:
    // Constructor and destructor
    NestedIntersectionProof(pairing_t pairing);
    NestedIntersectionProof(const NestedIntersectionProof& other);                // 拷贝构造函数
    NestedIntersectionProof& operator=(const NestedIntersectionProof& other);    // 赋值运算符
    ~NestedIntersectionProof() noexcept override;
    

    // Verify the proof
    bool verify_nested(AccValue& lhs_acc, AccValue& rhs_acc, AccValue& res_acc, AccPublicKey& pk) override;
        
};


/**
 * @brief Proof for nested union operation
 */
class NestedUnionProof : public NestedProof {
public:
    element_t I_r;          // g^{\sum_{i∈I} r^i}
    element_t I_r_beta;     // g^{beta·\sum_{i∈I} r^i}
    element_t Q_s_r;        // g^{q(r,s)}
    element_t Q_s_r_delta;  // g^{delta·q(r,s)} (renamed from q_x_y_gamma)
    element_t L_r;          // g^{\sum_{i∈I} r^i/r}

    element_t Z_s_r;
    element_t U_s_r_gamma;
    element_t Z_r_s;
    element_t U_r_s_gamma;

public:
    // Constructor and destructor
    NestedUnionProof(pairing_t pairing);
    NestedUnionProof(const NestedUnionProof& other);                // 拷贝构造函数
    NestedUnionProof& operator=(const NestedUnionProof& other);    // 赋值运算符
    ~NestedUnionProof() noexcept override;

    size_t getSizeADS();
    
    // Verify the proof
    bool verify_nested(AccValue& lhs_acc, AccValue& rhs_acc, AccValue& res_acc, AccPublicKey& pk) override;
};

/**
 * @brief Proof for Difference operation
 */
class NestedDifferenceProof : public NestedProof {
public:
    element_t I_r;          // g^{\sum_{i∈I} r^i}
    element_t I_r_beta;     // g^{beta·\sum_{i∈I} r^i}
    element_t Q_s_r;        // g^{q(r,s)}
    element_t Q_s_r_delta;  // g^{delta·q(r,s)} (renamed from q_x_y_gamma)
    element_t L_r;          // g^{\sum_{i∈I} r^i/r}

    element_t Z_s_r;
    element_t D_s_r_gamma;
    element_t Z_r_s;
    element_t D_r_s_gamma;

public:
    // Constructor and destructor
    NestedDifferenceProof(pairing_t pairing);
    NestedDifferenceProof(const NestedDifferenceProof& other);                // 拷贝构造函数
    NestedDifferenceProof& operator=(const NestedDifferenceProof& other);    // 赋值运算符
    ~NestedDifferenceProof() noexcept override;

    size_t getSizeADS();
    
    // Verify the proof
    bool verify_nested(AccValue& lhs_acc, AccValue& rhs_acc, AccValue& res_acc, AccPublicKey& pk) override;
};

/**
 * @brief Proof for nested range query
 */
class NestedRangeProof : public Proof {
public:
    element_t B_s;              // g^B(s)
    element_t B_s_alpha;        // g^{alpha·B(s)}
    element_t B_s_r;            // g^B(s,r)
    element_t B_s_r_gamma;      // g^{gamma·B(s,r)}
    element_t D_s;              // g^D(s)
    element_t D_s_alpha;        // g^{alpha·D(s)}
    element_t Z_s_r;            // g^Z(s,r)

    IntersectionProof pi_BC;    // pi_BC
    IntersectionProof pi_BD;    // pi_BD
    IntersectionProof pi_CD;    // pi_CD
    MaxProof pi_1;              // pi_1
    MinProof pi_2;              // pi_2
    
    uint64_t max_B;             // max_B
    uint64_t min_D;             // min_D

    AccValue d_B;
    AccValue d_D;

    element_t ZC_s_r;
    element_t ZC_r_s;

    uint64_t max_C;
    MaxProof pi_C_max;
    uint64_t min_C;
    MinProof pi_C_min;

public:
    // Constructor and destructor
    NestedRangeProof(pairing_t pairing);
    NestedRangeProof(const NestedRangeProof& other);                // 拷贝构造函数
    NestedRangeProof& operator=(const NestedRangeProof& other);    // 赋值运算符
    ~NestedRangeProof() noexcept override;

    size_t getSizeADS();

    // Override verify method to handle numeric results
    bool verify(AccValue& lhs_acc, AccValue& rhs_acc, const Set& result, AccPublicKey& pk) override {
        // This method will never be called directly
        return false;
    }
    
    // New verify method for numeric results
    bool verify_nested_range(AccValue& acc, AccValue& res_acc, uint64_t l, uint64_t r, AccPublicKey& pk);
};

} // namespace acc

#endif // ACC_PROOF_H 