/**
 * @file proof.cpp
 * @brief Implementation of proof classes for the Expressive Set Accumulator
 */

#include <iostream>

#include "acc/proof.h"

namespace acc {

/*
* IntersectionProof ====================================================
*/

IntersectionProof::IntersectionProof(pairing_t pairing) {
    // Initialize elements
    element_init_G1(I_r, pairing);
    element_init_G1(I_r_beta, pairing);
    element_init_G1(Q_s_r, pairing);
    element_init_G1(Q_s_r_delta, pairing);
    element_init_G1(L_r, pairing);
    
    // Set to identity elements
    element_set1(I_r);
    element_set1(I_r_beta);
    element_set1(Q_s_r);
    element_set1(Q_s_r_delta);
    element_set1(L_r);
}

// Copy constructor
IntersectionProof::IntersectionProof(const IntersectionProof& other) {
    // Initialize elements in the same field
    element_init_same_as(I_r, const_cast<element_t&>(other.I_r));
    element_init_same_as(I_r_beta, const_cast<element_t&>(other.I_r_beta));
    element_init_same_as(Q_s_r, const_cast<element_t&>(other.Q_s_r));
    element_init_same_as(Q_s_r_delta, const_cast<element_t&>(other.Q_s_r_delta));
    element_init_same_as(L_r, const_cast<element_t&>(other.L_r));
    
    // Copy values
    element_set(I_r, const_cast<element_t&>(other.I_r));
    element_set(I_r_beta, const_cast<element_t&>(other.I_r_beta));
    element_set(Q_s_r, const_cast<element_t&>(other.Q_s_r));
    element_set(Q_s_r_delta, const_cast<element_t&>(other.Q_s_r_delta));
    element_set(L_r, const_cast<element_t&>(other.L_r));
}

// Assignment operator
IntersectionProof& IntersectionProof::operator=(const IntersectionProof& other) {
    if (this != &other) {
        // Clear existing elements
        element_clear(I_r);
        element_clear(I_r_beta);
        element_clear(Q_s_r);
        element_clear(Q_s_r_delta);
        element_clear(L_r);
        
        // Initialize elements in the same field
        element_init_same_as(I_r, const_cast<element_t&>(other.I_r));
        element_init_same_as(I_r_beta, const_cast<element_t&>(other.I_r_beta));
        element_init_same_as(Q_s_r, const_cast<element_t&>(other.Q_s_r));
        element_init_same_as(Q_s_r_delta, const_cast<element_t&>(other.Q_s_r_delta));
        element_init_same_as(L_r, const_cast<element_t&>(other.L_r));
        
        // Copy values
        element_set(I_r, const_cast<element_t&>(other.I_r));
        element_set(I_r_beta, const_cast<element_t&>(other.I_r_beta));
        element_set(Q_s_r, const_cast<element_t&>(other.Q_s_r));
        element_set(Q_s_r_delta, const_cast<element_t&>(other.Q_s_r_delta));
        element_set(L_r, const_cast<element_t&>(other.L_r));
    }
    return *this;
}

IntersectionProof::~IntersectionProof() noexcept {
    // Clear elements
    element_clear(I_r);
    element_clear(I_r_beta);
    element_clear(Q_s_r);
    element_clear(Q_s_r_delta);
    element_clear(L_r);
}

bool IntersectionProof::verify(AccValue& lhs_acc, AccValue& rhs_acc, 
                              const Set& result, AccPublicKey& pk) {    
    /*
    * Check 1: e(A_s, B_{r,s}) = e(I_r, g^{s^q}) × e(Q_{s,r}, g)
    */
    element_t left_1, right_1, right_1_1, right_1_2;
    element_init_GT(left_1, pk.pairing);
    element_init_GT(right_1, pk.pairing);
    element_init_GT(right_1_1, pk.pairing);
    element_init_GT(right_1_2, pk.pairing);
    
    // e(A_s, B_{r,s})
    element_pairing(left_1, lhs_acc.g_s, rhs_acc.g_r_s);
    // e(I_r, g^{s^q})
    element_pairing(right_1_1, I_r, pk.g_sq);
    // e(Q_{s,r}, g)
    element_pairing(right_1_2, Q_s_r, pk.g);
    // Combine: right_1 = right_1_1 × right_1_2
    element_mul(right_1, right_1_1, right_1_2);
    // Check if left_1 == right_1
    bool check_1 = (element_cmp(left_1, right_1) == 0);
    if (!check_1) {
        std::cout << "intersection check_1: " << check_1 << std::endl;
    }


    // Clean up
    element_clear(left_1);
    element_clear(right_1);
    element_clear(right_1_1);
    element_clear(right_1_2);
    
    /*
    * Check 2: e(I_r, g^alpha) = e(I_{r,alpha}, g)
    */
    element_t left_2, right_2;
    element_init_GT(left_2, pk.pairing);
    element_init_GT(right_2, pk.pairing);
    
    // e(I_r, g^alpha)
    element_pairing(left_2, I_r, pk.g_beta);
    // e(I_{r,alpha}, g)
    element_pairing(right_2, I_r_beta, pk.g);
    // Check if left_2 == right_2
    bool check_2 = (element_cmp(left_2, right_2) == 0);
    if (!check_2) {
        std::cout << "intersection check_2: " << check_2 << std::endl;
    }

    // Clean up
    element_clear(left_2);
    element_clear(right_2);
    
    /*
    * Check 3: e(Q_{s,r}, g^delta) = e(Q_{s,r,delta}, g)
    */
    element_t left_3, right_3;
    element_init_GT(left_3, pk.pairing);
    element_init_GT(right_3, pk.pairing);
    
    // e(Q_{s,r}, g^delta)
    element_pairing(left_3, Q_s_r, pk.g_delta);
    // e(Q_{s,r,delta}, g)
    element_pairing(right_3, Q_s_r_delta, pk.g);
    // Check if left_3 == right_3
    bool check_3 = (element_cmp(left_3, right_3) == 0);
    if (!check_3) {
        std::cout << "intersection check_3: " << check_3 << std::endl;
    }

    // Clean up
    element_clear(left_3);
    element_clear(right_3);
    
    /*
    * Check 4: e(I_r, g) = e(L_r, g^r)
    */
    element_t left_4, right_4;
    element_init_GT(left_4, pk.pairing);
    element_init_GT(right_4, pk.pairing);
    
    // e(I_r, g)
    element_pairing(left_4, I_r, pk.g);
    // e(L_r, g^r)
    element_pairing(right_4, L_r, pk.g_ri[0]);
    // Check if left_4 == right_4
    bool check_4 = (element_cmp(left_4, right_4) == 0);
    if (!check_4) {
        std::cout << "intersection check_4: " << check_4 << std::endl;
    }

    // Clean up
    element_clear(left_4);
    element_clear(right_4);
    
    /*
    * Check 5: Verify I_r = g^{\sum_{i∈I} r^i} (I_r is indeed the accumulator for the intersection set)
    */
    element_t computed_I_r;
    element_init_G1(computed_I_r, pk.pairing);
    element_set1(computed_I_r);
    
    // Compute g^{\sum_{i∈I} r^i}
    for (const auto& i : result) {
        if (i < 1 || i >= pk.q) continue;  // Skip invalid elements
        element_mul(computed_I_r, computed_I_r, pk.g_ri[i - 1]);  // Multiply by g^{r^i}
    }
    
    // Check if computed_I_r == I_r
    bool check_5 = (element_cmp(computed_I_r, I_r) == 0);
    if (!check_5) {
        std::cout << "intersection check_5: " << check_5 << std::endl;
    }
    
    // Clean up
    element_clear(computed_I_r);
    
    // All checks must pass
    return check_1 && check_2 && check_3 && check_4 && check_5;
}

/*
* UnionProof ====================================================
*/

UnionProof::UnionProof(pairing_t pairing) {
    // Initialize elements
    element_init_G1(I_r, pairing);
    element_init_G1(I_r_beta, pairing);
    element_init_G1(Q_s_r, pairing);
    element_init_G1(Q_s_r_delta, pairing);
    element_init_G1(L_r, pairing);
    element_init_G1(U_r, pairing);

    // Set to identity elements
    element_set1(I_r);
    element_set1(I_r_beta);
    element_set1(Q_s_r);
    element_set1(Q_s_r_delta);
    element_set1(L_r);
    element_set1(U_r);
}

// Copy constructor
UnionProof::UnionProof(const UnionProof& other) {
    // Initialize elements in the same field
    element_init_same_as(I_r, const_cast<element_t&>(other.I_r));
    element_init_same_as(I_r_beta, const_cast<element_t&>(other.I_r_beta));
    element_init_same_as(Q_s_r, const_cast<element_t&>(other.Q_s_r));
    element_init_same_as(Q_s_r_delta, const_cast<element_t&>(other.Q_s_r_delta));
    element_init_same_as(L_r, const_cast<element_t&>(other.L_r));
    element_init_same_as(U_r, const_cast<element_t&>(other.U_r));
    
    // Copy values
    element_set(I_r, const_cast<element_t&>(other.I_r));
    element_set(I_r_beta, const_cast<element_t&>(other.I_r_beta));
    element_set(Q_s_r, const_cast<element_t&>(other.Q_s_r));
    element_set(Q_s_r_delta, const_cast<element_t&>(other.Q_s_r_delta));
    element_set(L_r, const_cast<element_t&>(other.L_r));
    element_set(U_r, const_cast<element_t&>(other.U_r));
}

// Assignment operator
UnionProof& UnionProof::operator=(const UnionProof& other) {
    if (this != &other) {
        // Clear existing elements
        element_clear(I_r);
        element_clear(I_r_beta);
        element_clear(Q_s_r);
        element_clear(Q_s_r_delta);
        element_clear(L_r);
        element_clear(U_r);

        // Initialize elements in the same field
        element_init_same_as(I_r, const_cast<element_t&>(other.I_r));
        element_init_same_as(I_r_beta, const_cast<element_t&>(other.I_r_beta));
        element_init_same_as(Q_s_r, const_cast<element_t&>(other.Q_s_r));
        element_init_same_as(Q_s_r_delta, const_cast<element_t&>(other.Q_s_r_delta));
        element_init_same_as(L_r, const_cast<element_t&>(other.L_r));
        element_init_same_as(U_r, const_cast<element_t&>(other.U_r));
        
        // Copy values
        element_set(I_r, const_cast<element_t&>(other.I_r));
        element_set(I_r_beta, const_cast<element_t&>(other.I_r_beta));
        element_set(Q_s_r, const_cast<element_t&>(other.Q_s_r));
        element_set(Q_s_r_delta, const_cast<element_t&>(other.Q_s_r_delta));
        element_set(L_r, const_cast<element_t&>(other.L_r));
        element_set(U_r, const_cast<element_t&>(other.U_r));
    }
    return *this;
}

UnionProof::~UnionProof() noexcept {
    // Clear elements
    element_clear(I_r);
    element_clear(I_r_beta);
    element_clear(Q_s_r);
    element_clear(Q_s_r_delta);
    element_clear(L_r);
    element_clear(U_r);
}

bool UnionProof::verify(AccValue& lhs_acc, AccValue& rhs_acc, 
                             const Set& result, AccPublicKey& pk) {
    // First verify the intersection part of the proof
    /*
    * Check 1: e(A_s, B_{r,s}) = e(I_r, g^{s^q}) × e(Q_{s,r}, g)
    */
    element_t left_1, right_1, right_1_1, right_1_2;
    element_init_GT(left_1, pk.pairing);
    element_init_GT(right_1, pk.pairing);
    element_init_GT(right_1_1, pk.pairing);
    element_init_GT(right_1_2, pk.pairing);
    
    // e(A_s, B_{r,s})
    element_pairing(left_1, lhs_acc.g_s, rhs_acc.g_r_s);
    // e(I_r, g^{s^q})
    element_pairing(right_1_1, I_r, pk.g_sq);
    // e(Q_{s,r}, g)
    element_pairing(right_1_2, Q_s_r, pk.g);
    // Combine: right_1 = right_1_1 × right_1_2
    element_mul(right_1, right_1_1, right_1_2);
    // Check if left_1 == right_1
    bool check_1 = (element_cmp(left_1, right_1) == 0);
    if (!check_1) {
        std::cout << "union check_1: " << check_1 << std::endl;
    }

    // Clean up
    element_clear(left_1);
    element_clear(right_1);
    element_clear(right_1_1);
    element_clear(right_1_2);
    
    /*
    * Check 2: e(I_r, g^beta) = e(I_{r,beta}, g)
    */
    element_t left_2, right_2;
    element_init_GT(left_2, pk.pairing);
    element_init_GT(right_2, pk.pairing);
    
    // e(I_r, g^beta)
    element_pairing(left_2, I_r, pk.g_beta);
    // e(I_{r,beta}, g)
    element_pairing(right_2, I_r_beta, pk.g);
    // Check if left_2 == right_2
    bool check_2 = (element_cmp(left_2, right_2) == 0);
    if (!check_2) {
        std::cout << "union check_2: " << check_2 << std::endl;
    }

    // Clean up
    element_clear(left_2);
    element_clear(right_2);
    
    /*
    * Check 3: e(Q_{s,r}, g^delta) = e(Q_{s,r,delta}, g)
    */
    element_t left_3, right_3;
    element_init_GT(left_3, pk.pairing);
    element_init_GT(right_3, pk.pairing);
    
    // e(Q_{s,r}, g^delta)
    element_pairing(left_3, Q_s_r, pk.g_delta);
    // e(Q_{s,r,delta}, g)
    element_pairing(right_3, Q_s_r_delta, pk.g);
    // Check if left_3 == right_3
    bool check_3 = (element_cmp(left_3, right_3) == 0);
    if (!check_3) {
        std::cout << "union check_3: " << check_3 << std::endl;
    }

    // Clean up
    element_clear(left_3);
    element_clear(right_3);
    
    /*
    * Check 4: e(I_r, g) = e(L_r, g^r)
    */
    element_t left_4, right_4;
    element_init_GT(left_4, pk.pairing);
    element_init_GT(right_4, pk.pairing);
    
    // e(I_r, g)
    element_pairing(left_4, I_r, pk.g);
    // e(L_r, g^r)
    element_pairing(right_4, L_r, pk.g_ri[0]);
    // Check if left_4 == right_4
    bool check_4 = (element_cmp(left_4, right_4) == 0);
    if (!check_4) {
        std::cout << "union check_4: " << check_4 << std::endl;
    }

    // Clean up
    element_clear(left_4);
    element_clear(right_4);
    
    /*
    * Check 6: Verify U_r = A_r * B_r / I_r (accumulator for union = accumulator for A / accumulator for intersection)
    */
    element_t expected_U_r, I_r_inv;
    element_init_G1(expected_U_r, pk.pairing);
    element_init_G1(I_r_inv, pk.pairing);
    
    // Compute the inverse of I_r
    element_set(I_r_inv, I_r);
    element_invert(I_r_inv, I_r_inv);
    
    // expected_U_r = A_r * B_r / I_r = A_r * B_r * I_r^(-1)
    element_mul(expected_U_r, lhs_acc.g_r, rhs_acc.g_r);
    element_mul(expected_U_r, expected_U_r, I_r_inv);
    
    // Check if expected_U_r == U_r
    bool check_6 = (element_cmp(expected_U_r, U_r) == 0);
    if (!check_6) {
        std::cout << "union check_6: " << check_6 << std::endl;
    }
    
    // Clean up
    element_clear(expected_U_r);
    element_clear(I_r_inv);
    
    /*
    * Check 7: Verify U_r = g^{\sum_{i∈U} r^i} (U_r is indeed the accumulator for the union set)
    */
    element_t computed_U_r;
    element_init_G1(computed_U_r, pk.pairing);
    element_set1(computed_U_r);
    
    // Compute g^{\sum_{i∈U} r^i}
    for (const auto& i : result) {
        if (i < 1 || i >= pk.q) continue;  // Skip invalid elements
        element_mul(computed_U_r, computed_U_r, pk.g_ri[i - 1]);  // Multiply by g^{r^i}
    }
    
    // Check if computed_U_r == U_r
    bool check_7 = (element_cmp(computed_U_r, U_r) == 0);
    if (!check_7) {
        std::cout << "union check_7: " << check_7 << std::endl;
    }
    
    // Clean up
    element_clear(computed_U_r);
    
    // All checks must pass
    return check_1 && check_2 && check_3 && check_4 && check_6 && check_7;
}

/*
* DifferenceProof ====================================================
*/

DifferenceProof::DifferenceProof(pairing_t pairing) {
    // Initialize elements
    element_init_G1(I_r, pairing);
    element_init_G1(I_r_beta, pairing);
    element_init_G1(Q_s_r, pairing);
    element_init_G1(Q_s_r_delta, pairing);
    element_init_G1(L_r, pairing);
    element_init_G1(D_r, pairing);
    
    // Set to identity elements
    element_set1(I_r);
    element_set1(I_r_beta);
    element_set1(Q_s_r);
    element_set1(Q_s_r_delta);
    element_set1(L_r);
    element_set1(D_r);
}

// Copy constructor
DifferenceProof::DifferenceProof(const DifferenceProof& other) {
    // Initialize elements in the same field
    element_init_same_as(I_r, const_cast<element_t&>(other.I_r));
    element_init_same_as(I_r_beta, const_cast<element_t&>(other.I_r_beta));
    element_init_same_as(Q_s_r, const_cast<element_t&>(other.Q_s_r));
    element_init_same_as(Q_s_r_delta, const_cast<element_t&>(other.Q_s_r_delta));
    element_init_same_as(L_r, const_cast<element_t&>(other.L_r));
    element_init_same_as(D_r, const_cast<element_t&>(other.D_r));
    
    // Copy values
    element_set(I_r, const_cast<element_t&>(other.I_r));
    element_set(I_r_beta, const_cast<element_t&>(other.I_r_beta));
    element_set(Q_s_r, const_cast<element_t&>(other.Q_s_r));
    element_set(Q_s_r_delta, const_cast<element_t&>(other.Q_s_r_delta));
    element_set(L_r, const_cast<element_t&>(other.L_r));
    element_set(D_r, const_cast<element_t&>(other.D_r));
}

// Assignment operator
DifferenceProof& DifferenceProof::operator=(const DifferenceProof& other) {
    if (this != &other) {
        // Clear existing elements
        element_clear(I_r);
        element_clear(I_r_beta);
        element_clear(Q_s_r);
        element_clear(Q_s_r_delta);
        element_clear(L_r);
        element_clear(D_r);
        
        // Initialize elements in the same field
        element_init_same_as(I_r, const_cast<element_t&>(other.I_r));
        element_init_same_as(I_r_beta, const_cast<element_t&>(other.I_r_beta));
        element_init_same_as(Q_s_r, const_cast<element_t&>(other.Q_s_r));
        element_init_same_as(Q_s_r_delta, const_cast<element_t&>(other.Q_s_r_delta));
        element_init_same_as(L_r, const_cast<element_t&>(other.L_r));
        element_init_same_as(D_r, const_cast<element_t&>(other.D_r));
        
        // Copy values
        element_set(I_r, const_cast<element_t&>(other.I_r));
        element_set(I_r_beta, const_cast<element_t&>(other.I_r_beta));
        element_set(Q_s_r, const_cast<element_t&>(other.Q_s_r));
        element_set(Q_s_r_delta, const_cast<element_t&>(other.Q_s_r_delta));
        element_set(L_r, const_cast<element_t&>(other.L_r));
        element_set(D_r, const_cast<element_t&>(other.D_r));
    }
    return *this;
}

DifferenceProof::~DifferenceProof() noexcept {
    // Clear elements
    element_clear(I_r);
    element_clear(I_r_beta);
    element_clear(Q_s_r);
    element_clear(Q_s_r_delta);
    element_clear(L_r);
    element_clear(D_r);
}

bool DifferenceProof::verify(AccValue& lhs_acc, AccValue& rhs_acc, 
                             const Set& result, AccPublicKey& pk) {
    // First verify the intersection part of the proof
    /*
    * Check 1: e(A_s, B_{r,s}) = e(I_r, g^{s^q}) × e(Q_{s,r}, g)
    */
    element_t left_1, right_1, right_1_1, right_1_2;
    element_init_GT(left_1, pk.pairing);
    element_init_GT(right_1, pk.pairing);
    element_init_GT(right_1_1, pk.pairing);
    element_init_GT(right_1_2, pk.pairing);
    
    // e(A_s, B_{r,s})
    element_pairing(left_1, lhs_acc.g_s, rhs_acc.g_r_s);
    // e(I_r, g^{s^q})
    element_pairing(right_1_1, I_r, pk.g_sq);
    // e(Q_{s,r}, g)
    element_pairing(right_1_2, Q_s_r, pk.g);
    // Combine: right_1 = right_1_1 × right_1_2
    element_mul(right_1, right_1_1, right_1_2);
    // Check if left_1 == right_1
    bool check_1 = (element_cmp(left_1, right_1) == 0);
    if (!check_1) {
        std::cout << "difference check_1: " << check_1 << std::endl;
    }

    // Clean up
    element_clear(left_1);
    element_clear(right_1);
    element_clear(right_1_1);
    element_clear(right_1_2);
    
    /*
    * Check 2: e(I_r, g^beta) = e(I_{r,beta}, g)
    */
    element_t left_2, right_2;
    element_init_GT(left_2, pk.pairing);
    element_init_GT(right_2, pk.pairing);
    
    // e(I_r, g^beta)
    element_pairing(left_2, I_r, pk.g_beta);
    // e(I_{r,beta}, g)
    element_pairing(right_2, I_r_beta, pk.g);
    // Check if left_2 == right_2
    bool check_2 = (element_cmp(left_2, right_2) == 0);
    if (!check_2) {
        std::cout << "difference check_2: " << check_2 << std::endl;
    }

    // Clean up
    element_clear(left_2);
    element_clear(right_2);
    
    /*
    * Check 3: e(Q_{s,r}, g^delta) = e(Q_{s,r,delta}, g)
    */
    element_t left_3, right_3;
    element_init_GT(left_3, pk.pairing);
    element_init_GT(right_3, pk.pairing);
    
    // e(Q_{s,r}, g^delta)
    element_pairing(left_3, Q_s_r, pk.g_delta);
    // e(Q_{s,r,delta}, g)
    element_pairing(right_3, Q_s_r_delta, pk.g);
    // Check if left_3 == right_3
    bool check_3 = (element_cmp(left_3, right_3) == 0);
    if (!check_3) {
        std::cout << "difference check_3: " << check_3 << std::endl;
    }

    // Clean up
    element_clear(left_3);
    element_clear(right_3);
    
    /*
    * Check 4: e(I_r, g) = e(L_r, g^r)
    */
    element_t left_4, right_4;
    element_init_GT(left_4, pk.pairing);
    element_init_GT(right_4, pk.pairing);
    
    // e(I_r, g)
    element_pairing(left_4, I_r, pk.g);
    // e(L_r, g^r)
    element_pairing(right_4, L_r, pk.g_ri[0]);
    // Check if left_4 == right_4
    bool check_4 = (element_cmp(left_4, right_4) == 0);
    if (!check_4) {
        std::cout << "difference check_4: " << check_4 << std::endl;
    }

    // Clean up
    element_clear(left_4);
    element_clear(right_4);
    
    // /*
    // * Check 5: Verify I_r = g^{\sum_{i∈I} r^i} (I_r is indeed the accumulator for the intersection set)
    // */
    // element_t computed_I_r;
    // element_init_G1(computed_I_r, pk.pairing);
    // element_set1(computed_I_r);
    
    // // Compute g^{\sum_{i∈I} r^i}
    // for (const auto& i : intersection) {
    //     if (i < 1 || i >= pk.q) continue;  // Skip invalid elements
    //     element_mul(computed_I_r, computed_I_r, pk.g_ri[i - 1]);  // Multiply by g^{r^i}
    // }
    
    // // Check if computed_I_r == I_r
    // bool check_5 = (element_cmp(computed_I_r, I_r) == 0);
    // std::cout << "difference check_5: " << check_5 << std::endl;
    
    // // Clean up
    // element_clear(computed_I_r);
    
    /*
    * Check 6: Verify D_r = A_r / I_r (accumulator for difference = accumulator for A / accumulator for intersection)
    */
    element_t expected_D_r, I_r_inv;
    element_init_G1(expected_D_r, pk.pairing);
    element_init_G1(I_r_inv, pk.pairing);
    
    // Compute the inverse of I_r
    element_set(I_r_inv, I_r);
    element_invert(I_r_inv, I_r_inv);
    
    // expected_D_r = A_r / I_r = A_r * I_r^(-1)
    element_set(expected_D_r, lhs_acc.g_r);
    element_mul(expected_D_r, expected_D_r, I_r_inv);
    
    // Check if expected_D_r == D_r
    bool check_6 = (element_cmp(expected_D_r, D_r) == 0);
    if (!check_6) {
        std::cout << "difference check_6: " << check_6 << std::endl;
    }
    
    // Clean up
    element_clear(expected_D_r);
    element_clear(I_r_inv);
    
    /*
    * Check 7: Verify D_r = g^{\sum_{i∈D} r^i} (D_r is indeed the accumulator for the difference set)
    */
    element_t computed_D_r;
    element_init_G1(computed_D_r, pk.pairing);
    element_set1(computed_D_r);
    
    // Compute g^{\sum_{i∈D} r^i}
    for (const auto& i : result) {
        if (i < 1 || i >= pk.q) continue;  // Skip invalid elements
        element_mul(computed_D_r, computed_D_r, pk.g_ri[i - 1]);  // Multiply by g^{r^i}
    }
    
    // Check if computed_D_r == D_r
    bool check_7 = (element_cmp(computed_D_r, D_r) == 0);
    if (!check_7) {
        std::cout << "difference check_7: " << check_7 << std::endl;
    }
    
    // Clean up
    element_clear(computed_D_r);
    
    // All checks must pass
    return check_1 && check_2 && check_3 && check_4 && check_6 && check_7;
}

/*
* CountProof ====================================================
*/

// CountProof implementation
CountProof::CountProof(pairing_t pairing) {
    // Initialize elements
    element_init_G1(g_a_s, pairing);
    element_set1(g_a_s);
}

// Copy constructor
CountProof::CountProof(const CountProof& other) {
    // Initialize elements in the same field
    element_init_same_as(g_a_s, const_cast<element_t&>(other.g_a_s));
    
    // Copy values
    element_set(g_a_s, const_cast<element_t&>(other.g_a_s));
}

// Assignment operator
CountProof& CountProof::operator=(const CountProof& other) {
    if (this != &other) {
        // Clear existing elements
        element_clear(g_a_s);
        
        // Initialize elements in the same field
        element_init_same_as(g_a_s, const_cast<element_t&>(other.g_a_s));
        
        // Copy values
        element_set(g_a_s, const_cast<element_t&>(other.g_a_s));
    }
    return *this;
}

CountProof::~CountProof() noexcept {
    // Clear elements
    element_clear(g_a_s);
}

bool CountProof::verify_aggr(AccValue& acc, uint64_t result, AccPublicKey& pk) {
    /*
    * Check: e(A_s / g^v, g) = e(g^{s-1}, pi)
    */
    mpz_t v;
    element_t left_tmp, left, right_tmp, right;
    mpz_init(v);
    element_init_G1(left_tmp, pk.pairing);
    element_init_GT(left, pk.pairing);
    element_init_G1(right_tmp, pk.pairing);
    element_init_GT(right, pk.pairing);
    
    // g^v
    mpz_import(v, 1, -1, sizeof(result), 0, 0, &result);
    element_pow_mpz(left_tmp, pk.g, v);
    // A_s / g^v
    element_div(left_tmp, acc.g_s, left_tmp);
    // e(A_s / g^v, g)
    element_pairing(left, left_tmp, pk.g);
    // g^{s-1}
    element_div(right_tmp, pk.g_si[0], pk.g);
    // e(g^{s-1}, pi)
    element_pairing(right, right_tmp, g_a_s);

    bool check = (element_cmp(left, right) == 0);
    if (!check) {
        std::cout << "count check: " << check << std::endl;
    }

    // Clean up
    mpz_clear(v);
    element_clear(left_tmp);
    element_clear(left);
    element_clear(right_tmp);
    element_clear(right);
    
    // All checks must pass
    return check;
}

/*
* SumProof ====================================================
*/

// SumProof implementation
SumProof::SumProof(pairing_t pairing) : count(0) {
    // Initialize elements
    element_init_G1(g_b_s, pairing);
    element_set1(g_b_s);
}

// Copy constructor
SumProof::SumProof(const SumProof& other): count(other.count) {
    // Initialize elements in the same field
    element_init_same_as(g_b_s, const_cast<element_t&>(other.g_b_s));
    
    // Copy values
    element_set(g_b_s, const_cast<element_t&>(other.g_b_s));
}

// Assignment operator
SumProof& SumProof::operator=(const SumProof& other) {
    if (this != &other) {
        // Clear existing elements
        element_clear(g_b_s);
        
        // Initialize elements in the same field
        element_init_same_as(g_b_s, const_cast<element_t&>(other.g_b_s));
        
        // Copy values
        element_set(g_b_s, const_cast<element_t&>(other.g_b_s));
        count = other.count;
    }
    return *this;
}

SumProof::~SumProof() noexcept {
    // Clear elements
    element_clear(g_b_s);
}

bool SumProof::verify_aggr(AccValue& acc, uint64_t result, AccPublicKey& pk) {
    /*
    * Check: e(A_s, g) = e(g^{(s-1)^2}, pi_1) * e(g^{v(s-1)+pi_2}, g)
    */
    mpz_t v, pi_2;
    element_t left, right, right_1, right_2, right_1_tmp, right_2_tmp, right_2_tmp2;
    mpz_init(v);
    mpz_init(pi_2);
    element_init_GT(left, pk.pairing);
    element_init_GT(right, pk.pairing);
    element_init_GT(right_1, pk.pairing);
    element_init_GT(right_2, pk.pairing);
    element_init_G1(right_1_tmp, pk.pairing);
    element_init_G1(right_2_tmp, pk.pairing);
    element_init_G1(right_2_tmp2, pk.pairing);

    // left: e(A_s, g)
    element_pairing(left, acc.g_s, pk.g);
    
    // right_1_tmp: g^{(s-1)^2} = g^{s^2} / (g^{s})^2 * g
    element_set1(right_1_tmp);
    element_mul(right_1_tmp, right_1_tmp, pk.g_si[1]);
    element_div(right_1_tmp, right_1_tmp, pk.g_si[0]);
    element_div(right_1_tmp, right_1_tmp, pk.g_si[0]);
    element_mul(right_1_tmp, right_1_tmp, pk.g);
    // right_1: e(g^{(s-1)^2}, pi_1)
    element_pairing(right_1, right_1_tmp, g_b_s);

    // right_2_tmp: g^{v(s-1)+pi_2}
    // g^{s-1}
    element_set1(right_2_tmp);
    element_mul(right_2_tmp, right_2_tmp, pk.g_si[0]);
    element_div(right_2_tmp, right_2_tmp, pk.g);
    // g^(s-1)v
    mpz_import(v, 1, -1, sizeof(result), 0, 0, &result);
    element_pow_mpz(right_2_tmp, right_2_tmp, v);
    // g^{pi_2}
    mpz_import(pi_2, 1, -1, sizeof(count), 0, 0, &count);
    element_pow_mpz(right_2_tmp2, pk.g, pi_2);
    // g^{v(s-1)+pi_2}
    element_mul(right_2_tmp, right_2_tmp, right_2_tmp2);
    // right_2: e(g^{v(s-1)+pi_2}, g)
    element_pairing(right_2, right_2_tmp, pk.g);

    // right: e(g^{(s-1)^2}, pi_1) * e(g^{v(s-1)+pi_2}, g)
    element_mul(right, right_1, right_2);

    bool check = (element_cmp(left, right) == 0);
    if (!check) {
        std::cout << "sum check: " << check << std::endl;
    }

    // Clean up
    mpz_clear(v);
    mpz_clear(pi_2);
    element_clear(left);
    element_clear(right);
    element_clear(right_1);
    element_clear(right_2);
    element_clear(right_1_tmp);
    element_clear(right_2_tmp);
    element_clear(right_2_tmp2);

    // All checks must pass
    return check;
}

/*
* MinProof ====================================================
*/

// MinProof implementation
MinProof::MinProof(pairing_t pairing) {
    // Initialize elements
    element_init_G1(pi_min, pairing);
    element_set1(pi_min);
}

// Copy constructor
MinProof::MinProof(const MinProof& other) {
    // Initialize elements in the same field
    element_init_same_as(pi_min, const_cast<element_t&>(other.pi_min));
    
    // Copy values
    element_set(pi_min, const_cast<element_t&>(other.pi_min));
}

// Assignment operator
MinProof& MinProof::operator=(const MinProof& other) {
    if (this != &other) {
        // Clear existing elements
        element_clear(pi_min);
        
        // Initialize elements in the same field
        element_init_same_as(pi_min, const_cast<element_t&>(other.pi_min));
        
        // Copy values
        element_set(pi_min, const_cast<element_t&>(other.pi_min));
    }
    return *this;
}

MinProof::~MinProof() noexcept {
    // Clear elements
    element_clear(pi_min);
}

bool MinProof::verify_aggr(AccValue& acc, uint64_t result, AccPublicKey& pk) {
    if (result < 1 || result >= pk.q) {
        // Minimum value is outside the valid range
        return false;
    }
    
    element_t left, right, right_1, right_2;
    element_init_GT(left, pk.pairing);
    element_init_GT(right, pk.pairing);
    element_init_GT(right_1, pk.pairing);
    element_init_GT(right_2, pk.pairing);

    // left: e(A_s, g)
    element_pairing(left, acc.g_s, pk.g);

    // right_1: e(g^{s^v}, g)
    element_pairing(right_1, pk.g_si[result - 1], pk.g);

    // right_2: e(pi_min, g^{s^{v+1}})
    element_pairing(right_2, pi_min, pk.g_si[result]);

    // right: e(A_s, g) == e(g^{s^v}, g) * e(pi_min, g^{s^{v+1}})
    element_mul(right, right_1, right_2);

    bool check = (element_cmp(left, right) == 0);
    if (!check) {
        std::cout << "min check: " << check << std::endl;
    }

    // Clean up
    element_clear(left);
    element_clear(right);
    element_clear(right_1);
    element_clear(right_2);
    
    return check;
}

/*
* MaxProof ====================================================
*/

// MaxProof implementation
MaxProof::MaxProof(pairing_t pairing) {
    // Initialize elements
    element_init_G1(pi_max, pairing);
    element_set1(pi_max);
}

// Copy constructor
MaxProof::MaxProof(const MaxProof& other) {
    // Initialize elements in the same field
    element_init_same_as(pi_max, const_cast<element_t&>(other.pi_max));
    
    // Copy values
    element_set(pi_max, const_cast<element_t&>(other.pi_max));
}

// Assignment operator
MaxProof& MaxProof::operator=(const MaxProof& other) {
    if (this != &other) {
        // Clear existing elements
        element_clear(pi_max);
        
        // Initialize elements in the same field
        element_init_same_as(pi_max, const_cast<element_t&>(other.pi_max));
        
        // Copy values
        element_set(pi_max, const_cast<element_t&>(other.pi_max));
    }
    return *this;
}

MaxProof::~MaxProof() noexcept {
    // Clear elements
    element_clear(pi_max);
}

bool MaxProof::verify_aggr(AccValue& acc, uint64_t result, AccPublicKey& pk) {
    if (result < 1 || result >= pk.q) {
        // Maximum value is outside the valid range
        return false;
    }
    
    element_t left, right, right_1, right_2;
    element_init_GT(left, pk.pairing);
    element_init_GT(right, pk.pairing);
    element_init_GT(right_1, pk.pairing);
    element_init_GT(right_2, pk.pairing);

    // left: e(A_{r,s}, g)
    element_pairing(left, acc.g_r_s, pk.g);

    // right_1: e(g^{r^v*s^{q-v}}, g)
    size_t idx = AccPublicKey::map_i_j_to_index(result, pk.q - result, pk.q);
    element_pairing(right_1, pk.g_ri_sj[idx], pk.g);

    // right_2: e(pi_max, g^{s^{q-v+1}})
    element_pairing(right_2, pi_max, pk.g_si[pk.q - result]);

    // right: e(A_{r,s}, g) == e(g^{r^v*s^{q-v}}, g) * e(pi_max, g^{s^{q-v+1}})
    element_mul(right, right_1, right_2);

    bool check = (element_cmp(left, right) == 0);
    if (!check) {
        std::cout << "max check: " << check << std::endl;
    }

    // Clean up
    element_clear(left);
    element_clear(right);
    element_clear(right_1);
    element_clear(right_2);
    
    return check;
}

/*
* RangeProof ====================================================
*/

RangeProof::RangeProof(pairing_t pairing) : d_B(pairing), d_D(pairing), max_B(0), min_D(0),
    pi_BC(pairing), pi_BD(pairing), pi_CD(pairing), pi_1(pairing), pi_2(pairing) {
    // Initialize elements
    element_init_G1(B_s, pairing);
    element_init_G1(B_s_alpha, pairing);
    element_init_G1(B_s_r, pairing);
    element_init_G1(B_s_r_gamma, pairing);
    element_init_G1(D_s, pairing);
    element_init_G1(D_s_alpha, pairing);
    element_init_G1(Z_s_r, pairing);
}

// Copy constructor
RangeProof::RangeProof(const RangeProof& other) : d_B(other.d_B), d_D(other.d_D), max_B(other.max_B), min_D(other.min_D),
    pi_BC(other.pi_BC), pi_BD(other.pi_BD), pi_CD(other.pi_CD), pi_1(other.pi_1), pi_2(other.pi_2) {
    // Initialize elements in the same field
    element_init_same_as(B_s, const_cast<element_t&>(other.B_s));
    element_init_same_as(B_s_alpha, const_cast<element_t&>(other.B_s_alpha));
    element_init_same_as(B_s_r, const_cast<element_t&>(other.B_s_r));
    element_init_same_as(B_s_r_gamma, const_cast<element_t&>(other.B_s_r_gamma));
    element_init_same_as(D_s, const_cast<element_t&>(other.D_s));
    element_init_same_as(D_s_alpha, const_cast<element_t&>(other.D_s_alpha));
    element_init_same_as(Z_s_r, const_cast<element_t&>(other.Z_s_r));
    
    // Copy values
    element_set(B_s, const_cast<element_t&>(other.B_s));
    element_set(B_s_alpha, const_cast<element_t&>(other.B_s_alpha));
    element_set(B_s_r, const_cast<element_t&>(other.B_s_r));
    element_set(B_s_r_gamma, const_cast<element_t&>(other.B_s_r_gamma));
    element_set(D_s, const_cast<element_t&>(other.D_s));
    element_set(D_s_alpha, const_cast<element_t&>(other.D_s_alpha));
    element_set(Z_s_r, const_cast<element_t&>(other.Z_s_r));
}

// Assignment operator
RangeProof& RangeProof::operator=(const RangeProof& other) {
    if (this != &other) {
        // Clear existing elements
        element_clear(B_s);
        element_clear(B_s_alpha);
        element_clear(B_s_r);
        element_clear(B_s_r_gamma);
        element_clear(D_s);
        element_clear(D_s_alpha);
        element_clear(Z_s_r);
        
        // Initialize elements in the same field
        element_init_same_as(B_s, const_cast<element_t&>(other.B_s));
        element_init_same_as(B_s_alpha, const_cast<element_t&>(other.B_s_alpha));
        element_init_same_as(B_s_r, const_cast<element_t&>(other.B_s_r));
        element_init_same_as(B_s_r_gamma, const_cast<element_t&>(other.B_s_r_gamma));
        element_init_same_as(D_s, const_cast<element_t&>(other.D_s));
        element_init_same_as(D_s_alpha, const_cast<element_t&>(other.D_s_alpha));
        element_init_same_as(Z_s_r, const_cast<element_t&>(other.Z_s_r));
        
        // Copy values
        element_set(B_s, const_cast<element_t&>(other.B_s));
        element_set(B_s_alpha, const_cast<element_t&>(other.B_s_alpha));
        element_set(B_s_r, const_cast<element_t&>(other.B_s_r));
        element_set(B_s_r_gamma, const_cast<element_t&>(other.B_s_r_gamma));
        element_set(D_s, const_cast<element_t&>(other.D_s));
        element_set(D_s_alpha, const_cast<element_t&>(other.D_s_alpha));
        element_set(Z_s_r, const_cast<element_t&>(other.Z_s_r));
        
        max_B = other.max_B;
        min_D = other.min_D;

        pi_BC = other.pi_BC;
        pi_BD = other.pi_BD;
        pi_CD = other.pi_CD;
        pi_1 = other.pi_1;
        pi_2 = other.pi_2;
    }
    return *this;
}

RangeProof::~RangeProof() noexcept {
    // Clear elements
    element_clear(B_s);
    element_clear(B_s_alpha);
    element_clear(B_s_r);
    element_clear(B_s_r_gamma);
    element_clear(D_s);
    element_clear(D_s_alpha);
    element_clear(Z_s_r);
}

bool RangeProof::verify_range(AccValue& acc, const Set& result, uint64_t l, uint64_t r, AccPublicKey& pk) {
    /*
    * Check 1.1: e(B_s, g^alpha) = e(B_s_alpha, g)
    */
    element_t left, right;
    element_init_GT(left, pk.pairing);
    element_init_GT(right, pk.pairing);

    // left: e(B_s, g^alpha)
    element_pairing(left, B_s, pk.g_alpha);
    // right: e(B_s_alpha, g)
    element_pairing(right, B_s_alpha, pk.g);

    bool check_1_1 = (element_cmp(left, right) == 0);
    if (!check_1_1) {
        std::cout << "range check_1_1: " << check_1_1 << std::endl;
    }
    
    /*
    * Check 1.2: e(B_s_r, g^gamma) = e(B_s_r_gamma, g)
    */
    // left: e(B_s_r, g^gamma)
    element_pairing(left, B_s_r, pk.g_gamma);
    // right: e(B_s_r_gamma, g)
    element_pairing(right, B_s_r_gamma, pk.g);

    bool check_1_2 = (element_cmp(left, right) == 0);
    if (!check_1_2) {
        std::cout << "range check_1_2: " << check_1_2 << std::endl;
    }

    /*
    * Check 1.3: e(D_s, g^alpha) = e(D_s_alpha, g)
    */
    // left: e(D_s, g^alpha)
    element_pairing(left, D_s, pk.g_alpha);
    // right: e(D_s_alpha, g)
    element_pairing(right, D_s_alpha, pk.g);

    bool check_1_3 = (element_cmp(left, right) == 0);
    if (!check_1_3) {
        std::cout << "range check_1_3: " << check_1_3 << std::endl;
    }

    /*
    * Check 2: e(B_s / B_s_r, g) = e(g^{r-1}, Z_s_r)
    */
    element_t left_tmp, right_tmp;
    element_init_G1(left_tmp, pk.pairing);
    element_init_G1(right_tmp, pk.pairing);

    // left_tmp: B_s / B_s_r
    element_div(left_tmp, B_s, B_s_r);
    // left: e(B_s / B_s_r, g)
    element_pairing(left, left_tmp, pk.g);
   
    // right_tmp: g^{r-1}
    element_div(right_tmp, pk.g_ri[0], pk.g);
    // right: e(g^{r-1}, Z_s_r)
    element_pairing(right, right_tmp, Z_s_r);

    bool check_2 = (element_cmp(left, right) == 0);
    if (!check_2) {
        std::cout << "range check_2: " << check_2 << std::endl;
    }

    /*
    * Check 3: verify pi_BC
    */
    AccValue d_C = AccValue::setup(result, pk);
    // empty set
    Set empty_set;
    bool check_3 = pi_BC.verify(d_B, d_C, empty_set, pk);
    if (!check_3) {
        std::cout << "range check_3: " << check_3 << std::endl;
    }

    /*
    * Check 4: verify pi_BD
    */
    bool check_4 = pi_BD.verify(d_B, d_D, empty_set, pk);
    if (!check_4) {
        std::cout << "range check_4: " << check_4 << std::endl;
    }

    /*
    * Check 5: verify pi_CD
    */
    bool check_5 = pi_CD.verify(d_C, d_D, empty_set, pk);
    if (!check_5) {
        std::cout << "range check_5: " << check_5 << std::endl;
    }
    
    /*
    * Check 6: B_s * C_s * D_s = A_s
    */
    // left: B_s * C_s
    element_mul(left_tmp, B_s, d_C.g_s);
    // left: B_s * C_s * D_s
    element_mul(left_tmp, left_tmp, d_D.g_s);
    bool check_6 = (element_cmp(left_tmp, acc.g_s) == 0);
    if (!check_6) {
        std::cout << "range check_6: " << check_6 << std::endl;
    }

    /*
    * Check 7: verify pi_1
    */
   bool check_7;
    if (!element_is1(pi_1.pi_max)) {
        check_7 = pi_1.verify_aggr(d_B, max_B, pk);
    } else {
        check_7 = max_B == 0;
    }
    if (!check_7) {
        std::cout << "range check_7: " << check_7 << std::endl;
    }
    
    /*
    * Check 8: verify pi_2
    */
    bool check_8;
    if (!element_is1(pi_2.pi_min)) {
        check_8 = pi_2.verify_aggr(d_D, min_D, pk);
    } else {
        check_8 = min_D == 0;
    }
    if (!check_8) {
        std::cout << "range check_8: " << check_8 << std::endl;
    }

    /*
    * Check 9: max_B < l, min_C >= l, max_C <= r, min_D > r
    */
    uint64_t min_C = UINT64_MAX, max_C = 0;
    for (const auto& i : result) {
        if (i < min_C) min_C = i;
        if (i > max_C) max_C = i;
    }
    
    bool check_9 = (max_B == 0 || max_B < l) && (min_C >= l) && (max_C <= r) && (min_D == 0 || min_D > r);
    if (!check_9) {
        std::cout << "range check_9: " << check_9 << std::endl;
    }


    element_clear(left);
    element_clear(right);
    element_clear(left_tmp);
    element_clear(right_tmp);

    return check_1_1 && check_1_2 && check_1_3 && check_2 && check_3 
    && check_4 && check_5 && check_6 && check_7 && check_8 && check_9;
}

/*
* NestedIntersectionProof ====================================================
*/

NestedIntersectionProof::NestedIntersectionProof(pairing_t pairing) {
    // Initialize elements
    element_init_G1(I_r_beta, pairing);
    element_init_G1(Q_s_r, pairing);
    element_init_G1(Q_s_r_delta, pairing);
    element_init_G1(L_r, pairing);

    element_init_G1(Z_s_r, pairing);
    element_init_G1(I_s_r_gamma, pairing);

    element_init_G1(I_s_alpha, pairing);
    element_init_G1(Q_r_s, pairing);
    element_init_G1(Q_r_s_delta, pairing);
    element_init_G1(L_s, pairing);

    element_init_G1(Z_r_s, pairing);
    element_init_G1(I_r_s_gamma, pairing);
}

// Copy constructor
NestedIntersectionProof::NestedIntersectionProof(const NestedIntersectionProof& other) {
    // Initialize elements in the same field
    element_init_same_as(I_r_beta, const_cast<element_t&>(other.I_r_beta));
    element_init_same_as(Q_s_r, const_cast<element_t&>(other.Q_s_r));
    element_init_same_as(Q_s_r_delta, const_cast<element_t&>(other.Q_s_r_delta));
    element_init_same_as(L_r, const_cast<element_t&>(other.L_r));

    element_init_same_as(Z_s_r, const_cast<element_t&>(other.Z_s_r));
    element_init_same_as(I_s_r_gamma, const_cast<element_t&>(other.I_s_r_gamma));   
    element_init_same_as(I_s_alpha, const_cast<element_t&>(other.I_s_alpha));
    element_init_same_as(Q_r_s, const_cast<element_t&>(other.Q_r_s));
    element_init_same_as(Q_r_s_delta, const_cast<element_t&>(other.Q_r_s_delta));
    element_init_same_as(L_s, const_cast<element_t&>(other.L_s));
    element_init_same_as(Z_r_s, const_cast<element_t&>(other.Z_r_s));
    element_init_same_as(I_r_s_gamma, const_cast<element_t&>(other.I_r_s_gamma));   
    
    // Copy values
    element_set(I_r_beta, const_cast<element_t&>(other.I_r_beta));
    element_set(Q_s_r, const_cast<element_t&>(other.Q_s_r));
    element_set(Q_s_r_delta, const_cast<element_t&>(other.Q_s_r_delta));
    element_set(L_r, const_cast<element_t&>(other.L_r));
    element_set(Z_s_r, const_cast<element_t&>(other.Z_s_r));
    element_set(I_s_r_gamma, const_cast<element_t&>(other.I_s_r_gamma));
    element_set(I_s_alpha, const_cast<element_t&>(other.I_s_alpha));
    element_set(Q_r_s, const_cast<element_t&>(other.Q_r_s));
    element_set(Q_r_s_delta, const_cast<element_t&>(other.Q_r_s_delta));
    element_set(L_s, const_cast<element_t&>(other.L_s));
    element_set(Z_r_s, const_cast<element_t&>(other.Z_r_s));
    element_set(I_r_s_gamma, const_cast<element_t&>(other.I_r_s_gamma));
}

// Assignment operator
NestedIntersectionProof& NestedIntersectionProof::operator=(const NestedIntersectionProof& other) {
    if (this != &other) {
        // Clear existing elements
        element_clear(I_r_beta);
        element_clear(Q_s_r);
        element_clear(Q_s_r_delta);
        element_clear(L_r);
        element_clear(Z_s_r);
        element_clear(I_s_r_gamma);
        element_clear(I_s_alpha);
        element_clear(Q_r_s);
        element_clear(Q_r_s_delta);
        element_clear(L_s);
        element_clear(Z_r_s);
        element_clear(I_r_s_gamma);
        
        // Initialize elements in the same field
        element_init_same_as(I_r_beta, const_cast<element_t&>(other.I_r_beta));
        element_init_same_as(Q_s_r, const_cast<element_t&>(other.Q_s_r));
        element_init_same_as(Q_s_r_delta, const_cast<element_t&>(other.Q_s_r_delta));
        element_init_same_as(L_r, const_cast<element_t&>(other.L_r));
        element_init_same_as(Z_s_r, const_cast<element_t&>(other.Z_s_r));
        element_init_same_as(I_s_r_gamma, const_cast<element_t&>(other.I_s_r_gamma));
        element_init_same_as(I_s_alpha, const_cast<element_t&>(other.I_s_alpha));
        element_init_same_as(Q_r_s, const_cast<element_t&>(other.Q_r_s));
        element_init_same_as(Q_r_s_delta, const_cast<element_t&>(other.Q_r_s_delta));
        element_init_same_as(L_s, const_cast<element_t&>(other.L_s));
        element_init_same_as(Z_r_s, const_cast<element_t&>(other.Z_r_s));
        element_init_same_as(I_r_s_gamma, const_cast<element_t&>(other.I_r_s_gamma));
        
        // Copy values
        element_set(I_r_beta, const_cast<element_t&>(other.I_r_beta));
        element_set(Q_s_r, const_cast<element_t&>(other.Q_s_r));
        element_set(Q_s_r_delta, const_cast<element_t&>(other.Q_s_r_delta));
        element_set(L_r, const_cast<element_t&>(other.L_r));
        element_set(Z_s_r, const_cast<element_t&>(other.Z_s_r));
        element_set(I_s_r_gamma, const_cast<element_t&>(other.I_s_r_gamma));
        element_set(I_s_alpha, const_cast<element_t&>(other.I_s_alpha));
        element_set(Q_r_s, const_cast<element_t&>(other.Q_r_s));
        element_set(Q_r_s_delta, const_cast<element_t&>(other.Q_r_s_delta));
        element_set(L_s, const_cast<element_t&>(other.L_s));
        element_set(Z_r_s, const_cast<element_t&>(other.Z_r_s));
        element_set(I_r_s_gamma, const_cast<element_t&>(other.I_r_s_gamma));
    }
    return *this;
}


NestedIntersectionProof::~NestedIntersectionProof() noexcept {
    // Clear elements
    element_clear(I_r_beta);
    element_clear(Q_s_r);
    element_clear(Q_s_r_delta);
    element_clear(L_r);

    element_clear(Z_s_r);
    element_clear(I_s_r_gamma);

    element_clear(I_s_alpha);
    element_clear(Q_r_s);
    element_clear(Q_r_s_delta);
    element_clear(L_s);

    element_clear(Z_r_s);
    element_clear(I_r_s_gamma);
}

bool NestedIntersectionProof::verify_nested(AccValue& lhs_acc, AccValue& rhs_acc, AccValue& res_acc, AccPublicKey& pk) {
    /*
    * Check 1: e(A_s, B_{r,s}) = e(I_r, g^{s^q}) × e(Q_{s,r}, g)
    */
    element_t left_1, right_1, right_1_1, right_1_2;
    element_init_GT(left_1, pk.pairing);
    element_init_GT(right_1, pk.pairing);
    element_init_GT(right_1_1, pk.pairing);
    element_init_GT(right_1_2, pk.pairing);
    
    // e(A_s, B_{r,s})
    element_pairing(left_1, lhs_acc.g_s, rhs_acc.g_r_s);
    // e(I_r, g^{s^q})
    element_pairing(right_1_1, res_acc.g_r, pk.g_sq);
    // e(Q_{s,r}, g)
    element_pairing(right_1_2, Q_s_r, pk.g);
    // Combine: right_1 = right_1_1 × right_1_2
    element_mul(right_1, right_1_1, right_1_2);
    // Check if left_1 == right_1
    bool SR_check_1 = (element_cmp(left_1, right_1) == 0);
    if (!SR_check_1) {
        std::cout << "nested S,R intersection check_1: " << SR_check_1 << std::endl;
    }

    // Clean up
    element_clear(left_1);
    element_clear(right_1);
    element_clear(right_1_1);
    element_clear(right_1_2);
    
    /*
    * Check 2: e(I_r, g^alpha) = e(I_{r,alpha}, g)
    */
    element_t left_2, right_2;
    element_init_GT(left_2, pk.pairing);
    element_init_GT(right_2, pk.pairing);
    
    // e(I_r, g^alpha)
    element_pairing(left_2, res_acc.g_r, pk.g_beta);
    // e(I_{r,alpha}, g)
    element_pairing(right_2, I_r_beta, pk.g);
    // Check if left_2 == right_2
    bool SR_check_2 = (element_cmp(left_2, right_2) == 0);
    if (!SR_check_2) {
        std::cout << "nested S,R intersection check_2: " << SR_check_2 << std::endl;
    }

    // Clean up
    element_clear(left_2);
    element_clear(right_2);
    
    /*
    * Check 3: e(Q_{s,r}, g^delta) = e(Q_{s,r,delta}, g)
    */
    element_t left_3, right_3;
    element_init_GT(left_3, pk.pairing);
    element_init_GT(right_3, pk.pairing);
    
    // e(Q_{s,r}, g^delta)
    element_pairing(left_3, Q_s_r, pk.g_delta);
    // e(Q_{s,r,delta}, g)
    element_pairing(right_3, Q_s_r_delta, pk.g);
    // Check if left_3 == right_3
    bool SR_check_3 = (element_cmp(left_3, right_3) == 0);
    if (!SR_check_3) {
        std::cout << "nested S,R intersection check_3: " << SR_check_3 << std::endl;
    }

    // Clean up
    element_clear(left_3);
    element_clear(right_3);
    
    /*
    * Check 4: e(I_r, g) = e(L_r, g^r)
    */
    element_t left_4, right_4;
    element_init_GT(left_4, pk.pairing);
    element_init_GT(right_4, pk.pairing);
    
    // e(I_r, g)
    element_pairing(left_4, res_acc.g_r, pk.g);
    // e(L_r, g^r)
    element_pairing(right_4, L_r, pk.g_ri[0]);
    // Check if left_4 == right_4
    bool SR_check_4 = (element_cmp(left_4, right_4) == 0);
    if (!SR_check_4) {
        std::cout << "nested S,R intersection check_4: " << SR_check_4 << std::endl;
    }

    // Clean up
    element_clear(left_4);
    element_clear(right_4);
    
    /*
    * Check 5: e(I_r / I_r_s, g) = e(g^{s-1}, Z_s_r)
    */
    element_t left_5, right_5, left_5_1, right_5_1;
    element_init_GT(left_5, pk.pairing);
    element_init_GT(right_5, pk.pairing);
    element_init_G1(left_5_1, pk.pairing);
    element_init_G1(right_5_1, pk.pairing);

    // left_5_1 = I_r / I_r_s
    element_div(left_5_1, res_acc.g_r, res_acc.g_r_s);
    // left_5 = e(I_r / I_r_s, g)
    element_pairing(left_5, left_5_1, pk.g);

    // right_5_1 = g^{s-1}
    element_div(right_5_1, pk.g_si[0], pk.g);
    // right_5 = e(g^{s-1}, Z_s_r)
    element_pairing(right_5, right_5_1, Z_s_r);
    
    // Check if left_5 == right_5
    bool SR_check_5 = (element_cmp(left_5, right_5) == 0);
    if (!SR_check_5) {
        std::cout << "nested S,R intersection check_5: " << SR_check_5 << std::endl;
    }

    // Clean up
    element_clear(left_5);
    element_clear(right_5);
    element_clear(left_5_1);
    element_clear(right_5_1);
    
    /*
    * Check 6: e(I_s_r, g^gamma) = e(I_s_r_gamma, g)
    */
    element_t left_6, right_6;
    element_init_GT(left_6, pk.pairing);
    element_init_GT(right_6, pk.pairing);
    
    // e(I_s_r, g^gamma)
    element_pairing(left_6, res_acc.g_s_r, pk.g_gamma);
    // e(I_s_r_gamma, g)
    element_pairing(right_6, I_s_r_gamma, pk.g);
    
    // Check if left_6 == right_6
    bool SR_check_6 = (element_cmp(left_6, right_6) == 0);
    if (!SR_check_6) {
        std::cout << "nested S,R intersection check_6: " << SR_check_6 << std::endl;
    }

    // Clean up
    element_clear(left_6);
    element_clear(right_6);

    // R, S ======================================================================
    /*
    * Check 1: e(A_r, B_{s,r}) = e(I_s, g^{r^q}) × e(Q_{r,s}, g)
    */
    // element_t left_1, right_1, right_1_1, right_1_2;
    element_init_GT(left_1, pk.pairing);
    element_init_GT(right_1, pk.pairing);
    element_init_GT(right_1_1, pk.pairing);
    element_init_GT(right_1_2, pk.pairing);
    
    // e(A_r, B_{s,r})
    element_pairing(left_1, lhs_acc.g_r, rhs_acc.g_s_r);
    // e(I_s, g^{r^q})
    element_pairing(right_1_1, res_acc.g_s, pk.g_rq);
    // e(Q_{r,s}, g)
    element_pairing(right_1_2, Q_r_s, pk.g);
    // Combine: right_1 = right_1_1 × right_1_2
    element_mul(right_1, right_1_1, right_1_2);
    // Check if left_1 == right_1
    bool RS_check_1 = (element_cmp(left_1, right_1) == 0);
    if (!RS_check_1) {
        std::cout << "nested R,S intersection check_1: " << RS_check_1 << std::endl;
    }

    // Clean up
    element_clear(left_1);
    element_clear(right_1);
    element_clear(right_1_1);
    element_clear(right_1_2);
    
    /*
    * Check 2: e(I_s, g^alpha) = e(I_{s,alpha}, g)
    */
    // element_t left_2, right_2;
    element_init_GT(left_2, pk.pairing);
    element_init_GT(right_2, pk.pairing);
    
    // e(I_s, g^alpha)
    element_pairing(left_2, res_acc.g_s, pk.g_alpha);
    // e(I_{s,alpha}, g)
    element_pairing(right_2, I_s_alpha, pk.g);
    // Check if left_2 == right_2
    bool RS_check_2 = (element_cmp(left_2, right_2) == 0);
    if (!RS_check_2) {
        std::cout << "nested R,S intersection check_2: " << RS_check_2 << std::endl;
    }

    // Clean up
    element_clear(left_2);
    element_clear(right_2);
    
    /*
    * Check 3: e(Q_{r,s}, g^delta) = e(Q_{r,s,delta}, g)
    */
    // element_t left_3, right_3;
    element_init_GT(left_3, pk.pairing);
    element_init_GT(right_3, pk.pairing);
    
    // e(Q_{r,s}, g^delta)
    element_pairing(left_3, Q_r_s, pk.g_delta);
    // e(Q_{r,s,delta}, g)
    element_pairing(right_3, Q_r_s_delta, pk.g);
    // Check if left_3 == right_3
    bool RS_check_3 = (element_cmp(left_3, right_3) == 0);
    if (!RS_check_3) {
        std::cout << "nested R,S intersection check_3: " << RS_check_3 << std::endl;
    }

    // Clean up
    element_clear(left_3);
    element_clear(right_3);
    
    /*
    * Check 4: e(I_s, g) = e(L_s, g^s)
    */
    // element_t left_4, right_4;
    element_init_GT(left_4, pk.pairing);
    element_init_GT(right_4, pk.pairing);
    
    // e(I_s, g)
    element_pairing(left_4, res_acc.g_s, pk.g);
    // e(L_s, g^s)
    element_pairing(right_4, L_s, pk.g_si[0]);
    // Check if left_4 == right_4
    bool RS_check_4 = (element_cmp(left_4, right_4) == 0);
    if (!RS_check_4) {
        std::cout << "nested R,S intersection check_4: " << RS_check_4 << std::endl;
    }

    // Clean up
    element_clear(left_4);
    element_clear(right_4);
    
    /*
    * Check 5: e(I_s / I_s_r, g) = e(g^{r-1}, Z_r_s)
    */
    // element_t left_5, right_5, left_5_1, right_5_1;
    element_init_GT(left_5, pk.pairing);
    element_init_GT(right_5, pk.pairing);
    element_init_G1(left_5_1, pk.pairing);
    element_init_G1(right_5_1, pk.pairing);

    // left_5_1 = I_s / I_s_r
    element_div(left_5_1, res_acc.g_s, res_acc.g_s_r);
    // left_5 = e(I_s / I_s_r, g)
    element_pairing(left_5, left_5_1, pk.g);

    // right_5_1 = g^{r-1}
    element_div(right_5_1, pk.g_ri[0], pk.g);
    // right_5 = e(g^{r-1}, Z_r_s)
    element_pairing(right_5, right_5_1, Z_r_s);
    
    // Check if left_5 == right_5
    bool RS_check_5 = (element_cmp(left_5, right_5) == 0);
    if (!RS_check_5) {
        std::cout << "nested R,S intersection check_5: " << RS_check_5 << std::endl;
    }

    // Clean up
    element_clear(left_5);
    element_clear(right_5);
    element_clear(left_5_1);
    element_clear(right_5_1);
    
    /*
    * Check 6: e(I_r_s, g^gamma) = e(I_r_s_gamma, g)
    */
    // element_t left_6, right_6;
    element_init_GT(left_6, pk.pairing);
    element_init_GT(right_6, pk.pairing);
    
    // e(I_r_s, g^gamma)
    element_pairing(left_6, res_acc.g_r_s, pk.g_gamma);
    // e(I_r_s_gamma, g)
    element_pairing(right_6, I_r_s_gamma, pk.g);
    
    // Check if left_6 == right_6
    bool RS_check_6 = (element_cmp(left_6, right_6) == 0);
    if (!RS_check_6) {
        std::cout << "nested R,S intersection check_6: " << RS_check_6 << std::endl;
    }

    // Clean up
    element_clear(left_6);
    element_clear(right_6);
    
    // All checks must pass
    return SR_check_1 && SR_check_2 && SR_check_3 && SR_check_4 && SR_check_5 && SR_check_6
    && RS_check_1 && RS_check_2 && RS_check_3 && RS_check_4 && RS_check_5 && RS_check_6;
}


/*
* NestedUnionProof ====================================================
*/

NestedUnionProof::NestedUnionProof(pairing_t pairing) {
    element_init_G1(I_r, pairing);
    element_init_G1(I_r_beta, pairing);
    element_init_G1(Q_s_r, pairing);
    element_init_G1(Q_s_r_delta, pairing);
    element_init_G1(L_r, pairing);

    element_init_G1(Z_s_r, pairing);
    element_init_G1(U_s_r_gamma, pairing);
    element_init_G1(Z_r_s, pairing);
    element_init_G1(U_r_s_gamma, pairing);
}

size_t NestedUnionProof::getSizeADS() {
    return (element_length_in_bytes(I_r) + element_length_in_bytes(I_r_beta) + element_length_in_bytes(Q_s_r) + 
    element_length_in_bytes(Q_s_r_delta) + element_length_in_bytes(L_r) + element_length_in_bytes(Z_s_r) + 
    element_length_in_bytes(U_s_r_gamma) + element_length_in_bytes(Z_r_s) + element_length_in_bytes(U_r_s_gamma)) * sizeof(char);
}

// Copy constructor
NestedUnionProof::NestedUnionProof(const NestedUnionProof& other) {
    element_init_same_as(I_r, const_cast<element_t&>(other.I_r));
    element_init_same_as(I_r_beta, const_cast<element_t&>(other.I_r_beta));
    element_init_same_as(Q_s_r, const_cast<element_t&>(other.Q_s_r));
    element_init_same_as(Q_s_r_delta, const_cast<element_t&>(other.Q_s_r_delta));
    element_init_same_as(L_r, const_cast<element_t&>(other.L_r));
    
    element_init_same_as(Z_s_r, const_cast<element_t&>(other.Z_s_r));
    element_init_same_as(U_s_r_gamma, const_cast<element_t&>(other.U_s_r_gamma));
    element_init_same_as(Z_r_s, const_cast<element_t&>(other.Z_r_s));
    element_init_same_as(U_r_s_gamma, const_cast<element_t&>(other.U_r_s_gamma));
    
    // Copy values
    element_set(I_r, const_cast<element_t&>(other.I_r));
    element_set(I_r_beta, const_cast<element_t&>(other.I_r_beta));
    element_set(Q_s_r, const_cast<element_t&>(other.Q_s_r));
    element_set(Q_s_r_delta, const_cast<element_t&>(other.Q_s_r_delta));
    element_set(L_r, const_cast<element_t&>(other.L_r));
    element_set(Z_s_r, const_cast<element_t&>(other.Z_s_r));
    element_set(U_s_r_gamma, const_cast<element_t&>(other.U_s_r_gamma));
    element_set(Z_r_s, const_cast<element_t&>(other.Z_r_s));
    element_set(U_r_s_gamma, const_cast<element_t&>(other.U_r_s_gamma));
}

// Assignment operator
NestedUnionProof& NestedUnionProof::operator=(const NestedUnionProof& other) {
    if (this != &other) {
        // Clear existing elements
        element_clear(I_r);
        element_clear(I_r_beta);
        element_clear(Q_s_r);
        element_clear(Q_s_r_delta);
        element_clear(L_r);
        element_clear(Z_s_r);
        element_clear(U_s_r_gamma);
        element_clear(Z_r_s);
        element_clear(U_r_s_gamma);
        
        // Initialize elements in the same field
        element_init_same_as(I_r, const_cast<element_t&>(other.I_r));
        element_init_same_as(I_r_beta, const_cast<element_t&>(other.I_r_beta));
        element_init_same_as(Q_s_r, const_cast<element_t&>(other.Q_s_r));
        element_init_same_as(Q_s_r_delta, const_cast<element_t&>(other.Q_s_r_delta));
        element_init_same_as(L_r, const_cast<element_t&>(other.L_r));
        element_init_same_as(Z_s_r, const_cast<element_t&>(other.Z_s_r));
        element_init_same_as(U_s_r_gamma, const_cast<element_t&>(other.U_s_r_gamma));
        element_init_same_as(Z_r_s, const_cast<element_t&>(other.Z_r_s));
        element_init_same_as(U_r_s_gamma, const_cast<element_t&>(other.U_r_s_gamma));
        
        // Copy values
        element_set(I_r, const_cast<element_t&>(other.I_r));
        element_set(I_r_beta, const_cast<element_t&>(other.I_r_beta));
        element_set(Q_s_r, const_cast<element_t&>(other.Q_s_r));
        element_set(Q_s_r_delta, const_cast<element_t&>(other.Q_s_r_delta));
        element_set(L_r, const_cast<element_t&>(other.L_r));
        element_set(Z_s_r, const_cast<element_t&>(other.Z_s_r));
        element_set(U_s_r_gamma, const_cast<element_t&>(other.U_s_r_gamma));
        element_set(Z_r_s, const_cast<element_t&>(other.Z_r_s));
        element_set(U_r_s_gamma, const_cast<element_t&>(other.U_r_s_gamma));
    }
    return *this;
}

NestedUnionProof::~NestedUnionProof() noexcept {
    element_clear(I_r);
    element_clear(I_r_beta);
    element_clear(Q_s_r);
    element_clear(Q_s_r_delta);
    element_clear(L_r);

    element_clear(Z_s_r);
    element_clear(U_s_r_gamma);
    element_clear(Z_r_s);
    element_clear(U_r_s_gamma);
}

bool NestedUnionProof::verify_nested(AccValue& lhs_acc, AccValue& rhs_acc, AccValue& res_acc, AccPublicKey& pk) {
    // First verify the intersection part of the proof
    /*
    * Check 1: e(A_s, B_{r,s}) = e(I_r, g^{s^q}) × e(Q_{s,r}, g)
    */
    element_t left_1, right_1, right_1_1, right_1_2;
    element_init_GT(left_1, pk.pairing);
    element_init_GT(right_1, pk.pairing);
    element_init_GT(right_1_1, pk.pairing);
    element_init_GT(right_1_2, pk.pairing);
    
    // e(A_s, B_{r,s})
    element_pairing(left_1, lhs_acc.g_s, rhs_acc.g_r_s);
    // e(I_r, g^{s^q})
    element_pairing(right_1_1, I_r, pk.g_sq);
    // e(Q_{s,r}, g)
    element_pairing(right_1_2, Q_s_r, pk.g);
    // Combine: right_1 = right_1_1 × right_1_2
    element_mul(right_1, right_1_1, right_1_2);
    // Check if left_1 == right_1
    bool check_1 = (element_cmp(left_1, right_1) == 0);
    if (!check_1) {
        std::cout << "nested union check_1: " << check_1 << std::endl;
    }

    // Clean up
    element_clear(left_1);
    element_clear(right_1);
    element_clear(right_1_1);
    element_clear(right_1_2);
    
    /*
    * Check 2: e(I_r, g^beta) = e(I_{r,beta}, g)
    */
    element_t left_2, right_2;
    element_init_GT(left_2, pk.pairing);
    element_init_GT(right_2, pk.pairing);
    
    // e(I_r, g^beta)
    element_pairing(left_2, I_r, pk.g_beta);
    // e(I_{r,beta}, g)
    element_pairing(right_2, I_r_beta, pk.g);
    // Check if left_2 == right_2
    bool check_2 = (element_cmp(left_2, right_2) == 0);
    if (!check_2) {
        std::cout << "nested union check_2: " << check_2 << std::endl;
    }

    // Clean up
    element_clear(left_2);
    element_clear(right_2);
    
    /*
    * Check 3: e(Q_{s,r}, g^delta) = e(Q_{s,r,delta}, g)
    */
    element_t left_3, right_3;
    element_init_GT(left_3, pk.pairing);
    element_init_GT(right_3, pk.pairing);
    
    // e(Q_{s,r}, g^delta)
    element_pairing(left_3, Q_s_r, pk.g_delta);
    // e(Q_{s,r,delta}, g)
    element_pairing(right_3, Q_s_r_delta, pk.g);
    // Check if left_3 == right_3
    bool check_3 = (element_cmp(left_3, right_3) == 0);
    if (!check_3) {
        std::cout << "nested union check_3: " << check_3 << std::endl;
    }

    // Clean up
    element_clear(left_3);
    element_clear(right_3);
    
    /*
    * Check 4: e(I_r, g) = e(L_r, g^r)
    */
    element_t left_4, right_4;
    element_init_GT(left_4, pk.pairing);
    element_init_GT(right_4, pk.pairing);
    
    // e(I_r, g)
    element_pairing(left_4, I_r, pk.g);
    // e(L_r, g^r)
    element_pairing(right_4, L_r, pk.g_ri[0]);
    // Check if left_4 == right_4
    bool check_4 = (element_cmp(left_4, right_4) == 0);
    if (!check_4) {
        std::cout << "nested union check_4: " << check_4 << std::endl;
    }

    // Clean up
    element_clear(left_4);
    element_clear(right_4);
    
    // /*
    // * Check 5: Verify I_r = g^{\sum_{i∈I} r^i} (I_r is indeed the accumulator for the intersection set)
    // */
    // element_t computed_I_r;
    // element_init_G1(computed_I_r, pk.pairing);
    // element_set1(computed_I_r);
    
    // // Compute g^{\sum_{i∈I} r^i}
    // for (const auto& i : intersection) {
    //     if (i < 1 || i >= pk.q) continue;  // Skip invalid elements
    //     element_mul(computed_I_r, computed_I_r, pk.g_ri[i - 1]);  // Multiply by g^{r^i}
    // }
    
    // // Check if computed_I_r == I_r
    // bool check_5 = (element_cmp(computed_I_r, I_r) == 0);
    // std::cout << "difference check_5: " << check_5 << std::endl;
    
    // // Clean up
    // element_clear(computed_I_r);
    
    /*
    * Check 6: Verify U_r = A_r * B_r / I_r
    */
    element_t expected_U_r, I_r_inv;
    element_init_G1(expected_U_r, pk.pairing);
    element_init_G1(I_r_inv, pk.pairing);
    
    // Compute the inverse of I_r
    element_set(I_r_inv, I_r);
    element_invert(I_r_inv, I_r_inv);
    
    // expected_U_r = A_r * B_r / I_r = A_r * B_r * I_r^(-1)
    element_mul(expected_U_r, lhs_acc.g_r, rhs_acc.g_r);
    element_mul(expected_U_r, expected_U_r, I_r_inv);
    
    // Check if expected_U_r == U_r
    bool check_6 = (element_cmp(expected_U_r, res_acc.g_r) == 0);
    if (!check_6) {
        std::cout << "nested union check_6: " << check_6 << std::endl;
    }
    
    // Clean up
    element_clear(expected_U_r);
    element_clear(I_r_inv);
    
    // /*
    // * Check 7: Verify D_r = g^{\sum_{i∈D} r^i} (D_r is indeed the accumulator for the difference set)
    // */
    // element_t computed_D_r;
    // element_init_G1(computed_D_r, pk.pairing);
    // element_set1(computed_D_r);
    
    // // Compute g^{\sum_{i∈D} r^i}
    // for (const auto& i : result) {
    //     if (i < 1 || i >= pk.q) continue;  // Skip invalid elements
    //     element_mul(computed_D_r, computed_D_r, pk.g_ri[i - 1]);  // Multiply by g^{r^i}
    // }
    
    // // Check if computed_D_r == D_r
    // bool check_7 = (element_cmp(computed_D_r, D_r) == 0);
    // std::cout << "difference check_7: " << check_7 << std::endl;
    
    // // Clean up
    // element_clear(computed_D_r);

    /*
    * Check 8: e(U_r / U_r_s, g) = e(g^{s-1}, Z_r_s)
    */
    element_t left_8, right_8, left_8_1, right_8_1;
    element_init_GT(left_8, pk.pairing);
    element_init_GT(right_8, pk.pairing);
    element_init_G1(left_8_1, pk.pairing);
    element_init_G1(right_8_1, pk.pairing);

    // left_8_1 = U_r / U_r_s
    element_div(left_8_1, res_acc.g_r, res_acc.g_r_s);
    // left_8 = e(U_r / U_r_s, g)
    element_pairing(left_8, left_8_1, pk.g);

    // right_5_1 = g^{s-1}
    element_div(right_8_1, pk.g_si[0], pk.g);
    // right_8 = e(g^{s-1}, Z_r_s)
    element_pairing(right_8, right_8_1, Z_r_s);
    
    // Check if left_8 == right_8
    bool check_8 = (element_cmp(left_8, right_8) == 0);
    if (!check_8) {
        std::cout << "nested union check_8: " << check_8 << std::endl;
    }

    // Clean up
    element_clear(left_8);
    element_clear(right_8);
    element_clear(left_8_1);
    element_clear(right_8_1);
    
    /*
    * Check 9: e(U_s_r, g^gamma) = e(U_s_r_gamma, g)
    */
    element_t left_9, right_9;
    element_init_GT(left_9, pk.pairing);
    element_init_GT(right_9, pk.pairing);
    
    // e(U_s_r, g^gamma)
    element_pairing(left_9, res_acc.g_s_r, pk.g_gamma);
    // e(U_s_r_gamma, g)
    element_pairing(right_9, U_s_r_gamma, pk.g);
    
    // Check if left_9 == right_9
    bool check_9 = (element_cmp(left_9, right_9) == 0);
    if (!check_9) {
        std::cout << "nested union check_9: " << check_9 << std::endl;
    }

    // Clean up
    element_clear(left_9);
    element_clear(right_9);

    /*
    * Check 10: e(U_s / U_s_r, g) = e(g^{r-1}, Z_s_r)
    */
    element_t left_10, right_10, left_10_1, right_10_1;
    element_init_GT(left_10, pk.pairing);
    element_init_GT(right_10, pk.pairing);
    element_init_G1(left_10_1, pk.pairing);
    element_init_G1(right_10_1, pk.pairing);

    // left_10_1 = U_s / U_s_r
    element_div(left_10_1, res_acc.g_s, res_acc.g_s_r);
    // left_10 = e(U_s / U_s_r, g)
    element_pairing(left_10, left_10_1, pk.g);

    // right_10_1 = g^{r-1}
    element_div(right_10_1, pk.g_ri[0], pk.g);
    // right_10 = e(g^{r-1}, Z_s_r)
    element_pairing(right_10, right_10_1, Z_s_r);
    
    // Check if left_10 == right_10
    bool check_10 = (element_cmp(left_10, right_10) == 0);
    if (!check_10) {
        std::cout << "nested union check_10: " << check_10 << std::endl;
    }

    // Clean up
    element_clear(left_10);
    element_clear(right_10);
    element_clear(left_10_1);
    element_clear(right_10_1);
    
    /*
    * Check 11: e(U_r_s, g^gamma) = e(U_r_s_gamma, g)
    */
    element_t left_11, right_11;
    element_init_GT(left_11, pk.pairing);
    element_init_GT(right_11, pk.pairing);
    
    // e(U_r_s, g^gamma)
    element_pairing(left_11, res_acc.g_r_s, pk.g_gamma);
    // e(U_r_s_gamma, g)
    element_pairing(right_11, U_r_s_gamma, pk.g);
    
    // Check if left_11 == right_11
    bool check_11 = (element_cmp(left_11, right_11) == 0);
    if (!check_11) {
        std::cout << "nested union check_11: " << check_11 << std::endl;
    }

    // Clean up
    element_clear(left_11);
    element_clear(right_11);

    // All checks must pass
    return check_1 && check_2 && check_3 && check_4 && check_6 
    && check_8 && check_9 && check_10 && check_11;
}


/*
* NestedDifferenceProof ====================================================
*/

NestedDifferenceProof::NestedDifferenceProof(pairing_t pairing) {
    element_init_G1(I_r, pairing);
    element_init_G1(I_r_beta, pairing);
    element_init_G1(Q_s_r, pairing);
    element_init_G1(Q_s_r_delta, pairing);
    element_init_G1(L_r, pairing);

    element_init_G1(Z_s_r, pairing);
    element_init_G1(D_s_r_gamma, pairing);
    element_init_G1(Z_r_s, pairing);
    element_init_G1(D_r_s_gamma, pairing);
}

size_t NestedDifferenceProof::getSizeADS() {
    return (element_length_in_bytes(I_r) + element_length_in_bytes(I_r_beta) + element_length_in_bytes(Q_s_r) + 
    element_length_in_bytes(Q_s_r_delta) + element_length_in_bytes(L_r) + element_length_in_bytes(Z_s_r) + 
    element_length_in_bytes(D_s_r_gamma) + element_length_in_bytes(Z_r_s) + element_length_in_bytes(D_r_s_gamma)) * sizeof(char);
}

// Copy constructor
NestedDifferenceProof::NestedDifferenceProof(const NestedDifferenceProof& other) {
    element_init_same_as(I_r, const_cast<element_t&>(other.I_r));
    element_init_same_as(I_r_beta, const_cast<element_t&>(other.I_r_beta));
    element_init_same_as(Q_s_r, const_cast<element_t&>(other.Q_s_r));
    element_init_same_as(Q_s_r_delta, const_cast<element_t&>(other.Q_s_r_delta));
    element_init_same_as(L_r, const_cast<element_t&>(other.L_r));
    
    element_init_same_as(Z_s_r, const_cast<element_t&>(other.Z_s_r));
    element_init_same_as(D_s_r_gamma, const_cast<element_t&>(other.D_s_r_gamma));
    element_init_same_as(Z_r_s, const_cast<element_t&>(other.Z_r_s));
    element_init_same_as(D_r_s_gamma, const_cast<element_t&>(other.D_r_s_gamma));
    
    // Copy values
    element_set(I_r, const_cast<element_t&>(other.I_r));
    element_set(I_r_beta, const_cast<element_t&>(other.I_r_beta));
    element_set(Q_s_r, const_cast<element_t&>(other.Q_s_r));
    element_set(Q_s_r_delta, const_cast<element_t&>(other.Q_s_r_delta));
    element_set(L_r, const_cast<element_t&>(other.L_r));
    element_set(Z_s_r, const_cast<element_t&>(other.Z_s_r));
    element_set(D_s_r_gamma, const_cast<element_t&>(other.D_s_r_gamma));
    element_set(Z_r_s, const_cast<element_t&>(other.Z_r_s));
    element_set(D_r_s_gamma, const_cast<element_t&>(other.D_r_s_gamma));
}

// Assignment operator
NestedDifferenceProof& NestedDifferenceProof::operator=(const NestedDifferenceProof& other) {
    if (this != &other) {
        // Clear existing elements
        element_clear(I_r);
        element_clear(I_r_beta);
        element_clear(Q_s_r);
        element_clear(Q_s_r_delta);
        element_clear(L_r);
        element_clear(Z_s_r);
        element_clear(D_s_r_gamma);
        element_clear(Z_r_s);
        element_clear(D_r_s_gamma);
        
        // Initialize elements in the same field
        element_init_same_as(I_r, const_cast<element_t&>(other.I_r));
        element_init_same_as(I_r_beta, const_cast<element_t&>(other.I_r_beta));
        element_init_same_as(Q_s_r, const_cast<element_t&>(other.Q_s_r));
        element_init_same_as(Q_s_r_delta, const_cast<element_t&>(other.Q_s_r_delta));
        element_init_same_as(L_r, const_cast<element_t&>(other.L_r));
        element_init_same_as(Z_s_r, const_cast<element_t&>(other.Z_s_r));
        element_init_same_as(D_s_r_gamma, const_cast<element_t&>(other.D_s_r_gamma));
        element_init_same_as(Z_r_s, const_cast<element_t&>(other.Z_r_s));
        element_init_same_as(D_r_s_gamma, const_cast<element_t&>(other.D_r_s_gamma));
        
        // Copy values
        element_set(I_r, const_cast<element_t&>(other.I_r));
        element_set(I_r_beta, const_cast<element_t&>(other.I_r_beta));
        element_set(Q_s_r, const_cast<element_t&>(other.Q_s_r));
        element_set(Q_s_r_delta, const_cast<element_t&>(other.Q_s_r_delta));
        element_set(L_r, const_cast<element_t&>(other.L_r));
        element_set(Z_s_r, const_cast<element_t&>(other.Z_s_r));
        element_set(D_s_r_gamma, const_cast<element_t&>(other.D_s_r_gamma));
        element_set(Z_r_s, const_cast<element_t&>(other.Z_r_s));
        element_set(D_r_s_gamma, const_cast<element_t&>(other.D_r_s_gamma));
    }
    return *this;
}

NestedDifferenceProof::~NestedDifferenceProof() noexcept {
    element_clear(I_r);
    element_clear(I_r_beta);
    element_clear(Q_s_r);
    element_clear(Q_s_r_delta);
    element_clear(L_r);

    element_clear(Z_s_r);
    element_clear(D_s_r_gamma);
    element_clear(Z_r_s);
    element_clear(D_r_s_gamma);
}

bool NestedDifferenceProof::verify_nested(AccValue& lhs_acc, AccValue& rhs_acc, AccValue& res_acc, AccPublicKey& pk) {
    // First verify the intersection part of the proof
    /*
    * Check 1: e(A_s, B_{r,s}) = e(I_r, g^{s^q}) × e(Q_{s,r}, g)
    */
    element_t left_1, right_1, right_1_1, right_1_2;
    element_init_GT(left_1, pk.pairing);
    element_init_GT(right_1, pk.pairing);
    element_init_GT(right_1_1, pk.pairing);
    element_init_GT(right_1_2, pk.pairing);
    
    // e(A_s, B_{r,s})
    element_pairing(left_1, lhs_acc.g_s, rhs_acc.g_r_s);
    // e(I_r, g^{s^q})
    element_pairing(right_1_1, I_r, pk.g_sq);
    // e(Q_{s,r}, g)
    element_pairing(right_1_2, Q_s_r, pk.g);
    // Combine: right_1 = right_1_1 × right_1_2
    element_mul(right_1, right_1_1, right_1_2);
    // Check if left_1 == right_1
    bool check_1 = (element_cmp(left_1, right_1) == 0);
    if (!check_1) {
        std::cout << "nested difference check_1: " << check_1 << std::endl;
    }

    // Clean up
    element_clear(left_1);
    element_clear(right_1);
    element_clear(right_1_1);
    element_clear(right_1_2);
    
    /*
    * Check 2: e(I_r, g^beta) = e(I_{r,beta}, g)
    */
    element_t left_2, right_2;
    element_init_GT(left_2, pk.pairing);
    element_init_GT(right_2, pk.pairing);
    
    // e(I_r, g^beta)
    element_pairing(left_2, I_r, pk.g_beta);
    // e(I_{r,beta}, g)
    element_pairing(right_2, I_r_beta, pk.g);
    // Check if left_2 == right_2
    bool check_2 = (element_cmp(left_2, right_2) == 0);
    if (!check_2) {
        std::cout << "nested difference check_2: " << check_2 << std::endl;
    }

    // Clean up
    element_clear(left_2);
    element_clear(right_2);
    
    /*
    * Check 3: e(Q_{s,r}, g^delta) = e(Q_{s,r,delta}, g)
    */
    element_t left_3, right_3;
    element_init_GT(left_3, pk.pairing);
    element_init_GT(right_3, pk.pairing);
    
    // e(Q_{s,r}, g^delta)
    element_pairing(left_3, Q_s_r, pk.g_delta);
    // e(Q_{s,r,delta}, g)
    element_pairing(right_3, Q_s_r_delta, pk.g);
    // Check if left_3 == right_3
    bool check_3 = (element_cmp(left_3, right_3) == 0);
    if (!check_3) {
        std::cout << "nested difference check_3: " << check_3 << std::endl;
    }

    // Clean up
    element_clear(left_3);
    element_clear(right_3);
    
    /*
    * Check 4: e(I_r, g) = e(L_r, g^r)
    */
    element_t left_4, right_4;
    element_init_GT(left_4, pk.pairing);
    element_init_GT(right_4, pk.pairing);
    
    // e(I_r, g)
    element_pairing(left_4, I_r, pk.g);
    // e(L_r, g^r)
    element_pairing(right_4, L_r, pk.g_ri[0]);
    // Check if left_4 == right_4
    bool check_4 = (element_cmp(left_4, right_4) == 0);
    if (!check_4) {
        std::cout << "nested difference check_4: " << check_4 << std::endl;
    }

    // Clean up
    element_clear(left_4);
    element_clear(right_4);
    
    // /*
    // * Check 5: Verify I_r = g^{\sum_{i∈I} r^i} (I_r is indeed the accumulator for the intersection set)
    // */
    // element_t computed_I_r;
    // element_init_G1(computed_I_r, pk.pairing);
    // element_set1(computed_I_r);
    
    // // Compute g^{\sum_{i∈I} r^i}
    // for (const auto& i : intersection) {
    //     if (i < 1 || i >= pk.q) continue;  // Skip invalid elements
    //     element_mul(computed_I_r, computed_I_r, pk.g_ri[i - 1]);  // Multiply by g^{r^i}
    // }
    
    // // Check if computed_I_r == I_r
    // bool check_5 = (element_cmp(computed_I_r, I_r) == 0);
    // std::cout << "difference check_5: " << check_5 << std::endl;
    
    // // Clean up
    // element_clear(computed_I_r);
    
    /*
    * Check 6: Verify D_r = A_r / I_r (accumulator for difference = accumulator for A / accumulator for intersection)
    */
    element_t expected_D_r, I_r_inv;
    element_init_G1(expected_D_r, pk.pairing);
    element_init_G1(I_r_inv, pk.pairing);
    
    // Compute the inverse of I_r
    element_set(I_r_inv, I_r);
    element_invert(I_r_inv, I_r_inv);
    
    // expected_D_r = A_r / I_r = A_r * I_r^(-1)
    element_set(expected_D_r, lhs_acc.g_r);
    element_mul(expected_D_r, expected_D_r, I_r_inv);
    
    // Check if expected_D_r == D_r
    bool check_6 = (element_cmp(expected_D_r, res_acc.g_r) == 0);
    if (!check_6) {
        std::cout << "nested difference check_6: " << check_6 << std::endl;
    }
    
    // Clean up
    element_clear(expected_D_r);
    element_clear(I_r_inv);
    
    // /*
    // * Check 7: Verify D_r = g^{\sum_{i∈D} r^i} (D_r is indeed the accumulator for the difference set)
    // */
    // element_t computed_D_r;
    // element_init_G1(computed_D_r, pk.pairing);
    // element_set1(computed_D_r);
    
    // // Compute g^{\sum_{i∈D} r^i}
    // for (const auto& i : result) {
    //     if (i < 1 || i >= pk.q) continue;  // Skip invalid elements
    //     element_mul(computed_D_r, computed_D_r, pk.g_ri[i - 1]);  // Multiply by g^{r^i}
    // }
    
    // // Check if computed_D_r == D_r
    // bool check_7 = (element_cmp(computed_D_r, D_r) == 0);
    // std::cout << "difference check_7: " << check_7 << std::endl;
    
    // // Clean up
    // element_clear(computed_D_r);

    /*
    * Check 8: e(D_r / D_r_s, g) = e(g^{s-1}, Z_r_s)
    */
    element_t left_8, right_8, left_8_1, right_8_1;
    element_init_GT(left_8, pk.pairing);
    element_init_GT(right_8, pk.pairing);
    element_init_G1(left_8_1, pk.pairing);
    element_init_G1(right_8_1, pk.pairing);

    // left_8_1 = D_r / D_r_s
    element_div(left_8_1, res_acc.g_r, res_acc.g_r_s);
    // left_8 = e(D_r / D_r_s, g)
    element_pairing(left_8, left_8_1, pk.g);

    // right_5_1 = g^{s-1}
    element_div(right_8_1, pk.g_si[0], pk.g);
    // right_8 = e(g^{s-1}, Z_r_s)
    element_pairing(right_8, right_8_1, Z_r_s);
    
    // Check if left_8 == right_8
    bool check_8 = (element_cmp(left_8, right_8) == 0);
    if (!check_8) {
        std::cout << "nested difference check_8: " << check_8 << std::endl;
    }

    // Clean up
    element_clear(left_8);
    element_clear(right_8);
    element_clear(left_8_1);
    element_clear(right_8_1);
    
    /*
    * Check 9: e(D_s_r, g^gamma) = e(D_s_r_gamma, g)
    */
    element_t left_9, right_9;
    element_init_GT(left_9, pk.pairing);
    element_init_GT(right_9, pk.pairing);
    
    // e(D_s_r, g^gamma)
    element_pairing(left_9, res_acc.g_s_r, pk.g_gamma);
    // e(D_s_r_gamma, g)
    element_pairing(right_9, D_s_r_gamma, pk.g);
    
    // Check if left_9 == right_9
    bool check_9 = (element_cmp(left_9, right_9) == 0);
    if (!check_9) {
        std::cout << "nested difference check_9: " << check_9 << std::endl;
    }

    // Clean up
    element_clear(left_9);
    element_clear(right_9);

    /*
    * Check 10: e(D_s / D_s_r, g) = e(g^{r-1}, Z_s_r)
    */
    element_t left_10, right_10, left_10_1, right_10_1;
    element_init_GT(left_10, pk.pairing);
    element_init_GT(right_10, pk.pairing);
    element_init_G1(left_10_1, pk.pairing);
    element_init_G1(right_10_1, pk.pairing);

    // left_10_1 = D_s / D_s_r
    element_div(left_10_1, res_acc.g_s, res_acc.g_s_r);
    // left_10 = e(D_s / D_s_r, g)
    element_pairing(left_10, left_10_1, pk.g);

    // right_10_1 = g^{r-1}
    element_div(right_10_1, pk.g_ri[0], pk.g);
    // right_10 = e(g^{r-1}, Z_s_r)
    element_pairing(right_10, right_10_1, Z_s_r);
    
    // Check if left_10 == right_10
    bool check_10 = (element_cmp(left_10, right_10) == 0);
    if (!check_10) {
        std::cout << "nested difference check_10: " << check_10 << std::endl;
    }

    // Clean up
    element_clear(left_10);
    element_clear(right_10);
    element_clear(left_10_1);
    element_clear(right_10_1);
    
    /*
    * Check 11: e(D_r_s, g^gamma) = e(D_r_s_gamma, g)
    */
    element_t left_11, right_11;
    element_init_GT(left_11, pk.pairing);
    element_init_GT(right_11, pk.pairing);
    
    // e(D_s_r, g^gamma)
    element_pairing(left_11, res_acc.g_r_s, pk.g_gamma);
    // e(D_r_s_gamma, g)
    element_pairing(right_11, D_r_s_gamma, pk.g);
    
    // Check if left_11 == right_11
    bool check_11 = (element_cmp(left_11, right_11) == 0);
    if (!check_11) {
        std::cout << "nested difference check_11: " << check_11 << std::endl;
    }

    // Clean up
    element_clear(left_11);
    element_clear(right_11);

    // All checks must pass
    return check_1 && check_2 && check_3 && check_4 && check_6 
    && check_8 && check_9 && check_10 && check_11;
}

/*
* NestedRangeProof ====================================================
*/

NestedRangeProof::NestedRangeProof(pairing_t pairing) : d_B(pairing), d_D(pairing), max_B(0), min_D(0),
    pi_BC(pairing), pi_BD(pairing), pi_CD(pairing), pi_1(pairing), pi_2(pairing),
    pi_C_max(pairing), pi_C_min(pairing), max_C(0), min_C(0) {
    // Initialize elements
    element_init_G1(B_s, pairing);
    element_init_G1(B_s_alpha, pairing);
    element_init_G1(B_s_r, pairing);
    element_init_G1(B_s_r_gamma, pairing);
    element_init_G1(D_s, pairing);
    element_init_G1(D_s_alpha, pairing);
    element_init_G1(Z_s_r, pairing);

    element_init_G1(ZC_s_r, pairing);
    element_init_G1(ZC_r_s, pairing);
}

size_t NestedRangeProof::getSizeADS() {
    return (element_length_in_bytes(B_s) + element_length_in_bytes(B_s_alpha) + element_length_in_bytes(B_s_r) + 
    element_length_in_bytes(B_s_r_gamma) + element_length_in_bytes(D_s) + element_length_in_bytes(D_s_alpha) + 
    element_length_in_bytes(Z_s_r) + element_length_in_bytes(ZC_s_r) + element_length_in_bytes(ZC_r_s)) * sizeof(char);
}

// Copy constructor
NestedRangeProof::NestedRangeProof(const NestedRangeProof& other) : d_B(other.d_B), d_D(other.d_D), max_B(other.max_B), min_D(other.min_D),
    pi_BC(other.pi_BC), pi_BD(other.pi_BD), pi_CD(other.pi_CD), pi_1(other.pi_1), pi_2(other.pi_2),
    pi_C_max(other.pi_C_max), pi_C_min(other.pi_C_min), max_C(other.max_C), min_C(other.min_C) {
    element_init_same_as(B_s, const_cast<element_t&>(other.B_s));
    element_init_same_as(B_s_alpha, const_cast<element_t&>(other.B_s_alpha));
    element_init_same_as(B_s_r, const_cast<element_t&>(other.B_s_r));
    element_init_same_as(B_s_r_gamma, const_cast<element_t&>(other.B_s_r_gamma));
    element_init_same_as(D_s, const_cast<element_t&>(other.D_s));
    element_init_same_as(D_s_alpha, const_cast<element_t&>(other.D_s_alpha));
    element_init_same_as(Z_s_r, const_cast<element_t&>(other.Z_s_r));

    element_init_same_as(ZC_s_r, const_cast<element_t&>(other.ZC_s_r));
    element_init_same_as(ZC_r_s, const_cast<element_t&>(other.ZC_r_s));

    // Copy values
    element_set(B_s, const_cast<element_t&>(other.B_s));
    element_set(B_s_alpha, const_cast<element_t&>(other.B_s_alpha));
    element_set(B_s_r, const_cast<element_t&>(other.B_s_r));
    element_set(B_s_r_gamma, const_cast<element_t&>(other.B_s_r_gamma));
    element_set(D_s, const_cast<element_t&>(other.D_s));
    element_set(D_s_alpha, const_cast<element_t&>(other.D_s_alpha));
    element_set(Z_s_r, const_cast<element_t&>(other.Z_s_r));
    element_set(ZC_s_r, const_cast<element_t&>(other.ZC_s_r));
    element_set(ZC_r_s, const_cast<element_t&>(other.ZC_r_s));
}

// Assignment operator
NestedRangeProof& NestedRangeProof::operator=(const NestedRangeProof& other) {
    if (this != &other) {
        // Clear existing elements
        element_clear(B_s);
        element_clear(B_s_alpha);
        element_clear(B_s_r);
        element_clear(B_s_r_gamma);
        element_clear(D_s);
        element_clear(D_s_alpha);
        element_clear(Z_s_r);
        element_clear(ZC_s_r);
        element_clear(ZC_r_s);

        // Initialize elements in the same field
        element_init_same_as(B_s, const_cast<element_t&>(other.B_s));
        element_init_same_as(B_s_alpha, const_cast<element_t&>(other.B_s_alpha));
        element_init_same_as(B_s_r, const_cast<element_t&>(other.B_s_r));
        element_init_same_as(B_s_r_gamma, const_cast<element_t&>(other.B_s_r_gamma));
        element_init_same_as(D_s, const_cast<element_t&>(other.D_s));
        element_init_same_as(D_s_alpha, const_cast<element_t&>(other.D_s_alpha));
        element_init_same_as(Z_s_r, const_cast<element_t&>(other.Z_s_r));
        element_init_same_as(ZC_s_r, const_cast<element_t&>(other.ZC_s_r));
        element_init_same_as(ZC_r_s, const_cast<element_t&>(other.ZC_r_s));
        
        // Copy values
        element_set(B_s, const_cast<element_t&>(other.B_s));
        element_set(B_s_alpha, const_cast<element_t&>(other.B_s_alpha));
        element_set(B_s_r, const_cast<element_t&>(other.B_s_r));
        element_set(B_s_r_gamma, const_cast<element_t&>(other.B_s_r_gamma));
        element_set(D_s, const_cast<element_t&>(other.D_s));
        element_set(D_s_alpha, const_cast<element_t&>(other.D_s_alpha));
        element_set(Z_s_r, const_cast<element_t&>(other.Z_s_r));
        element_set(ZC_s_r, const_cast<element_t&>(other.ZC_s_r));
        element_set(ZC_r_s, const_cast<element_t&>(other.ZC_r_s));

        max_B = other.max_B;
        min_D = other.min_D;

        pi_BC = other.pi_BC;
        pi_BD = other.pi_BD;
        pi_CD = other.pi_CD;
        pi_1 = other.pi_1;
        pi_2 = other.pi_2;

        pi_C_max = other.pi_C_max;
        pi_C_min = other.pi_C_min;
        max_C = other.max_C;
        min_C = other.min_C;
    }
    return *this;
}

NestedRangeProof::~NestedRangeProof() noexcept {
    // Clear elements
    element_clear(B_s);
    element_clear(B_s_alpha);
    element_clear(B_s_r);
    element_clear(B_s_r_gamma);
    element_clear(D_s);
    element_clear(D_s_alpha);
    element_clear(Z_s_r);

    element_clear(ZC_s_r);
    element_clear(ZC_r_s);
}

bool NestedRangeProof::verify_nested_range(AccValue& acc, AccValue& res_acc, uint64_t l, uint64_t r, AccPublicKey& pk) {
    /*
    * Check 1.1: e(B_s, g^alpha) = e(B_s_alpha, g)
    */
    element_t left, right;
    element_init_GT(left, pk.pairing);
    element_init_GT(right, pk.pairing);

    // left: e(B_s, g^alpha)
    element_pairing(left, B_s, pk.g_alpha);
    // right: e(B_s_alpha, g)
    element_pairing(right, B_s_alpha, pk.g);

    bool check_1_1 = (element_cmp(left, right) == 0);
    if (!check_1_1) {
        std::cout << "nested range check_1_1: " << check_1_1 << std::endl;
    }
    
    /*
    * Check 1.2: e(B_s_r, g^gamma) = e(B_s_r_gamma, g)
    */
    // left: e(B_s_r, g^gamma)
    element_pairing(left, B_s_r, pk.g_gamma);
    // right: e(B_s_r_gamma, g)
    element_pairing(right, B_s_r_gamma, pk.g);

    bool check_1_2 = (element_cmp(left, right) == 0);
    if (!check_1_2) {
        std::cout << "nested range check_1_2: " << check_1_2 << std::endl;
    }

    /*
    * Check 1.3: e(D_s, g^alpha) = e(D_s_alpha, g)
    */
    // left: e(D_s, g^alpha)
    element_pairing(left, D_s, pk.g_alpha);
    // right: e(D_s_alpha, g)
    element_pairing(right, D_s_alpha, pk.g);

    bool check_1_3 = (element_cmp(left, right) == 0);
    if (!check_1_3) {
        std::cout << "nested range check_1_3: " << check_1_3 << std::endl;
    }

    /*
    * Check 2: e(B_s / B_s_r, g) = e(g^{r-1}, Z_s_r)
    */
    element_t left_tmp, right_tmp;
    element_init_G1(left_tmp, pk.pairing);
    element_init_G1(right_tmp, pk.pairing);

    // left_tmp: B_s / B_s_r
    element_div(left_tmp, B_s, B_s_r);
    // left: e(B_s / B_s_r, g)
    element_pairing(left, left_tmp, pk.g);
   
    // right_tmp: g^{r-1}
    element_div(right_tmp, pk.g_ri[0], pk.g);
    // right: e(g^{r-1}, Z_s_r)
    element_pairing(right, right_tmp, Z_s_r);

    bool check_2 = (element_cmp(left, right) == 0);
    if (!check_2) {
        std::cout << "nested range check_2: " << check_2 << std::endl;
    }

    /*
    * Check 3: verify pi_BC
    */
    // empty set
    Set empty_set;
    bool check_3 = pi_BC.verify(d_B, res_acc, empty_set, pk);
    if (!check_3) {
        std::cout << "nested range check_3: " << check_3 << std::endl;
    }

    /*
    * Check 4: verify pi_BD
    */
    bool check_4 = pi_BD.verify(d_B, d_D, empty_set, pk);
    if (!check_4) {
        std::cout << "nested range check_4: " << check_4 << std::endl;
    }

    /*
    * Check 5: verify pi_CD
    */
    bool check_5 = pi_CD.verify(res_acc, d_D, empty_set, pk);
    if (!check_5) {
        std::cout << "nested range check_5: " << check_5 << std::endl;
    }
    
    /*
    * Check 6: B_s * C_s * D_s = A_s
    */
    // left: B_s * C_s
    element_mul(left_tmp, B_s, res_acc.g_s);
    // left: B_s * C_s * D_s
    element_mul(left_tmp, left_tmp, d_D.g_s);
    bool check_6 = (element_cmp(left_tmp, acc.g_s) == 0);
    if (!check_6) {
        std::cout << "nested range check_6: " << check_6 << std::endl;
    }

    /*
    * Check 7: verify pi_1
    */
   bool check_7;
    if (!element_is1(pi_1.pi_max)) {
        check_7 = pi_1.verify_aggr(d_B, max_B, pk);
    } else {
        check_7 = max_B == 0;
    }
    if (!check_7) {
        std::cout << "nested range check_7: " << check_7 << std::endl;
    }
    
    /*
    * Check 8: verify pi_2
    */
    bool check_8;
    if (!element_is1(pi_2.pi_min)) {
        check_8 = pi_2.verify_aggr(d_D, min_D, pk);
    } else {
        check_8 = min_D == 0;
    }
    if (!check_8) {
        std::cout << "nested range check_8: " << check_8 << std::endl;
    }

    /*
    * Check 9: max_B < l, min_C >= l, max_C <= r, min_D > r
    */
    bool check_9 = (max_B == 0 || max_B < l) && (min_C == 0 || min_C >= l) && (max_C == 0 || max_C <= r) && (min_D == 0 || min_D > r);
    if (!check_9) {
        std::cout << "nested range check_9: " << check_9 << std::endl;
    }

    /*
    * Check 10: e(C_s / C_s_r, g) = e(g^{r-1}, ZC_s_r)
    */
    // left_tmp: C_s / C_s_r
    element_div(left_tmp, res_acc.g_s, res_acc.g_s_r);
    // left: e(C_s / C_s_r, g)
    element_pairing(left, left_tmp, pk.g);
   
    // right_tmp: g^{r-1}
    element_div(right_tmp, pk.g_ri[0], pk.g);
    // right: e(g^{r-1}, ZC_s_r)
    element_pairing(right, right_tmp, ZC_s_r);

    bool check_10 = (element_cmp(left, right) == 0);
    if (!check_10) {
        std::cout << "nested range check_10: " << check_10 << std::endl;
    }

    /*
    * Check 11: e(C_r / C_r_s, g) = e(g^{s-1}, ZC_r_s)
    */
    // left_tmp: C_r / C_r_s
    element_div(left_tmp, res_acc.g_r, res_acc.g_r_s);
    // left: e(C_r / C_r_s, g)
    element_pairing(left, left_tmp, pk.g);

    // right_tmp: g^{s-1}
    element_div(right_tmp, pk.g_si[0], pk.g);
    // right: e(g^{s-1}, ZC_r_s)
    element_pairing(right, right_tmp, ZC_r_s);

    bool check_11 = (element_cmp(left, right) == 0);
    if (!check_11) {
        std::cout << "nested range check_11: " << check_11 << std::endl;
    }

    /*
    * Check 12: pi_C_max
    */
    bool check_12;
    if (!element_is1(pi_C_max.pi_max)) {
        check_12 = pi_C_max.verify_aggr(res_acc, max_C, pk);
    } else {
        check_12 = max_C == 0;
    }
    if (!check_12) {
        std::cout << "nested range check_12: " << check_12 << std::endl;
    }

    /*
    * Check 13: pi_C_min
    */
    bool check_13;
    if (!element_is1(pi_C_min.pi_min)) {
        check_13 = pi_C_min.verify_aggr(res_acc, min_C, pk);
    } else {
        check_13 = min_C == 0;
    }
    if (!check_13) {
        std::cout << "nested range check_13: " << check_13 << std::endl;
    }

    element_clear(left);
    element_clear(right);
    element_clear(left_tmp);
    element_clear(right_tmp);

    return check_1_1 && check_1_2 && check_1_3 && check_2 && check_3 && check_4 && check_5 
    && check_6 && check_7 && check_8 && check_9 && check_10 && check_11 && check_12 && check_13;
}

} // namespace acc