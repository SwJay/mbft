/**
 * @file acc_value.cpp
 * @brief Implementation of accumulator value class for the Expressive Set Accumulator
 */

#include "acc/acc_value.h"
#include "mht/hash.h"
#include <iostream>

namespace acc {

// 默认构造函数
// AccValue::AccValue() {
//     // 不初始化任何内容，需要外部调用者确保在使用前正确初始化
//     // 这个构造函数仅用于支持作为成员变量时的默认初始化
//     // 注意：使用前必须通过pairing_t参数的构造函数或赋值运算符进行正确初始化
//     std::cout << "**********!!!!!!!!!!!!!!" << std::endl;
// }

// Constructor
AccValue::AccValue(pairing_t pairing) {
    // Initialize elements in G1
    element_init_G1(g_s, pairing);
    element_init_G1(g_r, pairing);
    element_init_G1(g_s_r, pairing);  // Changed from h_s_r (G2) to g_s_r (G1)
    element_init_G1(g_r_s, pairing);  // Changed from h_r_s (G2) to g_r_s (G1)
    
    // Set to identity elements
    element_set1(g_s);
    element_set1(g_r);
    element_set1(g_s_r);
    element_set1(g_r_s);
}

// Copy constructor
AccValue::AccValue(const AccValue& other) {
    // Initialize elements in the same field
    element_init_same_as(g_s, const_cast<element_t&>(other.g_s));
    element_init_same_as(g_r, const_cast<element_t&>(other.g_r));
    element_init_same_as(g_s_r, const_cast<element_t&>(other.g_s_r));  // Changed from h_s_r to g_s_r
    element_init_same_as(g_r_s, const_cast<element_t&>(other.g_r_s));  // Changed from h_r_s to g_r_s
    
    // Copy values
    element_set(g_s, const_cast<element_t&>(other.g_s));
    element_set(g_r, const_cast<element_t&>(other.g_r));
    element_set(g_s_r, const_cast<element_t&>(other.g_s_r));  // Changed from h_s_r to g_s_r
    element_set(g_r_s, const_cast<element_t&>(other.g_r_s));  // Changed from h_r_s to g_r_s
}

// Assignment operator
AccValue& AccValue::operator=(const AccValue& other) {
    if (this != &other) {
        // Check if we need to re-initialize
        element_clear(g_s);
        element_clear(g_r);
        element_clear(g_s_r);  // Changed from h_s_r to g_s_r
        element_clear(g_r_s);  // Changed from h_r_s to g_r_s
        
        // Initialize elements in the same field
        element_init_same_as(g_s, const_cast<element_t&>(other.g_s));
        element_init_same_as(g_r, const_cast<element_t&>(other.g_r));
        element_init_same_as(g_s_r, const_cast<element_t&>(other.g_s_r));  // Changed from h_s_r to g_s_r
        element_init_same_as(g_r_s, const_cast<element_t&>(other.g_r_s));  // Changed from h_r_s to g_r_s
        
        // Copy values
        element_set(g_s, const_cast<element_t&>(other.g_s));
        element_set(g_r, const_cast<element_t&>(other.g_r));
        element_set(g_s_r, const_cast<element_t&>(other.g_s_r));  // Changed from h_s_r to g_s_r
        element_set(g_r_s, const_cast<element_t&>(other.g_r_s));  // Changed from h_r_s to g_r_s
    }
    return *this;
}

// Destructor
AccValue::~AccValue() {
    // Clear elements
    element_clear(g_s);
    element_clear(g_r);
    element_clear(g_s_r);  // Changed from h_s_r to g_s_r
    element_clear(g_r_s);  // Changed from h_r_s to g_r_s
}

AccValue AccValue::setup(const Set& set, AccPublicKey& pk) {
    // Create a new accumulator value
    AccValue acc(pk.pairing);
    
    // For each element in the set, multiply the accumulator by the corresponding element in the public key
    for (const auto& element : set) {
        if (element < 1 || element >= pk.q) {
            // std::cout << "Element " << element << " is out of range for public key " << pk.q << std::endl;
            continue;
        }
        
        // g^{s^i}
        element_mul(acc.g_s, acc.g_s, pk.g_si[element - 1]);
        
        // g^{r^i}
        element_mul(acc.g_r, acc.g_r, pk.g_ri[element - 1]);
        
        // g^{s^i \cdot r^{q-i}}
        size_t idx = AccPublicKey::map_i_j_to_index(pk.q - element, element, pk.q);
        element_mul(acc.g_s_r, acc.g_s_r, pk.g_ri_sj[idx]);
        
        // g^{r^i \cdot s^{q-i}}
        idx = AccPublicKey::map_i_j_to_index(element, pk.q - element, pk.q);
        element_mul(acc.g_r_s, acc.g_r_s, pk.g_ri_sj[idx]);
    }
    
    return acc;
}

uint32_t AccValue::hash(const unsigned int nSeed) const {
    uint32_t hash;

    int size1 = element_length_in_bytes(const_cast<element_t&>(g_s));
    int size2 = element_length_in_bytes(const_cast<element_t&>(g_r));
    int size3 = element_length_in_bytes(const_cast<element_t&>(g_s_r));
    int size4 = element_length_in_bytes(const_cast<element_t&>(g_r_s));

    size_t total_size = size1 + size2 + size3 + size4;
    
    // Create a single vector with the total size
    std::vector<unsigned char> result(total_size);
    
    // Keep track of current position in the vector
    size_t current_pos = 0;
    
    // Fill the vector with bytes from each element
    element_to_bytes(result.data() + current_pos, const_cast<element_t&>(g_s));
    current_pos += size1;
    
    element_to_bytes(result.data() + current_pos, const_cast<element_t&>(g_r));
    current_pos += size2;
    
    element_to_bytes(result.data() + current_pos, const_cast<element_t&>(g_s_r));
    current_pos += size3;
    
    element_to_bytes(result.data() + current_pos, const_cast<element_t&>(g_r_s));

    MurmurHash3_x86_32(result.data(), result.size(), nSeed, &hash);
    
    return hash;
}

bool AccValue::equals(const AccValue& other) const {
    return element_cmp(const_cast<element_t&>(g_s), const_cast<element_t&>(other.g_s)) == 0 &&
           element_cmp(const_cast<element_t&>(g_r), const_cast<element_t&>(other.g_r)) == 0 &&
           element_cmp(const_cast<element_t&>(g_s_r), const_cast<element_t&>(other.g_s_r)) == 0 &&
           element_cmp(const_cast<element_t&>(g_r_s), const_cast<element_t&>(other.g_r_s)) == 0;
}

size_t AccValue::getSizeADS() {
    size_t sizeADS = element_length_in_bytes(g_s) + element_length_in_bytes(g_r) + element_length_in_bytes(g_s_r) + element_length_in_bytes(g_r_s);
    return sizeADS * sizeof(char);
}

} // namespace acc 