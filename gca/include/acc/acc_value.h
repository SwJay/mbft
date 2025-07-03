/**
 * @file acc_value.h
 * @brief Accumulator value definition for the Expressive Set Accumulator
 */

#ifndef ACC_VALUE_H
#define ACC_VALUE_H

#include <pbc/pbc.h>
#include "acc/set.h"
#include "acc/keys.h"

namespace acc {

// Forward declarations
class AccPublicKey;

/**
 * @brief Represents an accumulator value for a set
 */
class AccValue {
public:
    element_t g_s;    // g^{\sum s^i}
    element_t g_r;    // g^{\sum r^i}
    element_t g_s_r;  // g^{\sum s^i \cdot r^{q-i}} (replacing h_s_r)
    element_t g_r_s;  // g^{\sum r^i \cdot s^{q-i}} (replacing h_r_s)

public:
    // Constructors and destructor
    // AccValue(); // 默认构造函数
    AccValue(pairing_t pairing);
    AccValue(const AccValue& other); // 拷贝构造函数
    AccValue& operator=(const AccValue& other); // 赋值运算符
    ~AccValue();

    uint32_t hash(const unsigned int nSeed=0) const;

    bool equals(const AccValue& other) const;

    size_t getSizeADS();

    // setup
    static AccValue setup(const Set& set, AccPublicKey& pk);
};

} // namespace acc

#endif // ACC_VALUE_H 