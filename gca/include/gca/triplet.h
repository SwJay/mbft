/**
 * @file triplet.h
 * @brief Triplet implementation for the GCA²-tree
 */

#ifndef GCA_TRIPLET_H
#define GCA_TRIPLET_H

#include <cstdint>
#include <vector>

#include "acc/acc_value.h"

namespace gca {

/**
 * @brief Represents a triplet in the GCA²-tree
 */
class Triplet {
public:
    /**
     * @brief Constructor for Triplet
     * @param value The j-th value of the attribute
     * @param nextValue The (j+1)-th value of the attribute
     * @param digest The digest of the prefix set
     */
    Triplet(uint32_t value, uint32_t nextValue, const acc::AccValue& digest);

    /**
     * @brief Get the j-th value of the attribute
     * @return The j-th value
     */
    uint32_t getValue() const;

    /**
     * @brief Get the (j+1)-th value of the attribute
     * @return The (j+1)-th value
     */
    uint32_t getNextValue() const;

    /**
     * @brief Get the digest of the prefix set
     * @return The digest
     */
    const acc::AccValue& getDigest() const;

    // /**
    //  * @brief Serialize the triplet to a byte array
    //  * @return The serialized triplet
    //  */
    // std::vector<uint8_t> serialize() const;

    // /**
    //  * @brief Deserialize a triplet from a byte array
    //  * @param data The serialized triplet
    //  * @return The deserialized triplet
    //  */
    // static Triplet deserialize(const std::vector<uint8_t>& data);

    /**
     * @brief Calculate the hash of the triplet
     * @return The hash of the triplet
     */
    uint32_t hash(const unsigned int nSeed=0) const;

    /**
     * @brief Get the size of the triplet
     * @return The size of the triplet
     */
    size_t getSizeADS();

private:
    uint32_t value;
    uint32_t nextValue;
    acc::AccValue digest;
};

} // namespace gca

#endif // GCA_TRIPLET_H
