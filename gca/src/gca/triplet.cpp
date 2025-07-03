/**
 * @file triplet.cpp
 * @brief Implementation of the Triplet class
 */

#include <cstring>
#include <string>

#include "mht/hash.h"
#include "gca/triplet.h"

namespace gca {

Triplet::Triplet(uint32_t value, uint32_t nextValue, const acc::AccValue& digest)
    : value(value), nextValue(nextValue), digest(digest) {
}

size_t Triplet::getSizeADS() {
    return sizeof(value) + sizeof(nextValue) + digest.getSizeADS();
}

uint32_t Triplet::getValue() const {
    return value;
}

uint32_t Triplet::getNextValue() const {
    return nextValue;
}

const acc::AccValue& Triplet::getDigest() const {
    return digest;
}

// std::vector<uint8_t> Triplet::serialize() const {
//     // Serialize value and nextValue
//     std::vector<uint8_t> result;
    
//     // Reserve space for value, nextValue, and a digest placeholder
//     // In a real implementation, we would need to determine the actual size needed for the digest
//     const size_t digestPlaceholderSize = 4; // Just a placeholder size
//     result.resize(sizeof(uint32_t) * 2 + digestPlaceholderSize);
    
//     // Serialize value
//     std::memcpy(result.data(), &value, sizeof(uint32_t));
    
//     // Serialize nextValue
//     std::memcpy(result.data() + sizeof(uint32_t), &nextValue, sizeof(uint32_t));
    
//     // Add a placeholder for the digest
//     // In a real implementation, we would need to serialize the elements of the digest
//     // For now, we'll just add a placeholder to indicate that the digest is present
//     uint32_t digestPlaceholder = 0xDEADBEEF; // Magic number as a placeholder
//     std::memcpy(result.data() + sizeof(uint32_t) * 2, &digestPlaceholder, digestPlaceholderSize);
    
//     return result;
// }

// Triplet Triplet::deserialize(const std::vector<uint8_t>& data) {
//     // Check if we have enough data for value, nextValue, and the digest placeholder
//     const size_t digestPlaceholderSize = 4; // Just a placeholder size
//     if (data.size() < sizeof(uint32_t) * 2 + digestPlaceholderSize) {
//         // Not enough data
//         return Triplet(0, 0, acc::AccValue());
//     }
    
//     // Deserialize value
//     uint32_t value;
//     std::memcpy(&value, data.data(), sizeof(uint32_t));
    
//     // Deserialize nextValue
//     uint32_t nextValue;
//     std::memcpy(&nextValue, data.data() + sizeof(uint32_t), sizeof(uint32_t));
    
//     // Check the digest placeholder
//     uint32_t digestPlaceholder;
//     std::memcpy(&digestPlaceholder, data.data() + sizeof(uint32_t) * 2, digestPlaceholderSize);
    
//     // In a real implementation, we would deserialize the digest here
//     // For now, we'll just check if the placeholder matches our magic number
//     if (digestPlaceholder != 0xDEADBEEF) {
//         // Invalid digest placeholder
//         return Triplet(0, 0, acc::AccValue());
//     }
    
//     // Create a default digest
//     // In a real implementation, we would deserialize the digest from the data
//     return Triplet(value, nextValue, acc::AccValue());
// }

uint32_t Triplet::hash(const unsigned int nSeed) const {
    std::string v1 = std::to_string(value);
    std::string v2 = std::to_string(nextValue);
    std::string a_s = std::to_string(digest.hash(nSeed));
    std::string str = v1 + v2 + a_s;
    
    std::vector<unsigned char> vec(str.begin(), str.end());
	uint32_t hash;
    MurmurHash3_x86_32(vec.data(), vec.size(), nSeed, &hash);
    return hash;
}

} // namespace gca
