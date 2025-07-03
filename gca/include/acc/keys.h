/**
 * @file keys.h
 * @brief Secret and public key definitions for the Expressive Set Accumulator
 */

#ifndef ACC_KEYS_H
#define ACC_KEYS_H

#include <pbc/pbc.h>
#include <vector>  // Added for serialize method return type

namespace acc {

/**
 * @brief Secret key for the accumulator
 */
class AccSecretKey {
public:
    element_t s;       // Secret value s
    element_t r;       // Secret value r
    element_t alpha;   // Secret value alpha
    element_t beta;    // Secret value beta
    element_t gamma;   // Secret value gamma
    element_t delta;   // Secret value delta (added as per the documentation)

    char param[1024];
    size_t count;
    pairing_t pairing; // The pairing to use (added to enable proper copy/assignment)

public:
    // Constructors and destructor
    AccSecretKey(char* param, size_t count);
    AccSecretKey(const AccSecretKey& other);
    // AccSecretKey& operator=(const AccSecretKey& other) = delete; // Delete copy assignment to prevent issues
    ~AccSecretKey();
    
    // Generate a random secret key
    static AccSecretKey random(char* param, size_t count);
};

/**
 * @brief Public key for the accumulator
 */
class AccPublicKey {
public:
    uint64_t q;                    // Size of the universe
    element_t g;                   // Generator (only using one generator as per docs)

    element_t g_alpha;             // g^alpha
    element_t g_beta;              // g^beta
    element_t g_gamma;             // g^gamma
    element_t g_delta;             // g^delta

    element_t g_sq;               // g^{s^q}
    element_t g_rq;               // g^{r^q}
    
    element_t* g_si;  // g^{s^i} for i ∈ [q-1]
    element_t* g_alpha_si;  // g^{alpha·s^i} for i ∈ [q-1]
    element_t* g_ri;  // g^{r^i} for i ∈ [q-1]
    element_t* g_beta_ri;  // g^{beta·r^i} for i ∈ [q-1]
    element_t* g_gamma_ri_sqi;  // g^{gamma·r^i·s^{q-i}} for i ∈ [q-1]
    // We'll store r^i·s^j as a flattened 2D array
    element_t* g_ri_sj;  // g^{r^i·s^j} for (i,j) ∈ ([2q-1] \ {q}) × ([2q-1] \ {q})
    element_t* g_delta_ri_sj;  // g^{delta·r^i·s^j} for (i,j) ∈ ([2q-1] \ {q}) × ([2q-1] \ {q})
    
    char param[1024];
    size_t count;
    pairing_t pairing;  // The pairing to use

public:
    // Constructors and destructor
    // AccPublicKey() = default;
    AccPublicKey(char* param, size_t count, uint64_t q);
    AccPublicKey(const AccPublicKey& other);
    // AccPublicKey& operator=(const AccPublicKey& other) = delete; // Delete copy assignment to prevent issues
    ~AccPublicKey();
    
    // Generate public key from secret key
    static AccPublicKey generate(AccSecretKey& sk, char* param, size_t count, uint64_t q);
    static AccPublicKey generate_mul_thread(AccSecretKey& sk, char* param, size_t count, uint64_t q);
    static AccPublicKey generate_mul_thread(element_t g, AccSecretKey& sk, char* param, size_t count, uint64_t q);
    // helper function
    static size_t map_i_j_to_index(uint64_t i, uint64_t j, uint64_t q);
    
    // Serialization methods
    /**
     * @brief Serialize the public key to bytes
     * @return Vector of bytes containing the serialized public key
     */
    std::vector<unsigned char> serialize();
    
    /**
     * @brief Calculate the byte size needed to serialize this public key
     * @return Size in bytes
     */
    size_t serialized_size();
    
    /**
     * @brief Deserialize the public key from bytes
     * @param data Pointer to byte array containing serialized public key
     * @param size Size of the byte array
     * @return Number of bytes read from the array
     */
    size_t deserialize(const unsigned char* data, size_t size);
    size_t deserialize_mul_thread(const unsigned char* data, size_t size);
    
    /**
     * @brief Construct a public key by deserializing from bytes
     * @param data Pointer to byte array containing serialized public key
     * @param size Size of the byte array
     * @return Deserialized public key
     */
    static AccPublicKey from_bytes(const unsigned char* data, size_t size);
};

} // namespace acc

#endif // ACC_KEYS_H 