/**
 * @file accumulator.cpp
 * @brief Implementation of the main Accumulator class methods
 */

#include <stdexcept>
#include <cstring>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <cmath>
#include <vector>
#include <iterator>

#include "acc/accumulator.h"

namespace acc {

// Key generation
// std::pair<AccSecretKey, AccPublicKey> Accumulator::genkey(uint64_t universe_size) {
//     FILE *file = fopen("../param/a.param","r");
//     char param[1024];
//     size_t count;

//     count = fread(param, 1, 1024, file);
//     // std::cout << "count: " << count << std::endl;
//     if(!count)
//         pbc_die("input error for setup");


//     AccSecretKey sk = AccSecretKey::random(param, count);
//     AccPublicKey pk = AccPublicKey::generate(sk, param, count, universe_size);
    
//     // return std::pair<AccSecretKey, AccPublicKey>(std::move(sk), std::move(pk));
//     return std::pair<AccSecretKey, AccPublicKey>(sk, pk);
// }

AccPublicKey Accumulator::genkey(std::string key_dir, std::string pbc_param_path, uint64_t universe_size) {
    if (!std::filesystem::exists(key_dir)) {
        std::filesystem::create_directories(key_dir);
    }
    // 1. determine if the public key from the file exists
    uint64_t q_bit = std::ceil(std::log2(universe_size));
    std::string key_path = key_dir + "/pk-" + std::to_string(q_bit);
    std::cout << "key_path: " << key_path << std::endl;

    std::ifstream inkey_file(key_path, std::ios::binary);
    if (inkey_file.good()) {
        // 如果文件存在，则从文件读取
        std::cout << "file exists" << std::endl;
        std::vector<unsigned char> key_data;
        if (inkey_file) {
            // Get the current position (likely the beginning, but good practice)
            auto start_pos = inkey_file.tellg();
            std::cout << "start_pos: " << start_pos << std::endl;

            // Seek to the end to get the size
            inkey_file.seekg(0, std::ios::end);
            auto file_size = inkey_file.tellg() - start_pos; // Calculate size relative to start
            std::cout << "file_size: " << file_size << std::endl;

            if (file_size < 0) {
                 // Handle error: unable to determine size correctly
                 inkey_file.close(); // Close the file before throwing
                 throw std::runtime_error("Could not determine size of the key file.");
            }


            // Seek back to the original position
            inkey_file.seekg(start_pos);

            if (file_size > 0) {
                // Resize the vector to fit the file content
                key_data.resize(static_cast<size_t>(file_size));

                // Read the file content directly into the vector
                auto start_time = std::chrono::high_resolution_clock::now();
                inkey_file.read(reinterpret_cast<char*>(key_data.data()), file_size);
                auto end_time = std::chrono::high_resolution_clock::now();
                auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
                std::cout << "time for reading the file content: " << duration.count() / 1000.0 << " seconds" << std::endl;

                // Check if the read was successful
                if (!inkey_file) {
                    inkey_file.close(); // Close the file before throwing
                    throw std::runtime_error("Error reading key file. Read " + std::to_string(inkey_file.gcount()) + " bytes out of " + std::to_string(file_size));
                }
            }
            // If file_size is 0, key_data remains empty, which is correct.
        } else {
             // Handle error: input file stream was invalid to begin with
             throw std::runtime_error("Input key file stream is invalid.");
        }
        return AccPublicKey::from_bytes(key_data.data(), key_data.size());
    }

    FILE *file = fopen(pbc_param_path.c_str(), "r");
    if (!file) {
        std::string error = "Cannot open parameter file: " + pbc_param_path;
        pbc_die(error.c_str());
    }
    
    char param[1024];
    size_t count;

    count = fread(param, 1, 1024, file);
    if (!count) {
        fclose(file);
        pbc_die("Input error for setup: couldn't read parameter file");
    }

    // 成功读取参数后关闭文件
    fclose(file);

    AccSecretKey sk = AccSecretKey::random(param, count);
    std::cout << "sk generated" << std::endl;
    AccPublicKey pk = AccPublicKey::generate_mul_thread(sk, param, count, universe_size);

    // 2. save the public key to the file
    std::ofstream outkey_file(key_path, std::ios::binary);
    if (!outkey_file.good()) {
        pbc_die("Cannot open key file for writing");
    }
    outkey_file.write(reinterpret_cast<const char*>(pk.serialize().data()), pk.serialize().size());
    outkey_file.close();

    return pk;
}

std::pair<AccPublicKey, AccPublicKey> Accumulator::genkey_test(std::string key_dir, std::string pbc_param_path, uint64_t universe_size) {
    // 1. determine if the public key from the file exists
    uint64_t q_bit = std::ceil(std::log2(universe_size));
    std::string key_path1 = key_dir + "/pk-" + std::to_string(q_bit) + "-1";
    std::string key_path2 = key_dir + "/pk-" + std::to_string(q_bit) + "-2";
    std::cout << "key_path1: " << key_path1 << std::endl;
    std::cout << "key_path2: " << key_path2 << std::endl;

    std::ifstream inkey_file1(key_path1, std::ios::binary);
    std::ifstream inkey_file2(key_path2, std::ios::binary);
    if (inkey_file1.good() && inkey_file2.good()) {
        // 如果文件存在，则从文件读取
        std::vector<unsigned char> key_data1(std::istreambuf_iterator<char>(inkey_file1), {});
        std::vector<unsigned char> key_data2(std::istreambuf_iterator<char>(inkey_file2), {});
        return std::make_pair(AccPublicKey::from_bytes(key_data1.data(), key_data1.size()), AccPublicKey::from_bytes(key_data2.data(), key_data2.size()));
    }

    FILE *file = fopen(pbc_param_path.c_str(), "r");
    if (!file) {
        std::string error = "Cannot open parameter file: " + pbc_param_path;
        pbc_die(error.c_str());
    }
    
    char param[1024];
    size_t count;

    count = fread(param, 1, 1024, file);
    if (!count) {
        fclose(file);
        pbc_die("Input error for setup: couldn't read parameter file");
    }

    // 成功读取参数后关闭文件
    fclose(file);

    AccSecretKey sk = AccSecretKey::random(param, count);
    AccPublicKey pk1 = AccPublicKey::generate(sk, param, count, universe_size);
    AccPublicKey pk2 = AccPublicKey::generate_mul_thread(pk1.g, sk, param, count, universe_size);

    // 2. save the public key to the file
    std::ofstream outkey_file1(key_path1, std::ios::binary);
    std::ofstream outkey_file2(key_path2, std::ios::binary);
    if (!outkey_file1.good() || !outkey_file2.good()) {
        pbc_die("Cannot open key file for writing");
    }
    outkey_file1.write(reinterpret_cast<const char*>(pk1.serialize().data()), pk1.serialize().size());
    outkey_file2.write(reinterpret_cast<const char*>(pk2.serialize().data()), pk2.serialize().size());
    outkey_file1.close();
    outkey_file2.close();

    return std::make_pair(pk1, pk2);
}


// Digest setup
AccValue Accumulator::setup(const Set& set, AccPublicKey& pk) {
    // Create a new accumulator value
    return AccValue::setup(set, pk);
}

/*
// Update accumulator
AccValue Accumulator::update(const AccValue& acc, bool is_insert, uint64_t element, const AccPublicKey& pk) {
    // Create a new accumulator value
    AccValue result(pk.get_pairing());
    
    // Copy the existing accumulator
    result = acc;
    
    // If element is out of range, return the original accumulator
    if (element < 1 || element >= pk.get_q()) {
        return result;
    }
    
    // Update the accumulator based on insert or delete operation
    if (is_insert) {
        // Multiply with the appropriate elements from the public key
        element_mul(result.g_s, result.g_s, pk.get_g_s_i(element - 1));
        element_mul(result.g_r, result.g_r, pk.get_g_r_i(element - 1));
        element_mul(result.g_s_r, result.g_s_r, pk.get_g_s_r_i(element - 1));
        element_mul(result.g_r_s, result.g_r_s, pk.get_g_r_s_i(element - 1));
    } else {
        // For deletion, compute the inverse of the elements
        element_t g_s_inv, g_r_inv, g_s_r_inv, g_r_s_inv;
        
        element_init_G1(g_s_inv, pk.get_pairing());
        element_init_G1(g_r_inv, pk.get_pairing());
        element_init_G1(g_s_r_inv, pk.get_pairing());
        element_init_G1(g_r_s_inv, pk.get_pairing());
        
        // Get the elements from the public key
        element_set(g_s_inv, pk.get_g_s_i(element - 1));
        element_set(g_r_inv, pk.get_g_r_i(element - 1));
        element_set(g_s_r_inv, pk.get_g_s_r_i(element - 1));
        element_set(g_r_s_inv, pk.get_g_r_s_i(element - 1));
        
        // Compute the inverse of each element
        element_invert(g_s_inv, g_s_inv);
        element_invert(g_r_inv, g_r_inv);
        element_invert(g_s_r_inv, g_s_r_inv);
        element_invert(g_r_s_inv, g_r_s_inv);
        
        // Multiply the accumulator by the inverse
        element_mul(result.g_s, result.g_s, g_s_inv);
        element_mul(result.g_r, result.g_r, g_r_inv);
        element_mul(result.g_s_r, result.g_s_r, g_s_r_inv);
        element_mul(result.g_r_s, result.g_r_s, g_r_s_inv);
        
        // Clean up
        element_clear(g_s_inv);
        element_clear(g_r_inv);
        element_clear(g_s_r_inv);
        element_clear(g_r_s_inv);
    }
    
    return result;
}

*/

// Compute the intersection and proof
std::pair<Set, IntersectionProof> Accumulator::query_intersection(const Set& a, const Set& b, AccPublicKey& pk) {
    // Calculate the intersection
    Set intersection = a.intersection(b);
    
    // Create a new intersection proof
    IntersectionProof proof(pk.pairing);
    
    // Pre-compute the accumulator values for the intersection
    AccValue i_acc = setup(intersection, pk);
    
    // 1. Calculate I_r =  g^I(r)
    element_set(proof.I_r, i_acc.g_r);
    
    // 2. Calculate I_r_beta =  g^{beta * I(r)}
    element_t temp;
    element_init_G1(temp, pk.pairing);
    element_set1(temp);
    
    for (const auto& i : intersection) {
        if (i < 1 || i >= pk.q) continue;
        element_mul(temp, temp, pk.g_beta_ri[i - 1]);
    }
    element_set(proof.I_r_beta, temp);
    
    // 3. Calculate Q_s_r =  g^q(s, r)
    element_set1(temp);
    size_t idx;
    for (const auto& i : a) {
        for (const auto& j : b) {
            if (i != j) {
                idx = AccPublicKey::map_i_j_to_index(j, pk.q + i - j, pk.q);
                element_mul(temp, temp, pk.g_ri_sj[idx]);
            }
        }
    }
    element_set(proof.Q_s_r, temp);
    
    // 4. Calculate Q_s_r_delta (g^{delta·q(s,r)})
    element_set1(temp);
    for (const auto& i : a) {
        for (const auto& j : b) {
            if (i != j) {
                idx = AccPublicKey::map_i_j_to_index(j, pk.q + i - j, pk.q);
                element_mul(temp, temp, pk.g_delta_ri_sj[idx]);
            }
        }
    }
    element_set(proof.Q_s_r_delta, temp);
    
    // 5. Calculate L_r = g^{I(r) / r}
    element_set1(temp);
    for (const auto& i : intersection) {
        if (i < 1 || i >= pk.q) continue;
        if (i == 1) {
            element_mul(temp, temp, pk.g);
        } else {
            element_mul(temp, temp, pk.g_ri[i - 2]);  // g^{r^(i-1)}
        }
    }
    element_set(proof.L_r, temp);
    
    // Clean up
    element_clear(temp);
    
    return std::make_pair(intersection, proof);
}

// Compute the Union and proof
std::pair<Set, UnionProof> Accumulator::query_union(const Set& a, const Set& b, AccPublicKey& pk) {
    // Calculate the union U = A ∪ B
    Set union_set = a.union_with(b);
    
    // Calculate the intersection I = A ∩ B
    Set intersection = a.intersection(b);
    
    // Create a new union proof
    UnionProof proof(pk.pairing);
    
    // Pre-compute the accumulator values
    AccValue a_acc = setup(a, pk);
    AccValue i_acc = setup(intersection, pk);
    AccValue u_acc = setup(union_set, pk);
    
    // 1. Calculate I_r = g^I(r) (accumulator for intersection)
    element_set(proof.I_r, i_acc.g_r);
    
    // 2. Calculate I_r_beta = g^{beta * I(r)}
    element_t temp;
    element_init_G1(temp, pk.pairing);
    element_set1(temp);
    
    for (const auto& i : intersection) {
        if (i < 1 || i >= pk.q) continue;
        element_mul(temp, temp, pk.g_beta_ri[i - 1]);
    }
    element_set(proof.I_r_beta, temp);
    
    // 3. Calculate Q_s_r = g^q(s, r)
    element_set1(temp);
    size_t idx;
    for (const auto& i : a) {
        for (const auto& j : b) {
            if (i != j) {
                idx = AccPublicKey::map_i_j_to_index(j, pk.q + i - j, pk.q);
                element_mul(temp, temp, pk.g_ri_sj[idx]);
            }
        }
    }
    element_set(proof.Q_s_r, temp);
    
    // 4. Calculate Q_s_r_delta (g^{delta·q(s,r)})
    element_set1(temp);
    for (const auto& i : a) {
        for (const auto& j : b) {
            if (i != j) {
                idx = AccPublicKey::map_i_j_to_index(j, pk.q + i - j, pk.q);
                element_mul(temp, temp, pk.g_delta_ri_sj[idx]);
            }
        }
    }
    element_set(proof.Q_s_r_delta, temp);
    
    // 5. Calculate L_r = g^{I(r) / r}
    element_set1(temp);
    for (const auto& i : intersection) {
        if (i < 1 || i >= pk.q) continue;
        if (i == 1) {
            element_mul(temp, temp, pk.g);
        } else {
            element_mul(temp, temp, pk.g_ri[i - 2]);  // g^{r^(i-1)}
        }
    }
    element_set(proof.L_r, temp);
    
    // 6. Set u_r (accumulator for union set)
    element_set(proof.U_r, u_acc.g_r);
    
    // Clean up
    element_clear(temp);
    
    return std::make_pair(union_set, proof);
}
    

// Compute the difference and proof
std::pair<Set, DifferenceProof> Accumulator::query_difference(const Set& a, const Set& b, AccPublicKey& pk) {
    // Calculate the difference D = A \ B
    Set difference = a.difference(b);
    
    // Calculate the intersection I = A ∩ B
    Set intersection = a.intersection(b);
    
    // Create a new difference proof
    DifferenceProof proof(pk.pairing);
    
    // Pre-compute the accumulator values
    AccValue a_acc = setup(a, pk);
    AccValue i_acc = setup(intersection, pk);
    AccValue d_acc = setup(difference, pk);
    
    // 1. Calculate I_r = g^I(r) (accumulator for intersection)
    element_set(proof.I_r, i_acc.g_r);
    
    // 2. Calculate I_r_beta = g^{beta * I(r)}
    element_t temp;
    element_init_G1(temp, pk.pairing);
    element_set1(temp);
    
    for (const auto& i : intersection) {
        if (i < 1 || i >= pk.q) continue;
        element_mul(temp, temp, pk.g_beta_ri[i - 1]);
    }
    element_set(proof.I_r_beta, temp);
    
    // 3. Calculate Q_s_r = g^q(s, r)
    element_set1(temp);
    size_t idx;
    for (const auto& i : a) {
        for (const auto& j : b) {
            if (i != j) {
                idx = AccPublicKey::map_i_j_to_index(j, pk.q + i - j, pk.q);
                element_mul(temp, temp, pk.g_ri_sj[idx]);
            }
        }
    }
    element_set(proof.Q_s_r, temp);
    
    // 4. Calculate Q_s_r_delta (g^{delta·q(s,r)})
    element_set1(temp);
    for (const auto& i : a) {
        for (const auto& j : b) {
            if (i != j) {
                idx = AccPublicKey::map_i_j_to_index(j, pk.q + i - j, pk.q);
                element_mul(temp, temp, pk.g_delta_ri_sj[idx]);
            }
        }
    }
    element_set(proof.Q_s_r_delta, temp);
    
    // 5. Calculate L_r = g^{I(r) / r}
    element_set1(temp);
    for (const auto& i : intersection) {
        if (i < 1 || i >= pk.q) continue;
        if (i == 1) {
            element_mul(temp, temp, pk.g);
        } else {
            element_mul(temp, temp, pk.g_ri[i - 2]);  // g^{r^(i-1)}
        }
    }
    element_set(proof.L_r, temp);
    
    // 6. Set D_r (accumulator for difference set)
    element_set(proof.D_r, d_acc.g_r);
    
    // Clean up
    element_clear(temp);
    
    return std::make_pair(difference, proof);
}


// Compute count with proof
std::pair<uint64_t, CountProof> Accumulator::query_count(const Set& a, AccPublicKey& pk) {
    // Count is just the number of elements in the set
    uint64_t count = a.size();
    
    // Create a new count proof
    CountProof proof(pk.pairing);
    
    // Calculate g_a(s) = g^a(s) (accumulator for the set)
    element_set1(proof.g_a_s);
    for (const auto& i : a) {
        for (uint64_t j = 0; j < i; j++) {
            if (j == 0) {
                element_mul(proof.g_a_s, proof.g_a_s, pk.g);
            } else {
                element_mul(proof.g_a_s, proof.g_a_s, pk.g_si[j - 1]);
            }
        }
    }
    
    return std::make_pair(count, proof);
}

// Compute sum with proof
std::pair<uint64_t, SumProof> Accumulator::query_sum(const Set& a, AccPublicKey& pk) {
    // Calculate the sum of elements in the set
    uint64_t sum = 0;
    for (const auto& element : a) {
        sum += element;
    }
    
    // Create a new sum proof
    SumProof proof(pk.pairing);

    // Count is just the number of elements in the set
    proof.count = a.size();
    
    // Calculate g_a(s) = g^a(s) (accumulator for the set)
    element_set1(proof.g_b_s);
    for (const auto& i : a) {
        for (uint64_t j = 0; j < i; j++) {
            if (j == 0) continue;
            for (uint64_t k = 0; k < j; k++) {
                if (k == 0) {
                    element_mul(proof.g_b_s, proof.g_b_s, pk.g);
                } else {
                    element_mul(proof.g_b_s, proof.g_b_s, pk.g_si[k - 1]);
                }
            }
        }
    }
    
    return std::make_pair(sum, proof);
}

// Compute minimum with proof
std::pair<uint64_t, MinProof> Accumulator::query_min(const Set& a, AccPublicKey& pk) {
    // If set is empty, return a default value and null proof
    if (a.empty()) {
        // 因为不能返回nullptr，创建一个默认的证明对象
        MinProof empty_proof(pk.pairing);
        return std::make_pair(0, empty_proof);
    }
    
    // Find the minimum element in the set
    uint64_t min = UINT64_MAX;
    for (const auto& element : a) {
        if (element < min) {
            min = element;
        }
    }
    
    // Create a new minimum proof
    MinProof proof(pk.pairing);
    
    // Calculate g^{(A(s) - s^min)/s^{min+1}}
    element_set1(proof.pi_min);
    for (const auto& i : a) {
        if (i == min) continue;
        if (i == min + 1) {
            element_mul(proof.pi_min, proof.pi_min, pk.g);
        } else {
            element_mul(proof.pi_min, proof.pi_min, pk.g_si[i - min - 2]);
        }
    }

    return std::make_pair(min, proof);
}

// Compute maximum with proof
std::pair<uint64_t, MaxProof> Accumulator::query_max(const Set& a, AccPublicKey& pk) {
    // If set is empty, return a default value and null proof
    if (a.empty()) {
        // 因为不能返回nullptr，创建一个默认的证明对象
        MaxProof empty_proof(pk.pairing);
        return std::make_pair(0, empty_proof);
    }
    
    // Find the maximum element in the set
    uint64_t max = 0;
    for (const auto& element : a) {
        if (element > max) {
            max = element;
        }
    }
    
    // Create a new maximum proof
    MaxProof proof(pk.pairing);
    
    // Calculate g^{(A(r,s) - r^max*s^{q-max})/s^{q-max+1}}
    element_set1(proof.pi_max);
    for (const auto& i : a) {
        if (i == max) continue;
        if (i == max - 1) {
            element_mul(proof.pi_max, proof.pi_max, pk.g_ri[i - 1]);
        } else {
            size_t idx = AccPublicKey::map_i_j_to_index(i, max - i - 1, pk.q);
            element_mul(proof.pi_max, proof.pi_max, pk.g_ri_sj[idx]);
        }
    }

    return std::make_pair(max, proof);
}

// Compute range query with proof
std::pair<Set, RangeProof> Accumulator::query_range(const Set& a, uint64_t l, uint64_t r, AccPublicKey& pk) {
    // Filter the set to include only elements in the range [l, r]
    Set B, C, D;
    for (const auto& element : a) {
        if (element >= l && element <= r) {
            C.insert(element);
        } else if (element < l) {
            B.insert(element);
        } else {
            D.insert(element);
        }
    }
    
    // We reuse the intersection proof structure for simplicity
    // In a real implementation, this would use a specific RANGE proof
    RangeProof proof(pk.pairing);

    AccValue B_acc = setup(B, pk);
    AccValue C_acc = setup(C, pk);
    AccValue D_acc = setup(D, pk);

    proof.d_B = B_acc;
    proof.d_D = D_acc;

    // 1.1 B_s, B_s_alpha
    // Calculate B_s = g^B(s)
    element_set(proof.B_s, B_acc.g_s);
    // Calculate B_s_alpha = g^{alpha·B(s)}
    element_set1(proof.B_s_alpha);
    for (const auto& i : B) {
        element_mul(proof.B_s_alpha, proof.B_s_alpha, pk.g_alpha_si[i - 1]);
    }

    // 1.2 B_s_r, B_s_r_gamma
    // Calculate B_s_r = g^B(s,r)
    element_set(proof.B_s_r, B_acc.g_s_r);
    // Calculate B_s_r_gamma = g^{gamma·B(s,r)}
    element_set1(proof.B_s_r_gamma);
    for (const auto& i : B) {
        element_mul(proof.B_s_r_gamma, proof.B_s_r_gamma, pk.g_gamma_ri_sqi[pk.q - i - 1]);
    }

    // 1.3 D_s, D_s_alpha
    // Calculate D_s = g^D(s)
    element_set(proof.D_s, D_acc.g_s);
    // Calculate D_s_alpha = g^{alpha·D(s)}
    element_set1(proof.D_s_alpha);
    for (const auto& i : D) {
        element_mul(proof.D_s_alpha, proof.D_s_alpha, pk.g_alpha_si[i - 1]);
    }

    // 2. Z_s_r
    // Calculate Z_s_r = g^Z(s,r)
    element_set1(proof.Z_s_r);
    for (const auto& i : B) {
        for (uint64_t j = 0; j < pk.q - i; j++) {
            if (j == 0) {
                element_div(proof.Z_s_r, proof.Z_s_r, pk.g_si[i - 1]);
            } else {
                size_t idx = AccPublicKey::map_i_j_to_index(j, i, pk.q);
                element_div(proof.Z_s_r, proof.Z_s_r, pk.g_ri_sj[idx]);
            }
        }
    }

    // 3. pi_BC
    auto [BC, proof_BC] = query_intersection(B, C, pk);
    proof.pi_BC = IntersectionProof(proof_BC);

    // 4. pi_BD
    auto [BD, proof_BD] = query_intersection(B, D, pk);
    proof.pi_BD = IntersectionProof(proof_BD);

    // 5. pi_CD
    auto [CD, proof_CD] = query_intersection(C, D, pk);
    proof.pi_CD = IntersectionProof(proof_CD);

    // 6. pi_1
    auto [max, proof_max] = query_max(B, pk);
    proof.pi_1 = MaxProof(proof_max);
    proof.max_B = max;

    // 7. pi_2
    auto [min, proof_min] = query_min(D, pk);
    proof.pi_2 = MinProof(proof_min);
    proof.min_D = min;

    return std::make_pair(C, proof);
}

// Compute nested intersection with proof
std::pair<AccValue, NestedIntersectionProof> Accumulator::query_nested_intersection(const Set& a, const Set& b, AccPublicKey& pk) {
    // Calculate the intersection
    Set intersection = a.intersection(b);
    
    // Create a new intersection proof
    NestedIntersectionProof proof(pk.pairing);
    
    // Pre-compute the accumulator values for the intersection
    AccValue i_acc = setup(intersection, pk);
    /*
    * S, R
    */

    // 1. Calculate I_r =  g^I(r)
    // element_set(proof.I_r, i_acc.g_r);
    
    // 2. Calculate I_r_beta =  g^{beta * I(r)}
    element_t temp;
    element_init_G1(temp, pk.pairing);
    element_set1(temp);
    
    for (const auto& i : intersection) {
        if (i < 1 || i >= pk.q) continue;
        element_mul(temp, temp, pk.g_beta_ri[i - 1]);
    }
    element_set(proof.I_r_beta, temp);
    
    // 3. Calculate Q_s_r =  g^q(s, r)
    element_set1(temp);
    size_t idx;
    for (const auto& i : a) {
        for (const auto& j : b) {
            if (i != j) {
                idx = AccPublicKey::map_i_j_to_index(j, pk.q + i - j, pk.q);
                element_mul(temp, temp, pk.g_ri_sj[idx]);
            }
        }
    }
    element_set(proof.Q_s_r, temp);
    
    // 4. Calculate Q_s_r_delta (g^{delta·q(s,r)})
    element_set1(temp);
    for (const auto& i : a) {
        for (const auto& j : b) {
            if (i != j) {
                idx = AccPublicKey::map_i_j_to_index(j, pk.q + i - j, pk.q);
                element_mul(temp, temp, pk.g_delta_ri_sj[idx]);
            }
        }
    }
    element_set(proof.Q_s_r_delta, temp);
    
    // 5. Calculate L_r = g^{I(r) / r}
    element_set1(temp);
    for (const auto& i : intersection) {
        if (i < 1 || i >= pk.q) continue;
        if (i == 1) {
            element_mul(temp, temp, pk.g);
        } else {
            element_mul(temp, temp, pk.g_ri[i - 2]);  // g^{r^(i-1)}
        }
    }
    element_set(proof.L_r, temp);

    // Nested 1.1: Z_s_r
    element_set1(proof.Z_s_r);
    for (const auto& i : intersection) {
        for (uint64_t j = 0; j < pk.q - i; j++) {
            if (j == 0) {
                element_div(proof.Z_s_r, proof.Z_s_r, pk.g_ri[i - 1]);
            } else {
                size_t idx = AccPublicKey::map_i_j_to_index(i, j, pk.q);
                element_div(proof.Z_s_r, proof.Z_s_r, pk.g_ri_sj[idx]);
            }
        }
    }

    // Nested 1.2: I_s_r_gamma
    element_set1(proof.I_s_r_gamma);
    for (const auto& i : intersection) {
        element_mul(proof.I_s_r_gamma, proof.I_s_r_gamma, pk.g_gamma_ri_sqi[pk.q - i - 1]);
    }

    /*
    * R, S
    */
    
    // 2. Calculate I_s_alpha =  g^{alpha * I()}
    element_set1(temp);
    for (const auto& i : intersection) {
        if (i < 1 || i >= pk.q) continue;
        element_mul(temp, temp, pk.g_alpha_si[i - 1]);
    }
    element_set(proof.I_s_alpha, temp);
    
    // 3. Calculate Q_r_s =  g^q(r, s)
    element_set1(temp);
    for (const auto& i : a) {
        for (const auto& j : b) {
            if (i != j) {
                idx = AccPublicKey::map_i_j_to_index(pk.q + i - j, j, pk.q);
                element_mul(temp, temp, pk.g_ri_sj[idx]);
            }
        }
    }
    element_set(proof.Q_r_s, temp);
    
    // 4. Calculate Q_r_s_delta (g^{delta·q(r,s)})
    element_set1(temp);
    for (const auto& i : a) {
        for (const auto& j : b) {
            if (i != j) {
                idx = AccPublicKey::map_i_j_to_index(pk.q + i - j, j, pk.q);
                element_mul(temp, temp, pk.g_delta_ri_sj[idx]);
            }
        }
    }
    element_set(proof.Q_r_s_delta, temp);
    
    // 5. Calculate L_s = g^{I(s) / s}
    element_set1(temp);
    for (const auto& i : intersection) {
        if (i < 1 || i >= pk.q) continue;
        if (i == 1) {
            element_mul(temp, temp, pk.g);
        } else {
            element_mul(temp, temp, pk.g_si[i - 2]);  // g^{s^(i-1)}
        }
    }
    element_set(proof.L_s, temp);

    // Nested 2.1
    element_set1(proof.Z_r_s);
    for (const auto& i : intersection) {
        for (uint64_t j = 0; j < pk.q - i; j++) {
            if (j == 0) {
                element_div(proof.Z_r_s, proof.Z_r_s, pk.g_si[i - 1]);
            } else {
                size_t idx = AccPublicKey::map_i_j_to_index(j, i, pk.q);
                element_div(proof.Z_r_s, proof.Z_r_s, pk.g_ri_sj[idx]);
            }
        }
    }

    // Nested 2.2: I_r_s_gamma
    element_set1(proof.I_r_s_gamma);
    for (const auto& i : intersection) {
        element_mul(proof.I_r_s_gamma, proof.I_r_s_gamma, pk.g_gamma_ri_sqi[i - 1]);
    }
    
    // Clean up
    element_clear(temp);
    
    return std::make_pair(i_acc, proof);
}

// Compute the nested union and proof
std::pair<AccValue, NestedUnionProof> Accumulator::query_nested_union(const Set& a, const Set& b, AccPublicKey& pk) {
    // Calculate the union U = A ∪ B
    Set union_set = a.union_with(b);
    
    // Calculate the intersection I = A ∩ B
    Set intersection = a.intersection(b);
    
    // Create a new difference proof
    NestedUnionProof proof(pk.pairing);
    
    // Pre-compute the accumulator values
    AccValue a_acc = setup(a, pk);
    AccValue i_acc = setup(intersection, pk);
    AccValue u_acc = setup(union_set, pk);
    
    // 1. Calculate I_r = g^I(r) (accumulator for intersection)
    element_set(proof.I_r, i_acc.g_r);
    
    // 2. Calculate I_r_beta = g^{beta * I(r)}
    element_t temp;
    element_init_G1(temp, pk.pairing);
    element_set1(temp);
    
    for (const auto& i : intersection) {
        if (i < 1 || i >= pk.q) continue;
        element_mul(temp, temp, pk.g_beta_ri[i - 1]);
    }
    element_set(proof.I_r_beta, temp);
    
    // 3. Calculate Q_s_r = g^q(s, r)
    element_set1(temp);
    size_t idx;
    for (const auto& i : a) {
        for (const auto& j : b) {
            if (i != j) {
                idx = AccPublicKey::map_i_j_to_index(j, pk.q + i - j, pk.q);
                element_mul(temp, temp, pk.g_ri_sj[idx]);
            }
        }
    }
    element_set(proof.Q_s_r, temp);
    
    // 4. Calculate Q_s_r_delta (g^{delta·q(s,r)})
    element_set1(temp);
    for (const auto& i : a) {
        for (const auto& j : b) {
            if (i != j) {
                idx = AccPublicKey::map_i_j_to_index(j, pk.q + i - j, pk.q);
                element_mul(temp, temp, pk.g_delta_ri_sj[idx]);
            }
        }
    }
    element_set(proof.Q_s_r_delta, temp);
    
    // 5. Calculate L_r = g^{I(r) / r}
    element_set1(temp);
    for (const auto& i : intersection) {
        if (i < 1 || i >= pk.q) continue;
        if (i == 1) {
            element_mul(temp, temp, pk.g);
        } else {
            element_mul(temp, temp, pk.g_ri[i - 2]);  // g^{r^(i-1)}
        }
    }
    element_set(proof.L_r, temp);
    
    // 6. Z_s_r
    element_set1(proof.Z_s_r);
    for (const auto& i : union_set) {
        for (uint64_t j = 0; j < pk.q - i; j++) {
            if (j == 0) {
                element_div(proof.Z_s_r, proof.Z_s_r, pk.g_si[i - 1]);
            } else {
                size_t idx = AccPublicKey::map_i_j_to_index(j, i, pk.q);
                element_div(proof.Z_s_r, proof.Z_s_r, pk.g_ri_sj[idx]);
            }
        }
    }

    // 7. U_s_r_gamma
    element_set1(proof.U_s_r_gamma);
    for (const auto& i : union_set) {
        element_mul(proof.U_s_r_gamma, proof.U_s_r_gamma, pk.g_gamma_ri_sqi[pk.q - i - 1]);
    }

    // 8. Z_r_s
    element_set1(proof.Z_r_s);
    for (const auto& i : union_set) {
        for (uint64_t j = 0; j < pk.q - i; j++) {
            if (j == 0) {
                element_div(proof.Z_r_s, proof.Z_r_s, pk.g_ri[i - 1]);
            } else {
                size_t idx = AccPublicKey::map_i_j_to_index(i, j, pk.q);
                element_div(proof.Z_r_s, proof.Z_r_s, pk.g_ri_sj[idx]);
            }
        }
    }

    // 9. U_r_s_gamma
    element_set1(proof.U_r_s_gamma);
    for (const auto& i : union_set) {
        element_mul(proof.U_r_s_gamma, proof.U_r_s_gamma, pk.g_gamma_ri_sqi[i - 1]);
    }
    
    // Clean up
    element_clear(temp);
    
    return std::make_pair(u_acc, proof);
}


// Compute the nested difference and proof
std::pair<AccValue, NestedDifferenceProof> Accumulator::query_nested_difference(const Set& a, const Set& b, AccPublicKey& pk) {
    // Calculate the difference D = A \ B
    Set difference = a.difference(b);
    
    // Calculate the intersection I = A ∩ B
    Set intersection = a.intersection(b);
    
    // Create a new difference proof
    NestedDifferenceProof proof(pk.pairing);
    
    // Pre-compute the accumulator values
    AccValue a_acc = setup(a, pk);
    AccValue i_acc = setup(intersection, pk);
    AccValue d_acc = setup(difference, pk);
    
    // 1. Calculate I_r = g^I(r) (accumulator for intersection)
    element_set(proof.I_r, i_acc.g_r);
    
    // 2. Calculate I_r_beta = g^{beta * I(r)}
    element_t temp;
    element_init_G1(temp, pk.pairing);
    element_set1(temp);
    
    for (const auto& i : intersection) {
        if (i < 1 || i >= pk.q) continue;
        element_mul(temp, temp, pk.g_beta_ri[i - 1]);
    }
    element_set(proof.I_r_beta, temp);
    
    // 3. Calculate Q_s_r = g^q(s, r)
    element_set1(temp);
    size_t idx;
    for (const auto& i : a) {
        for (const auto& j : b) {
            if (i != j) {
                idx = AccPublicKey::map_i_j_to_index(j, pk.q + i - j, pk.q);
                element_mul(temp, temp, pk.g_ri_sj[idx]);
            }
        }
    }
    element_set(proof.Q_s_r, temp);
    
    // 4. Calculate Q_s_r_delta (g^{delta·q(s,r)})
    element_set1(temp);
    for (const auto& i : a) {
        for (const auto& j : b) {
            if (i != j) {
                idx = AccPublicKey::map_i_j_to_index(j, pk.q + i - j, pk.q);
                element_mul(temp, temp, pk.g_delta_ri_sj[idx]);
            }
        }
    }
    element_set(proof.Q_s_r_delta, temp);
    
    // 5. Calculate L_r = g^{I(r) / r}
    element_set1(temp);
    for (const auto& i : intersection) {
        if (i < 1 || i >= pk.q) continue;
        if (i == 1) {
            element_mul(temp, temp, pk.g);
        } else {
            element_mul(temp, temp, pk.g_ri[i - 2]);  // g^{r^(i-1)}
        }
    }
    element_set(proof.L_r, temp);
    
    // 6. Z_s_r
    element_set1(proof.Z_s_r);
    for (const auto& i : difference) {
        for (uint64_t j = 0; j < pk.q - i; j++) {
            if (j == 0) {
                element_div(proof.Z_s_r, proof.Z_s_r, pk.g_si[i - 1]);
            } else {
                size_t idx = AccPublicKey::map_i_j_to_index(j, i, pk.q);
                element_div(proof.Z_s_r, proof.Z_s_r, pk.g_ri_sj[idx]);
            }
        }
    }

    // 7. D_s_r_gamma
    element_set1(proof.D_s_r_gamma);
    for (const auto& i : difference) {
        element_mul(proof.D_s_r_gamma, proof.D_s_r_gamma, pk.g_gamma_ri_sqi[pk.q - i - 1]);
    }

    // 8. Z_r_s
    element_set1(proof.Z_r_s);
    for (const auto& i : difference) {
        for (uint64_t j = 0; j < pk.q - i; j++) {
            if (j == 0) {
                element_div(proof.Z_r_s, proof.Z_r_s, pk.g_ri[i - 1]);
            } else {
                size_t idx = AccPublicKey::map_i_j_to_index(i, j, pk.q);
                element_div(proof.Z_r_s, proof.Z_r_s, pk.g_ri_sj[idx]);
            }
        }
    }

    // 9. D_r_s_gamma
    element_set1(proof.D_r_s_gamma);
    for (const auto& i : difference) {
        element_mul(proof.D_r_s_gamma, proof.D_r_s_gamma, pk.g_gamma_ri_sqi[i - 1]);
    }
    
    // Clean up
    element_clear(temp);
    
    return std::make_pair(d_acc, proof);
}


// Compute range query with proof
std::pair<AccValue, NestedRangeProof> Accumulator::query_nested_range(const Set& a, uint64_t l, uint64_t r, AccPublicKey& pk) {
    // Filter the set to include only elements in the range [l, r]
    Set B, C, D;
    for (const auto& element : a) {
        if (element >= l && element <= r) {
            C.insert(element);
        } else if (element < l) {
            B.insert(element);
        } else {
            D.insert(element);
        }
    }
    
    // We reuse the intersection proof structure for simplicity
    // In a real implementation, this would use a specific RANGE proof
    NestedRangeProof proof(pk.pairing);

    AccValue B_acc = setup(B, pk);
    AccValue C_acc = setup(C, pk);
    AccValue D_acc = setup(D, pk);

    proof.d_B = B_acc;
    proof.d_D = D_acc;

    // 1.1 B_s, B_s_alpha
    // Calculate B_s = g^B(s)
    element_set(proof.B_s, B_acc.g_s);
    // Calculate B_s_alpha = g^{alpha·B(s)}
    element_set1(proof.B_s_alpha);
    for (const auto& i : B) {
        element_mul(proof.B_s_alpha, proof.B_s_alpha, pk.g_alpha_si[i - 1]);
    }

    // 1.2 B_s_r, B_s_r_gamma
    // Calculate B_s_r = g^B(s,r)
    element_set(proof.B_s_r, B_acc.g_s_r);
    // Calculate B_s_r_gamma = g^{gamma·B(s,r)}
    element_set1(proof.B_s_r_gamma);
    for (const auto& i : B) {
        element_mul(proof.B_s_r_gamma, proof.B_s_r_gamma, pk.g_gamma_ri_sqi[pk.q - i - 1]);
    }

    // 1.3 D_s, D_s_alpha
    // Calculate D_s = g^D(s)
    element_set(proof.D_s, D_acc.g_s);
    // Calculate D_s_alpha = g^{alpha·D(s)}
    element_set1(proof.D_s_alpha);
    for (const auto& i : D) {
        element_mul(proof.D_s_alpha, proof.D_s_alpha, pk.g_alpha_si[i - 1]);
    }

    // 2. Z_s_r
    // Calculate Z_s_r = g^Z(s,r)
    element_set1(proof.Z_s_r);
    for (const auto& i : B) {
        for (uint64_t j = 0; j < pk.q - i; j++) {
            if (j == 0) {
                element_div(proof.Z_s_r, proof.Z_s_r, pk.g_si[i - 1]);
            } else {
                size_t idx = AccPublicKey::map_i_j_to_index(j, i, pk.q);
                element_div(proof.Z_s_r, proof.Z_s_r, pk.g_ri_sj[idx]);
            }
        }
    }

    // 3. pi_BC
    auto [BC, proof_BC] = query_intersection(B, C, pk);
    proof.pi_BC = IntersectionProof(proof_BC);

    // 4. pi_BD
    auto [BD, proof_BD] = query_intersection(B, D, pk);
    proof.pi_BD = IntersectionProof(proof_BD);

    // 5. pi_CD
    auto [CD, proof_CD] = query_intersection(C, D, pk);
    proof.pi_CD = IntersectionProof(proof_CD);

    // 6. pi_1
    auto [max, proof_max] = query_max(B, pk);
    proof.pi_1 = MaxProof(proof_max);
    proof.max_B = max;

    // 7. pi_2
    auto [min, proof_min] = query_min(D, pk);
    proof.pi_2 = MinProof(proof_min);
    proof.min_D = min;

    // 8. ZC_s_r, ZC_r_s
    element_set1(proof.ZC_s_r);
    for (const auto& i : C) {
        for (uint64_t j = 0; j < pk.q - i; j++) {
            if (j == 0) {
                element_div(proof.ZC_s_r, proof.ZC_s_r, pk.g_si[i - 1]);
            } else {
                size_t idx = AccPublicKey::map_i_j_to_index(j, i, pk.q);
                element_div(proof.ZC_s_r, proof.ZC_s_r, pk.g_ri_sj[idx]);
            }
        }
    }

    element_set1(proof.ZC_r_s);
    for (const auto& i : C) {
        for (uint64_t j = 0; j < pk.q - i; j++) {
            if (j == 0) {
                element_div(proof.ZC_r_s, proof.ZC_r_s, pk.g_ri[i - 1]);
            } else {
                size_t idx = AccPublicKey::map_i_j_to_index(i, j, pk.q);
                element_div(proof.ZC_r_s, proof.ZC_r_s, pk.g_ri_sj[idx]);
            }
        }
    }
    
    // 9. pi_C_max
    auto [max_C, proof_max_C] = query_max(C, pk);
    proof.pi_C_max = MaxProof(proof_max_C);
    proof.max_C = max_C;

    // 10. pi_C_min
    auto [min_C, proof_min_C] = query_min(C, pk);
    proof.pi_C_min = MinProof(proof_min_C);
    proof.min_C = min_C;
    
    return std::make_pair(C_acc, proof);
}

// Verify a proof
bool Accumulator::verify(AccValue& a_acc, AccValue& b_acc, Proof& proof, const Set& result, AccPublicKey& pk) {
    return proof.verify(a_acc, b_acc, result, pk);
}

// Verify a proof for numeric result
bool Accumulator::verify_aggr(AccValue& acc, AggrProof& proof, uint64_t result, AccPublicKey& pk) {
    return proof.verify_aggr(acc, result, pk);
}

// Verify a range proof
bool Accumulator::verify_range(AccValue& acc, RangeProof& proof, const Set& result, uint64_t l, uint64_t r, AccPublicKey& pk) {
    return proof.verify_range(acc, result, l, r, pk);
}

// Verify a nested intersection proof
bool Accumulator::verify_nested(AccValue& lhs_acc, AccValue& rhs_acc, NestedProof& proof, AccValue& res_acc, AccPublicKey& pk) {
    return proof.verify_nested(lhs_acc, rhs_acc, res_acc, pk);
}

// Verify a nested range proof
bool Accumulator::verify_nested_range(AccValue& acc, NestedRangeProof& proof, AccValue& res_acc, uint64_t l, uint64_t r, AccPublicKey& pk) {
    return proof.verify_nested_range(acc, res_acc, l, r, pk);
}

} // namespace acc 