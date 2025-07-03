/**
 * @file keys.cpp
 * @brief Implementation of key classes for the Expressive Set Accumulator
 */

#include <cstdlib>  // For random number generation
#include <iostream>  // For debug output
#include <cstring>  // For memcpy
#include <thread>   // For multithreading support
#include <vector>   // For std::vector
#include <mutex>    // For std::mutex
#include <algorithm> // For std::min

#include "acc/keys.h"

namespace acc {

// AccSecretKey implementation
AccSecretKey::AccSecretKey(char* _param, size_t _count) {
    // 确保count不超过固定数组大小
    count = _count;
    memcpy(param, _param, count);
    
    // Store the pairing
    pairing_init_set_buf(pairing, param, count);
    
    // Initialize elements
    element_init_Zr(s, pairing);
    element_init_Zr(r, pairing);
    element_init_Zr(alpha, pairing);
    element_init_Zr(beta, pairing);
    element_init_Zr(gamma, pairing);
    element_init_Zr(delta, pairing);
}

AccSecretKey::AccSecretKey(const AccSecretKey& other) {
    // 复制count和param数组内容
    count = other.count;
    memcpy(param, other.param, count);
    
    // 初始化自己的pairing
    pairing_init_set_buf(pairing, param, count);
    
    // 初始化元素
    element_init_Zr(s, pairing);
    element_init_Zr(r, pairing);
    element_init_Zr(alpha, pairing);
    element_init_Zr(beta, pairing);
    element_init_Zr(gamma, pairing);
    element_init_Zr(delta, pairing);
    
    // PBC functions need non-const element_t
    element_t& non_const_s = const_cast<element_t&>(other.s);
    element_t& non_const_r = const_cast<element_t&>(other.r);
    element_t& non_const_alpha = const_cast<element_t&>(other.alpha);
    element_t& non_const_beta = const_cast<element_t&>(other.beta);
    element_t& non_const_gamma = const_cast<element_t&>(other.gamma);
    element_t& non_const_delta = const_cast<element_t&>(other.delta);
    
    // 复制元素值
    element_set(s, non_const_s);
    element_set(r, non_const_r);
    element_set(alpha, non_const_alpha);
    element_set(beta, non_const_beta);
    element_set(gamma, non_const_gamma);
    element_set(delta, non_const_delta);
}

AccSecretKey::~AccSecretKey() {
    // Clear elements
    element_clear(s);
    element_clear(r);
    element_clear(alpha);
    element_clear(beta);
    element_clear(gamma);
    element_clear(delta);
    
    // Clear pairing
    pairing_clear(pairing);
}

AccSecretKey AccSecretKey::random(char* param, size_t count) {
    AccSecretKey sk(param, count);
    
    // Set to random values
    element_random(sk.s);
    element_random(sk.r);
    element_random(sk.alpha);
    element_random(sk.beta);
    element_random(sk.gamma);
    element_random(sk.delta);
    
    return sk;
}

// AccPublicKey implementation
AccPublicKey::AccPublicKey(char* _param, size_t _count, uint64_t q) : q(q) {
    auto start_time = std::chrono::high_resolution_clock::now();

    // --- Sequential Setup ---
    count = _count;
    memcpy(param, _param, count);

    // Store the pairing
    pairing_init_set_buf(pairing, param, count);

    // Initialize the single generator and elements sequentially
    element_init_G1(g, pairing);
    element_init_G1(g_alpha, pairing);
    element_init_G1(g_beta, pairing);
    element_init_G1(g_gamma, pairing);
    element_init_G1(g_delta, pairing);
    element_init_G1(g_sq, pairing);
    element_init_G1(g_rq, pairing);

    // --- Parallel Initialization of 1D Arrays ---
    if (q > 1) {
        // Allocate memory sequentially
        g_si = new element_t[q - 1];
        g_alpha_si = new element_t[q - 1];
        g_ri = new element_t[q - 1];
        g_beta_ri = new element_t[q - 1];
        g_gamma_ri_sqi = new element_t[q - 1];

        // Determine thread configuration
        const uint64_t total_elements_1d = q - 1;
        // Use hardware_concurrency, ensuring at least 1 thread
        const unsigned int num_threads = std::max(1u, std::thread::hardware_concurrency());
        const uint64_t elements_per_thread = (total_elements_1d + num_threads - 1) / num_threads; // Ceiling division

        std::vector<std::thread> threads_1d;
        threads_1d.reserve(num_threads); // Reserve space for threads

        // Lambda for initializing a chunk of 1D arrays
        // Capture pairing by value/reference safely as it's initialized before threads start
        auto init_1d_chunk = [&](uint64_t start_idx, uint64_t end_idx) {
            for (uint64_t i = start_idx; i < end_idx; ++i) {
                element_init_G1(g_si[i], pairing);
                element_init_G1(g_alpha_si[i], pairing);
                element_init_G1(g_ri[i], pairing);
                element_init_G1(g_beta_ri[i], pairing);
                element_init_G1(g_gamma_ri_sqi[i], pairing);
            }
        };

        // Launch threads for 1D initialization
        for (unsigned int t = 0; t < num_threads; ++t) {
            uint64_t start = t * elements_per_thread;
            uint64_t end = std::min(start + elements_per_thread, total_elements_1d);
            if (start < end) { // Ensure we don't launch threads for empty ranges
                threads_1d.emplace_back(init_1d_chunk, start, end);
            }
        }

        // Join threads
        for (auto& thread : threads_1d) {
            if (thread.joinable()) {
                thread.join();
            }
        }
    } else {
        // Handle q <= 1 case (no 1D arrays needed)
        g_si = nullptr;
        g_alpha_si = nullptr;
        g_ri = nullptr;
        g_beta_ri = nullptr;
        g_gamma_ri_sqi = nullptr;
    }

    // --- Parallel Initialization of 2D Arrays ---
     if (q > 1) {
        size_t total_pairs = (2 * q - 2) * (2 * q - 2);
        // Allocate memory sequentially
        g_ri_sj = new element_t[total_pairs];
        g_delta_ri_sj = new element_t[total_pairs];

        // Determine thread configuration (can reuse num_threads or recalculate)
        const unsigned int num_threads = std::max(1u, std::thread::hardware_concurrency());
        const uint64_t pairs_per_thread = (total_pairs + num_threads - 1) / num_threads; // Ceiling division

        std::vector<std::thread> threads_2d;
        threads_2d.reserve(num_threads);

         // Lambda for initializing a chunk of 2D arrays
        auto init_2d_chunk = [&](uint64_t start_idx, uint64_t end_idx) {
            for (size_t idx = start_idx; idx < end_idx; ++idx) {
                element_init_G1(g_ri_sj[idx], pairing);
                element_init_G1(g_delta_ri_sj[idx], pairing);
            }
        };

        // Launch threads for 2D initialization
        for (unsigned int t = 0; t < num_threads; ++t) {
            uint64_t start = t * pairs_per_thread;
            uint64_t end = std::min(start + pairs_per_thread, total_pairs);
             if (start < end) { // Ensure we don't launch threads for empty ranges
                threads_2d.emplace_back(init_2d_chunk, start, end);
             }
        }

        // Join threads
        for (auto& thread : threads_2d) {
             if (thread.joinable()) {
                thread.join();
             }
        }
     } else {
        // Handle q <= 1 case (no 2D arrays needed)
        g_ri_sj = nullptr;
        g_delta_ri_sj = nullptr;
     }

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    // Update print statement to reflect multithreading
    std::cout << "Time for public key initialization (multi-threaded): " << duration.count() / 1000.0 << " seconds" << std::endl;
}

AccPublicKey::AccPublicKey(const AccPublicKey& other) {
    auto start_time = std::chrono::high_resolution_clock::now();

    // Copy elements - using const_cast because PBC functions don't accept const
    q = other.q;
    count = other.count;
    memcpy(param, other.param, count);
    
    // Store the pairing
    pairing_init_set_buf(pairing, param, count);
    
    // 初始化其他元素，使用自己的pairing
    element_init_G1(g, pairing);
    element_init_G1(g_alpha, pairing);
    element_init_G1(g_beta, pairing);
    element_init_G1(g_gamma, pairing);
    element_init_G1(g_delta, pairing);

    element_init_G1(g_sq, pairing);
    element_init_G1(g_rq, pairing);

    g_si = new element_t[q-1];
    g_alpha_si = new element_t[q-1];
    g_ri = new element_t[q-1];
    g_beta_ri = new element_t[q-1];
    g_gamma_ri_sqi = new element_t[q-1];
    for (uint64_t i = 0; i < q-1; i++) {
        element_init_G1(g_si[i], pairing);
        element_init_G1(g_alpha_si[i], pairing);
        element_init_G1(g_ri[i], pairing);
        element_init_G1(g_beta_ri[i], pairing);
        element_init_G1(g_gamma_ri_sqi[i], pairing);
    }

    // Initialize the flattened 2D arrays
    size_t total_pairs = (2 * q - 2) * (2 * q - 2);
    g_ri_sj = new element_t[total_pairs];
    g_delta_ri_sj = new element_t[total_pairs];
    for (size_t i = 0; i < total_pairs; i++) {
        element_init_G1(g_ri_sj[i], pairing);
        element_init_G1(g_delta_ri_sj[i], pairing);
    }

    // Create non-const references for PBC functions to copy elements
    element_t& non_const_g = const_cast<element_t&>(other.g);
    element_t& non_const_g_alpha = const_cast<element_t&>(other.g_alpha);
    element_t& non_const_g_beta = const_cast<element_t&>(other.g_beta);
    element_t& non_const_g_gamma = const_cast<element_t&>(other.g_gamma);
    element_t& non_const_g_delta = const_cast<element_t&>(other.g_delta);
    element_t& non_const_g_sq = const_cast<element_t&>(other.g_sq);
    element_t& non_const_g_rq = const_cast<element_t&>(other.g_rq);

    // Copy elements
    element_set(g, non_const_g);
    element_set(g_alpha, non_const_g_alpha);
    element_set(g_beta, non_const_g_beta);
    element_set(g_gamma, non_const_g_gamma);
    element_set(g_delta, non_const_g_delta);
    element_set(g_sq, non_const_g_sq);
    element_set(g_rq, non_const_g_rq);

    for (uint64_t i = 0; i < q-1; i++) {
        element_t& non_const_g_si = const_cast<element_t&>(other.g_si[i]);
        element_t& non_const_g_alpha_si = const_cast<element_t&>(other.g_alpha_si[i]);
        element_t& non_const_g_ri = const_cast<element_t&>(other.g_ri[i]);
        element_t& non_const_g_beta_ri = const_cast<element_t&>(other.g_beta_ri[i]);
        element_t& non_const_g_gamma_ri_sqi = const_cast<element_t&>(other.g_gamma_ri_sqi[i]);
        
        element_set(g_si[i], non_const_g_si);
        element_set(g_alpha_si[i], non_const_g_alpha_si);
        element_set(g_ri[i], non_const_g_ri);
        element_set(g_beta_ri[i], non_const_g_beta_ri);
        element_set(g_gamma_ri_sqi[i], non_const_g_gamma_ri_sqi);
    }

    for (size_t i = 0; i < total_pairs; i++) {
        element_t& non_const_g_ri_sj = const_cast<element_t&>(other.g_ri_sj[i]);
        element_t& non_const_g_delta_ri_sj = const_cast<element_t&>(other.g_delta_ri_sj[i]);
        
        element_set(g_ri_sj[i], non_const_g_ri_sj);
        element_set(g_delta_ri_sj[i], non_const_g_delta_ri_sj);
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    std::cout << "Time for public key copy: " << duration.count() / 1000.0 << " seconds" << std::endl;
}

AccPublicKey::~AccPublicKey() {
    // Clear all elements
    element_clear(g);
    element_clear(g_alpha);
    element_clear(g_beta);
    element_clear(g_gamma);
    element_clear(g_delta);

    element_clear(g_sq);
    element_clear(g_rq);
    
    for (uint64_t i = 0; i < q-1; i++) {
        element_clear(g_si[i]);
        element_clear(g_alpha_si[i]);
        element_clear(g_ri[i]);
        element_clear(g_beta_ri[i]);
        element_clear(g_gamma_ri_sqi[i]);
    }

    size_t total_pairs = (2 * q - 2) * (2 * q - 2);
    for (size_t i = 0; i < total_pairs; i++) {
        element_clear(g_ri_sj[i]);
        element_clear(g_delta_ri_sj[i]);
    }
    
    // Clear the pairing
    pairing_clear(pairing);
}

// Helper function to calculate mapping from i,j to flat index for 2D array
size_t AccPublicKey::map_i_j_to_index(uint64_t i, uint64_t j, uint64_t q) {
    // Adjust indices to skip q
    uint64_t adj_i = (i < q) ? i - 1 : i - 2;
    uint64_t adj_j = (j < q) ? j - 1 : j - 2;
    
    // For a valid mapping, indices should be in range [0, 2*q-3]
    if (adj_i >= (2 * q - 2) || adj_j >= (2 * q - 2)) {
        std::cout << "  WARNING: Index out of bounds, returning 0" << std::endl;
        // Return a safe index if out of bounds (this is a guard against invalid input)
        return 0;
    }
    
    // Calculate flat index
    size_t idx = adj_i * (2 * q - 2) + adj_j;
    return idx;
}

// Generate public key from secret key
AccPublicKey AccPublicKey::generate(AccSecretKey& sk, char* param, size_t count, uint64_t q) {
    AccPublicKey pk(param, count, q);

    auto start_time = std::chrono::high_resolution_clock::now();
    
    // Set generator to random value (only one generator now)
    element_random(pk.g);
    // Precompute the generator
    element_pp_t g_pp;
    element_pp_init(g_pp, pk.g);
    
    // Compute g^alpha, g^beta, g^gamma, g^delta
    element_pp_pow_zn(pk.g_alpha, sk.alpha, g_pp);
    element_pp_pow_zn(pk.g_beta, sk.beta, g_pp);
    element_pp_pow_zn(pk.g_gamma, sk.gamma, g_pp);
    element_pp_pow_zn(pk.g_delta, sk.delta, g_pp);
    
    // precompute s^i, r^i
    element_t *si, *ri;
    si = new element_t[2 * q - 1];
    ri = new element_t[2 * q - 1];
    
    // Initialize the elements in the arrays
    for (uint64_t i = 0; i < 2 * q - 1; i++) {
        element_init_Zr(si[i], pk.pairing);
        element_init_Zr(ri[i], pk.pairing);
    }

    // Temporary variables for computations
    element_t tmp_c_x; // For calculations
    element_t tmp_xi;  // For powers
    element_t tmp_rs;  // For r^i * s^j
    
    element_init_Zr(tmp_c_x, pk.pairing);
    element_init_Zr(tmp_xi, pk.pairing);
    element_init_Zr(tmp_rs, pk.pairing);
    
    // precompute s^i, and compute g^(s^i), g^(alpha*s^i)
    element_set1(tmp_xi);
    for (uint64_t i = 1; i < 2 * q; i++) {
        // s^i: i \in [1, 2*q-1]
        element_mul(tmp_xi, sk.s, tmp_xi);
        element_set(si[i - 1], tmp_xi);

        if (i < q) {
            // g^(s^i)
            element_pp_pow_zn(pk.g_si[i - 1], si[i - 1], g_pp);
            // alpha * s^i
            element_mul(tmp_c_x, sk.alpha, si[i - 1]);
            // g^(alpha*s^i)
            element_pp_pow_zn(pk.g_alpha_si[i - 1], tmp_c_x, g_pp);
        }
        else if (i == q) {
            // g^(s^q)
            element_pp_pow_zn(pk.g_sq, si[i - 1], g_pp);
        }
        // else: i \in [q+1, 2*q-1], precompute
    }
    
    // precompute r^i, and compute g^(r^i), g^(beta*r^i), g^(gamma*r^i*s^(q-i))
    element_set1(tmp_xi);
    for (uint64_t i = 1; i < 2 * q; i++) {
        // r^i
        element_mul(tmp_xi, sk.r, tmp_xi);
        element_set(ri[i - 1], tmp_xi);
        
        if (i < q) {
            // g^(r^i)
            element_pp_pow_zn(pk.g_ri[i - 1], ri[i - 1], g_pp);
            // beta * r^i
            element_mul(tmp_c_x, sk.beta, ri[i - 1]);
            // g^(beta*r^i)
            element_pp_pow_zn(pk.g_beta_ri[i - 1], tmp_c_x, g_pp);
            // gamma * r^i * s^(q-i)
            element_mul(tmp_c_x, sk.gamma, ri[i - 1]);
            element_mul(tmp_c_x, tmp_c_x, si[q - i - 1]);
            // g^(gamma*r^i*s^(q-i))
            element_pp_pow_zn(pk.g_gamma_ri_sqi[i - 1], tmp_c_x, g_pp);
        }
        else if (i == q) {
            // g^(r^q)
            element_pp_pow_zn(pk.g_rq, ri[i - 1], g_pp);
        }
        // else: i \in [q+1, 2*q-1], precompute
    }

    // Compute g^{r^i·s^j} for (i,j) ∈ ([2q-1] \ {q}) × ([2q-1] \ {q})
    // Initialize all 2D array elements to identity first
    size_t total_pairs = (2 * q - 2) * (2 * q - 2);
    // std::cout << "Initializing " << total_pairs << " 2D array elements to identity" << std::endl;
    for (size_t idx = 0; idx < total_pairs; idx++) {
        element_set1(pk.g_ri_sj[idx]);
        element_set1(pk.g_delta_ri_sj[idx]);
    }

    for (uint64_t i = 1; i < 2 * q; i++) {
        if (i == q) continue;  // Skip q
        
        for (uint64_t j = 1; j < 2 * q; j++) {
            if (j == q) continue;  // Skip q
            
            // r^i·s^j
            element_mul(tmp_rs, ri[i - 1], si[j - 1]);
            // Compute index in the flattened array
            size_t idx = map_i_j_to_index(i, j, q);
            
            // g^{r^i·s^j}
            element_pp_pow_zn(pk.g_ri_sj[idx], tmp_rs, g_pp);            
            // delta * r^i·s^j
            element_mul(tmp_c_x, sk.delta, tmp_rs);
            // g^{delta·r^i·s^j}
            element_pp_pow_zn(pk.g_delta_ri_sj[idx], tmp_c_x, g_pp);
        }
    }
    
    // Clean up
    element_pp_clear(g_pp);
    element_clear(tmp_c_x);
    element_clear(tmp_xi);
    element_clear(tmp_rs);
    
    // Free allocated memory
    for (uint64_t i = 0; i < 2 * q - 1; i++) {
        element_clear(si[i]);
        element_clear(ri[i]);
    }
    delete[] si;
    delete[] ri;
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    std::cout << "Time for public key generation: " << duration.count() / 1000.0 << " seconds" << std::endl;
    
    return pk;
}

// Helper function for parallel computation
void compute_powers(element_t* result, element_t& base, uint64_t start, uint64_t end, element_pp_t& pp, pairing_t pairing) {
    element_t tmp;
    element_init_Zr(tmp, pairing);
    element_set1(tmp);
    
    for (uint64_t i = start; i < end; i++) {
        element_mul(tmp, tmp, base);
        element_pp_pow_zn(result[i], tmp, pp);
    }
    
    element_clear(tmp);
}

AccPublicKey AccPublicKey::generate_mul_thread(AccSecretKey& sk, char* param, size_t count, uint64_t q) {
    AccPublicKey pk(param, count, q);
    auto start_time = std::chrono::high_resolution_clock::now();
    
    // Set generator to random value
    element_random(pk.g);
    element_pp_t g_pp;
    element_pp_init(g_pp, pk.g);
    
    // Initialize arrays for precomputation
    element_t *si, *ri;
    si = new element_t[2 * q - 1];
    ri = new element_t[2 * q - 1];
    
    // Initialize all elements first
    for (uint64_t i = 0; i < 2 * q - 1; i++) {
        element_init_Zr(si[i], pk.pairing);
        element_init_Zr(ri[i], pk.pairing);
    }

    // Temporary variables for computations
    element_t tmp_xi;
    element_init_Zr(tmp_xi, pk.pairing);
    
    // Step 1: Sequentially compute s^i and r^i
    // Compute s^i
    element_set1(tmp_xi);
    for (uint64_t i = 1; i < 2 * q; i++) {
        element_mul(tmp_xi, sk.s, tmp_xi);
        element_set(si[i - 1], tmp_xi);
    }
    
    // Compute r^i
    element_set1(tmp_xi);
    for (uint64_t i = 1; i < 2 * q; i++) {
        element_mul(tmp_xi, sk.r, tmp_xi);
        element_set(ri[i - 1], tmp_xi);
    }
    
    // Step 2: Parallel computation of basic elements
    std::vector<std::thread> threads;
    std::mutex mutex;
    
    // Compute g^alpha, g^beta, g^gamma, g^delta in parallel
    threads.emplace_back([&]() { element_pp_pow_zn(pk.g_alpha, sk.alpha, g_pp); });
    threads.emplace_back([&]() { element_pp_pow_zn(pk.g_beta, sk.beta, g_pp); });
    threads.emplace_back([&]() { element_pp_pow_zn(pk.g_gamma, sk.gamma, g_pp); });
    threads.emplace_back([&]() { element_pp_pow_zn(pk.g_delta, sk.delta, g_pp); });
    
    // Wait for basic elements to complete
    for (auto& thread : threads) {
        thread.join();
    }
    threads.clear();
    
    // Step 3: Parallel computation of g_si, g_alpha_si, g_ri, g_beta_ri, g_gamma_ri_sqi
    const int num_threads = std::thread::hardware_concurrency();
    const uint64_t elements_per_thread = (q - 1) / num_threads + 1;
    
    for (int t = 0; t < num_threads; t++) {
        uint64_t start = t * elements_per_thread;
        uint64_t end = std::min(start + elements_per_thread, q - 1);
        
        threads.emplace_back([&, start, end]() {
            element_t tmp_c_x;
            element_init_Zr(tmp_c_x, pk.pairing);
            
            for (uint64_t i = start; i < end; i++) {
                // g^(s^i)
                element_pp_pow_zn(pk.g_si[i], si[i], g_pp);
                
                // g^(alpha*s^i)
                element_mul(tmp_c_x, sk.alpha, si[i]);
                element_pp_pow_zn(pk.g_alpha_si[i], tmp_c_x, g_pp);
                
                // g^(r^i)
                element_pp_pow_zn(pk.g_ri[i], ri[i], g_pp);
                
                // g^(beta*r^i)
                element_mul(tmp_c_x, sk.beta, ri[i]);
                element_pp_pow_zn(pk.g_beta_ri[i], tmp_c_x, g_pp);
                
                // g^(gamma*r^i*s^(q-i-1))
                element_mul(tmp_c_x, sk.gamma, ri[i]);
                element_mul(tmp_c_x, tmp_c_x, si[q - i - 2]);
                element_pp_pow_zn(pk.g_gamma_ri_sqi[i], tmp_c_x, g_pp);
            }
            
            element_clear(tmp_c_x);
        });
    }
    
    // Compute g^(s^q) and g^(r^q)
    element_pp_pow_zn(pk.g_sq, si[q - 1], g_pp);
    element_pp_pow_zn(pk.g_rq, ri[q - 1], g_pp);
    
    // Wait for all threads to complete
    for (auto& thread : threads) {
        thread.join();
    }
    threads.clear();
    
    // Step 4: Parallel computation of g^{r^i·s^j} and g^{delta·r^i·s^j}
    const uint64_t total_pairs = (2 * q - 2) * (2 * q - 2);
    const uint64_t pairs_per_thread = total_pairs / num_threads + 1;
    
    // Initialize all 2D array elements to identity first, like in generate method
    for (size_t idx = 0; idx < total_pairs; idx++) {
        element_set1(pk.g_ri_sj[idx]);
        element_set1(pk.g_delta_ri_sj[idx]);
    }
    
    for (int t = 0; t < num_threads; t++) {
        uint64_t start = t * pairs_per_thread;
        uint64_t end = std::min(start + pairs_per_thread, total_pairs);
        
        threads.emplace_back([&, start, end]() {
            element_t tmp_rs, tmp_c_x;
            element_init_Zr(tmp_rs, pk.pairing);
            element_init_Zr(tmp_c_x, pk.pairing);
            
            for (uint64_t i = 1; i < 2 * q; i++) {
                if (i == q) continue;  // Skip q
                
                for (uint64_t j = 1; j < 2 * q; j++) {
                    if (j == q) continue;  // Skip q
                    
                    // 使用map_i_j_to_index计算正确的索引
                    size_t idx = map_i_j_to_index(i, j, q);
                    
                    // 检查这个索引是否在当前线程的处理范围内
                    if (idx < start || idx >= end) continue;
                    
                    // r^i·s^j
                    element_mul(tmp_rs, ri[i - 1], si[j - 1]);
                    
                    // g^{r^i·s^j}
                    element_pp_pow_zn(pk.g_ri_sj[idx], tmp_rs, g_pp);
                    
                    // delta * r^i·s^j
                    element_mul(tmp_c_x, sk.delta, tmp_rs);
                    // g^{delta·r^i·s^j}
                    element_pp_pow_zn(pk.g_delta_ri_sj[idx], tmp_c_x, g_pp);
                }
            }
            
            element_clear(tmp_rs);
            element_clear(tmp_c_x);
        });
    }
    
    // Wait for all threads to complete
    for (auto& thread : threads) {
        thread.join();
    }
    
    // Clean up
    element_pp_clear(g_pp);
    element_clear(tmp_xi);
    
    for (uint64_t i = 0; i < 2 * q - 1; i++) {
        element_clear(si[i]);
        element_clear(ri[i]);
    }
    delete[] si;
    delete[] ri;
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    std::cout << "Time for public key generation: " << duration.count() / 1000.0 << " seconds" << std::endl;
    
    return pk;
}


AccPublicKey AccPublicKey::generate_mul_thread(element_t g, AccSecretKey& sk, char* param, size_t count, uint64_t q) {
    AccPublicKey pk(param, count, q);
    auto start_time = std::chrono::high_resolution_clock::now();
    
    // Set generator to random value
    // element_random(pk.g);
    element_set(pk.g, g);
    element_pp_t g_pp;
    element_pp_init(g_pp, pk.g);
    
    // Initialize arrays for precomputation
    element_t *si, *ri;
    si = new element_t[2 * q - 1];
    ri = new element_t[2 * q - 1];
    
    // Initialize all elements first
    for (uint64_t i = 0; i < 2 * q - 1; i++) {
        element_init_Zr(si[i], pk.pairing);
        element_init_Zr(ri[i], pk.pairing);
    }

    // Temporary variables for computations
    element_t tmp_xi;
    element_init_Zr(tmp_xi, pk.pairing);
    
    // Step 1: Sequentially compute s^i and r^i
    // Compute s^i
    element_set1(tmp_xi);
    for (uint64_t i = 1; i < 2 * q; i++) {
        element_mul(tmp_xi, sk.s, tmp_xi);
        element_set(si[i - 1], tmp_xi);
    }
    
    // Compute r^i
    element_set1(tmp_xi);
    for (uint64_t i = 1; i < 2 * q; i++) {
        element_mul(tmp_xi, sk.r, tmp_xi);
        element_set(ri[i - 1], tmp_xi);
    }
    
    // Step 2: Parallel computation of basic elements
    std::vector<std::thread> threads;
    std::mutex mutex;
    
    // Compute g^alpha, g^beta, g^gamma, g^delta in parallel
    threads.emplace_back([&]() { element_pp_pow_zn(pk.g_alpha, sk.alpha, g_pp); });
    threads.emplace_back([&]() { element_pp_pow_zn(pk.g_beta, sk.beta, g_pp); });
    threads.emplace_back([&]() { element_pp_pow_zn(pk.g_gamma, sk.gamma, g_pp); });
    threads.emplace_back([&]() { element_pp_pow_zn(pk.g_delta, sk.delta, g_pp); });
    
    // Wait for basic elements to complete
    for (auto& thread : threads) {
        thread.join();
    }
    threads.clear();
    
    // Step 3: Parallel computation of g_si, g_alpha_si, g_ri, g_beta_ri, g_gamma_ri_sqi
    const int num_threads = std::thread::hardware_concurrency();
    const uint64_t elements_per_thread = (q - 1) / num_threads + 1;
    
    for (int t = 0; t < num_threads; t++) {
        uint64_t start = t * elements_per_thread;
        uint64_t end = std::min(start + elements_per_thread, q - 1);
        
        threads.emplace_back([&, start, end]() {
            element_t tmp_c_x;
            element_init_Zr(tmp_c_x, pk.pairing);
            
            for (uint64_t i = start; i < end; i++) {
                // g^(s^i)
                element_pp_pow_zn(pk.g_si[i], si[i], g_pp);
                
                // g^(alpha*s^i)
                element_mul(tmp_c_x, sk.alpha, si[i]);
                element_pp_pow_zn(pk.g_alpha_si[i], tmp_c_x, g_pp);
                
                // g^(r^i)
                element_pp_pow_zn(pk.g_ri[i], ri[i], g_pp);
                
                // g^(beta*r^i)
                element_mul(tmp_c_x, sk.beta, ri[i]);
                element_pp_pow_zn(pk.g_beta_ri[i], tmp_c_x, g_pp);
                
                // g^(gamma*r^i*s^(q-i-1))
                element_mul(tmp_c_x, sk.gamma, ri[i]);
                element_mul(tmp_c_x, tmp_c_x, si[q - i - 2]);
                element_pp_pow_zn(pk.g_gamma_ri_sqi[i], tmp_c_x, g_pp);
            }
            
            element_clear(tmp_c_x);
        });
    }
    
    // Compute g^(s^q) and g^(r^q)
    element_pp_pow_zn(pk.g_sq, si[q - 1], g_pp);
    element_pp_pow_zn(pk.g_rq, ri[q - 1], g_pp);
    
    // Wait for all threads to complete
    for (auto& thread : threads) {
        thread.join();
    }
    threads.clear();
    
    // Step 4: Parallel computation of g^{r^i·s^j} and g^{delta·r^i·s^j}
    const uint64_t total_pairs = (2 * q - 2) * (2 * q - 2);
    const uint64_t pairs_per_thread = total_pairs / num_threads + 1;
    
    // Initialize all 2D array elements to identity first, like in generate method
    for (size_t idx = 0; idx < total_pairs; idx++) {
        element_set1(pk.g_ri_sj[idx]);
        element_set1(pk.g_delta_ri_sj[idx]);
    }
    
    for (int t = 0; t < num_threads; t++) {
        uint64_t start = t * pairs_per_thread;
        uint64_t end = std::min(start + pairs_per_thread, total_pairs);
        
        threads.emplace_back([&, start, end]() {
            element_t tmp_rs, tmp_c_x;
            element_init_Zr(tmp_rs, pk.pairing);
            element_init_Zr(tmp_c_x, pk.pairing);
            
            for (uint64_t i = 1; i < 2 * q; i++) {
                if (i == q) continue;  // Skip q
                
                for (uint64_t j = 1; j < 2 * q; j++) {
                    if (j == q) continue;  // Skip q
                    
                    // 使用map_i_j_to_index计算正确的索引
                    size_t idx = map_i_j_to_index(i, j, q);
                    
                    // 检查这个索引是否在当前线程的处理范围内
                    if (idx < start || idx >= end) continue;
                    
                    // r^i·s^j
                    element_mul(tmp_rs, ri[i - 1], si[j - 1]);
                    
                    // g^{r^i·s^j}
                    element_pp_pow_zn(pk.g_ri_sj[idx], tmp_rs, g_pp);
                    
                    // delta * r^i·s^j
                    element_mul(tmp_c_x, sk.delta, tmp_rs);
                    // g^{delta·r^i·s^j}
                    element_pp_pow_zn(pk.g_delta_ri_sj[idx], tmp_c_x, g_pp);
                }
            }
            
            element_clear(tmp_rs);
            element_clear(tmp_c_x);
        });
    }
    
    // Wait for all threads to complete
    for (auto& thread : threads) {
        thread.join();
    }
    
    // Clean up
    element_pp_clear(g_pp);
    element_clear(tmp_xi);
    
    for (uint64_t i = 0; i < 2 * q - 1; i++) {
        element_clear(si[i]);
        element_clear(ri[i]);
    }
    delete[] si;
    delete[] ri;
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    std::cout << "Time for public key generation: " << duration.count() / 1000.0 << " seconds" << std::endl;
    
    return pk;
}

} // namespace acc 