/**
 * @file keys_io.cpp
 * @brief Implementation of serialization and deserialization for key classes
 */

#include "acc/keys.h"
#include <cstring>  // For memcpy
#include <iostream> // For debug messages
#include <stdexcept> // For exceptions
#include <chrono>  // For timing
#include <thread>  // For multithreading
#include <vector>  // For std::vector
#include <mutex>   // For std::mutex
#include <algorithm> // For std::min

namespace acc {

// Helper function to write size_t to bytes
void write_size_t(std::vector<unsigned char>& data, size_t value) {
    const size_t size_bytes = sizeof(size_t);
    size_t pos = data.size();
    data.resize(pos + size_bytes);
    memcpy(data.data() + pos, &value, size_bytes);
}

// Helper function to write uint64_t to bytes
void write_uint64(std::vector<unsigned char>& data, uint64_t value) {
    const size_t size_bytes = sizeof(uint64_t);
    size_t pos = data.size();
    data.resize(pos + size_bytes);
    memcpy(data.data() + pos, &value, size_bytes);
}

// Helper function to read size_t from bytes
size_t read_size_t(const unsigned char* data, size_t& offset, size_t max_size) {
    if (offset + sizeof(size_t) > max_size) {
        throw std::runtime_error("Buffer overflow while reading size_t");
    }
    
    size_t value;
    memcpy(&value, data + offset, sizeof(size_t));
    offset += sizeof(size_t);
    return value;
}

// Helper function to read uint64_t from bytes
uint64_t read_uint64(const unsigned char* data, size_t& offset, size_t max_size) {
    if (offset + sizeof(uint64_t) > max_size) {
        throw std::runtime_error("Buffer overflow while reading uint64_t");
    }
    
    uint64_t value;
    memcpy(&value, data + offset, sizeof(uint64_t));
    offset += sizeof(uint64_t);
    return value;
}

// Calculate the serialized size of the public key
size_t AccPublicKey::serialized_size() {
    size_t total_size = 0;
    
    // Size of q and count
    total_size += sizeof(uint64_t) + sizeof(size_t);
    
    // Size of param
    total_size += count;
    
    // Size of single elements: g, g_alpha, g_beta, g_gamma, g_delta, g_sq, g_rq
    total_size += element_length_in_bytes(g);
    total_size += element_length_in_bytes(g_alpha);
    total_size += element_length_in_bytes(g_beta);
    total_size += element_length_in_bytes(g_gamma);
    total_size += element_length_in_bytes(g_delta);
    total_size += element_length_in_bytes(g_sq);
    total_size += element_length_in_bytes(g_rq);
    
    // Size of array elements: g_si, g_alpha_si, g_ri, g_beta_ri, g_gamma_ri_sqi
    for (uint64_t i = 0; i < q - 1; i++) {
        total_size += element_length_in_bytes(g_si[i]);
        total_size += element_length_in_bytes(g_alpha_si[i]);
        total_size += element_length_in_bytes(g_ri[i]);
        total_size += element_length_in_bytes(g_beta_ri[i]);
        total_size += element_length_in_bytes(g_gamma_ri_sqi[i]);
    }
    
    // Size of 2D array elements: g_ri_sj, g_delta_ri_sj
    size_t total_pairs = (2 * q - 2) * (2 * q - 2);
    for (size_t idx = 0; idx < total_pairs; idx++) {
        total_size += element_length_in_bytes(g_ri_sj[idx]);
        total_size += element_length_in_bytes(g_delta_ri_sj[idx]);
    }
    
    std::cout << "serialized_size: " << total_size << std::endl;
    return total_size;
}

// Serialize the public key to bytes
std::vector<unsigned char> AccPublicKey::serialize() {
    auto start_time = std::chrono::high_resolution_clock::now();
    
    std::vector<unsigned char> result;
    size_t total_size = serialized_size();
    result.reserve(total_size);  // Reserve memory to avoid reallocations
    
    // Write universe size q
    write_uint64(result, q);
    
    // Write pairing parameters
    write_size_t(result, count);
    size_t pos = result.size();
    result.resize(pos + count);
    memcpy(result.data() + pos, param, count);
    
    // Write single elements
    size_t bytes_written = 0;
    
    // g
    size_t len = element_length_in_bytes(g);
    pos = result.size();
    result.resize(pos + len);
    bytes_written = element_to_bytes(result.data() + pos, g);
    
    // g_alpha
    len = element_length_in_bytes(g_alpha);
    pos = result.size();
    result.resize(pos + len);
    bytes_written = element_to_bytes(result.data() + pos, g_alpha);
    
    // g_beta
    len = element_length_in_bytes(g_beta);
    pos = result.size();
    result.resize(pos + len);
    bytes_written = element_to_bytes(result.data() + pos, g_beta);
    
    // g_gamma
    len = element_length_in_bytes(g_gamma);
    pos = result.size();
    result.resize(pos + len);
    bytes_written = element_to_bytes(result.data() + pos, g_gamma);
    
    // g_delta
    len = element_length_in_bytes(g_delta);
    pos = result.size();
    result.resize(pos + len);
    bytes_written = element_to_bytes(result.data() + pos, g_delta);
    
    // g_sq
    len = element_length_in_bytes(g_sq);
    pos = result.size();
    result.resize(pos + len);
    bytes_written = element_to_bytes(result.data() + pos, g_sq);
    
    // g_rq
    len = element_length_in_bytes(g_rq);
    pos = result.size();
    result.resize(pos + len);
    bytes_written = element_to_bytes(result.data() + pos, g_rq);
    
    // Write array elements
    for (uint64_t i = 0; i < q - 1; i++) {
        // g_si[i]
        len = element_length_in_bytes(g_si[i]);
        pos = result.size();
        result.resize(pos + len);
        bytes_written = element_to_bytes(result.data() + pos, g_si[i]);
        
        // g_alpha_si[i]
        len = element_length_in_bytes(g_alpha_si[i]);
        pos = result.size();
        result.resize(pos + len);
        bytes_written = element_to_bytes(result.data() + pos, g_alpha_si[i]);
        
        // g_ri[i]
        len = element_length_in_bytes(g_ri[i]);
        pos = result.size();
        result.resize(pos + len);
        bytes_written = element_to_bytes(result.data() + pos, g_ri[i]);
        
        // g_beta_ri[i]
        len = element_length_in_bytes(g_beta_ri[i]);
        pos = result.size();
        result.resize(pos + len);
        bytes_written = element_to_bytes(result.data() + pos, g_beta_ri[i]);
        
        // g_gamma_ri_sqi[i]
        len = element_length_in_bytes(g_gamma_ri_sqi[i]);
        pos = result.size();
        result.resize(pos + len);
        bytes_written = element_to_bytes(result.data() + pos, g_gamma_ri_sqi[i]);
    }
    
    // Write 2D array elements
    size_t total_pairs = (2 * q - 2) * (2 * q - 2);
    for (size_t idx = 0; idx < total_pairs; idx++) {
        // g_ri_sj[idx]
        len = element_length_in_bytes(g_ri_sj[idx]);
        pos = result.size();
        result.resize(pos + len);
        bytes_written = element_to_bytes(result.data() + pos, g_ri_sj[idx]);
        
        // g_delta_ri_sj[idx]
        len = element_length_in_bytes(g_delta_ri_sj[idx]);
        pos = result.size();
        result.resize(pos + len);
        bytes_written = element_to_bytes(result.data() + pos, g_delta_ri_sj[idx]);
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    std::cout << "Time for serialization: " << duration.count() / 1000.0 << " seconds" << std::endl;
    
    return result;
}

// Deserialize the public key from bytes
size_t AccPublicKey::deserialize(const unsigned char* data, size_t size) {
    auto start_time = std::chrono::high_resolution_clock::now();
    std::cout << "Starting deserialization. Total size: " << size << std::endl;

    size_t offset = 0;

    // Read universe size q (need to update q if different from current value)
    std::cout << "Reading q. Current offset: " << offset << std::endl;
    uint64_t new_q = read_uint64(data, offset, size);
    std::cout << "Read q = " << new_q << ". New offset: " << offset << std::endl;

    // Read pairing parameters
    std::cout << "Reading param_count. Current offset: " << offset << std::endl;
    size_t param_count = read_size_t(data, offset, size);
    std::cout << "Read param_count = " << param_count << ". New offset: " << offset << std::endl;

    if (offset + param_count > size) {
        throw std::runtime_error("Buffer overflow while reading pairing parameters");
    }

    // If we need to reinitialize with new parameters or q value
    bool reinitialized = false;
    if (new_q != q || param_count != count || memcmp(param, data + offset, param_count) != 0) {
        std::cout << "Reinitializing public key..." << std::endl;
        reinitialized = true;
        // Clean up existing arrays if necessary
        this->~AccPublicKey();

        // Allocate memory for new parameters if needed (assuming param is a fixed-size buffer for now)
        if (param_count > sizeof(this->param)) {
             throw std::runtime_error("Parameter size too large for internal buffer");
        }
        memcpy(param, data + offset, param_count);
        count = param_count;
        q = new_q;

        // Re-initialize pairing and elements
        pairing_init_set_buf(pairing, param, count);

        // Initialize elements
        element_init_G1(g, pairing);
        element_init_G1(g_alpha, pairing);
        element_init_G1(g_beta, pairing);
        element_init_G1(g_gamma, pairing);
        element_init_G1(g_delta, pairing);
        element_init_G1(g_sq, pairing);
        element_init_G1(g_rq, pairing);

        // Initialize arrays
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

        // Initialize 2D arrays
        size_t total_pairs = (2 * q - 2) * (2 * q - 2);
        g_ri_sj = new element_t[total_pairs];
        g_delta_ri_sj = new element_t[total_pairs];

        for (size_t idx = 0; idx < total_pairs; idx++) {
            element_init_G1(g_ri_sj[idx], pairing);
            element_init_G1(g_delta_ri_sj[idx], pairing);
        }
         std::cout << "Reinitialization complete." << std::endl;
    }
    // Important: Advance offset past the param buffer *after* potential memcpy
    // std::cout << "Reading param data. Current offset: " << offset << std::endl;
    offset += param_count;
    // std::cout << "Finished reading param data. New offset: " << offset << std::endl;


    // Read single elements
    int bytes_read = 0;

    bytes_read = element_from_bytes(g, const_cast<unsigned char*>(data + offset));
    if (bytes_read <= 0) throw std::runtime_error("Failed to read element g or zero bytes read");
    offset += bytes_read;

    bytes_read = element_from_bytes(g_alpha, const_cast<unsigned char*>(data + offset));
     if (bytes_read <= 0) throw std::runtime_error("Failed to read element g_alpha or zero bytes read");
    offset += bytes_read;

    bytes_read = element_from_bytes(g_beta, const_cast<unsigned char*>(data + offset));
     if (bytes_read <= 0) throw std::runtime_error("Failed to read element g_beta or zero bytes read");
    offset += bytes_read;

    bytes_read = element_from_bytes(g_gamma, const_cast<unsigned char*>(data + offset));
    if (bytes_read <= 0) throw std::runtime_error("Failed to read element g_gamma or zero bytes read");
    offset += bytes_read;

    bytes_read = element_from_bytes(g_delta, const_cast<unsigned char*>(data + offset));
    if (bytes_read <= 0) throw std::runtime_error("Failed to read element g_delta or zero bytes read");
    offset += bytes_read;

    bytes_read = element_from_bytes(g_sq, const_cast<unsigned char*>(data + offset));
     if (bytes_read <= 0) throw std::runtime_error("Failed to read element g_sq or zero bytes read");
    offset += bytes_read;

    bytes_read = element_from_bytes(g_rq, const_cast<unsigned char*>(data + offset));
    if (bytes_read <= 0) throw std::runtime_error("Failed to read element g_rq or zero bytes read");
    offset += bytes_read;


    // Read array elements
    std::cout << "Reading array elements. Current offset: " << offset << std::endl;
    for (uint64_t i = 0; i < q - 1; i++) {
        bytes_read = element_from_bytes(g_si[i], const_cast<unsigned char*>(data + offset));
        if (bytes_read <= 0) throw std::runtime_error("Failed to read element g_si or zero bytes read");
        offset += bytes_read;

        bytes_read = element_from_bytes(g_alpha_si[i], const_cast<unsigned char*>(data + offset));
        if (bytes_read <= 0) throw std::runtime_error("Failed to read element g_alpha_si or zero bytes read");
        offset += bytes_read;

        bytes_read = element_from_bytes(g_ri[i], const_cast<unsigned char*>(data + offset));
        if (bytes_read <= 0) throw std::runtime_error("Failed to read element g_ri or zero bytes read");
        offset += bytes_read;

        bytes_read = element_from_bytes(g_beta_ri[i], const_cast<unsigned char*>(data + offset));
        if (bytes_read <= 0) throw std::runtime_error("Failed to read element g_beta_ri or zero bytes read");
        offset += bytes_read;

        bytes_read = element_from_bytes(g_gamma_ri_sqi[i], const_cast<unsigned char*>(data + offset));
        if (bytes_read <= 0) throw std::runtime_error("Failed to read element g_gamma_ri_sqi or zero bytes read");
        offset += bytes_read;
    }
    std::cout << "Finished reading array elements. New offset: " << offset << std::endl;

    // Read 2D array elements
    size_t total_pairs = (2 * q - 2) * (2 * q - 2);
    std::cout << "Reading 2D array elements (" << total_pairs << " pairs). Current offset: " << offset << std::endl;
    for (size_t idx = 0; idx < total_pairs; idx++) {
        bytes_read = element_from_bytes(g_ri_sj[idx], const_cast<unsigned char*>(data + offset));
        if (bytes_read <= 0) throw std::runtime_error("Failed to read element g_ri_sj or zero bytes read");
        offset += bytes_read;

        bytes_read = element_from_bytes(g_delta_ri_sj[idx], const_cast<unsigned char*>(data + offset));
        if (bytes_read <= 0) throw std::runtime_error("Failed to read element g_delta_ri_sj or zero bytes read");
        offset += bytes_read;
    }
    // std::cout << "Finished reading 2D array elements. Final offset: " << offset << std::endl;


    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    std::cout << "Time for deserialization: " << duration.count() / 1000.0 << " seconds" << std::endl;

    // Optional: Check if the final offset matches the total size (if known/expected)
     std::cout << "Deserialization finished. Total bytes read: " << offset << std::endl;
    if (offset != size) { // This check might be too strict if size includes padding etc.
        std::cerr << "Warning: Final offset (" << offset << ") does not match input size (" << size << ")" << std::endl;
    }


    return offset;
}

// Deserialize the public key from bytes (Multithreaded Version)
size_t AccPublicKey::deserialize_mul_thread(const unsigned char* data, size_t size) {
    auto start_time = std::chrono::high_resolution_clock::now();
    std::cout << "Starting deserialization (multi-thread). Total size: " << size << std::endl;

    size_t offset = 0;

    // === Sequential Part: Read q, param, reinitialize if needed ===
    std::cout << "Reading q. Current offset: " << offset << std::endl;
    uint64_t new_q = read_uint64(data, offset, size);
    std::cout << "Read q = " << new_q << ". New offset: " << offset << std::endl;

    std::cout << "Reading param_count. Current offset: " << offset << std::endl;
    size_t param_count = read_size_t(data, offset, size);
    std::cout << "Read param_count = " << param_count << ". New offset: " << offset << std::endl;

    if (offset + param_count > size) {
        throw std::runtime_error("Buffer overflow while reading pairing parameters");
    }

    bool reinitialized = false;
    if (new_q != q || param_count != count || memcmp(param, data + offset, param_count) != 0) {
        std::cout << "Reinitializing public key..." << std::endl;
        reinitialized = true;
        // Clean up existing resources first
        this->~AccPublicKey();

        // Allocate memory for new parameters if needed
        if (param_count > sizeof(this->param)) {
             throw std::runtime_error("Parameter size too large for internal buffer");
        }
        memcpy(param, data + offset, param_count);
        count = param_count;
        q = new_q; // Set the new q value

        // Re-initialize pairing and elements
        pairing_init_set_buf(pairing, param, count);

        // Initialize single elements
        element_init_G1(g, pairing);
        element_init_G1(g_alpha, pairing);
        element_init_G1(g_beta, pairing);
        element_init_G1(g_gamma, pairing);
        element_init_G1(g_delta, pairing);
        element_init_G1(g_sq, pairing);
        element_init_G1(g_rq, pairing);

        // Initialize 1D arrays (only if q > 1)
        if (q > 1) {
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
        } else {
             g_si = nullptr; g_alpha_si = nullptr; g_ri = nullptr; g_beta_ri = nullptr; g_gamma_ri_sqi = nullptr;
        }


        // Initialize 2D arrays (only if q > 1)
         if (q > 1) {
            size_t total_pairs = (2 * q - 2) * (2 * q - 2);
            g_ri_sj = new element_t[total_pairs];
            g_delta_ri_sj = new element_t[total_pairs];
            for (size_t idx = 0; idx < total_pairs; idx++) {
                element_init_G1(g_ri_sj[idx], pairing);
                element_init_G1(g_delta_ri_sj[idx], pairing);
            }
         } else {
              g_ri_sj = nullptr; g_delta_ri_sj = nullptr;
         }

         std::cout << "Reinitialization complete." << std::endl;
    }
    // Advance offset past the param buffer
    offset += param_count;
    std::cout << "Finished reading param data. New offset: " << offset << std::endl;

    // === Sequential Part: Read single elements ===
    int bytes_read = 0;
    std::cout << "Reading single elements sequentially. Current offset: " << offset << std::endl;

    bytes_read = element_from_bytes(g, const_cast<unsigned char*>(data + offset));
    if (bytes_read <= 0) throw std::runtime_error("Failed to read element g or zero bytes read");
    offset += bytes_read;

    bytes_read = element_from_bytes(g_alpha, const_cast<unsigned char*>(data + offset));
     if (bytes_read <= 0) throw std::runtime_error("Failed to read element g_alpha or zero bytes read");
    offset += bytes_read;

    bytes_read = element_from_bytes(g_beta, const_cast<unsigned char*>(data + offset));
     if (bytes_read <= 0) throw std::runtime_error("Failed to read element g_beta or zero bytes read");
    offset += bytes_read;

    bytes_read = element_from_bytes(g_gamma, const_cast<unsigned char*>(data + offset));
    if (bytes_read <= 0) throw std::runtime_error("Failed to read element g_gamma or zero bytes read");
    offset += bytes_read;

    bytes_read = element_from_bytes(g_delta, const_cast<unsigned char*>(data + offset));
    if (bytes_read <= 0) throw std::runtime_error("Failed to read element g_delta or zero bytes read");
    offset += bytes_read;

    bytes_read = element_from_bytes(g_sq, const_cast<unsigned char*>(data + offset));
     if (bytes_read <= 0) throw std::runtime_error("Failed to read element g_sq or zero bytes read");
    offset += bytes_read;

    bytes_read = element_from_bytes(g_rq, const_cast<unsigned char*>(data + offset));
    if (bytes_read <= 0) throw std::runtime_error("Failed to read element g_rq or zero bytes read");
    offset += bytes_read;
    std::cout << "Finished reading single elements. New offset: " << offset << std::endl;


    // === Parallel Part - 1D Arrays ===
    if (q > 1) {
        std::cout << "Reading 1D array elements (multithreaded)... Current offset: " << offset << std::endl;
        size_t offset_1d_start = offset; // Store offset before this block

        // Determine fixed element size (using the 'g' we just read)
        size_t element_size_g1 = element_length_in_bytes(g);
        if (element_size_g1 == 0) {
             throw std::runtime_error("Calculated G1 element size is zero. Cannot proceed with parallel read.");
        }
        std::cout << "Determined G1 element size: " << element_size_g1 << " bytes." << std::endl;

        // Calculate total size expected for this block
        const uint64_t total_elements_1d = q - 1;
        size_t expected_1d_block_size = total_elements_1d * 5 * element_size_g1;
        if (offset_1d_start + expected_1d_block_size > size) {
            throw std::runtime_error("Buffer overflow detected before reading 1D arrays.");
        }

        const int num_threads = std::max(1u, std::thread::hardware_concurrency()); // Use at least 1 thread
        const uint64_t elements_per_thread = (total_elements_1d + num_threads - 1) / num_threads; // Ceiling division

        std::vector<std::thread> threads_1d;
        std::vector<std::exception_ptr> thread_exceptions_1d(num_threads, nullptr); // To capture exceptions from threads

        // Lambda to read a chunk of the 1D arrays
        auto read_1d_chunk = [&](uint64_t thread_id, uint64_t start_idx, uint64_t end_idx) {
            try {
                int local_bytes_read = 0;
                for (uint64_t i = start_idx; i < end_idx; ++i) {
                    size_t current_element_base_offset = offset_1d_start + (i * 5 * element_size_g1);

                    // g_si[i]
                    local_bytes_read = element_from_bytes(g_si[i], const_cast<unsigned char*>(data + current_element_base_offset));
                    if (local_bytes_read <= 0 || (size_t)local_bytes_read != element_size_g1) throw std::runtime_error("Read error/size mismatch for g_si[" + std::to_string(i) + "]");

                    // g_alpha_si[i]
                    local_bytes_read = element_from_bytes(g_alpha_si[i], const_cast<unsigned char*>(data + current_element_base_offset + element_size_g1));
                     if (local_bytes_read <= 0 || (size_t)local_bytes_read != element_size_g1) throw std::runtime_error("Read error/size mismatch for g_alpha_si[" + std::to_string(i) + "]");

                    // g_ri[i]
                    local_bytes_read = element_from_bytes(g_ri[i], const_cast<unsigned char*>(data + current_element_base_offset + 2 * element_size_g1));
                     if (local_bytes_read <= 0 || (size_t)local_bytes_read != element_size_g1) throw std::runtime_error("Read error/size mismatch for g_ri[" + std::to_string(i) + "]");

                    // g_beta_ri[i]
                    local_bytes_read = element_from_bytes(g_beta_ri[i], const_cast<unsigned char*>(data + current_element_base_offset + 3 * element_size_g1));
                     if (local_bytes_read <= 0 || (size_t)local_bytes_read != element_size_g1) throw std::runtime_error("Read error/size mismatch for g_beta_ri[" + std::to_string(i) + "]");

                    // g_gamma_ri_sqi[i]
                    local_bytes_read = element_from_bytes(g_gamma_ri_sqi[i], const_cast<unsigned char*>(data + current_element_base_offset + 4 * element_size_g1));
                     if (local_bytes_read <= 0 || (size_t)local_bytes_read != element_size_g1) throw std::runtime_error("Read error/size mismatch for g_gamma_ri_sqi[" + std::to_string(i) + "]");
                }
            } catch (...) {
                thread_exceptions_1d[thread_id] = std::current_exception(); // Capture exception
            }
        };

        // Launch threads
        for (int t = 0; t < num_threads; ++t) {
            uint64_t start = t * elements_per_thread;
            uint64_t end = std::min(start + elements_per_thread, total_elements_1d);
            if (start < end) {
                threads_1d.emplace_back(read_1d_chunk, t, start, end);
            }
        }

        // Join threads and check for exceptions
        for (size_t t = 0; t < threads_1d.size(); ++t) {
             if (threads_1d[t].joinable()) {
                threads_1d[t].join();
                if (thread_exceptions_1d[t]) {
                    std::rethrow_exception(thread_exceptions_1d[t]); // Rethrow the first captured exception
                }
            }
        }


        offset = offset_1d_start + expected_1d_block_size; // Advance offset by the calculated block size
        std::cout << "Finished reading 1D array elements. New offset: " << offset << std::endl;
    } else {
         std::cout << "Skipping 1D array elements (q <= 1)." << std::endl;
    }


    // === Parallel Part - 2D Arrays ===
    if (q > 1) {
         std::cout << "Reading 2D array elements (multithreaded)... Current offset: " << offset << std::endl;
        size_t offset_2d_start = offset;
        const size_t total_pairs = (2 * q - 2) * (2 * q - 2);

        // Reuse element_size_g1 calculated before
        size_t element_size_g1 = element_length_in_bytes(g); // Ensure it's recalculated or stored if needed
         if (element_size_g1 == 0) {
             throw std::runtime_error("Calculated G1 element size is zero before 2D array read.");
         }

        size_t expected_2d_block_size = total_pairs * 2 * element_size_g1;
         if (offset_2d_start + expected_2d_block_size > size) {
             throw std::runtime_error("Buffer overflow detected before reading 2D arrays.");
         }

        const int num_threads = std::max(1u, std::thread::hardware_concurrency());
        const uint64_t pairs_per_thread = (total_pairs + num_threads - 1) / num_threads; // Ceiling division
        std::vector<std::thread> threads_2d;
        std::vector<std::exception_ptr> thread_exceptions_2d(num_threads, nullptr);

        // Lambda to read a chunk of the 2D arrays
        auto read_2d_chunk = [&](uint64_t thread_id, uint64_t start_idx, uint64_t end_idx) {
             try {
                int local_bytes_read = 0;
                for (size_t idx = start_idx; idx < end_idx; ++idx) {
                     size_t current_element_base_offset = offset_2d_start + (idx * 2 * element_size_g1);

                    // g_ri_sj[idx]
                    local_bytes_read = element_from_bytes(g_ri_sj[idx], const_cast<unsigned char*>(data + current_element_base_offset));
                     if (local_bytes_read <= 0 || (size_t)local_bytes_read != element_size_g1) throw std::runtime_error("Read error/size mismatch for g_ri_sj[" + std::to_string(idx) + "]");

                    // g_delta_ri_sj[idx]
                    local_bytes_read = element_from_bytes(g_delta_ri_sj[idx], const_cast<unsigned char*>(data + current_element_base_offset + element_size_g1));
                    if (local_bytes_read <= 0 || (size_t)local_bytes_read != element_size_g1) throw std::runtime_error("Read error/size mismatch for g_delta_ri_sj[" + std::to_string(idx) + "]");
                }
             } catch (...) {
                 thread_exceptions_2d[thread_id] = std::current_exception();
             }
        };

        // Launch threads
         for (int t = 0; t < num_threads; ++t) {
            uint64_t start = t * pairs_per_thread;
            uint64_t end = std::min(start + pairs_per_thread, total_pairs);
            if (start < end) {
                threads_2d.emplace_back(read_2d_chunk, t, start, end);
            }
        }

        // Join threads and check exceptions
        for (size_t t = 0; t < threads_2d.size(); ++t) {
             if (threads_2d[t].joinable()) {
                threads_2d[t].join();
                 if (thread_exceptions_2d[t]) {
                    std::rethrow_exception(thread_exceptions_2d[t]);
                }
            }
        }

        offset = offset_2d_start + expected_2d_block_size; // Advance offset
        std::cout << "Finished reading 2D array elements. Final offset: " << offset << std::endl;
    } else {
         std::cout << "Skipping 2D array elements (q <= 1). Final offset: " << offset << std::endl;
    }


    // === Finalization ===
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    std::cout << "Time for deserialization (multi-thread): " << duration.count() / 1000.0 << " seconds" << std::endl;

    std::cout << "Deserialization finished. Total bytes read: " << offset << std::endl;
    // Optional: Check if the final offset matches the total size (with caveats about padding)
    // if (offset != size) {
    //     std::cerr << "Warning: Final offset (" << offset << ") does not match input size (" << size << ")" << std::endl;
    // }

    return offset; // Return total bytes read
}

// Static method to create a public key from serialized bytes
AccPublicKey AccPublicKey::from_bytes(const unsigned char* data, size_t size) {
    std::cout << "start to deserialize the public key" << std::endl;
    size_t offset = 0;
    
    // Read q value
    uint64_t new_q = read_uint64(data, offset, size);
    
    // Read pairing parameters
    size_t param_count = read_size_t(data, offset, size);
    if (offset + param_count > size) {
        throw std::runtime_error("Buffer overflow while reading pairing parameters");
    }
    
    // Copy parameters to a temp buffer
    char param_copy[1024];
    if (param_count > 1024) {
        throw std::runtime_error("Parameter size too large");
    }
    memcpy(param_copy, data + offset, param_count);
    
    // Create a new public key
    std::cout << "create a new pk key" << std::endl;
    AccPublicKey pk(param_copy, param_count, new_q);
    std::cout << "pk created" << std::endl;
    
    // Deserialize the rest of the data
    std::cout << "deserialize the rest of the data" << std::endl;
    // pk.deserialize(data, size);
    pk.deserialize_mul_thread(data, size);
    std::cout << "data deserialized" << std::endl;
    
    return pk;
}

} // namespace acc
