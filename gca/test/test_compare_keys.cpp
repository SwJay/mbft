/**
 * @file test_compare_keys.cpp
 * @brief Test to compare public keys from serialized files
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>
#include <chrono>
#include <cmath>

#include "acc/accumulator.h"

// Function to read file into buffer
std::vector<unsigned char> read_file_to_buffer(const std::string& filename) {
    std::ifstream file(filename, std::ios::binary);
    
    if (!file) {
        std::cerr << "无法打开文件: " << filename << std::endl;
        return {};
    }
    
    // 获取文件大小
    file.seekg(0, std::ios::end);
    std::streamsize file_size = file.tellg();
    file.seekg(0, std::ios::beg);
    
    // 分配缓冲区
    std::vector<unsigned char> buffer(file_size);
    
    // 读取文件内容
    if (!file.read(reinterpret_cast<char*>(buffer.data()), file_size)) {
        std::cerr << "读取文件失败: " << filename << std::endl;
        return {};
    }
    
    return buffer;
}

// 比较两个元素是否相等
bool compare_elements(const element_t& e1, const element_t& e2) {
    // 需要使用const_cast因为PBC库的函数接受非const参数
    return element_cmp(const_cast<element_t&>(e1), const_cast<element_t&>(e2)) == 0;
}

int main(int argc, char** argv) {
    auto [pk1, pk2] = acc::Accumulator::genkey_test("../data/key", "../data/param/a.param", 32);

    // // 默认文件路径，可通过命令行参数覆盖
    // std::string pk1_file = "data/key/pk-5-1";
    // std::string pk2_file = "data/key/pk-5-2";
    
    // if (argc > 2) {
    //     pk1_file = argv[1];
    //     pk2_file = argv[2];
    // }
    
    // std::cout << "比较公钥文件：" << std::endl;
    // std::cout << "  文件1: " << pk1_file << std::endl;
    // std::cout << "  文件2: " << pk2_file << std::endl;
    
    // // 读取文件内容
    // auto start_time = std::chrono::high_resolution_clock::now();
    // std::vector<unsigned char> pk1_data = read_file_to_buffer(pk1_file);
    // std::vector<unsigned char> pk2_data = read_file_to_buffer(pk2_file);
    
    // if (pk1_data.empty() || pk2_data.empty()) {
    //     std::cerr << "文件读取失败，退出测试" << std::endl;
    //     return 1;
    // }
    
    // std::cout << "文件大小：" << std::endl;
    // std::cout << "  " << pk1_file << ": " << pk1_data.size() << " 字节" << std::endl;
    // std::cout << "  " << pk2_file << ": " << pk2_data.size() << " 字节" << std::endl;
    
    // // 从二进制数据创建公钥对象
    // std::cout << "从二进制数据创建公钥..." << std::endl;
    
    // // 由于AccPublicKey没有默认构造函数，使用静态方法创建公钥对象
    // acc::AccPublicKey pk1 = acc::AccPublicKey::from_bytes(pk1_data.data(), pk1_data.size());
    // std::cout << "成功加载第一个公钥" << std::endl;
    
    // acc::AccPublicKey pk2 = acc::AccPublicKey::from_bytes(pk2_data.data(), pk2_data.size());
    // std::cout << "成功加载第二个公钥" << std::endl;
    
    // auto end_time = std::chrono::high_resolution_clock::now();
    // auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    // std::cout << "加载公钥耗时: " << duration.count() / 1000.0 << " 秒" << std::endl;
    
    // 比较两个公钥的基本属性
    std::cout << "比较基本属性:" << std::endl;
    std::cout << "  Universe大小 q: " << (pk1.q == pk2.q ? "相同" : "不同") 
              << " (" << pk1.q << " vs " << pk2.q << ")" << std::endl;
    std::cout << "  参数大小: " << (pk1.count == pk2.count ? "相同" : "不同") 
              << " (" << pk1.count << " vs " << pk2.count << ")" << std::endl;
    
    // 比较单个元素
    std::cout << "比较单个元素:" << std::endl;
    
    bool g_equal = compare_elements(pk1.g, pk2.g);
    std::cout << "  g: " << (g_equal ? "相同" : "不同") << std::endl;
    
    bool g_alpha_equal = compare_elements(pk1.g_alpha, pk2.g_alpha);
    std::cout << "  g_alpha: " << (g_alpha_equal ? "相同" : "不同") << std::endl;
    
    bool g_beta_equal = compare_elements(pk1.g_beta, pk2.g_beta);
    std::cout << "  g_beta: " << (g_beta_equal ? "相同" : "不同") << std::endl;
    
    bool g_gamma_equal = compare_elements(pk1.g_gamma, pk2.g_gamma);
    std::cout << "  g_gamma: " << (g_gamma_equal ? "相同" : "不同") << std::endl;
    
    bool g_delta_equal = compare_elements(pk1.g_delta, pk2.g_delta);
    std::cout << "  g_delta: " << (g_delta_equal ? "相同" : "不同") << std::endl;
    
    bool g_sq_equal = compare_elements(pk1.g_sq, pk2.g_sq);
    std::cout << "  g_sq: " << (g_sq_equal ? "相同" : "不同") << std::endl;
    
    bool g_rq_equal = compare_elements(pk1.g_rq, pk2.g_rq);
    std::cout << "  g_rq: " << (g_rq_equal ? "相同" : "不同") << std::endl;
    
    // 比较数组元素（仅检查前几个和随机选择的元素）
    std::cout << "比较数组元素(采样):" << std::endl;
    
    uint64_t sample_size = std::min(static_cast<uint64_t>(5), pk1.q - 1);
    
    int g_si_diff_count = 0;
    for (uint64_t i = 0; i < sample_size; i++) {
        if (!compare_elements(pk1.g_si[i], pk2.g_si[i])) {
            g_si_diff_count++;
        }
    }
    std::cout << "  g_si: " << g_si_diff_count << "/" << sample_size << " 不同" << std::endl;
    
    int g_alpha_si_diff_count = 0;
    for (uint64_t i = 0; i < sample_size; i++) {
        if (!compare_elements(pk1.g_alpha_si[i], pk2.g_alpha_si[i])) {
            g_alpha_si_diff_count++;
        }
    }
    std::cout << "  g_alpha_si: " << g_alpha_si_diff_count << "/" << sample_size << " 不同" << std::endl;
    
    int g_ri_diff_count = 0;
    for (uint64_t i = 0; i < sample_size; i++) {
        if (!compare_elements(pk1.g_ri[i], pk2.g_ri[i])) {
            g_ri_diff_count++;
        }
    }
    std::cout << "  g_ri: " << g_ri_diff_count << "/" << sample_size << " 不同" << std::endl;
    
    int g_beta_ri_diff_count = 0;
    for (uint64_t i = 0; i < sample_size; i++) {
        if (!compare_elements(pk1.g_beta_ri[i], pk2.g_beta_ri[i])) {
            g_beta_ri_diff_count++;
        }
    }
    std::cout << "  g_beta_ri: " << g_beta_ri_diff_count << "/" << sample_size << " 不同" << std::endl;
    
    int g_gamma_ri_sqi_diff_count = 0;
    for (uint64_t i = 0; i < sample_size; i++) {
        if (!compare_elements(pk1.g_gamma_ri_sqi[i], pk2.g_gamma_ri_sqi[i])) {
            g_gamma_ri_sqi_diff_count++;
        }
    }
    std::cout << "  g_gamma_ri_sqi: " << g_gamma_ri_sqi_diff_count << "/" << sample_size << " 不同" << std::endl;
    
    // 比较 2D 数组 (采样)
    size_t total_pairs = (2 * pk1.q - 2) * (2 * pk1.q - 2);
    size_t sample_pairs = std::min(static_cast<size_t>(10), total_pairs);
    
    int g_ri_sj_diff_count = 0;
    int g_delta_ri_sj_diff_count = 0;
    
    for (size_t idx = 0; idx < sample_pairs; idx++) {
        size_t sample_idx = idx * (total_pairs / sample_pairs);
        if (sample_idx >= total_pairs) {
            sample_idx = total_pairs - 1;
        }
        
        if (!compare_elements(pk1.g_ri_sj[sample_idx], pk2.g_ri_sj[sample_idx])) {
            g_ri_sj_diff_count++;
        }
        
        if (!compare_elements(pk1.g_delta_ri_sj[sample_idx], pk2.g_delta_ri_sj[sample_idx])) {
            g_delta_ri_sj_diff_count++;
        }
    }
    
    std::cout << "  g_ri_sj: " << g_ri_sj_diff_count << "/" << sample_pairs << " 不同" << std::endl;
    std::cout << "  g_delta_ri_sj: " << g_delta_ri_sj_diff_count << "/" << sample_pairs << " 不同" << std::endl;
    
    // 总结比较结果
    bool all_equal = g_equal && g_alpha_equal && g_beta_equal && g_gamma_equal && 
                    g_delta_equal && g_sq_equal && g_rq_equal &&
                    g_si_diff_count == 0 && g_alpha_si_diff_count == 0 &&
                    g_ri_diff_count == 0 && g_beta_ri_diff_count == 0 &&
                    g_gamma_ri_sqi_diff_count == 0 &&
                    g_ri_sj_diff_count == 0 && g_delta_ri_sj_diff_count == 0;
    
    std::cout << "总结: 公钥" << (all_equal ? "完全相同" : "存在差异") << std::endl;
    
    return 0;
} 