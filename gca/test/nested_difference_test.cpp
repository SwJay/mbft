/**
 * @file nested_difference_test.cpp
 * @brief Test for nested difference operation in GCA
 */

#include <gtest/gtest.h>
#include "acc/accumulator.h"
#include <vector>
#include <string>

// 定义一个测试固件用于Nested Difference测试
class NestedDifferenceTest : public ::testing::Test {
protected:
    // acc::AccPublicKey pk;
    std::string key_path = "../../data/key/pk.bin"; // 公钥文件路径
    std::string pbc_param_path = "../../data/param/a.param"; // PBC参数文件路径
    uint64_t universe_size = 100; // 设置适当的宇宙大小

    // void SetUp() override {}

    // void TearDown() override {}
};

// 测试基本的嵌套差集操作
TEST_F(NestedDifferenceTest, BasicNestedDifference) {
    // 创建两个测试集合
    acc::AccPublicKey pk = acc::Accumulator::genkey(key_path, pbc_param_path, universe_size);
    
    acc::Set set_a, set_b;
    
    // 向集合A添加元素
    set_a.insert(1);
    set_a.insert(2);
    set_a.insert(3);
    set_a.insert(4);
    set_a.insert(5);
    
    // 向集合B添加元素
    set_b.insert(3);
    set_b.insert(4);
    set_b.insert(5);
    set_b.insert(6);
    
    // 计算累加器
    acc::AccValue acc_a = acc::Accumulator::setup(set_a, pk);
    acc::AccValue acc_b = acc::Accumulator::setup(set_b, pk);
    
    // 计算嵌套差集及其证明
    auto [diff_acc, proof] = acc::Accumulator::query_nested_difference(set_a, set_b, pk);
    
    // 预计的差集结果应该是 {1, 2}
    acc::Set expected_diff;
    expected_diff.insert(1);
    expected_diff.insert(2);
    
    // 验证嵌套差集操作的正确性
    acc::AccValue expected_diff_acc = acc::Accumulator::setup(expected_diff, pk); 
    
    // 验证proof是否有效
    bool is_valid = acc::Accumulator::verify_nested(acc_a, acc_b, proof, diff_acc, pk);
    EXPECT_TRUE(is_valid) << "嵌套差集证明验证失败";
    
    // 验证差集累加器是否与预期相符
    EXPECT_TRUE(diff_acc.equals(expected_diff_acc)) << "嵌套差集结果与预期不符";
}

// 测试空集差集
TEST_F(NestedDifferenceTest, EmptySetDifference) {
    // 测试情况1：A为空集
    acc::AccPublicKey pk = acc::Accumulator::genkey(key_path, pbc_param_path, universe_size);
    acc::Set empty_set, set_b;
    
    // 向集合B添加元素
    set_b.insert(1);
    set_b.insert(2);
    
    // 计算累加器
    acc::AccValue acc_empty = acc::Accumulator::setup(empty_set, pk);
    acc::AccValue acc_b = acc::Accumulator::setup(set_b, pk);
    
    // 计算嵌套差集及其证明 (空集 - B = 空集)
    auto [diff_acc1, proof1] = acc::Accumulator::query_nested_difference(empty_set, set_b, pk);
    
    // 验证证明
    bool is_valid1 = acc::Accumulator::verify_nested(acc_empty, acc_b, proof1, diff_acc1, pk);
    EXPECT_TRUE(is_valid1) << "空集与非空集的嵌套差集证明验证失败";
    EXPECT_TRUE(diff_acc1.equals(acc_empty)) << "空集与非空集的嵌套差集结果应为空集";
    
    // 测试情况2：B为空集
    acc::Set set_a;
    set_a.insert(1);
    set_a.insert(2);
    
    // 计算累加器
    acc::AccValue acc_a = acc::Accumulator::setup(set_a, pk);
    
    // 计算嵌套差集及其证明 (A - 空集 = A)
    auto [diff_acc2, proof2] = acc::Accumulator::query_nested_difference(set_a, empty_set, pk);
    
    // 验证证明
    bool is_valid2 = acc::Accumulator::verify_nested(acc_a, acc_empty, proof2, diff_acc2, pk);
    EXPECT_TRUE(is_valid2) << "非空集与空集的嵌套差集证明验证失败";
    EXPECT_TRUE(diff_acc2.equals(acc_a)) << "非空集与空集的嵌套差集结果应为原集合";
}

// 测试相同集合的差集
TEST_F(NestedDifferenceTest, SameSetDifference) {
    // 创建测试集合
    acc::AccPublicKey pk = acc::Accumulator::genkey(key_path, pbc_param_path, universe_size);
    acc::Set set_a;
    set_a.insert(1);
    set_a.insert(2);
    set_a.insert(3);
    
    // 计算累加器
    acc::AccValue acc_a = acc::Accumulator::setup(set_a, pk);
    
    // 计算嵌套差集及其证明 (A - A = 空集)
    auto [diff_acc, proof] = acc::Accumulator::query_nested_difference(set_a, set_a, pk);
    
    // 预计的差集结果应该是空集
    acc::Set empty_set;
    acc::AccValue empty_acc = acc::Accumulator::setup(empty_set, pk);
    
    // 验证证明
    bool is_valid = acc::Accumulator::verify_nested(acc_a, acc_a, proof, diff_acc, pk);
    EXPECT_TRUE(is_valid) << "相同集合的嵌套差集证明验证失败";
    EXPECT_TRUE(diff_acc.equals(empty_acc)) << "相同集合的嵌套差集结果应为空集";
}

// 测试大集合的差集
TEST_F(NestedDifferenceTest, LargeSetDifference) {
    // 创建两个测试集合
    acc::AccPublicKey pk = acc::Accumulator::genkey(key_path, pbc_param_path, universe_size);
    acc::Set set_a, set_b;
    
    // 向集合A添加多个元素
    for (uint64_t i = 1; i <= 20; i++) {
        set_a.insert(i);
    }
    
    // 向集合B添加部分元素
    for (uint64_t i = 10; i <= 30; i++) {
        set_b.insert(i);
    }
    
    // 计算累加器
    acc::AccValue acc_a = acc::Accumulator::setup(set_a, pk);
    acc::AccValue acc_b = acc::Accumulator::setup(set_b, pk);
    
    // 计算嵌套差集及其证明
    auto [diff_acc, proof] = acc::Accumulator::query_nested_difference(set_a, set_b, pk);
    
    // 预计的差集结果应该是 {1, 2, ..., 9}
    acc::Set expected_diff;
    for (uint64_t i = 1; i <= 9; i++) {
        expected_diff.insert(i);
    }
    
    // 验证嵌套差集操作的正确性
    acc::AccValue expected_diff_acc = acc::Accumulator::setup(expected_diff, pk);
    
    // 验证证明
    bool is_valid = acc::Accumulator::verify_nested(acc_a, acc_b, proof, diff_acc, pk);
    EXPECT_TRUE(is_valid) << "大集合的嵌套差集证明验证失败";
    EXPECT_TRUE(diff_acc.equals(expected_diff_acc)) << "大集合的嵌套差集结果与预期不符";
}

// 测试无交集的两个集合的差集
TEST_F(NestedDifferenceTest, DisjointSetDifference) {
    // 创建两个无交集的测试集合
    acc::AccPublicKey pk = acc::Accumulator::genkey(key_path, pbc_param_path, universe_size);
    acc::Set set_a, set_b;
    
    // 向集合A添加元素
    for (uint64_t i = 1; i <= 5; i++) {
        set_a.insert(i);
    }
    
    // 向集合B添加无交集的元素
    for (uint64_t i = 6; i <= 10; i++) {
        set_b.insert(i);
    }
    
    // 计算累加器
    acc::AccValue acc_a = acc::Accumulator::setup(set_a, pk);
    acc::AccValue acc_b = acc::Accumulator::setup(set_b, pk);
    
    // 计算嵌套差集及其证明
    auto [diff_acc, proof] = acc::Accumulator::query_nested_difference(set_a, set_b, pk);
    
    // 预计的差集结果应该是集合A本身
    
    // 验证嵌套差集操作的正确性
    bool is_valid = acc::Accumulator::verify_nested(acc_a, acc_b, proof, diff_acc, pk);
    EXPECT_TRUE(is_valid) << "无交集集合的嵌套差集证明验证失败";
    EXPECT_TRUE(diff_acc.equals(acc_a)) << "无交集集合的嵌套差集结果应为集合A本身";
}

// 主测试入口
int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
} 