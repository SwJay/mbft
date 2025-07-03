#ifndef MERKLE_H
#define MERKLE_H

//#include <vector>
//#include <stack>
//#include <tuple>
//#include <cstdint>

#include "mht/hash.h"
#include "mht/transaction.h"
#include <vector>
#include <stdexcept>

namespace mht {

// 新增 MerklePathElement 结构体，包含哈希值和位置信息
struct MerklePathNode {
    uint32_t hash;        // 兄弟节点的哈希值
    bool isRightSibling;  // 表示该节点是否为右兄弟节点

    MerklePathNode(uint32_t h, bool isRight) 
        : hash(h), isRightSibling(isRight) {}

    size_t getSizeADS() const {
        return sizeof(hash) + sizeof(isRightSibling);
    }
};

class MerkleNode
{
public:
	MerkleNode() = default;
	MerkleNode(uint32_t hash, unsigned int nSeed);
	MerkleNode(uint32_t hash, unsigned int nSeed, MerkleNode * pLeftChild, MerkleNode * pRightChild);
	MerkleNode(const vector<uint32_t>& vHashes, unsigned int nSeed);

	~MerkleNode();

	uint32_t Hash(uint32_t hash1, uint32_t hash2) const;

	uint32_t getHash() const;

    size_t getSizeADS() const;

    // 修改 getMerklePath 返回类型为 vector<MerklePathElement>
    vector<MerklePathNode> getMerklePath(int index) const;

    // 修改 VerifyMerklePath 为静态函数，移除 leafIndex 参数
    static bool VerifyMerklePath(uint32_t leafHash, const vector<MerklePathNode>& merklePath, unsigned int nSeed, uint32_t expectedRootHash);

    // 原有的验证方法保留，但设为私有
    // bool VerifyMerklePath(uint32_t leafHash, int leafIndex, const vector<MerkleNode>& merklePath, unsigned int nSeed, uint32_t expectedRootHash);

    MerkleNode* getLeftChild() const;
    MerkleNode* getRightChild() const;
	MerkleNode* getParent() const;

    void setChildrenNull();

    void printBFS() const;

    MerkleNode* getLeaf(int index) const;

private:

    // Property
	uint32_t hash = 0;
	unsigned int nSeed = 0;

	MerkleNode* pLeftChild = nullptr;
	MerkleNode* pRightChild = nullptr;
	MerkleNode* pParent = nullptr;

    std::vector<MerkleNode*> leaves;
};

} // namespace mht

#endif