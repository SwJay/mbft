#include <string>
#include <queue>
#include <sstream>
#include <iostream>
#include <stdexcept> // For std::out_of_range
#include <vector>
#include <assert.h>
#include <algorithm>
#include <cmath> // For std::ceil, std::log2

#include "mht/merkle.h"
#include "mht/hash.h"
#include "mht/transaction.h" // Ensure Transaction definition is available

namespace mht {

// Merkle Node -------------------------------------------------------------

uint32_t MerkleNode::Hash(const uint32_t hash1, const uint32_t hash2) const {
	std::string str1 = std::to_string(hash1);
	std::string str2 = std::to_string(hash2);
	std::string str = str1 + str2;
	std::vector<unsigned char> vec(str.begin(), str.end());
	uint32_t hash;
    MurmurHash3_x86_32(vec.data(), static_cast<int>(vec.size()), nSeed, &hash);
    return hash;
}

// --- Constructors ---

// Basic constructor for leaf/intermediate nodes (used internally mostly)
MerkleNode::MerkleNode(const uint32_t hashIn, const unsigned int nSeedIn):
	hash(hashIn),
	nSeed(nSeedIn),
	pLeftChild(nullptr),
	pRightChild(nullptr),
	pParent(nullptr) {}

// Basic constructor for intermediate nodes setting children (used internally)
MerkleNode::MerkleNode(const uint32_t hashIn, const unsigned int nSeedIn, MerkleNode * pLeftChildIn, MerkleNode * pRightChildIn):
	hash(hashIn),
	nSeed(nSeedIn),
	pLeftChild(pLeftChildIn),
	pRightChild(pRightChildIn),
	pParent(nullptr){
        // Set the parent pointer for the children
		if (pLeftChild) {
            pLeftChild->pParent = this;
        }
		if (pRightChild) {
		    pRightChild->pParent = this;
        }
	}

// Constructor to build the tree from hashes
MerkleNode::MerkleNode(const std::vector<uint32_t>& vHashes, const unsigned int nSeedIn):
    hash(0),
    nSeed(nSeedIn),
    pLeftChild(nullptr),
    pRightChild(nullptr),
    pParent(nullptr)
{
    if (vHashes.empty()) {
        return;
    }

    // 1. Create leaf nodes
    std::vector<MerkleNode*> currentLevelNodes;
    currentLevelNodes.reserve(vHashes.size());
    leaves.clear();
    leaves.reserve(vHashes.size());

    for (const uint32_t h : vHashes) {
        MerkleNode* leafNode = new MerkleNode(h, nSeedIn);
        currentLevelNodes.push_back(leafNode);
        leaves.push_back(leafNode);
    }

    if (currentLevelNodes.size() == 1) {
        this->hash = currentLevelNodes[0]->hash;
        this->pLeftChild = currentLevelNodes[0];
        currentLevelNodes[0]->pParent = this;
        return;
    }

    // 2. Build internal nodes level by level
    while (currentLevelNodes.size() > 1) {
        std::vector<MerkleNode*> nextLevelNodes;
        nextLevelNodes.reserve((currentLevelNodes.size() + 1) / 2);

        for (size_t i = 0; i < currentLevelNodes.size(); i += 2) {
            MerkleNode* leftChild = currentLevelNodes[i];
            // 如果是最后一个节点且没有右兄弟，则复制自身作为右兄弟
            MerkleNode* rightChild = (i + 1 < currentLevelNodes.size()) ? currentLevelNodes[i + 1] : leftChild;
            
            // 计算父节点的哈希值
            uint32_t parentHash = Hash(leftChild->getHash(), rightChild->getHash());
            MerkleNode* parentNode = new MerkleNode(parentHash, nSeedIn);
            
            // 设置父子关系
            parentNode->pLeftChild = leftChild;
            parentNode->pRightChild = rightChild;
            leftChild->pParent = parentNode;
            if (rightChild != leftChild) {
                rightChild->pParent = parentNode;
            }
            
            nextLevelNodes.push_back(parentNode);
        }
        
        currentLevelNodes = std::move(nextLevelNodes);
    }

    // 3. Set root properties
    MerkleNode* root = currentLevelNodes[0];
    this->hash = root->hash;
    this->pLeftChild = root->pLeftChild;
    this->pRightChild = root->pRightChild;
    
    // 更新子节点的父指针
    if (this->pLeftChild) {
        this->pLeftChild->pParent = this;
    }
    if (this->pRightChild && this->pRightChild != this->pLeftChild) {
        this->pRightChild->pParent = this;
    }

    // 清理临时根节点
    root->pLeftChild = nullptr;
    root->pRightChild = nullptr;
    delete root;
}


// --- Destructor ---
MerkleNode::~MerkleNode()
{
    // Recursively delete children.
    // The check prevents double deletion if a node was duplicated for an odd level.
	if (pLeftChild != nullptr) {
		delete pLeftChild;
    }
	if (pRightChild != nullptr && pRightChild != pLeftChild) { // Check they are not the same pointer
		delete pRightChild;
    }
    // The 'leaves' vector contains pointers to nodes that are part of the tree structure
    // managed by pLeftChild/pRightChild pointers. Deleting them here would be a double delete.
    // Just clear the vector itself.
    leaves.clear();
}

// --- Accessors ---

uint32_t MerkleNode::getHash() const {
	return hash;
}

MerkleNode* MerkleNode::getLeftChild() const {
	return pLeftChild;
}

MerkleNode* MerkleNode::getRightChild() const {
	return pRightChild;
}

MerkleNode* MerkleNode::getParent() const {
	return pParent;
}

MerkleNode* MerkleNode::getLeaf(int index) const {
    if (index < 0 || index >= leaves.size()) {
         throw std::out_of_range("Leaf index out of range");
    }
    return leaves[index];
}

// --- Utility Functions ---

void MerkleNode::setChildrenNull() {
	pLeftChild = nullptr;
	pRightChild = nullptr;
}

size_t MerkleNode::getSizeADS() const {
    size_t sizeADS = sizeof(hash) + sizeof(nSeed);

    if (pLeftChild) sizeADS += pLeftChild->getSizeADS();
    if (pRightChild) sizeADS += pRightChild->getSizeADS();
    
    return sizeADS;
}

void MerkleNode::printBFS() const
{
    if (!this) return; // Handle null pointer case

    std::queue<const MerkleNode*> q; // Use const MerkleNode*
	q.push(this); // push root node

	int levelNodes = 0;
    int level = 0;

    std::cout << "Level 0: ";
	while (!q.empty()) {
        levelNodes = q.size();
        if (levelNodes == 0) break;

        while(levelNodes > 0) {
            const MerkleNode* current = q.front();
            q.pop();
            std::cout << current->getHash() << " ";

            if (current->pLeftChild != nullptr) {
                q.push(current->pLeftChild);
            }
            // Only push right child if it's different from left to avoid duplicates in print
            if (current->pRightChild != nullptr && current->pRightChild != current->pLeftChild) {
                q.push(current->pRightChild);
            }
             // If left exists but right doesn't (shouldn't happen in full tree except maybe root of single node tree)
             // Or if right is same as left (duplicated node case) - don't push again.

            levelNodes--;
        }
        level++;
        if (!q.empty()) { // Don't print "Level X:" if the queue is now empty
            std::cout << std::endl << "Level " << level << ": ";
        }
	}
    std::cout << std::endl;
}


// --- Merkle Path Logic ---

std::vector<MerklePathNode> MerkleNode::getMerklePath(int leafIndex) const {
    std::vector<MerklePathNode> path;

    // 检查叶子节点索引是否有效
    if (leafIndex < 0 || static_cast<size_t>(leafIndex) >= leaves.size()) {
        std::cout << "Leaf index out of range: " << leafIndex << " (valid range: 0 to " << leaves.size() - 1 << ")" << std::endl;
        return path;
    }

    // 获取目标叶子节点
    MerkleNode* targetNode = leaves[leafIndex];
    if (!targetNode) {
        std::cout << "Target node is nullptr at index " << leafIndex << std::endl;
        return path;
    }

    // 如果树只有一个节点，返回空路径
    if (leaves.size() == 1) {
        std::cout << "Single leaf tree, returning empty path" << std::endl;
        return path;
    }

    // 从叶子节点开始向上遍历
    MerkleNode* currentNode = targetNode;
    int currentIndex = leafIndex;

    while (currentNode && currentNode->pParent) {
        MerkleNode* parent = currentNode->pParent;
        
        if (parent->pLeftChild == currentNode) {
            // 当前节点是左子节点，添加右兄弟节点
            // 如果右兄弟不存在或者是自己（奇数个节点的情况），则使用左节点
            MerkleNode* siblingNode = parent->pRightChild ? parent->pRightChild : parent->pLeftChild;
            // 标记为右兄弟
            path.push_back(MerklePathNode(siblingNode->getHash(), true));
        } else {
            // 当前节点是右子节点，添加左兄弟节点
            path.push_back(MerklePathNode(parent->pLeftChild->getHash(), false));
        }
        
        currentNode = parent;
    }

    return path;
}

// bool MerkleNode::VerifyMerklePath(uint32_t leafHash, int leafIndex, const std::vector<MerkleNode>& merklePath, unsigned int nSeed, uint32_t expectedRootHash) {
//     // 检查种子是否匹配
//     if (nSeed != this->nSeed) {
//         std::cout << "Seed mismatch: expected " << this->nSeed << ", got " << nSeed << std::endl;
//         return false;
//     }

//     // 检查叶子节点索引是否有效
//     if (leafIndex < 0 || static_cast<size_t>(leafIndex) >= leaves.size()) {
//         std::cout << "Leaf index out of range: " << leafIndex << " (valid range: 0 to " << leaves.size() - 1 << ")" << std::endl;
//         return false;
//     }

//     // 如果树只有一个节点，路径应该为空
//     if (leaves.size() == 1) {
//         bool result = leafHash == expectedRootHash && merklePath.empty();
//         if (!result) {
//             std::cout << "Single leaf verification failed: leafHash=" << leafHash 
//                       << ", expectedRootHash=" << expectedRootHash 
//                       << ", path.empty()=" << merklePath.empty() << std::endl;
//         }
//         return result;
//     }

//     uint32_t currentHash = leafHash;
//     int currentIndex = leafIndex;

//     for (const MerkleNode& sibling : merklePath) {
//         if (currentIndex % 2 == 0) {
//             // 当前节点是左子节点，与右兄弟节点计算哈希
//             currentHash = Hash(currentHash, sibling.getHash());
//         } else {
//             // 当前节点是右子节点，与左兄弟节点计算哈希
//             currentHash = Hash(sibling.getHash(), currentHash);
//         }
//         currentIndex /= 2; // 移动到父节点的索引
//     }

//     bool result = currentHash == expectedRootHash;
//     if (!result) {
//         std::cout << "Hash verification failed: currentHash=" << currentHash 
//                   << ", expectedRootHash=" << expectedRootHash << std::endl;
//     }
//     return result;
// }

// 新的静态验证方法，不需要知道树的大小和叶子索引
bool MerkleNode::VerifyMerklePath(uint32_t leafHash, const std::vector<MerklePathNode>& merklePath, unsigned int nSeed, uint32_t expectedRootHash) {
    // 如果路径为空，直接比较叶子哈希和根哈希
    if (merklePath.empty()) {
        return leafHash == expectedRootHash;
    }

    // 计算从叶子到根的哈希路径
    uint32_t currentHash = leafHash;

    for (const MerklePathNode& element : merklePath) {
        if (element.isRightSibling) {
            // 兄弟节点在右边，当前节点在左边
            currentHash = MerkleNode(0, nSeed).Hash(currentHash, element.hash);
        } else {
            // 兄弟节点在左边，当前节点在右边
            currentHash = MerkleNode(0, nSeed).Hash(element.hash, currentHash);
        }
    }

    return currentHash == expectedRootHash;
}

} // namespace mht
