/*
 * kdtree.h
 *
 *  Created on: Aug 2, 2016
 *      Author: hansondg
 */

#ifndef KDTREE_H_
#define KDTREE_H_

#include <functional>
#include <map>
#include <iterator>
#include <stack>

namespace dhlib {

	template<unsigned K, typename T>
	using entry = std::pair<std::array<double,K>, T>;

	template<unsigned K, typename T, typename Ealloc>
	struct KDNode;

	template<unsigned K, typename T, typename Ealloc>
	struct ParentNode;

	template<unsigned K, typename T, typename Ealloc>
	struct LeafNode;

	template<
		unsigned K, typename T,
		typename Ealloc = std::allocator<std::pair<double, entry<K,T> > >,
		typename Palloc = std::allocator<ParentNode<K,T,Ealloc> >,
		typename Lalloc = std::allocator<LeafNode<K,T,Ealloc> >
	>
	class KDTree;

	template<unsigned K, typename T>
	using kd_iter = std::iterator<std::bidirectional_iterator_tag, entry<K,T>, std::ptrdiff_t, entry<K,T>*, entry<K,T>& >;

	template<unsigned K, typename T, typename Ealloc>
	using leaf_map = std::multimap<double, entry<K,T>, std::less<double>, Ealloc>;

	template<unsigned K, typename T, typename Ealloc>
	struct KDNode {
		const bool isLeaf;
		KDNode(bool isLeaf) :
				isLeaf{isLeaf} {
		}
	};

	template<unsigned K, typename T, typename Ealloc>
	struct ParentNode : KDNode<K,T,Ealloc> {
		double value;
		const KDNode<K,T,Ealloc>* left, right;
		ParentNode(double value, KDNode<K, T, Ealloc>& left, KDNode<K, T, Ealloc>& right) :
				KDNode<K, T, Ealloc>(false), value{value}, left(left), right(right) { }
	};

	template<unsigned K, typename T, typename Ealloc>
	struct LeafNode : KDNode<K,T, Ealloc> {
		const leaf_map<K,T,Ealloc> values;
		template<typename InputIterator>
		LeafNode(InputIterator start, InputIterator end) : values(start, end) { }
		LeafNode() : values() { }
	};

	template<unsigned K, typename T, typename Ealloc>
	struct LeafData {
		LeafNode<K,T,Ealloc>& leaf;
		KDNode<K,T,Ealloc>*& childPtr;
		unsigned i;
	};

	template<
		unsigned K, typename T,
		typename Ealloc,
		typename Palloc,
		typename Lalloc
	>
	class KDTree {
	using leaf_iter = typename leaf_map<K,T,Ealloc>::iterator;
	private:

		// used to iterate over the elements in no particular order
		class StdIter : kd_iter<K,T> {
			bool isEnd;
			const KDTree<K,T,Ealloc,Palloc,Lalloc>* const tree;
			const std::stack<const ParentNode<K,T,Ealloc>* const, std::deque<const ParentNode<K,T,Ealloc>,Palloc>> parents;
			const LeafNode<K,T,Ealloc>* leaf;
			leaf_iter iter;

			LeafNode<K,T,Ealloc>& seekLeft(const KDNode<K,T,Ealloc>& parent){
				KDNode<K,T,Ealloc>* node = &parent;
				while(!node->isLeaf){
					auto parent = (ParentNode<K,T,Ealloc>*) node;
					this->parents.push(parent);
					node = parent->left;
				}
				return *((LeafNode<K,T,Ealloc>*) node);
			}

			LeafNode<K,T,Ealloc>& seekRight(const KDNode<K,T,Ealloc>& parent){
				KDNode<K,T,Ealloc>* node = &parent;
				while(!node->isLeaf){
					auto parent = (ParentNode<K,T,Ealloc>*) node;
					this->parents.push(parent);
					node = parent->right;
				}
				return *((LeafNode<K,T,Ealloc>*) node);
			}

		public:

			StdIter() : tree(nullptr), parents(), leaf(nullptr), iter(), isEnd(true) { }
			StdIter(const KDTree& tree, bool isEnd) : tree(tree), parents(), iter() {
				if(isEnd){
					this->isEnd = true;
					return;
				}
				auto leaf = seekLeft(*tree.root);
				this->map = leaf->values;
				this->iter = leaf->values.begin();
				this->isEnd = this->iter == leaf->values.end();
			}

			StdIter(const KDTree& tree) : StdIter(tree, false) { }

			StdIter& operator++(){
				if(isEnd){ // iterator is at end and/or uninitialized
					return *this;
				}
				this->iter++;
				while(this->iter == this->map.end()){
					KDNode<K,T,Ealloc>* child, *parent = this->leaf;
					do {
						if(this->stack.size() == 0){ // reached end
							isEnd = true;
							return *this;
						}
						child = parent;
						parent = this->stack.pop();
					} while(parent->right == child);
					this->leaf = seekLeft(*(parent->right));
					this->iter = this->leaf.values.begin();
				}
				return *this;
			}

			StdIter& operator--(){
				if(isEnd){
					if(this->tree == nullptr){ // iterator is uninitialized
						return *this;
					}
					if(this->leaf == nullptr){ // we need to calculate the last leaf
						this->leaf = seekRight(8(tree->root));
						this->iter = this->leaf.end(); // we decrement iter in the rest of this function
					}
				}
				while(this->iter == this->leaf->map.begin()){
					KDNode<K,T,Ealloc>* child, parent = this->leaf;
					do {
						if(this->stack.size() == 0){ // decremented passed begin
							this->leaf = seekRight(*(tree->root)); // move iterator back to front
							this->iter = this->leaf.begin();
							return *this;
						}
						child = parent;
						parent = this->stack.pop();
					} while(parent->left == child);
					this->leaf = seekRight(*(parent->left));
					this->iter = leaf->values.end();
				}
				this->iter--;
				return *this;
			}

			StdIter& operator=(const StdIter& iter){
				this(iter);
				return *this;
			}

			StdIter& operator=(StdIter&& iter){
				this(forward(iter));
				return *this;
			}

			bool operator==(const StdIter& iter){
				if(this->isEnd && iter.isEnd && this->tree == iter.tree){
					return true;
				}
				return this->iter == iter.iter;
			}

			bool operator!=(const StdIter& iter){
				return *this == iter;
			}

			entry<K,T>& operator*(){
				return *(this->iter);
			}
		};

		KDNode<K,T,Ealloc>* root;
		size_t leafSize;

		// gets the leaf the specified point goes to along with a reference to the pointer that points to the leaf
		template<typename coord>
		LeafData<K,T,Ealloc> getLeaf(const coord& point);
	public:
		using iterator = kd_iter<K,T>;

		KDTree(size_t leafSize) {
			this->leafSize = leafSize;
			KDNode<K,T,Ealloc>* rootPtr = Lalloc::allocate(1);
			Lalloc::construct(rootPtr);
			StdIter x();
			StdIter y(x);
			root = *rootPtr;
		}

		KDTree() : KDTree(32) { }

		template<typename coord>
		void insert(const coord&, T&& elem);

		template<typename coord>
		void insert(const coord&, const T& elem){
			insert(move(copy(elem)));
		}

		template<typename coord>
		unsigned erase(const coord&);

		iterator begin(){
			return new StdIter(*this);
		}

		iterator end(){
			return new StdIter(*this, true);
		}
	};
}
#endif /* KDTREE_H_ */
