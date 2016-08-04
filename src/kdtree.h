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
#include <vector>

namespace dhlib {

	template<unsigned K, typename T, unsigned C>
	struct KDNode;

	template<unsigned K, typename T, unsigned C>
	struct ParentNode;

	template<unsigned K, typename T, unsigned C>
	struct LeafNode;

	template<
		unsigned K, typename T,
		unsigned C = 8,
		typename Palloc = std::allocator<ParentNode<K,T,C> >,
		typename Lalloc = std::allocator<LeafNode<K,T,C> >
	>
	class KDTree;

	template<unsigned K>
	using coord = std::array<double,K>;

	template<unsigned K, typename T>
	using entry = std::pair<const coord<K>&, T&>;

	template<unsigned K, typename T>
	using kd_iter = std::iterator<std::bidirectional_iterator_tag, entry<K,T>, std::ptrdiff_t, entry<K,T>*, entry<K,T>& >;

	template<unsigned K, typename T, unsigned C>
	struct KDNode {
		const bool isLeaf;
		KDNode(bool isLeaf) :
				isLeaf{isLeaf} {
		}
	};

	template<unsigned K, typename T, unsigned C>
	struct ParentNode : KDNode<K,T,C> {
		double value;
		const KDNode<K,T,C>* left, right;
		ParentNode(double value, KDNode<K, T, C>& left, KDNode<K, T, C>& right) :
				KDNode<K, T, C>(false), value(value), left(left), right(right) { }
	};

	template<unsigned K, typename T, unsigned C>
	struct LeafNode : KDNode<K,T,C> {
		const std::array<bool, C> isFilled;
		const std::array<coord<K>,C> coords;
		const std::array<T,C> values;
		template<typename InputIterator>
		LeafNode(InputIterator start, InputIterator end) : values(start, end) { }
		LeafNode() : values() { }
		class iterator : kd_iter<K,T> {
		private:
			unsigned i;
			const LeafNode* leaf;
		public:
			iterator();
			iterator(const LeafNode& leaf);
			iterator& operator++();
			iterator& operator--();
			iterator& operator=(const iterator& iter);
			iterator& operator=(iterator&& iter);
			bool operator==(const iterator& iter);
			bool operator!=(const iterator& iter);
			entry<K,T> operator*();

		};
		iterator begin();
		iterator end();
	};

	template<unsigned K, typename T, unsigned C>
	struct GetLeafResult {
		LeafNode<K,T,C>& leaf;
		KDNode<K,T,C>*& childPtr;
		unsigned i;
	};

	template<
		unsigned K, typename T, unsigned C,
		typename Palloc,
		typename Lalloc
	>
	class KDTree {
	using leaf_iter = typename LeafNode<K,T,C>::iterator;
	private:

		// used to iterate over the elements in no particular order
		class StdIter : kd_iter<K,T> {
			bool isEnd;
			const KDTree<K,T,C,Palloc,Lalloc>* const tree;
			const std::stack<const ParentNode<K,T,C>* const, std::deque<const ParentNode<K,T,C>,Palloc>> parents;
			const LeafNode<K,T,C>* leaf;
			leaf_iter iter;

			LeafNode<K,T,C>& seekLeft(const KDNode<K,T,C>& parent);
			LeafNode<K,T,C>& seekRight(const KDNode<K,T,C>& parent);
		public:

			StdIter() : tree(nullptr), parents(), leaf(nullptr), isEnd(true) { }
			StdIter(const KDTree& tree, bool isEnd);
			StdIter(const KDTree& tree) : StdIter(tree, false) { }

			StdIter& operator++();
			StdIter& operator--();

			StdIter& operator=(const StdIter& iter);
			StdIter& operator=(StdIter&& iter);

			bool operator==(const StdIter& iter);
			bool operator!=(const StdIter& iter);

			entry<K,T> operator*();
		};
	protected:
		KDNode<K,T,C>* root;

		// gets the leaf the specified point goes to along with a reference to the pointer that points to the leaf
		template<typename coord>
		GetLeafResult<K,T,C> getLeaf(const coord& point);
	public:
		using iterator = kd_iter<K,T>;

		KDTree() {
			root = Lalloc::allocate(1);
			Lalloc::construct(root);
		}

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
