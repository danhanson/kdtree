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
#include <cmath>
#include <cstdarg>

namespace dhlib {

	template<int K, typename T, int C>
	struct KDNode;

	template<int K, typename T, int C>
	struct ParentNode;

	template<int K, typename T, int C>
	struct LeafNode;

	template<
		int K, typename T,
		int C = 8,
		typename Palloc = std::allocator<ParentNode<K,T,C> >,
		typename Lalloc = std::allocator<LeafNode<K,T,C> >
	>
	class KDTree;

	template<int K>
	using coord = std::array<double,K>;

	template<int K, typename T>
	using entry = std::pair<const coord<K>&, T&>;

	template<int K>
	double distanceSquared(const coord<K>& c1, const coord<K>& c2) {
		double ans = 0;
		for(int i = 0; i < K; i++){
			ans += pow(c1[i] - c2[i], 2);
		}
		return ans;
	}

	template<int K, typename T, int C>
	class KDNode {
	public:
		const bool isLeaf;
		KDNode(bool isLeaf) :
				isLeaf{isLeaf} {
		}
	};

	template<int K, typename T, int C>
	class ParentNode : public KDNode<K,T,C> {
	public:
		double value;
		KDNode<K,T,C> *left, *right;
		ParentNode(double value, KDNode<K, T, C>& left, KDNode<K, T, C>& right) :
				KDNode<K,T,C>(false), value(value), left(&left), right(&right) { }
	};

	template<int K, typename T, int C>
	class LeafNode : public KDNode<K,T,C> {
	public:
		std::array<bool, C> isFilled;
		std::array<coord<K>,C> coords;
		std::array<T,C> values;
		LeafNode() : KDNode<K,T,C>(true), isFilled(), coords(), values() { }

		template<typename InputIterator>
		LeafNode(InputIterator start, InputIterator end);

		class iterator {
			const LeafNode* leaf;
			int i;
		public:
			iterator();
			iterator(const LeafNode& leaf, bool isEnd);
			iterator(const iterator&);
			iterator& operator++();
			iterator operator++(int);
			iterator& operator--();
			iterator operator--(int);

			iterator& operator=(const iterator& iter);
			//iterator& operator=(iterator&& iter);
			bool operator==(const iterator& iter) const;
			bool operator!=(const iterator& iter) const;
			entry<K,T> operator*() const;
		};

		template<typename Palloc, typename Lalloc>
		ParentNode<K,T,C>* split(int dim);
		iterator begin() const;
		iterator end() const;
		void clear();
	};

	template<int K, typename T, int C>
	struct GetLeafResult {
		LeafNode<K,T,C>* leaf;
		KDNode<K,T,C>** childPtr;
		int i;
	};

	template<
		int K, typename T, int C,
		typename Palloc,
		typename Lalloc
	>
	class KDTree {
		using leaf_iter = typename LeafNode<K,T,C>::iterator;
		Lalloc lalloc;
		Palloc palloc;

		// used to iterate over the elements, do not expect them to be ordered
		class StdIter {
			const LeafNode<K,T,C>* leaf;
			std::stack<const ParentNode<K,T,C>*> parents;
			typename LeafNode<K,T,C>::iterator iter;

			const LeafNode<K,T,C>& seekLeft(const KDNode<K,T,C>& parent);
			const LeafNode<K,T,C>& seekRight(const KDNode<K,T,C>& parent);
		public:

			StdIter() : leaf(nullptr) { }
			StdIter(const KDTree& tree);
			StdIter(const StdIter&);
			StdIter(StdIter&&);

			StdIter& operator++();
			StdIter operator++(int);

			StdIter& operator=(const StdIter& iter);
			StdIter& operator=(StdIter&& iter);

			bool operator==(const StdIter& iter) const;
			bool operator!=(const StdIter& iter) const;

			entry<K,T> operator*() const;
		};

		struct Values {
			const KDTree<K,T,C,Palloc,Lalloc>* const tree;
			class iterator : public StdIter {
			public:
				iterator() : StdIter() { }
				iterator(const KDTree& tree) : StdIter(tree) { }
				iterator(const iterator& iter) : StdIter(iter) { }
				iterator(const iterator&& iter) : StdIter(iter) { }
				T& operator*() const {
					return StdIter::operator*().second;
				}
			};
			iterator begin(){
				return iterator(*tree);
			}
			iterator end(){
				return iterator();
			}
		};

		// a box that contains a volume of space covered by this tree
		class BBox {
		public:
			coord<K> lowerLeft; // corner of box with all minimum values for each dimension
			coord<K> upperRight; // corner of box with all maximum values for each dimension
			template<typename Point>
			bool contains(Point p) const;
		};

		class BoxIter {
			std::stack<const ParentNode<K,T,C>*> parents;
			const LeafNode<K,T,C>* leaf;
			leaf_iter iter;
			BBox bbox;
			void seek(const KDNode<K,T,C>&, int dim);
			void nextLeaf();
		public:
			BoxIter();
			BoxIter(const KDTree& tree, const BBox&);
			BoxIter(const BoxIter&);
			BoxIter(BoxIter&&);

			BoxIter& operator++();
			BoxIter operator++(int);

			BoxIter& operator=(const BoxIter& iter);
			BoxIter& operator=(BoxIter&& iter);

			bool operator==(const BoxIter& iter) const;
			bool operator!=(const BoxIter& iter) const;

			entry<K,T> operator*() const;
		};
	protected:
		KDNode<K,T,C>* root;

		// gets the leaf the specified point goes to along with a reference to the pointer that points to the leaf
		template<typename Point>
		GetLeafResult<K,T,C> getLeaf(const Point& point);
	public:
		using iterator = StdIter;
		using bbox_iterator = BoxIter;

		Values values { this };

		KDTree() {
			LeafNode<K,T,C> *leaf = lalloc.allocate(1);
			lalloc.construct(leaf);
			leaf->clear();
			root = leaf;
		}
		~KDTree();

		template<typename Point>
		T& operator[](const Point&);

		template<typename ...Numbers>
		T& operator()(Numbers... dims){
			std::array<double, sizeof...(Numbers)> x = {{ dims... }};
			return (*this)[x];
		}

		template<typename Point>
		int erase(const Point&);

		iterator begin() const {
			StdIter iter(*this);
			return iter;
		}

		iterator end() const {
			StdIter iter { };
			return iter;
		}

		bbox_iterator get(const BBox&) const;
	};
}
#endif /* KDTREE_H_ */
