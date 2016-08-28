//============================================================================
// Name        : kdtree.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================
#include "kdtree.h"

using namespace std;
using namespace dhlib;

template<int Y, typename Int>
inline Int inc(Int& x){
	return x = (x == Y - 1) ? 0 : x + 1;
}

template<int Y, typename Int>
inline Int dec(Int& x){
	return x = (x == 0) ? Y - 1 : x - 1;
}

template<int K, typename T, int C, typename Palloc, typename Lalloc>
KDTree<K,T,C,Palloc,Lalloc>::~KDTree(){
	stack<ParentNode<K,T,C>*> parents;
	KDNode<K,T,C>* current = root;

// since we need the parents to identify their children, we remove the children first before removing
// the parents by preorder traversal
SEEK:;
	while(!current->isLeaf){ // go to left most leaf below current
		ParentNode<K,T,C>* parent = (ParentNode<K,T,C>*) current;
		parents.push(parent);
		current = parent->left;
	}
	LeafNode<K,T,C>* leaf = (LeafNode<K,T,C>*) current;
	lalloc.deallocate(leaf, 1);
	const KDNode<K,T,C> *deleted = current;
	while(parents.size() > 0){ // back up from the right, deleting any parents whose right child is deleted
		ParentNode<K,T,C>* parent = parents.top();
		parents.pop();
		if(deleted == parent->right){
			palloc.deallocate(parent, 1);
			deleted = parent;
		} else {
			current = parent->right; // we deleted this parent's left children, now we delete its right children
			goto SEEK; // deal with it
		}
	}
}

template<int K, typename T,int C,typename Palloc,typename Lalloc>
template<typename Point>
GetLeafResult<K,T,C> KDTree<K,T,C,Palloc,Lalloc>::getLeaf(const Point& point){
	int i = 0;
	KDNode<K,T,C> **childPtr = &root, *current = root;
	while(!current->isLeaf){
		double eVal = point[i++];
		ParentNode<K,T,C>* parent = (ParentNode<K,T,C>*) current;
		if(eVal > parent->value){
			current = parent->right;
			childPtr = &parent->right;
		} else {
			current = parent->left;
			childPtr = &parent->left;
		}
		if(i == K){
			i = 0;
		}
	}
	LeafNode<K,T,C>* leaf = (LeafNode<K,T,C>*) current;
	GetLeafResult<K,T,C> ret {leaf, childPtr, i};
	return ret;
}


template<int K,typename T,int C,typename Palloc,typename Lalloc>
template<typename Point>
T& KDTree<K,T,C,Palloc,Lalloc>::operator[](const Point& point){
	GetLeafResult<K,T,C> leafData = getLeaf(point);
	LeafNode<K,T,C> *leaf = leafData.leaf;
	KDNode<K,T,C>** childPtr = leafData.childPtr;

	int i = 0;
	while(leaf->isFilled[i]){
		if(i == C){ // if leaf is full, split leaf
			ParentNode<K,T,C> *parent = leaf->template split<Palloc,Lalloc>(leafData.i);
			if(point[leafData.i] > parent->value){
				leaf = (LeafNode<K,T,C>*) parent->right;
				i = (C + 1) / 2;
			} else {
				leaf = (LeafNode<K,T,C>*) parent->left;
				i = C / 2;
			}
			*childPtr = parent; // move parent's pointer from deleted leaf to new parent
		} else {
			++i;
		}
	}
	leaf->isFilled[i] = true;
	for(int n = 0; n < K; n++){ // copy coordinates into tree
		leaf->coords[i][n] = point[n];
	}
	return leaf->values[i];
}



template<int K, typename T,int C,typename Palloc,typename Lalloc>
template<typename Point>
int KDTree<K,T,C,Palloc,Lalloc>::erase(const Point& point){
	GetLeafResult<K,T,C> leafData = getLeaf(point);
	LeafNode<K,T,C>& leaf = leafData.leaf;
	KDNode<K,T,C>& childPtr = leafData.childPtr;
	int count = 0;
	for(int i = 0; i < C; i++){
		if(leaf.isFilled[i]){
			for(int j = 0; j < K; j++){
				if(leaf.coords[i][j] == point[j]){
					count++;
					leaf.isFilled[i] = false;
				}
			}
		}
	}
	return count;
}

template<int K, typename T, int C>
template<typename InputIterator>
LeafNode<K,T,C>::LeafNode(InputIterator start, InputIterator end) : LeafNode<K,T,C>() {
	auto filledIter = isFilled.begin();
	auto coordsIter = coords.begin();
	auto valuesIter = values.begin();
	while(start != end){
		*filledIter = true;
		*coordsIter = *start.first;
		*valuesIter = *start.second;
		++start; ++filledIter; ++coordsIter; ++valuesIter;
	}
	while(filledIter != isFilled.end()){
		*(filledIter++) = false;
	}
}

template<int K, typename T, int C>
void LeafNode<K,T,C>::clear() {
	fill(isFilled.begin(), isFilled.end(), false);
}

template<int K, typename T, int C>
template<typename Palloc, typename Lalloc>
ParentNode<K,T,C>* LeafNode<K,T,C>::split(int dim) {
	Palloc palloc;
	Lalloc lalloc;

	// quick search to find median element
	double pivot;
	int min = 0;
	int max = K-1;
	while(min < max){
		int above = max;
		int below = min;
		pivot = coords[min][dim];
		while(above > below){
			while(coords[below][dim] <= pivot){
				below++;
			}
			while(coords[above][dim] >= pivot){
				above--;
			}
			swap(coords[above], coords[below]);
			swap(values[above], values[below]);
		}
		if(below < K / 2){ // we will choose the left median when there is even number of values
			min = below;
		} else {
			max = below;
		}
	}

	// create replacement node
	LeafNode<K,T,C>* newLeft = lalloc.allocate(1);
	lalloc.construct(newLeft);

	LeafNode<K,T,C>* newRight = lalloc.allocate(1);
	lalloc.construct(newRight);

	ParentNode<K,T,C>* parent = palloc.allocate(1);
	palloc.construct(parent, coords[min][dim], *newLeft, *newRight); // copy the right elements into the right child

	// partitioning of the values using the median
	coord<K>* left = newLeft->coords.begin();
	coord<K>* right = newRight->coords.begin();
	for(coord<K>* iter = coords.begin(); iter != coords.end(); ++iter){
		if(iter[0][dim] > parent->value){
			*right = *iter;
			++right;
		} else {
			*left = *iter;
			++left;
		}
	}

	fill(newLeft->isFilled.begin(), newLeft->isFilled.begin() + C/2, true);
	fill(newLeft->isFilled.begin() + C/2, newLeft->isFilled.end(), false);

	fill(newRight->isFilled.begin(), newRight->isFilled.begin() + (C+1)/2, true);
	fill(newRight->isFilled.begin() + (C+1)/2, newRight->isFilled.end(), false);
	lalloc.deallocate(this,1);
	return parent;
}

template<int K, typename T, int C>
typename LeafNode<K,T,C>::iterator LeafNode<K,T,C>::begin() const {
	return iterator(*this,false);
}

template<int K, typename T, int C>
typename LeafNode<K,T,C>::iterator LeafNode<K,T,C>::end() const {
	return iterator(*this,true);
}

template<int K, typename T, int C>
LeafNode<K,T,C>::iterator::iterator() : leaf(nullptr), i() { }

template<int K, typename T, int C>
LeafNode<K,T,C>::iterator::iterator(const LeafNode& leaf, bool isEnd) : leaf(&leaf), i(0) {
	if(isEnd){
		i = C;
	} else {
		i = 0;
		while(i < C && !leaf.isFilled[i]){
			++i;
		}
	}
}

template<int K, typename T, int C>
LeafNode<K,T,C>::iterator::iterator(const iterator& it) : leaf(it.leaf), i(it.i) { }

template<int K, typename T, int C>
typename LeafNode<K,T,C>::iterator& LeafNode<K,T,C>::iterator::operator++(){
	do {
		i++;
	} while(i < C && !leaf->isFilled[i]);
	return *this;
}

template<int K, typename T, int C>
typename LeafNode<K,T,C>::iterator LeafNode<K,T,C>::iterator::operator++(int){
	iterator it(*this);
	++(*this);
	return it;
}

template<int K, typename T, int C>
typename LeafNode<K,T,C>::iterator& LeafNode<K,T,C>::iterator::operator--(){
	do {
		i--;
	} while(i > -1 && !leaf->isFilled(i));
	return *this;
}

template<int K, typename T, int C>
typename LeafNode<K,T,C>::iterator LeafNode<K,T,C>::iterator::operator--(int){
	iterator ret (*this);
	--(*this);
	return ret;
}

template<int K, typename T, int C>
typename LeafNode<K,T,C>::iterator& LeafNode<K,T,C>::iterator::operator=(const iterator& iter){
	leaf = iter.leaf;
	i = iter.i;
	return *this;
}

template<int K, typename T, int C>
entry<K,T> LeafNode<K,T,C>::iterator::operator*() const {
	const coord<K>& point = leaf->coords[i];
	T& value = const_cast<T&>(leaf->values[i]); // XXX WHY IS leaf->values[i] CONST?
	pair<const coord<K>*, T*> pair (&point, &value);
	return pair;
}

template<int K, typename T, int C>
bool LeafNode<K,T,C>::iterator::operator==(const iterator& other) const {
	return leaf == other.leaf && i == other.i;
}


template<int K, typename T, int C>
bool LeafNode<K,T,C>::iterator::operator!=(const iterator& other) const {
	return !(*this == other);
}

template<int K, typename T,int C,typename Palloc,typename Lalloc>
KDTree<K,T,C,Palloc,Lalloc>::StdIter::StdIter(const KDTree& tree) {
	leaf = &seekLeft(*tree.root);
	iter = leaf->begin();
	if(iter == leaf->end()){
		leaf = nullptr;
	}
}

template<int K, typename T,int C,typename Palloc,typename Lalloc>
KDTree<K,T,C,Palloc,Lalloc>::StdIter::StdIter(const StdIter& iter)
		: leaf(iter.leaf), parents(iter.parents), iter(iter.iter) { }

template<int K, typename T,int C,typename Palloc,typename Lalloc>
KDTree<K,T,C,Palloc,Lalloc>::StdIter::StdIter(StdIter&& iter)
		: leaf(iter.leaf), parents(move(iter.parents)), iter(move(iter.iter)) { }

template<int K, typename T,int C,typename Palloc,typename Lalloc>
const LeafNode<K,T,C>& KDTree<K,T,C,Palloc,Lalloc>::StdIter::seekLeft(const KDNode<K,T,C>& parent){
	const KDNode<K,T,C>* node = &parent;
	while(!node->isLeaf){
		auto parent = (ParentNode<K,T,C>*) node;
		this->parents.push(parent);
		node = parent->left;
	}
	return *((LeafNode<K,T,C>*) node);
}

template<int K, typename T,int C,typename Palloc,typename Lalloc>
const LeafNode<K,T,C>& KDTree<K,T,C,Palloc,Lalloc>::StdIter::seekRight(const KDNode<K,T,C>& parent){
	const KDNode<K,T,C>* node = &parent;
	while(!node->isLeaf){
		ParentNode<K,T,C>* parent = (ParentNode<K,T,C>*) node;
		this->parents.push(parent);
		node = parent->right;
	}
	return *((LeafNode<K,T,C>*) node);
}

template<int K, typename T,int C,typename Palloc,typename Lalloc>
typename KDTree<K,T,C,Palloc,Lalloc>::StdIter& KDTree<K,T,C,Palloc,Lalloc>::StdIter::operator++(){
	if(this->leaf == nullptr){ // iterator is at end and/or uninitialized
		return *this;
	}
	++this->iter;
	while(this->iter == this->leaf->end()){
		const KDNode<K,T,C> *child = this->leaf;
		while(parents.size() > 0){ // back track to parent with right child
			const ParentNode<K,T,C> *parent = parents.top();
			if(parent->left == child){
				this->leaf = &seekLeft(*(parent->right)); // get left most child of right child
				this->iter = this->leaf->begin();
				goto CHECK_ITER;
			}
			child = parent;
		}
		// backtracked up to root, no nodes left, set leaf to null to indicate that this is an end
		this->leaf = nullptr;
		break;
	CHECK_ITER:;
	}
	return *this;
}

template<int K, typename T,int C,typename Palloc,typename Lalloc>
typename KDTree<K,T,C,Palloc,Lalloc>::StdIter KDTree<K,T,C,Palloc,Lalloc>::StdIter::operator++(int){
	StdIter ret (*this);
	++(*this);
	return ret;
}

template<int K, typename T,int C,typename Palloc,typename Lalloc>
typename KDTree<K,T,C,Palloc,Lalloc>::StdIter& KDTree<K,T,C,Palloc,Lalloc>::StdIter::operator=(const StdIter& iter){
	leaf = iter.leaf;
	parents = iter.parents;
	this->iter = iter.iter;
	return *this;
}

template<int K, typename T,int C,typename Palloc,typename Lalloc>
typename KDTree<K,T,C,Palloc,Lalloc>::StdIter& KDTree<K,T,C,Palloc,Lalloc>::StdIter::operator=(StdIter&& iter){
	leaf = iter.leaf;
	parents = move(iter.parents);
	this->iter = move(iter.iter);
	return *this;
}

template<int K, typename T,int C,typename Palloc,typename Lalloc>
bool KDTree<K,T,C,Palloc,Lalloc>::StdIter::operator==(const StdIter& iter) const {
	if(this->leaf == nullptr && iter.leaf == nullptr){
		return true;
	}
	return this->iter == iter.iter;
}

template<int K, typename T,int C,typename Palloc,typename Lalloc>
bool KDTree<K,T,C,Palloc,Lalloc>::StdIter::operator!=(const StdIter& iter) const {
	return !(*this == iter);
}

template<int K, typename T,int C,typename Palloc,typename Lalloc>
entry<K,T> KDTree<K,T,C,Palloc,Lalloc>::StdIter::operator*() const {
	return *iter;
}

template<int K, typename T,int C,typename Palloc,typename Lalloc>
typename KDTree<K,T,C,Palloc,Lalloc>::bbox_iterator KDTree<K,T,C,Palloc,Lalloc>::get(const BBox& bbox) const {
	return BoxIter(*this, bbox);
}

template<int K, typename T, int C, typename Palloc, typename Lalloc>
template<typename Point>
bool KDTree<K,T,C,Palloc,Lalloc>::BBox::contains(Point p) const {
	for(int i = 0; i < K; i++){
		if(p[i] > upperRight[i] || p[i] < lowerLeft[i]){
			return false;
		}
	}
	return true;
}

template<int K, typename T,int C,typename Palloc,typename Lalloc>
KDTree<K,T,C,Palloc,Lalloc>::BoxIter::BoxIter() : leaf(nullptr), bbox() { }

template<int K, typename T,int C,typename Palloc,typename Lalloc>
KDTree<K,T,C,Palloc,Lalloc>::BoxIter::BoxIter(const KDTree& tree, const BBox& box) : bbox(box) {
	seek(tree.root,0);
}

template<int K, typename T,int C,typename Palloc,typename Lalloc>
KDTree<K,T,C,Palloc,Lalloc>::BoxIter::BoxIter(const BoxIter& other)
		: parents(other.parents), leaf(other.leaf), iter(other.iter), bbox(other.bbox) { }


template<int K, typename T,int C,typename Palloc,typename Lalloc>
KDTree<K,T,C,Palloc,Lalloc>::BoxIter::BoxIter(BoxIter&& other)
		: parents(move(other.parents)), leaf(other.leaf), iter(move(other.iter)), bbox(other.bbox) { }

template<int K, typename T,int C,typename Palloc,typename Lalloc>
template<typename Point>
KDTree<K,T,C,Palloc,Lalloc>::BoxIter::BoxIter(const KDTree& tree, const Point& lowerLeft, const Point& upperRight) {
	for(int i = 0; i < K; ++i){
		bbox.lowerLeft[i] = lowerLeft[i];
		bbox.upperRight[i] = upperRight[i];
	}
	seek(*tree.root,0);
}

// starting from a leaf, goes to the next leaf in the bounded box, or nullptr if no next leaf exists
template<int K, typename T,int C,typename Palloc,typename Lalloc>
void KDTree<K,T,C,Palloc,Lalloc>::BoxIter::nextLeaf() {
	const KDNode<K,T,C> *child = leaf;
	while(parents.size() > 0){
		const ParentNode<K,T,C> *parent = parents.top();
		if(parent->left == child){
			int dim = parents.size() % K;
			dec<K>(dim);
			if(parent->value < bbox.upperRight[dim]){
				seek(*parent->right, dim == K - 1 ? 0 : dim + 1);
				return;
			}
		}
		parents.pop();
		child = parent;
	}
	leaf = nullptr; // reached end
}

// goes to the 'left' most leaf that intersects the bounding box below start
template<int K, typename T,int C,typename Palloc,typename Lalloc>
void KDTree<K,T,C,Palloc,Lalloc>::BoxIter::seek(const KDNode<K,T,C>& start, int dim) {
	const KDNode<K,T,C>* node = &start;
	do {
		while(!node->isLeaf){
			ParentNode<K,T,C>* parent = (ParentNode<K,T,C>*) node;
			if(parent->value > bbox.lowerLeft[dim]){
				parents.push(parent);
				node = parent->left;
				inc<K>(dim);
			} else {
				node = parents.top()->right;
			}
		}
		leaf = (LeafNode<K,T,C>*) node;
		iter = leaf->begin();
		while(iter != leaf->end()){
			if(bbox.contains(*((*iter).first))){
				return;
			}
			++iter;
		}
		nextLeaf();
		if(node == nullptr){
			return;
		}
		continue;
	} while(node != nullptr);
}

template<int K, typename T,int C,typename Palloc,typename Lalloc>
typename KDTree<K,T,C,Palloc,Lalloc>::BoxIter& KDTree<K,T,C,Palloc,Lalloc>::BoxIter::operator++() {
	if(leaf == nullptr){
		return *this;
	}
	do {
		if(iter == leaf->end()){
			nextLeaf();
			if(leaf == nullptr){
				return *this;
			}
		} else {
			++iter;
		}
	} while(!bbox.contains(*(*iter).first));
	return *this;
}


template<int K, typename T,int C,typename Palloc,typename Lalloc>
typename KDTree<K,T,C,Palloc,Lalloc>::BoxIter KDTree<K,T,C,Palloc,Lalloc>::BoxIter::operator++(int) {
	BoxIter old (*this);
	++(*this);
	return old;
}

template<int K, typename T,int C,typename Palloc,typename Lalloc>
entry<K,T> KDTree<K,T,C,Palloc,Lalloc>::BoxIter::operator*() const {
	return *iter;
}

template<int K, typename T,int C,typename Palloc,typename Lalloc>
typename KDTree<K,T,C,Palloc,Lalloc>::BoxIter& KDTree<K,T,C,Palloc,Lalloc>::BoxIter::operator=(const BoxIter& other) {
	parents = other.parents;
	leaf = other.leaf;
	iter = other.iter;
	bbox = other.bbox;
	return *this;
}

template<int K, typename T,int C,typename Palloc,typename Lalloc>
typename KDTree<K,T,C,Palloc,Lalloc>::BoxIter& KDTree<K,T,C,Palloc,Lalloc>::BoxIter::operator=(BoxIter&& other) {
	parents = move(other.parents);
	leaf = other.leaf;
	iter = move(other.iter);
	bbox = other.bbox;
	return *this;
}

template<int K, typename T,int C,typename Palloc,typename Lalloc>
bool KDTree<K,T,C,Palloc,Lalloc>::BoxIter::operator==(const BoxIter& iter) const {
	return (this->leaf == nullptr && iter.leaf == nullptr) || this->iter == iter.iter;
}

template<int K, typename T,int C,typename Palloc,typename Lalloc>
bool KDTree<K,T,C,Palloc,Lalloc>::BoxIter::operator!=(const BoxIter& iter) const {
	return !(*this == iter);
}

