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
// the parents
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
		int x = 0;
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
	KDNode<K,T,C> **childPtr = root, *current = root;
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
	LeafNode<K,T,C>& leaf = (ParentNode<K,T,C>&) current;
	return GetLeafResult<K,T,C>(leaf, *childPtr, i);
}


template<int K, typename T,int C,typename Palloc,typename Lalloc>
template<typename Point>
void KDTree<K,T,C,Palloc,Lalloc>::insert(const Point& point, T&& elem){
	GetLeafResult<K,T,C> leafData = getLeaf(point);
	LeafNode<K,T,C>& leaf = leafData.leaf;
	KDNode<K,T,C>*& childPtr = leafData.childPtr;

	int i = 0;
	while(leaf.isFilled[i++]){
		if(i == C){
			ParentNode<K,T,C>& parent = leaf.split(leafData.i);
			if(point[leafData.i] > parent.value){
				&leaf = parent.right;
				i = (C + 1) / 2;
			} else {
				&leaf = parent.left;
				i = C / 2;
			}
		}
	}
	leaf.isFilled[i] = true;
	for(int n = 0; n < K; n++){ // copy coordinates into tree
		leaf.coords[i][n] = point[n];
	}
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
template<typename Palloc, typename Lalloc>
ParentNode<K,T,C> LeafNode<K,T,C>::split(int dim) {
	Palloc palloc;
	Lalloc lalloc;

	// quick search to find median element
	int i = 0;
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
			swap(values[above], coords[below]);
		}
		if(below < values + K / 2 - 1){ // we will choose the left median when there is even values
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

	ParentNode<K,T,C>* parent = palloc.alocate(1);
	palloc.construct(parent, min[0][dim], *newLeft, *newRight); // copy the right elements into the right child

	// partitioning of the values using the median
	coord<K>* left = newLeft->values;
	coord<K>* right = newRight->values;
	for(coord<K>* iter = values.begin(); iter != values.end(); ++iter){
		if(*iter > parent->value){
			*right = *iter;
			++right;
		} else {
			*left = *iter;
			++left;
		}
	}

	fill(newLeft->isFilled, newLeft->isFilled.begin() + C/2, true);
	fill(newLeft->isFilled + C/2, newLeft->end(), false);

	fill(newRight->isFilled, newRight->isFilled + (C+1)/2, true);
	fill(newRight->isFilled + (C+1)/2, newRight->end(), false);
	lalloc.deallocate(this,1);
	return *parent;
}

template<int K, typename T, int C>
typename LeafNode<K,T,C>::iterator LeafNode<K,T,C>::begin() const {
	return iterator(*this,false);
}

template<int K, typename T, int C>
typename LeafNode<K,T,C>::iterator LeafNode<K,T,C>::end() const {
	return iterator(*this,true);
}

// simple insertion sort, performs well for small arrays
template<int K, typename T, int C>
void LeafNode<K,T,C>::sort(int dim){
	for(int i = 1; i < K; i++) {
		for(int j = i; j > 0; j--){
			if(coords[j][dim] < coords[j-1][dim]){
				swap(coords[j], coords[j-1]);
				swap(values[j], values[j-1]);
				swap(isFilled[j], isFilled[j-1]);
			} else {
				break;
			}
		}
	}
}

template<int K, typename T, int C>
LeafNode<K,T,C>::iterator::iterator() : leaf(nullptr), i() { }

template<int K, typename T, int C>
LeafNode<K,T,C>::iterator::iterator(const LeafNode& leaf, bool isEnd) : leaf(&leaf), i(0) {
	if(isEnd){
		i = K;
	} else {
		do {
			i++;
		} while(i < K && !leaf.isFilled[i]);
	}
}

template<int K, typename T, int C>
LeafNode<K,T,C>::iterator::iterator(const iterator& it) : leaf(it.leaf), i(it.i) { }

template<int K, typename T, int C>
typename LeafNode<K,T,C>::iterator& LeafNode<K,T,C>::iterator::operator++(){
	do {
		i++;
	} while(i < K && !leaf->isFilled[i]);
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
/*
template<int K, typename T, int C>
typename LeafNode<K,T,C>::iterator& LeafNode<K,T,C>::iterator::operator=(iterator&& iter){
	return (*this)(iter);
}
*/
template<int K, typename T, int C>
entry<K,T> LeafNode<K,T,C>::iterator::operator*() const {
	const coord<K>& point = leaf->coords[i];
	T& value = const_cast<T&>(leaf->values[i]); // XXX WHY IS leaf->values[i] CONST?
	pair<const coord<K>&, T&> pair (point, value);
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
	this->iter++;
	while(this->iter == this->leaf->end()){
		const KDNode<K,T,C> *child = this->leaf;
		while(parents.size() > 0){
			const ParentNode<K,T,C> *parent = parents.top();
			if(parent->left == child){
				this->leaf = &seekLeft(*(parent->right));
				this->iter = this->leaf->begin();
				goto outer;
			}
			child = parent;
		}
outer:;
	}
	this->leaf = nullptr;
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
	return *(this->iter);
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
	while(false) {
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
		if(iter == leaf->end()){
			nextLeaf();
			if(node == nullptr){
				return;
			}
			continue;
		}
	}
}

template<int K, typename T,int C,typename Palloc,typename Lalloc>
typename KDTree<K,T,C,Palloc,Lalloc>::BoxIter& KDTree<K,T,C,Palloc,Lalloc>::BoxIter::operator++() {
	if(leaf == nullptr){
		return *this;
	}
	do {
		++iter;
		if(iter == leaf->end()){
			nextLeaf();
		}
	} while(bbox.contains((*iter).first));
	return *this;
}


template<int K, typename T,int C,typename Palloc,typename Lalloc>
typename KDTree<K,T,C,Palloc,Lalloc>::BoxIter KDTree<K,T,C,Palloc,Lalloc>::BoxIter::operator++(int) {
	BoxIter old (*this);
	++(*this);
	return old;
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
	return this->iter == iter.iter;
}

template<int K, typename T,int C,typename Palloc,typename Lalloc>
bool KDTree<K,T,C,Palloc,Lalloc>::BoxIter::operator!=(const BoxIter& iter) const {
	return !(*this == iter);
}

template<int K, typename T,int C,typename Palloc,typename Lalloc>
entry<K,T> KDTree<K,T,C,Palloc,Lalloc>::BoxIter::operator*() const {
	return *(this->iter);
}
