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

template<unsigned K, typename T,unsigned C,typename Palloc,typename Lalloc>
template<typename coord>
GetLeafResult<K,T,C> KDTree<K,T,C,Palloc,Lalloc>::getLeaf(const coord& point){
	unsigned i (0);
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

template<unsigned K, typename T,unsigned C,typename Palloc,typename Lalloc>
template<typename coord>
void KDTree<K,T,C,Palloc,Lalloc>::insert(const coord& point, T&& elem){
	GetLeafResult<K,T,C> leafData = getLeaf(point);
	LeafNode<K,T,C>& leaf = leafData.leaf;
	KDNode<K,T,C>*& childPtr = leafData.childPtr;

	entry<K,T>& entry = *leaf.value.insert(point[leafData.i]);
	array<double,K>& addedPoint = entry.first;

	entry.second(forward(elem)); // move element into tree
	for(unsigned i = 0; i < K; i++){ // copy coordinates into tree
		addedPoint[i] = point[i];
	}

	if(leaf.values.size > C){ // leaf is to big, we will split it now
		LeafNode<K,T,C>* newLeft, newRight;
		double value;
		newLeft = Lalloc::allocate(1);
		Lalloc::construct(newLeft);
		auto iter = leaf.values.begin();
		unsigned n = 0;
		while(n++ < leaf.values.size() / 2){ // copy the left elements into the left child
			newLeft->values.insert(*iter++);
		}

		value = iter->first; // new node's value

		newRight = Lalloc::allocate(1);
		Lalloc::construct(newRight, iter, leaf.values.end());
		ParentNode<K,T,C>* parent = Palloc::alocate(1);
		Palloc::construct(parent, value, *left, *right); // copy the right elements into the right child
		*childPtr = parent; // set the parent's pointer to the new child
	}
}



template<unsigned K, typename T,unsigned C,typename Palloc,typename Lalloc>
template<typename coord>
unsigned KDTree<K,T,C,Palloc,Lalloc>::erase(const coord& point){
	GetLeafResult<K,T,C> leafData = getLeaf(point);
	LeafNode<K,T,C>& leaf = leafData.leaf;
	KDNode<K,T,C>& childPtr = leafData.childPtr;
	return *leaf.values.erase(point[leafData.i]);
}

template<unsigned K, typename T, unsigned C>
LeafNode<K,T,C>::iterator::iterator() : i(0), leaf(nullptr) { }

template<unsigned K, typename T, unsigned C>
LeafNode<K,T,C>::iterator::iterator(const LeafNode& leaf) : i(0), leaf(leaf) { }

template<unsigned K, typename T, unsigned C>
typename LeafNode<K,T,C>::iterator& LeafNode<K,T,C>::iterator::operator++(){
	do {
		i++;
	} while(i < K && !leaf->isFilled(i));
	return *this;
}

template<unsigned K, typename T, unsigned C>
typename LeafNode<K,T,C>::iterator& LeafNode<K,T,C>::iterator::operator--(){
	do {
		i--;
	} while(i < K && !leaf->isFilled(i));
	return *this;
}

template<unsigned K, typename T, unsigned C>
typename LeafNode<K,T,C>::iterator& LeafNode<K,T,C>::iterator::operator=(const iterator& iter){
	return *this(iter);
}

template<unsigned K, typename T, unsigned C>
typename LeafNode<K,T,C>::iterator& LeafNode<K,T,C>::iterator::operator=(iterator&& iter){
	return *this(move(iter));
}

template<unsigned K, typename T, unsigned C>
entry<K,T> LeafNode<K,T,C>::iterator::operator*(){
	const coord<K>& point = leaf->coords[i];
	T& value = leaf->values[i];
	return make_pair(point, value);
}

template<unsigned K, typename T,unsigned C,typename Palloc,typename Lalloc>
KDTree<K,T,C,Palloc,Lalloc>::StdIter::StdIter(const KDTree& tree, bool isEnd)
		: tree(tree) {
	if(isEnd){
		this->isEnd = true;
		return;
	}
	auto leaf = seekLeft(*tree.root);
	this->map = leaf->values;
	this->iter = leaf->values.begin();
	this->isEnd = this->iter == leaf->values.end();
}

template<unsigned K, typename T,unsigned C,typename Palloc,typename Lalloc>
LeafNode<K,T,C>& KDTree<K,T,C,Palloc,Lalloc>::StdIter::seekLeft(const KDNode<K,T,C>& parent){
	KDNode<K,T,C>* node = &parent;
	while(!node->isLeaf){
		auto parent = (ParentNode<K,T,C>*) node;
		this->parents.push(parent);
		node = parent->left;
	}
	return *((LeafNode<K,T,C>*) node);
}

template<unsigned K, typename T,unsigned C,typename Palloc,typename Lalloc>
LeafNode<K,T,C>& KDTree<K,T,C,Palloc,Lalloc>::StdIter::seekRight(const KDNode<K,T,C>& parent){
	KDNode<K,T,C>* node = &parent;
	while(!node->isLeaf){
		auto parent = (ParentNode<K,T,C>*) node;
		this->parents.push(parent);
		node = parent->right;
	}
	return *((LeafNode<K,T,C>*) node);
}

template<unsigned K, typename T,unsigned C,typename Palloc,typename Lalloc>
typename KDTree<K,T,C,Palloc,Lalloc>::StdIter& KDTree<K,T,C,Palloc,Lalloc>::StdIter::operator++(){
	if(isEnd){ // iterator is at end and/or uninitialized
		return *this;
	}
	this->iter++;
	while(this->iter == this->leaf.end()){
		KDNode<K,T,C>* child, *parent = this->leaf;
		do {
			if(this->stack.size() == 0){ // reached end
				isEnd = true;
				return *this;
			}
			child = parent;
			parent = this->stack.pop();
		} while(parent->right == child);
		this->leaf = seekLeft(*(parent->right));
		this->iter = this->leaf->begin();
	}
	return *this;
}

template<unsigned K, typename T,unsigned C,typename Palloc,typename Lalloc>
typename KDTree<K,T,C,Palloc,Lalloc>::StdIter& KDTree<K,T,C,Palloc,Lalloc>::StdIter::operator--(){
	if(isEnd){
		if(this->tree == nullptr){ // iterator is uninitialized
			return *this;
		}
		if(this->leaf == nullptr){ // we need to calculate the last leaf
			this->leaf = seekRight(*(tree->root));
			this->iter = this->leaf->end(); // we decrement iter in the rest of this function
		}
	}
	while(this->iter == this->leaf->begin()){
		KDNode<K,T,C>* child, parent = this->leaf;
		do {
			if(this->stack.size() == 0){ // decremented passed begin
				this->leaf = seekRight(*(tree->root)); // move iterator back to front
				this->iter = this->leaf->begin();
				return *this;
			}
			child = parent;
			parent = this->stack.pop();
		} while(parent->left == child);
		this->leaf = seekRight(*(parent->left));
		this->iter = leaf->values->end();
	}
	this->iter--;
	return *this;
}

template<unsigned K, typename T,unsigned C,typename Palloc,typename Lalloc>
typename KDTree<K,T,C,Palloc,Lalloc>::StdIter& KDTree<K,T,C,Palloc,Lalloc>::StdIter::operator=(const StdIter& iter){
	this(iter);
	return *this;
}

template<unsigned K, typename T,unsigned C,typename Palloc,typename Lalloc>
typename KDTree<K,T,C,Palloc,Lalloc>::StdIter& KDTree<K,T,C,Palloc,Lalloc>::StdIter::operator=(StdIter&& iter){
	return *this(move(iter));
}

template<unsigned K, typename T,unsigned C,typename Palloc,typename Lalloc>
bool KDTree<K,T,C,Palloc,Lalloc>::StdIter::operator==(const StdIter& iter){
	if(this->isEnd && iter.isEnd && this->tree == iter.tree){
		return true;
	}
	return this->iter == iter.iter;
}

template<unsigned K, typename T,unsigned C,typename Palloc,typename Lalloc>
bool KDTree<K,T,C,Palloc,Lalloc>::StdIter::operator!=(const StdIter& iter){
	return !(*this == iter);
}

template<unsigned K, typename T,unsigned C,typename Palloc,typename Lalloc>
entry<K,T> KDTree<K,T,C,Palloc,Lalloc>::StdIter::operator*(){
	return *(this->iter);
}


