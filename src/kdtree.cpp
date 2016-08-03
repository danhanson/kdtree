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

template<unsigned K, typename T,typename Ealloc,typename Palloc,typename Lalloc>
template<typename coord>
LeafData<K,T,Ealloc> KDTree<K,T,Ealloc,Palloc,Lalloc>::getLeaf(const coord& point){
	unsigned i (0);
	KDNode<K,T,Ealloc> **childPtr = root, *current = root;
	while(!current->isLeaf){
		double eVal = point[i++];
		ParentNode<K,T,Ealloc>* parent = (ParentNode<K,T,Ealloc>*) current;
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
	LeafNode<K,T,Ealloc>& leaf = (ParentNode<K,T,Ealloc>&) current;
	return std::make_pair(leaf, *childPtr);
}

template<unsigned K, typename T,typename Ealloc,typename Palloc,typename Lalloc>
template<typename coord>
void KDTree<K,T,Ealloc,Palloc,Lalloc>::insert(const coord& point, T&& elem){
	LeafData<K,T,Ealloc> leafData = getLeaf(point);
	LeafNode<K,T,Ealloc>& leaf = leafData.leaf;
	KDNode<K,T,Ealloc>*& childPtr = leafData.childPtr;

	entry<K,T>& entry = *leaf.value.insert(point[leafData.i]);
	array<double,K>& addedPoint = entry.first;

	entry.second(forward(elem)); // move element into tree
	for(unsigned i = 0; i < K; i++){ // copy coordinates into tree
		addedPoint[i] = point[i];
	}

	if(leaf.values.size > leafSize){ // leaf is to big, we will split it now
		LeafNode<K,T,Ealloc>* newLeft, newRight;
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
		ParentNode<K,T,Ealloc>* parent = Palloc::alocate(1);
		Palloc::construct(parent, value, *left, *right); // copy the right elements into the right child
		*childPtr = parent; // set the parent's pointer to the new child
	}
}



template<unsigned K, typename T,typename Ealloc,typename Palloc,typename Lalloc>
template<typename coord>
unsigned KDTree<K,T,Ealloc,Palloc,Lalloc>::erase(const coord& point){
	LeafData<K,T,Ealloc> leafData = getLeaf(point);
	LeafNode<K,T,Ealloc>& leaf = leafData.leaf;
	KDNode<K,T,Ealloc>& childPtr = leafData.childPtr;
	return *leaf.values.erase(point[leafData.i]);
}
