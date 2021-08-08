#ifndef MAXHEAP_HPP_INCLUDED
#define MAXHEAP_HPP_INCLUDED

#include <vector>
#include <iostream>

template<typename DataType>
class MaxHeap {
	size_t size;
	std::vector<DataType> *arr;
public:
	MaxHeap(size_t size, std::vector<DataType> *arr) {
		this->size = size;
		this->arr = arr;
		build();
	}
	void heapify(size_t i) {
		while (i < size / 2) {
			size_t left = 2 * i + 1;
			size_t right = 2 * i + 2;
			size_t root = (*arr)[i] < (*arr)[left] ? left : i;
			if (right < size)
				root = (*arr)[root] < (*arr)[right] ? right : root;
			if (root != i) {
				std::swap((*arr)[root], (*arr)[i]);
				i = root;
			}
			else {
				return;
			}
		}
	}
	void build() {
		if (size >= 2) {
			for (size_t i = (size_t)(size / 2) - 1;; i--) {
				heapify(i);
				if (i == 0) return;
			}
		}
	}
	void update(DataType d) {
		if ((*arr)[0] < d) return;
		(*arr)[0] = d;
		heapify(0);
	}
};
#endif