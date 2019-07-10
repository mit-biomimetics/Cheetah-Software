/*!
 * @file: SparseMatrixMath.h
 *
 *
 */

#ifndef PROJECT_SPARSEMATRIXMATH_H
#define PROJECT_SPARSEMATRIXMATH_H

#include "types.h"

/*!
 * Compressed Sparse Column format for matrix
 */
template<typename T>
struct MatCSC {
  u32 nnz, m, n;
  T* values = nullptr;
  u32* colPtrs = nullptr;
  u32* rowIdx = nullptr;

  void alloc(u32 _n, u32 _nnz) {
    nnz = _nnz;
    n = _n;
    m = _n;
    values = new T[nnz];
    colPtrs = new u32[n + 1];
    rowIdx = new u32[nnz];
  }

  void freeAll() {
    delete[] values;
    delete[] colPtrs;
    delete[] rowIdx;
  }
};

/*!
 * Triplet format of sparse matrix
 */
template<typename T>
struct SparseTriple {
  T value;
  u32 r, c;
};

/*!
 * Check if a sorted list of triples contains entries for the same spot in the matrix
 */
template<typename T>
bool checkSortedTripleDuplicates(std::vector<SparseTriple<T>>& triples) {
  for(u64 i = 0; i < triples.size() - 1; i++) {
    auto& a = triples[i];
    auto& b = triples[i + 1];
    if(a.r == b.r && a.c == b.c) {
      return true;
    }
  }
  return false;
}

/*!
 * Sort triplets in column major.  If check duplicates is set, and there are duplicates, a
 * std::runtime_error exception is thrown.
 */
template<typename T>
void sortTriples(std::vector<SparseTriple<T>>& triples, bool checkDuplicates) {
  std::sort(triples.begin(), triples.end(), [](SparseTriple<T>& a, SparseTriple<T>& b){
    if(a.c == b.c) {
      return a.r < b.r;
    } else {
      return a.c < b.c;
    }
  });

  if(checkDuplicates && checkSortedTripleDuplicates(triples)) {
    throw std::runtime_error("duplicate found");
  }
}

/*!
 * Take a sorted list of triples and merges them to create a sorted list of unique triples
 */
template<typename T>
void sumSortedTriples(std::vector<SparseTriple<T>>& triples) {
  std::vector<SparseTriple<T>> temp = triples;
  u64 oldSize = triples.size();
  triples.clear();
  triples.reserve(oldSize);

  u32 lastRow = UINT32_MAX, lastCol = UINT32_MAX;

  for(auto& triple : temp) {
    if(triple.value == 0.0) continue;
    if(triple.r == lastRow && triple.c == lastCol) {
      triples.back().value += triple.value;
    } else {
      triples.push_back(triple);
    }
    lastRow = triple.r;
    lastCol = triple.c;
  }
}

/*!
 * Takes a list of unsorted, possibly non-unique triples, sorts and merges as needed to get
 * a list of sorted, unique triples.
 */
template<typename T>
void sortAndSumTriples(std::vector<SparseTriple<T>>& triples) {
  sortTriples(triples, false);
  sumSortedTriples(triples);
}


#endif //PROJECT_SPARSEMATRIXMATH_H
