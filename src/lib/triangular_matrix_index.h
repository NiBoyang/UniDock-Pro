#ifndef VINA_TRIANGULAR_MATRIX_INDEX_H
#define VINA_TRIANGULAR_MATRIX_INDEX_H

#include "common.h"

inline sz triangular_matrix_index(sz n, sz i, sz j) {
    assert(j < n);
    assert(i <= j);

    return i + j * (j + 1) / 2;
}

inline sz triangular_matrix_index_permissive(sz n, sz i, sz j) {
    return (i <= j) ? triangular_matrix_index(n, i, j) : triangular_matrix_index(n, j, i);
}

#endif
