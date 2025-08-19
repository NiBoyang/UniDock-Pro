#ifndef VINA_INT_POW_H
#define VINA_INT_POW_H

#include "common.h"

template<int n>
inline fl int_pow(fl x) {
    return int_pow<n-1>(x) * x;
}

template<>
inline fl int_pow<0>(fl x) {
    return 1;
}

#endif