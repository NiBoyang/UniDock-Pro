#ifndef VINA_BRICK_H
#define VINA_BRICK_H

#include "common.h"

inline fl closest_between(fl begin, fl end, fl x) {
    assert(begin <= end);
    if (x <= begin)
        return begin;
    else if (x >= end)
        return end;
    return x;
}

inline vec brick_closest(const vec& begin, const vec& end, const vec& v) {
    vec tmp;
    VINA_FOR_IN(i, tmp)
    tmp[i] = closest_between(begin[i], end[i], v[i]);
    return tmp;
}

inline fl brick_distance_sqr(const vec& begin, const vec& end, const vec& v) {
    vec closest;
    closest = brick_closest(begin, end, v);
    return vec_distance_sqr(closest, v);
}

#endif
