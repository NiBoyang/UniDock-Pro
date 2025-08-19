#ifndef VINA_CURL_H
#define VINA_CURL_H

#include "common.h"

#if 1                  // prefer this to "hard curl"?
template <typename T>  // T = fl or vec
void curl(fl& e, T& deriv, fl v) {
    if (e > 0 && not_max(v)) {  // FIXME authentic_v can be gotten rid of everywhere now
        fl tmp = (v < epsilon_fl) ? 0 : (v / (v + e));
        e *= tmp;
        deriv *= sqr(tmp);
    }
}

inline void curl(fl& e, fl v) {
    if (e > 0 && not_max(v)) {
        fl tmp = (v < epsilon_fl) ? 0 : (v / (v + e));
        e *= tmp;
    }
}

#else

template <typename T>  // T = fl or vec
void curl(fl& e, T& deriv, fl v) {
    if (e > v) {
        e = v;
        deriv = 0;
    }
}

inline void curl(fl& e, fl v) {
    if (e > v) {
        e = v;
    }
}
#endif

#endif
