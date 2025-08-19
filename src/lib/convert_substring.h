#ifndef VINA_CONVERT_SUBSTRING_H
#define VINA_CONVERT_SUBSTRING_H

#include <cctype>  // for isspace
#include <boost/lexical_cast.hpp>
#include "common.h"

struct bad_conversion {};

template <typename T>
T convert_substring(const std::string& str, sz i,
                    sz j) {  // indexes are 1-based, the substring should be non-null
    if (i < 1 || i > j + 1 || j > str.size()) throw bad_conversion();

    // omit leading whitespace
    while (i <= j && std::isspace(str[i - 1])) ++i;

    T tmp;
    try {
        tmp = boost::lexical_cast<T>(str.substr(i - 1, j - i + 1));
    } catch (...) {
        throw bad_conversion();
    }
    return tmp;
}

inline bool substring_is_blank(const std::string& str, sz i,
                               sz j) {  // indexes are 1-based, the substring should be non-null
    if (i < 1 || i > j + 1 || j > str.size()) throw bad_conversion();
    VINA_RANGE(k, i - 1, j)
    if (!std::isspace(str[k])) return false;
    return true;
}

// when this was written, lexical cast to unsigned didn't work very well with "-123", etc.
template <> inline unsigned convert_substring<unsigned>(
    const std::string& str, sz i, sz j) {  // indexes are 1-based, the substring should be non-null
    int tmp = convert_substring<int>(str, i, j);
    if (tmp < 0) throw bad_conversion();
    return static_cast<unsigned>(tmp);
}

#endif
