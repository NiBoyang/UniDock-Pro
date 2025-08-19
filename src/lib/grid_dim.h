#ifndef VINA_GRID_DIM_H
#define VINA_GRID_DIM_H

#include <boost/array.hpp>

#include "common.h"

struct grid_dim {
    fl begin;
    fl end;
    sz n_voxels;  // number of intervals == number of sample points - 1
    grid_dim() : begin(0), end(0), n_voxels(0) {}
    fl span() const { return end - begin; }
    bool enabled() const { return (n_voxels > 0); }

private:
    friend class boost::serialization::access;
    template <class Archive> void serialize(Archive& ar, const unsigned version) {
        ar & begin;
        ar & end;
        ar & n_voxels;
    }
};

inline bool eq(const grid_dim& a, const grid_dim& b) {
    return a.n_voxels == b.n_voxels && eq(a.begin, b.begin) && eq(a.end, b.end);
}

typedef boost::array<grid_dim, 3> grid_dims;

inline bool eq(const grid_dims& a, const grid_dims& b) {
    return eq(a[0], b[0]) && eq(a[1], b[1]) && eq(a[2], b[2]);
}

inline void print(const grid_dims& gd, std::ostream& out = std::cout) {
    VINA_FOR_IN(i, gd)
    std::cout << gd[i].n_voxels << " [" << gd[i].begin << " .. " << gd[i].end << "]\n";
}

inline vec grid_dims_begin(const grid_dims& gd) {
    vec tmp;
    VINA_FOR_IN(i, gd)
    tmp[i] = gd[i].begin;
    return tmp;
}

inline vec grid_dims_end(const grid_dims& gd) {
    vec tmp;
    VINA_FOR_IN(i, gd)
    tmp[i] = gd[i].end;
    return tmp;
}

#endif
