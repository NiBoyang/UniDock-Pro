#ifndef VINA_SZV_GRID_H
#define VINA_SZV_GRID_H

#include "model.h"
#include "grid_dim.h"
#include "array3d.h"

struct szv_grid {
    szv_grid() {}
    szv_grid(const model& m, const grid_dims& gd, fl cutoff_sqr);
    const szv& possibilities(const vec& coords) const;

private:
    array3d<szv> m_data;
    vec m_init;
    vec m_range;
    vec index_to_coord(sz i, sz j, sz k) const;
};

grid_dims szv_grid_dims(const grid_dims& gd);

#endif
