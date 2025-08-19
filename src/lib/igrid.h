#ifndef VINA_IGRID_H
#define VINA_IGRID_H

#include "common.h"
#include "grid.h"

struct model;  // forward declaration

struct igrid {  // grids interface (that cache, etc. conform to)
    virtual fl eval(const model& m, fl v) const = 0;  // needs m.coords // clean up
    virtual fl eval_intra(model& m, fl v) const = 0;  // only flexres-grids
    virtual fl eval_deriv(model& m,
                          fl v) const
        = 0;  // needs m.coords, sets m.minus_forces // clean up
    virtual int get_atu() const = 0;
    virtual float get_slope() const = 0;
    virtual grid_dims get_gd() const = 0;
    virtual std::vector<grid> get_grids() const = 0;
};

#endif
