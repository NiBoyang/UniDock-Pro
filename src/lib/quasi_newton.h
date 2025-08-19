#ifndef VINA_QUASI_NEWTON_H
#define VINA_QUASI_NEWTON_H

#include "model.h"

struct quasi_newton {
    unsigned max_steps;
    fl average_required_improvement;
    quasi_newton() : max_steps(1000), average_required_improvement(0.0) {}
    // clean up
    void operator()(model& m, const precalculate_byatom& p, const igrid& ig, output_type& out,
                    change& g, const vec& v, int& evalcount) const;  // g must have correct size
};

#endif
