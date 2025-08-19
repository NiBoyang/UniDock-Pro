#include "quasi_newton.h"
#include "bfgs.h"

struct quasi_newton_aux {
    model* m;
    const precalculate_byatom* p;
    const igrid* ig;
    const vec v;

    quasi_newton_aux(model* m_, const precalculate_byatom* p_, const igrid* ig_, const vec& v_)
        : m(m_), p(p_), ig(ig_), v(v_) {}

    fl operator()(const conf& c, change& g) {
        // Before evaluating conf, we have to update model
        m->set(c);
        const fl tmp = m->eval_deriv(*p, *ig, v, g);
        return tmp;
    }
};

void quasi_newton::operator()(model& m, const precalculate_byatom& p, const igrid& ig,
                              output_type& out, change& g, const vec& v,
                              int& evalcount) const {  // g must have correct size
    quasi_newton_aux aux(&m, &p, &ig, v);

    fl res = bfgs(aux, out.c, g, max_steps, average_required_improvement, 10, evalcount);

    // Update model a last time after optimization
    m.set(out.c);
    out.e = res;
}
