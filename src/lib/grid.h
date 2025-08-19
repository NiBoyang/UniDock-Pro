#ifndef VINA_GRID_H
#define VINA_GRID_H

#include "array3d.h"
#include "grid_dim.h"
#include "curl.h"
#include "file.h"

class grid {  // FIXME rm 'm_', consistent with my new style
public:
    vec m_init;          // DSM was private
    vec m_range;         // DSM was private
    vec m_factor_inv;    // DSM was private
    array3d<fl> m_data;  // FIXME? - make cache a friend, and convert this back to private?
    grid()
        : m_init(0, 0, 0),
          m_range(1, 1, 1),
          m_factor(1, 1, 1),
          m_dim_fl_minus_1(-1, -1, -1),
          m_factor_inv(1, 1, 1) {}  // not private
    grid(const grid_dims& gd) { init(gd); }
    void init(const grid_dims& gd);
    vec index_to_argument(sz x, sz y, sz z) const {
        return vec(m_init[0] + m_factor_inv[0] * x, m_init[1] + m_factor_inv[1] * y,
                   m_init[2] + m_factor_inv[2] * z);
    }
    bool initialized() const { return m_data.dim0() > 0 && m_data.dim1() > 0 && m_data.dim2() > 0; }
    fl evaluate(const vec& location, fl slope, fl c) const {
        return evaluate_aux(location, slope, c, NULL);
    }
    fl evaluate(const vec& location, fl slope, fl c, vec& deriv) const {
        return evaluate_aux(location, slope, c, &deriv);
    }  // sets deriv
    // add to public
    vec m_factor;
    vec m_dim_fl_minus_1;

private:
    fl evaluate_aux(const vec& location, fl slope, fl v,
                    vec* deriv) const;  // sets *deriv if not NULL
    friend class boost::serialization::access;
    template <class Archive> void serialize(Archive& ar, const unsigned version) {
        ar & m_init;
        ar & m_data;
        ar & m_range;
        ar & m_factor;
        ar & m_dim_fl_minus_1;
        ar & m_factor_inv;
    }
};

#endif
