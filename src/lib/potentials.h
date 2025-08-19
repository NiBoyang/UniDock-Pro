#ifndef VINA_POTENTIALS_H
#define VINA_POTENTIALS_H

#include "atom.h"
#include "int_pow.h"

// Vina common functions
inline fl slope_step(fl x_bad, fl x_good, fl x) {
    if (x_bad < x_good) {
        if (x <= x_bad) return 0;
        if (x >= x_good) return 1;
    } else {
        if (x >= x_bad) return 0;
        if (x <= x_good) return 1;
    }
    return (x - x_bad) / (x_good - x_bad);
}

inline bool is_glue_type(sz xs_t) {
    if ((xs_t == XS_TYPE_G0) || (xs_t == XS_TYPE_G1) || (xs_t == XS_TYPE_G2)
        || (xs_t == XS_TYPE_G3))
        return true;
    return false;
}

inline fl optimal_distance(sz xs_t1, sz xs_t2) {
    if (is_glue_type(xs_t1) || is_glue_type(xs_t2)) return 0.0;  // G0, G1, G2 or G3
    return xs_radius(xs_t1) + xs_radius(xs_t2);
}

inline fl smooth_div(fl x, fl y) {
    if (std::abs(x) < epsilon_fl) return 0;
    if (std::abs(y) < epsilon_fl)
        return ((x * y > 0) ? max_fl : -max_fl);  // FIXME I hope -max_fl does not become NaN
    return x / y;
}

// Vinardo common functions
inline fl optimal_distance_vinardo(sz xs_t1, sz xs_t2) {
    if (is_glue_type(xs_t1) || is_glue_type(xs_t2)) return 0.0;  // G0, G1, G2 or G3
    return xs_vinardo_radius(xs_t1) + xs_vinardo_radius(xs_t2);
}

// Macrocycle - Vina
inline bool is_glued(sz xs_t1, sz xs_t2) {
    return (xs_t1 == XS_TYPE_G0 && xs_t2 == XS_TYPE_C_H_CG0)
           || (xs_t1 == XS_TYPE_G0 && xs_t2 == XS_TYPE_C_P_CG0)
           || (xs_t2 == XS_TYPE_G0 && xs_t1 == XS_TYPE_C_H_CG0)
           || (xs_t2 == XS_TYPE_G0 && xs_t1 == XS_TYPE_C_P_CG0) ||

           (xs_t1 == XS_TYPE_G1 && xs_t2 == XS_TYPE_C_H_CG1)
           || (xs_t1 == XS_TYPE_G1 && xs_t2 == XS_TYPE_C_P_CG1)
           || (xs_t2 == XS_TYPE_G1 && xs_t1 == XS_TYPE_C_H_CG1)
           || (xs_t2 == XS_TYPE_G1 && xs_t1 == XS_TYPE_C_P_CG1) ||

           (xs_t1 == XS_TYPE_G2 && xs_t2 == XS_TYPE_C_H_CG2)
           || (xs_t1 == XS_TYPE_G2 && xs_t2 == XS_TYPE_C_P_CG2)
           || (xs_t2 == XS_TYPE_G2 && xs_t1 == XS_TYPE_C_H_CG2)
           || (xs_t2 == XS_TYPE_G2 && xs_t1 == XS_TYPE_C_P_CG2) ||

           (xs_t1 == XS_TYPE_G3 && xs_t2 == XS_TYPE_C_H_CG3)
           || (xs_t1 == XS_TYPE_G3 && xs_t2 == XS_TYPE_C_P_CG3)
           || (xs_t2 == XS_TYPE_G3 && xs_t1 == XS_TYPE_C_H_CG3)
           || (xs_t2 == XS_TYPE_G3 && xs_t1 == XS_TYPE_C_P_CG3);
}

class Potential {
public:
    virtual ~Potential() {}
    virtual fl eval(const atom& a, const atom& b, fl r) { return 0; };
    virtual fl eval(sz t1, sz t2, fl r) { return 0; };
    virtual fl get_cutoff() { return 0; }
};

// Vina
class vina_gaussian : public Potential {
public:
    vina_gaussian(fl offset_, fl width_, fl cutoff_)
        : offset(offset_), width(width_), cutoff(cutoff_) {}
    //~vina_gaussian() { }
    fl eval(const atom& a, const atom& b, fl r) {
        if (r >= cutoff) return 0.0;
        if ((a.xs >= XS_TYPE_SIZE) || (b.xs >= XS_TYPE_SIZE)) return 0.0;
        return gauss(r - (optimal_distance(a.xs, b.xs) + offset));  // hard-coded to XS atom type
    };
    fl eval(sz t1, sz t2, fl r) {
        if (r >= cutoff) return 0.0;
        return gauss(r - (optimal_distance(t1, t2) + offset));  // hard-coded to XS atom type
    };
    fl get_cutoff() { return cutoff; }

    fl offset;  // added to optimal distance
    fl width;
    fl cutoff;

    fl gauss(fl x) { return std::exp(-sqr(x / width)); };
};

class vina_repulsion : public Potential {
public:
    vina_repulsion(fl offset_, fl cutoff_) : offset(offset_), cutoff(cutoff_) {}
    //~vina_repulsion() { }
    fl eval(const atom& a, const atom& b, fl r) {
        if (r >= cutoff) return 0.0;
        if ((a.xs >= XS_TYPE_SIZE) || (b.xs >= XS_TYPE_SIZE)) return 0.0;
        fl d = r - (optimal_distance(a.xs, b.xs) + offset);  // hard-coded to XS atom type
        if (d > 0.0) return 0.0;
        return d * d;
    };
    fl eval(sz t1, sz t2, fl r) {
        if (r >= cutoff) return 0.0;
        fl d = r - (optimal_distance(t1, t2) + offset);  // hard-coded to XS atom type
        if (d > 0.0) return 0.0;
        return d * d;
    };
    fl get_cutoff() { return cutoff; }

    fl offset;  // added to vdw
    fl cutoff;
};

class vina_hydrophobic : public Potential {
public:
    vina_hydrophobic(fl good_, fl bad_, fl cutoff_) : good(good_), bad(bad_), cutoff(cutoff_) {}
    //~vina_hydrophobic() { }
    fl eval(const atom& a, const atom& b, fl r) {
        if (r >= cutoff) return 0.0;
        if ((a.xs >= XS_TYPE_SIZE) || (b.xs >= XS_TYPE_SIZE)) return 0.0;
        if (xs_is_hydrophobic(a.xs) && xs_is_hydrophobic(b.xs))
            return slope_step(bad, good, r - optimal_distance(a.xs, b.xs));
        else
            return 0.0;
    };
    fl eval(sz t1, sz t2, fl r) {
        if (r >= cutoff) return 0.0;
        if (xs_is_hydrophobic(t1) && xs_is_hydrophobic(t2))
            return slope_step(bad, good, r - optimal_distance(t1, t2));
        else
            return 0.0;
    };
    fl get_cutoff() { return cutoff; }

    fl good;
    fl bad;
    fl cutoff;
};

class vina_non_dir_h_bond : public Potential {
public:
    vina_non_dir_h_bond(fl good_, fl bad_, fl cutoff_) : good(good_), bad(bad_), cutoff(cutoff_) {}
    //~vina_non_dir_h_bond() { }
    fl eval(const atom& a, const atom& b, fl r) {
        if (r >= cutoff) return 0.0;
        if ((a.xs >= XS_TYPE_SIZE) || (b.xs >= XS_TYPE_SIZE)) return 0.0;
        if (xs_h_bond_possible(a.xs, b.xs)) {
            if ((a.xs >= 32 && a.xs <= 35) || (b.xs >= 32 && b.xs <= 35)) {
                return 10.0 * slope_step(bad, good, r - optimal_distance(a.xs, b.xs));
            } else
                return slope_step(bad, good, r - optimal_distance(a.xs, b.xs));
        }

        return 0.0;
    };
    fl eval(sz t1, sz t2, fl r) {
        if (r >= cutoff) return 0.0;
        if (xs_h_bond_possible(t1, t2)) {
            if (t1 == XS_TYPE_O_XA || t1 == XS_TYPE_N_XA || t1 == XS_TYPE_O_XD || t1 == XS_TYPE_N_XD
                || t2 == XS_TYPE_O_XA || t2 == XS_TYPE_N_XA || t2 == XS_TYPE_O_XD
                || t2 == XS_TYPE_N_XD) {
                return 10.0 * slope_step(bad, good, r - optimal_distance(t1, t2));
            } else {
                return slope_step(bad, good, r - optimal_distance(t1, t2));
            }
        }
        return 0.0;
    };
    fl get_cutoff() { return cutoff; }

    fl good;
    fl bad;
    fl cutoff;
};

// Vinardo
class vinardo_gaussian : public Potential {
public:
    vinardo_gaussian(fl offset_, fl width_, fl cutoff_)
        : offset(offset_), width(width_), cutoff(cutoff_) {}
    //~vina_gaussian() { }
    fl eval(const atom& a, const atom& b, fl r) {
        if (r >= cutoff) return 0.0;
        if ((a.xs >= XS_TYPE_SIZE) || (b.xs >= XS_TYPE_SIZE)) return 0.0;
        return gauss(
            r - (optimal_distance_vinardo(a.xs, b.xs) + offset));  // hard-coded to XS atom type
    };
    fl eval(sz t1, sz t2, fl r) {
        if (r >= cutoff) return 0.0;
        return gauss(r
                     - (optimal_distance_vinardo(t1, t2) + offset));  // hard-coded to XS atom type
    };
    fl get_cutoff() { return cutoff; }

    fl offset;  // added to optimal distance
    fl width;
    fl cutoff;

    fl gauss(fl x) { return std::exp(-sqr(x / width)); };
};

class vinardo_repulsion : public Potential {
public:
    vinardo_repulsion(fl offset_, fl cutoff_) : offset(offset_), cutoff(cutoff_) {}
    //~vina_repulsion() { }
    fl eval(const atom& a, const atom& b, fl r) {
        if (r >= cutoff) return 0.0;
        if ((a.xs >= XS_TYPE_SIZE) || (b.xs >= XS_TYPE_SIZE)) return 0.0;
        fl d = r - (optimal_distance_vinardo(a.xs, b.xs) + offset);  // hard-coded to XS atom type
        if (d > 0.0) return 0.0;
        return d * d;
    };
    fl eval(sz t1, sz t2, fl r) {
        if (r >= cutoff) return 0.0;
        fl d = r - (optimal_distance_vinardo(t1, t2) + offset);  // hard-coded to XS atom type
        if (d > 0.0) return 0.0;
        return d * d;
    };
    fl get_cutoff() { return cutoff; }

    fl offset;  // added to vdw
    fl cutoff;
};

class vinardo_hydrophobic : public Potential {
public:
    vinardo_hydrophobic(fl good_, fl bad_, fl cutoff_) : good(good_), bad(bad_), cutoff(cutoff_) {}
    //~vina_hydrophobic() { }
    fl eval(const atom& a, const atom& b, fl r) {
        if (r >= cutoff) return 0.0;
        if ((a.xs >= XS_TYPE_SIZE) || (b.xs >= XS_TYPE_SIZE)) return 0.0;
        if (xs_is_hydrophobic(a.xs) && xs_is_hydrophobic(b.xs))
            return slope_step(bad, good, r - optimal_distance_vinardo(a.xs, b.xs));
        else
            return 0.0;
    };
    fl eval(sz t1, sz t2, fl r) {
        if (r >= cutoff) return 0.0;
        if (xs_is_hydrophobic(t1) && xs_is_hydrophobic(t2))
            return slope_step(bad, good, r - optimal_distance_vinardo(t1, t2));
        else
            return 0.0;
    };
    fl get_cutoff() { return cutoff; }

    fl good;
    fl bad;
    fl cutoff;
};

class vinardo_non_dir_h_bond : public Potential {
public:
    vinardo_non_dir_h_bond(fl good_, fl bad_, fl cutoff_)
        : good(good_), bad(bad_), cutoff(cutoff_) {}
    //~vina_non_dir_h_bond() { }
    fl eval(const atom& a, const atom& b, fl r) {
        if (r >= cutoff) return 0.0;
        if ((a.xs >= XS_TYPE_SIZE) || (b.xs >= XS_TYPE_SIZE)) return 0.0;
        if (xs_h_bond_possible(a.xs, b.xs))
            return slope_step(bad, good, r - optimal_distance_vinardo(a.xs, b.xs));
        return 0.0;
    };
    fl eval(sz t1, sz t2, fl r) {
        if (r >= cutoff) return 0.0;
        if (xs_h_bond_possible(t1, t2))
            return slope_step(bad, good, r - optimal_distance_vinardo(t1, t2));
        return 0.0;
    };
    fl get_cutoff() { return cutoff; }

    fl good;
    fl bad;
    fl cutoff;
};

// Macrocycle - Vina
class linearattraction : public Potential {
public:
    linearattraction(fl cutoff_) : cutoff(cutoff_) {}
    fl eval(const atom& a, const atom& b, fl r) {
        if (r >= cutoff) return 0.0;
        if (is_glued(a.xs, b.xs))
            return r;
        else
            return 0.0;
    };
    fl eval(sz t1, sz t2, fl r) {
        if (r >= cutoff) return 0.0;
        if (is_glued(t1, t2))
            return r;
        else
            return 0.0;
    };
    fl get_cutoff() { return cutoff; }

    fl cutoff;
};

#endif
