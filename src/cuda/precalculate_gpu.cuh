#include "math.h"

// Define GPU precalculate structures

/* atom related start */
#include "atom_constants_gpu.cuh"

/* atom related end */

/* potential related start */

__device__ __forceinline__ fl sqr(fl x) { return x * x; }

// Vina common functions
__device__ __forceinline__ fl slope_step_gpu(fl x_bad, fl x_good, fl x) {
    if (x_bad < x_good) {
        if (x <= x_bad) return 0;
        if (x >= x_good) return 1;
    } else {
        if (x >= x_bad) return 0;
        if (x <= x_good) return 1;
    }
    return (x - x_bad) / (x_good - x_bad);
}

__device__ __forceinline__ bool is_glue_type_gpu(sz xs_t) {
    if ((xs_t == XS_TYPE_G0) || (xs_t == XS_TYPE_G1) || (xs_t == XS_TYPE_G2)
        || (xs_t == XS_TYPE_G3))
        return true;
    return false;
}

__device__ __forceinline__ fl optimal_distance_gpu(sz xs_t1, sz xs_t2) {
    if (is_glue_type_gpu(xs_t1) || is_glue_type_gpu(xs_t2)) return 0.0;  // G0, G1, G2 or G3
    return xs_radius_gpu(xs_t1) + xs_radius_gpu(xs_t2);
}

__device__ __forceinline__ fl smooth_div_gpu(fl x, fl y) {
    if (std::abs(x) < epsilon_fl) return 0;
    if (std::abs(y) < epsilon_fl)
        return ((x * y > 0) ? max_fl : -max_fl);  // FIXME I hope -max_fl does not become NaN
    return x / y;
}

// Vinardo common functions
__device__ __forceinline__ fl optimal_distance_vinardo_gpu(sz xs_t1, sz xs_t2) {
    if (is_glue_type_gpu(xs_t1) || is_glue_type_gpu(xs_t2)) return 0.0;  // G0, G1, G2 or G3
    return xs_vinardo_radius_gpu(xs_t1) + xs_vinardo_radius_gpu(xs_t2);
}

// smoothing helper
__device__ __forceinline__ fl smoothen_gpu(fl r, fl rij, fl smoothing) {
    fl out;
    smoothing *= 0.5;

    if (r > rij + smoothing)
        out = r - smoothing;
    else if (r < rij - smoothing)
        out = r + smoothing;
    else
        out = rij;

    return out;
}


// Macrocycle - Vina
__device__ __forceinline__ bool is_glued_gpu(sz xs_t1, sz xs_t2) {
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

// Vina

__device__ __forceinline__ fl gauss_gpu(fl x, fl width) { return exp(-sqr(x / width)); };

__device__ __forceinline__ fl vina_gaussian_cuda_eval(sz t1, sz t2, fl r, fl cutoff, fl offset,
                                                      fl width) {
    if (r >= cutoff) return 0.0;
    if ((t1 >= XS_TYPE_SIZE) || (t2 >= XS_TYPE_SIZE)) return 0.0;
    return gauss_gpu(r - (optimal_distance_gpu(t1, t2) + offset),
                     width);  // hard-coded to XS atom type
};

__device__ __forceinline__ fl vina_repulsion_cuda_eval(sz t1, sz t2, fl r, fl cutoff, fl offset) {
    if (r >= cutoff) return 0.0;
    if ((t1 >= XS_TYPE_SIZE) || (t2 >= XS_TYPE_SIZE)) return 0.0;
    fl d = r - (optimal_distance_gpu(t1, t2) + offset);  // hard-coded to XS atom type
    if (d > 0.0) return 0.0;
    return d * d;
};

__device__ __forceinline__ fl vina_hydrophobic_cuda_eval(sz t1, sz t2, fl r, fl good, fl bad,
                                                         fl cutoff) {
    if (r >= cutoff) return 0.0;
    if ((t1 >= XS_TYPE_SIZE) || (t2 >= XS_TYPE_SIZE)) return 0.0;
    if (xs_is_hydrophobic_gpu(t1) && xs_is_hydrophobic_gpu(t2))
        return slope_step_gpu(bad, good, r - optimal_distance_gpu(t1, t2));
    else
        return 0.0;
};

__device__ __forceinline__ fl vina_non_dir_h_bond_cuda_eval(sz t1, sz t2, fl r, fl good, fl bad,
                                                            fl cutoff) {
    if (r >= cutoff) return 0.0;
    if ((t1 >= XS_TYPE_SIZE) || (t2 >= XS_TYPE_SIZE)) return 0.0;
    if (xs_h_bond_possible_gpu(t1, t2)){        
        if  ((t1 >= 32 && t1 <= 35) || (t2 >= 32 && t2 <= 35))
            return 10.0*slope_step_gpu(bad, good, r - optimal_distance_gpu(t1, t2));
        else 
            return slope_step_gpu(bad, good, r - optimal_distance_gpu(t1, t2));}
    return 0.0;
};

// Vinardo
__device__ __forceinline__ fl vinardo_gaussian_eval(sz t1, sz t2, fl r, fl offset, fl width,
                                                    fl cutoff) {
    if (r >= cutoff) return 0.0;
    if ((t1 >= XS_TYPE_SIZE) || (t2 >= XS_TYPE_SIZE)) return 0.0;
    return gauss_gpu(r - (optimal_distance_vinardo_gpu(t1, t2) + offset),
                     width);  // hard-coded to XS atom type
};

__device__ __forceinline__ fl vinardo_repulsion_eval(sz t1, sz t2, fl r, fl cutoff, fl offset) {
    if (r >= cutoff) return 0.0;
    if ((t1 >= XS_TYPE_SIZE) || (t2 >= XS_TYPE_SIZE)) return 0.0;
    fl d = r - (optimal_distance_vinardo_gpu(t1, t2) + offset);  // hard-coded to XS atom type
    if (d > 0.0) return 0.0;
    return d * d;
};

__device__ __forceinline__ fl vinardo_hydrophobic_eval(sz t1, sz t2, fl r, fl good, fl bad,
                                                       fl cutoff) {
    if (r >= cutoff) return 0.0;
    if ((t1 >= XS_TYPE_SIZE) || (t2 >= XS_TYPE_SIZE)) return 0.0;
    if (xs_is_hydrophobic_gpu(t1) && xs_is_hydrophobic_gpu(t2))
        return slope_step_gpu(bad, good, r - optimal_distance_vinardo_gpu(t1, t2));
    else
        return 0.0;
};

__device__ __forceinline__ fl vinardo_non_dir_h_bond_eval(sz t1, sz t2, fl r, fl good, fl bad,
                                                          fl cutoff) {
    if (r >= cutoff) return 0.0;
    if ((t1 >= XS_TYPE_SIZE) || (t2 >= XS_TYPE_SIZE)) return 0.0;
    if (xs_h_bond_possible_gpu(t1, t2))
        return slope_step_gpu(bad, good, r - optimal_distance_vinardo_gpu(t1, t2));
    return 0.0;
};


// Macrocycle - Vina
// Cutoff is large (20.0), may be biased if max_cutoff equals 8.0
__device__ __forceinline__ fl linearattraction_eval(sz t1, sz t2, fl r, fl cutoff) {
    if (r >= cutoff) return 0.0;
    if (is_glued_gpu(t1, t2))
        return r;
    else
        return 0.0;
};

/* potential related end */

/* scoring function related start */

typedef struct {
    int m_num_potentials;
    fl m_weights[6];
    int m_sf_choice;  // 0:vina, 1:vinardo
    // constants used in potential terms
    fl vina_gaussian_offset_1, vina_gaussian_width_1, vina_gaussian_cutoff_1;
    fl vina_gaussian_offset_2, vina_gaussian_width_2, vina_gaussian_cutoff_2;
    fl vina_repulsion_offset, vina_repulsion_cutoff;
    fl vina_hydrophobic_good, vina_hydrophobic_bad, vina_hydrophobic_cutoff;
    fl vina_non_dir_h_bond_good, vina_non_dir_h_bond_bad, vina_non_dir_h_bond_cutoff;
    fl vinardo_gaussian_offset, vinardo_gaussian_width, vinardo_gaussian_cutoff;
    fl vinardo_repulsion_offset, vinardo_repulsion_cutoff;
    fl vinardo_hydrophobic_good, vinardo_hydrophobic_bad, vinardo_hydrophobic_cutoff;
    fl vinardo_non_dir_h_bond_good, vinardo_non_dir_h_bond_bad, vinardo_non_dir_h_bond_cutoff;
    fl linearattraction_cutoff;  // shared by all scoring functions

} scoring_function_cuda_t;

/* scoring function related end */
