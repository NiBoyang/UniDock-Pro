#ifndef VINA_NON_CACHE_H
#define VINA_NON_CACHE_H

#include "igrid.h"
#include "szv_grid.h"
#include "precalculate.h"
#include <unordered_map>

struct non_cache : public igrid {
    non_cache() {}
    non_cache(const model& m, const grid_dims& gd_, const precalculate* p_, fl slope_);
    virtual fl eval(const model& m, fl v) const;  // needs m.coords // clean up
    virtual fl eval_intra(model& m, fl v) const;
    virtual fl eval_deriv(model& m, fl v) const;  // needs m.coords, sets m.minus_forces // clean up
    std::vector<grid> get_grids() const;
    int get_atu() const;
    float get_slope() const;
    bool within(const model& m, fl margin = -0.0001) const;
    fl slope;
    grid_dims get_gd() const { return gd; }
    
    // Reference ligand support (multiple reference ligands)
    void set_reference_ligands(const std::vector<model>& ref_ligs);
    void set_mode(bool pure_docking, bool similarity_searching, bool hybrid_mode);
    void set_reference_ligand_scale(double scale) { m_reference_ligand_scale = scale; }

    // Evaluate only reference ligand contribution (for hybrid mode energy decomposition)
    fl eval_reflig_only(const model& m, fl v) const;

private:
    // Applies the reference ligand bias to a given atom (averaged over all reference ligands)
    fl apply_reference_ligand_bias(const atom& a, const vec& a_coords, vec& deriv, bool calc_deriv) const;

    // LJ-like model for similarity bias
    fl lj_soft_cutoff(double d, double radius, double LJ_A, double LJ_B) const;
    fl lj_soft_cutoff_deriv(double d, double radius, double LJ_A, double LJ_B) const;

    szv_grid sgrid;
    grid_dims gd;
    const precalculate* p;

    // Reference ligand data (multiple reference ligands)
    bool has_reference_ligand = false;
    size_t num_reference_ligands = 0;
    // Stores coordinates per reference ligand: map[atom_type] -> vector of (ref_lig_idx, coords)
    std::unordered_map<sz, std::vector<std::pair<size_t, vec>>> ref_typed_atom_coords;

    // Docking mode flags
    bool m_pure_docking = true;
    bool m_similarity_searching = false;
    bool m_hybrid_mode = false;

    // Scale for reference ligand LJ parameters
    double m_reference_ligand_scale = 1.0;
};

#endif
