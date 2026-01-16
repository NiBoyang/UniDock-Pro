#include "non_cache.h"
#include "curl.h"

non_cache::non_cache(const model& m, const grid_dims& gd_, const precalculate* p_, fl slope_)
    : sgrid(m, szv_grid_dims(gd_), p_->cutoff_sqr()),
      gd(gd_),
      p(p_),
      slope(slope_) {}

// Initialize the mode flags
void non_cache::set_mode(bool pure_docking, bool similarity_searching, bool hybrid_mode) {
    m_pure_docking = pure_docking;
    m_similarity_searching = similarity_searching;
    m_hybrid_mode = hybrid_mode;
}

// Set up the reference ligand data (multiple reference ligands)
void non_cache::set_reference_ligands(const std::vector<model>& ref_ligs) {
    if (ref_ligs.empty()) {
        has_reference_ligand = false;
        num_reference_ligands = 0;
        return;
    }

    has_reference_ligand = true;
    num_reference_ligands = ref_ligs.size();
    ref_typed_atom_coords.clear();

    // Gather coordinates of reference ligand atoms by type, with reference ligand index
    atom_type::t atype = atom_type::XS;

    for (size_t ref_idx = 0; ref_idx < ref_ligs.size(); ++ref_idx) {
        const atomv& ligatoms = ref_ligs[ref_idx].get_atoms();
        for (const auto& a : ligatoms) {
            sz t = a.get(atype);
            ref_typed_atom_coords[t].push_back(std::make_pair(ref_idx, a.coords));
        }
    }
}

// LJ-like energy calculation (same as in vina.cpp)
fl non_cache::lj_soft_cutoff(double d, double radius, double LJ_A, double LJ_B) const {
    if (d <= radius) {
        double d2 = d * d + 1.5;  // Adjustment for a perfect slope
        double d6 = d2 * d2 * d2;
        double d12 = d6 * d6;
        return -((LJ_A / d12) + (LJ_B / d6));
    }
    else {
        return 0.0;
    }
}

// Derivative of the LJ-like energy
fl non_cache::lj_soft_cutoff_deriv(double d, double radius, double LJ_A, double LJ_B) const {
    if (d <= radius) {
        double d2 = d * d + 1.5;
        double d6 = d2 * d2 * d2;
        double d12 = d6 * d6;
        
        // Calculate dE/dd
        return (12.0 * LJ_A * d / (d12 * d2) + 6.0 * LJ_B * d / (d6 * d2));
    }
    else {
        return 0.0;
    }
}

// Apply reference ligand bias for a given atom (averaged over all reference ligands)
fl non_cache::apply_reference_ligand_bias(const atom& a, const vec& a_coords, vec& deriv, bool calc_deriv) const {
    if (!has_reference_ligand || num_reference_ligands == 0) return 0.0;

    const fl cutoff_sqr = p->cutoff_sqr();
    double radius = 1.54;

    // Get atom type
    sz t1 = a.get(atom_type::XS);
    switch (t1) {
        case XS_TYPE_G0:
        case XS_TYPE_G1:
        case XS_TYPE_G2:
        case XS_TYPE_G3:
            return 0.0;
        case XS_TYPE_C_H_CG0:
        case XS_TYPE_C_H_CG1:
        case XS_TYPE_C_H_CG2:
        case XS_TYPE_C_H_CG3:
            t1 = XS_TYPE_C_H;
            break;
        case XS_TYPE_C_P_CG0:
        case XS_TYPE_C_P_CG1:
        case XS_TYPE_C_P_CG2:
        case XS_TYPE_C_P_CG3:
            t1 = XS_TYPE_C_P;
            break;
    }

    // Get LJ parameters for this atom type
    double LJ_A = xs_lj(t1).LJ_A * 3.0 * m_reference_ligand_scale;
    double LJ_B = xs_lj(t1).LJ_B * 3.0 * m_reference_ligand_scale;

    // Find corresponding atom type in reference ligands
    auto it = ref_typed_atom_coords.find(t1);
    if (it == ref_typed_atom_coords.end()) {
        return 0.0;
    }

    // Calculate energy and derivative contribution from each reference ligand separately
    std::vector<fl> E_per_reflig(num_reference_ligands, 0.0);
    std::vector<vec> deriv_per_reflig(num_reference_ligands, vec(0, 0, 0));

    for (const auto& coord_pair : it->second) {
        size_t ref_idx = coord_pair.first;
        const vec& ref_coord = coord_pair.second;

        double d2 = vec_distance_sqr(ref_coord, a_coords);
        if (d2 > cutoff_sqr) continue;

        double d = std::sqrt(d2);
        fl dE = lj_soft_cutoff(d, radius, LJ_A, LJ_B);

        E_per_reflig[ref_idx] += dE;

        // Calculate and add derivatives if needed
        if (calc_deriv && dE != 0) {
            fl deriv_magnitude = lj_soft_cutoff_deriv(d, radius, LJ_A, LJ_B);
            vec deriv_direction = (a_coords - ref_coord) / d;  // Unit vector
            deriv_per_reflig[ref_idx] += deriv_magnitude * deriv_direction;
        }
    }

    // Average over all reference ligands
    fl bias_e = 0.0;
    vec avg_deriv(0, 0, 0);
    for (size_t i = 0; i < num_reference_ligands; ++i) {
        bias_e += E_per_reflig[i];
        avg_deriv += deriv_per_reflig[i];
    }
    fl inv_num_ref = 1.0 / static_cast<fl>(num_reference_ligands);
    bias_e *= inv_num_ref;
    avg_deriv = inv_num_ref * avg_deriv;

    if (calc_deriv) {
        deriv += avg_deriv;
    }

    return bias_e;
}

fl non_cache::eval(const model& m, fl v) const {  // clean up
    fl e = 0;
    const fl cutoff_sqr = p->cutoff_sqr();

    sz n = num_atom_types(atom_type::XS);

    if (m_similarity_searching) {
        VINA_FOR(i, m.num_movable_atoms()) {
            fl this_e = 0;
            fl out_of_bounds_penalty = 0;
            const atom& a = m.atoms[i];
            sz t1 = a.get(atom_type::XS);
            if (t1 >= n) continue;
            // Skip dummy atoms
            switch (t1) {
                case XS_TYPE_G0:
                case XS_TYPE_G1:
                case XS_TYPE_G2:
                case XS_TYPE_G3:
                    continue;
                case XS_TYPE_C_H_CG0:
                case XS_TYPE_C_H_CG1:
                case XS_TYPE_C_H_CG2:
                case XS_TYPE_C_H_CG3:
                    t1 = XS_TYPE_C_H;
                    break;
                case XS_TYPE_C_P_CG0:
                case XS_TYPE_C_P_CG1:
                case XS_TYPE_C_P_CG2:
                case XS_TYPE_C_P_CG3:
                    t1 = XS_TYPE_C_P;
                    break;
            }

            const vec& a_coords = m.coords[i];
            vec adjusted_a_coords;
            adjusted_a_coords = a_coords;
            
            // Handle out of bounds
            VINA_FOR_IN(j, gd) {
                if (gd[j].n_voxels > 0) {
                    if (a_coords[j] < gd[j].begin) {
                        adjusted_a_coords[j] = gd[j].begin;
                        out_of_bounds_penalty += std::abs(a_coords[j] - gd[j].begin);
                    } else if (a_coords[j] > gd[j].end) {
                        adjusted_a_coords[j] = gd[j].end;
                        out_of_bounds_penalty += std::abs(a_coords[j] - gd[j].end);
                    }
                }
            }
            out_of_bounds_penalty *= slope;

            // Apply reference ligand bias
            vec dummy_deriv; // Not used in eval
            this_e += apply_reference_ligand_bias(a, adjusted_a_coords, dummy_deriv, false);

            curl(this_e, v);
            e += this_e + out_of_bounds_penalty;
        }
        return e;
    }

    VINA_FOR(i, m.num_movable_atoms()) {
        fl this_e = 0;
        fl out_of_bounds_penalty = 0;
        const atom& a = m.atoms[i];
        sz t1 = a.get(atom_type::XS);
        if (t1 >= n) continue;
        switch (t1) {
            case XS_TYPE_G0:
            case XS_TYPE_G1:
            case XS_TYPE_G2:
            case XS_TYPE_G3:
                continue;
            case XS_TYPE_C_H_CG0:
            case XS_TYPE_C_H_CG1:
            case XS_TYPE_C_H_CG2:
            case XS_TYPE_C_H_CG3:
                t1 = XS_TYPE_C_H;
                break;
            case XS_TYPE_C_P_CG0:
            case XS_TYPE_C_P_CG1:
            case XS_TYPE_C_P_CG2:
            case XS_TYPE_C_P_CG3:
                t1 = XS_TYPE_C_P;
                break;
        }

        const vec& a_coords = m.coords[i];
        vec adjusted_a_coords;
        adjusted_a_coords = a_coords;
        VINA_FOR_IN(j, gd) {
            if (gd[j].n_voxels > 0) {
                if (a_coords[j] < gd[j].begin) {
                    adjusted_a_coords[j] = gd[j].begin;
                    out_of_bounds_penalty += std::abs(a_coords[j] - gd[j].begin);
                } else if (a_coords[j] > gd[j].end) {
                    adjusted_a_coords[j] = gd[j].end;
                    out_of_bounds_penalty += std::abs(a_coords[j] - gd[j].end);
                }
            }
        }
        out_of_bounds_penalty *= slope;

        const szv& possibilities = sgrid.possibilities(adjusted_a_coords);
        fl receptor_e = 0;

        VINA_FOR_IN(possibilities_j, possibilities) {
            const sz j = possibilities[possibilities_j];
            const atom& b = m.grid_atoms[j];
            sz t2 = b.get(atom_type::XS);
            if (t2 >= n) continue;
            vec r_ba;
            r_ba = adjusted_a_coords - b.coords;  // FIXME why b-a and not a-b ?
            fl r2 = sqr(r_ba);
            if (r2 < cutoff_sqr) {
                sz type_pair_index = get_type_pair_index(atom_type::XS, a, b);
                receptor_e += p->eval_fast(type_pair_index, r2);
            }
        }

        if (m_hybrid_mode) {
            this_e += receptor_e;

            // Add reference ligand bias
            vec dummy_deriv; // Not used in eval
            this_e += apply_reference_ligand_bias(a, adjusted_a_coords, dummy_deriv, false);
        } else {
            this_e += receptor_e;
        }

        // Bias terms removed

        curl(this_e, v);

        e += this_e + out_of_bounds_penalty;
    }
    return e;
}

bool non_cache::within(const model& m, fl margin) const {
    VINA_FOR(i, m.num_movable_atoms()) {
        if (m.atoms[i].is_hydrogen()) continue;
        const vec& a_coords = m.coords[i];
        VINA_FOR_IN(j, gd)
        if (gd[j].n_voxels > 0)
            if (a_coords[j] < gd[j].begin - margin || a_coords[j] > gd[j].end + margin)
                return false;
    }
    return true;
}

fl non_cache::eval_deriv(model& m, fl v) const {  // clean up
    fl e = 0;
    const fl cutoff_sqr = p->cutoff_sqr();

    sz n = num_atom_types(atom_type::XS);

    if (m_similarity_searching) {
        VINA_FOR(i, m.num_movable_atoms()) {
            fl this_e = 0;
            vec deriv(0, 0, 0);
            vec out_of_bounds_deriv(0, 0, 0);
            fl out_of_bounds_penalty = 0;
            const atom& a = m.atoms[i];
            sz t1 = a.get(atom_type::XS);
            if (t1 >= n) {
                m.minus_forces[i].assign(0);
                continue;
            }
            // Skip dummy atoms
            switch (t1) {
                case XS_TYPE_G0:
                case XS_TYPE_G1:
                case XS_TYPE_G2:
                case XS_TYPE_G3:
                    continue;
                case XS_TYPE_C_H_CG0:
                case XS_TYPE_C_H_CG1:
                case XS_TYPE_C_H_CG2:
                case XS_TYPE_C_H_CG3:
                    t1 = XS_TYPE_C_H;
                    break;
                case XS_TYPE_C_P_CG0:
                case XS_TYPE_C_P_CG1:
                case XS_TYPE_C_P_CG2:
                case XS_TYPE_C_P_CG3:
                    t1 = XS_TYPE_C_P;
                    break;
            }
            const vec& a_coords = m.coords[i];
            vec adjusted_a_coords;
            adjusted_a_coords = a_coords;
            
            // Handle out of bounds
            VINA_FOR_IN(j, gd) {
                if (gd[j].n_voxels > 0) {
                    if (a_coords[j] < gd[j].begin) {
                        adjusted_a_coords[j] = gd[j].begin;
                        out_of_bounds_deriv[j] = -1;
                        out_of_bounds_penalty += std::abs(a_coords[j] - gd[j].begin);
                    } else if (a_coords[j] > gd[j].end) {
                        adjusted_a_coords[j] = gd[j].end;
                        out_of_bounds_deriv[j] = 1;
                        out_of_bounds_penalty += std::abs(a_coords[j] - gd[j].end);
                    }
                }
            }
            out_of_bounds_penalty *= slope;
            out_of_bounds_deriv *= slope;

            // Apply reference ligand bias with derivatives
            this_e += apply_reference_ligand_bias(a, adjusted_a_coords, deriv, true);

            curl(this_e, deriv, v);
            m.minus_forces[i] = deriv + out_of_bounds_deriv;

            e += this_e + out_of_bounds_penalty;
        }
        return e;
    }

    VINA_FOR(i, m.num_movable_atoms()) {
        fl this_e = 0;
        vec deriv(0, 0, 0);
        vec out_of_bounds_deriv(0, 0, 0);
        fl out_of_bounds_penalty = 0;
        const atom& a = m.atoms[i];
        sz t1 = a.get(atom_type::XS);
        if (t1 >= n) {
            m.minus_forces[i].assign(0);
            continue;
        }
        switch (t1) {
            case XS_TYPE_G0:
            case XS_TYPE_G1:
            case XS_TYPE_G2:
            case XS_TYPE_G3:
                continue;
            case XS_TYPE_C_H_CG0:
            case XS_TYPE_C_H_CG1:
            case XS_TYPE_C_H_CG2:
            case XS_TYPE_C_H_CG3:
                t1 = XS_TYPE_C_H;
                break;
            case XS_TYPE_C_P_CG0:
            case XS_TYPE_C_P_CG1:
            case XS_TYPE_C_P_CG2:
            case XS_TYPE_C_P_CG3:
                t1 = XS_TYPE_C_P;
                break;
        }
        const vec& a_coords = m.coords[i];
        vec adjusted_a_coords;
        adjusted_a_coords = a_coords;
        VINA_FOR_IN(j, gd) {
            if (gd[j].n_voxels > 0) {
                if (a_coords[j] < gd[j].begin) {
                    adjusted_a_coords[j] = gd[j].begin;
                    out_of_bounds_deriv[j] = -1;
                    out_of_bounds_penalty += std::abs(a_coords[j] - gd[j].begin);
                } else if (a_coords[j] > gd[j].end) {
                    adjusted_a_coords[j] = gd[j].end;
                    out_of_bounds_deriv[j] = 1;
                    out_of_bounds_penalty += std::abs(a_coords[j] - gd[j].end);
                }
            }
        }
        out_of_bounds_penalty *= slope;
        out_of_bounds_deriv *= slope;

        const szv& possibilities = sgrid.possibilities(adjusted_a_coords);
        fl receptor_e = 0;
        vec receptor_deriv(0, 0, 0);

        VINA_FOR_IN(possibilities_j, possibilities) {
            const sz j = possibilities[possibilities_j];
            const atom& b = m.grid_atoms[j];
            sz t2 = b.get(atom_type::XS);
            if (t2 >= n) continue;
            vec r_ba;
            r_ba = adjusted_a_coords - b.coords;  // FIXME why b-a and not a-b ?
            fl r2 = sqr(r_ba);
            if (r2 < cutoff_sqr) {
                sz type_pair_index = get_type_pair_index(atom_type::XS, a, b);
                pr e_dor = p->eval_deriv(type_pair_index, r2);
                receptor_e += e_dor.first;
                receptor_deriv += e_dor.second * r_ba;
            }
        }

        if (m_hybrid_mode) {
            this_e += receptor_e;
            deriv += receptor_deriv;

            // Add reference ligand bias with derivatives
            vec ref_deriv(0, 0, 0);
            fl ref_e = apply_reference_ligand_bias(a, adjusted_a_coords, ref_deriv, true);
            this_e += ref_e;
            deriv += ref_deriv;
        } else {
            this_e += receptor_e;
            deriv += receptor_deriv;
        }
        curl(this_e, deriv, v);
        m.minus_forces[i] = deriv + out_of_bounds_deriv;
        e += this_e + out_of_bounds_penalty;
    }
    return e;
}

fl non_cache::eval_intra(model& m, fl v) const {  // clean up
    fl e = 0;
    const fl cutoff_sqr = p->cutoff_sqr();

    sz n = num_atom_types(atom_type::XS);

    VINA_FOR(i, m.num_movable_atoms()) {
        if (m.is_atom_in_ligand(i)) continue;  // we only want flex-rigid interactions
        fl this_e = 0;
        fl out_of_bounds_penalty = 0;
        const atom& a = m.atoms[i];
        sz t1 = a.get(atom_type::XS);
        if (t1 >= n) continue;
        switch (t1) {
            case XS_TYPE_G0:
            case XS_TYPE_G1:
            case XS_TYPE_G2:
            case XS_TYPE_G3:
                continue;
            case XS_TYPE_C_H_CG0:
            case XS_TYPE_C_H_CG1:
            case XS_TYPE_C_H_CG2:
            case XS_TYPE_C_H_CG3:
                t1 = XS_TYPE_C_H;
                break;
            case XS_TYPE_C_P_CG0:
            case XS_TYPE_C_P_CG1:
            case XS_TYPE_C_P_CG2:
            case XS_TYPE_C_P_CG3:
                t1 = XS_TYPE_C_P;
                break;
        }

        const vec& a_coords = m.coords[i];
        vec adjusted_a_coords;
        adjusted_a_coords = a_coords;
        VINA_FOR_IN(j, gd) {
            if (gd[j].n_voxels > 0) {
                if (a_coords[j] < gd[j].begin) {
                    adjusted_a_coords[j] = gd[j].begin;
                    out_of_bounds_penalty += std::abs(a_coords[j] - gd[j].begin);
                } else if (a_coords[j] > gd[j].end) {
                    adjusted_a_coords[j] = gd[j].end;
                    out_of_bounds_penalty += std::abs(a_coords[j] - gd[j].end);
                }
            }
        }
        out_of_bounds_penalty *= slope;

        const szv& possibilities = sgrid.possibilities(adjusted_a_coords);

        VINA_FOR_IN(possibilities_j, possibilities) {
            const sz j = possibilities[possibilities_j];
            const atom& b = m.grid_atoms[j];
            sz t2 = b.get(atom_type::XS);
            if (t2 >= n) continue;
            vec r_ba;
            r_ba = adjusted_a_coords - b.coords;  // FIXME why b-a and not a-b ?
            fl r2 = sqr(r_ba);
            if (r2 < cutoff_sqr) {
                sz type_pair_index = get_type_pair_index(atom_type::XS, a, b);
                this_e += p->eval_fast(type_pair_index, r2);
            }
        }
        curl(this_e, v);
        e += this_e + out_of_bounds_penalty;
    }
    return e;
}

// Evaluate only reference ligand contribution (for hybrid mode energy decomposition)
fl non_cache::eval_reflig_only(const model& m, fl v) const {
    if (!has_reference_ligand || !m_hybrid_mode) return 0.0;

    fl e = 0;
    sz n = num_atom_types(atom_type::XS);

    VINA_FOR(i, m.num_movable_atoms()) {
        fl this_e = 0;
        const atom& a = m.atoms[i];
        sz t1 = a.get(atom_type::XS);
        if (t1 >= n) continue;
        switch (t1) {
            case XS_TYPE_G0:
            case XS_TYPE_G1:
            case XS_TYPE_G2:
            case XS_TYPE_G3:
                continue;
            case XS_TYPE_C_H_CG0:
            case XS_TYPE_C_H_CG1:
            case XS_TYPE_C_H_CG2:
            case XS_TYPE_C_H_CG3:
                t1 = XS_TYPE_C_H;
                break;
            case XS_TYPE_C_P_CG0:
            case XS_TYPE_C_P_CG1:
            case XS_TYPE_C_P_CG2:
            case XS_TYPE_C_P_CG3:
                t1 = XS_TYPE_C_P;
                break;
        }

        const vec& a_coords = m.coords[i];
        vec adjusted_a_coords = a_coords;

        // Handle out of bounds
        VINA_FOR_IN(j, gd) {
            if (gd[j].n_voxels > 0) {
                if (a_coords[j] < gd[j].begin) {
                    adjusted_a_coords[j] = gd[j].begin;
                } else if (a_coords[j] > gd[j].end) {
                    adjusted_a_coords[j] = gd[j].end;
                }
            }
        }

        // Calculate reference ligand bias only
        vec dummy_deriv; // Not used in eval
        fl ref_e = apply_reference_ligand_bias(a, adjusted_a_coords, dummy_deriv, false);
        this_e = ref_e;

        curl(this_e, v);
        e += this_e;
    }
    return e;
}

// add to make get_grids() work
std::vector<grid> non_cache::get_grids() const {
    assert(false);  // This function should not be called!
    std::vector<grid> g;
    return g;
};

int non_cache::get_atu() const {
    assert(false);  // This function should not be called!
    return 0;
}

float non_cache::get_slope() const {
    assert(false);  // This function should not be called!
    return 0;
}