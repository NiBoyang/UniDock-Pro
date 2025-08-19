#include "conf_independent.h"
#include "model.h"

conf_independent_inputs::operator flv() const {
    flv tmp;
    tmp.push_back(num_tors);
    tmp.push_back(num_rotors);
    tmp.push_back(num_heavy_atoms);
    tmp.push_back(num_hydrophobic_atoms);
    tmp.push_back(ligand_max_num_h_bonds);
    tmp.push_back(num_ligands);
    tmp.push_back(ligand_lengths_sum);
    return tmp;
}

unsigned conf_independent_inputs::num_bonded_heavy_atoms(const model& m, const atom_index& i)
    const {  // FIXME? - could be static, but I don't feel like declaring function friends
    unsigned acc = 0;
    const std::vector<bond>& bonds = m.get_atom(i).bonds;

    VINA_FOR_IN(j, bonds) {
        const bond& b = bonds[j];
        const atom& a = m.get_atom(b.connected_atom_index);
        if (!a.is_hydrogen()) ++acc;
    }

    return acc;
}

// FIXME? - could be static, but I don't feel like declaring function friends
unsigned conf_independent_inputs::atom_rotors(const model& m, const atom_index& i)
    const {  // the number of rotatable bonds to heavy ligand atoms
    unsigned acc = 0;
    const std::vector<bond>& bonds = m.get_atom(i).bonds;

    VINA_FOR_IN(j, bonds) {
        const bond& b = bonds[j];
        const atom& a = m.get_atom(b.connected_atom_index);
        if (b.rotatable && !a.is_hydrogen()
            && num_bonded_heavy_atoms(m, b.connected_atom_index) > 1) {  // not counting CH_3, etc
            ++acc;
        }
    }

    return acc;
}

conf_independent_inputs::conf_independent_inputs()
    : num_tors(0),
      num_rotors(0),
      num_heavy_atoms(0),
      num_hydrophobic_atoms(0),
      ligand_max_num_h_bonds(0),
      num_ligands(0),
      ligand_lengths_sum(0) {}

conf_independent_inputs::conf_independent_inputs(const model& m) {
    torsdof = 0;
    num_tors = 0;
    num_rotors = 0;
    num_heavy_atoms = 0;
    num_hydrophobic_atoms = 0;
    ligand_max_num_h_bonds = 0;
    num_ligands = m.num_ligands();
    ligand_lengths_sum = 0;

    VINA_FOR(i, num_ligands) {
        const ligand& lig = m.get_ligand(i);
        ligand_lengths_sum += m.ligand_length(i);
        torsdof += lig.degrees_of_freedom;

        VINA_RANGE(j, lig.begin, lig.end) {
            const atom& a = m.get_atom(j);

            if (a.el != EL_TYPE_H) {
                unsigned ar = atom_rotors(m, atom_index(j, false));
                num_tors += 0.5 * ar;
                if (ar > 2)
                    num_rotors += 0.5;
                else
                    num_rotors += 0.5 * ar;

                if (xs_is_hydrophobic(a.xs)) ++num_hydrophobic_atoms;

                if (xs_is_acceptor(a.xs) || xs_is_donor(a.xs)) ++ligand_max_num_h_bonds;

                ++num_heavy_atoms;
            }
        }
    }
}


inline fl read_iterator(flv::const_iterator& i) {
    fl x = *i;
    ++i;
    return x;
}

fl conf_smooth_div(fl x, fl y) {
    if (std::abs(x) < epsilon_fl) return 0;
    if (std::abs(y) < epsilon_fl)
        return ((x * y > 0) ? max_fl : -max_fl);  // FIXME I hope -max_fl does not become NaN
    return x / y;
}

// Vina conf_independent functions
fl num_tors_div::eval(const conf_independent_inputs& in, fl x, flv::const_iterator& i) {
    fl weight = 0.1 * (read_iterator(i) + 1);  // weight is in [0..0.2]
    return conf_smooth_div(x, 1 + weight * in.num_tors / 5.0);
}

