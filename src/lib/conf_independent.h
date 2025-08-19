#ifndef VINA_CONF_INDEPENDENT_H
#define VINA_CONF_INDEPENDENT_H

#include <stdlib.h>
#include "atom.h"
#include "common.h"

// Forward declaration
struct model;

class conf_independent_inputs {
public:
    fl torsdof;  // from TORSDOF keyword in pdbqt file
    fl num_tors;
    fl num_rotors;
    fl num_heavy_atoms;
    fl num_hydrophobic_atoms;
    fl ligand_max_num_h_bonds;
    fl num_ligands;
    fl ligand_lengths_sum;

    operator flv() const;
    conf_independent_inputs();
    conf_independent_inputs(const model& m);

private:
    unsigned num_bonded_heavy_atoms(const model& m, const atom_index& i)
        const;  // FIXME? - could be static, but I don't feel like declaring function friends
    unsigned atom_rotors(const model& m, const atom_index& i)
        const;  // the number of rotatable bonds to heavy ligand atoms
};

// Conf independent
class ConfIndependent {
public:
    virtual ~ConfIndependent() {}
    virtual fl eval(const conf_independent_inputs& in, fl x, flv::const_iterator& i) { return 0; };
};

// Vina
class num_tors_div : public ConfIndependent {
public:
    num_tors_div() {}
    sz size() const { return 1; }
    fl eval(const conf_independent_inputs& in, fl x, flv::const_iterator& i);
};

#endif
