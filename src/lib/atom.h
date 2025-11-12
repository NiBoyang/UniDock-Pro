#ifndef VINA_ATOM_H
#define VINA_ATOM_H

#include "atom_base.h"

struct atom_index {
    sz i;
    bool in_grid;
    atom_index() : i(max_sz), in_grid(false) {}
    atom_index(sz i_, bool in_grid_) : i(i_), in_grid(in_grid_) {}

private:
    friend class boost::serialization::access;
    template <class Archive> void serialize(Archive& ar, const unsigned version) {
        ar & i;
        ar & in_grid;
    }
};

inline bool operator==(const atom_index& i, const atom_index& j) {
    return i.i == j.i && i.in_grid == j.in_grid;
}

struct bond {
    atom_index connected_atom_index;
    fl length;
    bool rotatable;
    bond() : length(0), rotatable(false) {}
    bond(const atom_index& connected_atom_index_, fl length_, bool rotatable_)
        : connected_atom_index(connected_atom_index_), length(length_), rotatable(rotatable_) {}

private:
    friend class boost::serialization::access;
    template <class Archive> void serialize(Archive& ar, const unsigned version) {
        ar & connected_atom_index;
        ar & length;
        ar & rotatable;
    }
};

struct atom : public atom_base {
    vec coords;
    std::vector<bond> bonds;
    atom() : coords(max_vec) {}

private:
    friend class boost::serialization::access;
    template <class Archive> void serialize(Archive& ar, const unsigned version) {
        ar& boost::serialization::base_object<atom_base>(*this);
        ar & coords;
        ar & bonds;
    }
};

typedef std::vector<atom> atomv;

#endif
