#ifndef VINA_ATOM_BASE_H
#define VINA_ATOM_BASE_H

#include "atom_type.h"

struct atom_base : public atom_type {
    fl charge;
    atom_base() : charge(0) {}

private:
    friend class boost::serialization::access;
    template <class Archive> void serialize(Archive& ar, const unsigned version) {
        ar& boost::serialization::base_object<atom_type>(*this);
        ar & charge;
    }
};

#endif
