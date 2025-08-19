#ifndef VINA_COORDS_H
#define VINA_COORDS_H

#include "conf.h"
#include "atom.h"  // for atomv

fl rmsd_upper_bound(const vecv& a, const vecv& b);
std::pair<sz, fl> find_closest(const vecv& a, const output_container& b);
void add_to_output_container(output_container& out, const output_type& t, fl min_rmsd, sz max_size);

#endif
