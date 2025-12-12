#include "coords.h"

fl rmsd_upper_bound(const vecv& a, const vecv& b) {
    VINA_CHECK(a.size() == b.size());
    fl acc = 0;
    VINA_FOR_IN(i, a)
    acc += vec_distance_sqr(a[i], b[i]);
    return (a.size() > 0) ? std::sqrt(acc / a.size()) : 0;
}

std::pair<sz, fl> find_closest(const vecv& a, const output_container& b) {
    std::pair<sz, fl> tmp(b.size(), max_fl);
    VINA_FOR_IN(i, b) {
        fl res = rmsd_upper_bound(a, b[i].coords);
        if (i == 0 || res < tmp.second) tmp = std::pair<sz, fl>(i, res);
    }
    return tmp;
}

void add_to_output_container(output_container& out, const output_type& t, fl min_rmsd,
                             sz max_size) {
    if (out.size() > 0 && t.coords.size() != out[0].coords.size()) {
        printf(
            "WARNING: in add_to_output_container, adding the %luth ligand\nt.coords.size()=%lu, "
            "out[0].coords.size()=%lu\n",
            out.size(), t.coords.size(), out[0].coords.size());
        return;
    }
    if (t.e < epsilon_fl && t.e > -epsilon_fl) {  // the energy is precisely 0, error
        return;
    }
    std::pair<sz, fl> closest_rmsd = find_closest(t.coords, out);
    if (closest_rmsd.first < out.size()
        && closest_rmsd.second < min_rmsd) {    // have a very similar one
        if (t.e < out[closest_rmsd.first].e) {  // the new one is better, apparently
            out[closest_rmsd.first] = t;        // FIXME? slow
        }
    } else {  // nothing similar
        if (out.size() < max_size)
            out.push_back(new output_type(t));  // the last one had the worst energy - replacing
        else if (!out.empty() && t.e < out.back().e)  // FIXME? - just changed
            out.back() = t;                           // FIXME? slow
    }
    out.sort();
}
