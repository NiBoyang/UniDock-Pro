#ifndef VINA_RANDOM_H
#define VINA_RANDOM_H

#include <random>
#include <boost/random.hpp>
#include "common.h"

typedef boost::mt19937 rng;

fl random_fl(fl a, fl b, rng& generator);             // expects a < b, returns rand in [a, b]
fl random_normal(fl mean, fl sigma, rng& generator);  // expects sigma >= 0
int random_int(int a, int b, rng& generator);         // expects a <= b, returns rand in [a, b]
sz random_sz(sz a, sz b, rng& generator);             // expects a <= b, returns rand in [a, b]
vec random_inside_sphere(
    rng& generator);  // returns a random vec inside the sphere centered at 0 with radius 1
vec random_in_box(const vec& corner1, const vec& corner2,
                  rng& generator);  // expects corner1[i] < corner2[i]
int auto_seed();                    // make seed from PID and time

#endif
