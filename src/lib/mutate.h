#ifndef VINA_MUTATE_H
#define VINA_MUTATE_H

#include "model.h"

// does not set model
void mutate_conf(conf& c, const model& m, fl amplitude, rng& generator);

#endif
