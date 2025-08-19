#ifndef VINA_MACROS_H
#define VINA_MACROS_H

#define VINA_FOR_IN(i, v) \
    for (std::size_t VINA_MACROS_TMP = (v).size(), (i) = 0; (i) < VINA_MACROS_TMP; ++(i))
#define VINA_FOR(i, n) \
    for (std::size_t VINA_MACROS_TMP = (n), (i) = 0; (i) < VINA_MACROS_TMP; ++(i))
#define VINA_U_FOR(i, n) for (unsigned VINA_MACROS_TMP = (n), (i) = 0; (i) < VINA_MACROS_TMP; ++(i))

#define VINA_RANGE(i, a, b) \
    for (std::size_t VINA_MACROS_TMP = (b), (i) = (a); (i) < VINA_MACROS_TMP; ++(i))
#define VINA_U_RANGE(i, a, b) \
    for (unsigned VINA_MACROS_TMP = (b), (i) = (a); (i) < VINA_MACROS_TMP; ++(i))
#define VINA_I_RANGE(i, a, b) \
    for (int VINA_MACROS_TMP = (b), (i) = (a); (i) < VINA_MACROS_TMP; ++(i))

#define VINA_LOOP_CONST(t, i, v) for (t::const_iterator(i) = (v).begin(); (i) != (v).end(); ++(i))
#define VINA_LOOP(t, i, v) for (t::iterator(i) = (v).begin(); (i) != (v).end(); ++(i))

#define VINA_SHOW(x)                                  \
    do {                                              \
        std::cout << #x << " = " << (x) << std::endl; \
    } while (false)
#define VINA_SHOW_FAST(x)                        \
    do {                                         \
        std::cout << #x << " = " << (x) << '\n'; \
    } while (false)
#define VINA_ESHOW(x)                            \
    do {                                         \
        std::cerr << #x << " = " << (x) << '\n'; \
    } while (false)

#endif
