#pragma once

#include <stdexcept>
template <typename T>
void check(T result, char const *const func, const char *const file, int const line) {
    if (result) {
        printf("CUDA error at %s:%d code=%d(%s) \"%s\" \n", file, line,
               static_cast<unsigned int>(result), cudaGetErrorName(result), func);
        throw std::runtime_error("CUDA Runtime Error");
    }
}
#define checkCUDA(val) check((val), #val, __FILE__, __LINE__)

template <size_t MAX_NUM_OF_LIG_TORSION_VALUE,
          size_t MAX_NUM_OF_FLEX_TORSION_VALUE,
          size_t MAX_NUM_OF_RIGID_VALUE,
          size_t MAX_NUM_OF_ATOMS_VALUE,
          size_t MAX_NUM_OF_LIG_PAIRS_VALUE,
          size_t MAX_P_DATA_M_DATA_SIZE_VALUE =
              MAX_NUM_OF_ATOMS_VALUE * (MAX_NUM_OF_ATOMS_VALUE + 1) / 2>
struct BaseConfig {
    static constexpr float TOLERANCE = 1e-16f;
    // Provide underscored aliases so templated code expecting Config::VALUE_
    // continues to compile without modification.
    static constexpr float TOLERANCE_ = TOLERANCE;
    static constexpr size_t MAX_NUM_OF_EVERY_M_DATA_ELEMENT = 512;
    static constexpr size_t MAX_NUM_OF_EVERY_M_DATA_ELEMENT_ =
        MAX_NUM_OF_EVERY_M_DATA_ELEMENT;
    static constexpr size_t MAX_M_DATA_MI = 16;
    static constexpr size_t MAX_M_DATA_MI_ = MAX_M_DATA_MI;
    static constexpr size_t MAX_M_DATA_MJ = 16;
    static constexpr size_t MAX_M_DATA_MJ_ = MAX_M_DATA_MJ;
    static constexpr size_t MAX_M_DATA_MK = 16;
    static constexpr size_t MAX_M_DATA_MK_ = MAX_M_DATA_MK;
    static constexpr size_t MAX_NUM_OF_TOTAL_M_DATA =
        MAX_M_DATA_MI * MAX_M_DATA_MJ * MAX_M_DATA_MK *
        MAX_NUM_OF_EVERY_M_DATA_ELEMENT;
    static constexpr size_t MAX_NUM_OF_TOTAL_M_DATA_ = MAX_NUM_OF_TOTAL_M_DATA;

    static constexpr size_t MAX_NUM_OF_LIG_TORSION =
        MAX_NUM_OF_LIG_TORSION_VALUE;
    static constexpr size_t MAX_NUM_OF_LIG_TORSION_ = MAX_NUM_OF_LIG_TORSION;
    static constexpr size_t MAX_NUM_OF_FLEX_TORSION =
        MAX_NUM_OF_FLEX_TORSION_VALUE;
    static constexpr size_t MAX_NUM_OF_FLEX_TORSION_ = MAX_NUM_OF_FLEX_TORSION;
    static constexpr size_t MAX_NUM_OF_RIGID = MAX_NUM_OF_RIGID_VALUE;
    static constexpr size_t MAX_NUM_OF_RIGID_ = MAX_NUM_OF_RIGID;
    static constexpr size_t MAX_NUM_OF_ATOMS = MAX_NUM_OF_ATOMS_VALUE;
    static constexpr size_t MAX_NUM_OF_ATOMS_ = MAX_NUM_OF_ATOMS;
    static constexpr size_t SIZE_OF_MOLEC_STRUC =
        ((3 + 4 + MAX_NUM_OF_LIG_TORSION + MAX_NUM_OF_FLEX_TORSION + 1) *
         sizeof(float));
    static constexpr size_t SIZE_OF_MOLEC_STRUC_ = SIZE_OF_MOLEC_STRUC;
    static constexpr size_t SIZE_OF_CHANGE_STRUC =
        ((3 + 3 + MAX_NUM_OF_LIG_TORSION + MAX_NUM_OF_FLEX_TORSION + 1) *
         sizeof(float));
    static constexpr size_t SIZE_OF_CHANGE_STRUC_ = SIZE_OF_CHANGE_STRUC;
    static constexpr size_t MAX_HESSIAN_MATRIX_D_SIZE =
        ((6 + MAX_NUM_OF_LIG_TORSION + MAX_NUM_OF_FLEX_TORSION) *
         (6 + MAX_NUM_OF_LIG_TORSION + MAX_NUM_OF_FLEX_TORSION + 1) / 2);
    static constexpr size_t MAX_HESSIAN_MATRIX_D_SIZE_ =
        MAX_HESSIAN_MATRIX_D_SIZE;
    static constexpr size_t MAX_NUM_OF_LIG_PAIRS = MAX_NUM_OF_LIG_PAIRS_VALUE;
    static constexpr size_t MAX_NUM_OF_LIG_PAIRS_ = MAX_NUM_OF_LIG_PAIRS;
    static constexpr size_t MAX_NUM_OF_BFGS_STEPS = 64;
    static constexpr size_t MAX_NUM_OF_BFGS_STEPS_ = MAX_NUM_OF_BFGS_STEPS;
    static constexpr size_t MAX_NUM_OF_RANDOM_MAP = 1000;  // not too large
    static constexpr size_t MAX_NUM_OF_RANDOM_MAP_ = MAX_NUM_OF_RANDOM_MAP;
    static constexpr size_t GRIDS_SIZE = 37;  // larger than vina1.1
    static constexpr size_t GRIDS_SIZE_ = GRIDS_SIZE;

    static constexpr size_t MAX_NUM_OF_GRID_MI = 128;  // 55
    static constexpr size_t MAX_NUM_OF_GRID_MI_ = MAX_NUM_OF_GRID_MI;
    static constexpr size_t MAX_NUM_OF_GRID_MJ = 128;  // 55
    static constexpr size_t MAX_NUM_OF_GRID_MJ_ = MAX_NUM_OF_GRID_MJ;
    static constexpr size_t MAX_NUM_OF_GRID_MK = 128;  // 81
    static constexpr size_t MAX_NUM_OF_GRID_MK_ = MAX_NUM_OF_GRID_MK;
    static constexpr size_t MAX_NUM_OF_GRID_POINT = 729000;
    static constexpr size_t MAX_NUM_OF_GRID_POINT_ = MAX_NUM_OF_GRID_POINT;

    static constexpr size_t GRID_MI = 65;  // 55
    static constexpr size_t GRID_MJ = 71;  // 55
    static constexpr size_t GRID_MK = 61;  // 81
    static constexpr size_t MAX_PRECAL_NUM_ATOM = 30;
    static constexpr size_t MAX_PRECAL_NUM_ATOM_ = MAX_PRECAL_NUM_ATOM;
    static constexpr size_t MAX_P_DATA_M_DATA_SIZE = MAX_P_DATA_M_DATA_SIZE_VALUE;
    static constexpr size_t MAX_P_DATA_M_DATA_SIZE_ = MAX_P_DATA_M_DATA_SIZE;
    static constexpr size_t MAX_NUM_OF_GRID_ATOMS = 150;
    static constexpr size_t MAX_NUM_OF_GRID_ATOMS_ = MAX_NUM_OF_GRID_ATOMS;
    static constexpr size_t FAST_SIZE = 2051;  // modified for vina1.2
    static constexpr size_t FAST_SIZE_ = FAST_SIZE;
    static constexpr size_t SMOOTH_SIZE = 2051;
    static constexpr size_t SMOOTH_SIZE_ = SMOOTH_SIZE;
    static constexpr size_t MAX_CONTAINER_SIZE_EVERY_WI = 5;
    static constexpr size_t MAX_CONTAINER_SIZE_EVERY_WI_ =
        MAX_CONTAINER_SIZE_EVERY_WI;

    static constexpr size_t MAX_THREAD = 41700000;  // modified for vina1.2
    static constexpr size_t MAX_THREAD_ = MAX_THREAD;
    static constexpr size_t MAX_LIGAND_NUM = 10250;  // modified for vina1.2
    static constexpr size_t MAX_LIGAND_NUM_ = MAX_LIGAND_NUM;
};

using DefaultConfig =
    BaseConfig<48, 1, 24, 150, 1024, 45150>;
using SizeConfig = BaseConfig<48, 1, 128, 300, 4096>;
using SmallConfig = BaseConfig<8, 1, 12, 40, 300>;
using MediumConfig = BaseConfig<16, 1, 18, 80, 600>;
using LargeConfig = BaseConfig<24, 1, 36, 120, 1024>;
using ExtraLargeConfig = BaseConfig<36, 1, 64, 160, 2048>;
using MaxConfig = BaseConfig<48, 1, 128, 300, 4096>;

template <typename Config>
struct matrix_d_ {
    float data[Config::MAX_HESSIAN_MATRIX_D_SIZE_];
    int dim;
} ;

template <typename Config>
struct affinities_cuda_t_ {
    float data[Config::GRIDS_SIZE_];
} ;

template <typename Config>
struct atom_cuda_t_ {
    int types[4];
    float coords[3];
};

template <typename Config>
struct grid_atoms_cuda_t_ {
    atom_cuda_t_<Config> atoms[Config::MAX_NUM_OF_ATOMS_];
} ;
template <typename Config>
struct m_coords_cuda_t_ {
    float coords[Config::MAX_NUM_OF_ATOMS_][3];
} ;
template <typename Config>
struct m_minus_forces_t_ {
    float coords[Config::MAX_NUM_OF_ATOMS_][3];
} ;
template <typename Config>
struct output_type_cuda_t_ {  // namely molec_struc
    float position[3];
    float orientation[4];
    float lig_torsion[Config::MAX_NUM_OF_LIG_TORSION_];
    float flex_torsion[Config::MAX_NUM_OF_FLEX_TORSION_];
    float lig_torsion_size;
    float coords[Config::MAX_NUM_OF_ATOMS_][3];
    float e;
} ;


template <typename Config>
struct change_cuda_t_{  // namely change_struc
    float position[3];
    float orientation[3];
    float lig_torsion[Config::MAX_NUM_OF_LIG_TORSION_];
    float flex_torsion[Config::MAX_NUM_OF_FLEX_TORSION_];
    float lig_torsion_size;
} ;
template <typename Config>
struct rigid_cuda_t_{  // depth-first order
    int atom_range[Config::MAX_NUM_OF_RIGID_][2];
    float origin[Config::MAX_NUM_OF_RIGID_][3];
    float orientation_m[Config::MAX_NUM_OF_RIGID_][9];  // This matrix is fixed to 3*3
    float orientation_q[Config::MAX_NUM_OF_RIGID_][4];

    float axis[Config::MAX_NUM_OF_RIGID_][3];             // 1st column is root node, all 0s
    float relative_axis[Config::MAX_NUM_OF_RIGID_][3];    // 1st column is root node, all 0s
    float relative_origin[Config::MAX_NUM_OF_RIGID_][3];  // 1st column is root node, all 0s

    int parent[Config::MAX_NUM_OF_RIGID_] ={0};  // every node has only 1 parent node
    bool children_map[Config::MAX_NUM_OF_RIGID_]
                     [Config::MAX_NUM_OF_RIGID_];  // chidren_map[i][j] = true if node i's child is node j
    bool descendant_map[Config::MAX_NUM_OF_RIGID_][Config::MAX_NUM_OF_RIGID_];
    int num_children;

} ;

template <typename Config>
struct lig_pairs_cuda_t_{
    int type_pair_index[Config::MAX_NUM_OF_LIG_PAIRS_];
    int a[Config::MAX_NUM_OF_LIG_PAIRS_];
    int b[Config::MAX_NUM_OF_LIG_PAIRS_];
    int num_pairs;
} ;
template <typename Config>
struct ligand_cuda_t_ {
    lig_pairs_cuda_t_<Config> pairs;
    rigid_cuda_t_<Config> rigid;
    int begin;
    int end;
} ;
template <typename Config>
struct random_maps_t_{
    int int_map[Config::MAX_NUM_OF_RANDOM_MAP_];
    float pi_map[Config::MAX_NUM_OF_RANDOM_MAP_];
    float sphere_map[Config::MAX_NUM_OF_RANDOM_MAP_][3];
} ;

template <typename Config>
struct m_cuda_t_{
    atom_cuda_t_<Config> atoms[Config::MAX_NUM_OF_ATOMS_];
    m_coords_cuda_t_<Config> m_coords;
    m_minus_forces_t_<Config> minus_forces;
    ligand_cuda_t_<Config> ligand;
    int m_num_movable_atoms;  // will be -1 if ligand parsing failed
} ;
template <typename Config>
struct grid_cuda_t_{
    float m_init[3];
    float m_range[3];
    float m_factor[3];
    float m_dim_fl_minus_1[3];
    float m_factor_inv[3];
    int m_i;
    int m_j;
    int m_k;
    float m_data[Config::MAX_NUM_OF_GRID_POINT_];
} ;
template <typename Config>
struct ig_cuda_t_ {
    int atu;
    float slope;
    grid_cuda_t_<Config> grids[Config::GRIDS_SIZE_];
} ;
template <typename Config>
struct p_m_data_cuda_t_{
    float fast[Config::FAST_SIZE_];
    float smooth[Config::SMOOTH_SIZE_][2];
    float factor;
} ;
template <typename Config>
struct p_cuda_t_ {
    float m_cutoff_sqr;
    int n;
    float factor;
    int m_data_size;
    p_m_data_cuda_t_<Config> *m_data;
} ;
template <typename Config>
struct p_cuda_t_cpu_ {
    float m_cutoff_sqr;
    int n;
    float factor;
    int m_data_size;
    p_m_data_cuda_t_<Config> m_data[Config::MAX_P_DATA_M_DATA_SIZE_];
} ;

template <typename Config>
struct output_container_cuda_t_{
    output_type_cuda_t_<Config> container[Config::MAX_CONTAINER_SIZE_EVERY_WI_];
    int current_size;
} ;
template <typename Config>
struct precalculate_element_cuda_t_ {
    float fast[Config::FAST_SIZE_];
    float smooth[Config::SMOOTH_SIZE_][2];  // smooth
    float factor;
} ;
template <typename Config>
struct pot_cuda_t_{
    float ptmp[Config::MAX_NUM_OF_RIGID_][3];
    float p[Config::MAX_NUM_OF_RIGID_][3];
    float o[Config::MAX_NUM_OF_RIGID_][3];
} ;

