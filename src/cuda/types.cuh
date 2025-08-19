#pragma once
#include "kernel.h"

using matrix_d = matrix_d_<DefaultConfig>;
using affinities_cuda_t = affinities_cuda_t_<DefaultConfig>;
using atom_cuda_t = atom_cuda_t_<DefaultConfig>;
using grid_atoms_cuda_t = grid_atoms_cuda_t_<DefaultConfig>;
using m_coords_cuda_t = m_coords_cuda_t_<DefaultConfig>;
using m_minus_forces_t = m_minus_forces_t_<DefaultConfig>;
using output_type_cuda_t = output_type_cuda_t_<DefaultConfig>;
using change_cuda_t = change_cuda_t_<DefaultConfig>;
using rigid_cuda_t = rigid_cuda_t_<DefaultConfig>;
using lig_pairs_cuda_t = lig_pairs_cuda_t_<DefaultConfig>;
using ligand_cuda_t = ligand_cuda_t_<DefaultConfig>;
using random_maps_t = random_maps_t_<DefaultConfig>;
using m_cuda_t = m_cuda_t_<DefaultConfig>;
using grid_cuda_t = grid_cuda_t_<DefaultConfig>;
using ig_cuda_t = ig_cuda_t_<DefaultConfig>;
using p_m_data_cuda_t = p_m_data_cuda_t_<DefaultConfig>;
using p_cuda_t = p_cuda_t_<DefaultConfig>;
using p_cuda_t_cpu = p_cuda_t_cpu_<DefaultConfig>;
using output_container_cuda_t = output_container_cuda_t_<DefaultConfig>;
using precalculate_element_cuda_t = precalculate_element_cuda_t_<DefaultConfig>;
using pot_cuda_t = pot_cuda_t_<DefaultConfig>;
