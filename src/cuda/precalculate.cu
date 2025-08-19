#include "kernel.h"
#include "types.cuh"
#include "cuda_memory.hpp"
#include "math.h"
#include "model.h"
#include <vector>

#include "precalculate.h"
#include "precalculate_gpu.cuh"

__constant__ scoring_function_cuda_t scoring_cuda_gpu_const;
__constant__ fl common_rs_gpu_const[DefaultConfig::FAST_SIZE];


#define PRECALC_BLOCK_SIZE 256

__global__ void precalculate_gpu(triangular_matrix_cuda_t *m_data_gpu_list,
                                 sz *atom_xs_gpu, sz *atom_ad_gpu,
                                 fl *atom_charge_gpu, int *atom_num_gpu,
                                 fl factor, fl max_fl, int ligand_count,
                                 int max_atom_num) {
    int ligand_idx = blockIdx.x;
    if (ligand_idx >= ligand_count) return;

    sz *xs_ptr = atom_xs_gpu + ligand_idx * max_atom_num;
    sz *ad_ptr = atom_ad_gpu + ligand_idx * max_atom_num;
    fl *charge_ptr = atom_charge_gpu + ligand_idx * max_atom_num;
    (void)ad_ptr;
    (void)charge_ptr;

    int atom_num = atom_num_gpu[ligand_idx];
    precalculate_element_cuda_t *p_data_gpu = m_data_gpu_list[ligand_idx].p_data;

    int pair_count = atom_num * (atom_num + 1) / 2;

    for (int pair_idx = threadIdx.x; pair_idx < pair_count; pair_idx += blockDim.x) {
        int j = int((sqrtf(8.0f * pair_idx + 1.0f) - 1.0f) * 0.5f);
        int i = pair_idx - j * (j + 1) / 2;
        int offset = pair_idx;
        int n = DefaultConfig::SMOOTH_SIZE;
        p_data_gpu[offset].factor = 32.0f;
        switch (scoring_cuda_gpu_const.m_sf_choice) {
            case SF_VINA: {
                for (int k = 0; k < n; ++k) {
                    fl sum = 0.0f;
                    sum += scoring_cuda_gpu_const.m_weights[0] *
                           vina_gaussian_cuda_eval(xs_ptr[i], xs_ptr[j], common_rs_gpu_const[k],
                                                   scoring_cuda_gpu_const.vina_gaussian_cutoff_1,
                                                   scoring_cuda_gpu_const.vina_gaussian_offset_1,
                                                   scoring_cuda_gpu_const.vina_gaussian_width_1);
                    sum += scoring_cuda_gpu_const.m_weights[1] *
                           vina_gaussian_cuda_eval(xs_ptr[i], xs_ptr[j], common_rs_gpu_const[k],
                                                   scoring_cuda_gpu_const.vina_gaussian_cutoff_2,
                                                   scoring_cuda_gpu_const.vina_gaussian_offset_2,
                                                   scoring_cuda_gpu_const.vina_gaussian_width_2);
                    sum += scoring_cuda_gpu_const.m_weights[2] *
                           vina_repulsion_cuda_eval(xs_ptr[i], xs_ptr[j], common_rs_gpu_const[k],
                                                    scoring_cuda_gpu_const.vina_repulsion_cutoff,
                                                    scoring_cuda_gpu_const.vina_repulsion_offset);
                    sum += scoring_cuda_gpu_const.m_weights[3] *
                           vina_hydrophobic_cuda_eval(xs_ptr[i], xs_ptr[j], common_rs_gpu_const[k],
                                                     scoring_cuda_gpu_const.vina_hydrophobic_good,
                                                     scoring_cuda_gpu_const.vina_hydrophobic_bad,
                                                     scoring_cuda_gpu_const.vina_hydrophobic_cutoff);
                    sum += scoring_cuda_gpu_const.m_weights[4] *
                           vina_non_dir_h_bond_cuda_eval(xs_ptr[i], xs_ptr[j], common_rs_gpu_const[k],
                                                         scoring_cuda_gpu_const.vina_non_dir_h_bond_good,
                                                         scoring_cuda_gpu_const.vina_non_dir_h_bond_bad,
                                                         scoring_cuda_gpu_const.vina_non_dir_h_bond_cutoff);
                    sum += scoring_cuda_gpu_const.m_weights[5] *
                           linearattraction_eval(xs_ptr[i], xs_ptr[j], common_rs_gpu_const[k],
                                                scoring_cuda_gpu_const.linearattraction_cutoff);
                    p_data_gpu[offset].smooth[k][0] = sum;
                }
                break;
            }
            case SF_VINARDO: {
                for (int k = 0; k < n; ++k) {
                    fl sum = 0.0f;
                    sum += scoring_cuda_gpu_const.m_weights[0] *
                           vinardo_gaussian_eval(xs_ptr[i], xs_ptr[j], common_rs_gpu_const[k],
                                                scoring_cuda_gpu_const.vinardo_gaussian_offset,
                                                scoring_cuda_gpu_const.vinardo_gaussian_width,
                                                scoring_cuda_gpu_const.vinardo_gaussian_cutoff);
                    sum += scoring_cuda_gpu_const.m_weights[1] *
                           vinardo_repulsion_eval(xs_ptr[i], xs_ptr[j], common_rs_gpu_const[k],
                                                  scoring_cuda_gpu_const.vinardo_repulsion_cutoff,
                                                  scoring_cuda_gpu_const.vinardo_repulsion_offset);
                    sum += scoring_cuda_gpu_const.m_weights[2] *
                           vinardo_hydrophobic_eval(xs_ptr[i], xs_ptr[j], common_rs_gpu_const[k],
                                                   scoring_cuda_gpu_const.vinardo_hydrophobic_good,
                                                   scoring_cuda_gpu_const.vinardo_hydrophobic_bad,
                                                   scoring_cuda_gpu_const.vinardo_hydrophobic_cutoff);
                    sum += scoring_cuda_gpu_const.m_weights[3] *
                           vinardo_non_dir_h_bond_eval(xs_ptr[i], xs_ptr[j], common_rs_gpu_const[k],
                                                       scoring_cuda_gpu_const.vinardo_non_dir_h_bond_good,
                                                       scoring_cuda_gpu_const.vinardo_non_dir_h_bond_bad,
                                                       scoring_cuda_gpu_const.vinardo_non_dir_h_bond_cutoff);
                    sum += scoring_cuda_gpu_const.m_weights[4] *
                           linearattraction_eval(xs_ptr[i], xs_ptr[j], common_rs_gpu_const[k],
                                                 scoring_cuda_gpu_const.linearattraction_cutoff);
                    p_data_gpu[offset].smooth[k][0] = sum;
                }
                break;
            }
            default:
                break;
        }
        for (int k = 0; k < n; ++k) {
            fl dor;
            if (k == 0 || k == n - 1) {
                dor = 0.0f;
            } else {
                fl delta = common_rs_gpu_const[k + 1] - common_rs_gpu_const[k - 1];
                fl r = common_rs_gpu_const[k];
                dor = (p_data_gpu[offset].smooth[k + 1][0] -
                       p_data_gpu[offset].smooth[k - 1][0]) /
                      (delta * r);
            }
            p_data_gpu[offset].smooth[k][1] = dor;
            fl f1 = p_data_gpu[offset].smooth[k][0];
            fl f2 = (k + 1 >= n) ? 0.0f : p_data_gpu[offset].smooth[k + 1][0];
            p_data_gpu[offset].fast[k] = (f2 + f1) / 2.0f;
        }
    }
}

void precalculate_parallel(triangular_matrix_cuda_t *m_data_list_cpu,
                           std::vector<precalculate_byatom> &m_precalculated_byatom_gpu,
                           const ScoringFunction &m_scoring_function,
                           std::vector<model> &m_model_gpu, const flv &common_rs, int thread) {
    // TODO: copy and transfer data to gpu array
    int max_atom_num = 0;
    for (int i = 0; i < thread; ++i) {
        if ((int)m_model_gpu[i].num_atoms() > max_atom_num)
            max_atom_num = (int)m_model_gpu[i].num_atoms();
    }
    DEBUG_PRINTF("max_atom_num = %d\n", max_atom_num);

    // copy atomv from m_model_gpu array and put into atom array

    // using cudaMallocHost may be better
    assert(DefaultConfig::MAX_LIGAND_NUM >= thread);
    sz atom_xs[thread * max_atom_num];  // maybe size_t is too large
    sz atom_ad[thread * max_atom_num];
    fl atom_charge[thread * max_atom_num];
    int atom_num[thread];
    int precalculate_matrix_size[thread];

    CudaMemory<sz> atom_xs_gpu(thread * max_atom_num);
    CudaMemory<sz> atom_ad_gpu(thread * max_atom_num);
    CudaMemory<fl> atom_charge_gpu(thread * max_atom_num);
    CudaMemory<int> atom_num_gpu(thread);

    for (int i = 0; i < thread; ++i) {
        atomv atoms = m_model_gpu[i].get_atoms();
        atom_num[i] = atoms.size();
        precalculate_matrix_size[i] = atom_num[i] * (atom_num[i] + 1) / 2;
        for (int j = 0; j < atoms.size(); ++j) {
            atom_xs[i * max_atom_num + j] = atoms[j].xs;
            DEBUG_PRINTF("atom[%d] on CPU: xs=%lu %lu\n", j, atoms[j].xs,
                         atom_xs[i * max_atom_num + j]);
            atom_ad[i * max_atom_num + j] = atoms[j].ad;
            atom_charge[i * max_atom_num + j] = atoms[j].charge;
        }
    }

    checkCUDA(cudaMemcpy(atom_xs_gpu, atom_xs, thread * max_atom_num * sizeof(sz),
                         cudaMemcpyHostToDevice));
    // // debug
    // sz atom_xs_check[max_atom_num];
    // checkCUDA(cudaMemcpy(atom_xs_check, atom_xs_gpu, max_atom_num * sizeof(sz),
    // cudaMemcpyDeviceToHost)); for (int i = 0;i < max_atom_num;++i){
    //     DEBUG_PRINTF("atom[%d] on gpu check: xs=%lu\n", i, atom_xs_check[i]);
    // }
    checkCUDA(cudaMemcpy(atom_ad_gpu, atom_ad, thread * max_atom_num * sizeof(sz),
                         cudaMemcpyHostToDevice));
    checkCUDA(cudaMemcpy(atom_charge_gpu, atom_charge, thread * max_atom_num * sizeof(fl),
                         cudaMemcpyHostToDevice));
    checkCUDA(cudaMemcpy(atom_num_gpu, atom_num, thread * sizeof(int), cudaMemcpyHostToDevice));

    // copy scoring function parameters
    scoring_function_cuda_t scoring_cuda;
    scoring_cuda.m_num_potentials = m_scoring_function.m_potentials.size();
    for (int w = 0; w < scoring_cuda.m_num_potentials; ++w) {
        scoring_cuda.m_weights[w] = m_scoring_function.m_weights[w];
    }
    scoring_cuda.m_sf_choice = m_scoring_function.m_sf_choice;
    switch (scoring_cuda.m_sf_choice) {
        case SF_VINA:  // vina
        {
            scoring_cuda.vina_gaussian_offset_1 = 0;
            scoring_cuda.vina_gaussian_width_1 = 0.5;
            scoring_cuda.vina_gaussian_cutoff_1 = 8;
            scoring_cuda.vina_gaussian_offset_2 = 3;
            scoring_cuda.vina_gaussian_width_2 = 2.0;
            scoring_cuda.vina_gaussian_cutoff_2 = 8;
            scoring_cuda.vina_repulsion_offset = 0.0;
            scoring_cuda.vina_repulsion_cutoff = 8.0;
            scoring_cuda.vina_hydrophobic_good = 0.5;
            scoring_cuda.vina_hydrophobic_bad = 1.5;
            scoring_cuda.vina_hydrophobic_cutoff = 8;
            scoring_cuda.vina_non_dir_h_bond_good = -0.7;
            scoring_cuda.vina_non_dir_h_bond_bad = 0;
            scoring_cuda.vina_non_dir_h_bond_cutoff = 8.0;
            scoring_cuda.linearattraction_cutoff = 20.0;
            break;
        }
        case SF_VINARDO:  // vinardo
        {
            scoring_cuda.vinardo_gaussian_offset = 0;
            scoring_cuda.vinardo_gaussian_width = 0.8;
            scoring_cuda.vinardo_gaussian_cutoff = 8.0;
            scoring_cuda.vinardo_repulsion_offset = 0;
            scoring_cuda.vinardo_repulsion_cutoff = 8.0;
            scoring_cuda.vinardo_hydrophobic_good = 0;
            scoring_cuda.vinardo_hydrophobic_bad = 2.5;
            scoring_cuda.vinardo_hydrophobic_cutoff = 8.0;
            scoring_cuda.vinardo_non_dir_h_bond_good = -0.6;
            scoring_cuda.vinardo_non_dir_h_bond_bad = 0;
            scoring_cuda.vinardo_non_dir_h_bond_cutoff = 8.0;
            scoring_cuda.linearattraction_cutoff = 20.0;
            break;
        }
    }

    // transfer common_rs to gpu
    checkCUDA(cudaMemcpyToSymbol(common_rs_gpu_const, common_rs.data(), DefaultConfig::FAST_SIZE * sizeof(fl)));

    // malloc output buffer for m_data, array of precalculate_element
    CudaMemory<triangular_matrix_cuda_t> m_data_gpu_list(thread);
    triangular_matrix_cuda_t m_data_cpu_list[thread];
    triangular_matrix_cuda_t m_data_gpu;
    precalculate_element_cuda_t *precalculate_element_list_ptr_gpu;
    for (int i = 0; i < thread; ++i) {
        checkCUDA(cudaMalloc(&precalculate_element_list_ptr_gpu,
                             sizeof(precalculate_element_cuda_t) * precalculate_matrix_size[i]));
        m_data_gpu.p_data = precalculate_element_list_ptr_gpu;
        m_data_cpu_list[i].p_data = precalculate_element_list_ptr_gpu;
        checkCUDA(cudaMemcpy(m_data_gpu_list + i, &m_data_gpu, sizeof(triangular_matrix_cuda_t),
                             cudaMemcpyHostToDevice));
    }

    // TODO: launch kernel
    DEBUG_PRINTF("launch kernel precalculate_gpu, thread=%d\n", thread);
    checkCUDA(cudaMemcpyToSymbol(scoring_cuda_gpu_const, &scoring_cuda, sizeof(scoring_function_cuda_t)));
    precalculate_gpu<<<thread, PRECALC_BLOCK_SIZE>>>(m_data_gpu_list, atom_xs_gpu,
                                                    atom_ad_gpu, atom_charge_gpu,
                                                    atom_num_gpu, 32, max_fl,
                                                    thread, max_atom_num);

    checkCUDA(cudaDeviceSynchronize());

    DEBUG_PRINTF("kernel exited\n");

    memcpy(m_data_list_cpu, m_data_cpu_list, sizeof(m_data_cpu_list));

    // device memory is released automatically when going out of scope
}