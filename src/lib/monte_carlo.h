#ifndef VINA_MONTE_CARLO_H
#define VINA_MONTE_CARLO_H

#include "model.h"
#include "kernel.h"
#include "grid.h"
#include "precalculate.h"

struct monte_carlo {
    unsigned max_evals;
    unsigned global_steps;
    fl temperature;
    vec hunt_cap;
    fl min_rmsd;
    sz num_saved_mins;
    fl mutation_amplitude;
    unsigned local_steps;
    unsigned threads_per_ligand;
    unsigned num_of_ligands;
    bool local_only;
    unsigned thread = 2048;  // for CUDA parallel option, num_of_ligands * threads_per_ligand
    
    // T = 600K, R = 2cal/(K*mol) -> temperature = RT = 1.2;  global_steps = 50*lig_atoms = 2500
    monte_carlo()
        : max_evals(0),
          global_steps(2500),
          threads_per_ligand(2048),
          temperature(1.2),
          hunt_cap(10, 1.5, 10),
          min_rmsd(0.5),
          num_saved_mins(50),
          mutation_amplitude(2) {}

    output_type operator()(model& m, const precalculate_byatom& p, const igrid& ig,
                           const vec& corner1, const vec& corner2, rng& generator) const;
    // out is sorted
    void operator()(model& m, output_container& out, const precalculate_byatom& p, const igrid& ig,
                    const vec& corner1, const vec& corner2, rng& generator) const;
    void operator()(std::vector<model>& m, std::vector<output_container>& out,
                    std::vector<precalculate_byatom>& p, triangular_matrix_cuda_t* m_data_list_gpu,
                    const igrid& ig, const vec& corner1, const vec& corner2, rng& generator,
                    int verbosity, unsigned long long seed) const;
    void mc_stream(std::vector<model>& m, std::vector<output_container>& out,
                   std::vector<precalculate_byatom>& p, triangular_matrix_cuda_t* m_data_list_gpu,
                   const igrid& ig, const vec& corner1, const vec& corner2, rng& generator,
                   int verbosity, unsigned long long seed) const;

}; 
struct monte_carlo_template {
    unsigned max_evals;
    unsigned global_steps;
    fl temperature;
    vec hunt_cap;
    fl min_rmsd;
    sz num_saved_mins;
    fl mutation_amplitude;
    unsigned local_steps;
    unsigned threads_per_ligand;
    unsigned num_of_ligands;
    bool local_only;
    unsigned thread = 2048;  // for CUDA parallel option, num_of_ligands * threads_per_ligand
    
    // T = 600K, R = 2cal/(K*mol) -> temperature = RT = 1.2;  global_steps = 50*lig_atoms = 2500
    monte_carlo_template()
        : max_evals(0),
          global_steps(2500),
          threads_per_ligand(2048),
          temperature(1.2),
          hunt_cap(10, 1.5, 10),
          min_rmsd(0.5),
          num_saved_mins(50),
          mutation_amplitude(2) {}

    output_type operator()(model& m, const precalculate_byatom& p, const igrid& ig,
                           const vec& corner1, const vec& corner2, rng& generator) const;
    // out is sorted
    void operator()(model& m, output_container& out, const precalculate_byatom& p, const igrid& ig,
                    const vec& corner1, const vec& corner2, rng& generator) const;
    void operator()(std::vector<model>& m, std::vector<output_container>& out,
                    std::vector<precalculate_byatom>& p, triangular_matrix_cuda_t* m_data_list_gpu,
                    const igrid& ig, const vec& corner1, const vec& corner2, rng& generator,
                    int verbosity, unsigned long long seed) const;
    template <typename Config>
    void do_docking();
    template <typename Config>
    void run_search(std::vector<model>& m, std::vector<output_container>& out,
                    std::vector<precalculate_byatom>& p, triangular_matrix_cuda_t* m_data_list_gpu,
                    const igrid& ig, const vec& corner1, const vec& corner2, rng& generator,
                    int verbosity, unsigned long long seed) const;
    template <typename Config>
    void do_docking(std::vector<model>& m, std::vector<output_container>& out,
                    std::vector<precalculate_byatom>& p, triangular_matrix_cuda_t* m_data_list_gpu,
                    const igrid& ig, const vec& corner1, const vec& corner2, rng& generator,
                    int verbosity, unsigned long long seed) const;
};


#endif
