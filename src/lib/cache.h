#ifndef VINA_CACHE_H
#define VINA_CACHE_H

#include <iostream>
#include <string>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <boost/serialization/split_member.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/static_assert.hpp>
#include "igrid.h"
#include "grid.h"
#include "model.h"
#include "file.h"
#include "szv_grid.h"

struct precalculate;
struct model;

struct cache : public igrid {
public:
    cache(fl slope = 1e6) : m_slope(slope), m_grids(XS_TYPE_SIZE) {}
    cache(const grid_dims& gd, fl slope = 1e6) : m_gd(gd), m_slope(slope), m_grids(XS_TYPE_SIZE) {}
    fl eval(const model& m, fl v) const;  // needs m.coords // clean up
    fl eval_intra(model& m, fl v) const;  // needs m.coords // clean up
    fl eval_deriv(model& m, fl v) const;  // needs m.coords, sets m.minus_forces // clean up
    grid_dims get_gd() const { return m_gd; }
    vec corner1() const {
        vec corner(m_gd[0].begin, m_gd[1].begin, m_gd[2].begin);
        return corner;
    }
    vec corner2() const {
        vec corner(m_gd[0].end, m_gd[1].end, m_gd[2].end);
        return corner;
    }
    bool is_in_grid(const model& m, fl margin = 0.0001) const;
    bool is_atom_type_grid_initialized(sz t) const { return m_grids[t].initialized(); }
    bool are_atom_types_grid_initialized(szv atom_types) const;
    void read(const std::string& str);
    void write(const std::string& out_prefix, const szv& atom_types,
               const std::string& gpf_filename = "NULL", const std::string& fld_filename = "NULL",
               const std::string& receptor_filename = "NULL");
    // add for gpu
    float get_slope() const;
    std::vector<grid> get_grids() const;
    int get_atu() const;
    std::vector<grid> m_grids;
    // std::vector<grid> grids;
    void populate_no_bias(const model& m, const precalculate& p, const szv& atom_types_needed);


    // 1) 返回网格在某个方向的 sample 点总数: n_voxels + 1
    //    e.g. dim(0) = m_gd[0].n_voxels+1
    inline int dim(int d) const {
        if(d < 0 || d > 2) return 0;
        return (m_gd[d].n_voxels + 1);
    }

    // 2) 将 (x,y,z) 索引转换到实际坐标
    //    原本 grid.cpp 里就有： grid::index_to_argument(...) 
    //    这里做一个 wrapper, 以便外部 code 直接调
    inline vec index_to_argument(sz x, sz y, sz z) const {
        // 见 grid::index_to_argument for reference
        vec ret;
        ret[0] = m_gd[0].begin + ( (fl)x / (fl)(dim(0)-1) ) * (m_gd[0].end - m_gd[0].begin );
        ret[1] = m_gd[1].begin + ( (fl)y / (fl)(dim(1)-1) ) * (m_gd[1].end - m_gd[1].begin );
        ret[2] = m_gd[2].begin + ( (fl)z / (fl)(dim(2)-1) ) * (m_gd[2].end - m_gd[2].begin );
        return ret;
    }

    // 3) allocate_type: 若还没初始化 m_grids[t], 给它 init
    //    这个函数通常只在需要时调用，如 compute_ligand_maps
    inline void allocate_type(sz t, const grid_dims& gd) {
        if(!m_grids[t].initialized()) {
            m_grids[t].init(gd);
        }
    }

    // 4) 访问 m_grids[t].m_data 的引用, 以对 x,y,z 做 get/set
    //    array3d<fl>& get_array3d_ref(sz t);
    //    const array3d<fl>& get_array3d_const_ref(sz t) const;
    // 
    //    需要注意: m_grids 是 vector<grid>, 
    //              grid 里有 array3d<fl> m_data
    //    你可以直接返回 reference
    array3d<fl>& get_array3d_ref(sz t) {
        return m_grids[t].m_data;
    }
    const array3d<fl>& get_array3d_const_ref(sz t) const {
        return m_grids[t].m_data;
    }

    // 5) get_value / set_value
    //    供 c++ 端直接读写 
    inline fl get_value(sz t, sz x, sz y, sz z) const {
        return m_grids[t].m_data(x, y, z);
    }
    inline void set_value(sz t, sz x, sz y, sz z, fl val) {
        m_grids[t].m_data(x, y, z) = val;
    }
private:
    grid_dims m_gd;
    fl m_slope;  // does not get (de-)serialized
};

#endif
