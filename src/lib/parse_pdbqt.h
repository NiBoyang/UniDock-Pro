#ifndef VINA_PARSE_PDBQT_H
#define VINA_PARSE_PDBQT_H

#include <string>
#include "model.h"
#include <set>
struct rigid {
    atomv atoms;
};
model parse_receptor_pdbqt(const std::string &rigid = std::string(),
                           const std::string &flex = std::string(),
                           atom_type::t atype = atom_type::XS);  // can throw struct_parse_error
model parse_receptor_pdb(const std::string &rigid = std::string(),
                         const std::string &flex = std::string(),
                         atom_type::t atype = atom_type::XS);  // can throw struct_parse_error

model parse_ligand_pdbqt_from_file(const std::string &name, atom_type::t atype,
                                   bool keep_H = false, bool skip_recenter = false);  // can throw struct_parse_error
model parse_ligand_pdbqt_from_file_no_failure(const std::string &name, atom_type::t atype,
                                              bool keep_H = false);  // can throw struct_parse_error
model parse_ligand_from_file_no_failure(const std::string &name, atom_type::t atype,
                                        bool keep_H = false);  // can throw struct_parse_error
model parse_ligand_sdf_from_file_no_failure(const std::string &name, atom_type::t atype,
                                            bool keep_H = false);  // can throw struct_parse_error

model parse_ligand_pdbqt_from_string(const std::string &string_name,
                                     atom_type::t atype);  // can exit with code EXIT_FAILURE
model parse_ligand_pdbqt_from_string_no_failure(
    const std::string &string_name, atom_type::t atype);  // can return empty model as failure
void parse_pdbqt_rigid(const path& name, rigid& r);

#endif
