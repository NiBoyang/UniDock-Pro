#ifndef VINA_UTILS_H
#define VINA_UTILS_H

#include <string>
#include <sys/types.h>
#include <sys/stat.h>
#include "file.h"

inline char separator();
path make_path(const std::string& str);
void doing(const std::string& str, int verbosity, int level = 0);
void done(int verbosity, int level = 0);
std::string default_output(const std::string& input_name);
std::string default_output(const std::string& input_name, const std::string& directory_pathname);
std::string default_score_output(const std::string& input_name);
bool is_directory(const std::string& directory_pathname);
std::string get_filename(const std::string& s);
std::string get_file_contents(const std::string& filename);

#endif