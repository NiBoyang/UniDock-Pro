#include "utils.h"

inline char separator() {

#ifdef _WIN32
    return '\\';
#else
    return '/';
#endif
}

path make_path(const std::string& str) {
    boost::filesystem::path p(str);
    return p;
}

void doing(const std::string& str, int verbosity, int level) {
    if (verbosity > level) {
        std::cout << str << std::string(" ... ") << std::flush;
    }
}

void done(int verbosity, int level) {
    if (verbosity > level) {
        std::cout << "done.\n" << std::flush;
    }
}

std::string default_output(const std::string& input_name) {
    std::string tmp = input_name;
    if (tmp.size() >= 6 && tmp.substr(tmp.size() - 6, 6) == ".pdbqt") {
        tmp.resize(tmp.size() - 6);  // FIXME?
        return tmp + "_out.pdbqt";
    } else if (tmp.size() >= 4 && tmp.substr(tmp.size() - 4, 4) == ".sdf") {
        tmp.resize(tmp.size() - 4);  // FIXME?
        return tmp + "_out.sdf";
    }
    return tmp + "_out.pdbqt";
}

std::string default_score_output(const std::string& input_name) {
    std::string tmp = input_name;
    if (tmp.size() >= 6 && tmp.substr(tmp.size() - 6, 6) == ".pdbqt") {
        tmp.resize(tmp.size() - 6);  // FIXME?
        return tmp + "_score.txt";
    } else if (tmp.size() >= 4 && tmp.substr(tmp.size() - 4, 4) == ".sdf") {
        tmp.resize(tmp.size() - 4);  // FIXME?
        return tmp + "_score.txt";
    }
    return tmp + "_score.txt";
}

std::string default_output(const std::string& input_name, const std::string& directory_pathname) {
    return directory_pathname + separator() + default_output(input_name);
}

bool is_directory(const std::string& directory_pathname) {
    // Source:
    // https://stackoverflow.com/questions/18100097/portable-way-to-check-if-directory-exists-windows-linux-c
    struct stat info;

    if (stat(directory_pathname.c_str(), &info) != 0)
        return false;
    else if (info.st_mode & S_IFDIR)  // S_ISDIR() doesn't exist on my windows
        return true;
    else
        return false;
}

std::string get_filename(const std::string& s) {
    size_t i = s.rfind(separator(), s.length());

    if (i != std::string::npos) {
        return (s.substr(i + 1, s.length() - i));
    }

    return (s);
}


std::string get_file_contents(const std::string& filename) {
    ifile in(make_path(filename));

    if (in) {
        std::string contents;
        in.seekg(0, std::ios::end);
        contents.resize(in.tellg());
        in.seekg(0, std::ios::beg);
        in.read(&contents[0], contents.size());
        in.close();
        return (contents);
    }
    throw file_error(filename, true);
}
