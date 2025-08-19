#ifndef VINA_PARSE_ERROR_H
#define VINA_PARSE_ERROR_H

#include "common.h"

class struct_parse_error : public std::exception {
public:
    explicit struct_parse_error(const std::string& message)
        : m_message("\n\nStructure parsing error: " + message + "\n") {}
    explicit struct_parse_error(const std::string& message, const std::string& pdbqt_line)
        : m_message("\n\nStructure parsing error: " + message + "\n > " + pdbqt_line + "\n") {}

    virtual const char* what() const throw() { return m_message.c_str(); }

private:
    const std::string m_message;
};

#endif
