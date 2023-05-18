#ifndef ROOT_TEST_PARSER_HPP
#define ROOT_TEST_PARSER_HPP

#include <algorithm>
#include <fmt/core.h>
#include <fmt/ranges.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>

class Parser {
private:
    std::unordered_map<std::string, std::string> map;

public:
    explicit Parser(const std::string &file) {
        std::ifstream cFile(file);
        if (cFile.is_open()) {
            std::string line;
            while (getline(cFile, line)) {
                line.erase(std::remove_if(line.begin(), line.end(), isspace), line.end());
                if (line.empty() or (line[0] == '#')) continue;
                std::string line_wo_comment = line.substr(0, line.find('#'));
                unsigned long delimiter_eq  = line_wo_comment.find('=');
                std::string name            = line_wo_comment.substr(0, delimiter_eq);
                std::string value           = line_wo_comment.substr(delimiter_eq + 1);
                map.insert({name, value});
            }
        } else {
            std::cerr << "Couldn't open config file for reading." << std::endl;
        }
    }

    template<typename T>
    auto get(const std::string &s) -> T {
        T v;
        std::stringstream(map.at(s)) >> v;
        return v;
    }
};

struct ParserOutput {
    std::string output;
    double EBene, PHmom;
    int target, process, model, Nb_event, NB_cells, samples, seed;
    double y_min, y_max, Q2_min, Q2_max, W2_min, theta_e_max, t_min, t_range, e_helicity;
};

#endif//ROOT_TEST_PARSER_HPP
