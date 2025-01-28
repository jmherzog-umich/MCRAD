#include <random>
#include <string>
#include <algorithm>
#include <cctype>
#include <iostream>

#define FP_FAST_FMA

#ifndef _UTILITY_H_
#define _UTILITY_H_

namespace util {

    //Random values
    namespace {
        std::default_random_engine _GEN;
        std::uniform_real_distribution<double> _DIST;
        std::uniform_real_distribution<float> _DISTF;
        std::exponential_distribution<double> _LOG;
    }

    //Helper functions
    bool isvec(const std::string& s);
    std::string makefilename(const std::string& name, const std::string& ext, int id);
    void writeheader(std::ostream& oout, std::string label);
    std::vector<std::string> evalrange(const std::string& t, int id, int max);
    std::string readfile(const std::string& name);
    std::vector<std::string> splitcomma(const std::string& s);

    //Methods
    void rand_init();
    double roll();
    float rollf();
    double logroll();
    float erfinvf(float a);
    std::string bitstring(uint_least16_t v);
}

#endif
