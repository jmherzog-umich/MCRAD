#include <cmath>
#include <random>
#include <ctime>
#include <string>
#include <cstdint>
#include <algorithm>
#include <cctype>
#include <fstream>
#include <iostream>

#include "utility.h"

#define FP_FAST_FMA

#define CONST_PI            3.1415926535897932
#define CONST_C             299.792458
#define CONST_HBAR          1.054571817679489e-22
#define CONST_HK            47.992430733662212
#define CONST_APERY         1.2020569031595942
#define CONST_WIEN          1.5936242600400401
#define CONST_EPS           1e-10
#define CONST_PLANCKMAX     0.6476102378919149

//Helper functions
bool util::isvec(const std::string& s) {
    for (char x : s)
        if (x == ',' or x == ';')
            return true;
    return false;
}

std::string util::makefilename(const std::string& name, const std::string& ext, int id) {
    std::string base = name.substr(0, name.find('.'));
    if (id >= 0)
        return base + "." + std::to_string(id) + "." + ext;
    else
        return base + "." + ext;
}

void util::writeheader(std::ostream& oout, std::string label) {
    oout << "==================================================================\n";
    oout << label << "\n";
    oout << "==================================================================" << std::endl;
}

std::vector<std::string> util::evalrange(const std::string& t, int id, int max) {
    //First, exit quick if this value isn't a range
    std::vector<std::string> out;
    unsigned long int i0 = t.find(":");
    if (i0 == std::string::npos) {
        out.push_back(t);
        return out;
    }
    double x0, xf, dx;
    unsigned long int i2 = t.find(":", i0+1);
    
    //Now get the values
    bool log = false;
    x0 = std::stod(t.substr(0, i0));
    if (i2 == std::string::npos) {
        xf = std::stod(t.substr(i0+1));
        dx = 1;
    } else {
        
        //Check for logscale
        if (t.at(i2+1) == '*')
            log = true;
            
        //Now get values
        xf = std::stod(t.substr(i0+1, i2-i0-1));
        if (log)
            dx = std::stod(t.substr(i2+2));
        else
            dx = std::stod(t.substr(i2+1));
    }
    
    //Now loop through the range
    int i = 0;
    for (double x = x0; x <= xf; (log) ? x*=dx : x+=dx) {
        if (i % max == id)
            out.push_back(std::to_string(x));
        i++;
    }
    return out;
}

std::string util::readfile(const std::string& name) {

    //Allocate string and load file
    std::string out;
    std::ifstream file = std::ifstream(name);
        
    //If file is bad, return empty string
    if (not file)
        return out;
    
    //Allocate a buffer string to read line by line
    std::string buf;
    while (not file.eof()) {
        std::getline(file, buf);
        out.append(buf);
        out.append("\n");
    }
    return out;
}

std::vector<std::string> util::splitcomma(const std::string& s) {
    std::vector<std::string> out;
    std::string tmp = "";
    for (char x : s+',') {
        if (x == ',' or x == ';') {
            if (tmp.length() > 0) {
                out.push_back(tmp);
                tmp = "";
            }
        } else
            tmp.push_back(x);
    }
    return out;
}

void util::rand_init() {
    util::_GEN.seed(time(0));
    util::_DIST = std::uniform_real_distribution<double>(0.0,1.0);
    util::_DISTF = std::uniform_real_distribution<float>(0.0,1.0);
    util::_LOG = std::exponential_distribution<double>(1.0);
}

double util::roll() {
    double eps = util::_DIST(util::_GEN);
    while (eps == 1.0)
        eps = util::_DIST(util::_GEN);
    return eps;
}

float util::rollf() {
    float eps = util::_DISTF(util::_GEN);
    while (eps == 1.0)
        eps = util::_DISTF(util::_GEN);
    return eps;
}

double util::logroll() {
    return util::_LOG(util::_GEN);
}

//Code adapted from:
// URL: https://stackoverflow.com/questions/27229371/inverse-error-function-in-c
// Author: njuffa
// Originally licensed under: CC BY-SA 4.0
// Modifications: Adapted from original code to use built-in log function
float util::erfinvf(float a) {
    float p, r, t;
    t = fmaf (a, 0.0f - a, 1.0f);
    t = logf (t);
    if (fabsf(t) > 6.125f) { 
        p =              3.03697567e-10f; //  0x1.4deb44p-32 
        p = fmaf (p, t,  2.93243101e-8f); //  0x1.f7c9aep-26 
        p = fmaf (p, t,  1.22150334e-6f); //  0x1.47e512p-20 
        p = fmaf (p, t,  2.84108955e-5f); //  0x1.dca7dep-16 
        p = fmaf (p, t,  3.93552968e-4f); //  0x1.9cab92p-12 
        p = fmaf (p, t,  3.02698812e-3f); //  0x1.8cc0dep-9 
        p = fmaf (p, t,  4.83185798e-3f); //  0x1.3ca920p-8 
        p = fmaf (p, t, -2.64646143e-1f); // -0x1.0eff66p-2 
        p = fmaf (p, t,  8.40016484e-1f); //  0x1.ae16a4p-1 
    } else { 
        p =              5.43877832e-9f;  //  0x1.75c000p-28 
        p = fmaf (p, t,  1.43285448e-7f); //  0x1.33b402p-23 
        p = fmaf (p, t,  1.22774793e-6f); //  0x1.499232p-20 
        p = fmaf (p, t,  1.12963626e-7f); //  0x1.e52cd2p-24 
        p = fmaf (p, t, -5.61530760e-5f); // -0x1.d70bd0p-15 
        p = fmaf (p, t, -1.47697632e-4f); // -0x1.35be90p-13 
        p = fmaf (p, t,  2.31468678e-3f); //  0x1.2f6400p-9 
        p = fmaf (p, t,  1.15392581e-2f); //  0x1.7a1e50p-7 
        p = fmaf (p, t, -2.32015476e-1f); // -0x1.db2aeep-3 
        p = fmaf (p, t,  8.86226892e-1f); //  0x1.c5bf88p-1 
    }
    r = a * p;
    return r;
}

std::string util::bitstring(uint_least16_t v) {
    std::string out;
    for (int i = 0; i < 16; i ++)
        if ( v & (1 << i) )
            out.push_back('1');
        else
            out.push_back('0');
    return out;
}


