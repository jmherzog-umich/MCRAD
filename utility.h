#include <cmath>
#include <random>
#include <ctime>
#include <string>
#include <cstdint>
#include <algorithm>
#include <cctype>
#include <fstream>
#include <iostream>

#define FP_FAST_FMA

#define CONST_PI            3.1415926535897932
#define CONST_C             299.792458
#define CONST_HBAR          1.054571817679489e-22
#define CONST_HK            47.992430733662212
#define CONST_APERY         1.2020569031595942
#define CONST_WIEN          1.5936242600400401
#define CONST_EPS           1e-10
#define CONST_PLANCKMAX     0.6476102378919149

#ifndef _UTILITY_H_
#define _UTILITY_H_

using namespace std;

//Helper functions
bool isvec(const string& s) {
    for (char x : s)
        if (x == ',' or x == ';')
            return true;
    return false;
}

string makefilename(const string& name, const string& ext, int id) {
    string base = name.substr(0, name.find('.'));
    if (id >= 0)
        return base + "." + to_string(id) + "." + ext;
    else
        return base + "." + ext;
}

void writeheader(ostream& oout, string label) {
    oout << "==================================================================" << endl;
    oout << label << endl;
    oout << "==================================================================" << endl;
}

vector<string> evalrange(const string& t, int id, int max) {
    //First, exit quick if this value isn't a range
    vector<string> out;
    unsigned long int i0 = t.find(":");
    if (i0 == string::npos) {
        out.push_back(t);
        return out;
    }
    double x0, xf, dx;
    unsigned long int i2 = t.find(":", i0+1);
    
    //Now get the values
    bool log = false;
    x0 = stod(t.substr(0, i0));
    if (i2 == string::npos) {
        xf = stod(t.substr(i0+1));
        dx = 1;
    } else {
        
        //Check for logscale
        if (t.at(i2+1) == '*')
            log = true;
            
        //Now get values
        xf = stod(t.substr(i0+1, i2-i0-1));
        if (log)
            dx = stod(t.substr(i2+2));
        else
            dx = stod(t.substr(i2+1));
    }
    
    //Now loop through the range
    int i = 0;
    for (double x = x0; x <= xf; (log) ? x*=dx : x+=dx) {
        if (i % max == id)
            out.push_back(to_string(x));
        i++;
    }
    return out;
}

string readfile(const string& name) {

    //Allocate string and load file
    string out;
    ifstream file = ifstream(name);
        
    //If file is bad, return empty string
    if (not file)
        return out;
    
    //Allocate a buffer string to read line by line
    string buf;
    while (not file.eof()) {
        getline(file, buf);
        out.append(buf);
        out.append("\n");
    }
    return out;
}

vector<string> splitcomma(const string& s) {
    vector<string> out;
    string tmp = "";
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

//Random values
namespace {
    default_random_engine _GEN;
    uniform_real_distribution<double> _DIST;
    uniform_real_distribution<float> _DISTF;
    exponential_distribution<double> _LOG;
}

//Methods
void rand_init() {
    _GEN.seed(time(0));
    _DIST = uniform_real_distribution<double>(0.0,1.0);
    _DISTF = uniform_real_distribution<float>(0.0,1.0);
    _LOG = exponential_distribution<double>(1.0);
}

double roll() {
    double eps = _DIST(_GEN);
    while (eps == 1.0)
        eps = _DIST(_GEN);
    return eps;
}

float rollf() {
    float eps = _DISTF(_GEN);
    while (eps == 1.0)
        eps = _DISTF(_GEN);
    return eps;
}

double logroll() {
    return _LOG(_GEN);
}

float erfinvf (float a)
{
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

string bitstring(uint_least16_t v) {
    string out;
    for (int i = 0; i < 16; i ++)
        if ( v & (1 << i) )
            out.push_back('1');
        else
            out.push_back('0');
    return out;
}

#endif
