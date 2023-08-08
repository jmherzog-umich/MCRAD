#include <cmath>
#include <vector>
#include <fstream>
#include <string>

#include "vec.h"

#ifndef _RAYPATH_H_
#define _RAYPATH_H_

using namespace std;

///
//Raypath class
///
struct RayPath {
    double f;       //frequency
    unsigned long int i;
    vector<double> x, y, z, t, I;
    
    static string ofbasename;
    
    void collide(const vec& xi, double t, double I);
    void print(ofstream& FILE) const;
    
    RayPath();
    RayPath(double v);
    RayPath(double v, double I0);
    RayPath(double v, double I0, const vec& x0, double t0);

    private:
        static unsigned long int i0;    //Number of actual instances
};

void RayPath::collide(const vec& xi, double ti, double Ii) {
    x.push_back(xi.X);
    y.push_back(xi.Y);
    z.push_back(xi.Z);
    t.push_back(ti);
    I.push_back(Ii);
}

void RayPath::print(ofstream& FILE) const {
    char dlm = ',';
    FILE << "#Raypath " << i << ": frequency = " << f << " THz" << endl;
    FILE << "#x[um], y[um], z[um], t[ps], I[photons]" << endl;
    for (unsigned int j = 0; j < x.size(); j ++)
        FILE << x.at(j) << dlm << y.at(j) << dlm << z.at(j) << dlm << t.at(j) << dlm << I.at(j) << endl;
    FILE << endl << endl;
}
    
RayPath::RayPath() {
    i = i0; i0 += 1; f = 0;
}

RayPath::RayPath(double v) {
    i = i0; i0 += 1; f = v;
}

RayPath::RayPath(double v, double I0) {
    i = i0; i0 += 1; f = v;
    I.push_back(I0);
    x.push_back(0);
    y.push_back(0);
    z.push_back(0);
    t.push_back(0);
}

RayPath::RayPath(double v, double I0, const vec& x0, double t0) {
    i = i0; i0 += 1;
    I.push_back(I0);
    x.push_back(x0.X);
    y.push_back(x0.Y);
    z.push_back(x0.Z);
    t.push_back(t0);
}

string RayPath::ofbasename = "raypath";
unsigned long int RayPath::i0 = 0;

#endif
