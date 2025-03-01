#include <cmath>
#include <vector>
#include <fstream>
#include <string>

#include "vec.h"
#include "raypath.h"

std::string RayPath::ofbasename = "raypath";
unsigned long int RayPath::i0 = 0;

void RayPath::collide(const vec& xi, double ti, double Ii) {
    x.push_back(xi.X);
    y.push_back(xi.Y);
    z.push_back(xi.Z);
    t.push_back(ti);
    I.push_back(Ii);
}

void RayPath::print(std::ofstream& FILE) const {
    char dlm = ',';
    FILE << "#Raypath " << i << ": frequency = " << f << " THz" << std::endl;
    FILE << "#x[um], y[um], z[um], t[ps], I[photons]" << std::endl;
    for (unsigned int j = 0; j < x.size(); j ++)
        FILE << x.at(j) << dlm << y.at(j) << dlm << z.at(j) << dlm << t.at(j) << dlm << I.at(j) << std::endl;
    FILE << std::endl << std::endl;
}
    
RayPath::RayPath() {
    i = i0; i0 += 1; f = 0; isFluorescence = false;
}

RayPath::RayPath(double v) {
    i = i0; i0 += 1; f = v; isFluorescence = false;
}

RayPath::RayPath(double v, double I0) {
    i = i0; i0 += 1; f = v; isFluorescence = false;
    I.push_back(I0);
    x.push_back(0);
    y.push_back(0);
    z.push_back(0);
    t.push_back(0);
}

RayPath::RayPath(double v, double I0, const vec& x0, double t0) {
    i = i0; i0 += 1; isFluorescence = false;
    I.push_back(I0);
    x.push_back(x0.X);
    y.push_back(x0.Y);
    z.push_back(x0.Z);
    t.push_back(t0);
}

