#include <vector>
#include <fstream>
#include <string>

#include "vec.h"

#ifndef _RAYPATH_H_
#define _RAYPATH_H_

///
//Raypath class
///
struct RayPath {
    double f;       //frequency
    unsigned long int i;
    std::vector<double> x, y, z, t, I;
    bool isFluorescence;
    
    static std::string ofbasename;
    
    void collide(const vec& xi, double t, double I);
    void print(std::ofstream& FILE) const;
    
    RayPath();
    RayPath(double v);
    RayPath(double v, double I0);
    RayPath(double v, double I0, const vec& x0, double t0);

    private:
        static unsigned long int i0;    //Number of actual instances
};

#endif
