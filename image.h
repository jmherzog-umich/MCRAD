#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>

#include "vec.h"
#include "photon.h"

#ifndef __IMAGE_H
#define __IMAGE_H

using namespace std;

struct Image {

    double Lpx;         //Pixel pitch in image (sensor) plane
    int nx;             //Number of pixels in grid (horizontal)
    int ny;             //Number of pixels in grid (vertical)
    
    bool interfere;     //Whether to include interference in the calculation or not
    
    vector<double> I;   //Intensity data
    
    void image(double x, double y, double f, double phi=1);     //Take photon and find where it strikes lens > sensor, then increment
    void clear();
    
    Image();
    Image(int n, double Lpx, bool interference = false);
    Image(int nx, int ny, double Lpx, bool interference = false);
    
    void print() const;
    void printGrid() const;
};

void Image::image(double x, double y, double f, double phase) {
    int binx = (int) floor((x/Lpx) + nx);
    int biny = (int) floor((y/Lpx) + ny);
    if ((binx < 0) or (binx > 2*nx) or (biny < 0) or (biny > 2*nx))
        return;
    if (interfere)
        I.at((2*nx+1)*biny + binx) += sqrt(f) * phase;
    else
        I.at((2*nx+1)*biny + binx) += f;
}

void Image::clear() {
    I = vector<double>((2*nx+1)*(2*ny+1),0);
}

Image::Image() {
    nx = 64;
    ny = 64;
    Lpx = 100;
    interfere = false;
    clear();
}

Image::Image(int n, double Lpx, bool interference) {
    nx = n;
    ny = n;
    this->Lpx = Lpx;
    interfere = interference;
    clear();
}

Image::Image(int nx, int ny, double Lpx, bool interference) {
    this->nx = nx;
    this->ny = ny;
    this->Lpx = Lpx;
    interfere = interference;
    clear();
}

void Image::print() const {
    for (int i = 0; i < 2*ny+1; i ++) {
        for (int j = 0; j < 2*nx + 1; j ++)
            cout << scientific << setw(18) << I.at((2*nx+1)*i + j);
        cout << endl;
    }
    cout << endl;
}

void Image::printGrid() const {
    cout << "Image domain (W,H):   " << (2*nx+1)*Lpx << " x " << (2*ny+1)*Lpx << " um";
    cout << endl << "  X:     ";
    for (int i = -nx; i <= nx; i ++) cout << scientific << setw(18) << i*Lpx;
    cout << endl << "  Y:     ";
    for (int i = -ny; i <= ny; i ++) cout << scientific << setw(18) << i*Lpx;
    cout << endl << endl;
}

#endif
