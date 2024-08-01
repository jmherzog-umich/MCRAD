#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>

#include "vec.h"
#include "photon.h"

#ifndef __CAMERA_H
#define __CAMERA_H

#define CONST_PI            3.1415926535897932
#define CONST_C             299.792458
#define CONST_HBAR          1.054571817679489e-22
#define CONST_HK            47.992430733662212
#define CONST_APERY         1.2020569031595942
#define CONST_WIEN          1.5936242600400401
#define CONST_EPS           1e-10
#define CONST_PLANCKMAX     0.6476102378919149

using namespace std;

struct Camera {

    public:
        double f;           //Lens f/number
        double F;           //Lens focal length in simulation coordinates
        double So;          //Object plane focal distance
        
        double M() const;         //Magnification
        double cosTheta() const;  //Cosine of collection half angle
        double Omega() const;     //Theoretical collection solid angle
        double Si() const;        //Image plane distance
        double D() const;         //Lens diameter
        double t() const;         //Depth of focus
        
        vec xo() const;     //Object plane center point

        double Lpx;         //Pixel pitch in image (sensor) plane
        int nx;             //Number of pixels in grid (horizontal)
        int ny;             //Number of pixels in grid (vertical)

        vec pos;            //Location of lens
        vec mu;             //Direction of lens
        
        vector<double> I;   //Intensity data
        
        void image(Photon& p, double n=1);     //Take photon and find where it strikes lens > sensor, then increment
        
        void setup();
        
        Camera();
        Camera(int n, double Lpx, bool interference = false);
        Camera(int n, double Lpx, const vec& pos, const vec& u, bool interference = false);
        Camera(int n, double Lpx, double M, double fn, double fl, bool interference = false);
        Camera(int n, double Lpx, double M, double fn, double fl, const vec& pos, const vec& mu, bool interference = false);
        
        void print(ostream& oout) const;
        void printGrid(ostream& oout) const;
        void printSetup(ostream& oout) const;
    
        ///TODO: IMPLEMENT THESE
        double fmin, fmax;  //Cutoff bands (in THz!) for filter. If both are zero, the filter is off.
        double ton, toff;   //Gate time of the camera (only include photons that strike within this time range)
        bool diffraction;   //If on, include random spreading at aperture
        bool interfere;     //If on, include phase in sums
        bool subjective;    //If accounting for interference, detect the subjective (interfere at sensor) rather than objective (interfere at object plane) pattern
    
    //Cached values
    private:
        vec _xo;
        double _M;
        double _cost;
        double _omega;
        double _Si;
        double _D;
        double _t;
        
        vec right;
        vec up;
};

void Camera::image(Photon& p, double n) {
    
    //Screen by ray direction (if ray is moving away from lens, nothing can happen)
    double ct = p.mu.dot(mu);
    if (ct >= 0)
        return;
    
    //Figure out how far it is to the object and lens planes
    vec dx = xo() - p.x;            //Distance from photon to object plane
    double s = dx.dot(mu) / ct;     //Travel distance from photon to object plane
    double sl = s - So / ct;        //Travel distance from photon to lens plane
    double r;
    
    //Figure out how far from the lens center the ray passes - if r > D/2, we're outside the aperture
    vec vl = (p.x + p.mu * sl) - pos;
    //cerr << "I = " << p.W << "     " << vl << endl;
    r = vl.norm();
    if (r > D()/2)
        return;
    
    //Figure out position projection in object plane
    double xi, yi, M, si, ri;
    vec vo = p.x + p.mu*s;
    M = this->M();
    xi = vo.dot(right) * M;
    yi = vo.dot(up) * M;
    
    //Store the Camera on the sensor
    int binx = (int)(xi/Lpx) + nx;
    int biny = (int)(yi/Lpx) + ny;
    if ((binx < 0) or (binx > 2*nx) or (biny < 0) or (biny > 2*nx))
        return;
    if (interfere) {
        if (subjective) {
            ri = hypot(xi,yi);
            si = sqrt(Si()*Si() + pow(ri*(1.0-M)/M - So * sqrt(1-ct*ct)/ct, 2));
            p.t += (sl + si) / CONST_C * n;
        }
        I.at((2*nx+1)*biny + binx) += sqrt(p.W) * p.phi();
    } else
        I.at((2*nx+1)*biny + binx) += p.W;
}

void Camera::setup() {
    _xo = pos + mu * So;
    _D = F/f;
    _M = F/(F-So);
    _cost = 2*So/sqrt(_D*_D + 4*So*So);
    _omega = 2*CONST_PI * (1.0 - _cost);
    _Si = 1.0/(1.0/F - 1.0/So);
    _t = 2 * f * Lpx * _Si / F;
    I = vector<double>((2*nx+1)*(2*ny+1),0);
    right = mu.perp(CONST_PI/2);
    up = mu.perp(0);
};

vec Camera::xo() const {
    return _xo;
}

double Camera::t() const {
    return _t;
}

double Camera::D() const {
    return _D;
}

double Camera::M() const {
    return _M;
}

double Camera::Omega() const {
    return _omega;
}

double Camera::Si() const {
    return _Si;
}

double Camera::cosTheta() const {
    return _cost;
}

Camera::Camera() {
    F = 50000;
    f = 1.2;
    So = 100000;
    Lpx = 30;
    nx = 5;
    ny = 5;
    mu = vec(-1,0,0);
    pos = vec(100000, 0, 0);
    interfere = false;
    this->setup();
}

Camera::Camera(int n, double Lpx, bool interference) {
    F = 50000;
    f = 1.2;
    So = 100000;
    this->Lpx = Lpx;
    nx = n;
    ny = n;
    mu = vec(-1,0,0);
    pos = vec(100000, 0, 0);
    interfere = interference;
    this->setup();
}

Camera::Camera(int n, double Lpx, const vec& pos, const vec& mu, bool interference) {
    F = 50000;
    f = 1.2;
    So = 100000;
    this->Lpx = Lpx;
    nx = n;
    ny = n;
    this->mu = mu;
    this->pos = pos;
    interfere = interference;
    this->setup();
}

Camera::Camera(int n, double Lpx, double M, double fn, double fl, bool interference) {
    F = fl;
    f = fn;
    this->So = (M-1.0)/M * F;
    this->Lpx = Lpx;
    nx = n;
    ny = n;
    mu = vec(-1,0,0);
    pos = vec(So, 0, 0);
    interfere = interference;
    this->setup();
}

Camera::Camera(int n, double Lpx, double M, double fn, double fl, const vec& pos, const vec& mu, bool interference) {
    F = fl;
    f = fn;
    So = (M-1.0)/M * F;
    this->Lpx = Lpx;
    nx = n;
    ny = n;
    this->mu = mu;
    this->pos = pos;
    interfere = interference;
    this->setup();
}

void Camera::print(ostream& oout) const {
    for (int i = 0; i < 2*ny+1; i ++) {
        for (int j = 0; j < 2*nx + 1; j ++)
            oout << scientific << setw(18) << I.at((2*nx+1)*i + j);
        oout << endl;
    }
    oout << endl;
}

void Camera::printGrid(ostream& oout) const {
    oout << "NOTE: CAMERA SIMULATION MODULE IS NOT COMPLETE. DO NOT RELY ON THESE RESULTS. " << endl;
    oout << "Image domain (W,H):   " << (2*nx+1)*Lpx/abs(M()) << " x " << (2*ny+1)*Lpx/abs(M()) << " um";
    oout << endl << "  X:     ";
    for (int i = -nx; i <= nx; i ++) oout << scientific << setw(18) << i*Lpx/abs(M());
    oout << endl << "  Y:     ";
    for (int i = -ny; i <= ny; i ++) oout << scientific << setw(18) << i*Lpx/abs(M());
    oout << endl << endl;
}

void Camera::printSetup(ostream& oout) const {
    vec right = mu.perp(CONST_PI/2);
    vec up = mu.perp(0);
    
    oout << "NOTE: CAMERA SIMULATION MODULE IS NOT COMPLETE. DO NOT RELY ON THESE RESULTS. " << endl;
    oout << "Camera lens focal length: " << fixed << F/1000 << " mm" << endl;
    oout << "Camera lens f-number: " << fixed << f << endl;
    oout << "Camera lens magnification: " << fixed << M() << endl;
    oout << "Object-plane pixel size: " << fixed << Lpx / abs(M()) << " um" << endl;
    oout << "Pixel size: " << fixed << Lpx << " um" << endl;
    oout << "Pixel count: " << fixed << (2*nx+1) << " x " << (2*ny+1) << endl;
    oout << "Detector size: " << fixed << (2*nx+1)*Lpx/1000 << " x " << (2*ny+1)*Lpx/1000 << " mm" << endl;
    oout << "ROI size: " << fixed << (2*nx+1)*Lpx/1000/abs(M()) << " x " << (2*ny+1)*Lpx/1000/abs(M()) << " mm" << endl;
    oout << "Camera lens diameter: " << fixed << D()/1000 << " mm" << endl;
    oout << "Collection half-angle: " << fixed << acos(cosTheta()) << ", " << acos(cosTheta())*180/CONST_PI << " rad, deg" << endl;
    oout << "Collection solid-angle: " << fixed << Omega() << " sr" << endl;
    oout << "Collection fraction: " << fixed << Omega()/4/CONST_PI << endl;
    oout << "Object distance: " << fixed << So / 1000 << " mm" << endl;
    oout << "Image distance: " << fixed << -So * M() / 1000 << " mm" << endl;
    oout << "Depth of focus: " << fixed << t()/1000 << " mm" << endl;
    oout << "Lens position: " << "<" << fixed << pos.X/1000 << ", " << pos.Y/1000 << ", " << pos.Z/1000 << "> mm" << endl;
    oout << "Lens normal: " << "<" << fixed << mu.X << ", " << mu.Y << ", " << mu.Z << ">" << endl;
    oout << "Camera x-direction: " << "<" << fixed << right.X << ", " << right.Y << ", " << right.Z << ">" << endl;
    oout << "Camera y-direction: " << "<" << fixed << up.X << ", " << up.Y << ", " << up.Z << ">" << endl;
    oout << "<NOT IMPLEMENTED YET> Filter bands: " << fmin << "   to   " << fmax << " THz" << endl;
    oout << "<NOT IMPLEMENTED YET> Gate time: " << ton << "   to   " << toff << " ps" << endl;
    oout << "<NOT IMPLEMENTED YET> Diffraction: " << ((diffraction) ? "True" : "False") << endl;
    oout << "<NOT IMPLEMENTED YET> Interference: " << ((interfere) ? "True" : "False") << endl;
    if (interfere)
        oout << "      " << ((subjective) ? "Subjective" : "Objective") << " interference pattern" << endl;
    oout << endl;
}

#endif
