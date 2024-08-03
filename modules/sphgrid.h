#include <vector>
#include <cmath>
#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <utility>

#include "../grid.h"
#include "../vec.h"
#include "../utility.h"

#define CONST_PI            3.1415926535897932
#define CONST_C             299.792458
#define CONST_HBAR          1.054571817679489e-22
#define CONST_HK            47.992430733662212
#define CONST_APERY         1.2020569031595942
#define CONST_WIEN          1.5936242600400401
#define CONST_EPS           1e-10
#define CONST_PLANCKMAX     0.6476102378919149

#ifndef __SPHGRID_H
#define __SPHGRID_H

using namespace std;

struct SphGrid : public Grid {
    //Data accessors
    virtual long unsigned int ind1(const vec& x) const {
        long unsigned int ii = (long unsigned int) floor(x.norm() / Lx * nx);
        if (ii >= nx)
            ii = nx;
        return ii;
    }
    
    virtual long unsigned int ind2(const vec& x) const {
        long unsigned int jj = (long unsigned int) floor(fmod(x.theta(), 2*CONST_PI)*ny/2/CONST_PI);
        if (jj >= ny)
            jj = ny-1;
        return jj;
    }
    
    virtual long unsigned int ind3(const vec& x) const {
        long unsigned int kk = (long unsigned int) floor(fmod(x.phi(), CONST_PI)*nz/CONST_PI);
        if (kk >= nz)
            kk = nz-1;
        return kk;
    }
    
    virtual long unsigned int ind2sub1(long unsigned int i) const {
        return (i % ((nx+1)*ny)) % (nx+1);
    }
    
    virtual long unsigned int ind2sub2(long unsigned int i) const {
        return (i % ((nx+1)*ny)) / (nx+1);
    }
    
    virtual long unsigned int ind2sub3(long unsigned int i) const {
        return i / ((nx+1)*ny);
    }
    
    virtual long unsigned int sub2ind(long unsigned int i, long unsigned int j, long unsigned int k) const {
        return (nx+1)*ny*k + (nx+1)*j + i;
    }
    
    virtual double Volume(long unsigned int i, long unsigned int j, long unsigned int k) const {
        cerr << "SphGrid::Volume not implemented yet!" << endl;
        return 0;
    }
    
    virtual vec rand(long unsigned int i, long unsigned int j, long unsigned int k) const {
        vec out;
        double r = cbrt(roll() * (3*i*i+3*i+1) + i*i*i)*Lx/nx;
        double eps = cos(CONST_PI/nz*(k+roll()));
        double sth = sin(2*CONST_PI*(j+roll())/ny);
        out.X = r * sqrt(1-sth*sth) * eps;
        out.Y = r * sth * eps;
        out.Z = r * sqrt(1-eps*eps);
        return out;
    }
    
    virtual bool contains(const vec& x) const {
        if (x.r2() > Lx*Lx)
            return false;
        return true;
    }
    
    virtual unsigned long int geti0() const {
        cerr << "SphGrid::geti0 not implemented yet!" << endl;
        return 0;
    }
    
    virtual unsigned long int getj0() const {
        cerr << "SphGrid::getj0 not implemented yet!" << endl;
        return 0;
    }
    
    virtual unsigned long int getk0() const {
        cerr << "SphGrid::getk0 not implemented yet!" << endl;
        return 0;
    }
    
    virtual long unsigned int ntotal() const {
        return (nx+1)*ny*nz;
    }
    
    virtual double intersect(const vec& x, const vec& mu) const {
        //Get current index
        double dx1, dx2, dx3;
        long unsigned int ii = ind1(x);
        long unsigned int jj = ind2(x);
        long unsigned int kk = ind3(x);
        
        double a = x.X * mu.X + x.Y * mu.Y + x.Z * mu.Z;
        double b, c, d, s;
        
        //Check if we're outside
        if (ii > nx)
            return -1;
        
        //Calculate intersections
        //-Radial
        if (a > 0)                                      //Moving outwards
            c = pow((ii+1)* Lx/nx, 2) - x.r2();         // C > 0
        else                                            //Moving inwards
            if (ii == 0) {
                c = pow(Lx/nx, 2) - x.r2();             // C > 0, moving through center
            } else
                c = pow(ii* Lx/nx, 2) - x.r2();         // C < 0   
        if (c < -a*a)                                   //Only get here if we don't intersect the circle
            dx1 = -1;
        else
            if ((ii == 0) and (a < 0))
                dx1 = -a * (sqrt(1 + c/a/a) + 1);
            else
                dx1 = a * (sqrt(1 + c/a/a) - 1);
        
        //Azimuthal
        if (ny == 1)                                    //Only a single shell, so no intersection
            dx2 = -1;
        else {                                          //Each phi intersection is a plane
            b = x.X * mu.Y - x.Y * mu.X;
            if (b > 0) {                                //Moving in +theta direction
                d = cos( 2*CONST_PI / ny * ((jj+1)%ny) );
                c = sqrt(1 - d*d);
            } else if (b < 0) {                         //Moving in -theta direction
                d = cos( 2*CONST_PI / ny * jj );
                c = sqrt(1 - d*d);
            } else
                d = 0;
            if (b == 0)                                 // perp to theta
                dx2 = -1;
            else
                dx2 = - (x.X * d + x.Y * c) / (mu.X * d + mu.Y * c);
        }
            
        //Polar
        if (nz == 1)
            dx3 = -1;
        else {
            //Angles of the position vector, and component of motion in polar direction (d)
            c = x.cosphi2();
            s = 1-c;
            d = x.costheta();
            d = sqrt(c) * (mu.X * d + mu.Y * sqrt(1-d*d)) - mu.Z * sqrt(s);
            
            //Now find the cos(theta)^2 of the cone we need to intersect
            if (d > 0)                                  //Increasing theta
                c = pow(cos(CONST_PI / nz * ((kk+1)%nz)),2);
            else
                c = pow(cos(CONST_PI / nz * kk),2);
            
            //Now calculate some parameters
            s = 1-c;
            a = x.Z * mu.Z * s - (x.X * mu.X + x.Y * mu.Y) * c;
            b = (mu.X*mu.X + mu.Y*mu.Y)*c - mu.Z*mu.Z*s;
            c = (x.X*x.X + x.Y*x.Y)*c - x.Z*x.Z*s;
            
            //And evaluate the quadratic equation
            if (c*b/a/a > 1)                            //No [real] solution
                dx3 = -1;
            else {                                      //Two real solutions
                //Evalulate both roots and take the smaller [positive] one
                s = a/b * (1 - sqrt(1 - c*b/a/a));
                dx3 = a/b * (1 + sqrt(1 - c*b/a/a));
                if ((s > 0) and (s < dx3))
                    dx3 = s;
                if (dx3 < 0)
                    dx3 = -1;
            }
        }
        
        //Check boundaries
        if ((ii == nx-1) and (a > 0))
            dx1 = -1;
        
        //Return nearest intersection
        double dx = dx1;
        if ((dx < 0) or (dx2 < dx))
            dx = dx2;
        if ((dx < 0) or (dx3 < dx))
            dx = dx3;
        return dx;
    }
    
    virtual void printOutside(ostream& oout, int id) const {
        //Declare some variables needed later
        stringstream ss;    
        double tmp = 0;
    
        //Outer surface
        tmp = 0; ss.str("");
        for (unsigned long int j = 0; j < ny; j ++) {
            for (unsigned long int k = 0; k < nz; k ++) {
                tmp += m.at(id).at(sub2ind(nx+1, j, k));
                ss << scientific << setw(18) << m.at(id).at(sub2ind(nx+1, j, k));
            }
            ss << endl;
        }
        oout << "    [Theta, Phi; R > Lx]: " << tmp << endl;
        if (tmp > CONST_EPS)
            oout << ss.str() << endl;
        oout << endl << endl;
    }
    
    virtual void printGrid(ostream& oout) const {
        //Print type, domain, and number of points
        oout << "Spherical simulation Domain (R,Theta,Phi): " << Lx << " um, " << "2pi rad, pi rad" << endl;
        oout << "Grid points (nR,nTheta,nPhi): "; 
        oout << nx << ", " << ny << ", " << nz << "    Total internal elements: " << nx*ny*nz << endl;
        
        //Print grid vectors
        if (nx > 1) {
            oout << "R_center [um]:     ";
            for (long unsigned int i = 0; i < nx; i ++) oout << scientific << setw(18) << (i+0.5) * Lx/nx;
            oout << endl;
        } if (ny > 1) {
            oout << "Theta_center [rad]: ";
            for (long unsigned int i = 0; i < ny; i ++) oout << scientific << setw(18) << (i+0.5) * Ly;
            oout << endl;
        } if (nz > 1) {
            oout << "Phi_center [rad]:   ";
            for (long unsigned int i = 0; i < nz; i ++) oout << scientific << setw(18) << (i+0.5) * Lz;
            oout << endl;
        }
        oout << endl;
    }
    
    virtual void printKeys(ostream& oout) const {
        if (ny > 1)
            oout << "[R, Theta, Phi]";
        else
            oout << "[R, Phi]";
    }
    
    SphGrid() {
        Lx = 1;
        nx = ny = nz = 20;
        Ly = 2.0*CONST_PI/(double)ny;
        Lz = CONST_PI/(double)ny;
    }
    
    SphGrid(double r, unsigned int nr=20, unsigned int ntheta=0, unsigned int nphi=0) {
        Lx = r;
        this->nx = nr;
        this->ny = (ntheta > 0) ? ntheta : nr;
        this->nz = (nphi > 0) ? nphi : nr;
        Ly = 2.0*CONST_PI/(double)ny;
        Lz = CONST_PI/(double)ny;
    }
};

#endif
