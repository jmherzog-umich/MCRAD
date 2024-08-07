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

#ifndef __CYLGRID_H
#define __CYLGRID_H

using namespace std;

struct CylGrid : public Grid {
    //Data accessors
    virtual long unsigned int ind1(const vec& x) const {
        long unsigned int ii = (long unsigned int) floor(x.r() / Lx * nx);
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
        long unsigned int kk = (x.Z > 0) ? (long unsigned int) floor(x.Z / Lz * nz)+1 : 0;
        if (kk > nz)
            kk = nz+1;
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
        return ((double)i+0.5)*pow(Lx/nx,2) * 2*CONST_PI/ny * Lz/nz;
    }
    
    virtual vec rand(long unsigned int i, long unsigned int j, long unsigned int k) const {
        vec out;
        double r = sqrt(roll() * (2*i+1) + i*i)*Lx/nx;
        double sth = sin(2*CONST_PI*(j+roll())/ny);
        out.X = sqrt(1-sth*sth) * r;
        out.Y = sth * r;
        out.Z = (k-1 + roll())*Lz/nz;
        return out;
    }
    
    virtual bool contains(const vec& x) const {
        if (x.r2() > Lx*Lx)
            return false;
        if (x.Z < 0 or x.Z > Lz)
            return false;
        return true;
    }
    
    virtual unsigned long int geti0() const {
        return 0;
    }
    
    virtual unsigned long int getj0() const {
        return 0;
    }
    
    virtual unsigned long int getk0() const {
        return 1;
    }
    
    virtual long unsigned int ntotal() const {
        return (nx+1)*ny*(nz+2);
    }
    
    virtual double intersect(const vec& x, const vec& mu) const {
        //Get current index
        double dx1, dx2, dx3;
        long unsigned int ii = ind1(x);
        long unsigned int jj = ind2(x);
        long unsigned int kk = ind3(x);
        
        double a = x.X * mu.X + x.Y * mu.Y;
        double b = mu.X*mu.X + mu.Y*mu.Y;
        double c, d;
        
        //Check if we're outside
        if ((kk == 0) or (ii > nx) or (kk > nz+1))
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
        if (c*b < -a*a)                                 //Only get here if we don't intersect the circle
            dx1 = -1;
        else
            if ((ii == 0) and (a < 0))
                dx1 = -a/b * (sqrt(1 + c*b/a/a) + 1);
            else
                dx1 = a/b * (sqrt(1 + c*b/a/a) - 1);
        //Azimuthal
        if (ny == 1)                                    //Only a single shell, so no intersection
            dx2 = -1;
        else {                                          //Each phi intersection is a plane
            b = x.X * mu.Y - x.Y * mu.X;                //r x mu
            if (b > 0) {                                //Moving in +theta direction
                d = cos( 2*CONST_PI / ny * ((jj+1)%ny) );
                c = sqrt(1 - d*d);
            } else if (b < 0) {                         //Moving in -theta direction
                d = cos( 2*CONST_PI / ny * jj );
                c = sqrt(1 - d*d);
            } else
                d = 0;
            if (b == 0)                                 //perp to theta
                dx2 = -1;
            else
                dx2 = - (x.X * d + x.Y * c) / (mu.X * d + mu.Y * c);
        }
        //Z-direction
        dx3 = (((mu.Z > 0) ? (kk+1) : kk) * Lz/nz - x.Z)/mu.Z;
        
        //Check boundaries
        if ((ii == nx-1) and (a > 0))
            dx1 = -1;
        if (((kk == nz) and (mu.Z > 0)) or ((kk == 1) and (mu.Z < 0)))
            dx3 = -1;
        
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
    
        //Front surface
        tmp = 0; ss.str("");
        for (unsigned long int j = 0; j < ny; j ++) {
            for (unsigned long int k = 0; k < nx+1; k ++) {
                tmp += m.at(id).at(sub2ind(k, j, 0));
                ss << scientific << setw(18) << m.at(id).at(sub2ind(k, j, 0));
            }
            ss << endl;
        }
        oout << "    [R, Theta; Z <  0]: " << tmp << endl;
        if (tmp > CONST_EPS)
            oout << ss.str() << endl;
        
        //Back surface
        tmp = 0; ss.str("");
        for (unsigned long int j = 0; j < ny; j ++) {
            for (unsigned long int k = 0; k < nx+1; k ++) {
                tmp += m.at(id).at(sub2ind(k, j, nz+1));
                ss << scientific << setw(18) << m.at(id).at(sub2ind(k, j, nz+1));
            }
            ss << endl;
        }
        oout << "    [R, Theta; Z > Lz]: " << tmp << endl;
        if (tmp > CONST_EPS)
            oout << ss.str() << endl;
        
        //Top surface
        tmp = 0; ss.str("");
        for (unsigned long int j = 0; j < ny; j ++) {
            for (unsigned long int k = 1; k < nz+1; k ++) {
                tmp += m.at(id).at(sub2ind(nx+1, j, k));
                ss << scientific << setw(18) << m.at(id).at(sub2ind(nx+1, j, k));
            }
            ss << endl;
        }
        oout << "    [Z, Theta; R > L]: " << tmp << endl;
        if (tmp > CONST_EPS)
            oout << ss.str() << endl;
        
        //No corners here, because we let k go up to r > R above for front/back sides
        oout << endl << endl;
    }
    
    virtual void printGrid(ostream& oout) const {
        //Print type, domain, and number of points
        oout << "Cylindrical simulation Domain (R,Theta,Z): " << Lx << " um, " << "2pi rad, " << Lz << " um" << endl;
        oout << "Grid points (nR,nTheta,nZ): ";
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
            oout << "Z_center [um]:     ";
            for (long unsigned int i = 0; i < nz; i ++) oout << scientific << setw(18) << (i+0.5) * Lz/nz;
            oout << endl;
        }
        oout << endl;
    }
    
    virtual void printKeys(ostream& oout) const {
        if (ny > 1)
            oout << "[R, Theta, Z]";
        else
            oout << "[R, Z]";
    }
    
    CylGrid() {
        Lx = Lz = 1;
        nx = ny = nz = 20;
        Ly = 2.0*CONST_PI/(double)ny;
    }
    
    CylGrid(double r, double z, unsigned int nr=20, unsigned int ntheta=0, unsigned int nz=0) {
        Lx = r; Lz = z;
        this->nx = nr;
        this->ny = (ntheta > 0) ? ntheta : nr;
        this->nz = (nz > 0) ? nz : nr;
        Ly = 2.0*CONST_PI/(double)ny;
    }
    
    virtual void collideSide(const vec& x, const vec& mu, double& ds, int& reflect) const {
        //Some edge cases
        if (abs(mu.r2()) <= 0)
            return;

        //Some preliminary variables
        double a = mu.X*mu.X + mu.Y*mu.Y;
        double b = (mu.X*x.X + mu.Y*x.Y)/a;
        double c = (x.r2() - Lx*Lx)/a;
        if (c/a > b*b)
            return;
                
        //Calculate the roots
        double r1 = -sqrt(b*b - c) - b;
        double r2 = sqrt(b*b - c) - b;
        
        //Check the two intersections
        if (r1 < ds and r1 > 0) {
            reflect = 3;
            ds = r1;
        }
        if (r2 < ds and r2 > 0) {
            reflect = 3;
            ds = r2;
        }
    }
    
    virtual void boundXY(vec& x) const {
        if (x.r2() < Lx*Lx)
            return;
        vec u(x.X, x.Y, 0);
        u = u / u.norm();
        while ( x.r2() > Lx*Lx )
            x = x - u * 2 * Lx;
    }
};

#endif
