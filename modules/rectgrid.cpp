#include <vector>
#include <cmath>
#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <utility>

#include "rectgrid.h"
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

    //Data accessors
long unsigned int RectGrid::ind1(const vec& x) const {
    double xx = x.X + Lx/2;
    long unsigned int ii = (xx > 0) ? (long unsigned int) std::floor(xx / Lx * nx)+1 : 0;
    if (ii > nx)
        ii = nx+1;
    return ii;
}

long unsigned int RectGrid::ind2(const vec& x) const {
    double yy = x.Y + Ly/2;
    long unsigned int jj = (yy > 0) ? (long unsigned int) std::floor(yy / Ly * ny)+1 : 0;
    if (jj > ny)
        jj = ny+1;
    return jj;
}

long unsigned int RectGrid::ind3(const vec& x) const {
    long unsigned int kk = (x.Z > 0) ? (long unsigned int) std::floor(x.Z / Lz * nz)+1 : 0;
    if (kk > nz)
        kk = nz+1;
    return kk;
}

long unsigned int RectGrid::ind2sub1(long unsigned int i) const {
    return (i % ((nx+2)*(ny+2))) % (nx+2);
}

long unsigned int RectGrid::ind2sub2(long unsigned int i) const {
    return (i % ((nx+2)*(ny+2))) / (nx+2);
}

long unsigned int RectGrid::ind2sub3(long unsigned int i) const {
    return i / ((nx+2)*(ny+2));
}

long unsigned int RectGrid::sub2ind(long unsigned int i, long unsigned int j, long unsigned int k) const {
    return (nx+2)*(ny+2)*k + (nx+2)*j + i;
}

double RectGrid::Volume(long unsigned int i, long unsigned int j, long unsigned int k) const {
    return (Lx/nx)*(Ly/ny)*(Lz/nz);
}

vec RectGrid::rand(long unsigned int i, long unsigned int j, long unsigned int k) const {
    vec out;
    out.X = (i-1 + util::roll())*Lx/nx - Lx/2;
    out.Y = (j-1 + util::roll())*Ly/ny - Ly/2;
    out.Z = (k-1 + util::roll())*Lz/nz;
    return out;
}

bool RectGrid::contains(const vec& x) const {
    if (std::abs(x.X) > Lx/2)
        return false;
    if (std::abs(x.Y) > Ly/2)
        return false;
    if (x.Z < 0 or x.Z > Lz)
        return false;
    return true;
}

unsigned long int RectGrid::geti0() const {
    return 1;
}

unsigned long int RectGrid::getj0() const {
    return 1;
}

unsigned long int RectGrid::getk0() const {
    return 1;
}

long unsigned int RectGrid::ntotal() const {
    return (nx+2)*(ny+2)*(nz+2);
}

double RectGrid::intersect(const vec& x, const vec& mu) const {
    //Get current index
    double dx1, dx2, dx3;
    long unsigned int ii = ind1(x);
    long unsigned int jj = ind2(x);
    long unsigned int kk = ind3(x);
    
    //If we're outside the grid, give up
    if ((ii == 0) or (jj == 0) or (kk == 0) or (ii > nx+1) or (jj > ny+1) or (kk > nz+1))
        return -1;
        
    //Calculate the three intersections (x,y,z)
    dx1 = (((mu.X > 0) ? (ii+1) : ii) * Lx/nx - x.X)/mu.X;
    dx2 = (((mu.Y > 0) ? (jj+1) : jj) * Ly/ny - x.Y)/mu.Y;
    dx3 = (((mu.Z > 0) ? (kk+1) : kk) * Lz/nz - x.Z)/mu.Z;
    
    //Check boundaries - if only one element, we have no collision
    if (((ii == nx) and (mu.X > 0)) or ((ii == 1) and (mu.X < 0)) or (nx==1))
        dx1 = -1;
    if (((jj == ny) and (mu.Y > 0)) or ((jj == 1) and (mu.Y < 0)) or (ny==1))
        dx2 = -1;
    if (((kk == nz) and (mu.Z > 0)) or ((kk == 1) and (mu.Z < 0)) or (nz==1))
        dx3 = -1;
    
    //Return nearest intersection    
    double dx = dx1;
    if ((dx < 0) or (dx2 < dx))
        dx = dx2;
    if ((dx < 0) or (dx3 < dx))
        dx = dx3;
    return dx;
}

void RectGrid::printOutside(std::ostream& oout, int id) const {
    //Declare some variables needed later
    std::stringstream ss;    
    double tmp = 0;

    //Front surface
    tmp = 0; ss.str("");
    for (unsigned long int j = 0; j < ny+2; j ++) {
        for (unsigned long int k = 0; k < nx+2; k ++) {
            tmp += m.at(id).at(sub2ind(k, j, 0));
            ss << std::scientific << std::setw(18) << m.at(id).at(sub2ind(k, j, 0));
        }
        ss << std::endl;
    }
    oout << "    [X, Y; Z <  0]: " << tmp << std::endl;
    if (tmp > CONST_EPS)
        oout << ss.str() << std::endl;
    
    //Back surface
    tmp = 0; ss.str("");
    for (unsigned long int j = 0; j < ny+2; j ++) {
        for (unsigned long int k = 0; k < nx+2; k ++) {
            tmp += m.at(id).at(sub2ind(k, j, nz+1));
            ss << std::scientific << std::setw(18) << m.at(id).at(sub2ind(k, j, nz+1));
        }
        ss << std::endl;
    }
    oout << "    [X, Y; Z > Lz]: " << tmp << std::endl;
    if (tmp > CONST_EPS)
        oout << ss.str() << std::endl;
    
    //Left surface
    tmp = 0; ss.str("");
    for (unsigned long int j = 1; j < nz+1; j ++) {
        for (unsigned long int k = 0; k < nx+2; k ++) {
            tmp += m.at(id).at(sub2ind(k, 0, j));
            ss << std::scientific << std::setw(18) << m.at(id).at(sub2ind(k, 0, j));
        }
        ss << std::endl;
    }
    oout << "    [X, Z; Y < -Ly/2]: " << tmp << std::endl;
    if (tmp > CONST_EPS)
        oout << ss.str() << std::endl;
    
    //Right surface
    tmp = 0; ss.str("");
    for (unsigned long int j = 1; j < nz+1; j ++) {
        for (unsigned long int k = 0; k < nx+2; k ++) {
            tmp += m.at(id).at(sub2ind(k, ny+1, j));
            ss << std::scientific << std::setw(18) << m.at(id).at(sub2ind(k, ny+1, j));
        }
        ss << std::endl;
    }
    oout << "    [X, Z; Y > +Ly/2]: " << tmp << std::endl;
    if (tmp > CONST_EPS)
        oout << ss.str() << std::endl;
    
    
    //Bottom surface
    tmp = 0; ss.str("");
    for (unsigned long int j = 1; j < ny+1; j ++) {
        for (unsigned long int k = 1; k < nz+1; k ++) {
            tmp += m.at(id).at(sub2ind(0, j, k));
            ss << std::scientific << std::setw(18) << m.at(id).at(sub2ind(0, j, k));
        }
        ss << std::endl;
    }
    oout << "    [Z, Y; X < -Lx/2]: " << tmp << std::endl;
    if (tmp > CONST_EPS)
        oout << ss.str() << std::endl;
    
    //Top surface
    tmp = 0; ss.str("");
    for (unsigned long int j = 1; j < ny+1; j ++) {
        for (unsigned long int k = 1; k < nz+1; k ++) {
            tmp += m.at(id).at(sub2ind(nx+1, j, k));
            ss << std::scientific << std::setw(18) << m.at(id).at(sub2ind(nx+1, j, k));
        }
        ss << std::endl;
    }
    oout << "    [Z, Y; X < +Lx/2]: " << tmp << std::endl;
    if (tmp > CONST_EPS)
        oout << ss.str() << std::endl;
    oout << std::endl << std::endl;
}

void RectGrid::printGrid(std::ostream& oout) const {
    //Print type, domain, and number of points
    oout << "Cartesian simulation Domain (X,Y,Z): " << Lx << " um, " << Ly << " um," << Lz << " um" << std::endl;
    oout << "Grid points (nX,nY,nZ): ";
    oout << nx << ", " << ny << ", " << nz << "    Total internal elements: " << nx*ny*nz << std::endl;
    
    //Print grid vectors
    if (nx > 1) {
        oout << "X_center [um]:     ";
        for (long unsigned int i = 0; i < nx; i ++) oout << std::scientific << std::setw(18) << (i+0.5) * Lx/nx - Lx/2;
        oout << std::endl;
    } if (ny > 1) {
        oout << "Y_center [um]:     ";
        for (long unsigned int i = 0; i < ny; i ++) oout << std::scientific << std::setw(18) << (i+0.5) * Ly/ny - Ly/2;
        oout << std::endl;
    } if (nz > 1) {
        oout << "Z_center [um]:     ";
        for (long unsigned int i = 0; i < nz; i ++) oout << std::scientific << std::setw(18) << (i+0.5) * Lz/nz;
        oout << std::endl;
    }
    oout << std::endl;
}

void RectGrid::printKeys(std::ostream& oout) const {
    if (ny > 1)
        oout << "[X, Y, Z]";
    else
        oout << "[X, Z]";
}

vec RectGrid::getNormal(int reflect, const vec& x, const vec& mu, vec& mur) {
    mur = mu;
    switch (reflect) {
        default:
        case 1:
            mur.Z = -mur.Z;
            return vec(0,0,-1);
        case 2:
            mur.Z = -mur.Z;
            return vec(0,0,1);
        case 3:
            mur.X = -mur.X;
            return vec(std::copysign(1,-mu.X),0,0);
        case 4:
            mur.Y = -mur.Y;
            return vec(0,std::copysign(1,-mu.Y),0);
    }
}

RectGrid::RectGrid() {
    Lx = Ly = Lz = 1;
    nx = ny = nz = 20;
}

RectGrid::RectGrid(double x, double y, double z, unsigned int nx, unsigned int ny, unsigned int nz) {
    Lx = x; Ly = y; Lz = z;
    this->nx = nx;
    this->ny = (ny > 0) ? ny : nx;
    this->nz = (nz > 0) ? nz : nx;
}
