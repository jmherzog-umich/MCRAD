#include <vector>
#include <cmath>
#include <string>
#include <iostream>
#include <iomanip>

#include "vec.h"

#ifndef __GRID_H
#define __GRID_H

using namespace std;

struct Grid {

    //Settings
    enum class CoordinateSystem { Cartesian = 0, Cylindrical = 1, Spherical = 2 };
  
    //Data  members
    vector<vector<double> > m;        //Data members, use a map to store multiple things
    vector<string> labels;            //Data member labels
    long unsigned int nx, ny, nz;     //Number of grid elements in each dimension
    double Lx, Ly, Lz;                //Dimensions of gridded region
    CoordinateSystem sys;

    //Data accessors
    long unsigned int ind1(const vec& x) const;
    long unsigned int ind2(const vec& x) const;
    long unsigned int ind3(const vec& x) const;
    long unsigned int sub2ind(long unsigned int i, long unsigned int j, long unsigned int k) const;
    double& at(int id, long unsigned int i, long unsigned int j, long unsigned int k);
    double& at(int id, const vec& x);
    double& at(const vec& x);
    double norm(int id, const vec& x) const;
    
    double Volume(long unsigned int i, long unsigned int j, long unsigned int k) const;
    double Volume(const vec& x) const;
    
    void newval(string label);
    void clear();
    
    double intersect(const vec& x, const vec& mu) const;
    
    void printGrid() const;
    void print() const;
    
    Grid();
    Grid(double x, double y, double z, int nx, int ny, int nz);
    Grid(double r, double z, int nr, int ntheta, int nz);
    Grid(double r, int nr, int ntheta, int nphi);
    
};

Grid::Grid() {
    Lx = Ly = Lz = 1;
    sys = CoordinateSystem::Cartesian;
    nx = ny = nz = 20;
}

Grid::Grid(double x, double y, double z, int nx=20, int ny=0, int nz=0) {
    //Initialize a spherical grid
    Lx = x; Ly = y; Lz = z;
    sys = CoordinateSystem::Cartesian;
    this->nx = nx;
    this->ny = (ny > 0) ? ny : nx;
    this->nz = (nz > 0) ? nz : nx;
}

Grid::Grid(double r, double z, int nr=20, int ntheta=0, int nz=0) {
    //Initialize a spherical grid
    Lx = r; Lz = z; Ly = 0;
    sys = CoordinateSystem::Cylindrical;
    this->nx = nr;
    this->ny = (ntheta > 0) ? ntheta : nr;
    this->nz = (nz > 0) ? nz : nr;
}

Grid::Grid(double r, int nr=20, int ntheta=0, int nphi=0) {
    //Initialize a spherical grid
    Lx = r; Ly = Lz = 0;
    sys = CoordinateSystem::Spherical;
    this->nx = nr;
    this->ny = (ntheta > 0) ? ntheta : nr;
    this->nz = (nphi > 0) ? nphi : nr;
}

double Grid::intersect(const vec& x, const vec& mu) const {
    //Get current coordinates
    long unsigned int ii, jj, kk;
    double dx1, dx2, dx3, a, b, c, d, s;
    a = b = c = d = s = 0;
    dx1 = dx2 = dx3 = -1;
    ii = ind1(x); jj = ind2(x); kk = ind3(x);
    
    //If we're outside the volume, exit
    switch (sys) {
        case CoordinateSystem::Cartesian:
            if ((ii == 0) or (jj == 0) or (kk == 0) or (ii > nx+1) or (jj > ny+1) or (kk > nz+1))
                return -1;
            break;
        case CoordinateSystem::Cylindrical:
            if ((kk == 0) or (ii > nx) or (kk > nz+1))
                return -1;
            break;
        case CoordinateSystem::Spherical:
            if (ii > nx)
                return -1;
            break;
    }
    
    //Now calculate distance to next volume element
    switch (sys) {
        case CoordinateSystem::Cartesian:
            dx1 = (((mu.X > 0) ? (ii+1) : ii) * Lx/nx - x.X)/mu.X;
            if (ny == 1)    //only a single giant slice, so no intersection
                dx2 = -1;
            else            //Intersect on each plane
                dx2 = (((mu.Y > 0) ? (jj+1) : jj) * Ly/ny - x.Y)/mu.Y;
            dx3 = (((mu.Z > 0) ? (kk+1) : kk) * Lz/nz - x.Z)/mu.Z;
            break;
        case CoordinateSystem::Cylindrical:
            //Some temporary variables
            a = x.X * mu.X + x.Y * mu.Y;
            b = mu.X*mu.X + mu.Y*mu.Y;
            
            //Calculate intersections
            //-Radial
            if (a > 0)                                  //Moving outwards
                c = pow((ii+1)* Lx/nx, 2) - x.r2();         // C > 0
            else                                        //Moving inwards
                if (ii == 0) {
                    c = pow(Lx/nx, 2) - x.r2();             // C > 0, moving through center
                } else
                    c = pow(ii* Lx/nx, 2) - x.r2();         // C < 0   
            if (c*b < -a*a)                               //Only get here if we don't intersect the circle
                dx1 = -1;
            else
                if ((ii == 0) and (a < 0))
                    dx1 = -a/b * (sqrt(1 + c*b/a/a) + 1);
                else
                    dx1 = a/b * (sqrt(1 + c*b/a/a) - 1);
            
            //Azimuthal
            if (ny == 1)    //Only a single shell, so no intersection
                dx2 = -1;
            else {          //Each phi intersection is a plane
                b = x.X * mu.Y - x.Y * mu.X;    //r x mu
                if (b > 0) {                //Moving in +theta direction
                    d = cos( 2*CONST_PI / ny * ((jj+1)%ny) );
                    c = sqrt(1 - d*d);
                } else if (b < 0) {         //Moving in -theta direction
                    d = cos( 2*CONST_PI / ny * jj );
                    c = sqrt(1 - d*d);
                } else
                    d = 0;
                if (b == 0)                 // perp to theta
                    dx2 = -1;
                else
                    dx2 = - (x.X * d + x.Y * c) / (mu.X * d + mu.Y * c);
            }
                
            //Z-direction
            dx3 = (((mu.Z > 0) ? (kk+1) : kk) * Lz/nz - x.Z)/mu.Z;
            break;
        case CoordinateSystem::Spherical:
            //Some temporary variables
            a = x.X * mu.X + x.Y * mu.Y + x.Z * mu.Z;
            
            //Calculate intersections
            //-Radial
            if (a > 0)                                  //Moving outwards
                c = pow((ii+1)* Lx/nx, 2) - x.r2();         // C > 0
            else                                        //Moving inwards
                if (ii == 0) {
                    c = pow(Lx/nx, 2) - x.r2();             // C > 0, moving through center
                } else
                    c = pow(ii* Lx/nx, 2) - x.r2();         // C < 0   
            if (c < -a*a)                               //Only get here if we don't intersect the circle
                dx1 = -1;
            else
                if ((ii == 0) and (a < 0))
                    dx1 = -a * (sqrt(1 + c/a/a) + 1);
                else
                    dx1 = a * (sqrt(1 + c/a/a) - 1);
            
            //Azimuthal
            if (ny == 1)    //Only a single shell, so no intersection
                dx2 = -1;
            else {          //Each phi intersection is a plane
                b = x.X * mu.Y - x.Y * mu.X;
                if (b > 0) {                //Moving in +theta direction
                    d = cos( 2*CONST_PI / ny * ((jj+1)%ny) );
                    c = sqrt(1 - d*d);
                } else if (b < 0) {         //Moving in -theta direction
                    d = cos( 2*CONST_PI / ny * jj );
                    c = sqrt(1 - d*d);
                } else
                    d = 0;
                if (b == 0)                 // perp to theta
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
                if (d > 0)              //Increasing theta
                    c = pow(cos(CONST_PI / nz * ((kk+1)%nz)),2);
                else
                    c = pow(cos(CONST_PI / nz * kk),2);
                
                //Now calculate some parameters
                s = 1-c;
                a = x.Z * mu.Z * s - (x.X * mu.X + x.Y * mu.Y) * c;
                b = (mu.X*mu.X + mu.Y*mu.Y)*c - mu.Z*mu.Z*s;
                c = (x.X*x.X + x.Y*x.Y)*c - x.Z*x.Z*s;
                
                //And evaluate the quadratic equation
                if (c*b/a/a > 1)        //No [real] solution
                    dx3 = -1;
                else {                  //Two real solutions
                    //Evalulate both roots and take the smaller [positive] one
                    s = a/b * (1 - sqrt(1 - c*b/a/a));
                    dx3 = a/b * (1 + sqrt(1 - c*b/a/a));
                    if ((s > 0) and (s < dx3))
                        dx3 = s;
                    if (dx3 < 0)
                        dx3 = -1;
                }
            }
            break;
    }
    
    //Now check boundary values
    switch (sys) {
        case CoordinateSystem::Cartesian:
            if (((ii == nx) and (mu.X > 0)) or ((ii == 1) and (mu.X < 0)))
                dx1 = -1;
            if (((jj == ny) and (mu.Y > 0)) or ((jj == 1) and (mu.Y < 0)))
                dx2 = -1;
            if (((kk == nz) and (mu.Z > 0)) or ((kk == 1) and (mu.Z < 0)))
                dx3 = -1;
            break;
        case CoordinateSystem::Cylindrical:
            if ((ii == nx-1) and (a > 0))
                dx1 = -1;
            if (((kk == nz) and (mu.Z > 0)) or ((kk == 1) and (mu.Z < 0)))
                dx3 = -1;
            break;
        case CoordinateSystem::Spherical:
            if ((ii == nx-1) and (a > 0))
                dx1 = -1;
            break;
    }
    
    //Now get the smallest dx value
    double dx = dx1;
    if ((dx < 0) or (dx2 < dx))
        dx = dx2;
    if ((dx < 0) or (dx3 < dx))
        dx = dx3;
    return dx;
}

void Grid::newval(string label) {
    m.push_back(vector<double>());
    labels.push_back(label);
}

void Grid::clear() {
    int N = 1;
    switch (sys) {
        case CoordinateSystem::Cartesian:
            N = (nx+2)*(ny+2)*(nz+2);
            break;
        case CoordinateSystem::Cylindrical:
            N = (nx+1)*ny*(nz+2);
            break;
        case CoordinateSystem::Spherical:
            N = (nx+1)*ny*nz;
            break;
    }
    for (unsigned long int i = 0; i < m.size(); i ++)
        m.at(i) = vector<double>(N, 0);
}

double& Grid::at(const vec& x) {
    return this->at(0, x);
}

double& Grid::at(int id, long unsigned int i, long unsigned int j, long unsigned int k) {
    return m.at(id).at(sub2ind(i,j,k));
}

double& Grid::at(int id, const vec& x) {
    long unsigned int ii, jj, kk;
    ii = ind1(x); jj = ind2(x); kk = ind3(x);
    return this->at(id, ii, jj, kk);
}

long unsigned int Grid::sub2ind(long unsigned int i, long unsigned int j, long unsigned int k) const {
    long unsigned int I = 0;
    switch (sys) {
        case CoordinateSystem::Cartesian:
            I = (nx+2)*(ny+2)*k + (nx+2)*j + i;
            break;
        case CoordinateSystem::Cylindrical:
            I = (nx+1)*ny*k + (nx+1)*j + i;
            break;
        case CoordinateSystem::Spherical:
            I = (nx+1)*ny*k + (nx+1)*j + i;
            break;
    }
    return I;
}

void Grid::printGrid() const {
    //Print the domain and number of data points
    switch (sys) {
        case CoordinateSystem::Cartesian:
            cout << "Cartesian simulation Domain (X,Y,Z): " << Lx << " um, " << Ly << " um," << Lz << " um" << endl;
            cout << "Grid points (nX,nY,nZ): ";
            break;
        case CoordinateSystem::Cylindrical:
            cout << "Cylindrical simulation Domain (R,Theta,Z): " << Lx << " um, " << "2pi rad, " << Lz << " um" << endl;
            cout << "Grid points (nR,nTheta,nZ): ";
            break;
        case CoordinateSystem::Spherical:
            cout << "Spherical simulation Domain (R,Theta,Phi): " << Lx << " um, " << "2pi rad, pi rad" << endl;
            cout << "Grid points (nR,nTheta,nPhi): "; 
            break;
    }
    cout << nx << ", " << ny << ", " << nz << "    Total elements: " << nx*ny*nz << endl;
    
    //Print the grid vectors
    switch (sys) {
        case CoordinateSystem::Cartesian:
            if (nx > 1) {
                cout << "X [um]:     ";
                for (long unsigned int i = 0; i < nx; i ++) cout << scientific << setw(18) << i * Lx/nx;
                cout << endl;
            } if (ny > 1) {
                cout << "Y [um]:     ";
                for (long unsigned int i = 0; i < ny; i ++) cout << scientific << setw(18) << i * Ly/ny;
                cout << endl;
            } if (nz > 1) {
                cout << "Z [um]:     ";
                for (long unsigned int i = 0; i < nz; i ++) cout << scientific << setw(18) << i * Lz/nz;
                cout << endl;
            }
            break;
        case CoordinateSystem::Cylindrical:
            if (nx > 1) {
                cout << "R [um]:     ";
                for (long unsigned int i = 0; i < nx; i ++) cout << scientific << setw(18) << i * Lx/nx;
                cout << endl;
            } if (ny > 1) {
                cout << "Theta [um]: ";
                for (long unsigned int i = 0; i < ny; i ++) cout << scientific << setw(18) << i * Ly/2/CONST_PI;
                cout << endl;
            } if (nz > 1) {
                cout << "Z [um]:     ";
                for (long unsigned int i = 0; i < nz; i ++) cout << scientific << setw(18) << i * Lz/nz;
                cout << endl;
            }
            break;
        case CoordinateSystem::Spherical:
            if (nx > 1) {
                cout << "R [um]:     ";
                for (long unsigned int i = 0; i < nx; i ++) cout << scientific << setw(18) << i * Lx/nx;
                cout << endl;
            } if (ny > 1) {
                cout << "Theta [um]: ";
                for (long unsigned int i = 0; i < ny; i ++) cout << scientific << setw(18) << i * Ly/2/CONST_PI;
                cout << endl;
            } if (nz > 1) {
                cout << "Phi [um]:   ";
                for (long unsigned int i = 0; i < nz; i ++) cout << scientific << setw(18) << i * Lz/CONST_PI;
                cout << endl;
            }
            break;
    }
    cout << endl;
}

void Grid::print() const {
    //Some initial parameters for the loops
    unsigned long int k0, j0, l0;
    double val = 0;
    k0 = j0 = l0 = 0;
    
    //Now loop through and print the grid
    for (unsigned long int i = 0; i < m.size(); i ++) {
        //print Key
        switch (sys) {
            case CoordinateSystem::Cartesian:
                if (Ly > 1)
                    cout << "[X, Y, Z]";
                else
                    cout << "[X, Z]";
                k0 = 1; j0 = 1; l0 = 1;
                break;
            case CoordinateSystem::Cylindrical:
                if (Ly > 1)
                    cout << "[R, Theta, Z]";
                else
                    cout << "[R, Z]";
                k0 = 1; j0 = 0; l0 = 0;
                break;
            case CoordinateSystem::Spherical:
                if (Ly > 1)
                    cout << "[R, Theta, Phi]";
                else
                    cout << "[R, Phi]";
                k0 = 0; j0 = 0; l0 = 0;
                break;
        }
    
        //Print the key name
        val = 0;
        for (unsigned long int z = 0; z < m.at(i).size(); z ++) val += m.at(i).at(z);
        cout << "   " << labels.at(i) << "   SUM=" << val << endl;
        
        //Now print the data
        for (unsigned long int k = k0; k < (k0+nz); k ++) {                //Blocks    
            for (unsigned long int l = l0; l < (l0+ny); l ++) {            //Rows
                for (unsigned long int j = j0; j < (j0+nx); j ++)          //Columns
                    cout << scientific << setw(18) << m.at(i).at(sub2ind(j, l, k));
            //Skip line after each row
            cout << endl;
            }
        //Skip an extra line after each block ONLY IF WE HAVE MULTIPLE LINES PER BLOCK
        if (ny > 1)
            cout << endl;
        }
    //Skip an extra line for fun
    cout << endl;   
    }
}

long unsigned int Grid::ind1(const vec& x) const {
    long unsigned int ii = 0;
    double xx;
    switch (sys) {
        //Clip domain to 1-n range, then use 0 and n+1 for out of bounds (total vector size is n+2) for x,y,z
        //Clip to 0 to n-1 range for r, and n for out of bounds (total size n+1)
        //Force theta and phi to 0-(n-1) range (total size n)
        case (CoordinateSystem::Cartesian):
            xx = x.X + Lx/2;
            ii = (xx > 0) ? (long unsigned int)(xx / Lx * nx)+1 : 0; if (ii > nx) ii = nx+1;
            break;
        case (CoordinateSystem::Cylindrical):
            ii = (long unsigned int)(x.r() / Lx * nx); if (ii >= nx) ii = nx;
            break;
        case (CoordinateSystem::Spherical):
            ii = (long unsigned int)(x.norm() / Lx * nx); if (ii >= nx) ii = nx;
            break;
    }
    return ii;
}

long unsigned int Grid::ind2(const vec& x) const {
    long unsigned int jj = 0;
    double yy;
    switch (sys) {
        //Clip domain to 1-n range, then use 0 and n+1 for out of bounds (total vector size is n+2) for x,y,z
        //Clip to 0 to n-1 range for r, and n for out of bounds (total size n+1)
        //Force theta and phi to 0-(n-1) range (total size n)
        case (CoordinateSystem::Cartesian):
            yy = x.Y + Ly/2;
            jj = (yy > 0) ? (long unsigned int)(yy / Ly * ny)+1 : 0; if (jj > ny) jj = ny+1;
            break;
        case (CoordinateSystem::Cylindrical):
            jj = (long unsigned int)(fmod(x.theta(), 2*CONST_PI)*ny/2/CONST_PI); if (jj >= ny) jj = ny-1;
            break;
        case (CoordinateSystem::Spherical):
            jj = (long unsigned int)(fmod(x.theta(), 2*CONST_PI)*ny/2/CONST_PI); if (jj >= ny) jj = ny-1;
            break;
    }
    return jj;
}

long unsigned int Grid::ind3(const vec& x) const {
    long unsigned int kk = 0;
    switch (sys) {
        //Clip domain to 1-n range, then use 0 and n+1 for out of bounds (total vector size is n+2) for x,y,z
        //Clip to 0 to n-1 range for r, and n for out of bounds (total size n+1)
        //Force theta and phi to 0-(n-1) range (total size n)
        case (CoordinateSystem::Cartesian):
            kk = (x.Z > 0) ? (long unsigned int)(x.Z / Lz * nz)+1 : 0; if (kk > nz) kk = nz+1;
            break;
        case (CoordinateSystem::Cylindrical):
            kk = (x.Z > 0) ? (long unsigned int)(x.Z / Lz * nz)+1 : 0; if (kk > nz) kk = nz+1;
            break;
        case (CoordinateSystem::Spherical):
            kk = (long unsigned int)(fmod(x.phi(), CONST_PI)*nz/CONST_PI); if (kk >= nz) kk = nz-1;
            break;
    }
    return kk;
}

double Grid::norm(int id, const vec& x) const {
    long unsigned int ii, jj, kk;
    ii = ind1(x); jj = ind2(x); kk = ind3(x);
    return m.at(id).at(sub2ind(ii, jj, kk)) / Volume(ii, jj, kk);
}

double Grid::Volume(const vec& x) const {
    long unsigned int ii, jj, kk;
    ii = ind1(x); jj = ind2(x); kk = ind3(x);
    return Volume(ii, jj, kk);
}

double Grid::Volume(long unsigned int i, long unsigned int j, long unsigned int k) const {
    double V = 1;
    switch (sys) {
        case (CoordinateSystem::Cartesian):
            V = (Lx/nx)*(Ly/ny)*(Lz/nz);
            break;
        case (CoordinateSystem::Cylindrical):
            V = 0;
            break;
        case (CoordinateSystem::Spherical):
            V = 0;
            break;
    }
    return V;
}
#endif
