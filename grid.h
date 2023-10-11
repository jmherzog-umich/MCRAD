#include <vector>
#include <cmath>
#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <utility>

#include "vec.h"
#include "utility.h"

#define CONST_EPS 1e-10

#ifndef __GRID_H
#define __GRID_H

using namespace std;

struct Grid {

    //Grid iterator
    struct iterator {
        //private:
            unsigned long int ix, iy, iz;
            const Grid& g;
            bool DONE;
        public:
            iterator(const Grid& p);
            iterator& operator++(int);
            unsigned long int index() const;
            bool end() const;
    };

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
    long unsigned int ind2sub1(long unsigned int i) const;
    long unsigned int ind2sub2(long unsigned int i) const;
    long unsigned int ind2sub3(long unsigned int i) const;
    long unsigned int sub2ind(long unsigned int i, long unsigned int j, long unsigned int k) const;

    //Simple iterator alternative
    long unsigned int ncell() const;
    iterator begin() const;

    double& at(int id, long unsigned int i, long unsigned int j, long unsigned int k);
    double& at(int id, long unsigned int i);
    double& at(int id, const vec& x);
    double& at(const vec& x);
    
    double norm(int id, const vec& x) const;
    
    double Volume(long unsigned int i) const;
    double Volume(long unsigned int i, long unsigned int j, long unsigned int k) const;
    double Volume(const vec& x) const;
    
    vec rand(long unsigned int i, long unsigned int j, long unsigned int k) const;
    vec rand(long unsigned int i) const;
    bool contains(const vec& x) const;
    
    bool isempty(int id) const;
    double sum(int id) const;
    bool less(int id, double val) const;
    
    vector<double> cum(int id, double& sum) const;
    vector<double> cum(int id) const;
    
    void newval(string label);
    void clear();
    
    double intersect(const vec& x, const vec& mu) const;
    
    void printOutside(int id) const;
    void printGrid() const;
    void print() const;
    
    Grid();
    Grid(double x, double y, double z, unsigned int nx, unsigned int ny, unsigned int nz);
    Grid(double r, double z, unsigned int nr, unsigned int ntheta, unsigned int nz);
    Grid(double r, unsigned int nr, unsigned int ntheta, unsigned int nphi);
    
};

Grid::Grid() {
    Lx = Ly = Lz = 1;
    sys = CoordinateSystem::Cartesian;
    nx = ny = nz = 20;
}

Grid::Grid(double x, double y, double z, unsigned int nx=20, unsigned int ny=0, unsigned int nz=0) {
    //Initialize a spherical grid
    Lx = x; Ly = y; Lz = z;
    sys = CoordinateSystem::Cartesian;
    this->nx = nx;
    this->ny = (ny > 0) ? ny : nx;
    this->nz = (nz > 0) ? nz : nx;
}

Grid::Grid(double r, double z, unsigned int nr=20, unsigned int ntheta=0, unsigned int nz=0) {
    //Initialize a spherical grid
    Lx = r; Lz = z; Ly = 0;
    sys = CoordinateSystem::Cylindrical;
    this->nx = nr;
    this->ny = (ntheta > 0) ? ntheta : nr;
    this->nz = (nz > 0) ? nz : nr;
}

Grid::Grid(double r, unsigned int nr=20, unsigned int ntheta=0, unsigned int nphi=0) {
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

double& Grid::at(int id, long unsigned int i) {
    return this->m.at(id).at(i);
}

double& Grid::at(const vec& x) {
    return this->at(0, x);
}

double& Grid::at(int id, long unsigned int i, long unsigned int j, long unsigned int k) {
    return this->m.at(id).at(sub2ind(i,j,k));
}

double& Grid::at(int id, const vec& x) {
    long unsigned int ii, jj, kk;
    ii = ind1(x); jj = ind2(x); kk = ind3(x);
    return this->at(id, ii, jj, kk);
}

long unsigned int Grid::ind2sub1(long unsigned int i) const {
    switch (sys) {
        case (CoordinateSystem::Cartesian):
            return (i % ((nx+2)*(ny+2))) % (nx+2);
        case (CoordinateSystem::Cylindrical):
            return 0;
        case (CoordinateSystem::Spherical):
            return 0;
    }
    return 0;
}

long unsigned int Grid::ind2sub2(long unsigned int i) const {
    switch (sys) {
        case (CoordinateSystem::Cartesian):
            return (i % ((nx+2)*(ny+2))) / (nx+2);
        case (CoordinateSystem::Cylindrical):
            return 0;
        case (CoordinateSystem::Spherical):
            return 0;
    }
    return 0;
}

long unsigned int Grid::ind2sub3(long unsigned int i) const {
    switch (sys) {
        case (CoordinateSystem::Cartesian):
            return i / ((nx+2)*(ny+2));
        case (CoordinateSystem::Cylindrical):
            return 0;
        case (CoordinateSystem::Spherical):
            return 0;
    }
    return 0;
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
                cout << "X_center [um]:     ";
                for (long unsigned int i = 0; i < nx; i ++) cout << scientific << setw(18) << (i+0.5) * Lx/nx - Lx/2;
                cout << endl;
            } if (ny > 1) {
                cout << "Y_center [um]:     ";
                for (long unsigned int i = 0; i < ny; i ++) cout << scientific << setw(18) << (i+0.5) * Ly/ny - Ly/2;
                cout << endl;
            } if (nz > 1) {
                cout << "Z_center [um]:     ";
                for (long unsigned int i = 0; i < nz; i ++) cout << scientific << setw(18) << (i+0.5) * Lz/nz;
                cout << endl;
            }
            break;
        case CoordinateSystem::Cylindrical:
            if (nx > 1) {
                cout << "R_center [um]:     ";
                for (long unsigned int i = 0; i < nx; i ++) cout << scientific << setw(18) << (i+0.5) * Lx/nx;
                cout << endl;
            } if (ny > 1) {
                cout << "Theta_center [um]: ";
                for (long unsigned int i = 0; i < ny; i ++) cout << scientific << setw(18) << (i+0.5) * Ly/2/CONST_PI;
                cout << endl;
            } if (nz > 1) {
                cout << "Z_center [um]:     ";
                for (long unsigned int i = 0; i < nz; i ++) cout << scientific << setw(18) << (i+0.5) * Lz/nz;
                cout << endl;
            }
            break;
        case CoordinateSystem::Spherical:
            if (nx > 1) {
                cout << "R_center [um]:     ";
                for (long unsigned int i = 0; i < nx; i ++) cout << scientific << setw(18) << (i+0.5) * Lx/nx;
                cout << endl;
            } if (ny > 1) {
                cout << "Theta_center [um]: ";
                for (long unsigned int i = 0; i < ny; i ++) cout << scientific << setw(18) << (i+0.5) * Ly/2/CONST_PI;
                cout << endl;
            } if (nz > 1) {
                cout << "Phi_center [um]:   ";
                for (long unsigned int i = 0; i < nz; i ++) cout << scientific << setw(18) << (i+0.5) * Lz/CONST_PI;
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
        //Don't print if the grid is empty
        if (isempty(i))
            continue;
    
        //print Key
        switch (sys) {
            case CoordinateSystem::Cartesian:
                if (ny > 1)
                    cout << "[X, Y, Z]";
                else
                    cout << "[X, Z]";
                k0 = 1; j0 = 1; l0 = 1;
                break;
            case CoordinateSystem::Cylindrical:
                if (ny > 1)
                    cout << "[R, Theta, Z]";
                else
                    cout << "[R, Z]";
                k0 = 1; j0 = 0; l0 = 0;
                break;
            case CoordinateSystem::Spherical:
                if (ny > 1)
                    cout << "[R, Theta, Phi]";
                else
                    cout << "[R, Phi]";
                k0 = 0; j0 = 0; l0 = 0;
                break;
        }
    
        //Print the key name
        val = 0;
        for (unsigned long int z = 0; z < m.at(i).size(); z ++) val += m.at(i).at(z);
        cout << "   " << labels.at(i) << "   SUM=" << val << "   (" << sum(i) << ")" << endl;
        
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
    //Now print the elements outside the grid
    printOutside(i);
    
    //Skip an extra line for fun
    cout << endl;   
    }
}

void Grid::printOutside(int id) const {
    //Initialize some stuff
    stringstream ss;    
    double tmp = 0;
    
    //Prepare to print
    cout << endl;
    switch (sys) {
        //Now switch based on CS type
        case CoordinateSystem::Cartesian:
        
            //Front surface
            tmp = 0; ss.str("");
            for (unsigned long int j = 1; j < ny+1; j ++) {
                for (unsigned long int k = 1; k < nx+1; k ++) {
                    tmp += m.at(id).at(sub2ind(k, j, 0));
                    ss << scientific << setw(18) << m.at(id).at(sub2ind(k, j, 0));
                }
                ss << endl;
            }
            cout << "    [X, Y; Z <  0]: " << tmp << endl;
            if (tmp > CONST_EPS)
                cout << ss.str() << endl;
            
            //Back surface
            tmp = 0; ss.str("");
            for (unsigned long int j = 1; j < ny+1; j ++) {
                for (unsigned long int k = 1; k < nx+1; k ++) {
                    tmp += m.at(id).at(sub2ind(k, j, nz+1));
                    ss << scientific << setw(18) << m.at(id).at(sub2ind(k, j, nz+1));
                }
                ss << endl;
            }
            cout << "    [X, Y; Z > Lz]: " << tmp << endl;
            if (tmp > CONST_EPS)
                cout << ss.str() << endl;
            
            //Left surface
            tmp = 0; ss.str("");
            for (unsigned long int j = 1; j < nz+1; j ++) {
                for (unsigned long int k = 1; k < nx+1; k ++) {
                    tmp += m.at(id).at(sub2ind(k, 0, j));
                    ss << scientific << setw(18) << m.at(id).at(sub2ind(k, 0, j));
                }
                ss << endl;
            }
            cout << "    [X, Z; Y < -Ly/2]: " << tmp << endl;
            if (tmp > CONST_EPS)
                cout << ss.str() << endl;
            
            //Right surface
            tmp = 0; ss.str("");
            for (unsigned long int j = 1; j < nz+1; j ++) {
                for (unsigned long int k = 1; k < nx+1; k ++) {
                    tmp += m.at(id).at(sub2ind(k, ny+1, j));
                    ss << scientific << setw(18) << m.at(id).at(sub2ind(k, ny+1, j));
                }
                ss << endl;
            }
            cout << "    [X, Z; Y > +Ly/2]: " << tmp << endl;
            if (tmp > CONST_EPS)
                cout << ss.str() << endl;
            
            
            //Bottom surface
            tmp = 0; ss.str("");
            for (unsigned long int j = 1; j < ny+1; j ++) {
                for (unsigned long int k = 1; k < nz+1; k ++) {
                    tmp += m.at(id).at(sub2ind(0, j, k));
                    ss << scientific << setw(18) << m.at(id).at(sub2ind(0, j, k));
                }
                ss << endl;
            }
            cout << "    [Z, Y; X < -Lx/2]: " << tmp << endl;
            if (tmp > CONST_EPS)
                cout << ss.str() << endl;
            
            //Top surface
            tmp = 0; ss.str("");
            for (unsigned long int j = 1; j < ny+1; j ++) {
                for (unsigned long int k = 1; k < nz+1; k ++) {
                    tmp += m.at(id).at(sub2ind(nx+1, j, k));
                    ss << scientific << setw(18) << m.at(id).at(sub2ind(nx+1, j, k));
                }
                ss << endl;
            }
            cout << "    [Z, Y; X < +Lx/2]: " << tmp << endl;
            if (tmp > CONST_EPS)
                cout << ss.str() << endl;
            
            //Now the corners
            tmp = m.at(id).at(sub2ind(0, 0, 0)) + m.at(id).at(sub2ind(0, 0, nz+1))
                + m.at(id).at(sub2ind(0, ny+1, 0)) + m.at(id).at(sub2ind(0, ny+1, nz+1))
                + m.at(id).at(sub2ind(nx+1, 0, 0)) + m.at(id).at(sub2ind(nx+1, 0, nz+1))
                + m.at(id).at(sub2ind(nx+1, ny+1, 0)) + m.at(id).at(sub2ind(nx+1, ny+1, nz+1));
            cout << endl << "[-,-,-], [-,-,+], [-,+,-], [-,+,+], [+,-,-], [+,-,+], [+,+,-], [+,+,+]: " << tmp << endl;
            if (tmp > CONST_EPS) {
                cout << scientific << setw(18) << m.at(id).at(sub2ind(0, 0, 0));
                cout << scientific << setw(18) << m.at(id).at(sub2ind(0, 0, nz+1));
                cout << scientific << setw(18) << m.at(id).at(sub2ind(0, ny+1, 0));
                cout << scientific << setw(18) << m.at(id).at(sub2ind(0, ny+1, nz+1));
                cout << scientific << setw(18) << m.at(id).at(sub2ind(nx+1, 0, 0));
                cout << scientific << setw(18) << m.at(id).at(sub2ind(nx+1, 0, nz+1));
                cout << scientific << setw(18) << m.at(id).at(sub2ind(nx+1, ny+1, 0));
                cout << scientific << setw(18) << m.at(id).at(sub2ind(nx+1, ny+1, nz+1));
                cout << endl;
            }
            break;
        case CoordinateSystem::Cylindrical:
            
            
            break;
        case CoordinateSystem::Spherical:
            for (unsigned long int j = 0; j < ny; j ++) {
                for (unsigned long int k = 0; k < nz; k ++) {
                    tmp += m.at(id).at(sub2ind(nx, j, k));
                    ss << scientific << setw(18) << m.at(id).at(sub2ind(nx, j, k));
                }
                ss << endl;
            }
            cout << "    [theta, phi; r > R]:  " << tmp << endl;
            if (tmp > CONST_EPS)
                cout << ss.str() << endl;
            break;
    }
    cout << endl << endl;
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

double Grid::Volume(long unsigned int i) const {
    return Volume(ind2sub1(i), ind2sub2(i), ind2sub3(i));
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

bool Grid::isempty(int id) const {
    unsigned long int k0, j0, i0;
    
    //Get search indices
    if (sys == CoordinateSystem::Cartesian)
            i0 = j0 = k0 = 1;
    else if (sys == CoordinateSystem::Cylindrical)
            k0 = 1;

    //Now loop through and add up values
    k0 = j0 = i0 = 0;
    for (unsigned long int i = i0; i < nx+i0; i++) {
        for (unsigned long int j = j0; j < ny+j0; j++) {
            for (unsigned long int k = k0; k < nz+k0; k++) {
                if (abs(m.at(id).at(sub2ind(i,j,k))) > CONST_EPS)
                    return false;
            }
        }
    }
    return true;
}

double Grid::sum(int id) const {
    unsigned long int k0, j0, i0;
    double out = 0;
    
    //Get search indices
    k0 = j0 = i0 = 0;
    if (sys == CoordinateSystem::Cartesian)
            i0 = j0 = k0 = 1;
    else if (sys == CoordinateSystem::Cylindrical)
            k0 = 1;

    //Now loop through and add up values    
    for (unsigned long int i = i0; i < nx+i0; i++) {
        for (unsigned long int j = j0; j < ny+j0; j++) {
            for (unsigned long int k = k0; k < nz+k0; k++) {
                out += m.at(id).at(sub2ind(i,j,k));
            }
        }
    }
    return out;
}

vec Grid::rand(long unsigned int i) const {
    return rand(ind2sub1(i), ind2sub2(i), ind2sub3(i));
}

vec Grid::rand(long unsigned int i, long unsigned int j, long unsigned int k) const {
    double eps, r, sth;
    vec out;
    
    switch (sys) {
        case (CoordinateSystem::Cartesian):
            out.X = (i-1 + roll())*Lx/nx - Lx/2;
            out.Y = (j-1 + roll())*Ly/ny - Ly/2;
            out.Z = (k-1 + roll())*Lz/nz;
            break;
        case (CoordinateSystem::Cylindrical):
            r = sqrt(roll() * (2*i+1) + i*i)*Lx/nx;
            sth = sin(2*CONST_PI*(j+roll())/ny);
            out.X = sqrt(1-sth*sth) * r;
            out.Y = sth * r;
            out.Z = (k-1 + roll())*Lz/nz;
            break;
        case (CoordinateSystem::Spherical):
            r = cbrt(roll() * (3*i*i+3*i+1) + i*i*i)*Lx/nx;
            eps = cos(CONST_PI/nz*(k+roll()));
            sth = sin(2*CONST_PI*(j+roll())/ny);
            out.X = r * sqrt(1-sth*sth) * eps;
            out.Y = r * sth * eps;
            out.Z = r * sqrt(1-eps*eps);
            break;
    }
    return out;
}

bool Grid::contains(const vec& x) const {
    switch (sys) {
       case (CoordinateSystem::Cartesian):
            if (abs(x.X) > Lx/2)
                return false;
            if (abs(x.Y) > Ly/2)
                return false;
            if (x.Z < 0 or x.Z > Lz)
                return false;
            return true;
        case (CoordinateSystem::Cylindrical):
            if (x.r2() > Lx*Lx)
                return false;
            if (x.Z < 0 or x.Z > Lz)
                return false;
            return true;
        case (CoordinateSystem::Spherical):
            if (x.r2() > Lx*Lx)
                return false;
            return true;
    }
    return true;
}

long unsigned int Grid::ncell() const {
    return nx * ny * nz;
}

Grid::iterator Grid::begin() const {
    return Grid::iterator(*this);
}

Grid::iterator& Grid::iterator::operator++(int) {
    //Get upper limits
    unsigned long int kmin = ((g.sys == CoordinateSystem::Spherical) ? 0 : 1);
    unsigned long int jmin = ((g.sys == CoordinateSystem::Cartesian) ? 1 : 0);
    unsigned long int lmin = ((g.sys == CoordinateSystem::Cartesian) ? 1 : 0);
    unsigned long int kmax = kmin+g.nz;
    unsigned long int jmax = jmin+g.ny;
    unsigned long int lmax = lmin+g.nx;
    
    //Increment x first
    ix ++;
    if (ix >= lmax) {
        ix = lmin;
        iy ++;
    }
    //If y was incremented
    if (iy >= jmax) {
        iy = jmin;
        iz ++;
    }
    //If z was incremented
    if (iz >= kmax) {
        ix = lmin;
        iy = jmin;
        iz = kmin;
        DONE = true;
    }
    return *this;
}

unsigned long int Grid::iterator::index() const {
    return g.sub2ind(ix, iy, iz);
}

Grid::iterator::iterator(const Grid& p) : g(p) {
    iz = ((g.sys == CoordinateSystem::Spherical) ? 0 : 1);
    iy = ((g.sys == CoordinateSystem::Cartesian) ? 1 : 0);
    ix = ((g.sys == CoordinateSystem::Cartesian) ? 1 : 0);
    DONE = false;
}

bool Grid::iterator::end() const {
    return DONE;
}

vector<double> Grid::cum(int id, double& sum) const {
    //Allocate some memory
    vector<double> out;
    out.push_back(sum);
    
    //Iterate through grid, create cumulant vector, and sum it up
    for (auto it = begin(); !it.end(); it ++) {
        sum += m.at(id).at(it.index());
        out.push_back(sum);
    }
    
    //Normalize by sum to get CDF
    for (long unsigned int i = 0; i < out.size(); i ++)
        out.at(i) /= sum;
    return out;
}

vector<double> Grid::cum(int id) const {
    //Allocate some memory
    double sum = 0;
    return cum(id, sum);
}

bool Grid::less(int id, double val) const {
    for (auto it = begin(); !it.end(); it ++ )
        if (m.at(id).at(it.index()) >= val)
            return false;
    return true;
}

#endif
