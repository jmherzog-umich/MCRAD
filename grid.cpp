#include <vector>
#include <cmath>
#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <utility>

#include "vec.h"
#include "utility.h"

#include "grid.h"

#define CONST_PI            3.1415926535897932
#define CONST_C             299.792458
#define CONST_HBAR          1.054571817679489e-22
#define CONST_HK            47.992430733662212
#define CONST_APERY         1.2020569031595942
#define CONST_WIEN          1.5936242600400401
#define CONST_EPS           1e-10
#define CONST_PLANCKMAX     0.6476102378919149

using namespace std;

void Grid::newval(string label) {
    m.push_back(vector<double>());
    labels.push_back(label);
}

void Grid::clear() {
    int N = ntotal();
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

void Grid::print(ostream& oout, bool ALL) const {
    //Some initial parameters for the loops
    unsigned long int k0, j0, l0;
    double val = 0;
    k0 = getk0();
    j0 = geti0();
    l0 = getj0();
    
    //Now loop through and print the grid
    for (unsigned long int i = 0; i < m.size(); i ++) {
        //Don't print if the grid is empty
        if (not ALL and isempty(i))
            continue;
    
        //print Key
        printKeys(oout);
    
        //Print the key name
        val = 0;
        for (unsigned long int z = 0; z < m.at(i).size(); z ++)
            val += m.at(i).at(z);
        oout << "   " << labels.at(i) << "   SUM=" << val << "   (" << sum(i) << ")" << endl;
        
        //Now print the data
        for (unsigned long int k = k0; k < (k0+nz); k ++) {                //Blocks    
            for (unsigned long int l = l0; l < (l0+ny); l ++) {            //Rows
                for (unsigned long int j = j0; j < (j0+nx); j ++)          //Columns
                    oout << scientific << setw(18) << m.at(i).at(sub2ind(j, l, k));
            
            //Skip line after each row
            oout << endl;
            }
        
        //Skip an extra line after each block ONLY IF WE HAVE MULTIPLE LINES PER BLOCK
        if (ny > 1)
            oout << endl;
        }

    //Now print the elements outside the grid
    printOutside(oout, i);
    
    //Skip an extra line for fun
    oout << endl;   
    }
}

double Grid::dens(int id, const vec& x) const {
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

bool Grid::isempty(int id) const {
    //Get starting indices for grid
    unsigned long int k0, j0, i0;
    i0 = geti0();
    j0 = getj0();
    k0 = getk0();

    //Now loop through and add up values
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
    //Get starting indices for grid
    unsigned long int k0, j0, i0;
    i0 = geti0();
    j0 = getj0();
    k0 = getk0();
    double out = 0;

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

long unsigned int Grid::ncell() const {
    return nx * ny * nz;
}

Grid::iterator Grid::begin() const {
    return Grid::iterator(*this);
}

Grid::iterator& Grid::iterator::operator++(int) {
    //Get upper limits
    unsigned long int kmin = g.getk0();
    unsigned long int jmin = g.getj0();
    unsigned long int lmin = g.geti0();
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
    iz = g.getk0();
    iy = g.getj0();
    ix = g.geti0();
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

void Grid::collideCell(const vec& x, const vec& mu, double& ds, int& reflect) const {
    double ds2 = intersect(x, mu);
    if (ds2 <= ds and ds > 0) {
        ds = ds2;
        reflect = 0;
    }
}

void Grid::collideFront(const vec& x, const vec& mu, double& ds, int& reflect) const {
    if (mu.Z >= 0)
        return;
        
    double ds2 = -x.Z/mu.Z;
    if (ds2 <= ds and ds2 >= 0) {
        ds = ds2;
        reflect = 2;
    }
}
    
void Grid::collideBack(const vec& x, const vec& mu, double& ds, int& reflect) const {
    if (mu.Z <= 0)
        return;
        
    double ds2 = (Lz - x.Z)/mu.Z;
    if (ds2 <= ds and ds2 >= 0) {
        ds = ds2;
        reflect = 1;
    }
}

void Grid::collideSide(const vec& x, const vec& mu, double& ds, int& reflect) const {
    //Allocate variables
    double sx, sy;
    
    //X direction
    if (abs(mu.X) <= 0)
        sx = -1;
    else if (mu.X > 0)
        sx = (Lx/2.0-x.X)/mu.X;
    else
        sx = (-Lx/2.0-x.X)/mu.X;
    
    //Y direction
    if (abs(mu.Y) < 0)
        sy = -1;
    else if (mu.Y > 0)
        sy = (Ly/2.0-x.Y)/mu.Y;
    else
        sy = (-Ly/2.0-x.Y)/mu.Y;
        
    //First check x
    if (sx < ds and sx >= 0) {
        reflect = 3;
        ds = sx;
    }
    
    //Now check y
    if (sy < ds and sy >= 0) {
        reflect = 4;
        ds = sy;
    }
        
}

void Grid::boundXY(vec& x) const {
    while (abs(x.X)>Lx)
        x.X -= copysign(2*Lx, x.X);
    while (abs(x.Y)>Ly)
        x.Y -= copysign(2*Ly, x.Y);
}

void Grid::boundZ(vec& x) const {
    while ( x.Z > Lz )
        x.Z -= Lz;
    while ( x.Z < 0 )
        x.Z += Lz;
}

