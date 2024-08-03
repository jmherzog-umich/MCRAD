#include <vector>
#include <cmath>
#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <utility>

#include "vec.h"
#include "utility.h"

#define CONST_PI            3.1415926535897932
#define CONST_C             299.792458
#define CONST_HBAR          1.054571817679489e-22
#define CONST_HK            47.992430733662212
#define CONST_APERY         1.2020569031595942
#define CONST_WIEN          1.5936242600400401
#define CONST_EPS           1e-10
#define CONST_PLANCKMAX     0.6476102378919149

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
  
    //Data  members
    vector<vector<double> > m;        //Data members, use a map to store multiple things
    vector<string> labels;            //Data member labels
    long unsigned int nx, ny, nz;     //Number of grid elements in each dimension
    double Lx, Ly, Lz;                //Dimensions of gridded region

    //Get first index of each dimension that corresponds to inside the grid
    virtual unsigned long int geti0() const = 0;
    virtual unsigned long int getj0() const = 0;
    virtual unsigned long int getk0() const = 0;
    virtual long unsigned int ntotal() const = 0;

    //Data accessors
    virtual long unsigned int ind1(const vec& x) const = 0;
    virtual long unsigned int ind2(const vec& x) const = 0;
    virtual long unsigned int ind3(const vec& x) const = 0;
    
    virtual long unsigned int ind2sub1(long unsigned int i) const = 0;
    virtual long unsigned int ind2sub2(long unsigned int i) const = 0;
    virtual long unsigned int ind2sub3(long unsigned int i) const = 0;
    
    virtual long unsigned int sub2ind(long unsigned int i, long unsigned int j, long unsigned int k) const = 0;

    //Geometry operations on cells
    virtual double Volume(long unsigned int i, long unsigned int j, long unsigned int k) const = 0;
    virtual vec rand(long unsigned int i, long unsigned int j, long unsigned int k) const = 0;
    virtual bool contains(const vec& x) const = 0;
    virtual double intersect(const vec& x, const vec& mu) const = 0;
    
    //Print grid specific info
    virtual void printOutside(ostream& oout, int id) const = 0;
    virtual void printGrid(ostream& oout) const = 0;
    virtual void printKeys(ostream& oout) const = 0;

    //Simple iterator alternative
    long unsigned int ncell() const;
    iterator begin() const;

    double& at(int id, long unsigned int i, long unsigned int j, long unsigned int k);
    double& at(int id, long unsigned int i);
    double& at(int id, const vec& x);
    double& at(const vec& x);
    
    double norm(int id, const vec& x) const;
    
    double Volume(long unsigned int i) const;
    double Volume(const vec& x) const;
    
    vec rand(long unsigned int i) const;
    
    bool isempty(int id) const;
    double sum(int id) const;
    bool less(int id, double val) const;
    
    vector<double> cum(int id, double& sum) const;
    vector<double> cum(int id) const;
    
    void newval(string label);
    void clear();
    
    void print(ostream& oout, bool ALL = false) const;
    
};

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

#endif
