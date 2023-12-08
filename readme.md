### MCRAD

MCRAD is a simple <emph>Monte Carlo</emph> photon transport code with a focus
on optical diagnostics and fluorescence analysis of dense slabs. Internal units
are picoseconds, terahertz, and micrometers.

## Compilation
The code requires no external dependencies and can in principle be compiled on
any platform that has a C++-17 compliant compiler and standard library. MCRAD
can be compiled with g++ using the following command

```
g++ -std=c++17 -static -Ofast -o mcrad.o mcrad.cpp photon.h stats.h vec.h
 simulation.h utility.h medium.h beam.h camera.h image.h raypath.h -Wall
```
 
MCRAD was compiled and tested on Windows 11 and Fedora 38 with GCC 13 and should
compile with no warnings or errors with the specified settings.

## Running the code
MCRAD accepts an input file and command line arguments, prints results to
stdout, and prints progress and errors to stderr. Typical usage is

mcrad.o inputfile.in > outputfile.out

where the input file has a single parameter name followed by value or values
(separated by whitespace) per line, e.g.,:

```
Ss 1
Sa 0.003
xca 0
xcs 0
N0 5000000
n0 1
nx 1
n 1
nr 1
Rres 20
Zres 20
THres 20
Tres 200
dens 0.0005
R 10000
L 10000
dT 0.01
Rb 1
Rbmax 1
phase 0
backwall 1
frontwall 1
sidewall 1
cartesian 1
periodicxy 0
fluorescence 1
timedependent 0
singlephoton 1
beamprofile 0
beamspread 1
dbfile test.csv
exportpaths 100
pathbasefilename test
exec
```

The parameter/command definitions are as follows:

| Parameter | Description |
| --- | --- |
| Ss | Scattering cross-section (square microns) |
| Sa | Absorption cross-section (square microns) |
| g | Scattering anisotropy parameter(s) (floating point) |
| FQY | Particle fluorescence quantum yield (floating point) |
| xca | Absorption cross-section model (enum value) |
| xcs | Scattering cross-section model (enum value) |
| xce | NOT YET IMPLEMENTED: Fluorescence spectrum model (enum value)
| femit | Emitted photon frequency parameter(s) (floating point) |
| dens | Particle density (per cubic micron) |
| phase | Scattering phase-function (enum value) |
|  |  |
| N0 | Initial number of photon packets |
| n0 | Front surface refractive index |
| nx | Rear surface refractive index |
| nr | Side surface refractive index |
| n | Bulk medium refractive index |
| R | Radial extent (or half-width) of simulation domain (micron) |
| L | Depth of simulation domain (micron) |
| dT | Temporal resolution of simulation (picoseconds) |
| Rres | Number of grid points in radial/XY directions |
| Zres | Number of grid points in depth direction |
| THres | Number of grid points in angular direction |
| Tres | Number of grid points in temporal direction |
|  |  |
| Wmin | Packet weight at which to begin roulette procedure |
| Wm | Minimum packet weight if packet survives roulette procedure |
| maxstep | Max number of steps before abandoning calculation |
| moments | Stores and calculates moments up to this order (integer, 1-4) |
|  |  |
| backwall | Enable back wall (boolean 1 or 0) |
| frontwall | Enable front wall (boolean 1 or 0) |
| sidewall | Enable side wall (boolean 1 or 0) |
| cartesian | Enable Cartesian wall geometry instead of cylindrical (boolean 1 or 0) |
| periodicxy | Enable periodic boundary conditions in R/XY direction (boolean 1 or 0) |
| periodicz | NOT YET IMPLEMENTED: Enables periodic boundary conditions in depth direction (boolean 1 or 0) |
| fluorescence | Enable fluorescence calculation (boolean 1 or 0) |
| timedependent | Enable time-resolved simulation (boolean 1 or 0) |
| singlephoton | Enable single-photon mode (boolean 1 or 0) |
| interference | NOT YET IMPLEMENTED: account for interference effects in images (boolean 1 or 0) |
| saturation | TESTING: account for saturation of absorbers based on grid resolution (boolean 1 or 0) |
|  |  |
| Rb | Source beam radius (micron) |
| Rbmax | Maximum radius at which photons will be generated (micron) |
| beamprofile | Beam profile shape (enum value) |
| beamspread | Beam spread model (enum value) |
| beamwidth | Beam temporal profile (enum value) |
| beamspec | Beam spectral profile (enum value) |
| E | Initial number of photons per packet/scaling factor (floating point) |
| sin0 | Sin of initial incidence angle of beam (floating point) |
| Pb | Angular width of beam (rad) |
| Sb | Semi-minor axis width or inner radius of beam (micron) |
| Zb | Focal depth of beam (micron) |
| Tb | Temporal width of beam (picoseconds) |
| wb | Peak frequency of beam (THz) |
| dwb | Frequency width of beam (THz) |
|  |  |
| trackphoton | DEBUG: prints to stderr the status of the photon stored at this index at each step (integer) |
| printsteps | Prints the grid data to the output file every 'printsteps' steps (integer) [Useful for time- or order-resolved data] |
| dbfile | Output filename (CSV) to store scalar results (string) |
| exportpaths | Number of paths (X-Y-Z-T values) to store to disk (integer) |
| pathbasefilename | Base filename of path.xyz file (string) |
|  |  |
| run | Runs simulation immediately |
| setup | Runs setup calculations based on currently set inputs |
| genBeam | Randomly samples the beam settings and generates initial photons |
| print | Prints statistical results to stdout |
| write | Writes any output paths, database files, etc., based on settings |
| exec | Reads command line input, then calls 'setup', 'genBeam', 'run', 'print', and 'write', sequentially |

Several parameters can take on multiple values. For example, the absorption
cross-section can be specified as parameters to a function depending on the
selected model. In this case, the parameters are specified as a sequence of values
separated by white space, e.g.,

```
Sa 0 700 0.003
xca 12
```

which uses a step-function model with a cutoff frequency of 700 THz such that the
absorption cross-section is 0.003 square microns for frequencies larger than 700 THz
and 0 for frequencies below 700 THz.

Currently, non-constant parameters are not fully supported so please refer to the
code for details on how to specify multi-parameter models.

## Future work
Several features are planned as indicated in the previous section:
- Better support for polychromatic simulations
- More and better models for cross-sections and spectra
- Image simulation
- Parallelization and automatic uncertainty analysis
