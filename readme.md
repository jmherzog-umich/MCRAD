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

```
mcrad.o inputfile.in > outputfile.out
```

where the input file has a single parameter name followed by value or values
(separated by whitespace) per line, e.g.,:

```
spectrum-scattering 1 5 800
spectrum-absorption 0.003 5 800
spectrum-fluorescence 1 1 550 50
scattering-anisotropy 0 0.98
lifetime-fluorescence 1000
density 0.5
index-medium 1
index-front 1
index-side 1
index-back 1
grid-xy-points 20
grid-z-points 20
grid-theta-points 20
grid-t-points 500
grid-freq-points 100
simulation-packets 5000000
simulation-frequency-range 400 1000
simulation-length 100
simulation-radius 10000
simulation-timestep 0.1
simulation-timestep-fluorescence 20
enable-backwall 1
enable-frontwall 1
enable-sidewall 1
enable-cartesian 1
enable-periodic-xy 0
enable-fluorescence 1
enable-time-dependent 1
enable-single-photon 1
beam-photons-per-packet 10000
beam-radius-cutoff 5000
beam-profile 1
database-filename test.csv
export-paths 100 test
execute
```

Parameters can be specified on the command line as well, but multiple arguments
for a single parameter must be comma-separated, e.g.,

```
mcrad.o inputfile.in simulation-length 300 spectrum-absorption 0.005,5,800 > outputfile.out
```

The parameter/command definitions are as follows:

| Parameter | Arguments | Description |
| --- | --- | --- |
| spectrum-scattering | value [model] [parameters] | Scattering cross-section distribution (square microns) |
| spectrum-absorption | value [model] [parameters] | Absorption cross-section distribution (square microns) |
| spectrum-fluorescence | quantum_yield [model] [parameters] | Fluorescence spectrum, FQY, and parameters|
| scattering-anisotropy | model value [parameters] | Scattering anisotropy model, value g, and any additional parameters |
| density | value | Particle density (per cubic micron) |
| lifetime-fluorescence | value | Single-exponential fluorescence lifetime |
| index-medium | value | Bulk medium refractive index |
| index-front | value | Front surface refractive index |
| index-back | value | Rear surface refractive index |
| index-side | value | Side surface refractive index |
| grid-xy-points | value | Number of x-y/r grid points |
| grid-z-points | value | Number of z grid points |
| grid-theta-points | value | Number of incidence/exit angle grid points |
| grid-t-points | value | Number of temporal grid points |
| grid-freq-points | value | Number of radiation frequency grid points |
| simulation-packets | value | Total number of photon packets |
| simulation-radius | value | Radial extent (or half-width) of simulation domain (micron) |
| simulation-length | value | Depth of simulation domain (micron) |
| simulation-timestep | value | Temporal resolution of simulation (picoseconds) |
| simulation-timestep-fluorescence | value | Temporal resolution of fluorescence emission and grid (picoseconds) |
| simulation-frequency-range | min max | Frequency range for gridding spectrally-resolved results, and generating fluorescence |
| roulette-minweight | value | Packet weight at which to begin roulette procedure |
| roulette-newweight | value | Minimum packet weight if packet survives roulette procedure |
| max-steps | value | Max number of steps before abandoning calculation |
| moments | value | Stores and calculates moments up to this order (integer, 1-4) |
| enable-backwall | value | Enable back wall (boolean 1 or 0) |
| enable-frontwall | value | Enable front wall (boolean 1 or 0) |
| enable-sidewall | value | Enable side wall (boolean 1 or 0) |
| enable-cartesian | value | Enable Cartesian wall geometry instead of cylindrical (boolean 1 or 0) |
| enable-periodic-xy | value | Enable periodic boundary conditions in R/XY direction (boolean 1 or 0) |
| enable-periodic-z | value | NOT YET IMPLEMENTED: Enables periodic boundary conditions in depth direction (boolean 1 or 0) |
| enable-fluorescence | value | Enable fluorescence calculation (boolean 1 or 0) |
| enable-fluorescence-trapping | value | Enable fluorescence calculation (boolean 1 or 0) |
| enable-time-dependent | value | Enable time-resolved simulation (boolean 1 or 0) |
| enable-single-photon | value | Enable single-photon mode (boolean 1 or 0) |
| enable-interference | value | NOT YET IMPLEMENTED: account for interference effects in images (boolean 1 or 0) |
| enable-saturation | value | TESTING: account for saturation of absorbers based on grid resolution (boolean 1 or 0) |
| beam-profile | radius [model] [parameter] | Source beam radius (micron), profile, and additional parameters |
| beam-divergence | width [model] [parameter] | Beam spread width (rad), model, and additional parameters |
| beam-pulse-duration | duration [model] | Beam temporal width (ps) and profile |
| beam-spectrum | [model] [parameter] | Beam peak frequency, spectrum model, and model parameters |
| beam-radius-cutoff | value | Maximum radius at which photons will be generated (micron) |
| beam-photons-per-packet | value | Initial number of photons per packet/scaling factor (photons) |
| beam-initial-angle | value | Sin of initial incidence angle of beam |
| beam-focal-depth | value | Focal depth of beam (micron) |
| track-photon | value | DEBUG: prints to stderr the status of the photon stored at this index at each step (integer) |
| print-steps | value | Prints the grid data to the output file every 'printsteps' steps (integer) [Useful for time- or order-resolved data] |
| database-filename | name | Output filename (CSV) to store scalar results (string) |
| export-paths | number filename | Number of paths (X-Y-Z-T values) to store to disk (integer) and basename of file |
|  |  |  |
| run | | Runs simulation immediately |
| setup | | Runs setup calculations based on currently set inputs |
| generate-beam | | Randomly samples the beam settings and generates initial photons |
| print | | Prints statistical results to stdout |
| write | | Writes any output paths, database files, etc., based on settings |
| execute | | Reads command line input, then sequentially performs setup, photon generation, simulation, printing, and writing tasks |

Models are specified by an integer value corresponding to the internal enum value,
and parameters are specific to each model. Currently, non-constant parameters are
not fully documented so please refer to the code for details on how to specify
multi-parameter models.

## Future work
Several features are planned as indicated in the previous section:
- More and better models for cross-sections and spectra
- Improved scattering anisotropy models
- Image simulation
- Parallelization and automatic uncertainty analysis
