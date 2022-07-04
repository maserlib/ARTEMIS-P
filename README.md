# ARTEMIS-P
Anisotropic Ray Tracer for Electromagnetism in Magnetospheres, Ionospheres and Solar wind, including Polarisation

A short description of the code is available in the [Proceedings of the 2013 URSI-France meeting](http://ursi-france.telecom-paristech.fr/fileadmin/journees_scient/docs_journees_2013/data/articles/000052.pdf).


## Requirement
### Compilation
The default compiler is [GNU Fortran Compiler](https://gcc.gnu.org/wiki/GFortranBinaries).  
It is also possible to compile with [Intel Fortran Compiler](https://www.intel.com/content/www/us/en/developer/tools/oneapi/fortran-compiler.html#gs.i844d1).

## How to run
### input files
There are 2 initialization files : 
* _init_environ.txt_ : contains initial values of the environement.
* _init_raytracing.txt_ : constains all parameters for the initialization of the program. 

### Defining Density and Magnetic Field
Currently, the density and magnetic field are defined by an analytical function only. In the file _environ.f90_, the routine _magneticf_ calculates the magnetic field and _density_ computes the electron density.

### Outputs
the program produces for each ray a file containing various data calculated at each time step and along the trajectory of the ray.

### Future improvements
It will be possible to load a complete environment from a simulation for instance.  
In addition, the format of the results will also be updated.


## Description of input files

You need to specify all the lines. You can find two simple examples on Toy1 and Toy2 branches.

### _init_environ.txt_ :
n0 in cm^-3 ; used in computation of Ne (see the routine _density_ in _environ.f90_)
z0 in km ; used in the computation of Ne (see the routine _density_ in _environ.f90_)
B0 in T ; used in computation of B (see the routine _magneticf_ in _environ.f90_)
z0B in km ; used in computation of B (see the routine _magneticf_ in _environ.f90_)

### _init_raytracing.txt_ :
<directory_name> : maximum length = 4
<integrator_name> : choose between RK4 and RK4F
frequency in kHz (real)
number of rays (integer)
number of iterations (integer)
coordinate system. Choose between cart or sphr
length units. Choose between km, m or rs
initial position of rays (x,y,z). Give the initial position of each ray on a different line.
k angle (kx, ky, kz). Give the initial position of each ray on a different line.
propagation mode. "X" for extraordinary mode or "O" for ordinary mode. Give the initial position of each ray on a different line.
maximum value of refractive index (integer)
initial time step (real)
maximum time step (real)
minimum time step (real)
