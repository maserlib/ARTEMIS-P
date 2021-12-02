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
* _init_raytracing_ : constains all parameters for the initialization of the program. 

### Defining Density and Magnetic Field
Currently, the density and magnetic field are defined by an analytical function only. In the file _environ.f90_, the routine _magneticf_ calculates the magnetic field and _density_ computes the electron density.

### Outputs
the program produces for each ray a file containing various data calculated at each time step and along the trajectory of the ray.

### Future improvements
It will be possible to load a complete environment from a simulation for instance.  
In addition, the format of the results will also be updated.
