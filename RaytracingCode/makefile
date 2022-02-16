#Makefile

F90 = gfortran
OPT = -O2 -fimplicit-none -fdefault-integer-8 -fdefault-double-8 -fdefault-real-8 -g -Wall
OPTD = -g -Wall -Wtabs -Wunderflow -fimplicit-none -Wconversion -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -Wuninitialized -Wextra -Wredundant-decls
#-O2 : compilation optimale
#-g :debug

##ifort
#OPT = -implicitnone -i8 -r8 -prec-div -prec-sqrt -fast -O5
#OPTD= -implicitnone -i8 -r8 -warn all -check all -fPIC -fpe0 -g -debug extended -traceback

CDF_DIR = /usr
INC_CDF = $(CDF_DIR)/include/
LIB_CDF = $(CDF_DIR)/lib/x86_64-linux-gnu/

all: module program clean

debug: moduleD programD clean

moduleD: #*.a	

	$(F90) -c $(OPTD) constantes.f90
	$(F90) -c $(OPTD) coord.f90	
	$(F90) -c $(OPTD) environ.f90
	$(F90) -c $(OPTD) subs_raytracing.f90
	$(F90) -c $(OPTD) RungeKuttaFehlberg.f90
	$(F90) -c $(OPTD) RungeKutta.f90
	$(F90) -c $(OPTD) params_raytracing.f90

	ar -r biblio.a *.o


module: #*.a	

	$(F90) -c $(OPT) constantes.f90 
	$(F90) -c $(OPT) coord.f90 
	$(F90) -c $(OPT) variables.f90 
	$(F90) -c $(OPT) defs_basic_cdf.f90 -I$(INC_CDF) 
	$(F90) -c $(OPT) environ.f90 -I$(INC_CDF) 
	$(F90) -c $(OPT) subs_raytracing.f90 
	$(F90) -c $(OPT) RungeKuttaFehlberg.f90 
	$(F90) -c $(OPT) RungeKutta.f90 
	$(F90) -c $(OPT) params_raytracing.f90 

	ar -r biblio.a *.o 

program:
	$(F90) -c $(OPT) raytracing_5.f90 -I$(INC_CDF) 
	$(F90) -o ARTEMIS-P raytracing_5.o biblio.a -I$(INC_CDF) -L$(LIB_CDF) -lnetcdf -lnetcdff
#*.o 
programD:
	$(F90) -c $(OPTD) raytracing_5.f90
	$(F90) -o ARTEMIS-P raytracing_5.o biblio.a
#*.o 

# -L chemin -lbiblio  si la bibliotheque n'est pas ds le meme dossier
clean:
	rm -f *~
	rm -f *.o
	rm -f *.mod
