#Makefile

F90 = ifort
OPT = -implicitnone -i8 -r8 -prec-div -prec-sqrt -fast -O5
OPTD= -implicitnone -i8 -r8 -warn all -check all -fPIC -fpe0 -g -debug extended -traceback
#-O2 -fimplicit-none -fdefault-integer-8 -fdefault-double-8 -fdefault-real-8
# -g :debug
# -O2 -fimplicit-none -fdefault-integer-8 -fdefault-double-8 -fdefault-real-8 :compilation optimale 
# OPT (debug) = -g -Wall -Wtabs -Wunderflow -fimplicit-none -Wconversion -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -Wuninitialized -Wextra -Wredundant-decls


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

#	ar -r biblio.a *.o


module: #*.a	

	$(F90) -c $(OPT) constantes.f90
	$(F90) -c $(OPT) coord.f90	
	$(F90) -c $(OPT) environ.f90
	$(F90) -c $(OPT) subs_raytracing.f90
	$(F90) -c $(OPT) RungeKuttaFehlberg.f90
	$(F90) -c $(OPT) RungeKutta.f90
	$(F90) -c $(OPT) params_raytracing.f90

#	ar -r biblio.a *.o

program:
	$(F90) -c $(OPT) raytracing_5.0.f90
	$(F90) -o raytracing raytracing_5.0.o *.o #biblio.a

programD:
	$(F90) -c $(OPTD) raytracing_5.0.f90
	$(F90) -o raytracing raytracing_5.0.o *.o #biblio.a

# -L chemin -lbiblio  si la bibliotheque n'est pas ds le meme dossier
clean:
	rm -f *~
	rm -f *.o
	rm -f *.mod


