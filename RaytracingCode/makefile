#Makefile

F90G = gfortran
OPTG = -O2 -fimplicit-none -fdefault-integer-8 -fdefault-double-8 -fdefault-real-8
OPTDG = -g -Wall -Wtabs -Wunderflow -fimplicit-none -Wconversion -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -Wuninitialized -Wextra -Wredundant-decls
MAING = raytracing_5.0_gfortran

F90I = ifort
OPTI = -implicitnone -i8 -r8 -prec-div -prec-sqrt -fast -O5
OPTDI= -implicitnone -i8 -r8 -warn all -check all -fPIC -fpe0 -g -debug extended -traceback
MAINI = raytracing_5.0_ifortran

#-O2 -fimplicit-none -fdefault-integer-8 -fdefault-double-8 -fdefault-real-8
# -g :debug
# -O2 -fimplicit-none -fdefault-integer-8 -fdefault-double-8 -fdefault-real-8 :compilation optimale 
# OPT (debug) = -g -Wall -Wtabs -Wunderflow -fimplicit-none -Wconversion -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -Wuninitialized -Wextra -Wredundant-decls


allG: moduleG programG clean

debugG: moduleDG programDG clean

allI: moduleI programI clean

debugI: moduleDI programDI clean

moduleDG: #*.a	

	$(F90G) -c $(OPTDG) constantes.f90
	$(F90G) -c $(OPTDG) coord.f90	
	$(F90G) -c $(OPTDG) environ.f90
	$(F90G) -c $(OPTDG) subs_raytracing.f90
	$(F90G) -c $(OPTDG) RungeKuttaFehlberg.f90
	$(F90G) -c $(OPTDG) RungeKutta.f90
	$(F90G) -c $(OPTDG) params_raytracing.f90

#	ar -r biblio.a *.o

moduleDI: #*.a	

	$(F90I) -c $(OPTDI) constantes.f90
	$(F90I) -c $(OPTDI) coord.f90	
	$(F90I) -c $(OPTDI) environ.f90
	$(F90I) -c $(OPTDI) subs_raytracing.f90
	$(F90I) -c $(OPTDI) RungeKuttaFehlberg.f90
	$(F90I) -c $(OPTDI) RungeKutta.f90
	$(F90I) -c $(OPTDI) params_raytracing.f90

#	ar -r biblio.a *.o


moduleG: #*.a	

	$(F90G) -c $(OPTI) constantes.f90
	$(F90G) -c $(OPTI) coord.f90	
	$(F90G) -c $(OPTI) environ.f90
	$(F90G) -c $(OPTI) subs_raytracing.f90
	$(F90G) -c $(OPTI) RungeKuttaFehlberg.f90
	$(F90G) -c $(OPTI) RungeKutta.f90
	$(F90G) -c $(OPTI) params_raytracing.f90

#	ar -r biblio.a *.o

moduleI: #*.a	

	$(F90I) -c $(OPTI) constantes.f90
	$(F90I) -c $(OPTI) coord.f90	
	$(F90I) -c $(OPTI) environ.f90
	$(F90I) -c $(OPTI) subs_raytracing.f90
	$(F90I) -c $(OPTI) RungeKuttaFehlberg.f90
	$(F90I) -c $(OPTI) RungeKutta.f90
	$(F90I) -c $(OPTI) params_raytracing.f90

#	ar -r biblio.a *.o

programG:
	$(F90G) -c $(OPTG) $(MAING).f90
	$(F90G) -o raytracing *.o #biblio.a

programDG:
	$(F90G) -c $(OPTDG) $(MAING).f90 
	$(F90G) -o raytracing *.o #biblio.a
	
programI:
	$(F90I) -c $(OPTI) $(MAINI).f90
	$(F90I) -o raytracing *.o #biblio.a

programDI:
	$(F90I) -c $(OPTDI) $(MAINI).f90 
	$(F90I) -o raytracing *.o #biblio.a	

# -L chemin -lbiblio  si la bibliotheque n'est pas ds le meme dossier
clean:
	rm -f *~
	rm -f *.o
	rm -f *.mod


