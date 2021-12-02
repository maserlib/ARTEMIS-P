MODULE constantes
	
	implicit none
	
    real(kind=8), parameter :: pi = acos(-1.d0)
	real(kind=8), parameter :: c = 	3.d8           !m/s
	real(kind=8), parameter :: me = 9.10000d-31    !kg
	real(kind=8), parameter :: e  = -1.60000d-19   !C
	real(kind=8), parameter :: eps0 = 8.8542d-12   !F/m
	real(kind=8), parameter :: RT = 6378.0d0       !km
	real(kind=8), parameter :: RS = 60300.0d0      !km
	real(kind=8), parameter :: RJ = 69911.0d0      !km
	real(kind=8), parameter :: RSol = 6.960d5      !km
	complex(kind=8), parameter :: ic = cmplx(0.0,1.0) 

        real(kind=8) :: z0,n0,B0,z0B
END MODULE constantes
