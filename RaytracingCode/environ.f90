MODULE environ
	
	use constantes
	implicit none
	
	contains

subroutine magneticf(V,B)
	!----------------------------------------------------------
	! Magnetic field computation.
	! CALLING SEQUENCE: call magneticf(V,B)
	! INPUTS:  V: position vector in cartesian coordinates
	! OUTPUTS: B: magnetic field vector in cartesian coordinate (Tesla)
	!----------------------------------------------------------
	real(kind=8), dimension(:,:),intent(in)  :: V
	real(kind=8), dimension(:,:),intent(out) :: B
	real(kind=8)							 :: B0,z0

	open(40, file='init_environ.txt')
	read(40,*) B0
	read(40,*) z0

	B(1,:)=0. 
	B(2,:)=0.!B0*1.d-4*(V(2,:)/z0)**(-3)
	B(3,:)=0.
	close(40)
end subroutine magneticf

subroutine density(V,Ne)
	!-----------------------------------------------------
	! Electron density computation
	! CALLING SEQUENCE: call density(V,Ne)
	! INPUTS:   V: position vector in cartesian coordinates
	! OUTPUTS: Ne: electron density (cm-3)
	!-----------------------------------------------------
	real(kind=8),dimension(:,:),intent(in)    :: V
	real(kind=8), dimension(:), intent(out)   :: Ne
	real(kind=8)                              :: ne1,ne2,x0,a

	open(40, file='init_environ.txt')
	read(40,*) 
	read(40,*) 
	read(40,*) ne1
	read(40,*) ne2
	read(40,*) x0
	read(40,*) a
	

	Ne(:)=ne1+(ne2-ne1)/2.d0*(1+tanh((V(1,:)-x0)/a))
	
	close(40)
end subroutine density

END MODULE environ
