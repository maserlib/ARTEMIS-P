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

	B(1,:)=0. 
	B(2,:)=0.
	B(3,:)=B0*(V(3,:)/z0B)**(-3)

end subroutine magneticf


subroutine read_environ()
        !----------------------------------------------------
        ! Lit le fichier init_environ.txt et range les donnees
        !  dans un fichier
        !----------------------------------------------------
        use constantes
        use netcdf

        open(40,file='init_environ.txt')
        read(40,*)n0
        read(40,*)z0
        read(40,*)B0
        read(40,*)z0B
        close(40)
        write (*,*)"n0 ",n0," z0 ",z0," B0 ",B0," z0B ",z0B


end subroutine read_environ

subroutine density(V,Ne)
	!-----------------------------------------------------
	! Electron density computation
	! CALLING SEQUENCE: call density(V,Ne)
	! INPUTS:   V: position vector in cartesian coordinates
	! OUTPUTS: Ne: electron density (cm-3)
	!-----------------------------------------------------
        use constantes
        real(kind=8),dimension(:,:),intent(in)    :: V
	real(kind=8), dimension(:), intent(out)   :: Ne



!	real(kind=8)                              :: n0,z0

!	open(40, file='init_environ.txt')
!	read(40,*) n0
!	read(40,*) z0

	Ne(:)=n0*(V(3,:)/z0)**(-2)
	
!	close(40)



end subroutine density

subroutine read_environ_netcdf()
  use netcdf
  integer(kind=4) :: ncid,stId,varid
  real(kind=8)::dens(3)

  stId = nf90_open("prod_jia_o2p_330.nc",nf90_nowrite,ncid)
  stId = nf90_inq_varid(ncid,"Prod",varid)
  stId = nf90_get_var(ncid,varid,dens)

  write(*,*)"dens = ",dens

  stId = nf90_close(ncid)
end subroutine read_environ_netcdf


END MODULE environ
