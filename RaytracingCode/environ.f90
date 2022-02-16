MODULE environ
  
  use constantes
  use variables
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
 ! B(3,:)=B0*(V(3,:)/z0B)**(-3)
 B(3,:)=0. !sqrt(2.0)*79e-9 ! Tesla from lathys env_ganymede

end subroutine magneticf

subroutine read_magnetic_field_netcdf()


end subroutine read_magnetic_field_netcdf


subroutine read_environ()
  !----------------------------------------------------
  ! Lit le fichier init_environ.txt et range les donnees
  !  dans un fichier
  !----------------------------------------------------

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
  real(kind=8),dimension(:,:),intent(in)    :: V
  real(kind=8), dimension(:), intent(out)   :: Ne
  real(kind=8),dimension(:),allocatable :: r,th,ph
  real(kind=8) :: y
  integer :: it,l,ii,jj,kk, N

!  do ii=1,size_r
!    do jj=1,size_th
!      do kk=1,size_ph
!        write(*,*)"dens : ",ii,jj,kk,Dn(ii,jj,kk)
!      enddo
!    end do
!  enddo
  !Ne(:)=50.0!n0*(V(3,:)/z0)**(-2)
  N = size(Ne)
  allocate(r(N),th(N),ph(N))
  r(:)=0 ; th(:)=0 ; ph(:) = 0


!  do ii = 1,Ni
  !  if(V(2,ii)<0)then 
  !    y = 0
  !  else 
  !    y = V(2,ii)
  !  endif
  !i  if(V(1,ii) .gt. 0)then
    where(V(1,:)==0)
      ph(:)=pi/2
    elsewhere
      ph(:) = mod(atan(V(2,:)/V(1,:))+4*pi,2*pi)
    end where
  !  else if (V(1,ii) .lt. 0) then 
  !    ph(:) = mod(atan(y/V(1,:)) + 5*pi,2*pi)
  !  else
  !    ph(:)=pi/2
  !  endif
 ! enddo
  r(:)  = sqrt(V(1,:)*V(1,:) + V(2,:)*V(2,:) + V(3,:)*V(3,:))
  th(:) = mod(acos(V(3,:)/r(:))+4*pi, 2*pi)

  write(*,*)"ici r_ ",r(:)
  write(*,*)"ici th ",th(:)
  write(*,*)"ici ph ",ph(:)
  do it=1,N
    ii=0; jj=0; kk=0
    do l=1,size_r
      if(r_low(l) <= r(it) .and. r(it) < r_upp(l))&
         ii = l
    enddo
    if(ii==0) write(*,*)"ii=0",r_low(l),r(it),r_upp(l)
    do l=1,size_th
      if(th_low(l) <= th(it) .and. th(it) < th_upp(l))&
        jj = l
    enddo
    if(jj==0) write(*,*)"jj=0",th_low(l),th(it),th_upp(l)

    do l=1,size_ph
      if(ph_low(l) <= ph(it) .and. ph(it) < ph_upp(l))&
         kk = l
    enddo
    if(kk==0) write(*,*)"kk=0",ph_low(l),ph(it),ph_upp(l)
    Ne(it) = Dn(ii,jj,kk)
   write(*,*)"Ne = ",Ne(it)
  enddo

  deallocate(r,th,ph)

end subroutine density

subroutine read_environ_netcdf()
  use netcdf
  use defs_basic_cdf
  integer(kind=4) :: ncid,stId
!  real(kind=8),allocatable,dimension(:,:,:) ::dens
!  real(kind=8),allocatable,dimension(:) :: r_low,r_upp,th_low,&
!    th_upp,ph_low,ph_upp
!  real(kind=8) :: s_min(3),s_max(3),gstep(3)
!  real(kind=8) :: r_planet,s_centr
!  integer(kind=4) :: size_ph,size_r,size_th

  stId = nf90_open("../prod_jia_o2p_330.nc",nf90_nowrite,ncid)
  call test_cdf(stId)
  call get_simple_dimens_cdf(ncid,"size_r",size_r)
  call get_simple_dimens_cdf(ncid,"size_ph",size_ph)
  call get_simple_dimens_cdf(ncid,"size_th",size_th)
  call get_simple_variable_cdf(ncid,"s_min",s_min)
  call get_simple_variable_cdf(ncid,"s_max",s_max)
  call get_simple_variable_cdf(ncid,"gstep",gstep)
  call get_simple_variable_cdf(ncid,"phys_length",phys_length)
  call get_simple_variable_cdf(ncid,"r_planet",r_planet)
  call get_simple_variable_cdf(ncid,"s_centr",s_centr)

  allocate(Dn(size_r,size_th,size_ph))
  allocate(r_low(size_r),r_upp(size_r))
  allocate(th_low(size_th),th_upp(size_th))
  allocate(ph_low(size_ph),ph_upp(size_ph))

  call get_simple_variable_cdf(ncid,"Prod",Dn)
  call get_simple_variable_cdf(ncid,"r_low",r_low)
  call get_simple_variable_cdf(ncid,"r_upp",r_upp)
  call get_simple_variable_cdf(ncid,"elev_low",th_low)
  call get_simple_variable_cdf(ncid,"elev_upp",th_upp)
  call get_simple_variable_cdf(ncid,"azimut_low",ph_low)
  call get_simple_variable_cdf(ncid,"azimut_upp",ph_upp)

  r_planet = r_planet*phys_length
  r_low  = r_low*phys_length
  r_upp  = r_upp*phys_length
!  th_low = th_low*180.0/pi
!  th_upp = th_upp*180.0/pi
!  ph_low = ph_low*180.0/pi
!  ph_upp = ph_upp*180.0/pi
  write(*,*)"s_min",s_min,"s_max",s_max
  write(*,*)"size_r",size_r,"size_th",size_th,"size_ph",size_ph
  write(*,*)"r_planet",r_planet,"s_centr",s_centr

  write(*,*)"r_low",r_low
  write(*,*)"r_upp",r_upp
  write(*,*)"th_low",th_low
  write(*,*)"th_upp",th_upp
  write(*,*)"ph_low",ph_low
  write(*,*)"ph_upp",ph_upp

!  deallocate(Dn)
!  deallocate(r_low,r_upp)
!  deallocate(th_low,th_upp)
!  deallocate(ph_low,ph_upp)
  stId = nf90_close(ncid)
end subroutine read_environ_netcdf


END MODULE environ
