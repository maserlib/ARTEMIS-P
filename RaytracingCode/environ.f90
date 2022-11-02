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
  integer :: it,l,ii,jj,kk, N

!  B(1,:)=0. 
!  B(2,:)=0.
! B(3,:)=B0*(V(3,:)/z0B)**(-3)
! B(3,:)=0. !sqrt(2.0)*79e-9 ! Tesla from lathys env_ganymede

!write(*,*)"N=",Nray,size(V(1,:))
  N=size(V(1,:))
  do it=1,N
  ii=0;jj=0;kk=0
    do l=1,size_x-1
      if(x_axis(l) <= V(1,it) .and. V(1,it) < x_axis(l+1))then
        ii=l
      endif
    enddo
    
    do l=1,size_y-1
      if(y_axis(l) <= V(2,it) .and. V(2,it) < y_axis(l+1))then
        jj=l
      endif
    enddo

    do l=1,size_z-1
      if(z_axis(l) <= V(3,it) .and. V(3,it) < z_axis(l+1))then
        kk=l
      endif
    enddo
    if(ii==0 .or. jj==0 .or. kk==0) then
      write(*,*)"error detected in magneticf"
      write(*,*)"xyz (curr) = ",V(:,it)
      write(*,*)"xyz (min) = ",x_axis(1),y_axis(1),z_axis(1)
      write(*,*)"xyz (min) = ",x_axis(size_x),y_axis(size_y),z_axis(size_z)
      stop
    else
      B(1,it)=Bx(ii,jj,kk)
      B(2,it)=By(ii,jj,kk)
      B(3,it)=Bz(ii,jj,kk)
    endif
!write(*,*)V(:,it)," B = ",B(:,it),ii,jj,kk!,Bx(75,105,50),By(75,105,50),Bz(75,105,50)!x_axis(ii),y_axis(jj),z_axis(kk)
  enddo
  
!  B(1,:)=0. 
!  B(2,:)=0. 
!  B(3,:)=0. 

end subroutine magneticf

!----------------------------------------------
! Charge l'environnement magnétique à partir d'un fichier netcdf
!----------------------------------------------
subroutine read_magnetic_field_netcdf()
  use netcdf
  use defs_basic_cdf
  integer(kind=4) :: ncid,stId
  !real(kind=8),dimension(:),allocatable :: x_axis,y_axis,z_axis
  !real(kind=8),dimension(:,:,:),allocatable :: Bx,By,Bz
  real(kind=8) :: gstep(3),s_centr(3),phys_length
  integer(kind=4) :: ii

  stId = nf90_open("../Magw_PG_19_03_16_t00250.nc",nf90_nowrite,ncid)
  call test_cdf(stId)
  
  call get_simple_dimens_cdf(ncid,"size_x",size_x)
  call get_simple_dimens_cdf(ncid,"size_y",size_y)
  call get_simple_dimens_cdf(ncid,"size_z",size_z)
  
!  call get_simple_variable_cdf(ncid,"s_min",s_min)
!  call get_simple_variable_cdf(ncid,"s_max",s_max)
  call get_simple_variable_cdf(ncid,"gstep",gstep)
  call get_simple_variable_cdf(ncid,"phys_length",phys_length)
  !call get_simple_variable_cdf(ncid,"r_planet",r_planet)
  call get_simple_variable_cdf(ncid,"s_centr",s_centr)

  allocate(Bx(size_x,size_y,size_z))
  allocate(By(size_x,size_y,size_z))
  allocate(Bz(size_x,size_y,size_z))

  allocate(x_axis(size_x),y_axis(size_y),z_axis(size_z))

  call get_simple_variable_cdf(ncid,"Bx",Bx)
  call get_simple_variable_cdf(ncid,"By",By)
  call get_simple_variable_cdf(ncid,"Bz",Bz)
  !call get_simple_variable_cdf(ncid,"X_axis",X_axis)
  !call get_simple_variable_cdf(ncid,"Y_axis",Y_axis)
  !call get_simple_variable_cdf(ncid,"Z_axis",Z_axis)

  stId = nf90_close(ncid)

  x_axis = 0
  y_axis = 0
  z_axis = 0

  do ii=1,size_x
    x_axis(ii)=((ii-1)*gstep(1)*phys_length - s_centr(1)*phys_length)!*1e3
  enddo

  do ii=1,size_y
    y_axis(ii)=((ii-1)*gstep(2)*phys_length - s_centr(2)*phys_length)!*1e3

  enddo


  do ii=1,size_x
    z_axis(ii)=((ii-1)*gstep(3)*phys_length - s_centr(3)*phys_length)!*1e3
  enddo
  !write(*,*)"x_axis = ",x_axis  
  !write(*,*)"y_axis = ",y_axis  
  !write(*,*)"z_axis = ",z_axis  
!   write(*,*)"Bx=",Bx
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
  integer :: it,l,ii,jj,kk, N

  N = size(Ne)
  allocate(r(N),th(N),ph(N))
  r(:)=0 ; th(:)=0 ; ph(:) = 0


 ! where(V(1,:)==0)
 !   ph(:)=pi/2
 ! elsewhere
 !   ph(:) = mod(atan(V(2,:)/V(1,:))+4*pi,2*pi)
 ! end where
 

  r(:)  = sqrt(V(1,:)*V(1,:) + V(2,:)*V(2,:) + V(3,:)*V(3,:))
  th(:) = mod(acos(V(3,:)/r(:))+4*pi, 2*pi)

!  write(*,*)"ici r_ ",r(:), V(:,:)
!  write(*,*)"ici th ",th(:)
!  write(*,*)"ici ph ",ph(:)
  do it=1,N

    if (V(1,it) .ne. 0.0) then
       if (V(1,it) > 0.0) then
         ph(it) = atan(V(2,it)/V(1,it))
       else if (V(2,it) >= 0.0) then
         ph(it) = atan(V(2,it)/V(1,it))+pi
       else 
         ph(it) = atan(V(2,it)/V(1,it))-pi
       endif
    else if(V(2,it) > 0.0) then
       ph(it) = pi
    else if(V(2,it) < 0.0) then
       ph(it) = -pi
    else
       write(*,*) "error in density, x=0 and y=0"
       stop
    endif
    ph(it) = mod(ph(it)+4*pi,2*pi)


    ii=0; jj=0; kk=0
    do l=1,size_r
      if(r_low(l) <= r(it) .and. r(it) < r_upp(l))&
         ii = l
    enddo
    do l=1,size_th
      if(th_low(l) <= th(it) .and. th(it) < th_upp(l))&
        jj = l
    enddo

    do l=1,size_ph
      if(ph_low(l) <= ph(it) .and. ph(it) < ph_upp(l))&
         kk = l
    enddo

    if (ii==0 .or. jj==0 .or. kk==0)then
      write(*,*)"error detected in density"
      write(*,*)"r (curr,min,max) ",r(it),r_low(1),r_upp(size_r)
      write(*,*)"th (curr,min,max) ",th(it),th_low(1),th_upp(size_th)
      write(*,*)"ph (curr,min,max) ",ph(it),ph_low(1),ph_upp(size_ph)
      stop
    endif
      
    Ne(it) = Dn(ii,jj,kk)
 !  write(*,*)"Ne = ",Ne(it)
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
  real(kind=8) :: s_min(3),s_max(3),gstep(3)
  real(kind=8) :: r_planet,s_centr,phys_length
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

  r_planet = r_planet*phys_length*1e3
  r_low  = r_low*phys_length!*1e3
  r_upp  = r_upp*phys_length!*1e3
!  th_low = th_low*180.0/pi
!  th_upp = th_upp*180.0/pi
!  ph_low = ph_low*180.0/pi
!  ph_upp = ph_upp*180.0/pi
!  write(*,*)"s_min",s_min,"s_max",s_max
!  write(*,*)"size_r",size_r,"size_th",size_th,"size_ph",size_ph
!  write(*,*)"r_planet",r_planet,"s_centr",s_centr

!  write(*,*)"r_low",r_low
!  write(*,*)"r_upp",r_upp
!  write(*,*)"th_low",th_low
!  write(*,*)"th_upp",th_upp
!  write(*,*)"ph_low",ph_low
!  write(*,*)"ph_upp",ph_upp

!  deallocate(Dn)
!  deallocate(r_low,r_upp)
!  deallocate(th_low,th_upp)
!  deallocate(ph_low,ph_upp)
  stId = nf90_close(ncid)
end subroutine read_environ_netcdf


END MODULE environ
