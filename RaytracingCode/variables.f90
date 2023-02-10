module variables
  implicit none

  real(kind=8) :: z0, n0 ! initial density and altitude
  real(kind=8) :: B0,z0B ! initial magnetic value and altitude for the magnetic analytic function

  real(kind=8),allocatable,dimension(:,:,:) :: Dn,Bx,By,Bz ! density in the environment input file (polar coordinates) in cm-3
  real(kind=8),allocatable,dimension(:) :: r_low,r_upp,&  ! limit low and upp of each cells in the environment input file in rad for
                                           th_low,th_upp,ph_low,ph_upp,& ! netcdf density file
                                           x_axis,y_axis,z_axis ! axis of magnetic field cube (netcdf Magw file)
!  real(kind=8) :: s_min(3),s_max(3),gstep(3)  ! minimum and maximum of the grid and step of a cell 
  real(kind=8) :: r_planet!, s_centr           ! radius and center of the planetary object
!  real(kind=8) :: phys_length,gstep                 ! physical length in km (used to normalize)
  integer(kind=4) :: size_r,size_th,size_ph,&   ! number of cells in each direction for netcdf density file
                    size_x,size_y,size_z        ! number of cells in each direction for netcdf magnetic file
  integer(kind=8) :: Nray

end module variables
