module defs_basic_cdf

use netcdf

implicit none
private

public  ::                        &
     test_cdf,                    &
      get_simple_dimens_cdf,       &
      get_simple_variable_cdf

 interface get_simple_variable_cdf
  module procedure get_simple_variable_cdf_int
!  module procedure get_simple_variable_cdf_int2_kind1
!  module procedure get_simple_variable_cdf_int2_kind2
!  module procedure get_simple_variable_cdf_int2
  module procedure get_simple_variable_cdf_dp
  module procedure get_simple_variable_cdf_dparr
!  module procedure get_simple_variable_cdf_dparr2
  module procedure get_simple_variable_cdf_dparr3
!  module procedure get_simple_variable_cdf_dparr4
!  module procedure get_simple_variable_cdf_dparr6
!  module procedure get_simple_variable_cdf_string
!  module procedure get_simple_variable_cdf_string2
 end interface
contains
 !!=============================================================
 !!subroutine: m_basic_cdf/test_cdf
 !! FUNCTION
 !!  Test the error flag of netcdf function
 !!
 !! OUTPUT
 !!  Only write
subroutine test_cdf(stId)

  integer(kind=4),intent(in) :: stId

  if (stId /= nf90_noerr) then
   write(*,*)  nf90_strerror(stId)
   stop
  endif

end subroutine test_cdf


 !!=============================================================
 !!subroutine: m_basic_cdf/get_simple_dimens_cdf
 !! FUNCTION
 !!  Get dimension from an open Netcdf file
 !!
 !! OUTPUT
 subroutine get_simple_dimens_cdf(ncid,namedim,dimens)

  integer(kind=4),intent(in)  :: ncid
  integer(kind=4),intent(out) :: dimens
  character(len=*),intent(in)  :: namedim

  integer(kind=4) :: stId, dimid

  stId = nf90_inq_dimid(ncid, namedim, dimid)
  call test_cdf(stId)
  stId = nf90_inquire_dimension(ncid,dimid,len=dimens)
  call test_cdf(stId)

 end subroutine get_simple_dimens_cdf


 !!=============================================================
 !!subroutine: m_basic_cdf/get_simple_variable_cdf_int2_kind2
 !! FUNCTION
 !!  Get variable from an open Netcdf file
 !!
 !! OUTPUT
 !!  Only write
! subroutine get_simple_variable_cdf_int2_kind2(ncid,namevar,values)
!
!  integer,intent(in)  :: ncid
!  integer(i2b),intent(out) :: values(:)
!  character(len=*),intent(in)  :: namevar
!
!  integer(kind=4) :: stId
!  integer :: varid
!
!  stId = nf90_inq_varid(ncid, namevar, varid)
!  call test_cdf(stId)
!  stId = nf90_get_var(ncid,varid,values)
!  call test_cdf(stId)
!
! end subroutine get_simple_variable_cdf_int2_kind2

 !!=============================================================
 !!subroutine: m_basic_cdf/get_simple_variable_cdf_int2_kind1
 !! FUNCTION
 !!  Get variable from an open Netcdf file
 !!
 !! OUTPUT
 !!  Only write
! subroutine get_simple_variable_cdf_int2_kind1(ncid,namevar,values)
!
!  integer,intent(in)  :: ncid
!  integer(i1b),intent(out) :: values(:)
!  character(len=*),intent(in)  :: namevar
!
!  integer :: stId,varid
!
!  stId = nf90_inq_varid(ncid, namevar, varid)
!  call test_cdf(stId)
!  stId = nf90_get_var(ncid,varid,values)
!
! end subroutine get_simple_variable_cdf_int2_kind1
!
!
 !!=============================================================
 !!subroutine: m_basic_cdf/get_simple_variable_cdf_int
 !! FUNCTION
 !!  Get variable from an open Netcdf file
 !!
 !! OUTPUT
 !!  Only write
 subroutine get_simple_variable_cdf_int(ncid,namevar,values)

  integer(kind=4),intent(in)  :: ncid
  integer,intent(out) :: values
  character(len=*),intent(in)  :: namevar

  integer(kind=4) :: stId, varid

  stId = nf90_inq_varid(ncid, namevar, varid)
  call test_cdf(stId)
  stId = nf90_get_var(ncid,varid,values)
  call test_cdf(stId)

 end subroutine get_simple_variable_cdf_int

! !!=============================================================
! !!subroutine: m_basic_cdf/get_simple_variable_cdf_int2
! !! FUNCTION
! !!  Get variable from an open Netcdf file
! !!
! !! OUTPUT
! !!  Only write
! subroutine get_simple_variable_cdf_int2(ncid,namevar,values)
! integer,intent(in)  :: ncid
!  integer,intent(out) :: values(:)
!  character(len=*),intent(in)  :: namevar
!
!  integer :: stId,varid
!
!  stId = nf90_inq_varid(ncid, namevar, varid)
!  call test_cdf(stId)
!  stId = nf90_get_var(ncid,varid,values)
!  call test_cdf(stId)
!
! end subroutine get_simple_variable_cdf_int2
!
 !!=============================================================
 !!subroutine: m_basic_cdf/get_simple_variable_cdf_dp
 !! FUNCTION
 !!  Get variable from an open Netcdf file
 !!
 !! OUTPUT
 !!  Only write
 subroutine get_simple_variable_cdf_dp(ncid,namevar,values)

  integer(kind=4),intent(in)   :: ncid
  real(kind=8),intent(out) :: values
  character(len=*),intent(in)  :: namevar

  integer(kind=4) :: stId,varid

  stId = nf90_inq_varid(ncid, namevar, varid)
  call test_cdf(stId)
  stId = nf90_get_var(ncid,varid,values)
  call test_cdf(stId)

 end subroutine get_simple_variable_cdf_dp

 !!=============================================================
 !!subroutine: m_basic_cdf/get_simple_variable_dparr
 !! FUNCTION
 !!  Get variable from an open Netcdf file
 !!
 !! OUTPUT
 !!  Only write
 subroutine get_simple_variable_cdf_dparr(ncid,namevar,values)

  integer(kind=4),intent(in)   :: ncid
  real(kind=8),intent(out) :: values(:)
  character(len=*),intent(in)  :: namevar

  integer(kind=4) :: stId,varid

  stId = nf90_inq_varid(ncid, namevar, varid)
  call test_cdf(stId)
  stId = nf90_get_var(ncid,varid,values)
  call test_cdf(stId)

 end subroutine get_simple_variable_cdf_dparr

 !!=============================================================
 !!subroutine: m_basic_cdf/get_simple_variable_cdf_dparr2
 !! FUNCTION
 !!  Get variable from an open Netcdf file
 !!
 !! OUTPUT
 !!  Only write
! subroutine get_simple_variable_cdf_dparr2(ncid,namevar,values)
!
! integer,intent(in)   :: ncid
!  real(dp),intent(out) :: values(:,:)
!  character(len=*),intent(in)  :: namevar
!
!  integer :: stId,varid
!
!  stId = nf90_inq_varid(ncid, namevar, varid)
!  call test_cdf(stId)
!  stId = nf90_get_var(ncid,varid,values)
!  call test_cdf(stId)
!
! end subroutine get_simple_variable_cdf_dparr2
!
 !!=============================================================
 !!subroutine: m_basic_cdf/get_simple_variable_cdf_dparr3
 !! FUNCTION
 !!  Get variable from an open Netcdf file
 !!
 !! OUTPUT
 !!  Only write
 subroutine get_simple_variable_cdf_dparr3(ncid,namevar,values)

  integer(kind=4),intent(in)   :: ncid
  real(kind=8),intent(out) :: values(:,:,:)
  character(len=*),intent(in)  :: namevar

  integer(kind=4) :: stId,varid

  stId = nf90_inq_varid(ncid, namevar, varid)
  call test_cdf(stId)
  stId = nf90_get_var(ncid,varid,values)
  call test_cdf(stId)

 end subroutine get_simple_variable_cdf_dparr3

  !!=============================================================
!  !!subroutine: m_basic_cdf/get_simple_variable_cdf_dparr4
!  !! FUNCTION
!  !!  Get variable from an open Netcdf file
!  !!
!  !! OUTPUT
!  !!  Only write
!  subroutine get_simple_variable_cdf_dparr4(ncid,namevar,values)
!
!   integer,intent(in)   :: ncid
!   real(dp),intent(out) :: values(:,:,:,:)
!   character(len=*),intent(in)  :: namevar
!
!   integer :: stId,varid
!
!   stId = nf90_inq_varid(ncid, namevar, varid)
!   call test_cdf(stId)
!   stId = nf90_get_var(ncid,varid,values)
!   call test_cdf(stId)
!
! end subroutine get_simple_variable_cdf_dparr4
!
!  !!=============================================================
!  !!subroutine: m_basic_cdf/get_simple_variable_cdf_dparr6
!  !! FUNCTION
!  !!  Get variable from an open Netcdf file
!  !!
!  !! OUTPUT
!  !!  Only write
!  subroutine get_simple_variable_cdf_dparr6(ncid,namevar,values)
!
!
!   integer,intent(in)   :: ncid
!   real(dp),intent(out) :: values(:,:,:,:,:,:)
!   character(len=*),intent(in)  :: namevar
!
!   integer :: stId,varid
!
!   stId = nf90_inq_varid(ncid, namevar, varid)
!   call test_cdf(stId)
!   stId = nf90_get_var(ncid,varid,values)
!   call test_cdf(stId)
!
! end subroutine get_simple_variable_cdf_dparr6
!
! !!=============================================================
! !!subroutine: m_basic_cdf/get_simple_variable_cdf_string
! !! FUNCTION
! !!  Get variable from an open Netcdf file
! !!
! !! OUTPUT
! !!  Only write
! subroutine get_simple_variable_cdf_string(ncid,namevar,values)
!
!  integer,intent(in)   :: ncid
!  character(len=*),intent(out) :: values
!  character(len=*),intent(in)  :: namevar
!
!  integer :: stId,varid
!
!  stId = nf90_inq_varid(ncid, namevar, varid)
!  call test_cdf(stId)
!  stId = nf90_get_var(ncid,varid,values)
!  call test_cdf(stId)
!
! end subroutine get_simple_variable_cdf_string
!
! !!=============================================================
! !!subroutine: m_basic_cdf/get_simple_variable_cdf_string2
! !! FUNCTION
! !!  Get variable from an open Netcdf file
! !!
! !! OUTPUT
! !!  Only write
! subroutine get_simple_variable_cdf_string2(ncid,namevar,values)
!
!  integer,intent(in)   :: ncid
!  character(len=*),intent(out) :: values(:)
!  character(len=*),intent(in)  :: namevar
!
!  integer :: stId,varid
!
!  stId = nf90_inq_varid(ncid, namevar, varid)
!  call test_cdf(stId)
!  stId = nf90_get_var(ncid,varid,values)
!  call test_cdf(stId)
!
! end subroutine get_simple_variable_cdf_string2
!
!
end module defs_basic_cdf
