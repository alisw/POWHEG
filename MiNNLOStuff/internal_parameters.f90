!========================================================
!--------------------------------------------------------
! Includes a module providing the resummation ingredients
! -------------------------------------------------------
!========================================================
module internal_parameters
  use types; use consts_dp
  use hoppet_v1, EvolvePDF_hoppet => EvolvePDF, InitPDF_hoppet => InitPDF 
  implicit none

  ! some useful parameters
  ! EW masses
  real(dp), public :: MZ=91.1876_dp
  ! top quark mass
  real(dp), public :: mt = 173.5
  ! set default number of flavours
  integer, public  :: nf_lcl=5
  
  public :: set_masses, set_nflav
contains

  subroutine set_nflav(nf_int)
    integer, intent(in) :: nf_int
    nf_lcl = nf_int
  end subroutine set_nflav
  
  subroutine set_masses(MZ_in, mt_in)
    real(dp), intent(in) :: MZ_in
    real(dp), optional, intent(in) :: mt_in
    MZ = MZ_in
    if (present(mt_in)) then
       mt = mt_in
    endif
  end subroutine set_masses
  
end module internal_parameters
