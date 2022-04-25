!========================================================
!--------------------------------------------------------
! Includes a module providing the resummation ingredients
! -------------------------------------------------------
!========================================================
module rad_tools
  use types; use consts_dp
  use hoppet_v1; use internal_parameters
  implicit none
  private
  

  ! a type to conveniently hold the coupling and information on 
  ! scale ratios
  type process_and_parameters
     character(len=5) :: collider, proc
     real(dp) :: rts
     real(dp) :: alphas_muR, as2pi
     real(dp) :: M, muR, muF, Q
     real(dp) :: ln_muR2_M2
     real(dp) :: ln_Q2_M2
     real(dp) :: ln_Q2_muR2
     real(dp) :: ln_muF2_M2
     real(dp) :: M2_rts2
  end type process_and_parameters


  public :: process_and_parameters, set_process_and_parameters, reset_dynamical_parameters
  public :: init_anom_dim, virtual_scale_dep, update_b2_proc_dep, compare_process_and_parameters, reset_profiled_scales_parameters
  public :: init_profiled_scales_parameters
  real(dp), public :: A(3), B(3), H(2), B2Rad, B3Rad
  real(dp), public :: Asud(3), Bsud(3), as_pow  
  real(dp), public :: CC, BB 
  type(process_and_parameters), public, save :: cs
  logical, save :: init_anom_dim_called = .false.
  real(dp), public, save :: Q0, npow
  logical,  public, save :: profiled_scales
  
contains

  subroutine init_profiled_scales_parameters(use_profiled_scales, Q0_scale, npower)
    real(dp), intent(in) :: Q0_scale, npower
    logical,  intent(in) :: use_profiled_scales

    profiled_scales = use_profiled_scales
    Q0 = Q0_scale
    npow = npower
    return    
  end subroutine init_profiled_scales_parameters

  subroutine reset_profiled_scales_parameters(Q0_scale)
    real(dp), intent(in) :: Q0_scale

    Q0 = Q0_scale
    return    
  end subroutine reset_profiled_scales_parameters


  ! Subroutines to compare process_and_parameters structures
  logical function compare_process_and_parameters(cs1,cs2) result(res)
    type(process_and_parameters) :: cs1,cs2
    res = cs1%rts == cs2%rts .and. &
         &  cs1%alphas_muR == cs2%alphas_muR .and. &
         &  cs1%as2pi == cs2%as2pi .and. &
         &  cs1%M == cs2%M .and. &
         &  cs1%muR == cs2%muR .and. &
         &  cs1%muF == cs2%muF .and. &
         &  cs1%Q == cs2%Q .and. &
         &  cs1%ln_muR2_M2 == cs2%ln_muR2_M2 .and. &
         &  cs1%ln_Q2_M2 == cs2%ln_Q2_M2 .and. &
         &  cs1%ln_Q2_MuR2 == cs2%ln_Q2_MuR2 .and. &
         &  cs1%ln_muF2_M2 == cs2%ln_muF2_M2 .and. &
         &  cs1%M2_rts2 == cs2%M2_rts2
  end function compare_process_and_parameters
  
    !======================================================================
  subroutine init_anom_dim
    implicit none
  !======================================================================
    real(dp) :: Lt, Gamma_4
    real(dp), parameter :: zeta4=pi**4/90._dp, zeta5=1.0369277551433699263313654864570d0
    logical, save :: flg_B3ON

    init_anom_dim_called = .true.
    flg_B3ON = .false.

    if(trim(cs%proc).ne.'qq'.and.trim(cs%proc).ne.'gg') then 
       call wae_error('ERROR: init_anom_dim, unrecognized cs%proc value '//cs%proc)
    endif
   
    select case(trim(cs%proc))
    case('gg')
       as_pow = two
       A(1)  = ca_def/pi
       A(2)  = cmw_K*A(1)/twopi 
       B(1)  = -two*twopi*beta0/twopi
       H(1)  = (ca_def*(5._dp+7._dp/6._dp*pisq)-3._dp*cf_def)/twopi
       B2Rad = (-two*(ca_def**2*(8._dp/3._dp+three*zeta3)-cf_def*tf_def-four/three*ca_def*tf_def) &
            & +twopi_beta0*zeta2*ca_def)/twopi**2
       B2Rad = B2Rad + two*A(1)**2*zeta3
       B(2)  = B2Rad + H(1)*beta0
       !B(2)  = B(2) - Hvirt + NewVirt ! add Hvirt
       A(3)  = (cmw_K2*A(1)*twopi+pi*beta0*ca_def*(ca_def*(808._dp/27._dp-28._dp*zeta3)-224._dp/27._dp*tf_def))/twopi**3
       B(3)  = zero
    case('qq')
       as_pow = zero
       A(1)  = cf_def/pi
       A(2)  = cmw_K*A(1)/twopi
       B(1)  = -three*cf_def/twopi
       H(1)  = cf_def*(-8._dp + 7._dp/6._dp*pisq)/twopi
       H(2)  = (-57433d0/972d0+281d0/162d0*pi**2+22d0/27d0*pi**4+1178d0/27d0*zeta3)/twopi**2
       ! add the additional term from our Sudakov (in direct space instead of b-space)
       H(2)  = H(2) + (8._dp*A(1)*twopi*pi*beta0/3._dp*zeta3)/twopi**2
       ! add the new O(as^2) term coming from the double emission integral (MiNNLO only)
       H(2)  = H(2) + (two*zeta3*A(1)*twopi*B(1)*twopi)/twopi**2
       B2Rad = (-two*(cf_def**2*(-half*pisq+3._dp/8._dp+6._dp*zeta3)&
            & + cf_def*ca_def*(11._dp/18._dp*pisq+17._dp/24._dp-three*zeta3) &
            & + cf_def*tf_def*(-one/6._dp-two/9._dp*pisq)) &
            & + twopi_beta0*zeta2*cf_def)/twopi**2
       B2Rad = B2Rad + two*A(1)**2*zeta3
       B(2)  = B2Rad + H(1)*beta0
       !B(2)  = B(2) - DYvirt + NewVirt ! add DYvirt
       A(3)  = (cmw_K2*A(1)*twopi+pi*beta0*cf_def*(ca_def*(808._dp/27._dp-28._dp*zeta3)-224._dp/27._dp*tf_def))/twopi**3
       B(3)  = zero
!!$       if(flg_B3ON) then
!!$          B3Rad = ((ca_def*cf_def*nf*(-47.51165980795611d0 + (5188d0*Pi**2)/243d0 + (44d0*Pi**4)/45d0 - (3856d0*zeta3)/27d0))/2d0 - &
!!$               & cf_def**2*nf*(63.370370370370374d0 - (8d0*Pi**4)/45d0 - (304d0*zeta3)/9d0) - &
!!$               & ca_def*cf_def*nf*(85.90672153635117d0 - (412d0*Pi**2)/243d0 + (2d0*Pi**4)/27d0 - (904d0*zeta3)/27d0) - cf_def*nf**2*(-2.5459533607681757d0 - (32d0*zeta3)/9d0) + &
!!$               & (cf_def*nf**2*(26.52400548696845d0 - (80d0*Pi**2)/27d0 - (64d0*zeta3)/27d0))/4d0 + &
!!$               & (cf_def**2*nf*(218.74074074074073d0 - (52d0*Pi**2)/9d0 - (56d0*Pi**4)/27d0 + (1024d0*zeta3)/9d0))/2d0 + &
!!$               & ca_def**2*cf_def*(-95.5727023319616d0 - (7163d0*Pi**2)/243d0 - (83d0*Pi**4)/45d0 + (7052d0*zeta3)/9d0 - (88d0*Pi**2*zeta3)/9d0 - 272d0*zeta5) + &
!!$               & ca_def*cf_def**2*(-75.5d0 + (410d0*Pi**2)/9d0 + (494d0*Pi**4)/135d0 - (1688d0*zeta3)/3d0 - (16d0*Pi**2*zeta3)/3d0 - 240d0*zeta5) - &
!!$               & ca_def**2*cf_def*(-407.4471879286694d0 + (3196d0*Pi**2)/243d0 + (77d0*Pi**4)/135d0 + (12328d0*zeta3)/27d0 - (88d0*Pi**2*zeta3)/9d0 - 192d0*zeta5) + &
!!$               & cf_def**3*(-29d0 - 6d0*Pi**2 - (16d0*Pi**4)/5d0 - 136d0*zeta3 + (32d0*Pi**2*zeta3)/3d0 + 480d0*zeta5))/8d0
!!$          B3Rad/twopi**3 ! correct normalization to MiNLO convention
!!$!          B3Rad = B3Rad + two*A(2)**2*zeta3 MW: add the correct zeta pieces here from radish test code!
!!$          B(3) = B3Rad + 2d0*H(2)*beta0 + 2d0*H(1)*beta1/twopi - H(1)**2*beta0 ! MW: from Mathematica notebook, not sure about factors of 2!
!!$       endif
       
    case default
       call wae_error("ERROR: init_proc cs%proc: "//cs%proc)
       stop
    end select

    
    ! convert to our convention for the anomalous dimensions from MiNLO
    A(1) = twopi*A(1); B(1) = twopi*B(1)
    H(1) = twopi*H(1)
    A(2) = twopi**2*A(2); B(2) = twopi**2*B(2)
    H(2) = twopi**2*H(2)
    A(3) = twopi**3*A(3); B(3) = twopi**3*B(3)
    
    !write(*,*) 'rad_tools: Resummation coefficients'
    !write(*,*) '{A1 -> ', A(1),', B1 ->',B(1),', A2 -> ', A(2),', B2 ->',B(2),', A3 -> ', A(3),', B3 ->',B(3),'}'
    !write(*,*) '{beta0 -> ', beta0,' , beta1 -> ',beta1,' ,beta2 -> ',beta2 ,'}'


!!$       B(3) = B(3) + two*twopi*beta0*B(2)*cs%ln_muR2_M2 &
!!$            & + twopi**2*beta1*B(1)*cs%ln_muR2_M2 + twopi**2*beta0**2*B(1)*cs%ln_muR2_M2**2

!!$       A(3) = A(3) + two*twopi*beta0*A(2)*cs%ln_muR2_M2 &
!!$            & + twopi**2*beta1*A(1)*cs%ln_muR2_M2 + twopi**2*beta0**2*A(1)*cs%ln_muR2_M2**2
      
!!$       B(3) = B(3) + twopi**3*beta0**3*as_pow*cs%ln_muR2_M2**2 &
!!$            & + three*twopi**3*beta0*beta1*as_pow*cs%ln_muR2_M2

    B(2) = B(2) + twopi*B(1)*beta0*cs%ln_muR2_M2 + twopi**2*as_pow*beta0**2*cs%ln_muR2_M2
    A(2) = A(2) + twopi*A(1)*beta0*cs%ln_muR2_M2

    ! add B(3) resummation scale dependence coming from H1 and H2
    B(3) = B(3) &
         & -       twopi   *beta0   *cs%ln_Q2_M2   *(-two*B(2) +           A(2)*cs%ln_Q2_M2) &
         & -       twopi**2*beta0**2*cs%ln_Q2_M2**2*(-    B(1) + one/three*A(1)*cs%ln_Q2_M2) &
         & - half *twopi**2*beta1   *cs%ln_Q2_M2   *(-two*B(1) +           A(1)*cs%ln_Q2_M2)
    ! uncomment to test higher-order dependence!

    ! add B(2) resummation scale dependence coming from H1
    B(2) = B(2) - pi*beta0*cs%ln_Q2_M2*(-two*B(1) + A(1)*cs%ln_Q2_M2)
    

    Asud(:) = A(:); Bsud(:) = B(:)
    ! add additional resummation scale dependence from redefining B(as) in the integrand
    Bsud(:) = Bsud(:) - Asud(:)*cs%ln_Q2_M2

  end subroutine init_anom_dim

  !======================================================================
  subroutine virtual_scale_dep(msqB, msqV1, msqV2)
    implicit none
    real(dp), intent(inout) :: msqB(-nf_lcl:nf_lcl,-nf_lcl:nf_lcl)
    real(dp), intent(inout) :: msqV1(-nf_lcl:nf_lcl,-nf_lcl:nf_lcl), msqV2(-nf_lcl:nf_lcl,-nf_lcl:nf_lcl)
  !======================================================================
    if(.not. init_anom_dim_called) then
       write(*,*) ' virtual_scale_dep: called before init_anom_dim, exiting...'
       call exit(-1)
    endif
    ! now handle the virtual corrections, and implement their scale dependence
    msqV2 = msqV2 &
         ! add the additional term from our Sudakov
         & + 8._dp*A(1)*pi*beta0/3._dp*zeta3*msqB &
         ! add the new O(as^2) term coming from the double emission integral (MiNNLO only)
         & - two*zeta3*A(1)*B(1)*msqB &
         ! add renormalisation scale depedence
         & - as_pow*four*(-half*(one+as_pow)*pisq*beta0**2*cs%ln_muR2_M2**2-pisq*beta1*cs%ln_muR2_M2)*msqB &
         & + two*(one+as_pow)*pi*beta0*cs%ln_muR2_M2*msqV1
    ! now add resummation scale dependence
    msqV2 = msqV2 + &
         & (B(1)*cs%ln_Q2_M2 - (A(1)*cs%ln_Q2_M2**2)/2._dp)*msqV1 + &
         & ((A(1)**2*cs%ln_Q2_M2**4)/8._dp + &
         & cs%ln_Q2_M2**3*(-(A(1)*B(1))/2._dp - (A(1)*beta0*pi)/3._dp) + &
         & cs%ln_Q2_M2**2*(-A(2)/2._dp + B(1)**2/2._dp + B(1)*beta0*pi - &
         & A(1)*as_pow*beta0*cs%ln_muR2_M2*pi) + &
         & cs%ln_Q2_M2*(2*as_pow*B(1)*beta0*cs%ln_muR2_M2*pi))*msqB + &
         ! finally add the B2 term which DOES NOT have to contain the Q dependence
         ! (subtract it explicitly from the definition in init_anom_dim) 
         & cs%ln_Q2_M2*(B(2) + pi*beta0*cs%ln_Q2_M2*(-two*B(1) + A(1)*cs%ln_Q2_M2))*msqB 

    
    msqV1 = msqV1 + &
         ! add renormalisation scale dependence
         & two*as_pow*pi*beta0*cs%ln_muR2_M2*msqB
    ! now add resummation scale dependence
    msqV1 = msqV1 + (-half*A(1)*cs%ln_Q2_M2+B(1))*cs%ln_Q2_M2*msqB

    ! additional term to compensate for the difference in the scale of the coupling w.r.t. RadISH
    ! (this comes from expanding as[KR*Q*Exp[-L]] about as[KR*M*Exp[-L]]
    msqV2 = msqV2 - twopi*beta0*cs%ln_Q2_M2*msqV1

  end subroutine virtual_scale_dep

  subroutine update_b2_proc_dep()
    implicit none
    write(*,*) 'error: update_b2_proc_dep is not supposed to be used'
    stop
  end subroutine update_b2_proc_dep

  subroutine set_process_and_parameters(collider, proc, rts, M, KR, KF, KQ)
    implicit none
    character(len=*),          intent(in)  :: collider, proc
    real(dp),                  intent(in)  :: rts, M, KR, KF, KQ
    !-------------------------------------------------------
    real(dp) :: muR, muF, Q
    cs%collider   = collider
    cs%proc       = proc
    cs%rts        = rts
    cs%M          = M
    cs%Q          = KQ*M
    cs%muR        = KR*M
    cs%muF        = KF*M
    cs%alphas_muR = RunningCoupling(cs%muR)
    cs%as2pi      = cs%alphas_muR/twopi
    cs%ln_muR2_M2 = 2*log(cs%muR/M)
    cs%ln_Q2_M2   = 2*log(cs%Q/M)
    cs%ln_Q2_muR2 = 2*log(cs%Q/cs%muR)
    cs%ln_muF2_M2 = 2*log(cs%muF/M)
    cs%M2_rts2    = M**2/rts**2
  end subroutine set_process_and_parameters

  !======================================================================
  ! 
  subroutine reset_dynamical_parameters(Msinglet, KR, KF, KQ)
    real(dp),                  intent(in)  :: Msinglet, KR, KF, KQ
    !-------------------------------------------------------
    cs%M          = Msinglet
    cs%muR        = KR*cs%M
    cs%muF        = KF*cs%M
    cs%Q          = KQ*cs%M 
    cs%alphas_muR = RunningCoupling(cs%muR)
    cs%as2pi      = cs%alphas_muR/twopi
    cs%ln_Q2_muR2 = 2*log( cs%Q   / cs%muR )
    cs%ln_muF2_M2 = 2*log( cs%muF / cs%M)
    cs%ln_muR2_M2 = 2*log( cs%muR / cs%M)
    cs%ln_Q2_M2   = 2*log( cs%Q   / cs%M)
    cs%M2_rts2    = cs%M**2/cs%rts**2
  end subroutine reset_dynamical_parameters

end module rad_tools
