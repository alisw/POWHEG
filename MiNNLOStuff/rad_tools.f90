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
  public :: init_anom_dim, virtual_scale_dep, update_b2_proc_dep, compare_process_and_parameters
  real(dp), public :: A(3), B(3), H(1), B2Rad
  real(dp), public :: Asud(3), Bsud(3), as_pow  
  real(dp), public :: CC, BB 
  type(process_and_parameters), public, save :: cs
  logical, save :: init_anom_dim_called = .false.
  
contains


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
    real(dp), parameter :: zeta4=pi**4/90._dp, zeta5=1.0369277551433699263313654864570

    init_anom_dim_called = .true.
    
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
       A(3)  = (cmw_K2*A(1)+pi*beta0*ca_def*(ca_def*(808._dp/27._dp-28._dp*zeta3)-224._dp/27._dp*tf_def))/twopi**3
       B(3)  = zero
    case('qq')
       as_pow = zero
       A(1)  = cf_def/pi
       A(2)  = cmw_K*A(1)/twopi
       B(1)  = -three*cf_def/twopi
       H(1)  = cf_def*(-8._dp + 7._dp/6._dp*pisq)/twopi
       B2Rad = (-two*(cf_def**2*(-half*pisq+3._dp/8._dp+6._dp*zeta3)&
            & + cf_def*ca_def*(11._dp/18._dp*pisq+17._dp/24._dp-three*zeta3) &
            & + cf_def*tf_def*(-one/6._dp-two/9._dp*pisq)) &
            & + twopi_beta0*zeta2*cf_def)/twopi**2
       B2Rad = B2Rad + two*A(1)**2*zeta3
       B(2)  = B2Rad + H(1)*beta0
       !B(2)  = B(2) - DYvirt + NewVirt ! add DYvirt
       A(3)  = (cmw_K2*A(1)+pi*beta0*cf_def*(ca_def*(808._dp/27._dp-28._dp*zeta3)-224._dp/27._dp*tf_def))/twopi**3
       B(3)  = zero
    case default
       call wae_error("ERROR: init_proc cs%proc: "//cs%proc)
       stop
    end select

    
    ! convert to our convention for the anomalous dimensions from MiNLO
    A(1) = twopi*A(1); B(1) = twopi*B(1)
    H(1) = twopi*H(1)    
    A(2) = twopi**2*A(2); B(2) = twopi**2*B(2)
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

    Asud(:) = A(:); Bsud(:) = B(:)

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
         & + two*zeta3*A(1)*B(1)*msqB &
         ! add renormalisation scale depedence
         & - as_pow*four*(-half*(one+as_pow)*pisq*beta0**2*cs%ln_muR2_M2**2-pisq*beta1*cs%ln_muR2_M2)*msqB &
         & + two*(one+as_pow)*pi*beta0*cs%ln_muR2_M2*msqV1

    msqV1 = msqV1 &
         ! add renormalisation scale dependence
         & + two*as_pow*pi*beta0*cs%ln_muR2_M2*msqB

  end subroutine virtual_scale_dep

  subroutine update_b2_proc_dep(H1proc)
    implicit none
    real(dp), intent(in) :: H1proc
    ! This for Higgs  and Z is not called. To be worked more
    ! for more complex processes (when H1proc becomes flavour and kinematic dependent)
    B(2) = B2Rad + beta0*H1proc
    B(2) = twopi**2*B(2)
    B(2) = B(2) + twopi*B(1)*beta0*cs%ln_muR2_M2 + twopi**2*as_pow*beta0**2*cs%ln_muR2_M2
    Bsud(2) = B(2)
  end subroutine update_b2_proc_dep

  subroutine set_process_and_parameters(collider, proc, rts, M, KR, KF)
    implicit none
    character(len=*),          intent(in)  :: collider, proc
    real(dp),                  intent(in)  :: rts, M, KR, KF
    !-------------------------------------------------------
    real(dp) :: muR, muF, Q
    cs%collider   = collider
    cs%proc       = proc
    cs%rts        = rts
    cs%M          = M
    cs%Q          = M
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
  subroutine reset_dynamical_parameters(Msinglet, KR, KF)
    real(dp),                  intent(in)  :: Msinglet, KR, KF
    !-------------------------------------------------------
    cs%M          = Msinglet
    cs%muR        = KR*cs%M
    cs%muF        = KF*cs%M
    cs%Q          = cs%M 
    cs%alphas_muR = RunningCoupling(cs%muR)
    cs%as2pi      = cs%alphas_muR/twopi
    cs%ln_Q2_muR2 = 2*log( cs%Q   / cs%muR )
    cs%ln_muF2_M2 = 2*log( cs%muF / cs%M)
    cs%ln_muR2_M2 = 2*log( cs%muR / cs%M)
    cs%ln_Q2_M2   = 2*log( cs%Q   / cs%M)
    cs%M2_rts2    = cs%M**2/cs%rts**2
  end subroutine reset_dynamical_parameters

end module rad_tools
