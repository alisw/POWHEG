module sudakov_radiators
  use types; use consts_dp
  use hoppet_v1
  use rad_tools
  use frcgauss_intrfc
  
  implicit none

  private
  
  public :: Sudakov_pt, Sudakov_pt_exact
  public :: expsudakov_pt, exp2sudakov_pt
  public :: Smearing, ptDlogSmearing

  logical  :: use_analytic_alphas
  real(dp) :: alphas_freezing_scale
contains

  !-----------------------------------------------------------------------
  ! numerical calculation of the Sudakov radiator
  ! Int_0^Lmod dL Sudakov_integrand(L)
  ! L_mod = modified log of Q/pt
  function Sudakov_pt_exact(L_mod, freezing_scale, analytic_alphas) result(res)
    real(dp), intent(in) :: L_mod, freezing_scale
    logical,  intent(in) :: analytic_alphas
    real(dp)             :: res
    real(dp), parameter  :: dgauss_accuracy = 1e-3

    ! set alphas for integrand
    use_analytic_alphas   = analytic_alphas
    alphas_freezing_scale = freezing_scale

    ! perform numerical integration with dgauss
    res = frcgauss(Sudakov_integrand, zero, L_mod, dgauss_accuracy)
    
    if (use_analytic_alphas) then
       !if (cs%alphas_muR*beta0*L_mod < half) then
       if ((cs%alphas_muR*beta0*L_mod < 0.43_dp).or.(profiled_scales)) then
          res = exp(-res)
       else
          res = zero ! set to zero below the Landau pole
       end if
    else
       res = exp(-res)
    end if
  end function Sudakov_pt_exact

  !-----------------------------------------------------------------------
  ! integrand of the Sudakov
  ! L = log(Q/pt)
  ! cs%alphas_muR = as(KR*Q)
  function Sudakov_integrand(L) result(res)
    real(dp), intent(in) :: L
    real(dp)             :: res, lambda, alphas2pi
    real(dp)             :: A, B, muR, pwhg_alphas
    
    if (use_analytic_alphas) then
       if (profiled_scales) then
          lambda = cs%alphas_muR*beta0 * log(cs%Q / (cs%Q*exp(-L) + Q0 / (one + (cs%Q/Q0*exp(-L))**npow)))
       else
          lambda = cs%alphas_muR*beta0*L
       endif
       alphas2pi = cs%alphas_muR/twopi/ (one-two*lambda) * (one &
            &        - cs%alphas_muR / (one-two*lambda) * beta1/beta0 * log(one-two*lambda))
    else
       if (profiled_scales) then
          !>> linear scaling
          !muR = cs%muR/cs%Q * (cs%Q*exp(-L) + Q0)
          !>> with extra suppression at large pt (not much of a difference)
!          muR = cs%muR/cs%Q * (cs%Q*exp(-L) + Q0 / (one + (cs%Q/Q0*exp(-L))**npow)) --> this was without the kappaQ implementation
          !>> change to implement the resummation scale (KR*pt)
          muR = cs%muR/cs%Q * cs%Q/cs%M * (cs%Q*exp(-L) + Q0 / (one + (cs%Q/Q0*exp(-L))**npow))
          alphas2pi = pwhg_alphas(muR**2, zero, zero)/twopi
          !alphas2pi = RunningCoupling(muR)/twopi          
       else
          !muR = cs%muR * exp(-L)
          !>> change to implement the resummation scale (KR*pt)
          muR = cs%muR * cs%Q/cs%M * exp(-L)
          !>> implement same freezing as in setlocalscales
!          if (muR < alphas_freezing_scale) then --> this was without the kappaQ implementation
          if (muR < alphas_freezing_scale) then
             muR = alphas_freezing_scale
          end if
          alphas2pi = pwhg_alphas(muR**2, zero, zero)/twopi
          !alphas2pi = RunningCoupling(muR)/twopi
       end if
    end if
       
    A   = Asud(1)*alphas2pi + Asud(2)*alphas2pi**2 + Asud(3)*alphas2pi**3
    B   = Bsud(1)*alphas2pi + Bsud(2)*alphas2pi**2 + Bsud(3)*alphas2pi**3
    res = two*(two*A*L + B)
    return
  end function Sudakov_integrand
  
  !-----------------------------------------------------------------------
  ! this is already for two legs
  ! scale dependence is implemented according to RadISH conventions [L = ln(M/pt), as = as(KR*M)]
  function Sudakov_pt(log_M_over_pt, alphas_in) result(res)
    real(dp), intent(in) :: log_M_over_pt
    real(dp), optional, intent(in) :: alphas_in
    !--------------------------
    real(dp) :: L, lambda, res, beta2_MiNLO, as
    !----------------------------------------
    !----------------------------------------   
    res = 0d0
    if(present(alphas_in)) then
       as = alphas_in
    else
       as = RunningCoupling(cs%muR)
    endif
    L = log_M_over_pt
    lambda=as*L*beta0

    if (lambda == zero) then
       res=one ! this is the Sudakov
       return
    end if
           
    if (lambda<half) then
       ! include the L*g1 contribution
       res = L*Asud(1)/pi/beta0*(one+half*log(one-two*lambda)/lambda)

       ! add g2
       res = res + Bsud(1)/twopi/beta0*log(one-two*lambda)&
            & -Asud(2)/four/pisq/beta0**2*(two*lambda/(one-two*lambda)+log(one-two*lambda))&
            & +Asud(1)/two*(pisq*beta1)/(pi*beta0)**3*(half*log(one-two*lambda)**2&
            & +(log(one-two*lambda)+two*lambda)/(one-two*lambda))

       ! add as/Pi*g3 (include terms consistently with MiNLO)
       beta2_MiNLO = beta2
       !beta2_MiNLO = zero
       res = res + as/pi * ( &
            & -half*Asud(3)/8._dp/pisq/beta0**2*(two*lambda/(one-two*lambda))**2 &
            & -Bsud(2)/four/pi/beta0*two*lambda/(one-two*lambda) &
            & +Asud(2)/four*(pisq*beta1)/(pi*beta0)**3*(lambda*(three*two*lambda-two)/(one-two*lambda)**2 &
            & -(one-four*lambda)/(one-two*lambda)**2*log(one-two*lambda)) &
            & +(Bsud(1))/two*(pisq*beta1)/pisq/beta0**2&
            &        *(two*lambda/(one-two*lambda)+log(one-two*lambda)/(one-two*lambda)) &
            & +Asud(1)/two*((pisq*beta1)**2/two/(pi*beta0)**4*(one-four*lambda)/&
            & (one-two*lambda)**2*log(one-two*lambda)**2 &
            & +log(one-two*lambda)*((pi*beta0*pi**3*beta2_MiNLO-(pisq*beta1)**2)/&
            &          (pi*beta0)**4+(pisq*beta1)**2/(pi*beta0)**4/(one-two*lambda))&
            & +one/(pi*beta0)**4*lambda/(one-two*lambda)**2&
            & *(pi*beta0*pi**3*beta2_MiNLO*(two-three*two*lambda)+(pisq*beta1)**2*two*lambda)))

       ! add Bsud(3) term (pure scale variation in MiNNLO)
       res = res - (as/pi)**2 * Bsud(3)/(two*twopi*beta0) * lambda*(one - lambda)/(one - two*lambda)**2

       res = exp(res)
    else
       res=zero ! this is the Sudakov
    end if
    return
  end function Sudakov_pt
  !-----------------------------------------------------------------------

  ! Expansion of the above Sudakov to O(as)
  ! scale dependence is implemented according to MiNLO conventions
  function expsudakov_pt(L) result(res)
    real(dp), intent(in) :: L
    real(dp) :: res

    res =  Asud(1)/pi * L**2 + Bsud(1)*L/pi
    return
  end function expsudakov_pt

  
  ! Expansion of the above Sudakov to O(as^2)
  ! scale dependence is implemented according to MiNLO conventions
  function exp2sudakov_pt(L) result(res)
    real(dp), intent(in) :: L
    real(dp) :: res

    res = (Bsud(2)*L)/(2.*pi**2) - (2*Asud(1)*beta0*L**3)/(3.*pi) &
         & + (L**2*(Asud(2) - 2*Bsud(1)*beta0*pi))/(2.*pi**2) &
         ! adding scale dependence pieces to expand about alphas(KR*M*e^-L) to compute Delta2
         & - Bsud(1)*L*beta0*cs%ln_Q2_M2/pi - Asud(1)*L**2*beta0*cs%ln_Q2_M2/pi
    return
  end function exp2sudakov_pt

  function Smearing(exp_minus_L,Qsmear,modlog_p) result(res)
    real(dp), intent(in) :: exp_minus_L,Qsmear,modlog_p
    real(dp) :: res, pt

    pt = cs%Q * exp_minus_L
    if (modlog_p .gt. 0d0) then
       res = exp(-abs(Qsmear**2/cs%Q**2 * ((1d0,0d0) - cmplx(pt/cs%Q,0d0)**(-modlog_p))**(2d0/modlog_p)))
    else
       res = exp(- Qsmear**2d0/pt**2d0)
    endif
  end function Smearing

  function ptDlogSmearing(exp_minus_L,Qsmear,modlog_p) result(res)
    real(dp), intent(in) :: exp_minus_L,Qsmear,modlog_p
    real(dp) :: res, pt

    pt = cs%Q * exp_minus_L
    if (modlog_p .gt. 0d0) then
       res = 2d0* (cs%Q/pt)**(modlog_p) * (Qsmear/cs%Q)**2 * &
            & abs(-(1d0,0d0) + cmplx(pt/cs%Q,0d0)**(-modlog_p))**(2d0/modlog_p-1d0) ! argument of abs is real, but taken for safety reasons
    else
       res = 2d0*Qsmear**2d0/pt**2d0
    endif
  end function ptDlogSmearing

end module sudakov_radiators
