module sudakov_radiators
  use types; use consts_dp
  use hoppet_v1
  use rad_tools
  
  implicit none

  private
  
  public :: Sudakov_pt
  public :: expsudakov_pt, exp2sudakov_pt
  
contains


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
         & + (L**2*(Asud(2) - 2*Bsud(1)*beta0*pi))/(2.*pi**2)
    return
  end function exp2sudakov_pt

end module sudakov_radiators
