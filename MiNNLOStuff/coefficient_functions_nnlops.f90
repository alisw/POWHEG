module coefficient_functions_VH
  use types; use consts_dp; use convolution_communicator; use convolution
  use dglap_objects
  use qcd
  use rad_tools; use internal_parameters
  implicit none
  private

  public :: InitCoeffMatrix
  

contains
   !----------------------------------------------------------------------
   function C1qq(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x

    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       res = CF*(1-x)
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       ! no virtual piece
    case(cc_DELTA)
       res = -pisq/12._dp*CF
    end select

    if (cc_piece /= cc_DELTA) res = res * x     
  end function C1qq

   !----------------------------------------------------------------------
   function C1gq(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x

    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       res = CF*x
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       ! no virtual piece
    case(cc_DELTA)
       res = zero
    end select

    if (cc_piece /= cc_DELTA) res = res * x     
  end function C1gq


   !----------------------------------------------------------------------
   function C1qg(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x

    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       res = 2*TR*x*(1-x)
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       ! no virtual piece
    case(cc_DELTA)
       res = zero
    end select

    if (cc_piece /= cc_DELTA) res = res * x     
  end function C1qg


  !----------------------------------------------------------------------
  function C1gg(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x

    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       res = 0.0_dp
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       ! no virtual piece
    case(cc_DELTA)
       res = -pisq/12._dp*CA
    end select

    if (cc_piece /= cc_DELTA) res = res * x     
  end function C1gg


  ! Below are the two-loop coefficient functions C2 and
  ! the spin-correlation coefficient functions G derived from
  ! Catani-Grazzini 1106.4652 for ptB resummation.
  ! The resummation scheme is set
  ! by separating the hard-collinear coefficient function
  ! from the hard-virtual corrections to the form factor
  ! (this is different from what is adopted in Catani et al. 1311.1654)

  ! G(z) coefficient functions, accounting for spin correlations in
  ! gluonic processes:
   !----------------------------------------------------------------------
   function G1qq(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x

    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       res = zero
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       ! no virtual piece
    case(cc_DELTA)
       res = zero
    end select

    if (cc_piece /= cc_DELTA) res = res * x     
  end function G1qq

  !----------------------------------------------------------------------
   function G1gq(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x

    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       res = two*CF*(1-x)/x
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       ! no virtual piece
    case(cc_DELTA)
       res = zero
    end select

    if (cc_piece /= cc_DELTA) res = res * x     
  end function G1gq


  !----------------------------------------------------------------------
   function G1gg(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x

    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       res = two*CA*(1-x)/x
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       ! no virtual piece
    case(cc_DELTA)
       res = zero
    end select

    if (cc_piece /= cc_DELTA) res = res * x     
  end function G1gg

  !----------------------------------------------------------------------
   function G1qg(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x

    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       res = zero
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       ! no virtual piece
    case(cc_DELTA)
       res = zero
    end select

    if (cc_piece /= cc_DELTA) res = res * x     
  end function G1qg
  

  ! Two-loop coefficient functions C2(z) from Catani-Grazzini et al.
  ! obtained as described in our paper 1705.09127
  !----------------------------------------------------------------------
  function C2qq(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    complex(dp) :: HPL2, HPL3
    complex(dp) :: HPL20, HPL30, HPL31
    real(dp), parameter    :: cutoff=1e-11 ! cutoff to regularise the x->1 singularity in the gg coefficient function
    !------------------------------------------
    ! The following declarations are necessary for hplog
    integer, parameter :: n1=-1
    integer, parameter :: n2= 1
    integer, parameter :: nw= 4
    real(dp)    :: Hr1(n1:n2),Hr2(n1:n2,n1:n2),Hr3(n1:n2,n1:n2,n1:n2),&
         & Hr4(n1:n2,n1:n2,n1:n2,n1:n2)
    real(dp)    :: Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2),&
         & Hi4(n1:n2,n1:n2,n1:n2,n1:n2)
    complex(dp) :: Hc1(n1:n2),Hc2(n1:n2,n1:n2),Hc3(n1:n2,n1:n2,n1:n2), &
         & Hc4(n1:n2,n1:n2,n1:n2,n1:n2)

    
    x = exp(-y)
    res = zero

    
    ! evaluate HPL's using hplog (much faster than Chaplin)
    call hplog(x,nw,Hc1,Hc2,Hc3,Hc4, &
         &     Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2)
    HPL20 = Hc2(0,1)
    HPL30 = Hc3(0,0,1)
    call hplog(one-x,nw,Hc1,Hc2,Hc3,Hc4, &
         &     Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2)
    HPL31 = Hc3(0,0,1)

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       if (one - x > cutoff) then
          res = (CF*(28*nf + CA*(-202 + 189*Zeta3)))/27./(one-x) + & ! this is the regular part of the plus distribution 1/(1-x)_+
               ! now add the regular terms
               & (CF*(-344 + 24*Pi**2 + 974*x - 1600*CA*x + 1188*CF*x - 432*CA*HPL30*x + 1080*CF*HPL30*x + 148*nf*x - 60*Pi**2*x + 54*CA*Pi**2*x - 54*CF*Pi**2*x - 1188*x**2 + 1584*CA*x**2 - 2376*CF*x**2 - 72*nf*x**2 + 72*Pi**2*x**2 - 108*CA*Pi**2*x**2 + 108*CF*Pi**2*x**2 + 830*x**3 + 16*CA*x**3 + 1188*CF*x**3 - 432*CA*HPL30*x**3 + 1080*CF*HPL30*x**3 - 76*nf*x**3 - 60*Pi**2*x**3 + 54*CA*Pi**2*x**3 - 54*CF*Pi**2*x**3 - 272*x**4 + 24*Pi**2*x**4 + 216*(CA - CF)*HPL31*x*(1 + x**2) + 1188*CA*x*zeta3 - 1080*CF*x*zeta3 - 324*CA*x**3*zeta3 - 1080*CF*x**3*zeta3 - 36*CA*Pi**2*x*log(1 - x) + 36*CF*Pi**2*x*log(1 - x) - 108*CA*x**2*log(1 - x) + 108*CF*x**2*log(1 - x) + 108*CA*x**3*log(1 - x) - 108*CF*x**3*log(1 - x) - 36*CA*Pi**2*x**3*log(1 - x) + 36*CF*Pi**2*x**3*log(1 - x) - 252*x*log(x) + 348*CA*x*log(x) - 540*CF*x*log(x) - 60*nf*x*log(x) + 612*x**2*log(x) - 432*CA*x**2*log(x) + 1404*CF*x**2*log(x) - 744*x**3*log(x) + 996*CA*x**3*log(x) - 1728*CF*x**3*log(x) - 60*nf*x**3*log(x) + 384*x**4*log(x) - 144*log(1 - x)*log(x) + 360*x*log(1 - x)*log(x) - 216*CA*x*log(1 - x)*log(x) + 648*CF*x*log(1 - x)*log(x) - 432*x**2*log(1 - x)*log(x) + 432*CA*x**2*log(1 - x)*log(x) - 1296*CF*x**2*log(1 - x)*log(x) + 360*x**3*log(1 - x)*log(x) - 216*CA*x**3*log(1 - x)*log(x) + 648*CF*x**3*log(1 - x)*log(x) - 144*x**4*log(1 - x)*log(x) + 216*CA*x*log(1 - x)**2*log(x) - 324*CF*x*log(1 - x)**2*log(x) + 216*CA*x**3*log(1 - x)**2*log(x) - 324*CF*x**3*log(1 - x)**2*log(x) + 27*x*log(x)**2 + 99*CA*x*log(x)**2 - 162*CF*x*log(x)**2 - 18*nf*x*log(x)**2 + 108*CA*x**2*log(x)**2 - 108*CF*x**2*log(x)**2 + 45*x**3*log(x)**2 - 9*CA*x**3*log(x)**2 + 108*CF*x**3*log(x)**2 - 18*nf*x**3*log(x)**2 - 72*x**4*log(x)**2 - 108*CF*x*log(1 - x)*log(x)**2 - 108*CF*x**3*log(1 - x)*log(x)**2 - 18*x*log(x)**3 + 18*CA*x*log(x)**3 - 18*CF*x*log(x)**3 + 18*x**3*log(x)**3 + 18*CA*x**3*log(x)**3 + 18*CF*x**3*log(x)**3 - 72*HPL20*((-1 + x)**2*(2 + (-1 + 3*CA - 6*CF)*x + 2*x**2) - 3*(CA - CF)*x*(1 + x**2)*log(1 - x) - 3*(CA - 3*CF)*x*(1 + x**2)*log(x))))/(216*(-1 + x)*x)
       else
          res = zero
       end if
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       if (one - x > cutoff) then 
          ! subtract the singular part of the 1/(1-x)_+ distribution
          res = res - (CF*(28*nf + CA*(-202 + 189*Zeta3)))/27./(one-x)
       else
          res = res
       end if
    case(cc_DELTA)
       res = (CF*(9*CF*Pi**4 + 2*CA*(4856 - 603*Pi**2 + 18*Pi**4 - 2772*Zeta3) &
            & + 4*nf*(-328 + 45*Pi**2 + 252*Zeta3)))/2592.
    end select

    if (cc_piece /= cc_DELTA) res = res * x     
  end function C2qq


  !----------------------------------------------------------------------
  function C2qqbar(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    complex(dp) :: HPL2, HPL3
    complex(dp) :: HPL20, HPL2_1, HPL30, HPL3_1, HPL31_1px
    real(dp), parameter    :: cutoff=1e-11 ! cutoff to regularise the x->1 singularity in the gg coefficient function
    !------------------------------------------
    ! The following declarations are necessary for hplog
    integer, parameter :: n1=-1
    integer, parameter :: n2= 1
    integer, parameter :: nw= 4
    real(dp)    :: Hr1(n1:n2),Hr2(n1:n2,n1:n2),Hr3(n1:n2,n1:n2,n1:n2),&
         & Hr4(n1:n2,n1:n2,n1:n2,n1:n2)
    real(dp)    :: Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2),&
         & Hi4(n1:n2,n1:n2,n1:n2,n1:n2)
    complex(dp) :: Hc1(n1:n2),Hc2(n1:n2,n1:n2),Hc3(n1:n2,n1:n2,n1:n2), &
         & Hc4(n1:n2,n1:n2,n1:n2,n1:n2)
    
    x = exp(-y)
    res = zero

    ! evaluate HPL's using hplog (much faster than Chaplin)
    call hplog(x,nw,Hc1,Hc2,Hc3,Hc4, &
         &     Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2)
    HPL20 = Hc2(0,1)
    HPL30 = Hc3(0,0,1)
    call hplog(-x,nw,Hc1,Hc2,Hc3,Hc4, &
         &     Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2)
    HPL2_1 = Hc2(0,1)
    HPL3_1 = Hc3(0,0,1)
    call hplog(one/(one+x),nw,Hc1,Hc2,Hc3,Hc4, &
         &     Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
    HPL31_1px = Hc3(0,0,1)

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       if (one - x > cutoff) then
          res = 4*(CF*(((1 - x)*(172 - 143*x + 136*x**2))/(432.*x) + ((21 - 30*x + 32*x**2)*Log(x))/72. - &
               &        ((3 + 3*x + 8*x**2)*Log(x)**2)/96. + ((1 + x)*Log(x)**3)/48. + &
               &        ((1 - x)*(2 - x + 2*x**2)*(HPL20 - Pi**2/6. + Log(1 - x)*Log(x)))/(12.*x)) + &
               &     CF*(-CA/2. + CF)*((Pi**2*(-3 + x))/24. + ((3 + 11*x)*Log(x))/8. + &
               &        (1 - x)*(1.875 + HPL20/2. + (Log(1 - x)*Log(x))/2.) - &
               &        ((1 + x)*(HPL2_1 + Log(x)*Log(1 + x)))/2. + &
               &        ((1 + x**2)*(HPL30 + HPL31_1px + (3*HPL3_1)/2. - (3*zeta3)/4. - (HPL20*Log(x))/2. - &
               &             (HPL2_1*Log(x))/2. - Log(x)**3/24. + (Pi**2*Log(1 + x))/12. + &
               &             (Log(x)**2*Log(1 + x))/4. - Log(1 + x)**3/6.))/(1 + x)))
       else
          res = zero
       end if
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       ! no virtual piece
    case(cc_DELTA)
       res = zero
    end select
    
    if (cc_piece /= cc_DELTA) res = res * x     
  end function C2qqbar

  !----------------------------------------------------------------------
  function C2qqp(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    complex(dp) :: HPL2, HPL3
    complex(dp) :: HPL20
    real(dp), parameter    :: cutoff=1e-11 ! cutoff to regularise the x->1 singularity in the gg coefficient function
    !------------------------------------------
    ! The following declarations are necessary for hplog
    integer, parameter :: n1=-1
    integer, parameter :: n2= 1
    integer, parameter :: nw= 4
    real(dp)    :: Hr1(n1:n2),Hr2(n1:n2,n1:n2),Hr3(n1:n2,n1:n2,n1:n2),&
         & Hr4(n1:n2,n1:n2,n1:n2,n1:n2)
    real(dp)    :: Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2),&
         & Hi4(n1:n2,n1:n2,n1:n2,n1:n2)
    complex(dp) :: Hc1(n1:n2),Hc2(n1:n2,n1:n2),Hc3(n1:n2,n1:n2,n1:n2), &
         & Hc4(n1:n2,n1:n2,n1:n2,n1:n2)
    
    x = exp(-y)
    res = zero

    ! evaluate HPL's using hplog (much faster than Chaplin)
    call hplog(x,nw,Hc1,Hc2,Hc3,Hc4, &
         &     Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2)
    HPL20 = Hc2(0,1)

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       if (one - x > cutoff) then          
          res = 4*CF*(((1 - x)*(172 - 143*x + 136*x**2))/(432.*x) + ((21 - 30*x + 32*x**2)*Log(x))/72. - &
               &    ((3 + 3*x + 8*x**2)*Log(x)**2)/96. + ((1 + x)*Log(x)**3)/48. + &
               &    ((1 - x)*(2 - x + 2*x**2)*(HPL20 - Pi**2/6. + Log(1 - x)*Log(x)))/(12.*x))
       else
          res = zero
       end if
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       ! no virtual piece
    case(cc_DELTA)
       res = zero
    end select
    
    if (cc_piece /= cc_DELTA) res = res * x     
  end function C2qqp

  
  ! Eq. (30) of 1106.4652
  !----------------------------------------------------------------------
  function C2gq(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    complex(dp) :: HPL2, HPL3
    complex(dp) :: HPL20, HPL2_1, HPL30, HPL31, HPL3_1
    real(dp), parameter    :: cutoff=1e-11 ! cutoff to regularise the x->1 singularity in the gg coefficient function
    !------------------------------------------
    ! The following declarations are necessary for hplog
    integer, parameter :: n1=-1
    integer, parameter :: n2= 1
    integer, parameter :: nw= 4
    real(dp)    :: Hr1(n1:n2),Hr2(n1:n2,n1:n2),Hr3(n1:n2,n1:n2,n1:n2),&
         & Hr4(n1:n2,n1:n2,n1:n2,n1:n2)
    real(dp)    :: Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2),&
         & Hi4(n1:n2,n1:n2,n1:n2,n1:n2)
    complex(dp) :: Hc1(n1:n2),Hc2(n1:n2,n1:n2),Hc3(n1:n2,n1:n2,n1:n2), &
         & Hc4(n1:n2,n1:n2,n1:n2,n1:n2)

    x = exp(-y)
    res = zero

    ! evaluate HPL's using hplog (much faster than Chaplin)
    !if (x.ne.one) then !PM debug
       call hplog(x,nw,Hc1,Hc2,Hc3,Hc4, &
            &     Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
       HPL20 = Hc2(0,1)
       HPL30 = Hc3(0,0,1)
       call hplog(-x,nw,Hc1,Hc2,Hc3,Hc4, &
            &     Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
       HPL2_1 = Hc2(0,1)
       HPL3_1 = Hc3(0,0,1)
       call hplog(one/(one+x),nw,Hc1,Hc2,Hc3,Hc4, &
            &     Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
       HPL31 = Hc3(0,0,1)
    !else
    !   ! these numbers are evaluated with HPL 2.0
    !   HPL20 = pisq/6._dp
    !   HPL30 = zeta3
    !   HPL2_1 = -pisq/12._dp
    !   HPL3_1 = -3._dp/4._dp*zeta3
    !   HPL31 = 0.5372131936_dp
    !end if

    ! evaluate HPL's using Chaplin
    !HPL20  = HPL2(0,1,x)
    !HPL2_1 = HPL2(0,1,-x)
    !HPL30  = HPL3(0,0,1,x)
    !HPL31  = HPL3(0,0,1,1/(1 + x))
    !HPL3_1 = HPL3(0,0,1,-x)
    
    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       if (one-x > cutoff) then
          res = 4._dp*(CF*nf*(-0.5185185185185185 + 14/(27.*x) + (13*x)/108. - (5*Log(1 - x))/18. + (5*Log(1 - x))/(18.*x) + (x*Log(1 - x))/18. - &
               &       Log(1 - x)**2/12. + Log(1 - x)**2/(12.*x) + (x*Log(1 - x)**2)/24.) + &
               &    CF**2*(0.625 - x/16. - 2*Log(1 - x) + (2*Log(1 - x))/x + (5*x*Log(1 - x))/8. - (3*Log(1 - x)**2)/8. + (3*Log(1 - x)**2)/(8.*x) + &
               &       (x*Log(1 - x)**2)/16. - Log(1 - x)**3/12. + Log(1 - x)**3/(12.*x) + (x*Log(1 - x)**3)/24. - (15*Log(x))/16. + (5*x*Log(x))/16. - &
               &       Log(x)**2/8. - (3*x*Log(x)**2)/32. + Log(x)**3/24. - (x*Log(x)**3)/48.) + &
               &    CA*CF*(7.324074074074074 - Pi**2/3. - 395/(54.*x) + (11*Pi**2)/(36.*x) - (67*x)/27. + (Pi**2*x)/16. + (38*x**2)/27. - (Pi**2*x**2)/18. - &
               &       (5*zeta3)/2. + (4*zeta3)/x + 2*x*zeta3 + (x*HPL2_1)/4. + 2*HPL20 - (11*HPL20)/(6.*x) - (x*HPL20)/2. + &
               &       (x**2*HPL20)/3. - (3*HPL3_1)/2. - (3*HPL3_1)/(2.*x) - (3*x*HPL3_1)/4. + HPL30/2. - &
               &       (5*HPL30)/(2.*x) - (5*x*HPL30)/4. - HPL31 - HPL31/x - (x*HPL31)/2. + &
               &       (19*Log(1 - x))/9. - (19*Log(1 - x))/(9.*x) - (43*x*Log(1 - x))/72. + (11*Log(1 - x)**2)/24. - (11*Log(1 - x)**2)/(24.*x) - &
               &       (5*x*Log(1 - x)**2)/48. + Log(1 - x)**3/12. - Log(1 - x)**3/(12.*x) - (x*Log(1 - x)**3)/24. - (83*Log(x))/24. + (x*Log(x))/12. - &
               &       (11*x**2*Log(x))/9. + (HPL2_1*Log(x))/2. + (HPL2_1*Log(x))/(2.*x) + (x*HPL2_1*Log(x))/4. - (HPL20*Log(x))/2. + &
               &       (3*HPL20*Log(x))/(2.*x) + (3*x*HPL20*Log(x))/4. + 2*Log(1 - x)*Log(x) - (11*Log(1 - x)*Log(x))/(6.*x) - &
               &       (3*x*Log(1 - x)*Log(x))/4. + (x**2*Log(1 - x)*Log(x))/3. - (Log(1 - x)**2*Log(x))/4. + (Log(1 - x)**2*Log(x))/(4.*x) + &
               &       (x*Log(1 - x)**2*Log(x))/8. + (3*Log(x)**2)/4. + (3*x*Log(x)**2)/16. + (x**2*Log(x)**2)/6. - (Log(1 - x)*Log(x)**2)/4. + &
               &       (Log(1 - x)*Log(x)**2)/(4.*x) + (x*Log(1 - x)*Log(x)**2)/8. - Log(x)**3/12. - (x*Log(x)**3)/24. - (Pi**2*Log(1 + x))/12. - &
               &       (Pi**2*Log(1 + x))/(12.*x) - (Pi**2*x*Log(1 + x))/24. + (x*Log(x)*Log(1 + x))/4. - (Log(x)**2*Log(1 + x))/4. - &
               &       (Log(x)**2*Log(1 + x))/(4.*x) - (x*Log(x)**2*Log(1 + x))/8. + Log(1 + x)**3/6. + Log(1 + x)**3/(6.*x) + (x*Log(1 + x)**3)/12.))
       else
          res = zero
       end if
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       ! no virtual piece
    case(cc_DELTA)
       res = zero
    end select

    if (cc_piece /= cc_DELTA) res = res * x     
  end function C2gq


  ! Derived from Eq. (32) of 1106.4652, using two-loop virtual corrections
  ! from Becher-Neubert 1205.3806 (obtained with nf=5, N=3)
  !----------------------------------------------------------------------
  function C2gg(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x, Lt
    complex(dp) :: HPL2, HPL3
    complex(dp) :: HPL20, HPL21, HPL2_1, HPL30, HPL31, HPL3_1
    real(dp), parameter    :: cutoff=1e-11 ! cutoff to regularise the x->1 singularity in the gg coefficient function
    !------------------------------------------
    ! The following declarations are necessary for hplog
    integer, parameter :: n1=-1
    integer, parameter :: n2= 1
    integer, parameter :: nw= 4
    real(dp)    :: Hr1(n1:n2),Hr2(n1:n2,n1:n2),Hr3(n1:n2,n1:n2,n1:n2),&
         & Hr4(n1:n2,n1:n2,n1:n2,n1:n2)
    real(dp)    :: Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2),&
         & Hi4(n1:n2,n1:n2,n1:n2,n1:n2)
    complex(dp) :: Hc1(n1:n2),Hc2(n1:n2,n1:n2),Hc3(n1:n2,n1:n2,n1:n2), &
         & Hc4(n1:n2,n1:n2,n1:n2,n1:n2)
    
    x = exp(-y)
    res = zero

    ! evaluate HPL's using hplog (much faster than Chaplin)
    !if (x.ne.one) then ! PM debug
       call hplog(x,nw,Hc1,Hc2,Hc3,Hc4, &
            &     Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2)
       HPL20 = Hc2(0,1)
       HPL30 = Hc3(0,0,1)
       call hplog(-x,nw,Hc1,Hc2,Hc3,Hc4, &
            &     Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2)
       HPL2_1 = Hc2(0,1)
       HPL3_1 = Hc3(0,0,1)
       call hplog(one-x,nw,Hc1,Hc2,Hc3,Hc4, &
            &     Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2)
       HPL21 = Hc2(0,1)
       call hplog(x/(one+x),nw,Hc1,Hc2,Hc3,Hc4, &
            &     Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
       HPL31 = Hc3(0,0,1)
    !else
    !   ! these numbers are evaluated with HPL 2.0
    !   HPL20 = pisq/6._dp
    !   HPL30 = zeta3
    !   HPL2_1 = -pisq/12._dp
    !   HPL3_1 = -3._dp/4._dp*zeta3
    !   HPL21 = zero
    !   HPL31 = 0.5372131936_dp
    !end if

    ! evaluate HPL's using Chaplin
    !HPL20  = HPL2(0,1,x)
    !HPL2_1 = HPL2(0,1,-x)
    !HPL21  = HPL2(0,1,1-x)
    !HPL30  = HPL3(0,0,1,x)
    !HPL31  = HPL3(0,0,1,x/(1 + x))
    !HPL3_1 = HPL3(0,0,1,-x)

    Lt = log(cs%M**2/mt**2)

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       if (one - x > cutoff) then
          res =4._dp*((7._dp*CA*nf)/27._dp + CA**2*(-1.8703703703703705_dp + (7._dp*Zeta3)/4._dp))/(one-x) +& ! this is the regular part of the plus distribution 1/(1-x)_+
               ! now add the regular terms
               &       4._dp*(CA*nf*(-0.7685185185185185 + 121/(216.*x) + (55*x)/108. - (139*x**2)/216. - (x*Log(1 - x))/24. + (13*Log(x))/72. + (5*x*Log(x))/36. + &
               &       Log(x)**2/24. + (x*Log(x)**2)/24.) + CF*nf*(2 - 1/(12.*x) - 2*x + x**2/12. + (3*Log(x))/4. + (3*x*Log(x))/4. + (3*Log(x)**2)/16. + &
               &       (x*Log(x)**2)/16. + Log(x)**3/24. + (x*Log(x)**3)/24.) + &
               &       CA**2*(8.574074074074074 - 395/(54.*x) - (190*x)/27. + (835*x**2)/108. - (17*zeta3)/(4.*(1 + x)) - zeta3/(2.*(1 - x)*(1 + x)) + &
               &       zeta3/(2.*x*(1 + x)) + (5*zeta3)/(2.*(1 - x)*x*(1 + x)) - (11*x*zeta3)/(2.*(1 + x)) + (5*x*zeta3)/(2.*(1 - x)*(1 + x)) - &
               &       (5*x**2*zeta3)/(2.*(1 + x)) + (x**2*zeta3)/(2.*(1 - x)*(1 + x)) - (3*x**3*zeta3)/(1 + x) - (5*x**3*zeta3)/(2.*(1 - x)*(1 + x)) + &
               &       (x**4*zeta3)/(2.*(1 - x)*(1 + x)) - 2*HPL21 + (11*HPL21)/(6.*x) + 2*x*HPL21 - (11*x**2*HPL21)/6. - &
               &       HPL3_1/(1 + x) - HPL3_1/(2.*x*(1 + x)) - (3*x*HPL3_1)/(2.*(1 + x)) - (x**2*HPL3_1)/(1 + x) - &
               &       (x**3*HPL3_1)/(2.*(1 + x)) + HPL30/(2.*(1 - x)*(1 + x)) - (5*HPL30)/(2.*(1 - x)*x*(1 + x)) - &
               &       (5*x*HPL30)/(2.*(1 - x)*(1 + x)) - (x**2*HPL30)/(2.*(1 - x)*(1 + x)) + (5*x**3*HPL30)/(2.*(1 - x)*(1 + x)) - &
               &       (x**4*HPL30)/(2.*(1 - x)*(1 + x)) + (2*HPL31)/(1 + x) + HPL31/(x*(1 + x)) + &
               &       (3*x*HPL31)/(1 + x) + (2*x**2*HPL31)/(1 + x) + (x**3*HPL31)/(1 + x) + (x*Log(1 - x))/24. - &
               &       (701*Log(x))/144. - (149*x*Log(x))/144. - (67*x**2*Log(x))/18. + (HPL2_1*Log(x))/(1 + x) + (HPL2_1*Log(x))/(2.*x*(1 + x)) + &
               &       (3*x*HPL2_1*Log(x))/(2.*(1 + x)) + (x**2*HPL2_1*Log(x))/(1 + x) + (x**3*HPL2_1*Log(x))/(2.*(1 + x)) - &
               &       (HPL20*Log(x))/(2.*(1 - x)*(1 + x)) + (3*HPL20*Log(x))/(2.*(1 - x)*x*(1 + x)) + (3*x*HPL20*Log(x))/(2.*(1 - x)*(1 + x)) + &
               &       (x**2*HPL20*Log(x))/(2.*(1 - x)*(1 + x)) - (3*x**3*HPL20*Log(x))/(2.*(1 - x)*(1 + x)) + &
               &       (x**4*HPL20*Log(x))/(2.*(1 - x)*(1 + x)) - (Log(1 - x)**2*Log(x))/(2.*(1 - x)) + (Log(1 - x)**2*Log(x))/(4.*(1 - x)*x) + &
               &       (3*x*Log(1 - x)**2*Log(x))/(4.*(1 - x)) - (x**2*Log(1 - x)**2*Log(x))/(2.*(1 - x)) + (x**3*Log(1 - x)**2*Log(x))/(4.*(1 - x)) + &
               &       (25*Log(x)**2)/48. - (11*x*Log(x)**2)/48. + (11*x**2*Log(x)**2)/12. - (Log(1 - x)*Log(x)**2)/(2.*(1 - x)) + &
               &       (Log(1 - x)*Log(x)**2)/(4.*(1 - x)*x) + (3*x*Log(1 - x)*Log(x)**2)/(4.*(1 - x)) - (x**2*Log(1 - x)*Log(x)**2)/(2.*(1 - x)) + &
               &       (x**3*Log(1 - x)*Log(x)**2)/(4.*(1 - x)) - Log(x)**3/(12.*(1 - x)*(1 + x)) - (x*Log(x)**3)/(6.*(1 - x)*(1 + x)) + &
               &       (x**2*Log(x)**3)/(12.*(1 - x)*(1 + x)) + (x**3*Log(x)**3)/(6.*(1 - x)*(1 + x)) - (x**4*Log(x)**3)/(12.*(1 - x)*(1 + x)) + &
               &       (Pi**2*Log(1 + x))/(6.*(1 + x)) + (Pi**2*Log(1 + x))/(12.*x*(1 + x)) + (Pi**2*x*Log(1 + x))/(4.*(1 + x)) + &
               &       (Pi**2*x**2*Log(1 + x))/(6.*(1 + x)) + (Pi**2*x**3*Log(1 + x))/(12.*(1 + x)) - (Log(x)**2*Log(1 + x))/(2.*(1 + x)) - &
               &       (Log(x)**2*Log(1 + x))/(4.*x*(1 + x)) - (3*x*Log(x)**2*Log(1 + x))/(4.*(1 + x)) - (x**2*Log(x)**2*Log(1 + x))/(2.*(1 + x)) - &
               &       (x**3*Log(x)**2*Log(1 + x))/(4.*(1 + x)) + (Log(x)*Log(1 + x)**2)/(1 + x) + (Log(x)*Log(1 + x)**2)/(2.*x*(1 + x)) + &
               &       (3*x*Log(x)*Log(1 + x)**2)/(2.*(1 + x)) + (x**2*Log(x)*Log(1 + x)**2)/(1 + x) + (x**3*Log(x)*Log(1 + x)**2)/(2.*(1 + x)) - &
               &       Log(1 + x)**3/(3.*(1 + x)) - Log(1 + x)**3/(6.*x*(1 + x)) - (x*Log(1 + x)**3)/(2.*(1 + x)) - (x**2*Log(1 + x)**3)/(3.*(1 + x)) - &
               &       (x**3*Log(1 + x)**3)/(6.*(1 + x))))
       else
          res = zero
       end if
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       if (one - x > cutoff) then 
          ! subtract the singular part of the 1/(1-x)_+ distribution
          res = res  - 4._dp*((7._dp*CA*nf)/27._dp + CA**2*(-1.8703703703703705_dp + (7._dp*Zeta3)/4._dp))/(one-x)
       else
          res = res
       end if
    case(cc_DELTA)
       res = 4._dp*(-12.405092592592593 - (5*CA)/192. - CF/24. + (9*CF**2)/8. - (137*Lt)/48. - (1679*Pi**2)/192. - (37*Pi**4)/64. + &
            &    CA*CF*(-3.0208333333333335 - (11*Lt)/16. - (7*Pi**2)/16.) + &
            &    CA**2*(5.532986111111111 + (7*Lt)/16. + (43*Pi**2)/36. + (79*Pi**4)/1152. - (55*Zeta3)/36.) + &
            &    CA*nf*(-0.9965277777777778 - (5*Pi**2)/72. - (2*Zeta3)/9.) + CF*nf*(-0.8541666666666666 + Lt/4. + Zeta3/2.) + (499*Zeta3)/48.)
    end select

    if (cc_piece /= cc_DELTA) res = res * x     
  end function C2gg

  
  !----------------------------------------------------------------------
  function C2qg(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    complex(dp) :: HPL2, HPL3
    complex(dp) :: HPL20, HPL21, HPL2_1, HPL30, HPL31, HPL3_1, HPL3x_1px, HPL31_1px
    real(dp), parameter    :: cutoff=1e-11 ! cutoff to regularise the x->1 singularity in the gg coefficient function
    !------------------------------------------
    ! The following declarations are necessary for hplog
    integer, parameter :: n1=-1
    integer, parameter :: n2= 1
    integer, parameter :: nw= 4
    real(dp)    :: Hr1(n1:n2),Hr2(n1:n2,n1:n2),Hr3(n1:n2,n1:n2,n1:n2),&
         & Hr4(n1:n2,n1:n2,n1:n2,n1:n2)
    real(dp)    :: Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2),&
         & Hi4(n1:n2,n1:n2,n1:n2,n1:n2)
    complex(dp) :: Hc1(n1:n2),Hc2(n1:n2,n1:n2),Hc3(n1:n2,n1:n2,n1:n2), &
         & Hc4(n1:n2,n1:n2,n1:n2,n1:n2)
    
    x = exp(-y)
    res = zero

    ! evaluate HPL's using hplog (much faster than Chaplin)
    call hplog(x,nw,Hc1,Hc2,Hc3,Hc4, &
         &     Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2)
    HPL20 = Hc2(0,1)
    HPL30 = Hc3(0,0,1)
    call hplog(-x,nw,Hc1,Hc2,Hc3,Hc4, &
         &     Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2)
    HPL2_1 = Hc2(0,1)
    HPL3_1 = Hc3(0,0,1)
    call hplog(one-x,nw,Hc1,Hc2,Hc3,Hc4, &
         &     Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2)
    HPL21 = Hc2(0,1)
    HPL31 = Hc3(0,0,1)
    call hplog(x/(one+x),nw,Hc1,Hc2,Hc3,Hc4, &
         &     Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
    HPL3x_1px = Hc3(0,0,1)
    call hplog(one/(one+x),nw,Hc1,Hc2,Hc3,Hc4, &
         &     Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
    HPL31_1px = Hc3(0,0,1)


    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       if (one-x > cutoff) then
          res = 4*((CF*(-78 - 24*HPL30 - 24*HPL31 + 450*x + 48*HPL30*x + 48*HPL31*x - 12*Pi**2*x - 432*x**2 - 48*HPL30*x**2 - 48*HPL31*x**2 + 12*Pi**2*x**2 + 192*zeta3 - 384*x*zeta3 + 384*x**2*zeta3 - 36*x*log(1 - x) + 48*x**2*log(1 - x) + 24*HPL21*(1 - 2*x + 2*x**2)*log(1 - x) + 24*x*log(1 - x)**2 - 24*x**2*log(1 - x)**2 - 4*log(1 - x)**3 + 8*x*log(1 - x)**3 - 8*x**2*log(1 - x)**3 + 48*log(x) + 90*x*log(x) - 48*x**2*log(x) + 24*HPL20*(1 - 2*x + 2*x**2)*log(x) - 48*x*log(1 - x)*log(x) + 48*x**2*log(1 - x)*log(x) + 12*log(1 - x)**2*log(x) - 24*x*log(1 - x)**2*log(x) + 24*x**2*log(1 - x)**2*log(x) + 3*log(x)**2 + 36*x*log(x)**2 - 24*x**2*log(x)**2 + 12*log(1 - x)*log(x)**2 - 24*x*log(1 - x)*log(x)**2 + 24*x**2*log(1 - x)*log(x)**2 - 2*log(x)**3 + 4*x*log(x)**3 - 8*x**2*log(x)**3))/192 + (CA*(344 - 630*x + 108*HPL31*x + 216*HPL31_1px*x + 324*HPL3_1*x + 774*x**2 + 864*HPL30*x**2 - 216*HPL31*x**2 + 432*HPL31_1px*x**2 + 648*HPL3_1*x**2 + 36*Pi**2*x**2 - 596*x**3 + 216*HPL31*x**3 + 432*HPL31_1px*x**3 + 648*HPL3_1*x**3 - 324*x*zeta3 - 648*x**3*zeta3 + 162*x**2*log(1 - x) - 216*x**3*log(1 - x) - 108*x**2*log(1 - x)**2 + 108*x**3*log(1 - x)**2 + 18*x*log(1 - x)**3 - 36*x**2*log(1 - x)**3 + 36*x**3*log(1 - x)**3 + 252*x*log(x) - 360*x**2*log(x) - 864*HPL20*x**2*log(x) + 72*Pi**2*x**2*log(x) + 816*x**3*log(x) - 27*x*log(x)**2 + 108*x**2*log(x)**2 - 396*x**3*log(x)**2 - 432*x**2*log(1 - x)*log(x)**2 + 18*x*log(x)**3 + 36*x**2*log(x)**3 - 36*HPL21*(4 - 6*x + 24*x**2 - 22*x**3 + 3*x*(1 - 2*x + 2*x**2)*log(1 - x) + 12*x**2*log(x)) - 108*HPL2_1*x*(-2*x*(1 + x) + (1 + 2*x + 2*x**2)*log(x)) + 18*Pi**2*x*log(1 + x) + 36*Pi**2*x**2*log(1 + x) + 36*Pi**2*x**3*log(1 + x) + 216*x**2*log(x)*log(1 + x) + 216*x**3*log(x)*log(1 + x) + 54*x*log(x)**2*log(1 + x) + 108*x**2*log(x)**2*log(1 + x) + 108*x**3*log(x)**2*log(1 + x) - 36*x*log(1 + x)**3 - 72*x**2*log(1 + x)**3 - 72*x**3*log(1 + x)**3))/(864*x))
       else
          res = zero
       end if
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       ! no virtual piece
    case(cc_DELTA)
       res = zero
    end select

    if (cc_piece /= cc_DELTA) res = res * x     
  end function C2qg
  
  !======================================================================
  !! Initialise a coefficient function matrix, repeating exactly
  !! what's done for a LO unpolarised splitting matrix, with the nf
  !! value that is current from the qcd module.
  subroutine InitCoeffMatrix(grid, C1, C2, G1)
    type(grid_def),  intent(in)    :: grid
    type(split_mat), intent(inout) :: C1
    type(split_mat), optional, intent(inout) :: C2,G1
    type(grid_conv) :: C2qqS_grid, C2qqbar_grid
    
    write(*,*) 'Initialising coefficient function grids ...'

    !-- O(as) coefficient functions
    !C1%loops  = 1
    C1%nf_int = nf_int
    
    call cobj_InitSplitLinks(C1)

    call InitGridConv(grid, C1%gg, C1gg)
    call InitGridConv(grid, C1%qq, C1qq)
    call InitGridConv(grid, C1%gq, C1gq)
    call InitGridConv(grid, C1%qg, C1qg)

    !-- now fix up pieces so that they can be directly used as a matrix
    call Multiply(C1%qg, 2*nf)

    !-- PqqV +- PqqbarV
    call InitGridConv(C1%NS_plus,  C1%qq)
    call InitGridConv(C1%NS_minus, C1%qq)

    !-- PNSminus + nf * (PqqS - PqqbarS)
    call InitGridConv(C1%NS_V, C1%NS_minus)


    !-- O(as^2) coefficient functions. The C2qq is non-trivial at this order and we need
    !-- to decompose it consistently in the Hoppet's evolution basis (singlet, non-singlet, valence)
    !C2%loops  = 1
    ! fix the number of active flavours
    C2%nf_int = nf_int

    call cobj_InitSplitLinks(C2)

    !-- Interpret C2qq and C2qqbar as non-singlet components, and
    !-- C2qq', C2qqbar' as singlet components. At as^2 C2qq'=C2qqbar' as
    !-- for the singlet splitting functions. Therefore it seems that all
    !-- the symmetries are preserved. This follows Hoppet's implementation
    !-- of the NLO splitting functions

    !-- start with the +- components of the non-singlet
    !-- PqqV +- PqqbarV, remember that C2qq'=C2qqbar' at as^2
    ! plus component
    call InitGridConv(grid, C2%NS_plus,  C2qq)
    call InitGridConv(grid, C2qqbar_grid,  C2qqbar)
    call AddWithCoeff(C2%NS_plus,  C2qqbar_grid, one)
    call InitGridConv(grid, C2qqS_grid,  C2qqp)
    call AddWithCoeff(C2%NS_plus,  C2qqS_grid, -two)
    ! minus component
    call InitGridConv(grid, C2%NS_minus, C2qq)
    call AddWithCoeff(C2%NS_minus, C2qqbar_grid, -one)

    !-- now the valence component
    !-- PNSminus + nf * (PqqS - PqqbarS), remember that C2qq'=C2qqbar' at as^2
    call InitGridConv(C2%NS_V, C2%NS_minus)

    !-- finally, build the singlet component (eq. 2.4 of Moch et al. 0404111)
    !-- Pqq = PNS_plus + nf*(PqqS + PqqbarS), remember that C2qq'=C2qqbar' at as^2
    call InitGridConv(C2%qq, C2%NS_plus)
    call AddWithCoeff(C2%qq, C2qqS_grid, two*nf)

    ! Clean up
    call Delete(C2qqS_grid)
    call Delete(C2qqbar_grid)

    !-- now build the rest of the singlet
    call InitGridConv(grid, C2%gg, C2gg)
    call InitGridConv(grid, C2%gq, C2gq)
    call InitGridConv(grid, C2%qg, C2qg)
    !-- add the factor of 2*nf in Cqg to sum over all species
    call Multiply(C2%qg, 2*nf)


    !-- O(as) coefficient functions describing spin correlation
    !G1%loops  = 1
    G1%nf_int = nf_int

    call cobj_InitSplitLinks(G1)

    call InitGridConv(grid, G1%gg, G1gg)
    call InitGridConv(grid, G1%qq, G1qq)
    call InitGridConv(grid, G1%gq, G1gq)
    call InitGridConv(grid, G1%qg, G1qg)

    !-- now fix up pieces so that they can be directly used as a matrix
    call Multiply(G1%qg, 2*nf)

    !-- PqqV +- PqqbarV
    call InitGridConv(G1%NS_plus,  G1%qq)
    call InitGridConv(G1%NS_minus, G1%qq)

    !-- PNSminus + nf * (PqqS - PqqbarS)
    call InitGridConv(G1%NS_V, G1%NS_minus)
  end subroutine InitCoeffMatrix

end module coefficient_functions_VH
