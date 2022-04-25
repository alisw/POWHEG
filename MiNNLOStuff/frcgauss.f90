module frcgauss_intrfc

  interface
     recursive function frcgauss(fun, a, b, eps)
       implicit none
       !--- arguments -------------------------------------------------
       real(kind(1d0)) :: frcgauss
       interface
          function fun(x)
            real(kind(1d0))              :: fun
            real(kind(1d0)), intent(in)  :: x
          end function fun
       end interface
       real(kind(1d0)), intent(in) :: a, b, eps
     end function frcgauss
  end interface
end module frcgauss_intrfc


!----------------------------------------------------------------------
! Routine for adaptive gaussian integration. The algorithm is taken
! straight from  the description of the CERNLIB routine GAUSS
! (D103). The only differences are that (a) it is recursive, and (b)
! it calls a subroutine which evaluates the function at an array of
! points -- x(:) (saves on multiple sub function calls -- important if
! the function being evaluated is 'light')
!
! Integration limits are [a,b], relative/absolute (better of the two)
! accuracy eps (for details see CERNLIB writeup D103).
!
! GPS 15/11/96 (CCN9 84)
!----------------------------------------------------------------------
recursive function frcgauss(fun, a, b, eps)
  implicit none
  !--- arguments ------------------------------------------------------
  real(kind(1d0)) :: frcgauss
  interface
     function fun(x)
       real(kind(1d0))              :: fun
       real(kind(1d0)), intent(in)  :: x
     end function fun
  end interface
  real(kind(1d0)), intent(in) :: a, b, eps
  !--- internal vars --------------------------------------------------
  integer, parameter :: nl = 8, nh = 16 ! lower and higher # of absicsae
  real(kind(1d0))    :: gl, fl(nl), xl(nl), wl(nl)
  real(kind(1d0))    :: gh, fh(nh), xh(nh), wh(nh)
  real(kind(1d0))    :: c0, c1, r, diffab, diff01
  integer i

  frcgauss = 0d0
!!  write(*,*) a, b
  diffab = b - a

  !-- first interval pair
  c0 = a
  c1 = b

  !-- calculations of absiscae could probably be optimised if
  !   incorporated into remainder of the structure here -- but that
  !   would add an extra level of complication...
  do 
     diff01 = c1 - c0
     !-- first check to see if interval being deal with is acceptable.
     if ((1d0 + diff01*0.005d0/diffab) == 1d0) then
        write(*,*) 'frcgauss: desired accuracy could not be achieved'
        write(*,*) 'stop'
        stop
     end if
     !-- calculate 8 and 16 point quadratures for this interval
     call dgset(c0,c1,nl,xl,wl)
     call dgset(c0,c1,nh,xh,wh)
     gl = 0d0; gh = 0d0
     !--- not the most efficient approach?
     do i = 1, nl
        gl = gl + fun(xl(i))*wl(i)
     end do
     do i = 1, nh
        gh = gh + fun(xh(i))*wh(i)
     end do

     !-- test to see if accuracy has been reached
     r = abs(gh - gl)/(1d0 + abs(gh))
!     write(*,*) c0, c1, r
     if (r > eps) then           ! if not then subdivide interval
        c1 = c0 + 0.5 * diff01   
     elseif (c1 /= b) then       ! else add this interval, move to next
        frcgauss = frcgauss + gh   
        c0 = c1 ; c1 = b         
     else                        ! or add this interval and finish
        frcgauss = frcgauss + gh
        exit
     end if
  end do
end function frcgauss

