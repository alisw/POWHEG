      subroutine setreal(p,rflav,amp2)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_st.h'
      include 'pwhg_math.h'
      real * 8 p(0:3,nlegreal)
      integer rflav(nlegreal)
      real * 8 amp2
      double precision pstd(0:3,nlegreal), fermfac
      integer stdflav(nlegreal), index(nlegreal)

c     compute standard ordering of flavours and momenta
      call stdorderedflav(p, rflav, nlegreal, pstd, stdflav, index,
     &     fermfac)
c     get real via crossing
      call stdreal(pstd, stdflav, index, fermfac, amp2)

      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccccc
ccc Routine to compute real MEs through        ccc
ccc crossing, where the crossing is determined ccc
ccc by an array index(nleg) specifying the     ccc
ccc permutation. The momenta and flavours have ccc
ccc to be given in standard ordering.          ccc
cccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine stdreal(p, flav, index, fermfac, amp)

      implicit none
      include 'nlegborn.h'
      include 'pwhg_math.h'
      double precision p(0:3, nlegreal), fermfac, amp
      integer flav(nlegreal), index(nlegreal)
      double precision preal(0:3, nlegreal)
      double precision Q, s12, sab, ta1, ta2, tb1, prefac
      double precision qqFac, qgFac, ggFac, qqAqq, qqbAgg, qqpAqqp,
     &     udAud
      logical isuptype
      integer i

c     initialize to zero
      amp = 0.d0

      if (flav(nlegreal).EQ.22) then
         if (flav(1).EQ.flav(2)) then ! -> qb qb q q A
c     reorder to q q -> A q q (no fermfac needed for two fermion flips)
            do i = 0, 3
               preal(i, 1) = -p(i, 1)
               preal(i, 2) = -p(i, 2)
               preal(i, 3) = p(i, 5)
               preal(i, 4) = p(i, 4)
               preal(i, 5) = p(i, 3)
            enddo
            call mandelstam23(preal, s12, sab, ta1, ta2, tb1)
            call setcharge(flav(3), Q)
c     find out which particles are initial, i.e. which have index = 1, 2
            if ((index(1)*index(2).EQ.2).OR.(index(3)*index(4).EQ.2))
     &           then
               prefac = 0.5d0 * qqFac(NC) ! final state symmetry
            else
               prefac = qqFac(NC)
            endif
            amp = fermfac * qqAqq(Q, s12, sab, ta1, ta2, tb1) * prefac

         elseif (flav(1).EQ.-flav(2)) then ! -> qb q g g A
            do i = 0, 3         ! reorder to q qb -> A g g
               preal(i, 1) = -p(i, 1)
               preal(i, 2) = -p(i, 2)
               preal(i, 3) = p(i, 5)
               preal(i, 4) = p(i, 4)
               preal(i, 5) = p(i, 3)
            enddo
            call mandelstam23(preal, s12, sab, ta1, ta2, tb1)
            call setcharge(flav(2), Q)
            if (index(1)*index(2).EQ.2) then
               prefac = qqFac(NC) * 0.5d0 ! final state symmetry
            elseif ((index(1)*index(3).EQ.2).OR.
     &              (index(1)*index(4).EQ.2).OR.
     &              (index(2)*index(3).EQ.2).OR.
     &              (index(2)*index(4).EQ.2)) then
               prefac = qgFac(NC)
            elseif ((index(3)*index(4).EQ.2)) then
               prefac = ggFac(NC)
            else
               stop 'stdreal: qbqggA: illegal initial particles'
            endif
            amp = fermfac * qqbAgg(Q, s12, sab, ta1, ta2, tb1) * prefac

         elseif (MOD(flav(1)-flav(2), 2).EQ.0) then ! -> qb qb' q' q A (same charge)
            do i = 0, 3         ! reorder to q q' -> A q q'
               preal(i, 1) = -p(i, 1)
               preal(i, 2) = -p(i, 2)
               preal(i, 3) = p(i, 5)
               preal(i, 4) = p(i, 4)
               preal(i, 5) = p(i, 3)
            enddo
            call mandelstam23(preal, s12, sab, ta1, ta2, tb1)
            call setcharge(flav(3), Q)
            amp = fermfac * qqpAqqp(Q, s12, sab, ta1, ta2, tb1) *
     &           qqFac(NC)

         elseif (isuptype(flav(1))) then ! -> ub db d u A
            do i = 0, 3         ! reorder to u d -> A u d
               preal(i, 1) = -p(i, 1)
               preal(i, 2) = -p(i, 2)
               preal(i, 3) = p(i, 5)
               preal(i, 4) = p(i, 4)
               preal(i, 5) = p(i, 3)
            enddo
            call mandelstam23(preal, s12, sab, ta1, ta2, tb1)
            amp = fermfac * udAud(s12, sab, ta1, ta2, tb1) * qqFac(NC)

         elseif (isuptype(flav(2))) then ! -> db ub u d A
            do i = 0, 3         ! reorder to u d -> A u d
               preal(i, 1) = -p(i, 2)
               preal(i, 2) = -p(i, 1)
               preal(i, 3) = p(i, 5)
               preal(i, 4) = p(i, 3)
               preal(i, 5) = p(i, 4)
            enddo
            call mandelstam23(preal, s12, sab, ta1, ta2, tb1)
            amp = fermfac * udAud(s12, sab, ta1, ta2, tb1) * qqFac(NC)
         endif

      else
         stop 'stdreal: Pure QCD not implemented!'
      endif

      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccccc
ccc Compute Mandelstam vars for 2 -> 3 process ccc
cccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine mandelstam23(p, s12, sab, ta1, ta2, tb1)

      implicit none
      include 'nlegborn.h'
      double precision p(0:3, nlegreal), s12, sab, ta1, ta2, tb1
      double precision dotp

      sab = 2.d0 * dotp(p(0,1),p(0,2))
      s12 = 2.d0 * dotp(p(0,3),p(0,4))
      ta1 = -2.d0 * dotp(p(0,1),p(0,3))
      ta2 = -2.d0 * dotp(p(0,1),p(0,4))
      tb1 = -2.d0 * dotp(p(0,2),p(0,3))

      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccccc
ccc Determine if particle is up-type.          ccc
cccccccccccccccccccccccccccccccccccccccccccccccccc

      logical function isuptype(flav)

      implicit none
      integer flav

      if ((flav.NE.0).AND.(MOD(flav,2).EQ.0)) then
         isuptype = .TRUE.
      else
         isuptype = .FALSE.
      endif

      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccccc
ccc The MEs have been obtained via FormCalc7.5 ccc
ccc Since in 2 -> 3 gamma production           ccc
ccc up-type and down-type photon emission      ccc
ccc interfere and FormCalc doesn't allow one   ccc
ccc to use charges as variables, one has to    ccc
ccc separate processes which involve up- and   ccc
ccc down-type quarks from processes where only ccc
ccc one charge is present.                     ccc
ccc That is one has three cases:               ccc
ccc - all quarks equal                         ccc
ccc - different quarks but same charge (type)  ccc
ccc - different quarks with different charge   ccc
cccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccc
ccc q q -> gamma q q
ccc Charge is given as argument!
ccc Symmetry factor 1/2 for identical particles
ccc in the final state has to be applied to this!

      double precision function
     &     qqAqq(Q, s12, sab, ta1, ta2, tb1)
      implicit none
      include 'pwhg_math.h'
      include 'pwhg_st.h'
      include 'pwhg_em.h'
      double precision Q, s12, sab, ta1, ta2, tb1

      qqAqq = (-4608*Q**2*em_alpha*st_alpha*(-1 + NC)*(1 + NC)*Pi**4*
     &    (sab*ta1**2 + s12**2*(sab + 2*ta1) +
     &      s12*(2*sab*ta1 + 3*ta1**2 + 2*ta1*ta2 +
     &         ta1*tb1 - 2*ta2*tb1) +
     &      (ta1 + tb1)*(ta1**2 + ta1*ta2 - ta2*tb1))*
     &    ((2*sab**2 + 2*sab*(ta1 + tb1) + (ta1 + tb1)**2)*
     &       (sab*(ta1 + 2*ta2) + s12*(sab + ta1 + 2*ta2) +
     &         (ta1 + ta2)*(ta1 + 2*ta2 + tb1)) +
     &      NC*(4*sab**4 + s12**3*(sab + ta1 + 2*ta2) +
     &         2*sab**3*(7*ta1 + 6*ta2 + 4*tb1) +
     &         sab**2*(19*ta1**2 + 34*ta1*ta2 + 16*ta2**2 +
     &            22*ta1*tb1 + 18*ta2*tb1 + 6*tb1**2) +
     &         3*s12**2*(sab**2 +
     &            (ta1 + ta2)*(ta1 + 2*ta2 + tb1) +
     &            sab*(2*ta1 + 2*ta2 + tb1)) +
     &         (ta1 + ta2)*
     &          (3*ta1**3 + 4*ta2**3 + 4*ta2**2*tb1 +
     &            5*ta2*tb1**2 + 2*tb1**3 +
     &            ta1**2*(9*ta2 + 7*tb1) +
     &            2*ta1*(4*ta2**2 + 6*ta2*tb1 + 3*tb1**2))
     &          + sab*(12*ta1**3 +
     &            2*(ta2 + tb1)*(2*ta2 + tb1)**2 +
     &            ta1**2*(32*ta2 + 21*tb1) +
     &            4*ta1*(7*ta2**2 + 9*ta2*tb1 + 3*tb1**2))
     &          + s12*(6*sab**3 + 5*ta1**3 +
     &            8*sab**2*(2*ta1 + 2*ta2 + tb1) +
     &            ta1**2*(17*ta2 + 8*tb1) +
     &            2*ta1*(9*ta2**2 + 8*ta2*tb1 + 2*tb1**2) +
     &            ta2*(8*ta2**2 + 6*ta2*tb1 + 5*tb1**2) +
     &            sab*(15*ta1**2 + 4*ta1*(7*ta2 + 4*tb1) +
     &               4*(3*ta2**2 + 4*ta2*tb1 + tb1**2))))))/
     &  (9.d0*NC*s12*ta1*ta2*(s12 + ta1 + ta2)*
     &    (sab + ta1 + ta2)*tb1*(s12 + ta1 + tb1)*
     &    (s12 + sab + ta1 + ta2 + tb1))

      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccc
ccc q qbar -> gamma g g
ccc Charge is given as argument!
ccc Symmetry factor 1/2 for identical particles
ccc in the final state has to be applied to this!

      double precision function
     &     qqbAgg(Q, s12, sab, ta1, ta2, tb1)
      implicit none
      include 'pwhg_math.h'
      include 'pwhg_st.h'
      include 'pwhg_em.h'
      double precision Q, s12, sab, ta1, ta2, tb1

      qqbAgg =(-4608*Q**2*em_alpha*st_alpha*(-1 + NC)*(1 + NC)*Pi**4*
     &    (-(sab*(sab + ta1 + tb1)) +
     &      NC**2*(sab**2 + s12*(sab + ta1 + 2*ta2) +
     &         (ta1 + ta2)*(ta1 + 2*ta2 + tb1) +
     &         sab*(2*ta1 + 2*ta2 + tb1)))*
     &    (2*ta1**4 + 9*ta1**3*ta2 + 15*ta1**2*ta2**2 +
     &      12*ta1*ta2**3 + 4*ta2**4 +
     &      sab**3*(ta1 + 2*ta2) +
     &      s12**3*(sab + ta1 + 2*ta2) - ta1**3*tb1 +
     &      3*ta1**2*ta2*tb1 + 6*ta1*ta2**2*tb1 +
     &      4*ta2**3*tb1 + 3*ta1*ta2*tb1**2 +
     &      3*ta2**2*tb1**2 - ta1*tb1**3 + ta2*tb1**3 +
     &      3*sab**2*(ta1**2 + 3*ta1*ta2 +
     &         ta2*(2*ta2 + tb1)) +
     &      3*s12**2*(ta1**2 + 3*ta1*ta2 +
     &         sab*(ta1 + 2*ta2) + ta2*(2*ta2 + tb1)) +
     &      sab*(4*ta1**3 + 15*ta1**2*ta2 +
     &         6*ta1*ta2*(3*ta2 + tb1) +
     &         ta2*(8*ta2**2 + 6*ta2*tb1 + 3*tb1**2)) +
     &      s12*(sab**3 + 4*ta1**3 + 15*ta1**2*ta2 +
     &         3*sab**2*(ta1 + 2*ta2) +
     &         6*ta1*ta2*(3*ta2 + tb1) +
     &         ta2*(8*ta2**2 + 6*ta2*tb1 + 3*tb1**2) +
     &         6*sab*(ta1**2 + 3*ta1*ta2 +
     &            ta2*(2*ta2 + tb1)))))/
     &  (9.d0*NC*ta1*ta2*(s12 + ta1 + ta2)*
     &    (sab + ta1 + ta2)*tb1*(sab + ta1 + tb1)*
     &    (s12 + sab + ta1 + ta2 + tb1))

      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccc
ccc u d -> gamma u d
ccc up-down-type mixing channel
ccc Charge is explicit in ME.

      double precision function
     &     udAud(s12, sab, ta1, ta2, tb1)
      implicit none
      include 'pwhg_math.h'
      include 'pwhg_st.h'
      include 'pwhg_em.h'
      double precision s12, sab, ta1, ta2, tb1

      udAud = (512*em_alpha*st_alpha*(-1 + NC)*(1 + NC)*Pi**4*
     &    (s12**4*(2*sab + ta1) +
     &      s12**3*(4*sab**2 + 5*ta1**2 + 3*ta1*ta2 +
     &         4*ta1*tb1 + 2*ta2*tb1 +
     &         2*sab*(5*ta1 + 2*(ta2 + tb1))) +
     &      s12**2*(8*sab**3 + 11*ta1**3 +
     &         8*sab**2*(3*ta1 + ta2 + tb1) +
     &         4*ta2*tb1*(ta2 + 2*tb1) +
     &         2*ta1**2*(7*ta2 + 8*tb1) +
     &         2*ta1*(2*ta2**2 + 9*ta2*tb1 + 3*tb1**2) +
     &         sab*(28*ta1**2 + 4*(ta2 + tb1)**2 +
     &            ta1*(22*ta2 + 24*tb1))) +
     &      2*(4*sab**3*ta1**2 +
     &         2*sab**2*(5*ta1**3 + 6*ta1*ta2*tb1 +
     &            4*ta2*tb1**2 + 4*ta1**2*(ta2 + tb1)) +
     &         (ta1 + tb1)*
     &          (3*ta1**4 + ta1**3*(7*ta2 + 4*tb1) +
     &            4*ta2*tb1*(ta2**2 + ta2*tb1 + tb1**2) +
     &            2*ta1**2*
     &             (3*ta2**2 + 6*ta2*tb1 + tb1**2) +
     &            2*ta1*ta2*(ta2**2 + 5*ta2*tb1 + 5*tb1**2))
     &           + sab*(9*ta1**4 + 14*ta1**3*(ta2 + tb1) +
     &            8*ta2*tb1**2*(ta2 + tb1) +
     &            12*ta1*ta2*tb1*(ta2 + 2*tb1) +
     &            ta1**2*(6*ta2**2 + 28*ta2*tb1 + 6*tb1**2))
     &         ) + s12*(16*sab**3*ta1 + 13*ta1**4 +
     &         ta1**3*(23*ta2 + 26*tb1) +
     &         4*sab**2*(10*ta1**2 + 5*ta1*ta2 +
     &            6*ta1*tb1 + 2*ta2*tb1) +
     &         4*ta2*tb1*(ta2**2 + 3*ta2*tb1 + 3*tb1**2) +
     &         2*ta1**2*
     &          (7*ta2**2 + 22*ta2*tb1 + 9*tb1**2) +
     &         2*ta1*(ta2**3 + 13*ta2**2*tb1 +
     &            17*ta2*tb1**2 + 2*tb1**3) +
     &         2*sab*(19*ta1**3 + 4*ta2*tb1*(ta2 + 2*tb1) +
     &            3*ta1**2*(7*ta2 + 8*tb1) +
     &            ta1*(6*ta2**2 + 22*ta2*tb1 + 8*tb1**2)))))
     &   /(9.d0*s12*ta1*ta2*(s12 + ta1 + ta2)*tb1*
     &    (s12 + ta1 + tb1))

      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccc
ccc q q' -> gamma q q'
ccc with q,q' of same charge
ccc Charge is given as argument!

      double precision function
     &     qqpAqqp(Q, s12, sab, ta1, ta2, tb1)
      implicit none
      include 'pwhg_math.h'
      include 'pwhg_st.h'
      include 'pwhg_em.h'
      double precision Q, s12, sab, ta1, ta2, tb1

      qqpAqqp = (-4608*Q**2*em_alpha*st_alpha*(-1 + NC)*(1 + NC)
     &     * Pi**4*
     &    (4*sab**3*ta1**2 + s12**4*(sab + 2*ta1) +
     &      2*sab**2*(5*ta1**3 - 2*ta2*tb1**2 +
     &         4*ta1**2*(ta2 + tb1)) +
     &      s12**3*(2*sab**2 + 7*ta1**2 + 6*ta1*ta2 +
     &         5*ta1*tb1 - 2*ta2*tb1 +
     &         2*sab*(4*ta1 + ta2 + tb1)) +
     &      (ta1 + tb1)*(3*ta1**4 +
     &         ta1**3*(7*ta2 + 4*tb1) +
     &         2*ta1*ta2*(ta2**2 - ta2*tb1 - tb1**2) -
     &         2*ta2*tb1*(ta2**2 + ta2*tb1 + tb1**2) +
     &         ta1**2*(6*ta2**2 + 3*ta2*tb1 + 2*tb1**2)) +
     &      sab*(9*ta1**4 - 6*ta1*ta2*tb1**2 +
     &         14*ta1**3*(ta2 + tb1) -
     &         4*ta2*tb1**2*(ta2 + tb1) +
     &         2*ta1**2*(3*ta2**2 + 5*ta2*tb1 + 3*tb1**2))
     &       + s12**2*(4*sab**3 + 13*ta1**3 -
     &         ta2*tb1*(4*ta2 + 5*tb1) +
     &         ta1**2*(19*ta2 + 17*tb1) +
     &         2*sab*(13*ta1**2 + 10*ta1*ta2 + ta2**2 +
     &            9*ta1*tb1 - ta2*tb1 + tb1**2) +
     &         ta1*(8*ta2**2 + 6*ta2*tb1 + 6*tb1**2) +
     &         2*sab**2*(9*ta1 + 2*(ta2 + tb1))) +
     &      s12*(8*sab**3*ta1 + 11*ta1**4 +
     &         ta1**3*(22*ta2 + 19*tb1) +
     &         2*sab**2*(13*ta1**2 + 8*ta1*ta2 +
     &            6*ta1*tb1 - 4*ta2*tb1) -
     &         2*ta2*tb1*
     &          (2*ta2**2 + 3*ta2*tb1 + 3*tb1**2) +
     &         4*ta1**2*(4*ta2**2 + 4*ta2*tb1 + 3*tb1**2) +
     &         2*ta1*(2*ta2**3 - ta2**2*tb1 -
     &            2*ta2*tb1**2 + tb1**3) +
     &         2*sab*(14*ta1**3 - ta2*tb1*(4*ta2 + 5*tb1) +
     &            3*ta1**2*(6*ta2 + 5*tb1) +
     &            2*ta1*(3*ta2**2 + ta2*tb1 + 2*tb1**2)))))/
     &  (9.d0*s12*ta1*ta2*(s12 + ta1 + ta2)*tb1*
     &    (s12 + ta1 + tb1))

      return
      end
