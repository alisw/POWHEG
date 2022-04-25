c     Wrapper function to call setborneps with eps=0
      subroutine setborn(p,bflav,born,bornjk,bmunu)

      implicit none
      include 'nlegborn.h'
      integer nlegs
      parameter (nlegs=nlegborn)
      real * 8 p(0:3,nlegs),bornjk(nlegs,nlegs)
      integer bflav(nlegs)
      real * 8 bmunu(0:3,0:3,nlegs),born

      call setborneps(p,bflav,0,born,bornjk,bmunu)

      return
      end


c     Compute Born coefficient of epsilon^eps in D = 4 - 2*epsilon
c     dimensions. 
      subroutine setborneps(p,bflav,eps,born,bornjk,bmunu)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_math.h'
      include 'pwhg_kn.h'
      integer nlegs
      parameter (nlegs=nlegborn)
      real * 8 p(0:3,nlegs),bornjk(nlegs,nlegs)
      integer bflav(nlegs), eps
      real * 8 bmunu(0:3,0:3,nlegs),born
      double precision pstd(0:3,nlegs), fermfac
      integer stdflav(nlegs), index(nlegs)
      integer j,k,mu,nu
      real * 8 dotp
      double precision colcor, eta(0:3), gmunu

c     compute standard ordering of flavours and momenta
      call stdorderedflav(p, bflav, nlegs, pstd, stdflav, index,
     &     fermfac)
c     get born via crossing
      call stdborn(pstd, stdflav, index, fermfac, eps, born)

c     Born with final state photon
      if (bflav(3).EQ.22) then
c     Colour correlation (factorizes for 3 coloured particles)
         do j = 1, nlegborn
            do k = 1, nlegborn
               if ((j.NE.3).AND.(k.NE.3)) then
                  bornjk(j,k) = - born * colcor(j,k,bflav(j),bflav(k))
               else
                  bornjk(j,k) = 0.d0
               endif
            enddo
         enddo

c     Spin correlation (turns out this factorizes for q qbar -> gamma g and crossings)
         do j = 1, nlegborn
            if (bflav(j).EQ.0) then
               call etaCnst(p(0,j),eta)
               do mu = 0, 3
                  do nu = 0, 3
                     bmunu(mu,nu,j) = born * 0.5d0 * (-gmunu(mu,nu)
     &                    + (p(mu,j) * eta(nu) + p(nu,j) * eta(mu))
     &                    /dotp(p(0,j),eta))
                  enddo
               enddo
            else
               do mu = 0, 3
                  do nu = 0, 3
                     bmunu(mu,nu,j) = 0.d0
                  enddo
               enddo
            endif
         enddo

      else
c     Colour- and spin-correlations shouldn't be needed,
c     since we don't want to radiate partons off of pure QCD Borns
         bornjk = 0.d0
         bmunu = 0.d0
      endif

      end


cccccccccccccccccccccccccccccccccccccccccccccccc
ccc Put the flavourstructure in a standard   ccc
ccc ascending order (PDG numbering scheme    ccc
ccc i.e. gluon = 21) with all particles      ccc
ccc outgoing. Collect factors of -1 for the  ccc
ccc reversing of fermionlines in fermfac.    ccc
ccc p and flav are the nleg-entries-arrays   ccc
ccc of momenta and flavour structure that    ccc
ccc are to be reordered in pstd and stdflav; ccc
ccc the permutation is saved in index.       ccc
cccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine stdorderedflav(p, flav, nleg, pstd, stdflav, index,
     &     fermfac)

      implicit none
      integer nleg
      double precision p(0:3,nleg)
      integer flav(nleg)
      integer outflav(nleg), stdflav(nleg), index(nleg),
     &     i, j
      double precision pout(0:3,nleg), pstd(0:3,nleg), fermfac

c     factor introduced for reversing fermion lines
      fermfac = 1.d0

      outflav = flav
c     and turn all particles outgoing
      do i = 1, nleg
c     set all gluons from 0 to 21 (as in PDG scheme)
         if (flav(i).EQ.0) then
            outflav(i) = 21
         endif
c     particle/antiparticle flip
         if (i.LT.3) then
            if ((outflav(i).GT.-19).AND.(outflav(i).LT.19)) then
               outflav(i) = -outflav(i)
               fermfac = fermfac * (-1.d0)
            endif
            do j = 0, 3
               pout(j, i) = -p(j, i)
            enddo
         else
            do j = 0, 3
               pout(j, i) = p(j, i)
            enddo
         endif
      enddo

c     put into ascending standard order
      call sortzv(outflav, index, nleg, -1, 0, 0)
      do i = 1, nleg
         stdflav(i) = outflav(index(i))
         do j = 0, 3
            pstd(j, i) = pout(j, index(i))
         enddo
      enddo

      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccc
ccc For standard ordered momentum (p) and    ccc
ccc flavour (flav) arrays return the         ccc
ccc corresponding Born amplitude.            ccc
ccc Select coefficient of epsilon^eps via    ccc
ccc eps, e.g. eps=0 gives the finite normal  ccc
ccc finite amplitude, eps=1 the coefficient  ccc
ccc of eps^1 in D=4-2eps dimensions.         ccc
cccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine stdborn(p, flav, index, fermfac, eps, born)

      implicit none
      include 'nlegborn.h'
      include 'pwhg_math.h'
      double precision p(0:3, nlegborn), fermfac, born
      integer flav(nlegborn), index(nlegborn), eps
      double precision pborn(0:3, nlegborn)
      double precision Q, S, T, prefac
      integer i
      double precision qqFac, qgFac, ggFac, qqbAg, gggg, qqqq, qqbgg,
     &     qqpqqp

      if ((eps.LT.0).OR.(eps.GT.2)) then
         stop 'stdborn: value of eps not implemented'
      endif

c     initialize to zero
      born = 0.d0
      
      if (flav(nlegborn).EQ.22) then ! -> qb q g A
c     reorder to q qb -> A g (no fermfac needed for two fermion flips)
         do i = 0, 3
            pborn(i, 1) = -p(i, 1)
            pborn(i, 2) = -p(i, 2)
            pborn(i, 3) = p(i, 4)
            pborn(i, 4) = p(i, 3)
         enddo
         call mandelstam22(pborn, S, T)
         call setcharge(flav(2), Q)
c     find out which particles are initial, i.e. which have index = 1, 2
         if (index(1)*index(2).EQ.2) then
            prefac = qqFac(NC)
         elseif ((index(1)*index(3).EQ.2).OR.(index(2)*index(3).EQ.2))
     &           then
            prefac = qgFac(NC)
c$$$            if (eps.EQ.1) then  ! construct coefficient of epsilon
c$$$               born = fermfac * qqbAg(Q, S, T, 0) * prefac
c$$$            elseif (eps.EQ.2) then
c$$$               born = fermfac * (qqbAg(Q, S, T, 0) + qqbAg(Q, S, T, 1))
c$$$     &              * prefac
c$$$            endif
c     Checking virtuals against MadLoop shows that the above is actually
c     not to be taken into account. No idea why.
         else
            stop 'stdborn: illegal initial particles'
         endif
         born = fermfac * qqbAg(Q, S, T, eps) * prefac + born

      elseif (flav(1).EQ.21) then ! -> g g g g
         do i = 0, 3            ! reorder to g g -> g g
            pborn(i, 1) = -p(i, 1)
            pborn(i, 2) = -p(i, 2)
            pborn(i, 3) = p(i, 3)
            pborn(i, 4) = p(i, 4)
         enddo
         call mandelstam22(pborn, S, T)
         prefac = ggFac(NC) * 0.5d0 ! final state symmetry
         born = fermfac * gggg(S, T) * prefac
c$$$         if (eps.EQ.1) then
c$$$            born = born(eps=1) + fermfac * gggg(S, T, eps=0) * prefac
c$$$         elseif (eps.EQ.2) then
c$$$            born = born(eps=2) + fermfac * gggg(S, T, eps=0) * prefac
c$$$     &           + fermfac * gggg(S, T, eps=1) * prefac
c$$$         endif
         if (eps.NE.0) then
            stop 'stdborn: eps for pure QCD Born not implemented'
c     See how it is done in comment above.
         endif

      elseif (flav(1).EQ.flav(2)) then ! -> qb qb q q
         do i = 0, 3            ! reorder to q q -> q q
            pborn(i, 1) = -p(i, 1)
            pborn(i, 2) = -p(i, 2)
            pborn(i, 3) = p(i, 3)
            pborn(i, 4) = p(i, 4)
         enddo
         call mandelstam22(pborn, S, T)
         if ((index(1)*index(2).EQ.2).OR.(index(3)*index(4).EQ.2)) then
            prefac = qqFac(NC) * 0.5d0 ! final state symmetry
         else
            prefac = qqFac(NC)
         endif
         born = fermfac * qqqq(S, T) * prefac
         if (eps.NE.0) then
            stop 'stdborn: eps for pure QCD Born not implemented'
         endif

      elseif (flav(1).EQ.-flav(2)) then ! -> qb q g g
         do i = 0, 3            ! reorder to q qb -> g g
            pborn(i, 1) = -p(i, 1)
            pborn(i, 2) = -p(i, 2)
            pborn(i, 3) = p(i, 3)
            pborn(i, 4) = p(i, 4)
         enddo
         call mandelstam22(pborn, S, T)
         if (index(1)*index(2).EQ.2) then
            prefac = qqFac(NC) * 0.5d0 ! final state symmetry
         elseif (index(3)*index(4).EQ.2) then
            prefac = ggFac(NC)
         else
            prefac = qgFac(NC)
         endif
         born = fermfac * qqbgg(S, T) * prefac
         if (eps.NE.0) then
            stop 'stdborn: eps for pure QCD Born not implemented'
         endif

      elseif ((flav(1).LT.0).AND.(flav(2).LT.0)) then ! -> qb qb' q' q
         do i = 0, 3            ! reorder to q q' -> q q'
            pborn(i, 1) = -p(i, 1)
            pborn(i, 2) = -p(i, 2)
            pborn(i, 3) = p(i, 4)
            pborn(i, 4) = p(i, 3)
         enddo
         call mandelstam22(pborn, S, T)
         prefac = qqFac(NC)
         born = fermfac * qqpqqp(S, T) * prefac
         if (eps.NE.0) then
            stop 'stdborn: eps for pure QCD Born not implemented'
         endif

      else
         stop 'stdborn: nonexistent flavour structure'
      endif

      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccccc
ccc Determine charge on every leg.             ccc
cccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine setcharge(flav, Q)

      implicit none
      include 'nlegborn.h'
      integer flav
      double precision Q

      if (flav.EQ.0) then
         Q = 0.d0
      elseif (MOD(flav,2).EQ.0) then
         Q = 2.d0/3.d0
      else
         Q = -1.d0/3.d0
      endif

      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccc
ccc Compute Mandelstam variables from momenta
ccc of 2->2 process with massless particles
ccc U = -S - T

      subroutine mandelstam22(p, S, T)

      implicit none
      double precision p(0:3,4), S, T
      real * 8 dotp

      S = 2.d0 * dotp(p(0,1),p(0,2))
      T = -2.d0 * dotp(p(0,1),p(0,3))

      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccc
ccc Averaging prefactor for
ccc (anti)quark/(anti)quark initial states

      double precision function qqFac(NC)
      implicit none
      integer NC

      qqFac = 1.d0/(4*NC**2)

      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccc
ccc Averaging prefactor for (anti)quark/gluon
ccc initial states

      double precision function qgFac(NC)
      implicit none
      integer NC

      qgFac = 1.d0/(4*NC*(NC**2-1))

      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccc
ccc Averaging prefactor for gluon/gluon initial
ccc states

      double precision function ggFac(NC)

      implicit none
      integer NC

      ggFac = 1.d0/(4.d0*(NC**2-1)**2)

      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccc
ccc q qbar -> gamma g

      double precision function qqbAg(Q, S, T, eps)
      implicit none
      include 'pwhg_math.h'
      include 'pwhg_st.h'
      include 'pwhg_em.h'
      double precision Q, S, T
      integer eps

      if (eps.EQ.0) then
         qqbAg = (128*em_alpha*st_alpha*CF*NC*Pi**2*Q**2*(T**2 +
     &        (-S-T)**2))/(T*(-S-T))
      elseif (eps.EQ.1) then
         qqbAg =(-128*em_alpha*st_alpha*CF*NC*Pi**2*Q**2*(-2*S**2 - 2*S
     &        *T - 2*T**2))/(T*(S + T))
      elseif (eps.EQ.2) then
         qqbAg =(-128*em_alpha*st_alpha*CF*NC*Pi**2*Q**2*S**2)/(T*(S +
     &        T))
      endif

      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccc
ccc q q -> q q

      double precision function qqqq(S, T)

      implicit none
      include 'pwhg_math.h'
      include 'pwhg_st.h'
      double precision S, T

      qqqq = (64*st_alpha**2*(NC**2 - 1)*Pi**2*
     &     (S**2*T*(S + T) +
     &     NC*(S**4 + 3*S**3*T + 4*S**2*T**2 + 2*S*T**3 +
     &     T**4)))/(NC*T**2*(S + T)**2)

      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccc
ccc q q' -> q q'

      double precision function qqpqqp(S, T)

      implicit none
      include 'pwhg_math.h'
      include 'pwhg_st.h'
      double precision S, T

      qqpqqp = (32*st_alpha**2*(NC**2 - 1)*Pi**2*
     &     (2*S**2 + 2*S*T + T**2))/(T**2)

      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccc
ccc q qb -> g g

      double precision function qqbgg(S, T)

      implicit none
      include 'pwhg_math.h'
      include 'pwhg_st.h'
      double precision S, T

      qqbgg = (-32*st_alpha**2*(-1 + NC**2)*Pi**2*(S**2 + 2*S*T + 2*T
     &     **2)*((-1 + NC**2)*S**2 + 2*NC**2*S*T + 2*NC**2*T**2))/(NC*S
     &     **2*T*(S + T))

      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccc
ccc g g -> g g

      double precision function gggg(S, T)

      implicit none
      include 'pwhg_math.h'
      include 'pwhg_st.h'
      double precision S, T

      gggg = (256*st_alpha**2*NC**2*Pi**2*(S**2 + S*T + T**2)**3
     &     *(-1 + NC**2))/
     &     (S**2*T**2*(S + T)**2)

      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccccc
ccc A function to calculate the colour         ccc
ccc correlators for parton a on leg j with     ccc
ccc parton b on leg k.                         ccc
ccc The return is the colour correlator        ccc
ccc Ta.Tb under the assumption of only         ccc
ccc three coloured particles q,q(b),g in a     ccc
ccc standard model ME.                         ccc
cccccccccccccccccccccccccccccccccccccccccccccccccc

      double precision function colcor(j,k,a,b)

      implicit none
      include 'pwhg_math.h'
      integer j, k, a, b, adum, bdum

c     colour does not care for flavours or charges
      adum = a
      bdum = b
      call deflavour(adum)
      call deflavour(bdum)
      if (j.EQ.k) then
         if ((adum.EQ.1).AND.(bdum.EQ.1)) then
            colcor = CF         ! Tq^2
         else if ((adum.EQ.0).AND.(bdum.EQ.0)) then
            colcor = CA         ! Tg^2
         else
            stop 'Colcor: Different flavours on same leg.'
         endif
      else
         if ((adum.EQ.1).AND.(bdum.EQ.1)) then
            colcor = (CA - 2.d0*CF)/2.d0 ! 1/2(Tg^2-2*Tq^2)
         else if ((adum.EQ.1).AND.(bdum.EQ.0)) then
            colcor = -1.d0/2.d0*CA ! 1/2(Tq^2-Tq^2-Tg^2)
         else if ((adum.EQ.0).AND.(bdum.EQ.1)) then
            colcor = -1.d0/2.d0*CA ! 1/2(Tq^2-Tg^2-Tq^2)
         else
            write(*,*) '[DEBUG] Colcor: Flavour combination',a,b
            stop 'Colcor: Combination does not exist in a b -> gamma c'
         endif
      endif

      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccccc
ccc Since the several functions are            ccc
ccc insensitive to the flavour and charge this ccc
ccc reduces the indices to 0 or 1 for gluon or ccc
ccc quark respectively                         ccc
cccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine deflavour(a)
      implicit none
      integer a

      a=ABS(a)
      if (a.GT.6) then
         stop 'LHA ID is not corresponding to a quark or gluon'
      else if (a.GT.0) then
         a=1
      endif

      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccccc
ccc   Compute eta-vector used in polarization  ccc
ccc   sum for particle with momentum p.        ccc
ccc   eta has to be the fourth vector besides  ccc
ccc   p and the two polarization vectors to    ccc
ccc   span Minkowski space.                    ccc
cccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine etaCnst(p, eta)

      implicit none
      double precision p(0:3), eta(0:3)

      eta(0) = p(0)
      eta(1) = -p(1)
      eta(2) = -p(2)
      eta(3) = -p(3)

      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccccc
ccc   Björken-Drell metric.                    ccc
cccccccccccccccccccccccccccccccccccccccccccccccccc

      double precision function gmunu(mu,nu)

      implicit none
      integer mu, nu

      if (mu.NE.nu) then
         gmunu = 0.d0
      elseif (mu.EQ.0) then
         gmunu = 1.d0
      elseif (mu.LE.3) then
         gmunu = -1.d0
      else
         stop 'gmunu: Dimension > 4'
      endif

      return
      end




      subroutine borncolour_lh
c Sets up the colour for the given flavour configuration
c already filled in the Les Houches interface.
c In case there are several colour structure, one
c should pick one with a probability proportional to
c the value of the corresponding cross section, for the
c kinematics defined in the Les Houches interface
      implicit none
      include 'nlegborn.h'
      include 'LesHouches.h'
      include 'pwhg_kn.h'
      integer cline
      integer iq,ia,iq1,iq2,ia1,ia2,ig1,ig2,j,itmp
      real * 8 s,t,u
      real * 8 dotp

      cline = 500

      if (idup(3).EQ.22) then
         icolup(1,3) = 0
         icolup(2,3) = 0
         if (idup(1).EQ.21) then
            if (idup(2).GT.0) then ! g q -> gamma q
               call qqbgcol(cline, icolup(1,4), icolup(1,2), icolup(1
     &              ,1))
            else                ! g qb -> gamma qb
               call qqbgcol(cline, icolup(1,2), icolup(1,4), icolup(1
     &              ,1))
            endif
         elseif (idup(2).EQ.21) then
            if (idup(1).GT.0) then ! q g -> gamma q
               call qqbgcol(cline, icolup(1,4), icolup(1,1), icolup(1
     &              ,2))
            else                ! qb g -> gamma qb
               call qqbgcol(cline, icolup(1,1), icolup(1,4), icolup(1
     &              ,2))
            endif
         else
            if (idup(1).GT.0) then ! q qb -> gamma g
               call qqbgcol(cline, icolup(1,1), icolup(1,2), icolup(1
     &              ,4))
            else                ! qb q -> gamma g
               call qqbgcol(cline, icolup(1,2), icolup(1,1), icolup(1
     &              ,4))
            endif
         endif

      else
c     Pure QCD colour assignment, taken from dijet (arXiv:1012.3380)
c     g g g g
         if(idup(1).eq.21.and.idup(2).eq.21
     1        .and.idup(3).eq.21.and.idup(4).eq.21) then
            s=2*dotp(kn_cmpborn(0,1),kn_cmpborn(0,2))
            t=-2*dotp(kn_cmpborn(0,2),kn_cmpborn(0,3))
            u=-s-t
            call borncolour4g(icolup(1,1),icolup(1,2),
     1           icolup(1,3),icolup(1,4),s,t,u)
c     q qb g g or permutations-crossing
         elseif(idup(1).eq.21.or.idup(2).eq.21
     1           .or.idup(3).eq.21.or.idup(4).eq.21) then
c     find the quarks and gluons
            ig1=-1
            do j=1,4
               if(idup(j).eq.21) then
                  if(ig1.lt.0) then
                     ig1=j
                  else
                     ig2=j
                  endif
               elseif(idup(j)*istup(j).gt.0) then
                  iq=j
               elseif(idup(j)*istup(j).lt.0) then
                  ia=j
               else
                  write(*,*) 'borncolour_lh: should not be here!'
                  call exit(1)
               endif
            enddo
c     using istup we reverse the sign of incoming particle momenta,
c     so that the choice of colour can be made independently of which
c     particle is incoming.
            s=istup(iq)*istup(ia)*2*dotp(kn_cmpborn(0,iq),kn_cmpborn(0
     &           ,ia))
            t=istup(ia)*istup(ig1)*2*
     1           dotp(kn_cmpborn(0,ia),kn_cmpborn(0,ig1))
            u=-s-t
            call borncolour2g(icolup(1,iq),icolup(1,ia),
     1           icolup(1,ig1),icolup(1,ig2),s,t,u)
c     q q qb qb, or q Q qb Qb, plus permutations-crossing
         else
            iq1=-1
            iq2=-1
            ia1=-1
            ia2=-1
            do j=1,4
               if(idup(j)*istup(j).gt.0) then
                  if(iq1.lt.0) then
                     iq1=j
                  else
                     iq2=j
                  endif
               else
                  if(ia1.lt.0) then
                     ia1=j
                  else
                     ia2=j
                  endif
               endif
            enddo
            if(idup(iq1)*istup(iq1).eq.idup(iq2)*istup(iq2)) then
c     q q qb qb
               s=istup(iq1)*istup(iq2)*2*
     1              dotp(kn_cmpborn(0,iq1),kn_cmpborn(0,iq2))
               t=istup(iq2)*istup(ia1)*2*
     1              dotp(kn_cmpborn(0,iq2),kn_cmpborn(0,ia1))
               u=-s-t
               call borncolour4q(icolup(1,iq1),icolup(1,iq2),
     1              icolup(1,ia1),icolup(1,ia2),s,t,u)
            else
c     q Q qb Qb
               if(istup(iq1)*idup(iq1).ne.-istup(ia1)*idup(ia1)) then
                  itmp=ia1
                  ia1=ia2
                  ia2=itmp
               endif
               call colourjoin4q(icolup(1,iq1),icolup(1,iq2),
     1              icolup(1,ia1),icolup(1,ia2))
            endif
         endif

      endif

c     turn outgoing colours into incoming
      call conjcolor(icolup(1,1))
      call conjcolor(icolup(1,2))

c$$$      write(*,*) '[DEBUG] ****'
c$$$      write(*,*) '[DEBUG] idup:', idup(1), idup(2), '->', idup(3),
c$$$     &     idup(4)
c$$$      write(*,*) '[DEBUG] icolup:', icolup(1,1), icolup(2,1)
c$$$      write(*,*) '[DEBUG] icolup:', icolup(1,2), icolup(2,2)
c$$$      write(*,*) '[DEBUG] icolup:', icolup(1,3), icolup(2,3)
c$$$      write(*,*) '[DEBUG] icolup:', icolup(1,4), icolup(2,4)

      end


cccccccccccccccccccccccccccccccccccccccccccccccccc
ccc Taken from hvq.                            ccc
ccc Turn outgoin into incoming.                ccc
cccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine conjcolor(cl)
      integer cl(2),i
      i=cl(1)
      cl(1)=cl(2)
      cl(2)=i
      end


cccccccccccccccccccccccccccccccccccccccccccccccccc
ccc Routine to assign colourlines to outgoing  ccc
ccc quark-antiquark-gluon.                     ccc
cccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine qqbgcol(cline, qcol, qbcol, gcol)

      implicit none
      integer cline, qcol(2), qbcol(2), gcol(2)

      qcol(1) = cline + 1
      qcol(2) = 0
      qbcol(1) = 0
      qbcol(2) = cline + 2
      gcol(1) = cline + 2
      gcol(2) = cline + 1

      cline = cline + 3

      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccccc
ccc Routines to generate colour for pure QCD   ccc
ccc Born. Taken from dijet (arXiv:1012.3380)!  ccc
cccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine borncolour4g(icol1,icol2,icol3,icol4,s,t,u)
c     g g g g
      implicit none
      integer icol1(2),icol2(2),icol3(2),icol4(2)
      real * 8 s,t,u
      real * 8 rst,rtu,rsu,r
      real * 8 random
c     planar results for st channel
      rst=(t/s+s/t+1)**2
c     su channel
      rsu=(u/s+s/u+1)**2
c     tu channel
      rtu=(t/u+u/t+1)**2
c     Obtained by maxima; check:
c     rst+rsu+rtu=2(3-us/t^2-ut/s^2-st/u^2)
      r=random()*(rst+rsu+rtu)
      if(r.lt.rst) then
         call colourjoin4g(icol1,icol2,icol3,icol4)
      elseif(r.lt.rst+rsu) then
         call colourjoin4g(icol1,icol2,icol4,icol3)
      else
         call colourjoin4g(icol1,icol3,icol2,icol4)
      endif
      end

      subroutine colourjoin4g(icol1,icol2,icol3,icol4)
c     perform a planar colour connection on the planar sequence
c     of gluons
      implicit none
      integer icol1(2),icol2(2),icol3(2),icol4(2)
      integer newcolor
      call getnewcolor(newcolor)
      icol1(2)=newcolor
      icol2(1)=newcolor
      call getnewcolor(newcolor)
      icol2(2)=newcolor
      icol3(1)=newcolor
      call getnewcolor(newcolor)
      icol3(2)=newcolor
      icol4(1)=newcolor
      call getnewcolor(newcolor)
      icol4(2)=newcolor
      icol1(1)=newcolor
      end

      subroutine borncolour2g(icol1,icol2,icol3,icol4,s,t,u)
c     q     qbar  g     g
c     q qb g g
      implicit none
      integer icol1(2),icol2(2),icol3(2),icol4(2)
      real * 8 s,t,u
      real * 8 rt,ru,r
      real * 8 random
c     rt=u/t*(u**2+t**2)/s**2
c     ru=t/u*(u**2+t**2)/s**2
c     obtained by maxima; check: rt+ru=(1/(tu)-2/s^2)*(u^2+t^2)
c     watch out! crossin a fermion line to get q g->q g needs an
c     extra - sign;
      rt=abs(u/t)
      ru=abs(t/u)
      r=random()*(rt+ru)
      if(r.lt.rt) then
         call colourjoin2g(icol1,icol2,icol3,icol4)
      else
         call colourjoin2g(icol1,icol2,icol4,icol3)
      endif
      end

      subroutine colourjoin2g(icol1,icol2,icol3,icol4)
c     q     qbar  g     g
c     perform a planar colour connection on the planar sequence
c     q qbar g g
      implicit none
      integer icol1(2),icol2(2),icol3(2),icol4(2)
      integer newcolor
      icol1(2)=0
      icol2(1)=0
      call getnewcolor(newcolor)
      icol1(1)=newcolor
      icol4(2)=newcolor
      call getnewcolor(newcolor)
      icol4(1)=newcolor
      icol3(2)=newcolor
      call getnewcolor(newcolor)
      icol3(1)=newcolor
      icol2(2)=newcolor
      end


      subroutine borncolour4q(icol1,icol2,icol3,icol4,s,t,u)
c     q     q     qbar  qbar
      implicit none
      integer icol1(2),icol2(2),icol3(2),icol4(2)
      real * 8 s,t,u
      real * 8 rt,ru,r
      real * 8 random
c     the following is from q Q -> Q q
      rt=(s**2+u**2)/t**2
c     q Q->q Q
      ru=(s**2+t**2)/u**2
      r=random()*(rt+ru)
      if(r.lt.rt) then
c     t channel gluon (23 channel; thus colour is exchanged
c     2->4 and 1->3
         call colourjoin4q(icol1,icol2,icol4,icol3)
      else
         call colourjoin4q(icol1,icol2,icol3,icol4)
      endif
      end

      subroutine colourjoin4q(icol1,icol2,icol3,icol4)
c     q     q     qbar  qbar
c     perform a planar colour connection on the planar sequence
c     q qbar g g
      implicit none
      integer icol1(2),icol2(2),icol3(2),icol4(2)
      integer newcolor
      icol1(2)=0
      icol2(2)=0
      icol3(1)=0
      icol4(1)=0
      call getnewcolor(newcolor)
      icol1(1)=newcolor
      icol4(2)=newcolor
      call getnewcolor(newcolor)
      icol2(1)=newcolor
      icol3(2)=newcolor
      end


      subroutine finalize_lh
c     Use this routine that gets called before an event is written
c     to set up the AQEDUP value in the LHEF standard, since it is
c     set to -1 in born_lh
      implicit none
      include 'LesHouches.h'
      include 'pwhg_em.h'

      aqedup = em_alpha

      end
