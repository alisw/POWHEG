c     returns 2 Re(M_B * M_V)/(as/(2pi)), 
c     where M_B is the Born amplitude and 
c     M_V is the finite part of the virtual amplitude
c     The as/(2pi) factor is attached at a later point
      subroutine setvirtual(p,vflav,virtual)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_math.h'
      include 'pwhg_st.h'
c      include 'constants.h'
      real * 8 p(0:3,nlegborn)
      integer vflav(nlegborn)
      real * 8 virtual
      double precision muf, mur
      double precision pstd(0:3,nlegborn), fermfac
      integer stdflav(nlegborn), index(nlegborn)
      integer i, j
      double precision c(-6:6), gamma(-6:6), a, cij, kij, dotp
      double precision born2, bornjk1(nlegborn, nlegborn), 
     &     born, bornjk(nlegborn, nlegborn), bmunu(0:3,0:3,nlegborn)

c     initialize LoopTools at first call
      logical ini
      data ini/.TRUE./
      save ini
      integer cc
      data cc/0/
      save cc

      if (ini) then
         call ltini
c     Delta has to be set to the following if CTs depend on it (see EOF)
c     then one has to uncomment the include constants.h above, too
c$$$  call setdelta(-DLOG(4 * pi) + EulerGamma)
         do j=-6,6
            if(j.eq.0) then
               c(j)=ca
               gamma(j)=(11*ca-4*tf*st_nlight)/6
            else
               c(j)=cf
               gamma(j)=3d0/2*cf
            endif
         enddo
         ini = .FALSE.
      endif

c     set the renormalization scale in LoopTools
      call set_fac_ren_scales(muf,mur)
      call setmudim(mur**2)

c     compute standard ordering of flavours and momenta
      call stdorderedflav(p, vflav, nlegborn, pstd, stdflav, index,
     &     fermfac)
c     get virtuals via crossing
      call stdvirt(pstd, stdflav, index, fermfac, virtual)

c     supply finite terms proportional to a and c_ij in eq. (2.11)
c     (arXiv:1002.2581) to get 
c     V_fin = myvirt - a B^(2) - sum_ij c_ij B_ij^(1)
c     where the superscripts denote the coefficients in an expansion in
c     epsilon
      if (vflav(3).eq.22) then
         call setborneps(p, vflav, 2, born2, bornjk, bmunu)
         call setborneps(p, vflav, 1, born, bornjk1, bmunu)
         a = 0.d0
         do i = 1, nlegborn
            a = a - c(vflav(i))
            do j = 1, nlegborn
               if (i.EQ.j) cycle
               kij = dotp(p(0,i),p(0,j))
               cij = -gamma(vflav(i))/c(vflav(i)) + log(2*kij/mur**2)
               virtual = virtual - cij * bornjk1(i,j)
            enddo
         enddo
         virtual = virtual - a * born2
      endif

c     prevent LoopTools from cluttering up the memory
      cc = cc + 1
      if (cc.GE.10000) then    ! works on my machine with 8GB RAM
         call clearcache
         cc = 0
      endif

      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccc
ccc For standard ordered momentum (p) and    ccc
ccc flavour (flav) arrays return the         ccc
ccc corresponding virtual*born amplitude.    ccc
cccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine stdvirt(p, flav, index, fermfac, virt)

      implicit none
      include 'nlegborn.h'
      include 'pwhg_math.h'
      double precision p(0:3, nlegborn), fermfac, virt
      integer flav(nlegborn), index(nlegborn)
      double precision pvirt(0:3, nlegborn)
      double precision Q, S, T, prefac
      integer i
      double precision qqFac, qgFac, qqbAgVIRT

c     initialize to zero
      virt = 0.d0

      if (flav(nlegborn).EQ.22) then ! -> qb q g A
c     reorder to q qb -> A g (no fermfac needed for two fermion flips)
         do i = 0, 3
            pvirt(i, 1) = -p(i, 1)
            pvirt(i, 2) = -p(i, 2)
            pvirt(i, 3) = p(i, 4)
            pvirt(i, 4) = p(i, 3)
         enddo
         call mandelstam22(pvirt, S, T)
         call setcharge(flav(2), Q)
c     find out which particles are initial, i.e. which have index = 1, 2
         if (index(1)*index(2).EQ.2) then
            prefac = qqFac(NC)
         elseif ((index(1)*index(3).EQ.2).OR.(index(2)*index(3).EQ.2))
     &           then
            prefac = qgFac(NC)
         else
            stop 'stdvirt: illegal initial particles'
         endif

         virt = fermfac * qqbAgVIRT(Q, S, T) * prefac
      else
c$$$         stop 'stdvirt: pure QCD with QED corrections not implemented'
         virt = 0.d0
      endif

      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccc
ccc q qbar -> gamma g
cccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccc


cccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccc
ccc This is the complete virtual correction
ccc to the process q qbar -> gamma g
ccc including UV counterterms (but IR divergent)
ccc All other 2->2 processes can be obtained
ccc from this via crossing

      double precision function qqbAgVIRT(Q, S, T)
      implicit none
      double precision Q, S, T, U
      double precision qqbAgExt, qqbAgInt, qqbAgVert, qqbAgBox,
     &     qqbAgCT

      U = -S -T
      qqbAgVIRT = qqbAgExt(Q, S, T, U) + qqbAgInt(Q, S, T, U)
     &     + qqbAgVert(Q, S, T, U) + qqbAgBox(Q, S, T, U)
     &     + qqbAgCT(Q, S, T, U)

      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccc
ccc External corrections (via LSZ)
ccc minus contributions proportional to
ccc 1/eps * Born (in D dimensions)

      double precision function qqbAgExt(Q, S, T, U)
      implicit none
      include 'pwhg_math.h'
      include 'pwhg_st.h'
      include 'pwhg_em.h'
      double precision Q, S, T, U
      double precision qqbAg

      qqbAgExt = (-64*em_alpha*st_alpha*CF*NC*(6*CF - 5*NC
     &     + 2*st_nlight)*Pi**2*Q**2*(T**2 + T*U + U**2))/(3.*T*U)

      qqbAgExt = 2 * qqbAgExt

      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccc
ccc Internal corrections
ccc Result corresponds to the sum of all
ccc interference diagrams between treelevel
ccc and one loop on an internal propagator

      double precision function qqbAgInt(Q, S, T, U)
      implicit none
      include 'pwhg_math.h'
      include 'pwhg_st.h'
      include 'pwhg_em.h'
      double precision Q, S, T, U
      double complex con0, con1, con2 ! for contributions to O(1/eps,1)
      double complex qqbAgIntC ! the complex result to be taken 2*Re(..)
      double complex b0xT00, b1xT00, b0xU00, b1xU00

      call setBeps(-2, S, T, U, b0xT00, b1xT00, b0xU00, b1xU00)
      con2 = (-128*em_alpha*st_alpha*CF**2*NC*Pi**2*Q**2*
     &    (U*(2*T + 3*U)*b0xT00 +
     &      T*(3*T + 2*U)*b0xU00 +
     &      U*(2*T + 3*U)*b1xT00 +
     &      T*(3*T + 2*U)*b1xU00))/(T*U)

      call setBeps(-1, S, T, U, b0xT00, b1xT00, b0xU00, b1xU00)
      con1 = (128*em_alpha*st_alpha*CF**2*NC*Pi**2*Q**2*
     &    (U*(T + 3*U)*b0xT00 +
     &      T*(3*T + U)*b0xU00 +
     &      U*(T + 3*U)*b1xT00 +
     &      T*(3*T + U)*b1xU00))/(T*U)

      call setBeps(0, S, T, U, b0xT00, b1xT00, b0xU00, b1xU00)
      con0 = (-128*em_alpha*st_alpha*CF**2*NC*Pi**2*Q**2*
     &    (U**2*(b0xT00 + b1xT00) +
     &      T**2*(b0xU00 + b1xU00)))/(T*U)

      qqbAgIntC = con0 + con1 + con2

      qqbAgInt = 2 * DBLE(qqbAgIntC)

      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccc
ccc subroutine to set the values of the PV functions
ccc no idea how to solve this nicer w.o. rewriting
ccc my MEs...

      subroutine setBeps(eps, S, T, U, b0xT00, b1xT00, b0xU00, b1xU00)
      implicit none
#include "looptools.h"
      double complex b0xT00, b1xT00, b0xU00, b1xU00
      double precision S, T, U, deps
      integer eps

      deps = DBLE(eps)
      call setlambda(deps)

      b0xT00 = B0i(bb0,T,0.d0,0.d0)
      b1xT00 = B0i(bb1,T,0.d0,0.d0)
      b0xU00 = B0i(bb0,U,0.d0,0.d0)
      b1xU00 = B0i(bb1,U,0.d0,0.d0)

      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccc
ccc Vertex corrections

      double precision function qqbAgVert(Q, S, T, U)
      implicit none
      include 'pwhg_math.h'
      include 'pwhg_st.h'
      include 'pwhg_em.h'
      double precision Q, S, T, U
      double complex con0, con1, con2 ! for contributions to O(1/eps,1)
      double complex qqbAgVertC ! the complex result to be taken 2*Re(..)
      double complex c0x00T000, c0x00U000, c0xT00000, c0xU00000,
     &     c00x00T000, c00x00U000, c00xT00000, c00xU00000, c1x00T000,
     &     c1x00U000, c1xT00000, c1xU00000, c11xT00000, c11xU00000,
     &     c12x00T000, c12x00U000, c12xT00000, c12xU00000, c2x00T000,
     &     c2x00U000, c2xT00000, c2xU00000, c22x00T000, c22x00U000

      call setCeps(-2, S, T, U,
     &     c0x00T000, c0x00U000, c0xT00000, c0xU00000,
     &     c00x00T000, c00x00U000, c00xT00000, c00xU00000, c1x00T000,
     &     c1x00U000, c1xT00000, c1xU00000, c11xT00000, c11xU00000,
     &     c12x00T000, c12x00U000, c12xT00000, c12xU00000, c2x00T000,
     &     c2x00U000, c2xT00000, c2xU00000, c22x00T000, c22x00U000)
      con2 = (32*em_alpha*st_alpha*CF*NC*Pi**2*Q**2*
     &    (4*CF*T*U*(T + U)*c0x00T000 +
     &      2*T*U*(-(NC*U) + 2*CF*(T + U))*
     &       c0x00U000 +
     &      4*CF*T**2*U*c0xT00000 -
     &      2*NC*T**2*U*c0xT00000 +
     &      4*CF*T*U**2*c0xT00000 -
     &      2*NC*T*U**2*c0xT00000 +
     &      4*CF*T**2*U*c0xU00000 -
     &      2*NC*T**2*U*c0xU00000 +
     &      4*CF*T*U**2*c0xU00000 +
     &      24*CF*T*U*c00x00T000 +
     &      48*CF*U**2*c00x00T000 +
     &      48*CF*T**2*c00x00U000 +
     &      24*CF*T*U*c00x00U000 +
     &      8*NC*T*U*c00x00U000 +
     &      24*CF*T*U*c00xT00000 +
     &      8*NC*T*U*c00xT00000 +
     &      48*CF*U**2*c00xT00000 +
     &      4*NC*U**2*c00xT00000 +
     &      48*CF*T**2*c00xU00000 +
     &      4*NC*T**2*c00xU00000 +
     &      24*CF*T*U*c00xU00000 +
     &      4*CF*T**2*U*c1x00T000 +
     &      4*CF*T*U**2*c1x00T000 +
     &      4*CF*T**2*U*c1x00U000 +
     &      4*CF*T*U**2*c1x00U000 +
     &      8*CF*T**2*U*c1xT00000 +
     &      NC*T**2*U*c1xT00000 +
     &      12*CF*T*U**2*c1xT00000 -
     &      NC*T*U**2*c1xT00000 +
     &      12*CF*T**2*U*c1xU00000 -
     &      2*NC*T**2*U*c1xU00000 +
     &      8*CF*T*U**2*c1xU00000 +
     &      3*NC*T*U**2*c1xU00000 +
     &      4*CF*T**2*U*c11xT00000 +
     &      2*NC*T**2*U*c11xT00000 +
     &      8*CF*T*U**2*c11xT00000 +
     &      8*CF*T**2*U*c11xU00000 +
     &      4*CF*T*U**2*c11xU00000 +
     &      4*NC*T*U**2*c11xU00000 +
     &      4*CF*T**2*U*c12x00T000 +
     &      8*CF*T*U**2*c12x00T000 +
     &      8*CF*T**2*U*c12x00U000 +
     &      4*NC*T**2*U*c12x00U000 +
     &      4*CF*T*U**2*c12x00U000 +
     &      2*NC*T*U**2*c12x00U000 +
     &      4*CF*T**2*U*c12xT00000 +
     &      2*NC*T**2*U*c12xT00000 +
     &      8*CF*T*U**2*c12xT00000 +
     &      8*CF*T**2*U*c12xU00000 -
     &      4*NC*T**2*U*c12xU00000 +
     &      4*CF*T*U**2*c12xU00000 +
     &      8*CF*T**2*U*c2x00T000 +
     &      12*CF*T*U**2*c2x00T000 +
     &      12*CF*T**2*U*c2x00U000 +
     &      NC*T**2*U*c2x00U000 +
     &      8*CF*T*U**2*c2x00U000 -
     &      2*NC*T*U**2*c2x00U000 +
     &      4*CF*T**2*U*c2xT00000 +
     &      4*CF*T*U**2*c2xT00000 +
     &      4*CF*T**2*U*c2xU00000 +
     &      4*CF*T*U**2*c2xU00000 +
     &      4*CF*T**2*U*c22x00T000 +
     &      8*CF*T*U**2*c22x00T000 +
     &      2*T*U*(-(NC*U) + 2*CF*(2*T + U))*
     &       c22x00U000))/(T*U)

      call setCeps(-1, S, T, U,
     &     c0x00T000, c0x00U000, c0xT00000, c0xU00000,
     &     c00x00T000, c00x00U000, c00xT00000, c00xU00000, c1x00T000,
     &     c1x00U000, c1xT00000, c1xU00000, c11xT00000, c11xU00000,
     &     c12x00T000, c12x00U000, c12xT00000, c12xU00000, c2x00T000,
     &     c2x00U000, c2xT00000, c2xU00000, c22x00T000, c22x00U000)
      con1 = (-32*em_alpha*st_alpha*CF*NC*Pi**2*Q**2*
     &    (4*CF*T*U**2*c0x00T000 +
     &      2*(2*CF + NC)*T**2*U*c0x00U000 +
     &      4*CF*T*U**2*c0xT00000 -
     &      2*NC*T*U**2*c0xT00000 +
     &      4*CF*T**2*U*c0xU00000 -
     &      4*NC*T**2*U*c0xU00000 +
     &      8*CF*T*U*c00x00T000 +
     &      32*CF*U**2*c00x00T000 +
     &      32*CF*T**2*c00x00U000 +
     &      8*CF*T*U*c00x00U000 +
     &      8*NC*T*U*c00x00U000 +
     &      8*CF*T*U*c00xT00000 +
     &      8*NC*T*U*c00xT00000 +
     &      32*CF*U**2*c00xT00000 +
     &      16*NC*U**2*c00xT00000 +
     &      32*CF*T**2*c00xU00000 +
     &      16*NC*T**2*c00xU00000 +
     &      8*CF*T*U*c00xU00000 +
     &      4*CF*T*U**2*c1x00T000 +
     &      4*CF*T**2*U*c1x00U000 +
     &      NC*T**2*U*c1x00U000 +
     &      NC*T*U**2*c1x00U000 -
     &      4*CF*T**2*U*c1xT00000 +
     &      7*NC*T**2*U*c1xT00000 +
     &      8*CF*T*U**2*c1xT00000 +
     &      4*NC*T*U**2*c1xT00000 +
     &      8*CF*T**2*U*c1xU00000 +
     &      2*NC*T**2*U*c1xU00000 -
     &      4*CF*T*U**2*c1xU00000 +
     &      6*NC*T*U**2*c1xU00000 -
     &      4*CF*T**2*U*c11xT00000 +
     &      8*NC*T**2*U*c11xT00000 +
     &      4*CF*T*U**2*c11xT00000 +
     &      6*NC*T*U**2*c11xT00000 +
     &      4*CF*T**2*U*c11xU00000 +
     &      6*NC*T**2*U*c11xU00000 -
     &      4*CF*T*U**2*c11xU00000 +
     &      8*NC*T*U**2*c11xU00000 -
     &      4*CF*T**2*U*c12x00T000 +
     &      4*CF*T*U**2*c12x00T000 +
     &      4*CF*T**2*U*c12x00U000 +
     &      8*NC*T**2*U*c12x00U000 -
     &      4*CF*T*U**2*c12x00U000 +
     &      8*NC*T*U**2*c12x00U000 -
     &      4*CF*T**2*U*c12xT00000 +
     &      8*NC*T**2*U*c12xT00000 +
     &      4*CF*T*U**2*c12xT00000 +
     &      6*NC*T*U**2*c12xT00000 +
     &      4*CF*T**2*U*c12xU00000 -
     &      2*NC*T**2*U*c12xU00000 -
     &      4*CF*T*U**2*c12xU00000 -
     &      4*CF*T**2*U*c2x00T000 +
     &      8*CF*T*U**2*c2x00T000 +
     &      8*CF*T**2*U*c2x00U000 +
     &      2*NC*T**2*U*c2x00U000 -
     &      4*CF*T*U**2*c2x00U000 +
     &      NC*T*U**2*c2x00U000 +
     &      NC*T**2*U*c2xT00000 +
     &      4*CF*T*U**2*c2xT00000 +
     &      NC*T*U**2*c2xT00000 +
     &      4*CF*T**2*U*c2xU00000 -
     &      4*CF*T**2*U*c22x00T000 +
     &      4*CF*T*U**2*c22x00T000 +
     &      4*CF*T*(T - U)*U*c22x00U000))/(T*U)

      call setCeps(0, S, T, U,
     &     c0x00T000, c0x00U000, c0xT00000, c0xU00000,
     &     c00x00T000, c00x00U000, c00xT00000, c00xU00000, c1x00T000,
     &     c1x00U000, c1xT00000, c1xU00000, c11xT00000, c11xU00000,
     &     c12x00T000, c12x00U000, c12xT00000, c12xU00000, c2x00T000,
     &     c2x00U000, c2xT00000, c2xU00000, c22x00T000, c22x00U000)
      con0 = (-32*em_alpha*st_alpha*CF*NC*Pi**2*Q**2*
     &    (4*CF*T**2*U*c0x00T000 -
     &      2*T*U*(-2*CF*U + NC*(T + U))*
     &       c0x00U000 +
     &      4*CF*T**2*U*c0xT00000 -
     &      2*NC*T**2*U*c0xT00000 +
     &      2*NC*T**2*U*c0xU00000 +
     &      4*CF*T*U**2*c0xU00000 -
     &      8*CF*U**2*c00x00T000 -
     &      8*CF*T**2*c00x00U000 -
     &      8*CF*U**2*c00xT00000 -
     &      8*NC*U**2*c00xT00000 -
     &      8*CF*T**2*c00xU00000 -
     &      8*NC*T**2*c00xU00000 +
     &      4*CF*T**2*U*c1x00T000 -
     &      NC*T**2*U*c1x00U000 +
     &      4*CF*T*U**2*c1x00U000 -
     &      NC*T*U**2*c1x00U000 +
     &      8*CF*T**2*U*c1xT00000 -
     &      4*NC*T**2*U*c1xT00000 -
     &      3*NC*T*U**2*c1xT00000 -
     &      2*NC*T**2*U*c1xU00000 +
     &      8*CF*T*U**2*c1xU00000 -
     &      3*NC*T*U**2*c1xU00000 +
     &      4*CF*T**2*U*c11xT00000 -
     &      4*NC*T**2*U*c11xT00000 -
     &      4*NC*T*U**2*c11xT00000 -
     &      4*NC*T**2*U*c11xU00000 +
     &      4*CF*T*U**2*c11xU00000 -
     &      4*NC*T*U**2*c11xU00000 +
     &      4*CF*T**2*U*c12x00T000 -
     &      4*NC*T**2*U*c12x00U000 +
     &      4*CF*T*U**2*c12x00U000 -
     &      4*NC*T*U**2*c12x00U000 +
     &      4*CF*T**2*U*c12xT00000 -
     &      4*NC*T**2*U*c12xT00000 -
     &      4*NC*T*U**2*c12xT00000 +
     &      4*CF*T*U**2*c12xU00000 +
     &      8*CF*T**2*U*c2x00T000 -
     &      NC*T**2*U*c2x00U000 +
     &      8*CF*T*U**2*c2x00U000 -
     &      NC*T*U**2*c2x00U000 +
     &      4*CF*T**2*U*c2xT00000 -
     &      NC*T**2*U*c2xT00000 -
     &      NC*T*U**2*c2xT00000 +
     &      4*CF*T*U**2*c2xU00000 +
     &      4*CF*T**2*U*c22x00T000 +
     &      4*CF*T*U**2*c22x00U000))/(T*U)

      qqbAgVertC = con0 + con1 + con2

      qqbAgVert = 2 * DBLE(qqbAgVertC)

      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccc
ccc subroutine to set the values of the PV functions

      subroutine setCeps(eps, S, T, U,
     &     c0x00T000, c0x00U000, c0xT00000, c0xU00000,
     &     c00x00T000, c00x00U000, c00xT00000, c00xU00000, c1x00T000,
     &     c1x00U000, c1xT00000, c1xU00000, c11xT00000, c11xU00000,
     &     c12x00T000, c12x00U000, c12xT00000, c12xU00000, c2x00T000,
     &     c2x00U000, c2xT00000, c2xU00000, c22x00T000, c22x00U000)
      implicit none
#include "looptools.h"
      double complex c0x00T000, c0x00U000, c0xT00000, c0xU00000,
     &     c00x00T000, c00x00U000, c00xT00000, c00xU00000, c1x00T000,
     &     c1x00U000, c1xT00000, c1xU00000, c11xT00000, c11xU00000,
     &     c12x00T000, c12x00U000, c12xT00000, c12xU00000, c2x00T000,
     &     c2x00U000, c2xT00000, c2xU00000, c22x00T000, c22x00U000
      double precision S, T, U, deps
      integer eps

      deps = DBLE(eps)
      call setlambda(deps)

      c0x00T000 = C0i(cc0,0.d0,0.d0,T,0.d0,0.d0,0.d0)
      c0x00U000 = C0i(cc0,0.d0,0.d0,U,0.d0,0.d0,0.d0)
      c0xT00000 = C0i(cc0,T,0.d0,0.d0,0.d0,0.d0,0.d0)
      c0xU00000 = C0i(cc0,U,0.d0,0.d0,0.d0,0.d0,0.d0)
      c00x00T000 = C0i(cc00,0.d0,0.d0,T,0.d0,0.d0,0.d0)
      c00x00U000 = C0i(cc00,0.d0,0.d0,U,0.d0,0.d0,0.d0)
      c00xT00000 = C0i(cc00,T,0.d0,0.d0,0.d0,0.d0,0.d0)
      c00xU00000 = C0i(cc00,U,0.d0,0.d0,0.d0,0.d0,0.d0)
      c1x00T000 = C0i(cc1,0.d0,0.d0,T,0.d0,0.d0,0.d0)
      c1x00U000 = C0i(cc1,0.d0,0.d0,U,0.d0,0.d0,0.d0)
      c1xT00000 = C0i(cc1,T,0.d0,0.d0,0.d0,0.d0,0.d0)
      c1xU00000 = C0i(cc1,U,0.d0,0.d0,0.d0,0.d0,0.d0)
      c11xT00000 = C0i(cc11,T,0.d0,0.d0,0.d0,0.d0,0.d0)
      c11xU00000 = C0i(cc11,U,0.d0,0.d0,0.d0,0.d0,0.d0)
      c12x00T000 = C0i(cc12,0.d0,0.d0,T,0.d0,0.d0,0.d0)
      c12x00U000 = C0i(cc12,0.d0,0.d0,U,0.d0,0.d0,0.d0)
      c12xT00000 = C0i(cc12,T,0.d0,0.d0,0.d0,0.d0,0.d0)
      c12xU00000 = C0i(cc12,U,0.d0,0.d0,0.d0,0.d0,0.d0)
      c2x00T000 = C0i(cc2,0.d0,0.d0,T,0.d0,0.d0,0.d0)
      c2x00U000 = C0i(cc2,0.d0,0.d0,U,0.d0,0.d0,0.d0)
      c2xT00000 = C0i(cc2,T,0.d0,0.d0,0.d0,0.d0,0.d0)
      c2xU00000 = C0i(cc2,U,0.d0,0.d0,0.d0,0.d0,0.d0)
      c22x00T000 = C0i(cc22,0.d0,0.d0,T,0.d0,0.d0,0.d0)
      c22x00U000 = C0i(cc22,0.d0,0.d0,U,0.d0,0.d0,0.d0)

      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccc
ccc Box corrections

      double precision function qqbAgBox(Q, S, T, U)
      implicit none
      include 'pwhg_math.h'
      include 'pwhg_st.h'
      include 'pwhg_em.h'
      double precision Q, S, T, U
      double complex con0, con1, con2 ! for contributions to O(1/eps,1)
      double complex qqbAgBoxC ! the complex result to be taken 2*Re(..)
      double complex d0x0000TS0000, d0x0000TU0000, d0x0000US0000,
     &     d00x0000TS0000, d00x0000TU0000, d00x0000US0000,
     &     d001x0000TS0000, d001x0000TU0000, d001x0000US0000,
     &     d002x0000TS0000, d002x0000TU0000, d002x0000US0000,
     &     d003x0000TS0000, d003x0000TU0000, d003x0000US0000,
     &     d1x0000TS0000, d1x0000TU0000, d1x0000US0000,
     &     d11x0000TS0000, d11x0000TU0000, d11x0000US0000,
     &     d112x0000TS0000, d112x0000TU0000, d112x0000US0000,
     &     d113x0000TS0000, d113x0000TU0000, d113x0000US0000,
     &     d12x0000TS0000, d12x0000TU0000, d12x0000US0000,
     &     d122x0000TS0000, d122x0000TU0000, d122x0000US0000,
     &     d123x0000TS0000, d123x0000TU0000, d123x0000US0000,
     &     d13x0000TS0000, d13x0000TU0000, d13x0000US0000,
     &     d133x0000TS0000, d133x0000US0000, d2x0000TS0000,
     &     d2x0000TU0000, d2x0000US0000, d22x0000TS0000,
     &     d22x0000TU0000, d22x0000US0000, d222x0000TS0000,
     &     d222x0000TU0000, d222x0000US0000, d223x0000TS0000,
     &     d223x0000TU0000, d223x0000US0000, d23x0000TS0000,
     &     d23x0000TU0000, d23x0000US0000, d233x0000TS0000,
     &     d233x0000US0000, d3x0000TS0000, d3x0000TU0000,
     &     d3x0000US0000, d33x0000TS0000, d33x0000TU0000,
     &     d33x0000US0000, d133x0000TU0000, d233x0000TU0000


      call setDeps(-2, S, T, U,
     &     d0x0000TS0000, d0x0000TU0000, d0x0000US0000,
     &     d00x0000TS0000, d00x0000TU0000, d00x0000US0000,
     &     d001x0000TS0000, d001x0000TU0000, d001x0000US0000,
     &     d002x0000TS0000, d002x0000TU0000, d002x0000US0000,
     &     d003x0000TS0000, d003x0000TU0000, d003x0000US0000,
     &     d1x0000TS0000, d1x0000TU0000, d1x0000US0000,
     &     d11x0000TS0000, d11x0000TU0000, d11x0000US0000,
     &     d112x0000TS0000, d112x0000TU0000, d112x0000US0000,
     &     d113x0000TS0000, d113x0000TU0000, d113x0000US0000,
     &     d12x0000TS0000, d12x0000TU0000, d12x0000US0000,
     &     d122x0000TS0000, d122x0000TU0000, d122x0000US0000,
     &     d123x0000TS0000, d123x0000TU0000, d123x0000US0000,
     &     d13x0000TS0000, d13x0000TU0000, d13x0000US0000,
     &     d133x0000TS0000, d133x0000US0000, d2x0000TS0000,
     &     d2x0000TU0000, d2x0000US0000, d22x0000TS0000,
     &     d22x0000TU0000, d22x0000US0000, d222x0000TS0000,
     &     d222x0000TU0000, d222x0000US0000, d223x0000TS0000,
     &     d223x0000TU0000, d223x0000US0000, d23x0000TS0000,
     &     d23x0000TU0000, d23x0000US0000, d233x0000TS0000,
     &     d233x0000US0000, d3x0000TS0000, d3x0000TU0000,
     &     d3x0000US0000, d33x0000TS0000, d33x0000TU0000,
     &     d33x0000US0000, d133x0000TU0000, d233x0000TU0000)
      con2 = 32*em_alpha*st_alpha*CF*NC*Pi**2*Q**2*
     &  (2*(2*CF - NC)*(T + U)**2*
     &     d0x0000TS0000 -
     &    2*NC*U*(T + U)*d0x0000TU0000 +
     &    4*CF*T**2*d0x0000US0000 -
     &    2*NC*T**2*d0x0000US0000 +
     &    8*CF*T*U*d0x0000US0000 -
     &    4*NC*T*U*d0x0000US0000 +
     &    4*CF*U**2*d0x0000US0000 -
     &    2*NC*U**2*d0x0000US0000 +
     &    8*CF*T*d00x0000TS0000 -
     &    4*NC*T*d00x0000TS0000 +
     &    32*CF*U*d00x0000TS0000 -
     &    16*NC*U*d00x0000TS0000 +
     &    18*NC*T*d00x0000TU0000 +
     &    16*NC*U*d00x0000TU0000 +
     &    32*CF*T*d00x0000US0000 -
     &    16*NC*T*d00x0000US0000 +
     &    8*CF*U*d00x0000US0000 -
     &    4*NC*U*d00x0000US0000 +
     &    24*CF*T*d001x0000TS0000 -
     &    12*NC*T*d001x0000TS0000 +
     &    24*CF*U*d001x0000TS0000 -
     &    12*NC*U*d001x0000TS0000 +
     &    2*NC*T*d001x0000TU0000 -
     &    2*NC*U*d001x0000TU0000 +
     &    24*CF*T*d001x0000US0000 -
     &    12*NC*T*d001x0000US0000 +
     &    24*CF*U*d001x0000US0000 -
     &    12*NC*U*d001x0000US0000 +
     &    56*CF*T*d002x0000TS0000 -
     &    28*NC*T*d002x0000TS0000 +
     &    32*CF*U*d002x0000TS0000 -
     &    16*NC*U*d002x0000TS0000 -
     &    6*NC*T*d002x0000TU0000 -
     &    10*NC*U*d002x0000TU0000 +
     &    32*CF*T*d002x0000US0000 -
     &    16*NC*T*d002x0000US0000 +
     &    56*CF*U*d002x0000US0000 -
     &    28*NC*U*d002x0000US0000 +
     &    24*CF*T*d003x0000TS0000 -
     &    12*NC*T*d003x0000TS0000 +
     &    24*CF*U*d003x0000TS0000 -
     &    12*NC*U*d003x0000TS0000 -
     &    8*NC*T*d003x0000TU0000 -
     &    8*NC*U*d003x0000TU0000 +
     &    24*CF*T*d003x0000US0000 -
     &    12*NC*T*d003x0000US0000 +
     &    24*CF*U*d003x0000US0000 -
     &    12*NC*U*d003x0000US0000 +
     &    4*CF*T**2*d1x0000TS0000 -
     &    2*NC*T**2*d1x0000TS0000 +
     &    8*CF*T*U*d1x0000TS0000 -
     &    4*NC*T*U*d1x0000TS0000 +
     &    4*CF*U**2*d1x0000TS0000 -
     &    2*NC*U**2*d1x0000TS0000 +
     &    NC*T**2*d1x0000TU0000 -
     &    2*NC*T*U*d1x0000TU0000 -
     &    3*NC*U**2*d1x0000TU0000 +
     &    4*CF*T**2*d1x0000US0000 -
     &    2*NC*T**2*d1x0000US0000 +
     &    8*CF*T*U*d1x0000US0000 -
     &    4*NC*T*U*d1x0000US0000 +
     &    4*CF*U**2*d1x0000US0000 -
     &    2*NC*U**2*d1x0000US0000 -
     &    NC*T**2*d11x0000TU0000 -
     &    2*NC*T*U*d11x0000TU0000 -
     &    NC*U**2*d11x0000TU0000 +
     &    4*CF*T**2*d112x0000TS0000 -
     &    2*NC*T**2*d112x0000TS0000 +
     &    4*CF*T*U*d112x0000TS0000 -
     &    2*NC*T*U*d112x0000TS0000 +
     &    4*CF*T*U*d112x0000US0000 -
     &    2*NC*T*U*d112x0000US0000 +
     &    4*CF*U**2*d112x0000US0000 -
     &    2*NC*U**2*d112x0000US0000 +
     &    4*CF*T**2*d113x0000TS0000 -
     &    2*NC*T**2*d113x0000TS0000 +
     &    8*CF*T*U*d113x0000TS0000 -
     &    4*NC*T*U*d113x0000TS0000 +
     &    4*CF*U**2*d113x0000TS0000 -
     &    2*NC*U**2*d113x0000TS0000 +
     &    4*CF*T**2*d113x0000US0000 -
     &    2*NC*T**2*d113x0000US0000 +
     &    8*CF*T*U*d113x0000US0000 -
     &    4*NC*T*U*d113x0000US0000 +
     &    4*CF*U**2*d113x0000US0000 -
     &    2*NC*U**2*d113x0000US0000 +
     &    4*CF*T**2*d12x0000TS0000 -
     &    2*NC*T**2*d12x0000TS0000 +
     &    4*CF*T*U*d12x0000TS0000 -
     &    2*NC*T*U*d12x0000TS0000 +
     &    3*NC*T**2*d12x0000TU0000 +
     &    2*NC*T*U*d12x0000TU0000 -
     &    NC*U**2*d12x0000TU0000 +
     &    4*CF*T*U*d12x0000US0000 -
     &    2*NC*T*U*d12x0000US0000 +
     &    4*CF*U**2*d12x0000US0000 -
     &    2*NC*U**2*d12x0000US0000 +
     &    12*CF*T**2*d122x0000TS0000 -
     &    6*NC*T**2*d122x0000TS0000 +
     &    8*CF*T*U*d122x0000TS0000 -
     &    4*NC*T*U*d122x0000TS0000 -
     &    2*NC*T**2*d122x0000TU0000 -
     &    2*NC*T*U*d122x0000TU0000 +
     &    8*CF*T*U*d122x0000US0000 -
     &    4*NC*T*U*d122x0000US0000 +
     &    12*CF*U**2*d122x0000US0000 -
     &    6*NC*U**2*d122x0000US0000 +
     &    16*CF*T**2*d123x0000TS0000 -
     &    8*NC*T**2*d123x0000TS0000 +
     &    20*CF*T*U*d123x0000TS0000 -
     &    10*NC*T*U*d123x0000TS0000 +
     &    4*CF*U**2*d123x0000TS0000 -
     &    2*NC*U**2*d123x0000TS0000 -
     &    2*NC*T**2*d123x0000TU0000 +
     &    2*NC*U**2*d123x0000TU0000 +
     &    4*CF*T**2*d123x0000US0000 -
     &    2*NC*T**2*d123x0000US0000 +
     &    20*CF*T*U*d123x0000US0000 -
     &    10*NC*T*U*d123x0000US0000 +
     &    16*CF*U**2*d123x0000US0000 -
     &    8*NC*U**2*d123x0000US0000 +
     &    4*CF*T*U*d13x0000TS0000 -
     &    2*NC*T*U*d13x0000TS0000 +
     &    4*CF*U**2*d13x0000TS0000 -
     &    2*NC*U**2*d13x0000TS0000 -
     &    NC*T**2*d13x0000TU0000 -
     &    6*NC*T*U*d13x0000TU0000 -
     &    5*NC*U**2*d13x0000TU0000 +
     &    4*CF*T**2*d13x0000US0000 -
     &    2*NC*T**2*d13x0000US0000 +
     &    4*CF*T*U*d13x0000US0000 -
     &    2*NC*T*U*d13x0000US0000 +
     &    4*CF*T**2*d133x0000TS0000 -
     &    2*NC*T**2*d133x0000TS0000 +
     &    8*CF*T*U*d133x0000TS0000 -
     &    4*NC*T*U*d133x0000TS0000 +
     &    4*CF*U**2*d133x0000TS0000 -
     &    2*NC*U**2*d133x0000TS0000 +
     &    2*NC*T*U*d133x0000TU0000 +
     &    2*NC*U**2*d133x0000TU0000 +
     &    4*CF*T**2*d133x0000US0000 -
     &    2*NC*T**2*d133x0000US0000 +
     &    8*CF*T*U*d133x0000US0000 -
     &    4*NC*T*U*d133x0000US0000 +
     &    4*CF*U**2*d133x0000US0000 -
     &    2*NC*U**2*d133x0000US0000 +
     &    4*CF*T**2*d2x0000TS0000 -
     &    2*NC*T**2*d2x0000TS0000 +
     &    8*CF*T*U*d2x0000TS0000 -
     &    4*NC*T*U*d2x0000TS0000 +
     &    4*CF*U**2*d2x0000TS0000 -
     &    2*NC*U**2*d2x0000TS0000 +
     &    4*NC*T**2*d2x0000TU0000 +
     &    2*NC*T*U*d2x0000TU0000 -
     &    2*NC*U**2*d2x0000TU0000 +
     &    4*CF*T**2*d2x0000US0000 -
     &    2*NC*T**2*d2x0000US0000 +
     &    8*CF*T*U*d2x0000US0000 -
     &    4*NC*T*U*d2x0000US0000 +
     &    4*CF*U**2*d2x0000US0000 -
     &    2*NC*U**2*d2x0000US0000 +
     &    8*CF*T**2*d22x0000TS0000 -
     &    4*NC*T**2*d22x0000TS0000 +
     &    4*CF*T*U*d22x0000TS0000 -
     &    2*NC*T*U*d22x0000TS0000 +
     &    2*NC*T**2*d22x0000TU0000 +
     &    2*NC*T*U*d22x0000TU0000 +
     &    4*CF*T*U*d22x0000US0000 -
     &    2*NC*T*U*d22x0000US0000 +
     &    8*CF*U**2*d22x0000US0000 -
     &    4*NC*U**2*d22x0000US0000 +
     &    8*CF*T**2*d222x0000TS0000 -
     &    4*NC*T**2*d222x0000TS0000 +
     &    4*CF*T*U*d222x0000TS0000 -
     &    2*NC*T*U*d222x0000TS0000 -
     &    2*NC*T**2*d222x0000TU0000 -
     &    2*NC*T*U*d222x0000TU0000 +
     &    4*CF*T*U*d222x0000US0000 -
     &    2*NC*T*U*d222x0000US0000 +
     &    8*CF*U**2*d222x0000US0000 -
     &    4*NC*U**2*d222x0000US0000 +
     &    12*CF*T**2*d223x0000TS0000 -
     &    6*NC*T**2*d223x0000TS0000 +
     &    8*CF*T*U*d223x0000TS0000 -
     &    4*NC*T*U*d223x0000TS0000 -
     &    4*NC*T**2*d223x0000TU0000 -
     &    4*NC*T*U*d223x0000TU0000 +
     &    8*CF*T*U*d223x0000US0000 -
     &    4*NC*T*U*d223x0000US0000 +
     &    12*CF*U**2*d223x0000US0000 -
     &    6*NC*U**2*d223x0000US0000 +
     &    4*CF*T**2*d23x0000TS0000 -
     &    2*NC*T**2*d23x0000TS0000 +
     &    4*CF*T*U*d23x0000TS0000 -
     &    2*NC*T*U*d23x0000TS0000 +
     &    2*NC*T**2*d23x0000TU0000 +
     &    2*NC*T*U*d23x0000TU0000 +
     &    4*CF*T*U*d23x0000US0000 -
     &    2*NC*T*U*d23x0000US0000 +
     &    4*CF*U**2*d23x0000US0000 -
     &    2*NC*U**2*d23x0000US0000 +
     &    4*CF*T**2*d233x0000TS0000 -
     &    2*NC*T**2*d233x0000TS0000 +
     &    4*CF*T*U*d233x0000TS0000 -
     &    2*NC*T*U*d233x0000TS0000 -
     &    2*NC*T**2*d233x0000TU0000 -
     &    2*NC*T*U*d233x0000TU0000 +
     &    4*CF*T*U*d233x0000US0000 -
     &    2*NC*T*U*d233x0000US0000 +
     &    4*CF*U**2*d233x0000US0000 -
     &    2*NC*U**2*d233x0000US0000 +
     &    4*CF*T**2*d3x0000TS0000 -
     &    2*NC*T**2*d3x0000TS0000 +
     &    8*CF*T*U*d3x0000TS0000 -
     &    4*NC*T*U*d3x0000TS0000 +
     &    4*CF*U**2*d3x0000TS0000 -
     &    2*NC*U**2*d3x0000TS0000 -
     &    2*NC*T*U*d3x0000TU0000 -
     &    2*NC*U**2*d3x0000TU0000 -
     &    2*(-2*CF + NC)*(T + U)**2*
     &     d3x0000US0000)

      call setDeps(-1, S, T, U,
     &     d0x0000TS0000, d0x0000TU0000, d0x0000US0000,
     &     d00x0000TS0000, d00x0000TU0000, d00x0000US0000,
     &     d001x0000TS0000, d001x0000TU0000, d001x0000US0000,
     &     d002x0000TS0000, d002x0000TU0000, d002x0000US0000,
     &     d003x0000TS0000, d003x0000TU0000, d003x0000US0000,
     &     d1x0000TS0000, d1x0000TU0000, d1x0000US0000,
     &     d11x0000TS0000, d11x0000TU0000, d11x0000US0000,
     &     d112x0000TS0000, d112x0000TU0000, d112x0000US0000,
     &     d113x0000TS0000, d113x0000TU0000, d113x0000US0000,
     &     d12x0000TS0000, d12x0000TU0000, d12x0000US0000,
     &     d122x0000TS0000, d122x0000TU0000, d122x0000US0000,
     &     d123x0000TS0000, d123x0000TU0000, d123x0000US0000,
     &     d13x0000TS0000, d13x0000TU0000, d13x0000US0000,
     &     d133x0000TS0000, d133x0000US0000, d2x0000TS0000,
     &     d2x0000TU0000, d2x0000US0000, d22x0000TS0000,
     &     d22x0000TU0000, d22x0000US0000, d222x0000TS0000,
     &     d222x0000TU0000, d222x0000US0000, d223x0000TS0000,
     &     d223x0000TU0000, d223x0000US0000, d23x0000TS0000,
     &     d23x0000TU0000, d23x0000US0000, d233x0000TS0000,
     &     d233x0000US0000, d3x0000TS0000, d3x0000TU0000,
     &     d3x0000US0000, d33x0000TS0000, d33x0000TU0000,
     &     d33x0000US0000, d133x0000TU0000, d233x0000TU0000)
      con1 = 32*em_alpha*st_alpha*CF*NC*Pi**2*Q**2*
     &  (-2*(2*CF - NC)*(T + U)*(T + 2*U)*
     &     d0x0000TS0000 +
     &    NC*(T**2 + 3*T*U + 4*U**2)*
     &     d0x0000TU0000 -
     &    8*CF*T**2*d0x0000US0000 +
     &    4*NC*T**2*d0x0000US0000 -
     &    12*CF*T*U*d0x0000US0000 +
     &    6*NC*T*U*d0x0000US0000 -
     &    4*CF*U**2*d0x0000US0000 +
     &    2*NC*U**2*d0x0000US0000 -
     &    24*CF*T*d00x0000TS0000 +
     &    12*NC*T*d00x0000TS0000 -
     &    48*CF*U*d00x0000TS0000 +
     &    24*NC*U*d00x0000TS0000 -
     &    12*NC*T*d00x0000TU0000 -
     &    12*NC*U*d00x0000TU0000 -
     &    48*CF*T*d00x0000US0000 +
     &    24*NC*T*d00x0000US0000 -
     &    24*CF*U*d00x0000US0000 +
     &    12*NC*U*d00x0000US0000 -
     &    24*CF*T*d001x0000TS0000 +
     &    12*NC*T*d001x0000TS0000 -
     &    24*CF*U*d001x0000TS0000 +
     &    12*NC*U*d001x0000TS0000 -
     &    24*CF*T*d001x0000US0000 +
     &    12*NC*T*d001x0000US0000 -
     &    24*CF*U*d001x0000US0000 +
     &    12*NC*U*d001x0000US0000 -
     &    56*CF*T*d002x0000TS0000 +
     &    28*NC*T*d002x0000TS0000 -
     &    32*CF*U*d002x0000TS0000 +
     &    16*NC*U*d002x0000TS0000 -
     &    4*NC*T*d002x0000TU0000 -
     &    4*NC*U*d002x0000TU0000 -
     &    32*CF*T*d002x0000US0000 +
     &    16*NC*T*d002x0000US0000 -
     &    56*CF*U*d002x0000US0000 +
     &    28*NC*U*d002x0000US0000 -
     &    24*CF*T*d003x0000TS0000 +
     &    12*NC*T*d003x0000TS0000 -
     &    24*CF*U*d003x0000TS0000 +
     &    12*NC*U*d003x0000TS0000 -
     &    4*NC*T*d003x0000TU0000 -
     &    4*NC*U*d003x0000TU0000 -
     &    24*CF*T*d003x0000US0000 +
     &    12*NC*T*d003x0000US0000 -
     &    24*CF*U*d003x0000US0000 +
     &    12*NC*U*d003x0000US0000 -
     &    8*CF*T**2*d1x0000TS0000 +
     &    4*NC*T**2*d1x0000TS0000 -
     &    20*CF*T*U*d1x0000TS0000 +
     &    10*NC*T*U*d1x0000TS0000 -
     &    12*CF*U**2*d1x0000TS0000 +
     &    6*NC*U**2*d1x0000TS0000 -
     &    2*NC*T**2*d1x0000TU0000 +
     &    4*NC*U**2*d1x0000TU0000 -
     &    12*CF*T**2*d1x0000US0000 +
     &    6*NC*T**2*d1x0000US0000 -
     &    20*CF*T*U*d1x0000US0000 +
     &    10*NC*T*U*d1x0000US0000 -
     &    8*CF*U**2*d1x0000US0000 +
     &    4*NC*U**2*d1x0000US0000 -
     &    4*CF*T**2*d11x0000TS0000 +
     &    2*NC*T**2*d11x0000TS0000 -
     &    8*CF*T*U*d11x0000TS0000 +
     &    4*NC*T*U*d11x0000TS0000 -
     &    4*CF*U**2*d11x0000TS0000 +
     &    2*NC*U**2*d11x0000TS0000 +
     &    NC*T**2*d11x0000TU0000 +
     &    NC*U**2*d11x0000TU0000 -
     &    4*CF*T**2*d11x0000US0000 +
     &    2*NC*T**2*d11x0000US0000 -
     &    8*CF*T*U*d11x0000US0000 +
     &    4*NC*T*U*d11x0000US0000 -
     &    4*CF*U**2*d11x0000US0000 +
     &    2*NC*U**2*d11x0000US0000 -
     &    8*CF*T**2*d112x0000TS0000 +
     &    4*NC*T**2*d112x0000TS0000 -
     &    8*CF*T*U*d112x0000TS0000 +
     &    4*NC*T*U*d112x0000TS0000 -
     &    NC*T**2*d112x0000TU0000 -
     &    3*NC*T*U*d112x0000TU0000 -
     &    8*CF*T*U*d112x0000US0000 +
     &    4*NC*T*U*d112x0000US0000 -
     &    8*CF*U**2*d112x0000US0000 +
     &    4*NC*U**2*d112x0000US0000 -
     &    8*CF*T**2*d113x0000TS0000 +
     &    4*NC*T**2*d113x0000TS0000 -
     &    16*CF*T*U*d113x0000TS0000 +
     &    8*NC*T*U*d113x0000TS0000 -
     &    8*CF*U**2*d113x0000TS0000 +
     &    4*NC*U**2*d113x0000TS0000 -
     &    3*NC*T*U*d113x0000TU0000 -
     &    NC*U**2*d113x0000TU0000 -
     &    8*CF*T**2*d113x0000US0000 +
     &    4*NC*T**2*d113x0000US0000 -
     &    16*CF*T*U*d113x0000US0000 +
     &    8*NC*T*U*d113x0000US0000 -
     &    8*CF*U**2*d113x0000US0000 +
     &    4*NC*U**2*d113x0000US0000 -
     &    20*CF*T**2*d12x0000TS0000 +
     &    10*NC*T**2*d12x0000TS0000 -
     &    24*CF*T*U*d12x0000TS0000 +
     &    12*NC*T*U*d12x0000TS0000 -
     &    4*CF*U**2*d12x0000TS0000 +
     &    2*NC*U**2*d12x0000TS0000 -
     &    2*NC*T**2*d12x0000TU0000 -
     &    8*NC*T*U*d12x0000TU0000 -
     &    4*CF*T**2*d12x0000US0000 +
     &    2*NC*T**2*d12x0000US0000 -
     &    24*CF*T*U*d12x0000US0000 +
     &    12*NC*T*U*d12x0000US0000 -
     &    20*CF*U**2*d12x0000US0000 +
     &    10*NC*U**2*d12x0000US0000 -
     &    20*CF*T**2*d122x0000TS0000 +
     &    10*NC*T**2*d122x0000TS0000 -
     &    12*CF*T*U*d122x0000TS0000 +
     &    6*NC*T*U*d122x0000TS0000 -
     &    4*NC*T*U*d122x0000TU0000 -
     &    12*CF*T*U*d122x0000US0000 +
     &    6*NC*T*U*d122x0000US0000 -
     &    20*CF*U**2*d122x0000US0000 +
     &    10*NC*U**2*d122x0000US0000 -
     &    28*CF*T**2*d123x0000TS0000 +
     &    14*NC*T**2*d123x0000TS0000 -
     &    40*CF*T*U*d123x0000TS0000 +
     &    20*NC*T*U*d123x0000TS0000 -
     &    12*CF*U**2*d123x0000TS0000 +
     &    6*NC*U**2*d123x0000TS0000 +
     &    NC*T**2*d123x0000TU0000 -
     &    6*NC*T*U*d123x0000TU0000 -
     &    3*NC*U**2*d123x0000TU0000 -
     &    12*CF*T**2*d123x0000US0000 +
     &    6*NC*T**2*d123x0000US0000 -
     &    40*CF*T*U*d123x0000US0000 +
     &    20*NC*T*U*d123x0000US0000 -
     &    28*CF*U**2*d123x0000US0000 +
     &    14*NC*U**2*d123x0000US0000 -
     &    12*CF*T**2*d13x0000TS0000 +
     &    6*NC*T**2*d13x0000TS0000 -
     &    32*CF*T*U*d13x0000TS0000 +
     &    16*NC*T*U*d13x0000TS0000 -
     &    20*CF*U**2*d13x0000TS0000 +
     &    10*NC*U**2*d13x0000TS0000 +
     &    2*NC*T**2*d13x0000TU0000 +
     &    2*NC*T*U*d13x0000TU0000 +
     &    2*NC*U**2*d13x0000TU0000 -
     &    20*CF*T**2*d13x0000US0000 +
     &    10*NC*T**2*d13x0000US0000 -
     &    32*CF*T*U*d13x0000US0000 +
     &    16*NC*T*U*d13x0000US0000 -
     &    12*CF*U**2*d13x0000US0000 +
     &    6*NC*U**2*d13x0000US0000 -
     &    8*CF*T**2*d133x0000TS0000 +
     &    4*NC*T**2*d133x0000TS0000 -
     &    16*CF*T*U*d133x0000TS0000 +
     &    8*NC*T*U*d133x0000TS0000 -
     &    8*CF*U**2*d133x0000TS0000 +
     &    4*NC*U**2*d133x0000TS0000 -
     &    2*NC*T*U*d133x0000TU0000 -
     &    2*NC*U**2*d133x0000TU0000 -
     &    8*CF*T**2*d133x0000US0000 +
     &    4*NC*T**2*d133x0000US0000 -
     &    16*CF*T*U*d133x0000US0000 +
     &    8*NC*T*U*d133x0000US0000 -
     &    8*CF*U**2*d133x0000US0000 +
     &    4*NC*U**2*d133x0000US0000 -
     &    12*CF*T**2*d2x0000TS0000 +
     &    6*NC*T**2*d2x0000TS0000 -
     &    20*CF*T*U*d2x0000TS0000 +
     &    10*NC*T*U*d2x0000TS0000 -
     &    8*CF*U**2*d2x0000TS0000 +
     &    4*NC*U**2*d2x0000TS0000 +
     &    NC*T**2*d2x0000TU0000 +
     &    NC*T*U*d2x0000TU0000 +
     &    4*NC*U**2*d2x0000TU0000 -
     &    8*CF*T**2*d2x0000US0000 +
     &    4*NC*T**2*d2x0000US0000 -
     &    20*CF*T*U*d2x0000US0000 +
     &    10*NC*T*U*d2x0000US0000 -
     &    12*CF*U**2*d2x0000US0000 +
     &    6*NC*U**2*d2x0000US0000 -
     &    20*CF*T**2*d22x0000TS0000 +
     &    10*NC*T**2*d22x0000TS0000 -
     &    12*CF*T*U*d22x0000TS0000 +
     &    6*NC*T*U*d22x0000TS0000 +
     &    NC*T**2*d22x0000TU0000 -
     &    3*NC*T*U*d22x0000TU0000 -
     &    12*CF*T*U*d22x0000US0000 +
     &    6*NC*T*U*d22x0000US0000 -
     &    20*CF*U**2*d22x0000US0000 +
     &    10*NC*U**2*d22x0000US0000 -
     &    12*CF*T**2*d222x0000TS0000 +
     &    6*NC*T**2*d222x0000TS0000 -
     &    4*CF*T*U*d222x0000TS0000 +
     &    2*NC*T*U*d222x0000TS0000 +
     &    NC*T**2*d222x0000TU0000 -
     &    NC*T*U*d222x0000TU0000 -
     &    4*CF*T*U*d222x0000US0000 +
     &    2*NC*T*U*d222x0000US0000 -
     &    12*CF*U**2*d222x0000US0000 +
     &    6*NC*U**2*d222x0000US0000 -
     &    20*CF*T**2*d223x0000TS0000 +
     &    10*NC*T**2*d223x0000TS0000 -
     &    12*CF*T*U*d223x0000TS0000 +
     &    6*NC*T*U*d223x0000TS0000 +
     &    3*NC*T**2*d223x0000TU0000 +
     &    NC*T*U*d223x0000TU0000 -
     &    12*CF*T*U*d223x0000US0000 +
     &    6*NC*T*U*d223x0000US0000 -
     &    20*CF*U**2*d223x0000US0000 +
     &    10*NC*U**2*d223x0000US0000 -
     &    20*CF*T**2*d23x0000TS0000 +
     &    10*NC*T**2*d23x0000TS0000 -
     &    24*CF*T*U*d23x0000TS0000 +
     &    12*NC*T*U*d23x0000TS0000 -
     &    4*CF*U**2*d23x0000TS0000 +
     &    2*NC*U**2*d23x0000TS0000 +
     &    3*NC*T**2*d23x0000TU0000 +
     &    3*NC*T*U*d23x0000TU0000 +
     &    2*NC*U**2*d23x0000TU0000 -
     &    4*CF*T**2*d23x0000US0000 +
     &    2*NC*T**2*d23x0000US0000 -
     &    24*CF*T*U*d23x0000US0000 +
     &    12*NC*T*U*d23x0000US0000 -
     &    20*CF*U**2*d23x0000US0000 +
     &    10*NC*U**2*d23x0000US0000 -
     &    8*CF*T**2*d233x0000TS0000 +
     &    4*NC*T**2*d233x0000TS0000 -
     &    8*CF*T*U*d233x0000TS0000 +
     &    4*NC*T*U*d233x0000TS0000 +
     &    2*NC*T**2*d233x0000TU0000 +
     &    2*NC*T*U*d233x0000TU0000 -
     &    8*CF*T*U*d233x0000US0000 +
     &    4*NC*T*U*d233x0000US0000 -
     &    8*CF*U**2*d233x0000US0000 +
     &    4*NC*U**2*d233x0000US0000 -
     &    8*CF*T**2*d3x0000TS0000 +
     &    4*NC*T**2*d3x0000TS0000 -
     &    20*CF*T*U*d3x0000TS0000 +
     &    10*NC*T*U*d3x0000TS0000 -
     &    12*CF*U**2*d3x0000TS0000 +
     &    6*NC*U**2*d3x0000TS0000 +
     &    2*NC*T**2*d3x0000TU0000 +
     &    6*NC*T*U*d3x0000TU0000 +
     &    6*NC*U**2*d3x0000TU0000 -
     &    12*CF*T**2*d3x0000US0000 +
     &    6*NC*T**2*d3x0000US0000 -
     &    20*CF*T*U*d3x0000US0000 +
     &    10*NC*T*U*d3x0000US0000 -
     &    8*CF*U**2*d3x0000US0000 +
     &    4*NC*U**2*d3x0000US0000 -
     &    4*CF*T**2*d33x0000TS0000 +
     &    2*NC*T**2*d33x0000TS0000 -
     &    8*CF*T*U*d33x0000TS0000 +
     &    4*NC*T*U*d33x0000TS0000 -
     &    4*CF*U**2*d33x0000TS0000 +
     &    2*NC*U**2*d33x0000TS0000 +
     &    NC*T**2*d33x0000TU0000 +
     &    3*NC*T*U*d33x0000TU0000 +
     &    2*NC*U**2*d33x0000TU0000 +
     &    2*(-2*CF + NC)*(T + U)**2*
     &     d33x0000US0000)

      call setDeps(0, S, T, U,
     &     d0x0000TS0000, d0x0000TU0000, d0x0000US0000,
     &     d00x0000TS0000, d00x0000TU0000, d00x0000US0000,
     &     d001x0000TS0000, d001x0000TU0000, d001x0000US0000,
     &     d002x0000TS0000, d002x0000TU0000, d002x0000US0000,
     &     d003x0000TS0000, d003x0000TU0000, d003x0000US0000,
     &     d1x0000TS0000, d1x0000TU0000, d1x0000US0000,
     &     d11x0000TS0000, d11x0000TU0000, d11x0000US0000,
     &     d112x0000TS0000, d112x0000TU0000, d112x0000US0000,
     &     d113x0000TS0000, d113x0000TU0000, d113x0000US0000,
     &     d12x0000TS0000, d12x0000TU0000, d12x0000US0000,
     &     d122x0000TS0000, d122x0000TU0000, d122x0000US0000,
     &     d123x0000TS0000, d123x0000TU0000, d123x0000US0000,
     &     d13x0000TS0000, d13x0000TU0000, d13x0000US0000,
     &     d133x0000TS0000, d133x0000US0000, d2x0000TS0000,
     &     d2x0000TU0000, d2x0000US0000, d22x0000TS0000,
     &     d22x0000TU0000, d22x0000US0000, d222x0000TS0000,
     &     d222x0000TU0000, d222x0000US0000, d223x0000TS0000,
     &     d223x0000TU0000, d223x0000US0000, d23x0000TS0000,
     &     d23x0000TU0000, d23x0000US0000, d233x0000TS0000,
     &     d233x0000US0000, d3x0000TS0000, d3x0000TU0000,
     &     d3x0000US0000, d33x0000TS0000, d33x0000TU0000,
     &     d33x0000US0000, d133x0000TU0000, d233x0000TU0000)
      con0 = 32*em_alpha*st_alpha*CF*NC*Pi**2*Q**2*
     &  (2*(2*CF - NC)*U*(T + U)*
     &     d0x0000TS0000 -
     &    NC*(T**2 + T*U + 2*U**2)*
     &     d0x0000TU0000 +
     &    4*CF*T**2*d0x0000US0000 -
     &    2*NC*T**2*d0x0000US0000 +
     &    4*CF*T*U*d0x0000US0000 -
     &    2*NC*T*U*d0x0000US0000 +
     &    8*CF*U*d00x0000TS0000 -
     &    4*NC*U*d00x0000TS0000 +
     &    8*NC*T*d00x0000TU0000 +
     &    10*NC*U*d00x0000TU0000 +
     &    8*CF*T*d00x0000US0000 -
     &    4*NC*T*d00x0000US0000 +
     &    8*CF*T*d001x0000TS0000 -
     &    4*NC*T*d001x0000TS0000 +
     &    8*CF*U*d001x0000TS0000 -
     &    4*NC*U*d001x0000TS0000 -
     &    2*NC*T*d001x0000TU0000 +
     &    2*NC*U*d001x0000TU0000 +
     &    8*CF*T*d001x0000US0000 -
     &    4*NC*T*d001x0000US0000 +
     &    8*CF*U*d001x0000US0000 -
     &    4*NC*U*d001x0000US0000 -
     &    8*CF*U*d002x0000TS0000 +
     &    4*NC*U*d002x0000TS0000 +
     &    6*NC*T*d002x0000TU0000 +
     &    10*NC*U*d002x0000TU0000 -
     &    8*CF*T*d002x0000US0000 +
     &    4*NC*T*d002x0000US0000 +
     &    8*CF*T*d003x0000TS0000 -
     &    4*NC*T*d003x0000TS0000 +
     &    8*CF*U*d003x0000TS0000 -
     &    4*NC*U*d003x0000TS0000 +
     &    8*NC*T*d003x0000TU0000 +
     &    8*NC*U*d003x0000TU0000 +
     &    8*CF*T*d003x0000US0000 -
     &    4*NC*T*d003x0000US0000 +
     &    8*CF*U*d003x0000US0000 -
     &    4*NC*U*d003x0000US0000 +
     &    4*CF*T**2*d1x0000TS0000 -
     &    2*NC*T**2*d1x0000TS0000 +
     &    12*CF*T*U*d1x0000TS0000 -
     &    6*NC*T*U*d1x0000TS0000 +
     &    8*CF*U**2*d1x0000TS0000 -
     &    4*NC*U**2*d1x0000TS0000 +
     &    NC*T**2*d1x0000TU0000 +
     &    2*NC*T*U*d1x0000TU0000 -
     &    NC*U**2*d1x0000TU0000 +
     &    8*CF*T**2*d1x0000US0000 -
     &    4*NC*T**2*d1x0000US0000 +
     &    12*CF*T*U*d1x0000US0000 -
     &    6*NC*T*U*d1x0000US0000 +
     &    4*CF*U**2*d1x0000US0000 -
     &    2*NC*U**2*d1x0000US0000 +
     &    4*CF*T**2*d11x0000TS0000 -
     &    2*NC*T**2*d11x0000TS0000 +
     &    8*CF*T*U*d11x0000TS0000 -
     &    4*NC*T*U*d11x0000TS0000 +
     &    4*CF*U**2*d11x0000TS0000 -
     &    2*NC*U**2*d11x0000TS0000 +
     &    2*NC*T*U*d11x0000TU0000 +
     &    4*CF*T**2*d11x0000US0000 -
     &    2*NC*T**2*d11x0000US0000 +
     &    8*CF*T*U*d11x0000US0000 -
     &    4*NC*T*U*d11x0000US0000 +
     &    4*CF*U**2*d11x0000US0000 -
     &    2*NC*U**2*d11x0000US0000 +
     &    4*CF*T**2*d112x0000TS0000 -
     &    2*NC*T**2*d112x0000TS0000 +
     &    4*CF*T*U*d112x0000TS0000 -
     &    2*NC*T*U*d112x0000TS0000 +
     &    NC*T**2*d112x0000TU0000 +
     &    3*NC*T*U*d112x0000TU0000 +
     &    4*CF*T*U*d112x0000US0000 -
     &    2*NC*T*U*d112x0000US0000 +
     &    4*CF*U**2*d112x0000US0000 -
     &    2*NC*U**2*d112x0000US0000 +
     &    4*CF*T**2*d113x0000TS0000 -
     &    2*NC*T**2*d113x0000TS0000 +
     &    8*CF*T*U*d113x0000TS0000 -
     &    4*NC*T*U*d113x0000TS0000 +
     &    4*CF*U**2*d113x0000TS0000 -
     &    2*NC*U**2*d113x0000TS0000 +
     &    3*NC*T*U*d113x0000TU0000 +
     &    NC*U**2*d113x0000TU0000 +
     &    4*CF*T**2*d113x0000US0000 -
     &    2*NC*T**2*d113x0000US0000 +
     &    8*CF*T*U*d113x0000US0000 -
     &    4*NC*T*U*d113x0000US0000 +
     &    4*CF*U**2*d113x0000US0000 -
     &    2*NC*U**2*d113x0000US0000 +
     &    8*CF*T**2*d12x0000TS0000 -
     &    4*NC*T**2*d12x0000TS0000 +
     &    12*CF*T*U*d12x0000TS0000 -
     &    6*NC*T*U*d12x0000TS0000 +
     &    4*CF*U**2*d12x0000TS0000 -
     &    2*NC*U**2*d12x0000TS0000 +
     &    3*NC*T**2*d12x0000TU0000 +
     &    10*NC*T*U*d12x0000TU0000 +
     &    NC*U**2*d12x0000TU0000 +
     &    4*CF*T**2*d12x0000US0000 -
     &    2*NC*T**2*d12x0000US0000 +
     &    12*CF*T*U*d12x0000US0000 -
     &    6*NC*T*U*d12x0000US0000 +
     &    8*CF*U**2*d12x0000US0000 -
     &    4*NC*U**2*d12x0000US0000 +
     &    4*CF*T**2*d122x0000TS0000 -
     &    2*NC*T**2*d122x0000TS0000 +
     &    2*NC*T**2*d122x0000TU0000 +
     &    6*NC*T*U*d122x0000TU0000 +
     &    4*CF*U**2*d122x0000US0000 -
     &    2*NC*U**2*d122x0000US0000 +
     &    8*CF*T**2*d123x0000TS0000 -
     &    4*NC*T**2*d123x0000TS0000 +
     &    12*CF*T*U*d123x0000TS0000 -
     &    6*NC*T*U*d123x0000TS0000 +
     &    4*CF*U**2*d123x0000TS0000 -
     &    2*NC*U**2*d123x0000TS0000 +
     &    NC*T**2*d123x0000TU0000 +
     &    6*NC*T*U*d123x0000TU0000 +
     &    NC*U**2*d123x0000TU0000 +
     &    4*CF*T**2*d123x0000US0000 -
     &    2*NC*T**2*d123x0000US0000 +
     &    12*CF*T*U*d123x0000US0000 -
     &    6*NC*T*U*d123x0000US0000 +
     &    8*CF*U**2*d123x0000US0000 -
     &    4*NC*U**2*d123x0000US0000 +
     &    8*CF*T**2*d13x0000TS0000 -
     &    4*NC*T**2*d13x0000TS0000 +
     &    20*CF*T*U*d13x0000TS0000 -
     &    10*NC*T*U*d13x0000TS0000 +
     &    12*CF*U**2*d13x0000TS0000 -
     &    6*NC*U**2*d13x0000TS0000 -
     &    NC*T**2*d13x0000TU0000 -
     &    NC*U**2*d13x0000TU0000 +
     &    12*CF*T**2*d13x0000US0000 -
     &    6*NC*T**2*d13x0000US0000 +
     &    20*CF*T*U*d13x0000US0000 -
     &    10*NC*T*U*d13x0000US0000 +
     &    8*CF*U**2*d13x0000US0000 -
     &    4*NC*U**2*d13x0000US0000 +
     &    4*CF*T**2*d133x0000TS0000 -
     &    2*NC*T**2*d133x0000TS0000 +
     &    8*CF*T*U*d133x0000TS0000 -
     &    4*NC*T*U*d133x0000TS0000 +
     &    4*CF*U**2*d133x0000TS0000 -
     &    2*NC*U**2*d133x0000TS0000 +
     &    4*CF*T**2*d133x0000US0000 -
     &    2*NC*T**2*d133x0000US0000 +
     &    8*CF*T*U*d133x0000US0000 -
     &    4*NC*T*U*d133x0000US0000 +
     &    4*CF*U**2*d133x0000US0000 -
     &    2*NC*U**2*d133x0000US0000 +
     &    4*CF*T*U*d2x0000TS0000 -
     &    2*NC*T*U*d2x0000TS0000 +
     &    4*CF*U**2*d2x0000TS0000 -
     &    2*NC*U**2*d2x0000TS0000 -
     &    NC*T**2*d2x0000TU0000 +
     &    NC*T*U*d2x0000TU0000 -
     &    2*NC*U**2*d2x0000TU0000 +
     &    4*CF*T**2*d2x0000US0000 -
     &    2*NC*T**2*d2x0000US0000 +
     &    4*CF*T*U*d2x0000US0000 -
     &    2*NC*T*U*d2x0000US0000 -
     &    4*CF*T*U*d22x0000TS0000 +
     &    2*NC*T*U*d22x0000TS0000 +
     &    NC*T**2*d22x0000TU0000 +
     &    5*NC*T*U*d22x0000TU0000 -
     &    4*CF*T*U*d22x0000US0000 +
     &    2*NC*T*U*d22x0000US0000 -
     &    4*CF*T*U*d222x0000TS0000 +
     &    2*NC*T*U*d222x0000TS0000 +
     &    NC*T**2*d222x0000TU0000 +
     &    3*NC*T*U*d222x0000TU0000 -
     &    4*CF*T*U*d222x0000US0000 +
     &    2*NC*T*U*d222x0000US0000 +
     &    4*CF*T**2*d223x0000TS0000 -
     &    2*NC*T**2*d223x0000TS0000 +
     &    NC*T**2*d223x0000TU0000 +
     &    3*NC*T*U*d223x0000TU0000 +
     &    4*CF*U**2*d223x0000US0000 -
     &    2*NC*U**2*d223x0000US0000 +
     &    8*CF*T**2*d23x0000TS0000 -
     &    4*NC*T**2*d23x0000TS0000 +
     &    12*CF*T*U*d23x0000TS0000 -
     &    6*NC*T*U*d23x0000TS0000 +
     &    4*CF*U**2*d23x0000TS0000 -
     &    2*NC*U**2*d23x0000TS0000 -
     &    NC*T**2*d23x0000TU0000 -
     &    NC*T*U*d23x0000TU0000 -
     &    2*NC*U**2*d23x0000TU0000 +
     &    4*CF*T**2*d23x0000US0000 -
     &    2*NC*T**2*d23x0000US0000 +
     &    12*CF*T*U*d23x0000US0000 -
     &    6*NC*T*U*d23x0000US0000 +
     &    8*CF*U**2*d23x0000US0000 -
     &    4*NC*U**2*d23x0000US0000 +
     &    4*CF*T**2*d233x0000TS0000 -
     &    2*NC*T**2*d233x0000TS0000 +
     &    4*CF*T*U*d233x0000TS0000 -
     &    2*NC*T*U*d233x0000TS0000 +
     &    4*CF*T*U*d233x0000US0000 -
     &    2*NC*T*U*d233x0000US0000 +
     &    4*CF*U**2*d233x0000US0000 -
     &    2*NC*U**2*d233x0000US0000 +
     &    4*CF*T**2*d3x0000TS0000 -
     &    2*NC*T**2*d3x0000TS0000 +
     &    12*CF*T*U*d3x0000TS0000 -
     &    6*NC*T*U*d3x0000TS0000 +
     &    8*CF*U**2*d3x0000TS0000 -
     &    4*NC*U**2*d3x0000TS0000 -
     &    2*NC*T**2*d3x0000TU0000 -
     &    4*NC*T*U*d3x0000TU0000 -
     &    4*NC*U**2*d3x0000TU0000 +
     &    8*CF*T**2*d3x0000US0000 -
     &    4*NC*T**2*d3x0000US0000 +
     &    12*CF*T*U*d3x0000US0000 -
     &    6*NC*T*U*d3x0000US0000 +
     &    4*CF*U**2*d3x0000US0000 -
     &    2*NC*U**2*d3x0000US0000 +
     &    4*CF*T**2*d33x0000TS0000 -
     &    2*NC*T**2*d33x0000TS0000 +
     &    8*CF*T*U*d33x0000TS0000 -
     &    4*NC*T*U*d33x0000TS0000 +
     &    4*CF*U**2*d33x0000TS0000 -
     &    2*NC*U**2*d33x0000TS0000 -
     &    NC*T**2*d33x0000TU0000 -
     &    3*NC*T*U*d33x0000TU0000 -
     &    2*NC*U**2*d33x0000TU0000 -
     &    2*(-2*CF + NC)*(T + U)**2*
     &     d33x0000US0000)

      qqbAgBoxC = con0 + con1 + con2

      qqbAgBox = 2 * DBLE(qqbAgBoxC)

      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccc
ccc subroutine to set the values of the PV functions
ccc no idea how to solve this nicer w.o. rewriting
ccc my MEs...

      subroutine setDeps(eps, S, T, U,
     &     d0x0000TS0000, d0x0000TU0000, d0x0000US0000,
     &     d00x0000TS0000, d00x0000TU0000, d00x0000US0000,
     &     d001x0000TS0000, d001x0000TU0000, d001x0000US0000,
     &     d002x0000TS0000, d002x0000TU0000, d002x0000US0000,
     &     d003x0000TS0000, d003x0000TU0000, d003x0000US0000,
     &     d1x0000TS0000, d1x0000TU0000, d1x0000US0000,
     &     d11x0000TS0000, d11x0000TU0000, d11x0000US0000,
     &     d112x0000TS0000, d112x0000TU0000, d112x0000US0000,
     &     d113x0000TS0000, d113x0000TU0000, d113x0000US0000,
     &     d12x0000TS0000, d12x0000TU0000, d12x0000US0000,
     &     d122x0000TS0000, d122x0000TU0000, d122x0000US0000,
     &     d123x0000TS0000, d123x0000TU0000, d123x0000US0000,
     &     d13x0000TS0000, d13x0000TU0000, d13x0000US0000,
     &     d133x0000TS0000, d133x0000US0000, d2x0000TS0000,
     &     d2x0000TU0000, d2x0000US0000, d22x0000TS0000,
     &     d22x0000TU0000, d22x0000US0000, d222x0000TS0000,
     &     d222x0000TU0000, d222x0000US0000, d223x0000TS0000,
     &     d223x0000TU0000, d223x0000US0000, d23x0000TS0000,
     &     d23x0000TU0000, d23x0000US0000, d233x0000TS0000,
     &     d233x0000US0000, d3x0000TS0000, d3x0000TU0000,
     &     d3x0000US0000, d33x0000TS0000, d33x0000TU0000,
     &     d33x0000US0000, d133x0000TU0000, d233x0000TU0000)
      implicit none
#include "looptools.h"
      double complex d0x0000TS0000, d0x0000TU0000, d0x0000US0000,
     &     d00x0000TS0000, d00x0000TU0000, d00x0000US0000,
     &     d001x0000TS0000, d001x0000TU0000, d001x0000US0000,
     &     d002x0000TS0000, d002x0000TU0000, d002x0000US0000,
     &     d003x0000TS0000, d003x0000TU0000, d003x0000US0000,
     &     d1x0000TS0000, d1x0000TU0000, d1x0000US0000,
     &     d11x0000TS0000, d11x0000TU0000, d11x0000US0000,
     &     d112x0000TS0000, d112x0000TU0000, d112x0000US0000,
     &     d113x0000TS0000, d113x0000TU0000, d113x0000US0000,
     &     d12x0000TS0000, d12x0000TU0000, d12x0000US0000,
     &     d122x0000TS0000, d122x0000TU0000, d122x0000US0000,
     &     d123x0000TS0000, d123x0000TU0000, d123x0000US0000,
     &     d13x0000TS0000, d13x0000TU0000, d13x0000US0000,
     &     d133x0000TS0000, d133x0000US0000, d2x0000TS0000,
     &     d2x0000TU0000, d2x0000US0000, d22x0000TS0000,
     &     d22x0000TU0000, d22x0000US0000, d222x0000TS0000,
     &     d222x0000TU0000, d222x0000US0000, d223x0000TS0000,
     &     d223x0000TU0000, d223x0000US0000, d23x0000TS0000,
     &     d23x0000TU0000, d23x0000US0000, d233x0000TS0000,
     &     d233x0000US0000, d3x0000TS0000, d3x0000TU0000,
     &     d3x0000US0000, d33x0000TS0000, d33x0000TU0000,
     &     d33x0000US0000, d133x0000TU0000, d233x0000TU0000
      double precision S, T, U, deps
      integer eps

      deps = DBLE(eps)
      call setlambda(deps)

      d0x0000TS0000 = D0i(dd0,0.d0,0.d0,0.d0,0.d0,T,S,
     &     0.d0,0.d0,0.d0,0.d0)
      d0x0000TU0000 = D0i(dd0,0.d0,0.d0,0.d0,0.d0,T,U,
     &     0.d0,0.d0,0.d0,0.d0)
      d0x0000US0000 = D0i(dd0,0.d0,0.d0,0.d0,0.d0,U,S,
     &     0.d0,0.d0,0.d0,0.d0)
      d00x0000TS0000 = D0i(dd00,0.d0,0.d0,0.d0,0.d0,T,S,
     &     0.d0,0.d0,0.d0,0.d0)
      d00x0000TU0000 = D0i(dd00,0.d0,0.d0,0.d0,0.d0,T,U,
     &     0.d0,0.d0,0.d0,0.d0)
      d00x0000US0000 = D0i(dd00,0.d0,0.d0,0.d0,0.d0,U,S,
     &     0.d0,0.d0,0.d0,0.d0)
      d001x0000TS0000 = D0i(dd001,0.d0,0.d0,0.d0,0.d0,T,S,
     &     0.d0,0.d0,0.d0,0.d0)
      d001x0000TU0000 = D0i(dd001,0.d0,0.d0,0.d0,0.d0,T,U,
     &     0.d0,0.d0,0.d0,0.d0)
      d001x0000US0000 = D0i(dd001,0.d0,0.d0,0.d0,0.d0,U,S,
     &     0.d0,0.d0,0.d0,0.d0)
      d002x0000TS0000 = D0i(dd002,0.d0,0.d0,0.d0,0.d0,T,S,
     &     0.d0,0.d0,0.d0,0.d0)
      d002x0000TU0000 = D0i(dd002,0.d0,0.d0,0.d0,0.d0,T,U,
     &     0.d0,0.d0,0.d0,0.d0)
      d002x0000US0000 = D0i(dd002,0.d0,0.d0,0.d0,0.d0,U,S,
     &     0.d0,0.d0,0.d0,0.d0)
      d003x0000TS0000 = D0i(dd003,0.d0,0.d0,0.d0,0.d0,T,S,
     &     0.d0,0.d0,0.d0,0.d0)
      d003x0000TU0000 = D0i(dd003,0.d0,0.d0,0.d0,0.d0,T,U,
     &     0.d0,0.d0,0.d0,0.d0)
      d003x0000US0000 = D0i(dd003,0.d0,0.d0,0.d0,0.d0,U,S,
     &     0.d0,0.d0,0.d0,0.d0)
      d1x0000TS0000 = D0i(dd1,0.d0,0.d0,0.d0,0.d0,T,S,
     &     0.d0,0.d0,0.d0,0.d0)
      d1x0000TU0000 = D0i(dd1,0.d0,0.d0,0.d0,0.d0,T,U,
     &     0.d0,0.d0,0.d0,0.d0)
      d1x0000US0000 = D0i(dd1,0.d0,0.d0,0.d0,0.d0,U,S,
     &     0.d0,0.d0,0.d0,0.d0)
      d11x0000TS0000 = D0i(dd11,0.d0,0.d0,0.d0,0.d0,T,S,
     &     0.d0,0.d0,0.d0,0.d0)
      d11x0000TU0000 = D0i(dd11,0.d0,0.d0,0.d0,0.d0,T,U,
     &     0.d0,0.d0,0.d0,0.d0)
      d11x0000US0000 = D0i(dd11,0.d0,0.d0,0.d0,0.d0,U,S,
     &     0.d0,0.d0,0.d0,0.d0)
      d112x0000TS0000 = D0i(dd112,0.d0,0.d0,0.d0,0.d0,T,S,
     &     0.d0,0.d0,0.d0,0.d0)
      d112x0000TU0000 = D0i(dd112,0.d0,0.d0,0.d0,0.d0,T,U,
     &     0.d0,0.d0,0.d0,0.d0)
      d112x0000US0000 = D0i(dd112,0.d0,0.d0,0.d0,0.d0,U,S,
     &     0.d0,0.d0,0.d0,0.d0)
      d113x0000TS0000 = D0i(dd113,0.d0,0.d0,0.d0,0.d0,T,S,
     &     0.d0,0.d0,0.d0,0.d0)
      d113x0000TU0000 = D0i(dd113,0.d0,0.d0,0.d0,0.d0,T,U,
     &     0.d0,0.d0,0.d0,0.d0)
      d113x0000US0000 = D0i(dd113,0.d0,0.d0,0.d0,0.d0,U,S,
     &     0.d0,0.d0,0.d0,0.d0)
      d12x0000TS0000 = D0i(dd12,0.d0,0.d0,0.d0,0.d0,T,S,
     &     0.d0,0.d0,0.d0,0.d0)
      d12x0000TU0000 = D0i(dd12,0.d0,0.d0,0.d0,0.d0,T,U,
     &     0.d0,0.d0,0.d0,0.d0)
      d12x0000US0000 = D0i(dd12,0.d0,0.d0,0.d0,0.d0,U,S,
     &     0.d0,0.d0,0.d0,0.d0)
      d122x0000TS0000 = D0i(dd122,0.d0,0.d0,0.d0,0.d0,T,S,
     &     0.d0,0.d0,0.d0,0.d0)
      d122x0000TU0000 = D0i(dd122,0.d0,0.d0,0.d0,0.d0,T,U,
     &     0.d0,0.d0,0.d0,0.d0)
      d122x0000US0000 = D0i(dd122,0.d0,0.d0,0.d0,0.d0,U,S,
     &     0.d0,0.d0,0.d0,0.d0)
      d123x0000TS0000 = D0i(dd123,0.d0,0.d0,0.d0,0.d0,T,S,
     &     0.d0,0.d0,0.d0,0.d0)
      d123x0000TU0000 = D0i(dd123,0.d0,0.d0,0.d0,0.d0,T,U,
     &     0.d0,0.d0,0.d0,0.d0)
      d123x0000US0000 = D0i(dd123,0.d0,0.d0,0.d0,0.d0,U,S,
     &     0.d0,0.d0,0.d0,0.d0)
      d13x0000TS0000 = D0i(dd13,0.d0,0.d0,0.d0,0.d0,T,S,
     &     0.d0,0.d0,0.d0,0.d0)
      d13x0000TU0000 = D0i(dd13,0.d0,0.d0,0.d0,0.d0,T,U,
     &     0.d0,0.d0,0.d0,0.d0)
      d13x0000US0000 = D0i(dd13,0.d0,0.d0,0.d0,0.d0,U,S,
     &     0.d0,0.d0,0.d0,0.d0)
      d133x0000TS0000 = D0i(dd133,0.d0,0.d0,0.d0,0.d0,T,S,
     &     0.d0,0.d0,0.d0,0.d0)
      d133x0000US0000 = D0i(dd133,0.d0,0.d0,0.d0,0.d0,U,S,
     &     0.d0,0.d0,0.d0,0.d0)
      d2x0000TS0000 = D0i(dd2,0.d0,0.d0,0.d0,0.d0,T,S,
     &     0.d0,0.d0,0.d0,0.d0)
      d2x0000TU0000 = D0i(dd2,0.d0,0.d0,0.d0,0.d0,T,U,
     &     0.d0,0.d0,0.d0,0.d0)
      d2x0000US0000 = D0i(dd2,0.d0,0.d0,0.d0,0.d0,U,S,
     &     0.d0,0.d0,0.d0,0.d0)
      d22x0000TS0000 = D0i(dd22,0.d0,0.d0,0.d0,0.d0,T,S,
     &     0.d0,0.d0,0.d0,0.d0)
      d22x0000TU0000 = D0i(dd22,0.d0,0.d0,0.d0,0.d0,T,U,
     &     0.d0,0.d0,0.d0,0.d0)
      d22x0000US0000 = D0i(dd22,0.d0,0.d0,0.d0,0.d0,U,S,
     &     0.d0,0.d0,0.d0,0.d0)
      d222x0000TS0000 = D0i(dd222,0.d0,0.d0,0.d0,0.d0,T,S,
     &     0.d0,0.d0,0.d0,0.d0)
      d222x0000TU0000 = D0i(dd222,0.d0,0.d0,0.d0,0.d0,T,U,
     &     0.d0,0.d0,0.d0,0.d0)
      d222x0000US0000 = D0i(dd222,0.d0,0.d0,0.d0,0.d0,U,S,
     &     0.d0,0.d0,0.d0,0.d0)
      d223x0000TS0000 = D0i(dd223,0.d0,0.d0,0.d0,0.d0,T,S,
     &     0.d0,0.d0,0.d0,0.d0)
      d223x0000TU0000 = D0i(dd223,0.d0,0.d0,0.d0,0.d0,T,U,
     &     0.d0,0.d0,0.d0,0.d0)
      d223x0000US0000 = D0i(dd223,0.d0,0.d0,0.d0,0.d0,U,S,
     &     0.d0,0.d0,0.d0,0.d0)
      d23x0000TS0000 = D0i(dd23,0.d0,0.d0,0.d0,0.d0,T,S,
     &     0.d0,0.d0,0.d0,0.d0)
      d23x0000TU0000 = D0i(dd23,0.d0,0.d0,0.d0,0.d0,T,U,
     &     0.d0,0.d0,0.d0,0.d0)
      d23x0000US0000 = D0i(dd23,0.d0,0.d0,0.d0,0.d0,U,S,
     &     0.d0,0.d0,0.d0,0.d0)
      d233x0000TS0000 = D0i(dd233,0.d0,0.d0,0.d0,0.d0,T,S,
     &     0.d0,0.d0,0.d0,0.d0)
      d233x0000US0000 = D0i(dd233,0.d0,0.d0,0.d0,0.d0,U,S,
     &     0.d0,0.d0,0.d0,0.d0)
      d3x0000TS0000 = D0i(dd3,0.d0,0.d0,0.d0,0.d0,T,S,
     &     0.d0,0.d0,0.d0,0.d0)
      d3x0000TU0000 = D0i(dd3,0.d0,0.d0,0.d0,0.d0,T,U,
     &     0.d0,0.d0,0.d0,0.d0)
      d3x0000US0000 = D0i(dd3,0.d0,0.d0,0.d0,0.d0,U,S,
     &     0.d0,0.d0,0.d0,0.d0)
      d33x0000TS0000 = D0i(dd33,0.d0,0.d0,0.d0,0.d0,T,S,
     &     0.d0,0.d0,0.d0,0.d0)
      d33x0000TU0000 = D0i(dd33,0.d0,0.d0,0.d0,0.d0,T,U,
     &     0.d0,0.d0,0.d0,0.d0)
      d33x0000US0000 = D0i(dd33,0.d0,0.d0,0.d0,0.d0,U,S,
     &     0.d0,0.d0,0.d0,0.d0)
      d133x0000TU0000 = D0i(dd133,0.d0,0.d0,0.d0,0.d0,T,U,
     &     0.d0,0.d0,0.d0,0.d0)
      d233x0000TU0000 = D0i(dd233,0.d0,0.d0,0.d0,0.d0,T,U,
     &     0.d0,0.d0,0.d0,0.d0)

      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccc
ccc Counterterms
ccc This finite piece has to be taken into
ccc account if finite terms in the virtual
ccc corrections are generated by cancellation
ccc of the 1/eps_UV pole.

      double precision function qqbAgCT(Q, S, T, U)
      implicit none
      include 'pwhg_math.h'
      include 'pwhg_st.h'
      include 'pwhg_em.h'
      double precision Q, S, T, U
c     if CTs depend on EulerGamma - DLOG(4*Pi) then setdelta has to be
c     used in setvirtual!
c$$$      include 'constants.h'
c$$$
c$$$      qqbAgCT = (-64*em_alpha*st_alpha*CF*NC*(CF + NC)*Pi**2*Q**2*
c$$$     &    (-2*T*U + T**2*
c$$$     &       (-2 + EulerGamma - DLOG(4*Pi)) +
c$$$     &      U**2*(-2 + EulerGamma - DLOG(4*Pi))
c$$$     &      ))/(T*U)
      qqbAgCT = (-64*em_alpha*st_alpha*CF*NC*(CF + NC)*Pi**2*Q**2*
     &    (-2*T*U + T**2*(-2) + U**2*(-2)))/(T*U)

      qqbAgCT = 2 * qqbAgCT

      return
      end
