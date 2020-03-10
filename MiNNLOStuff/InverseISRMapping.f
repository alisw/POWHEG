c     inverse integral of spreading, denominator of F_corr of paper
      subroutine evaluubjakob(res)
      implicit none
      include 'brinclude.h'
      include 'pwhg_kn.h'
      include 'pwhg_math.h'
      include 'pwhg_pdf.h'
      real * 8 res
      real * 8 amp2uub,ebeam 
      real * 8 pb(0:3,nlegborn-1),sh, phint, fourPtsqOsh,x1b,x2b,pt,
     1     ptmax,mh2,yh,e,e2,pbsyst(0:3)
      real * 8 pdf1(-pdf_nparton:pdf_nparton),
     1     pdf2(-pdf_nparton:pdf_nparton)
      integer i

      call invISRmap(kn_pborn,nlegborn,pb)
      e2 =  kn_sbeams
      e = sqrt(e2)
      ebeam = e/2
      x1b = pb(0,1)/ebeam
      x2b = pb(0,2)/ebeam
      sh = x1b*x2b*kn_sbeams 

      brkn_xb1 = x1b
      brkn_xb2 = x2b
      brkn_sborn = sh
      brkn_pborn(:,:) = pb(:,:)

C     compute Phase space integral and divide it out 

c     pbsyst(0:3) is the colourless "boson" momentum
      pbsyst(:)=0d0
      do i=3,nlegborn
c     the sequence of colourless particles is the same in all flst_born(i,X) (X does not matter), so we pick X=1
         if (abs(flst_born(i,1)).gt.6) then
            pbsyst(:)=pbsyst(:) + kn_pborn(:,i)
         endif
      enddo
      pt = sqrt(pbsyst(1)**2+pbsyst(2)**2)
      
c this should be the same as the pt of the radiated parton
      if (abs(pt-sqrt(kn_pborn(1,nlegborn)**2+kn_pborn(2,nlegborn)**2)) .gt. 1d-7) then
          write(*,*) 'ERROR: in computing pt of colourless system!!!'
          stop
      endif

      fourPtsqOsh = 4d0*pt**2/sh
      call borndenomint(x1b,x2b,fourPtsqOsh,phint)
      phint = phint * pt / (4d0*pi**2) 

C     include here flux factor 
      if (phint .ne. 0d0) then
         res = 1d0/phint / (2d0*sh )
      else
         res = -1d10
      endif
      end


      subroutine evaluubsigma(res)
      implicit none
      include 'brinclude.h'
      include 'pwhg_kn.h'
      include 'pwhg_math.h'
      include 'pwhg_pdf.h'
      real * 8 res
      real * 8 amp2uub,ebeam 
      real * 8 pb(0:3,nlegborn-1),sh, phint, fourPtsqOsh,x1b,x2b,pt,
     1     ptmax,mh2,yh,e,e2
      real * 8 pdf1(-pdf_nparton:pdf_nparton),
     1     pdf2(-pdf_nparton:pdf_nparton),pbsyst(0:3)
      integer i

      call invISRmap(kn_pborn,nlegborn,pb)
      call uub_for_minnlo(pb,2,amp2uub)
      e2 =  kn_sbeams
      e = sqrt(e2)
      ebeam = e/2
      x1b = pb(0,1)/ebeam
      x2b = pb(0,2)/ebeam
      sh = x1b*x2b*kn_sbeams 
      call pdfcall(1,x1b,pdf1)
      call pdfcall(2,x2b,pdf2)
      res = amp2uub * pdf1(0)*pdf2(0) 

C     compute Phase space integral and divide it out 

c     pbsyst(0:3) is the colourless "boson" momentum
      pbsyst(:)=0d0
      do i=3,nlegborn
c     the sequence of colourless particles is the same in all flst_born(i,X) (X does not matter), so we pick X=1
         if (abs(flst_born(i,1)).gt.6) then
            pbsyst(:)=pbsyst(:) + kn_pborn(:,i)
         endif
      enddo
      pt = sqrt(pbsyst(1)**2+pbsyst(2)**2)
      
      if (abs(pt-sqrt(kn_pborn(1,nlegborn)**2+kn_pborn(2,nlegborn)**2)) .gt. 1d-7) then
          write(*,*) 'ERROR: in computing pt of colourless system!!!'
          stop
      endif
!      pt = sqrt(kn_pborn(1,3)**2+kn_pborn(2,3)**2)

      fourPtsqOsh = 4d0*pt**2/sh
      call borndenomint(x1b,x2b,fourPtsqOsh,phint)
      phint = phint * pt / (4d0*pi**2)

      yh = 0.5d0*log(((pb(0,1)+pb(0,2)) + (pb(3,1)+pb(3,2)))/
     1     ((pb(0,1)+pb(0,2)) - (pb(3,1)+pb(3,2))))
      mh2 = kn_sbeams*x1b*x2b

      ptmax= sqrt((((e2-mh2)/2/e)**2+mh2)/cosh(yh)**2 - mh2)
C     include here flux factor 
      res = res/phint / (2d0*sh ) / ptmax
      
      end


c     Inverse mapping from p(0:3,n) -> pb(0:3,n-1), following      
c     Sec 5.1.1 of arXiv:0709.2092      
c     The n-th momentum is interpreted as the radiated particle in ISR
      subroutine invISRmap(p,n,pb)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_kn.h'
c      real * 8 p(0:3,nlegrealexternal)
      integer n
      real * 8 p(0:3,n)
      real * 8 pb(0:3,n-1)
      real * 8 ktot(0:3),ktotl(0:3),prad(0:3)
      real * 8 betalvec(3), betal
      real * 8 betatvec(3), betat,mod
      real * 8 pCM(0:3,n)
      integer i,mu
      real * 8 dotp
      external dotp
      real * 8 ecm, xplus, xminus, y, csi, xplusbar, xminusbar,v,ebeam
      real * 8 zero
      
      ktot = p(:,1) + p(:,2) - p(:,n)
c     longitudinal boost
      betalvec(1) = 0d0
      betalvec(2) = 0d0
      betalvec(3) = 1d0
      betal = - ktot(3)/ktot(0)
      v = ktot(3)/ktot(0)
      call mboost(1,betalvec,betal,ktot,ktot)

c     transverse boost
      mod = sqrt(ktot(1)**2+ktot(2)**2)
      betatvec(1) = ktot(1)/mod
      betatvec(2) = ktot(2)/mod
      betatvec(3) = 0d0
      betat = -mod/ktot(0)
      call mboost(1,betatvec,betat,ktot,ktot)

c     apply  BL^(-1) BT BL
      call mboost(n-3,betalvec,betal,p(:,3),pb(:,3))
      call mboost(n-3,betatvec,betat,pb(:,3),pb(:,3))
      call mboost(n-3,betalvec,-betal,pb(:,3),pb(:,3))
      ktot=0d0
      do i=3,n-1
         ktot=ktot+pb(:,i)
      enddo

      ebeam=sqrt(kn_sbeams)

c     use the fact that
c     xplusbar * Kplus + xminusbar * Kminus = ktot (after boosts)!
      xplusbar = (ktot(0)+ktot(3))/ebeam
      xminusbar = (ktot(0)-ktot(3))/ebeam

      pb(0,1) = ebeam/2*xplusbar
      pb(1,1) = 0
      pb(2,1) = 0
      pb(3,1) = ebeam/2*xplusbar

      pb(0,2) = ebeam/2*xminusbar
      pb(1,2) = 0
      pb(2,2) = 0
      pb(3,2) = -ebeam/2*xminusbar

      return
      
c     check momentum conservation
      ktot=ktot-pb(:,1)-pb(:,2)
      zero = 0d0
      do mu=0, 3
         zero = zero + ktot(mu)**2
      enddo
      if (zero.gt.1d-12) then
         write(*,*) 'VIOLATION MOM CONSER ', zero
      endif
      
      end
