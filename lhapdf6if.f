      subroutine genericpdfset(ndns)
      implicit none
      include 'pwhg_pdf.h'
      integer ndns
c wrap for pdfset; avoids subsequent
c calls to pdfset (you never know)
      integer, parameter:: maxsets = 10
      integer sets(maxsets),idset(maxsets),nsets
      real * 8 lam5set(maxsets),q2minset(maxsets)
      save lam5set,q2minset,sets,idset,nsets
      real * 8 genericxlambdL,genericxlambdNL,genericxlambdNNL,
     ,         alphasPDF,asmz,mz
      parameter (mz=91.1876d0)
      real * 8 lam5
      integer iord,iset
      common/cgenericpdf/lam5,iord,iset
      logical,save:: ini=.true.
      integer j
      if(ini) then
         nsets = 0
         iset = 0
         ini = .false.
      endif
      if(iset /= 0) then
         if(idset(iset) == ndns) return
      endif
      do j=1,nsets
         if(ndns == idset(j)) then
            iset = j
            exit
         endif
      enddo
      if(j == nsets+1) then
c... initalise set
         nsets = nsets+1
         if(nsets .gt. maxsets) then
            write(*,*) ' generipdfset: too many sets,'
            write(*,*) ' increase maxsets. Exiting ...'
            call exit(-1)
         endif
         iset = nsets
         idset(iset) = ndns
         call setlha6set(iset,ndns,iord,mz,asmz,pdf_q2min)
         write(*,*) ' check: alpha_s(Mz)=',asmz
         if(iord.eq.0) then
c better to use the NLO formula anyhow;
c we don't have LO alpha around
            lam5=genericxlambdNL(asmz,mz,5)
         elseif(iord.eq.1) then
            lam5=genericxlambdNL(asmz,mz,5)
         elseif(iord.eq.2) then
            lam5=genericxlambdNNL(asmz,mz,5)
         endif
         write(*,*) ' alpha_s order (0,1,2): ',iord
         write(*,*) ' Lambda 5 is ',lam5
         lam5set(iset) = lam5
         q2minset(iset) = pdf_q2min
      else
         lam5 = lam5set(iset)
         pdf_q2min = q2minset(iset)
      endif
      end

      function genericxlambdL(as,q,nf)
      implicit none
      real * 8 genericxlambdL,as,q
      integer nf
      real * 8 pi,b,t,xlt,ot,as0,as1
      parameter (pi=3.14159265358979312D0)
      b  = (33-2*nf)/pi/12
      t  = 1/b/as
    1 xlt = log(t)
      ot = t
c-----------------------------------------------------------
c Value and Derivative of alfa with respect to t
      as0  = 1/b/t
      as1  = - 1/b/t**2
      t  = (as-as0)/as1 + t
      if(abs(ot-t)/ot.gt.0.00000001d0)goto 1
      genericxlambdL = q/exp(t/2)
      return
      end

      function genericxlambdNL(as,q,nf)
      implicit none
      real * 8 genericxlambdNL,as,q
      integer nf
      real * 8 pi,b,bp,t,xlt,ot,as0,as1
      parameter (pi=3.14159265358979312D0)
      b  = (33-2*nf)/pi/12
      bp = (153 - 19*nf) / pi / 2 / (33 - 2*nf)
      t  = 1/b/as
    1 xlt = log(t)
      ot = t
c-----------------------------------------------------------
c Value and Derivative of alfa with respect to t
      as0  = 1/b/t - bp*xlt/(b*t)**2
      as1  = - 1/b/t**2 -bp/b**2*(1-2*xlt)/t**3
      t  = (as-as0)/as1 + t
      if(abs(ot-t)/ot.gt.0.00000001d0)goto 1
      genericxlambdNL = q/exp(t/2)
      return
      end

      function genericxlambdNNL(as,q,nf)
      implicit none
      real * 8 genericxlambdNNL,as,q
      integer nf
      real * 8 pi,b0,b1,b2,t,xlt,ot,as0,as1
      integer icount
      parameter (pi=3.14159265358979312D0)
      b0  = (33.d0-2.d0*nf)/(12.d0*pi)
      b1  = (153.d0 - 19.d0*nf) / (24.d0*pi**2)
      b2  = (2857.d0/2.d0-5033.d0/18.d0*nf+325.d0/54.d0*nf**2)
     #     /(4.d0*pi)**3
      t  = 1/b0/as
      icount=0
    1 xlt = log(t)
      if(icount.gt.10000) then
          write(*,*) ' xlambd: cannot converge '
          stop
      endif
      icount=icount+1
      ot = t
c-----------------------------------------------------------
c Value and Derivative of alfa with respect to t
      as0  =   1/(t*b0)*(1-b1/b0**2*log(t)/t
     #         +(b1/b0**2*log(t)/t)**2
     #       -(b1**2*(log(t)+1)-b0*b2)/b0**4/t**2)
      as1  =
     5(-2*b1**2*log(t)**2/(b0**4*t**3)+2*(b1**2*(log(t)+1)-b0*b2)/(b0**4
     1   *t**3)+b1*log(t)/(b0**2*t**2)+2*b1**2*log(t)/(b0**4*t**3)-b1/(b
     2   0**2*t**2)-b1**2/(b0**4*t**3))/(b0*t)-(b1**2*log(t)**2/(b0**4*t
     3   **2)-(b1**2*(log(t)+1)-b0*b2)/(b0**4*t**2)-b1*log(t)/(b0**2*t)+
     4   1)/(b0*t**2)
      t  = (as-as0)/as1 + t
      if(abs(ot-t)/ot.gt.0.00000001d0)goto 1
      genericxlambdNNL = q/exp(t/2)
      return
      end

      subroutine genericpdf(ndns,ih,xmu2,x,fx)
c Interface to lhapdf package.
      implicit none
      include 'nlegborn.h'
      include 'pwhg_pdf.h'

      real * 8 lam5
      integer iord,iset
      common/cgenericpdf/lam5,iord,iset


      integer ndns,ih
      real * 8 xmu2,x,fx(-pdf_nparton:pdf_nparton)
      real * 8 fxlha(-6:6)
      integer j
      real * 8 tmp
      real*8 photon
      call genericpdfset(ndns)

      call xfxq2(iset,x,xmu2,fxlha)
c photon induced work only with MRST2004QED (ndns = 20460)
      if (ndns.eq.20460) then
          call xfphoton(x,sqrt(xmu2),photon)
      else
          photon=0d0
      endif

      fx=0
      fx(-6:6)=fxlha/x
      if(pdf_nparton.ge.22) then
         fx(22)=photon/x
      endif
c 1 is proton, -1 is antiproton, 3 is pi+, -3 is pi-
      if(ih.eq.1) then
         return
      elseif(ih.eq.-1) then
         do j=1,6
            tmp=fx(j)
            fx(j)=fx(-j)
            fx(-j)=tmp
         enddo
      elseif(ih.eq.3) then
         tmp=fx(1)
         fx(1)=fx(-1)
         fx(-1)=tmp
      elseif(ih.eq.-3) then
         tmp=fx(2)
         fx(2)=fx(-2)
         fx(-2)=tmp
      elseif(ih.eq.0) then
c 0 is deuteron
        fx(1)  = 0.5 * ( fx(1)+fx(2) )
        fx(-1) = 0.5 * ( fx(-1)+fx(-2) )
        fx(2)  = fx(1)
        fx(-2) = fx(-1)         
      elseif(ih.eq.4) then
c photon pdf
         continue
      else
         write(*,*) ' genericpdf: unimplemented hadron type ',ih
         stop
      endif
c Bug fixes for version 5.3 of lhapdf
c      if(ndns.eq.363) then
c         do j=1,6
c            fx(j)=fx(j)/2
c            fx(-j)=fx(-j)/2
c         enddo
c      endif
      
      end

      subroutine genericpdfpar(ndns,ih,xlam,scheme,iorder,iret)
      implicit none

      integer iord,iset
      real * 8 lam5
      common/cgenericpdf/lam5,iord,iset

      integer ndns,ih,iorder
      real * 8 xlam
      character * 2 scheme
      integer iret

      call genericpdfset(ndns)
      scheme='MS'
      iret=0
      xlam=lam5
      iorder=iord
      end

      function whichpdfpk()
      character * 3 whichpdfpk
      whichpdfpk='lha'
      end
