      subroutine genericpdfset(ndns)
      implicit none
      include 'pwhg_pdf.h'
      integer ndns
c wrap for pdfset; avoids subsequent
c calls to pdfset (you never know)
      integer,allocatable :: seq(:),idset(:),nsets
      real * 8, allocatable :: lam5set(:),q2minset(:)
      save lam5set,q2minset,seq,idset,nsets
      real * 8 genericxlambdL,genericxlambdNL,genericxlambdNNL,
     ,         alphasPDF,asmz,mz
      parameter (mz=91.1876d0)
      real * 8 lam5
      integer iord,iset,maxsets
      common/cgenericpdf/lam5,iord,iset,maxsets
      logical,save:: ini=.true.
      integer j
      real * 8 powheginput
      if(ini) then
         nsets = 0
         iset = 0
         maxsets=nint(powheginput("#lhapdf6maxsets"))
         if(maxsets < 0) maxsets = 10
         allocate(seq(maxsets),idset(maxsets),
     1        lam5set(maxsets),q2minset(maxsets))
         call setlha6init(maxsets)
         ini = .false.
      endif
      if(iset /= 0) then
         if(idset(iset) == ndns) return
      endif
      do j=1,nsets
         if(ndns == idset(j)) then
            iset = j
            call ontop(iset)
            exit
         endif
      enddo
      if(j == nsets+1) then
c...  initalise set
         if(nsets == maxsets) then
c     delete the oldest
            iset = seq(1)
            call setlha6del(iset)
            call ontop(iset)
         else
            nsets = nsets+1
            iset = nsets
            seq(iset)=iset
         endif
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
      contains
      subroutine ontop(j)
c     look for j in seq
      integer j,k,l
      do k=1,nsets
         if(seq(k) == j) then
            do l=k,nsets-1
               seq(l)=seq(l+1)
            enddo
            seq(nsets)=j
            return
         endif
      enddo
      end subroutine
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
      integer iord,iset,maxsets
      common/cgenericpdf/lam5,iord,iset,maxsets


      integer ndns,ih
      real * 8 xmu2,x,fx(-pdf_nparton:pdf_nparton)
      real * 8 fxlha(-6:6)
      integer j,jid
      real * 8 tmp,mu2_loc

      logical generic_has_id
      logical,save:: ini=.true.
      
      if(ini) then
         write(*,*) '==============================='
         write(*,*) 'LHAPDF called in POWHEG' 
         write(*,*) 'pdf cutoff factor = ',pdf_cutoff_fact
         write(*,*) 'pdf cutoff [GeV] = ',
     1    sqrt(pdf_cutoff_fact**2 *pdf_q2min)
         write(*,*) '==============================='
         ini=.false.
      endif


      call genericpdfset(ndns)

      mu2_loc=xmu2
      if(mu2_loc.lt.(pdf_cutoff_fact**2 *pdf_q2min)) then
         mu2_loc=pdf_cutoff_fact**2 *pdf_q2min
      endif

      fx=0

      call xfxq2(iset,x,mu2_loc,fxlha)
      fx(-6:6)=fxlha/x

      do jid=11,pdf_nparton
         if(jid==11 .or. jid==13 .or. jid==15 .or. jid==22) then
            if(generic_has_id(iset,jid)) then
               call xf_pdgid(iset,jid,x,mu2_loc,fx(jid))
               fx(jid)=fx(jid)/x
               if(jid /= 22) then
                  fx(-jid)=fx(jid)
               endif
            endif
         endif
      enddo

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
c ndns is a photon pdf; nothing to do
         continue
      else
         write(*,*) ' genericpdf: unimplemented hadron type ',ih
         call exit(-1)
      endif
      end

      subroutine genericpdfpar(ndns,ih,xlam,scheme,iorder,iret)
      implicit none

      integer iord,iset,maxsets
      real * 8 lam5
      common/cgenericpdf/lam5,iord,iset,maxsets

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

      function alphasfrompdf(mu)
      implicit none
      real * 8 lam5
      integer iord,iset,maxsets
      common/cgenericpdf/lam5,iord,iset,maxsets
      real *8 mu,alphasfrompdf,tmp
      call alphasfrompdf0(iset,mu,tmp)
      alphasfrompdf=tmp
      return
      end
