      subroutine genericpdfset(ndns)
      use rad_tools
      use pdfs_tools
      use internal_parameters
      implicit none
      include 'PhysPars.h'
      include 'pwhg_pdf.h'
      include 'pwhg_st.h'
      include 'minnlo_flg.h'
      integer ndns
      character*100 pdf_name
      character*100 string
      integer stringlength
      real * 8 genericxlambdL,genericxlambdNL,genericxlambdNNL,
     ,         alphasPDF,asmz
      real * 8 lam5
      integer iord,iset,maxsets
      common/cgenericpdf/lam5,iord,iset,maxsets
      logical,save:: ini=.true.

      real * 8 pdf_cutoff_hoppet, powheginput
      real * 8 whichpdfconst, get_M_for_init_Dterms
      common/ccwhichpdfpk/whichpdfconst
      real * 8, parameter :: whichpdfconst0=3.141592653589793d0
      integer imem
      logical pdfs_to_zero, lhapdf_in_hoppet
      real * 8, parameter :: mz_input=91.1876d0
      

c signal that we entered here already
      whichpdfconst=whichpdfconst0
      if(ini) then
c     Here we deal with all inputs that will be needed by hoppet (name
c     of pdf, definition of z mass)
         string = " "
         call lhapdfname(ndns,string,imem)
c     imem is the member name for the given ndns
         pdf_name = trim(string(1:stringlength(string)-1))
     c        //trim(".LHgrid")
         ! set the global mass of the colour singlet (this has to happen before the pdf initialisation)
! at least the next two lines should be irrelevant; and maybe even the third (although it sets the top mass)
         cs%M = get_M_for_init_Dterms()
         cs%Q = get_M_for_init_Dterms()
         call init_masses_Dterms() ! initializes hard-coded mass values in NNLOPS_plugin
         call set_masses(mz_input)

         flg_use_NNLOPS_pdfs = .true.
         if(powheginput('#use_NNLOPS_pdfs').eq.0) flg_use_NNLOPS_pdfs = .false.
         pdfs_to_zero = powheginput('#negative_pdfs_zero').eq.1
         lhapdf_in_hoppet = powheginput('#lhapdf_in_hoppet').eq.1

         pdf_cutoff_fact = powheginput('#pdf_cutoff_fact')
         if(pdf_cutoff_fact.lt.0d0) then
            if(flg_use_NNLOPS_pdfs) then
c     If use_NNLOPS_pdfs is true, we set it to 1.1 by default, so that
c     we have an intermediate scale just above the cutoff of the PDF set.             
               pdf_cutoff_fact = 1.1d0
            else
c     Otherwise, we set it to 1. matching the cutoff of the PDF set.
               pdf_cutoff_fact = 1d0
            endif
         endif
               
c     Now we can initialize hoppet
         if(flg_use_NNLOPS_pdfs) then
c     use the hybrid PDF evolution for the NNLOPS through hoppet
            call init_pdfs_NNLOPS(pdf_name, imem, cutoff_high_fact=pdf_cutoff_fact,
     c           set_pdfs_to_zero_when_negative=pdfs_to_zero)
            call GetQ2min(imem,pdf_q2min)
            write(*,*) 'initialization of PDFs, hoppet interface, '//
     c           'using init_pdfs_NNLOPS: pdf_q2min = ',pdf_q2min
         else
c     use LHAPDF through hoppet
            call init_pdfs_from_LHAPDF(pdf_name, imem, pdf_cutoff_fact,
     c           set_pdfs_to_zero_when_negative=pdfs_to_zero,
     c           lhapdf_in_hoppet=lhapdf_in_hoppet)
            call GetQ2min(imem,pdf_q2min)

            if(powheginput('#profiledscales').ne.0) then
               Q0 = powheginput('#Q0')
               if(Q0.lt.0d0) Q0 = 0d0
               if(Q0.eq.0d0) Q0 = 1d-10
               if(Q0**2.lt.pdf_q2min*pdf_cutoff_fact**2) then
                  write(*,*) 'Error: pdf cutoff > Q0 (profiled scales)'
                  write(*,*) 'This is not consistent.'
                  call exit(-1)
               endif
            endif
            
            write(*,*) 'initialization of PDFs, hoppet interface, '//
     c           'using init_pdfs_from_LHAPDF: pdf_q2min = ',pdf_q2min
         endif
         flg_hoppet_initialized = .true.

         ini = .false.
      endif
      asmz=alphasPDF(mz_input)
      if(iord.eq.0) then
         lam5=genericxlambdNL(asmz,mz_input,5)
      elseif(iord.eq.1) then
         lam5=genericxlambdNL(asmz,mz_input,5)
      elseif(iord.eq.2) then
         lam5=genericxlambdNNL(asmz,mz_input,5)
      endif


      call plotpdf(0.6d0,9d0)
      
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
      use pdfs_tools
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
      integer j
      real * 8 tmp
      real*8 photon
      call genericpdfset(ndns)

      ! call to LHAPDF
      !call xfxq2(iset,x,xmu2,fxlha)
      ! call to HOPPET
      call xfxq2_hoppet(x, xmu2, fxlha)
c photon induced work only with MRST2004QED (ndns = 20460)
c      if (ndns.eq.20460) then
c          call xfphoton(x,sqrt(xmu2),photon)
c      else
c          photon=0d0
c      endif

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
      logical, save :: ini=.true.
c     Tempoorary hack to bypass the checks in init_phys         
      real * 8 whichpdfconst
      common/ccwhichpdfpk/whichpdfconst
      real * 8, parameter :: whichpdfconst0=3.141592653589793d0
      ! MW+ER: removed this here and called it lha to get through stage3
c$$$      if(whichpdfconst==whichpdfconst0) then
c$$$         whichpdfpk='hop'
c$$$      else
c$$$         whichpdfpk='lha'
c$$$      endif
      whichpdfpk='lha'
      end

      function alphasfrompdf(mu)
      use pdfs_tools
      implicit none
      real *8 mu,alphasfrompdf,powheginput
c      running from hoppet (seems to yield ub violations when generating
c      powheg radiation, cfr gen_radiation)
      alphasfrompdf=alphas_hoppet(mu,
     f     lhapdf_in_hoppet=powheginput('#lhapdf_in_hoppet').eq.1)

      end

      subroutine plotpdf(mu2min,mu2max)
      use pdfs_tools
      implicit none
      real * 8 mu2min,mu2max
      real * 8 xmu2,x,fxlha(-6:6)
      real * 8 xmin,xmax
      integer j,k,l,n,iun
      character * 20 file
      logical, save :: ini=.true.
      if(ini) then
         ini=.false.
      else
         return
      endif
      xmin=1d-5
      xmax=0.9999d0
      n=100
      call newunit(iun)
      do l=-5,5
         if(l>=0) then
            write(file,'(i1)') l
            file="plotpdf_"//trim(adjustl(file))//'.top'
         else
            write(file,'(i1)') -l
            file="plotpdf_min_"//trim(adjustl(file))//'.top'
         endif
         open(unit=iun,file=file,status='unknown')
         do k=1,n
            x=exp(log(xmin)+k*log(xmax/xmin)/n)
            write(iun,*) ' title " x=',x,'"'
            do j=1,100
               xmu2=mu2min+j*(mu2max-mu2min)/100
               call xfxq2_hoppet(x, xmu2, fxlha)
               if(fxlha(l) /= 0) then
                  write(iun,*) xmu2,fxlha(l)
               endif
            enddo
            write(iun,*) ' join'
            do j=1,100
               xmu2=mu2min+j*(mu2max-mu2min)/100
               call evolvepdf(x, sqrt(xmu2), fxlha)
               if(fxlha(l) /= 0) then
                  write(iun,*) xmu2,fxlha(l)
               endif
            enddo
            write(iun,*) ' set color red'
            write(iun,*) ' join dashes'
            write(iun,*) ' set color black'
            write(iun,*) ' newplot'
         enddo
         close(iun)
      enddo
      end
      
            
      
