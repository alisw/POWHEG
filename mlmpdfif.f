      subroutine genericpdf(ndns,ih,xmu2,x,fx)
c Interface to mlmpdf package.
      implicit none
      include 'nlegborn.h'
      include 'pwhg_pdf.h'
      integer ndns,ih
      real * 8 xmu2,x,fx(-pdf_nparton:pdf_nparton)
      real * 4 sxmu2,sx,sfx(-5:5)
      integer j
      sx=x
      sxmu2=xmu2
      fx=0
      call mlmpdf(ndns,ih,sxmu2,sx,sfx,5)
      do j=-5,5
         fx(j)=sfx(j)
      enddo
      fx(1)=sfx(2)
      fx(-1)=sfx(-2)
      fx(2)=sfx(1)
      fx(-2)=sfx(-1)
      end

c This subroutine is in LHAPDF, and is invoked in
c setstrongcoupl.f in case the lhapdfif.f (and LHAPDF) is linked in
c If not, like now, it should be set to a dummy function with the same (dummy) arguments as in LHAPDF, to avoid
c link errors. It is never invoked in the present case.
      subroutine getq2min(dum1,dum2)
      implicit none
      integer dum1
      real * 8 dum2
      end

      subroutine genericpdfpar(ndns,ih,xlam,scheme,iorder,iret)
      implicit none
      include 'pwhg_pdf.h'
      integer ndns,ih
      real * 8 xlam
      character * 2 scheme
      integer iret,iorder
      call pdfpar(ndns,ih,xlam,scheme,iret)
c set ad hoc value for q2min (not provided in mlmpdf)
      pdf_q2min = 2d0
c not yet implemented
      iorder=-1
      end

      function whichpdfpk()
      character * 3 whichpdfpk
      whichpdfpk='mlm'
      end

c     Dummy function to avoid compilation errors.
c     This function exists in LHAPDF v5, and is defined in lhapdf6if.f
c     when using LHAPDF v6. When the code is compiled with native pdf,
c     a dummy function is needed, as a call to alphasfrompdf is present in
c     setstrongcoupl.f (inside an if statement).
c     With native pdf, this function should never be used, hence we put an error
c     message.
      function alphasfrompdf(mu)
      implicit none
      real *8 alphasfrompdf,mu
      write(*,*)
      write(*,*) '***************************************************'
      write(*,*) ' error: alphasfrompdf called, but not provided'
      write(*,*) ' by native pdf'
      write(*,*) ' either switch to lhapdf, or remove the flag'
      write(*,*) ' alphas_from_pdf from the input card'
      write(*,*) '***************************************************'
      call pwhg_exit(-1)
      end
