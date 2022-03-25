      subroutine pdfcall(ih,x,pdf)
      implicit none
      include 'pwhg_pdf.h'
      integer ih
      real * 8 x,pdf(-pdf_nparton:pdf_nparton)
      include 'pwhg_st.h'
      integer npdf
      if(ih == 1 .and. pdf_ndns1 == -100) then
         pdf = 1
         return
      elseif(ih == 2 .and. pdf_ndns2 == -100) then
         pdf = 1
         return
      endif
      if(x .ge. 1) then
         if(x-1 .gt. 1d-4) then
            write(*,*) 'pdfcall: warning, x=',x
            write(*,*) 'returning pdf=0'
         endif
         pdf = 0
         return
      endif
      if(ih.eq.1) then
         call genericpdf0(pdf_ndns1,pdf_ih1,st_mufact2,x,pdf)
      elseif(ih.eq.2) then
         call genericpdf0(pdf_ndns2,pdf_ih2,st_mufact2,x,pdf)
      else
         write(*,*) ' pdfcall: invalid call, ih=',ih
         stop
      endif
      end


c Front end to genericpdf; it stores the arguments and return values of
c the nrec most recent calls to genericpdf. When invoked it looks in the
c stored calls; if a match its found, its return value is used.
c In this framework it is found that nrec=8 would be enough.
c This provides a remarkable increase in spead (better than a factor of 3)
c when cteq6 pdf are used.
      subroutine genericpdf0(ns,ih,xmu2,x,fx)
      implicit none
      include 'pwhg_pdf.h'
      integer maxparton
      parameter (maxparton=22)
      integer ns,ih
      real * 8 xmu2,x,fx(-pdf_nparton:pdf_nparton)
      integer nrec
      parameter (nrec=10)
      real * 8 oxmu2(nrec),ox(nrec),ofx(-maxparton:maxparton,nrec)
      integer ons(nrec),oih(nrec)
      integer irec
      save oxmu2,ox,ofx,ons,oih,irec
c set to impossible values to begin with
      data ox/nrec*-1d0/
      data irec/0/
      integer j,k
      real * 8 charmthr2,bottomthr2
      logical ini
      data ini/.true./
      save ini,charmthr2,bottomthr2
      real * 8 powheginput
      external powheginput
      real * 8 tmp
c$$$      if(ini) then
c$$$         charmthr2=powheginput('#charmthrpdf')
c$$$         bottomthr2=powheginput('#bottomthrpdf')
c$$$         if(charmthr2.lt.0) charmthr2=1.5
c$$$         if(bottomthr2.lt.0) bottomthr2=5
c$$$         charmthr2=charmthr2**2
c$$$         bottomthr2=bottomthr2**2
c$$$         ini=.false.
c$$$      endif
      do j=irec,1,-1
         if(x.eq.ox(j)) then
            if(xmu2.eq.oxmu2(j)) then
               if(ns.eq.ons(j).and.ih.eq.oih(j)) then
                  fx=ofx(-pdf_nparton:pdf_nparton,j)
                  return
               endif
            endif
         endif
      enddo
      do j=nrec,irec+1,-1
         if(x.eq.ox(j)) then
            if(xmu2.eq.oxmu2(j)) then
               if(ns.eq.ons(j).and.ih.eq.oih(j)) then
                  fx=ofx(-pdf_nparton:pdf_nparton,j)
                  return
               endif
            endif
         endif
      enddo
      irec=irec+1
      if(irec.gt.nrec) irec=1
      ons(irec)=ns
      oih(irec)=ih
      oxmu2(irec)=xmu2
      ox(irec)=x
      if(abs(ih) == 2) then
c     Do the neutron or antineutron
         call genericpdf(ns,ih/abs(ih),xmu2,x,
     1        ofx(-pdf_nparton:pdf_nparton,irec))
         tmp = ofx(1,irec)
         ofx(1,irec)=ofx(2,irec)
         ofx(2,irec)=tmp
         tmp = ofx(-1,irec)
         ofx(-1,irec)=ofx(-2,irec)
         ofx(-2,irec)=tmp
      else
         call genericpdf(ns,ih,xmu2,x,
     1        ofx(-pdf_nparton:pdf_nparton,irec))
      endif
c Flavour thresholds:
c$$$      if(xmu2.lt.bottomthr2) then
c$$$         ofx(5,irec)=0
c$$$         ofx(-5,irec)=0
c$$$      endif
c$$$      if(xmu2.lt.charmthr2) then
c$$$         ofx(4,irec)=0
c$$$         ofx(-4,irec)=0
c$$$      endif
      if(powheginput('#negative_pdfs_zero').eq.1) then
         do k=-pdf_nparton,pdf_nparton
            if (ofx(k,irec).lt.0) then
               call increasecnt("negative pdf values");
               ofx(k,irec)=0
            endif
         enddo
      endif
      fx=ofx(-pdf_nparton:pdf_nparton,irec)
      end

