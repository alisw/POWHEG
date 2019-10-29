      subroutine inc_nhist
c increase nhist
      include 'pwhg_bookhist-multi-new.h'
      type(histptr), pointer :: save_ptr(:)

      if(nhist == 0) then
         allocate(hist_ptr(20))
         nhist = 20
      else
c allocate larger array
         allocate(save_ptr(2*nhist))
c transfer the entries already present in hist_ptr
         save_ptr(1:nhist) = hist_ptr
c associate hist_ptr to the new array
         hist_ptr => save_ptr
c IMPORTANT REMARK: we do not deallocate hist_ptr before associating it
c with the save_ptr pointer. By doing it, the memory taken by its previous
c allocation remain inaccessible and unclaimed. However, this avoids subtle
c errors: if this subroutine was called by other subroutines that receive
c the hist_ptr pointer as an argument, the associated dummy argument of the
c calling subroutine becomes undefined as hist_ptr is deallocate.
         nhist = 2 * nhist
      endif
      end

      subroutine inihists
      implicit none
      include 'pwhg_bookhist-multi-new.h'
      jhist = 0
      nhist = 0
c nmulti defaults to 1;
c It should be set at the first call of the analysis.
      nmulti = 1
      end

      subroutine setupmulti(n)
      implicit none
      integer n
      include 'pwhg_bookhist-multi-new.h'
      real * 8 weirdnum
      common/c_setupmulti/weirdnum
      save /c_setupmulti/
c if this is set we are sure that setupmulti was called
      weirdnum=317d0/12345d0
      if (nmulti.ne.n) then
         nmulti=n
         call rebookhist
      endif
      end

      subroutine bookupeqbins(string,binsize,xlow,xhigh)
      implicit none
      character *(*) string
      real * 8 binsize,xlow,xhigh
      real * 8, allocatable :: x(:)
      real * 8 xx,xhighmt,tiny
      parameter (tiny=1d-4)
      integer j,nbins
      xx = xlow
      nbins = 0
      xhighmt = xhigh - binsize * tiny
      do while(xx < xhighmt)
         nbins = nbins + 1
         xx = xx + binsize
      enddo
      allocate(x(nbins+1))
      x(1) = xlow
      do j=1,nbins
         x(j+1) = x(j) + binsize
         if(abs(x(j+1))/binsize < tiny) then
c avoid funny numbers for bin edges near zero
            x(j+1) = 0
         endif
      enddo
      if((x(nbins+1)-xhigh)/binsize.gt.tiny) then
         write(*,*) 'upper limit incompatible with bin size'
         write(*,*) 'replacing ',xhigh,' with ',x(nbins+1)
         write(*,*) ' in histogram ',string
      endif

      call bookup(string,nbins,x)

      deallocate(x)

      end


      subroutine bookup(string,n,x)
c Books up a histogram characterized by the tag string <string>,
c with n bins. The array x(n+1) are the bins endpoints,
c x(i) is the low extreme of bin i.
      implicit none
      character *(*) string
      integer n
      real * 8 x(n+1)
      include 'pwhg_bookhist-multi-new.h'
      integer j,k
      integer indexhist
      real * 8 weirdnum
      common/c_setupmulti/weirdnum
      save /c_setupmulti/
      logical ini
      data ini/.true./
      save ini

c We assume that this routine is always called first when the package is used
      if(ini) then
         if(.not.weirdnum.eq.317d0/12345d0) then
c setupmulti was not called! setup default value!
            call setupmulti(1)
         endif
         ini=.false.
      endif

c indexhist(string) returns the histogram index if a histogram
c with tag string was already booked, otherwise it books a new histogram,
c and returns minus the value of its index
      j=-indexhist(string)
      if(j.lt.0) then
         write(*,*) 'Histogram ',string,' already booked'
         call exit(-1)
      endif

      allocate(hist_ptr(j)%xhistarr(n+1))
      allocate(hist_ptr(j)%yhistarr(nmulti,0:n+1))
      allocate(hist_ptr(j)%yhistarr1(nmulti,0:n+1))
      allocate(hist_ptr(j)%errhistarr1(nmulti,0:n+1))
      allocate(hist_ptr(j)%yhistarr2(nmulti,0:n+1))
      allocate(hist_ptr(j)%errhistarr2(nmulti,0:n+1))
      allocate(hist_ptr(j)%nhits(0:n+1))

      hist_ptr(j)%xhistarr=x

c y and err values go from 0 to n+1, 0 being the underflow and n+1
c the overflow.
      hist_ptr(j)%yhistarr=0
      hist_ptr(j)%yhistarr1=0
      hist_ptr(j)%errhistarr1=0
      hist_ptr(j)%yhistarr2=0
      hist_ptr(j)%errhistarr2=0
      hist_ptr(j)%nhits=0

      hist_ptr(j)%nbins = n
      hist_ptr(j)%ient1 = 0
      hist_ptr(j)%nmulti = nmulti

      end


      function indexhist(string)
      implicit none
      character * (*) string
      include 'pwhg_bookhist-multi-new.h'
      integer indexhist
      integer j,khist
      if(string.eq.' ') then
         write(*,*) ' indexhist: error, empty name'
         call exit(-1)
      endif
      khist = -1
      do j=1,jhist
         if(hist_ptr(j)%id.eq.' ') then
            khist = j
         endif
         if(hist_ptr(j)%id.eq.string) then
            indexhist=j
            return
         endif
      enddo

      if(khist.lt.0) then
         if(jhist.eq.nhist) then
            call inc_nhist
         endif
         jhist=jhist+1
         khist = jhist
      endif

      hist_ptr(khist)%id=trim(adjustl(string))
      if(hist_ptr(khist)%id.ne.trim(adjustl(string))) then
         write(*,*) ' Histogram string "',string,'" too long'
         call exit(-1)
      endif
c the negative sign indicates a new histogram
      indexhist=-khist

      end


      subroutine deletehist(string,iret)
      implicit none
      character * (*) string
      integer iret
      include 'pwhg_bookhist-multi-new.h'
      integer j
      if(string.eq.' ') then
         write(*,*) ' deletehist: error, empty name'
         call exit(-1)
      endif
      do j=1,jhist
         if(hist_ptr(j)%id.eq.string) then
            hist_ptr(j)%id=' '
            deallocate(hist_ptr(j)%xhistarr)
            deallocate(hist_ptr(j)%yhistarr)
            deallocate(hist_ptr(j)%yhistarr1)
            deallocate(hist_ptr(j)%errhistarr1)
            deallocate(hist_ptr(j)%yhistarr2)
            deallocate(hist_ptr(j)%errhistarr2)
            deallocate(hist_ptr(j)%nhits)
            return
            iret = 0
         endif
      enddo
      write(*,*) ' deletehist: histogram '//string
     1     //' not present'
      iret = -1
      end


      subroutine filld(string,xval,weight)
      implicit none
      character *(*) string
      include 'pwhg_bookhist-multi-new.h'
      real * 8 xval,weight(1:nmulti)
      integer j,k,indexhist
      j=indexhist(string)
      if(j.lt.0) then
         write(*,*) ' histogram "',string,'" was not booked'
         call exit(-1)
      endif
c Make sure nmulti hasn't changed; fix up for older code
      if(nmulti.gt.hist_ptr(j)%nmulti) then
         call rebookhist
      endif
c     underflow
      if(xval.lt.hist_ptr(j)%xhistarr(1)) then
         hist_ptr(j)%yhistarr(1:nmulti,0)=
     1        hist_ptr(j)%yhistarr(1:nmulti,0)+weight
         hist_ptr(j)%nhits(0)=hist_ptr(j)%nhits(0)+1
         return
      else
         do k=1,hist_ptr(j)%nbins
            if(xval.lt.hist_ptr(j)%xhistarr(k+1)) then
               hist_ptr(j)%yhistarr(1:nmulti,k)=
     1              hist_ptr(j)%yhistarr(1:nmulti,k)+weight/
     2              (hist_ptr(j)%xhistarr(k+1)-hist_ptr(j)%xhistarr(k))
               hist_ptr(j)%nhits(k)=hist_ptr(j)%nhits(k)+1
               return
            endif
         enddo
      endif
c overflow
      hist_ptr(j)%yhistarr(1:nmulti,hist_ptr(j)%nbins+1) =
     1     hist_ptr(j)%yhistarr(1:nmulti,hist_ptr(j)%nbins+1) + weight
      end

      subroutine rebookhist
      implicit none
      integer j
      include 'pwhg_bookhist-multi-new.h'
      integer n
      logical alreadycalledonce
      data alreadycalledonce/.false./
      save alreadycalledonce
c rebook all nmulti dependent arrays.
c It is better to book the histogram after the number of weights
c are made available. This happens as soon as the first
c event is read. Older code was not organized in this way.
c In order not to break it, if nmulti changes by the time
c one tries to fill the first histogram, we rebook all of them
c according to the current nmulti. This is just a fix, and should
c not happen more than once.

      if(alreadycalledonce) then
         write(*,*) ' filld: error: number of weights'//
     1        ' no longer consistent'
         write(*,*) '        with its previous vlaue'
         write(*,*) ' exiting ...'
         call exit(-1)
      endif

      do j=1,jhist
         if(hist_ptr(j)%ient1.eq.0) then

            deallocate(hist_ptr(j)%yhistarr)
            deallocate(hist_ptr(j)%yhistarr1)
            deallocate(hist_ptr(j)%errhistarr1)
            deallocate(hist_ptr(j)%yhistarr2)
            deallocate(hist_ptr(j)%errhistarr2)

            n = hist_ptr(j)%nbins

            allocate(hist_ptr(j)%yhistarr(nmulti,0:n+1))
            allocate(hist_ptr(j)%yhistarr1(nmulti,0:n+1))
            allocate(hist_ptr(j)%errhistarr1(nmulti,0:n+1))
            allocate(hist_ptr(j)%yhistarr2(nmulti,0:n+1))
            allocate(hist_ptr(j)%errhistarr2(nmulti,0:n+1))

            hist_ptr(j)%yhistarr = 0
            hist_ptr(j)%yhistarr1 = 0
            hist_ptr(j)%errhistarr1 = 0
            hist_ptr(j)%yhistarr2 = 0
            hist_ptr(j)%errhistarr2 = 0

            hist_ptr(j)%nmulti = nmulti
         else
            write(*,*) ' filld: error: number of weights'//
     1           ' no longer consistent'
            write(*,*) '        with its previous value'
            write(*,*) ' exiting ...'
            call exit(-1)
         endif
      enddo

      end



      subroutine pwhgtopout(filename)
      implicit none
      character * (*) filename
      include 'pwhg_bookhist-multi-new.h'
      integer k,j,iun,l
      character * 3 cl
      call newunit(iun)
      do l=1,nmulti
         if(nmulti.eq.1) then
            open(unit=iun,file=trim(adjustl(filename))//'.top',
     1           status='unknown')
         else
            write(cl,'(i3)') l
            open(unit=iun,file=trim(adjustl(filename))//'-W'//
     1           trim(adjustl(cl))//'.top',status='unknown')
         endif
         do j=1,jhist
            if(hist_ptr(j)%id .ne. ' ') then
               write(iun,'(a,i3)')'# '//trim(adjustl(hist_ptr(j)%id))//
     1           ' index ',j-1
               do k=1,hist_ptr(j)%nbins
                  write(iun,'(4(1x,e14.8))') hist_ptr(j)%xhistarr(k),
     1                 hist_ptr(j)%xhistarr(k+1),
     2                 hist_ptr(j)%yhistarr2(l,k),
     3                 hist_ptr(j)%errhistarr2(l,k)
               enddo
               write(iun,*)
               write(iun,*)
            endif
         enddo
         close(iun)
      enddo
      end


      subroutine pwhgaccumup
c values histogrammed so far are transferred to array yhistarr1,
c and the square of the values are transferred to array errhistarr1.
c yhistarr is zeroed. The index ient1 is increased by one unit.
      implicit none
      include 'pwhg_bookhist-multi-new.h'
      integer j,k
      do j=1,jhist
         do k=0,hist_ptr(j)%nbins + 1
            hist_ptr(j)%yhistarr1(1:nmulti,k)=
     1           hist_ptr(j)%yhistarr1(1:nmulti,k)
     2           +hist_ptr(j)%yhistarr(1:nmulti,k)
            hist_ptr(j)%errhistarr1(1:nmulti,k)=
     1           hist_ptr(j)%errhistarr1(1:nmulti,k)
     2           +hist_ptr(j)%yhistarr(1:nmulti,k)**2
            hist_ptr(j)%yhistarr(1:nmulti,k)=0
         enddo
         hist_ptr(j)%ient1 = hist_ptr(j)%ient1 + 1
      enddo
      end

      subroutine pwhgsetout
c provides a snapshot of the current result of the
c analysis, leaving the yhistarr1 and errhistarr1 unchanged.
      implicit none
      include 'pwhg_bookhist-multi-new.h'
      integer j,k
      real *8 xxx,sum(1:nmulti),sumsq(1:nmulti)
      do j=1,jhist
         xxx = 1d0/hist_ptr(j)%ient1
         do k = 0, hist_ptr(j)%nbins + 1
            sum = hist_ptr(j)%yhistarr1(1:nmulti,k)
            sumsq = hist_ptr(j)%errhistarr1(1:nmulti,k)
            hist_ptr(j)%yhistarr2(1:nmulti,k) = xxx*sum
            hist_ptr(j)%errhistarr2(1:nmulti,k) =
     1           sqrt(xxx**2*abs(sumsq-sum**2*xxx))
         enddo
      enddo
      end

      subroutine pwhgaddout
c accumulates the results obtained so far in yhistarr2 and errhistarr2.
c It zeroes yhistarr1 and errhistarr1. To be used if we compute
c a cross section with several contributions.
      implicit none
      include 'pwhg_bookhist-multi-new.h'
      integer j,k
      real *8 xxx,sum(1:nmulti),sumsq(1:nmulti)
      do j=1,jhist
         xxx=1d0/hist_ptr(j)%ient1
         do k=0,hist_ptr(j)%nbins+1
            sum=hist_ptr(j)%yhistarr1(1:nmulti,k)
            sumsq=hist_ptr(j)%errhistarr1(1:nmulti,k)
            hist_ptr(j)%yhistarr2(1:nmulti,k)=
     1           hist_ptr(j)%yhistarr2(1:nmulti,k)+xxx*sum
            hist_ptr(j)%errhistarr2(1:nmulti,k)=
     1           sqrt(hist_ptr(j)%errhistarr2(1:nmulti,k)**2+
     2           xxx**2*abs(sumsq-sum**2*xxx))
         enddo
      enddo
      do j=1,jhist
         hist_ptr(j)%yhistarr1 = 0
         hist_ptr(j)%errhistarr1 = 0
         hist_ptr(j)%ient1 = 0
      enddo
      end
