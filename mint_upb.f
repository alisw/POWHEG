c      implicit none
c      integer ndim
c      character * 20 pwgprefix
c      integer lprefix
c      common/cpwgprefix/pwgprefix,lprefix
c      include 'pwhg_rnd.h'
c      parameter (ndim=19)
c      real * 8 ymax(50,ndim),xint,ymmm,ymin
c      integer k,j
c      lprefix=3
c      pwgprefix='pwg'
c      rnd_cwhichseed='none'
c      xint=8.8d-3
c      call loadmintupb(ndim,'btildeupb',xint,ymax)
c      ymmm=0
c      ymin=1d30
c      do k=1,50
c         do j=1,ndim
c            ymmm=max(ymmm,ymax(k,j))
c            ymin=min(ymin,ymax(k,j))
c         enddo
c         write(*,'(19(d8.2,1x))') (ymax(k,j),j=1,ndim)
c      enddo
c      write(*,*) ymmm, ymin, xint**(1d0/ndim)
c      end


c initialize the storage of values for the determination of the
c upper bounding envelope in MINT
      subroutine startstoremintupb(filetag)
      implicit none
      include 'pwhg_rnd.h'
      include 'pwhg_flg.h'
      character * (*) filetag
      character * 20 pwgprefix
      integer lprefix
      common/cpwgprefix/pwgprefix,lprefix
      integer iunit
      logical active
      common/storeubc/iunit,active
      save /storeubc/
      data active/.false./
      integer iret
      active=.true.
      if(rnd_cwhichseed.eq.'none') then
         call pwhg_io_open_write(pwgprefix(1:lprefix)//
     1    trim(filetag)//'.dat',iunit,flg_compress_upb,iret)
      else
         call pwhg_io_open_write(pwgprefix(1:lprefix)//
     1        trim(filetag)//'-'//rnd_cwhichseed//'.dat',
     2        iunit,flg_compress_upb,iret)
      endif
      if(iret /= 0) then
         write(*,*) 'startstoremintupb: cannot open output file'
         write(*,*) 'exiting ...'
         call exit(-1)
      endif
      end

      subroutine storemintupb(ndim,ncell,imode,f,f0)
      implicit none
      include 'nlegborn.h'
      integer ndim,ncell(ndim),imode
      real * 8 f,f0
      integer iunit
      logical active
      common/storeubc/iunit,active
      character * 30 fmt
      character(len=2*ndim+22) string
      integer k
      logical ini
      data ini/.true./
      save ini,fmt
      if(active) then
         if(ini) then
            fmt='(   i2,2d11.5)'
            write(fmt(2:4),'(i2)') ndim
            ini=.false.
         endif
         if(imode.eq.0) then
            write(string,fmt) (ncell(k),k=1,ndim),f,f0
            call pwhg_io_write(iunit,string)
         else
            write(string,fmt) (ncell(k),k=1,ndim),f
            call pwhg_io_write(iunit,string(1:2*ndim+11))
         endif
      endif
      end

      subroutine stopstoremintupb
      implicit none
      integer iunit
      logical active
      common/storeubc/iunit,active
      call pwhg_io_close(iunit)
      active=.false.
      end

      subroutine getlinemintupb1(filetag,ndim,cells,f,f0,iret)
c iret = 0 normally, iret = 1 on end of data.
c Subsequent call restart from the beginning of the data set.
c A special call with iret = -10 causes the deallocation
c of the memory arrays used to store the data.
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flg.h'
      character *(*) filetag
      integer ndim,cells(ndiminteg),iret
      real * 8 f,f0
      character * 1, dimension(:,:), allocatable, save :: allcells
      real, dimension(:), allocatable, save :: allf
      real, dimension(:), allocatable, save :: allf0
      real * 8 upblimitsv1,upblimitsv2
      logical print_outliers
      common/cfind_outliers/upblimitsv1,upblimitsv2,print_outliers
      integer nlines,j
      integer status,index
c the internal flag status is
c 0     if the data has not been loaded into memory
c       (typically upon the first invocation)
c 1     if the data s in memory
c the integer index is the data line to be read.
c It is increased after each call, and upon the last
c call (the one returning iret = 1) it is reset to 1
      data status/0/
      save status,index,nlines
c a call with iret=-10 deallocate all arrays and returns ;
      if(iret.eq.-10) then
         write(*,*) ' deallocating mintupb arrays'
         deallocate(allcells,stat=j)
         deallocate(allf,stat=j)
         deallocate(allf0,stat=j)
         write(*,*) ' end deallocating '
         status=0
         return
      endif
      if(status.eq.0) then
         write(*,*) ' getlinemintupb1: loading file(s)'
c status=0 is the initial call;
         iret=0
         nlines=0
c     In this block, find maximum below outliers
         upblimitsv1=1d308  ! DBL_MAX
         upblimitsv2=1d308
         print_outliers = .false.
         if(filetag.eq.'btildeupb'.and.flg_fastbtlbound) then
            call find_outliers_limit('startdouble',f,f0)
         else
            call find_outliers_limit('startsimple',f,f0)
         endif
         do while(iret == 0)
            call getlinemintupb(filetag,ndim,cells,f,f0,iret)
            if(iret == 0) call find_outliers_limit('add',f,f0)
         enddo
         call find_outliers_limit('lim',upblimitsv1,upblimitsv2)
         call store_outliers_limit('put',trim(filetag),
     1        upblimitsv1,upblimitsv2)
         iret = 0
c     in this block count the lines in the file (nlines
         write(*,*) ' getlinemintupb1: counting lines in files'
c     this flag determines whether outliers found in the sequence are
c     printed in the log file and counted in the counters.
         print_outliers = .true.
 1       continue
         call getlinemintupb(filetag,ndim,cells,f,f0,iret)
         if(iret.eq.0) then
            nlines=nlines+1
            goto 1
         endif
c     in the following repeated reading of the stored calls
c     outliers are no longer counted and printed, since they
c     have already been counted.
         print_outliers = .false.
c lines counted
         write(*,*) ' getlinemintupb1: found ',nlines,' lines in files'
c allocates enough stuff for nlines lines
         deallocate(allcells,stat=iret)
         allocate(allcells(ndim,nlines),stat = iret)
         deallocate(allf,stat=iret)
         allocate(allf(nlines),stat=iret)
         deallocate(allf0,stat=iret)
         allocate(allf0(nlines),stat=iret)
c store file content in allocated array
         write(*,*) ' getlinemintupb1: load content from files'
         do j=1,nlines
            call getlinemintupb(filetag,ndim,cells,f,f0,iret)
            call inttochar(cells,allcells(:,j),ndim)
            allf(j)=f
            allf0(j)=f0
         enddo
         write(*,*) ' getlinemintupb1: content loaded'
c this call should return 1 in iret (no more lines)
         call getlinemintupb(filetag,ndim,cells,f,f0,iret)
c Set status to 1 (cells loaded), index to 1 (initiate reading)
         status=1
         index=1
         write(*,*) 'getlinemintupb1: loaded file(s)'
      endif
      if(index.le.nlines) then
         f=allf(index)
         f0=allf0(index)
         call chartoint(allcells(:,index),cells,ndim)
         index=index+1
         iret=0
      else
c end of data
         iret=1
         index=1
      endif
      end
         
      subroutine getlinemintupb(filetag,ndim,cells,f,f0,iret)
c reads a line from the upper bound file or files.
c iret=-1: failure
c iret=0 : success
c iret=1 : end of stream; next call restart reading
c          from the beginning
      implicit none
      character *(*) filetag
      integer ndim,cells(ndim),iret
      real * 8 f,f0
      include 'pwhg_rnd.h'
      include 'pwhg_flg.h'
      character * 20 pwgprefix
      character(len=2*ndim+22) string
      integer lprefix
      common/cpwgprefix/pwgprefix,lprefix
      integer ltag
      character * 100 fname
      integer jfile,k,iunit
      character * 20 fmt
      character * 4 chnum
      logical lpresent
      real * 8 upblimitsv1,upblimitsv2
      logical print_outliers
      common/cfind_outliers/upblimitsv1,upblimitsv2,print_outliers
      integer status
      data status/0/
      save status,iunit,jfile,fmt
c initial call
      if(status.eq.0) then
         fmt='(   i2,2d11.5)'
         ltag=len(filetag)
 33      if(filetag(ltag:ltag).eq.' ') then
            ltag=ltag-1
            goto 33
         endif
c The format for reading the bounds file, for ndim dimensions
         write(fmt(2:4),'(i2)') ndim
c see if there are files to load
         if(rnd_cwhichseed.eq.'none') then
            fname=pwgprefix(1:lprefix)//filetag(1:ltag)//'.dat'
         else
            do jfile=1,9999
               write(chnum,'(i4)') jfile
               if(chnum(1:1).eq.' ') chnum(1:1)='0'
               if(chnum(2:2).eq.' ') chnum(2:2)='0'
               if(chnum(3:3).eq.' ') chnum(3:3)='0'
               fname=pwgprefix(1:lprefix)//filetag(1:ltag)
     1              //'-'//chnum//'.dat'
               inquire(file=fname,exist=lpresent)
               if(lpresent) goto 10
            enddo
            goto 999
         endif
 10      continue
         call pwhg_io_open_read(trim(fname),iunit,iret)
         if(iret /= 0) goto 999
         write(*,*) ' opened ',trim(fname)
c file opened for reading
         status=1
      endif
 12   continue
      call pwhg_io_read(iunit,string,iret)
      if(iret < 0) goto 11
      if(filetag.eq.'btildeupb'.and.flg_fastbtlbound) then
         read(string,fmt=fmt) (cells(k),k=1,ndim),f,f0
         if(flg_storemintupb_nooutliers) then
            if(f > upblimitsv1) then
               if(print_outliers) then
                  call increasecnt(trim(filetag)//' f  outliers')
                  write(*,*) 'getlinemintupb: '//trim(filetag)//
     1                 ' f  outlier=',f
               endif
               goto 12
            elseif(f0 > upblimitsv2) then
               call increasecnt(trim(filetag)//' f0 outliers')
               write(*,*) 'getlinemintupb: '//trim(filetag)//
     1          ' f0  outlier=',f
               goto 12
            endif
         endif
      else
         read(string,fmt=fmt) (cells(k),k=1,ndim),f
         if(flg_storemintupb_nooutliers) then
            if(f > upblimitsv1) then
               if(print_outliers) then
                  call increasecnt(trim(filetag)//' f  outliers')
                  write(*,*) 'getlinemintupb: '//trim(filetag)//
     1                 ' f0  outlier=',f
               endif
               goto 12
            endif
         endif
         f0=-1
      endif
      iret=0
      return
 11   continue
      call pwhg_io_close(iunit)
      if(rnd_cwhichseed.ne.'none') then
 13      jfile=jfile+1
         if(jfile.lt.9999) then
            write(chnum,'(i4)') jfile
            if(chnum(1:1).eq.' ') chnum(1:1)='0'
            if(chnum(2:2).eq.' ') chnum(2:2)='0'
            if(chnum(3:3).eq.' ') chnum(3:3)='0'
            fname=pwgprefix(1:lprefix)//filetag//'-'//chnum//'.dat'
            inquire(file=fname,exist=lpresent)
            if(lpresent) then
               call pwhg_io_open_read(trim(fname),iunit,iret)
               if(iret /= 0) goto 999
               write(*,*) ' opened ',trim(fname)
               goto 12
            else
               goto 13
            endif
         endif
      endif
      iret=1
      status=0
      return
 999  iret=-1
      end
         





      subroutine loadmintupb(ndim,filetag,ymax,ymaxrat)
      implicit none
      include 'pwhg_flg.h'
      include 'pwhg_par.h'
      integer ndim
      character *(*) filetag
      real * 8 ymax(50,ndim),xint,xintrat
      real * 8 ymaxrat(50,ndim)
      integer cells(ndim)
      integer kdim,iret,j,iunit
      real * 8 f,f0,prod,prodrat,xless,xmore,xlessrat,xmorerat
      real * 8 fail,tot,ubtot,failrat,totrat,ubtotrat,ratlim
      integer ipoints,ndiscarded
      integer iterations
      logical ratflg
      iterations=1
      if(filetag.eq.'btildeupb'.and.flg_fastbtlbound) then
         ratflg=.true.
      else
         ratflg=.false.
      endif
      ratlim=par_mintupb_ratlim
c First compute total for f and f/f0 (f/b)
      iret=0
      xint=0
      xintrat=0
      ipoints=0
      ndiscarded = 0
      do
         call getlinemintupb1(filetag,ndim,cells,f,f0,iret)
c        iret = 1 is end of file.
         if(iret == 1) exit
         if(iret /= 0) goto 998
         if(f0>0) then
            if(f/f0 > ratlim) then
               ndiscarded = ndiscarded + 1
               cycle
            endif
         endif
         ipoints=ipoints+1
         xint=xint+f
         if(ratflg) then
            if(f0.gt.0) xintrat=xintrat+f/f0
         endif
      enddo
c     This was introduced in rev. 3205 to reproduce a minor old bug, just not to break binary compatibility
c     of previous versions. In the previous version the exit condition was faulty, the do loop was
c     was executed one more time after getlinemintupb1 returned an end of file. Thus the last
c     value was adde in twice. Here we reproduce the same behaviour if flg_mintupb_xless is not set.
c     The flg_mintupb_xless was introduced in 3205, so if flg_mintupb_xless is true there no binary
c     compatibility with previous output anyhow.
      if(.not.flg_mintupb_xless) then
         xint=xint+f
         if(ratflg) then
            if(f0.gt.0) xintrat=xintrat+f/f0
         endif
      endif
      write(*,*) 'loadmintupb: num. discarded/total:',
     1     ndiscarded,'/',ipoints
      xint=xint/ipoints
      write(*,*) 'loadmintupb: estimated integral for '
     1    //trim(filetag)//' (pos.+|neg.|)  is ',xint
      xintrat=xintrat/ipoints/10
      xmore=0.01
      xmorerat=xmore/5
      if(flg_mintupb_xless) then
         xless=xmore
         xlessrat=xmorerat
      else
         xless=0
         xlessrat=0
      endif
c The bound is initially set as if the function was uniform
c with the given integral
      do kdim=1,ndim
         do j=1,50
            ymax(j,kdim)=xint**(1d0/ndim)
            if(ratflg) then
               if(xintrat.gt.0) ymaxrat(j,kdim)=xintrat**(1d0/ndim)
            endif
         enddo
      enddo
c     Large loop for finding ymax,ymaxrat, etc. Exit when efficiency is OK
      do
c     loop for reading u-bound data
         do
            call getlinemintupb1(filetag,ndim,cells,f,f0,iret)
            if(iret == 1) exit
            if(iret /= 0) then
               write(*,*) ' error while loading bound files'
               call exit(-1)
            endif
            if(f0>0) then
               if(f/f0 > ratlim) then
                  cycle
               endif
            endif

            call evalprod

            if(f.gt.prod) then
               do kdim=1,ndim
                  ymax(cells(kdim),kdim)=
     1                 ymax(cells(kdim),kdim)*(f/prod+0.1)**(xmore/ndim)
               enddo
            else
               do kdim=1,ndim
                  ymax(cells(kdim),kdim)=
     1                 ymax(cells(kdim),kdim)/(2d0)**(xless/ndim)
               enddo
            endif
            if(ratflg) then
               if(f0.gt.0) then
                  if(f/f0.gt.prodrat) then
                     do kdim=1,ndim
                        ymaxrat(cells(kdim),kdim)=
     1                       ymaxrat(cells(kdim),kdim)*
     2                       (f/f0/prodrat+0.1)**(xmorerat/ndim)
                     enddo
                  else
                     do kdim=1,ndim
                        ymaxrat(cells(kdim),kdim)=
     1                       ymaxrat(cells(kdim),kdim)/
     2                       (2d0)**(xlessrat/ndim)
                     enddo
                  endif
               endif
            endif
         enddo
c     check if the failure rate is satisfactory
         fail=0
         tot=0
         ubtot=0
         if(ratflg) then
            failrat=0
            totrat=0
            ubtotrat=0
         endif
         do
            call getlinemintupb1(filetag,ndim,cells,f,f0,iret)
            if(iret.eq.1) exit
            if(iret.ne.0) goto 998
            if(f0>0) then
               if(f/f0 > ratlim) then
                  cycle
               endif
            endif

            call evalprod

            if(f.gt.prod) fail=fail+(f-prod)
            tot=tot+f
            ubtot=ubtot+prod
            if(ratflg) then
               if(f0.ne.0) then
                  if(f.gt.f0*prodrat) failrat=failrat+(f-f0*prodrat)
                  totrat=totrat+f
                  ubtotrat=ubtotrat+f0*prodrat
               endif
            endif
         enddo
         if(fail/tot.gt.1d-3.or.(ratflg.and.failrat/totrat.gt.1d-3))then
c     stop updating the rat grid, if satisfactory
            if(ratflg.and.failrat/totrat.lt.1d-3) then
               ratflg=.false.
               write(*,*) ' ratios envelope efficiency',totrat/ubtotrat
               write(*,*) ' ratios failure estimate',failrat/totrat
            endif
            if(iterations.lt.4) then
               write(*,*) ' iterating upper bounding envelope formation'
            elseif(iterations.lt.5) then
               write(*,*) ' more iterations needed'
               write(*,*) ' this can take a moment ...'
            endif
            write(*,*) ' failure estimate, efficiency ',
     1           fail/tot,(tot-fail)/ubtot
            if(ratflg) then
               write(*,*) ' ratios failure estimate, efficiency ',
     1              failrat/totrat,(totrat-failrat)/ubtotrat
            endif
            iterations=iterations+1
            xless = xless * 0.9
            xlessrat = xlessrat * 0.9
         else
            exit
         endif
      enddo
      if(filetag.eq.'btildeupb'.and.flg_fastbtlbound) then
         ratflg=.true.
      endif
      write(*,*) ' summary of bounds for '//trim(filetag)
      if(ratflg) then
         write(*,*)
     1 ' estimated efficiency without fastbtlbound ',
     2        tot/ubtot
         write(*,*) ' failure estimate',fail/tot
         write(*,*)
     1 ' estimated efficiency with fastbtlbound ',
     2        totrat/ubtotrat
         write(*,*) ' failure estimate',failrat/totrat
      else
         write(*,*) ' envelope efficiency=',tot/ubtot
         write(*,*) ' failure estimate',fail/tot
      endif
      write(*,*) 'processed ',ipoints,' points',iterations,'iterations'
      call newunit(iunit)
      if(filetag == 'btlupb' .and. flg_fastbtlbound) then
         open(unit=iunit,file='mint_upb_'//trim(filetag)//'_rat.top',
     1        status='unknown')
         do kdim=1,ndim
            write(iunit,*) 'set limits x 0 50 y 0 5'
            do j=1,50
               write(iunit,*) j, ymaxrat(j,kdim)
            enddo
            write(iunit,*) 'hist'
            write(iunit,*) 'newplot'
         enddo
         close(iunit)
      endif
      open(unit=iunit,file='mint_upb_'//trim(filetag)//'.top',
     1     status='unknown')
      do kdim=1,ndim
         write(iunit,*) 'set limits x 0 50 y 0 5'
         do j=1,50
            write(iunit,*) j, ymax(j,kdim)
         enddo
         write(iunit,*) 'hist'
         write(iunit,*) 'newplot'
      enddo
      close(iunit)
      call getlinemintupb1(filetag,ndim,cells,f,f0,-10)
      return
 998  continue
      write(*,*) ' error while loading bound files'
      call exit(-1)      

      contains
      subroutine evalprod
      prod=1
      prodrat=1
      do kdim=1,ndim
         prod=prod*ymax(cells(kdim),kdim)
         if(ratflg) prodrat=prodrat*ymaxrat(cells(kdim),kdim)
      enddo
      end subroutine evalprod

      end


      subroutine monitorubound(x,icalls)
      implicit none
      real * 8 x
      integer icalls
      include 'pwhg_rnd.h'
      include 'pwhg_flg.h'
      character * 20 pwgprefix
      integer lprefix
      common/cpwgprefix/pwgprefix,lprefix
      character * 80 file
      integer iun
      if(flg_monitorubound) then
         if(rnd_cwhichseed.eq.'none') then
            file=pwgprefix(1:lprefix)//'boundviolations.dat'
         else
            file=pwgprefix(1:lprefix)//'boundviolations-'//
     1           rnd_cwhichseed//'.dat'
         endif
         call newunit(iun)
         open(unit=iun,file=file,access='append')
         write(iun,*) 'calls=',icalls,'f/ubound',x
         close(iun)
      endif
      end

      subroutine inttochar(int_arr,ch_arr,ndim)
      integer ndim
      integer int_arr(ndim)
      character * 1 ch_arr(ndim)
      integer k
      do k=1,ndim
         ch_arr(k)=char(int_arr(k))
      enddo
      end

      subroutine chartoint(ch_arr,int_arr,ndim)
      integer ndim
      integer int_arr(ndim)
      character * 1 ch_arr(ndim)
      integer k
      do k=1,ndim
         int_arr(k)=ichar(ch_arr(k))
      enddo
      end




      subroutine find_outliers_limit(mode,v1,v2)
c     First call:
c     call find_outliers_limit('startsimple')   ! initializes the subroutine for single value
c     call find_outliers_limit('startdouble')   ! initializes the subroutine for double value
c     call find_outliers_limit('add',v1 (, v2) )  ! add a value (single or double)
c     call find_outliers_limit('lim',v1 (, v2) )  ! returns in v1 (v2) the largest value that is not an outlier.
c     the subroutine builds up a histograms of the number of entries that have been added in bins
c     according to their size. Each bin is a fator of 10 wide (in log scale). When called with
c     mode='lim' it returns the end of the last non-isolated bin.
      character * (*)  mode
      double precision v1,v2
      double precision size, v(2)
      integer logsize
      integer, save :: l
      double precision, save, allocatable :: hist(:,:)
      if(mode == 'startsimple') then
         if(allocated(hist)) then
            write(*,*) 'find_outliers_limit: called with mode='
     1       //trim(mode)//','
            write(*,*) 'but was not reset, exiting ...'
            call exit(-1)
         endif
         allocate(hist(1,-100:100))
         hist = 0
         l=1
      else if(mode == 'startdouble') then
         if(allocated(hist)) then
            write(*,*) 'find_outliers_limit: called with mode='
     1           //trim(mode)//','
            write(*,*) 'but was not reset, exiting ...'
            call exit(-1)
         endif
         allocate(hist(2,-100:100))
         hist = 0
         l=2
      elseif(mode == 'add') then
         do k=1,l
            select case(k)
            case(1)
               v(1) = v1
            case(2)
               v(2) = v2
            end select            
            size=abs(v(k))
            if(size == 0) then
               continue
            else
               logsize = min(max(-100,nint(log(size))),100)
               hist(k,logsize) = hist(k,logsize) + 1
            endif
         enddo
      elseif(mode == 'lim') then
         do k=1,l
            maxj=-100
            do j=-100,100
               if(hist(k,j)>hist(k,maxj)) then
                  maxj=j
               endif
            enddo
            do j=maxj,100
               if(hist(k,j) == 0) then
c Make it two orders larger than the largest value, for safety.                  
                  v(k) = exp(dble(j))*100
                  exit
               endif
            enddo
            select case(k)
            case(1)
               v1 = v(1)
            case(2)
               v2 = v(2)
            end select
         enddo
         deallocate(hist)
      else
         write(*,*) 'find_outliers_limit: called with mode='
     1    //trim(mode)//','
         write(*,*) 'not known, exiting ...'
         call exit(-1)
      endif
      end



      subroutine store_outliers_limit(mode,filetag,v1,v2)
      implicit none
      character *(*) mode, filetag
      real * 8 v1,v2
      character * 20, save :: tags(10)
      real * 8, save :: vals(2,10)
      integer, save :: j=0
      integer k
      if(mode == 'put') then
         if(j>=10) then
            write(*,*) 'store_outliers_limit: too many entries'
            write(*,*) 'exiting ...'
            call exit(-1)
         endif
         do k=1,j
            if(tags(k) == filetag) then
               exit
            endif
         enddo
         if(k == j+1) then
            j=k
            tags(j) = filetag
            vals(1,j) = v1
            vals(2,j) = v2
         endif
      elseif(mode == 'get') then
         v1 = 1d308
         v2 = 1d308
         do k=1,j
            if(tags(k) == filetag) then
               v1 = vals(1,k)
               v2 = vals(2,k)
            endif
         enddo
      else
         write(*,*) 'store_outliers_limit: invalid mode',mode
         write(*,*) 'exiting ...'
         call exit(-1)
      endif
      end
