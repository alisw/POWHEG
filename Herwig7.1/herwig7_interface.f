c     herwig7_interface.f - (c) Silvia Ferrario Ravasio, Tomas Jezo, 
c       Paolo Nason and Carlo Oleari

      subroutine herwig7_init(maxev,file)
      implicit none
      integer maxev
      character *(*) file

      character * 6 WHCPRG
      common/cWHCPRG/WHCPRG
      include 'herwigsettings.h'
      real * 8 powheginput
      external powheginput
      logical, save :: ini = .true.
      integer, save :: saved_maxev
      integer innlodec
      common/c_innlodec/innlodec
      integer itmp
      real * 8 ptmin,ptsqmin
      common/ptmin/ptmin


      if(ini) then
         write(*,*)'herwig7_init: trying to open file ',
     $        '<'//trim(file)//'>'
         call opencountlocal(file,maxev)
         write(*,*)'herwig7_init: found ',maxev,' events'
         itmp = powheginput("#maxev")
         if( itmp > 0) then
            maxev = min(maxev, itmp)
         endif
         saved_maxev = maxev
         ptsqmin = powheginput("#ptsqmin")
         if(ptsqmin .le. 0d0) ptsqmin = 0.8d0
         ptmin = sqrt(ptsqmin)
   
      else
c     Herwig calls open many times; we count the events
c     and initialize things only once         
         maxev=saved_maxev
         return
      endif
      ini = .false.

      WHCPRG='HERWIG'
      
      eventCounter = 0


      scalupfac=powheginput('#scalupfac')
      if(scalupfac.lt.0) scalupfac=1

c     read in btilde and remn corrections factors (used together with
c     ubexcess_correct at the generation stage)

      ub_btilde_corr = powheginput('#ub_btilde_corr')
      if (ub_btilde_corr < 0d0) then
        ub_btilde_corr = 1d0
      endif
      ub_remn_corr = powheginput('#ub_remn_corr')
      if (ub_remn_corr < 0d0) then
        ub_remn_corr = 1d0
      endif

      call init_hist

      contains
      subroutine opencountlocal(file,maxev)
      implicit none
      include 'pwhg_rnd.h'
      integer maxev,iun
      character * (*) file
      character * 20 pwgprefix
      integer lprefix
      common/cpwgprefix/pwgprefix,lprefix
      integer ios
      character * 7 string
      real * 8 powheginput
      external powheginput
      integer nev,j
      maxev=0
      call pwhg_io_open_read(trim(file),iun,ios)
c      open(unit=iun,file=file,status='old',iostat=ios)
      if(ios.ne.0) then
         write(*,*) 'cannot open; aborting ...'
         call exit(-1)
      endif
 1    continue
      call pwhg_io_read(iun,string,ios)
      if(ios /= 0) goto 2
      if(string.eq.'</event') then
         maxev=maxev+1
         goto 1
      endif
      goto 1
 2    continue
      write(*,*) ' Found ',maxev,' events in file ',file
      if (maxev.eq.0) then
         write(*,*) ' NO EVENTS!! Program exits'
         call exit(3)
      endif
      call pwhg_io_close(iun)
      end subroutine

      end

  

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C---- To call the analysis and fill the histrograms     
      subroutine herwiganalysis()
      implicit none
      include 'hepevt.h'
      include 'LesHouches.h'
      if (mod(nevhep,2000).eq.0) then
         write(*,*) 'analyzing nevhep', nevhep
         call flush(6) 
      endif
      !write(*,*)'ratio'
      !write(*,*) weight, xwgtup, weight/xwgtup
!call analysis(xwgtup)
      call analysis(xwgtup)
      call pwhgaccumup 
      end    

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C---  To write the histograms in a top file     
      subroutine herwig7_end(runNumber)
      implicit none
      integer runNumber, j
      character*4 chRunNumber
      character*30 prefix
      call pwhgsetout
      if (runNumber == 0) then
        prefix = 'pwg'
      else
        write(chRunNumber,'(i4)') runNumber
        do j=1,4
          if(chRunNumber(j:j).eq.' ') chRunNumber(j:j)='0'
        enddo
        prefix = 'pwg-'//trim(chRunNumber)//'-'
      endif      
      prefix=adjustl(prefix)      
      call pwhgtopout(trim(prefix)//'POWHEG+HERWIG7-output')
      flush(6)
      end
