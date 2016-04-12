      implicit none
      include 'LesHouches.h'
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_rad.h'
      include 'pwhg_st.h'
      include 'pwhg_kn.h'
      include 'pwhg_rnd.h'
      include 'pwhg_flg.h'
      include 'pwhg_par.h'
      include 'pwhg_weights.h'
      include 'pwhg_lhrwgt.h'
      integer j,iun,iunin,iunrwgt,nev,maxev
      character * 10 statlhe
      common/cnev/nev
      real * 8 weight,tmp
      real * 8 powheginput
      character * 20 pwgprefix
      integer lprefix
      common/cpwgprefix/pwgprefix,lprefix
      integer ios
      character * 6 WHCPRG
      character * 100 filename
      common/cWHCPRG/WHCPRG
      integer iseed,n1,n2
      logical testplots

c Print out svn information, if any
      iun = 6
      include 'svn.version'


      if (powheginput('#testplots').eq.1d0) then
         testplots=.true.
      else
         testplots=.false.
      endif
      nev=powheginput('numevts')

c whether to save btilde calls to set up upper bounding envelope
      if(powheginput('#storemintupb').eq.1d0) then
         flg_storemintupb = .true.
      else
         flg_storemintupb = .false.
      endif
c whether to save btilde calls to set up upper bounding envelope
      if(powheginput('#fastbtlbound').eq.1d0) then
         flg_fastbtlbound = .true.
      else
         flg_fastbtlbound = .false.
      endif

      call newunit(iun)

c The following allows to perform multiple runs with
c different random seeds in the same directory.
c If manyseeds is set to 1, the program asks for an integer j;
c The file 'pwgprefix'seeds.dat at line j is read, and the
c integer at line j is used to initialize the random
c sequence for the generation of the event.
c The event file is called 'pwgprefix'events-'j'.lhe
      if(powheginput("#manyseeds").eq.1) then

         par_maxseeds=powheginput("#maxseeds")
         if(par_maxseeds < 0) then
            par_maxseeds = 200
         endif

         open(unit=iun,status='old',iostat=ios,
     1        file=pwgprefix(1:lprefix)//'seeds.dat')
          if(ios.ne.0) then
             write(*,*) 'option manyseeds required but '
             write(*,*) 'file ',pwgprefix(1:lprefix)/
     $            /'seeds.dat not found'
            call exit(-1)
         endif 
         do j=1,1000000
            read(iun,*,iostat=ios)  rnd_initialseed
            if(ios.ne.0) goto 10
         enddo
 10      continue
         rnd_numseeds=j-1
         write(*,*) 'enter which seed'
         read(*,*) rnd_iwhichseed
         if(rnd_iwhichseed.gt.rnd_numseeds) then
            write(*,*) ' no more than ',rnd_numseeds, ' seeds in ',
     1           pwgprefix(1:lprefix)//'seeds.dat'
            call exit(-1)
         endif
         if(rnd_iwhichseed.gt.par_maxseeds) then
            write(*,*)
     1           ' maximum seed value exceeded ',
     2           rnd_iwhichseed, '>', par_maxseeds
            write(*,*) ' Add to the powheg.input file a line like'
            write(*,*) ' maxseeds <maximum seed you need>'
            call exit(-1)
         endif
         rewind(iun)
         do j=1,rnd_iwhichseed
c Commented line to be used instead, for testing that manyseed runs
c yield the same results as single seed runs, provided the total number
c of calls is the same.
c     read(iun,*) rnd_initialseed,rnd_i1,rnd_i2
            read(iun,*) rnd_initialseed
            rnd_i1=0
            rnd_i2=0
         enddo
         close(iun)
         write(rnd_cwhichseed,'(i4)') rnd_iwhichseed
         do j=1,4
            if(rnd_cwhichseed(j:j).eq.' ') rnd_cwhichseed(j:j)='0'
         enddo
      else
         rnd_cwhichseed='none'
      endif
c If multiple weights may be used in the analysis, set them
c initially to 0 (no multiple weights). These are normally used
c only when reading lh files containing multiple weights.
      weights_num=0
c
      if (testplots) WHCPRG='NLO   '
      call pwhginit
      if(nev.gt.0) then
         if(powheginput("#clobberlhe").eq.1) then
            statlhe = 'unknown'
         else
            statlhe = 'new'
         endif
         if(flg_newweight) then
            if (testplots) then 
               write(*,*) '-------> Warning: testplots has been reset to
     1 false since we are doing reweighting' 
               testplots = .false. 
            endif
            continue
         else
            if(rnd_cwhichseed.ne.'none') then
               write(*,*) pwgprefix(1:lprefix)//'events-'//
     1            rnd_cwhichseed//'.lhe', rnd_iwhichseed,rnd_initialseed
               open(unit=iun,status=statlhe,file=pwgprefix(1:lprefix)
     1              //'events-'//rnd_cwhichseed//'.lhe')
            else
               open(unit=iun,status=statlhe,
     1              file=pwgprefix(1:lprefix)//'events.lhe')
            endif
         endif
      else
         write(*,*) ' No events requested'
         goto 999
      endif

c Input the string variables for the standard xml reweight format
      call getreweightinput

      if(.not.flg_newweight) then
         call lhefwritehdr(iun)
      endif
      if (testplots) then
         call init_hist 
c     let the analysis subroutine know that it is run by this program
         WHCPRG='LHE   '
      endif
c if we are using manyseeds, and iseed is given, it means that we want
c to examine that event in particular
      if(rnd_cwhichseed.ne.'none') then
         iseed=powheginput('#iseed')
         n1=powheginput('#rand1')
         n2=powheginput('#rand2')
         if(iseed.ge.0.and.n1.ge.0.and.n2.ge.0)
     1        call setrandom(iseed,n1,n2)
      endif
      call resetcnt
     1       ('upper bound failure in inclusive cross section')
      call resetcnt
     1       ('vetoed calls in inclusive cross section')
      call resetcnt(
     1 'upper bound failures in generation of radiation')
      call resetcnt('vetoed radiation')
      write(*,*)
      write(*,*)' POWHEG: generating events'
      if(flg_newweight) then
         flg_fullrwgt = powheginput("#fullrwgt") .eq. 1
         call opencountunit(maxev,iunin)
         if(flg_fullrwgt) then
c the following reads the pdf used in the input .lhe file;
c this is needed for fullrwgt.
            call readpowheginputinfo(iunin)
         endif
         call openoutputrw(iunrwgt)
         if(lhrwgt_id.ne.' ') then
            call lhrwgt_copyheader(iunin,iunrwgt) 
         endif
         if(maxev.ne.nev) then
            write(*,*) ' Warning: powheg.input says ',nev,' events'
            write(*,*) ' the file contains ', maxev, ' events'
            write(*,*) ' Doing ',maxev,' events'
            nev = maxev 
         endif
      endif

      flg_noevents = powheginput("#noevents") .eq. 1

      if(flg_noevents) then
         testplots = .true.
         write(*,*) 
     1' Since noevents is specified, testplots will be produced'
         write(*,*) ' irrespective of the testplot flag setting.'
      endif

      do j=1,nev
         if(j==68) then
            write(*,*) ' one hundreth event'
         endif
         if(flg_newweight) then
            call pwhgnewweight(iunin,iunrwgt)
         else
            call pwhgevent
            if(nup.eq.0) then
               write(*,*) ' nup = 0 skipping event'
               goto 111
            endif
            call lhefwritev(iun)
         endif
         if(idwtup.eq.3) then
            weight=rad_totgen*xwgtup*rad_branching
         elseif(idwtup.eq.-4) then
            weight=xwgtup
         else
            write(*,*) ' only 3 and -4 are allowed for idwtup'
            call exit(-1)
         endif
         if(testplots) then
            call lhtohep
            call analysis(weight)
            call pwhgaccumup
            if (mod(j,5000).eq.0) then
               if(rnd_cwhichseed.eq.'none') then
                  filename=pwgprefix(1:lprefix)//
     1                 'pwhgalone-output'
               else
                  filename=pwgprefix(1:lprefix)//
     1                 'pwhgalone-output'//rnd_cwhichseed
               endif
               call pwhgsetout
               call pwhgtopout(filename)
            endif
         endif
 111     continue
      enddo
      if (testplots) then
         if(rnd_cwhichseed.eq.'none') then
            filename=pwgprefix(1:lprefix)//
     1           'pwhgalone-output'
         else
            filename=pwgprefix(1:lprefix)//
     1           'pwhgalone-output'//rnd_cwhichseed
         endif
         call pwhgsetout
         call pwhgtopout(filename)
      endif
      if(flg_newweight) then
         call lhefwritetrailernw(iunin,iunrwgt)
         close(iunin)
         close(iunrwgt)
      else
         call lhefwritetrailer(iun)
         close(iun)
      endif
 999  continue
      call write_counters
c this causes powheginput to print all unused keywords
c in the powheg.input file; useful to catch mispelled keywords
      tmp=powheginput('print unused tokens')
      end


      subroutine opencountunit(maxev,iun)
      implicit none
      include 'pwhg_rnd.h'
      integer maxev,iun
      character * 30 file
      character * 20 pwgprefix
      integer lprefix
      common/cpwgprefix/pwgprefix,lprefix
      integer ios
      character * 7 string
      real * 8 powheginput
      external powheginput
      integer nev,j
      call newunit(iun)
      maxev=0
      file='pwgevents.lhe'
      open(unit=iun,file=file,status='old',iostat=ios)
      if(ios.ne.0) then
         do j=30,0,-1
            if(file(j:j).ne.' ') exit
         enddo
         write(*,*)' file not found:',file(1:j)
         write(*,*)' enter name of event file'
         read(*,'(a)') file
         open(unit=iun,file=file,status='old',iostat=ios)
         if(ios.ne.0) then
            write(*,*) 'cannot open; aborting ...'
            call exit(-1)
         endif
c get the name prefix
         j=index(file,'events')
         if(j.gt.1) then
            pwgprefix=file(1:j-1)
            lprefix=j-1
         else
            lprefix=0
         endif
c or a sequence number for manyseeds usage
         if(file(j+6:j+6).eq.'-') then
            rnd_cwhichseed=file(j+7:j+10)
         else
            rnd_cwhichseed='none'
         endif
      else
         pwgprefix='pwg'
         lprefix=3
         rnd_cwhichseed='none'
      endif
      write(*,*) 'prefix: "'//pwgprefix(1:lprefix)//'",', ' seed: "'//
     1     rnd_cwhichseed//'"'
c
      write(*,*) ' Opened event file ',file
      write(*,*) ' Counting events in ', file
      write(*,*) ' This may take some time...'
 1    continue
      read(unit=iun,fmt='(a)',end=2) string
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
      rewind(iun)
      end
