      subroutine wrtpowheginput(nlf)
      implicit none
      integer nlf
      include 'pwhg_pwin.h'
      character * (maxlin) line
      character * 20 pwgprefix
      integer lprefix
      common/cpwgprefix/pwgprefix,lprefix
      integer ios,iun,iret
      logical inrwgt
      call newunit(iun)
      if(pwgprefix(1:lprefix).eq.'pwg') then
         open(unit=iun,file='powheg.input',
     1     status='old',iostat=ios)
      else
         open(unit=iun,file=pwgprefix(1:lprefix)
     1     //'powheg.input',status='old',iostat=ios)
      endif
      if(ios.ne.0) then
         write(*,*) ' cannot open powheg input file'
         call exit(-1)
      endif
      inrwgt = .false.
      do
         read(unit=iun,fmt='(a)',iostat=ios,end=999) line
         if(adjustl(line) == '</initrwgt>') then
            read(unit=iun,fmt='(a)',iostat=ios,end=999) line
            inrwgt = .false.
         elseif(adjustl(line) == '<initrwgt>') then
            inrwgt = .true.
         endif
         if(inrwgt) cycle
         if(ios.ne.0) then
            write(*,*) ' cannot read powheg input file'
            call exit(-1)
         endif
         call pwhg_io_write(nlf,trim(line))
      enddo
 999  continue
      close(iun)
      end

      subroutine powheginputfile(file)
      implicit none
      character * 20 pwgprefix
      integer lprefix
      common/cpwgprefix/pwgprefix,lprefix
      character(len=*) file
      if(pwgprefix(1:lprefix).eq.'pwg') then
         file = 'powheg.input'
      else
         file = pwgprefix(1:lprefix)//'powheg.input'
      endif
      end

      function powheginput(stringa)
      implicit none
      real * 8 powheginput
      character *(*) stringa
      include 'pwhg_pwin.h'
      character * (maxlin) line,line0
      character * (maxkey) string
      character * 20 pwgprefix
      integer lprefix
      common/cpwgprefix/pwgprefix,lprefix
      logical exist
      integer ios,iun,j,k,l,imode,iret
      integer ini
      data ini/0/      
      save ini
      if(ini.eq.0) then
         call newunit(iun)
         exist = .false.
         if(lprefix.gt.0.and.lprefix.le.len(pwgprefix)) then
            inquire(file=pwgprefix(1:lprefix)//'powheg.input',
     1           exist=exist)
            if(exist) then
                write(*,*) ' found input file ',
     1              pwgprefix(1:lprefix)//'powheg.input'
            else
               write(*,*) ' file '//pwgprefix(1:lprefix)//
     1              'powheg.input does not exist'
            endif
              
         endif
         if(exist) then
            open(unit=iun,file=pwgprefix(1:lprefix)//'powheg.input',
     1           status='old',iostat=ios)
            if(ios.ne.0) then         
               write(*,*) ' cannot open ',
     1              pwgprefix(1:lprefix)//'powheg.input'
               call exit(-1)
            else
               write(*,*) ' opened '//pwgprefix(1:lprefix)//
     1              'powheg.input'
            endif
         else
            inquire(file='powheg.input',exist=exist)
            if(exist) then
               write(*,*) ' found input file powheg.input'
               open(unit=iun,file='powheg.input',
     1              status='old',iostat=ios)
               if(ios.ne.0) then
                  write(*,*) ' cannot open powheg.input'
                  call exit(-1)
               else
                  write(*,*) ' opened powheg.input'
               endif
               lprefix=3
               pwgprefix='pwg'
            else
               write(*,*) ' file powheg.input does not exists'
               write(*,*)
     1          ' Enter the prefix for this run < 20 characters'
               read(*,'(a)') pwgprefix
               do lprefix=20,1,-1
                  if(pwgprefix(lprefix:lprefix).ne.' ') then
                     goto 11
                  endif
               enddo
 11            continue
               lprefix=lprefix+1
               if(lprefix.gt.20) lprefix=20
               pwgprefix(lprefix:lprefix)='-'
               open(unit=iun,file=pwgprefix(1:lprefix)//'powheg.input',
     1              status='old',iostat=ios)
               if(ios.ne.0) then            
                  write(*,*) ' cannot open ',
     1                 pwgprefix(1:lprefix)//'powheg.input'
                  call exit(-1)
               else
                  write(*,*) ' opened ',
     1                 pwgprefix(1:lprefix)//'powheg.input'
                  
               endif
            endif
         endif
         pwin_numvalues=0
         do l=1,1000000
            line0=' '
            read(unit=iun,fmt='(a)',iostat=ios) line0
            if(ios.ne.0.and.line0.eq.' ') goto 10
c skip completely a reweigting section
            if(adjustl(line0) == '<initrwgt>') then
               do while(adjustl(line0) /= '</initrwgt>')
c     we empty the line before reading, since in case of EOF the behaviour
c     is undefined
                  line0 = ' '
                  read(unit=iun,fmt='(a)',iostat=ios) line0
                  if(ios.ne.0.and.line0.eq.' ') goto 10
               enddo
               line0 = ' '
               read(unit=iun,fmt='(a)',iostat=ios) line0
               if(ios.ne.0.and.line0.eq.' ') goto 10
            endif
c end reweigting section
            line=line0
            do k=1,maxlin
               if(line(k:k).eq.'#'.or.line(k:k).eq.'!') then
                  line(k:)=' '
               endif
            enddo
            if(line.ne.' ') then
               if(pwin_numvalues.eq.maxnum) then
                  write(*,*) ' too many entries in powheginput.dat'
                  call exit(-1)
               endif
               pwin_numvalues=pwin_numvalues+1
c get first word in line
               call firststringword(line,pwin_keywords(pwin_numvalues),
     1              iret)
               if(iret.lt.0) then
                  write(*,*) ' powheginput: keyword too long:'
                  write(*,*) trim(line)
                  call exit(-1)
               endif
c See if the same keyword is already there: give error in this case
               if(pwin_numvalues.gt.1) then
                  do j=1,pwin_numvalues-1
                     if(pwin_keywords(j).eq.
     1                    pwin_keywords(pwin_numvalues)) then
                        write(*,*) 'powheginput: keyword '//
     1                       trim(pwin_keywords(pwin_numvalues))
     1                     //' appears more than once in powheg.input:'
                        write(*,*) trim(line)
                        write(*,*) ' appeared after '
                        write(*,*)  pwin_keywords(j), pwin_values(j)
                        write(*,*) 'Exiting'
                        call exit(-1)
                     endif
                  enddo
               endif
               call skipfirstword(line,line,iret)
c if the first character is a quote, it is a string
               if(line(1:1).eq."'".or.line(1:1).eq.'"') then
                  pwin_numstrings = pwin_numstrings + 1
                  if(pwin_numstrings.gt.maxstrings) then
                     write(*,*) ' powheginput: too many strings'
                     write(*,*) ' increase maxstrings'
                     call exit(-1)
                  endif 
                  pwin_stringptr(pwin_numvalues) = pwin_numstrings
                  call getquotedstring(line,
     1                 pwin_strings(pwin_numstrings),iret)
                  if(iret.lt.0) then
                     write(*,*) ' powheginput: cannot read string'
                     call exit(-1)
                  endif
                  pwin_values(pwin_numvalues)=-1d6
                  pwin_used(pwin_numvalues)=.false.
               else
                  pwin_stringptr(pwin_numvalues)=0
                  read(unit=line,fmt=*,iostat=ios)
     1                 pwin_values(pwin_numvalues)
                  pwin_used(pwin_numvalues)=.false.
                  if(ios.ne.0) then
                     write(*,*) ' powheginput error: cannot parse '
                     write(*,'(a)') line0
                     stop
                  endif
               endif
            endif
         enddo
 10      continue
         close(iun)
         ini=1
      endif
      if(stringa.eq.'initialize powheginput') then
         return
      endif
      if(stringa.eq.'print unused tokens') then
         do j=1,pwin_numvalues
            if(.not.pwin_used(j)) then
               write(*,*)'powheginput WARNING: unused variable ',
     1              pwin_keywords(j)
            endif
         enddo
         return
      endif

      call  assignstring(stringa,string,iret)
      if(iret.lt.0) then
         write(*,*) ' powheginput: input string too long:',trim(stringa)
         call exit(-1)
      endif

      if(string(1:1).eq.'#') then
         string=string(2:)
         imode=0
      else
         imode=1
      endif
      do j=1,pwin_numvalues
         if(string.eq.pwin_keywords(j)) then
            powheginput=pwin_values(j)
            if(.not.pwin_used(j)) then
               pwin_used(j)=.true.
               write(*,*) ' powheginput keyword ',pwin_keywords(j),
     1                    ' set to ',pwin_values(j)
            endif
            return
         endif
      enddo
      if(imode.eq.1) then
         write(*,*) ' powheginput: keyword ',string,' not found'
         call exit(-1)
      endif
c Not found; assign value -1d6; store the token anyhow
      if(pwin_numvalues.eq.maxnum) then
         write(*,*) ' too many entries in powheginput.dat'
         write(*,*) ' increase maxnum in powheginput.f'
         call exit(-1)
      endif
      pwin_numvalues=pwin_numvalues+1
      pwin_keywords(pwin_numvalues)=string
      pwin_values(pwin_numvalues)=-1d6
      pwin_used(pwin_numvalues)=.true.
      powheginput=-1d6
      write(*,*) ' powheginput keyword ',pwin_keywords(j),
     1     ' absent; set to ',pwin_values(j)
      end

      subroutine powheginputstring(stringa,stringout)
      implicit none
      character * (*) stringa,stringout
      include 'pwhg_pwin.h'
      character * (maxkey) string
      real * 8 tmp,powheginput
      integer j,imode,iret
      logical ini
      data ini/.true./
      save ini
      if(ini) then
c just force loading of the powheg.input file
         tmp = powheginput('initialize powheginput')
         ini = .false.
      endif

      if(stringa(1:1).eq.'#') then
         imode = 0
      else
         imode = 1
      endif

      call assignstring(stringa(2-imode:),string,iret)
      if(iret.lt.0) then
         write(*,*) ' powheginputstring: '
         write(*,*) ' keyword too long'
         call exit(-1)
      endif

      do j=1,pwin_numvalues
         if(string.eq.pwin_keywords(j)) then
            if(pwin_stringptr(j) == 0) then
               write(*,*) 'powheginputstring: error, keyword ',
     1              trim(pwin_keywords(j)),
     2              ' is not associated to a string'
               write(*,*) 'exiting ...'
               call exit(-1)
            endif
            call assignstring(pwin_strings(pwin_stringptr(j)),
     1           stringout,iret)
            if(iret.lt.0) then
               write(*,*) ' powheginputstring:'
               write(*,*) ' output string too short'
               call exit(-1)
            endif
            pwin_used(j)=.true.
            goto 999
         endif
      enddo
      if(imode.eq.0) then
         stringout = ' '
      else
         write(*,*) ' powheginputstring: '
         write(*,*) ' keyword '//trim(string)//' not present!'
         call exit(-1)
      endif
 999  continue
      write(*,*) ' powheginput keyword ',string,
     1                    ' set to ','"'//trim(stringout)//'"'

      end


c String utilities follow here. This may be useful in a separate
c file ...

      subroutine assignstring(stringin,stringout,iret)
c assign string checking for overflow. Leading and trailing blanks are ignored
c iret = 1: output string too short
c iret = 0 OK
      implicit none
      character * (*) stringin,stringout
      integer iret
      integer lin,lout
      integer j
      lin = len(trim(adjustl(stringin)))
      lout = len(stringout)
      if(lout.lt.lin) then
         write(*,*)
     1        'assignstring: input string does not fit in output string'
         iret = -1
         return
      endif
      stringout = adjustl(stringin)
      iret = 0
      end

      subroutine firststringword(stringin,word,iret)
      implicit none
      character * (*) stringin,word
      integer iret
      integer first,last
      integer j,l
      l=len(stringin)
      first = 0
      last = 0
      do j=1,l
         if(stringin(j:j).ne.' ') then
            if(first.eq.0) first = j
         else
            if(first.ne.0) then
               last = j-1
               exit
            endif
         endif
      enddo
      if(first.eq.0) then
         iret = 1
         word = ' '
         return
      endif
      if(last.eq.0) last = l
      call assignstring(stringin(first:last),word,iret)
      if(iret.lt.0) then
         write(*,*) ' firststringword: output word too short'
         return
      endif
      end
            
      subroutine skipfirstword(stringin,stringout,iret)
      implicit none
      character * (*) stringin,stringout
      integer iret
      integer first,last
      integer j,l
      l=len(stringin)
      first = 0
      last = 0
      do j=1,l
         if(stringin(j:j).ne.' ') then
            if(first.eq.0) first = j
         else
            if(first.ne.0) then
               last = j-1
               exit
            endif
         endif
      enddo
      if(first.gt.0.and.last.eq.0) last = l
      if(last.lt.l) then
         call assignstring(stringin(last+1:),stringout,iret)
         if(iret.lt.0) then
            write(*,*) 'skipfirstword: output string too short'
         endif
      else
         stringout = ' '
         iret = 1
      endif
      end

      subroutine getquotedstring(stringin,stringout,iret)
      character *(*) stringin,stringout
      integer iret
      integer lin,lout,lclose
      character * 1 quote 
      lin=len(stringin)
      lout=len(stringout)
      stringout = adjustl(stringin)
      quote = stringout(1:1)
      if(quote.ne."'".and.quote.ne.'"') then
         write(*,*) ' ************* ERROR ***********'
         write(*,*) ' getquotedstring: did not find a quote in string'
         iret = -1
         return
      endif
      lclose = index(stringout(2:),quote) + 1
      if(lclose.eq.0) then
         write(*,*) ' ************* ERROR ***********'
         write(*,*) ' getquotedstring: did not find end quote in string'
         iret = -2
         return
      endif
      stringout = stringout(2:lclose-1)
      iret = 0
      end

      function stringlength(string)
c length of string neglecting leading and trailing blanks
      integer stringlength
      character * (*) string
      stringlength = len(trim(adjustl(string)))
      end
