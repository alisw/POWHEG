      subroutine rwl_init
      implicit none
      include 'pwhg_rwl.h'
      rwl_initialized = rwl_initialized_const
      rwl_num_weights = 0
      rwl_num_groups = 0
      end

      subroutine rwl_loadweights(iun,iret)
      use, intrinsic :: ISO_C_BINDING
      implicit none
      include 'pwhg_rwl.h'
      integer iun,iret
      character(len=9) tag
      character(len=:), allocatable :: buf
      integer nskip,jpos,jw,i,j
      logical start_equal_strings
c make sure we have enough space in rwl_weights array
      call rwl_allocate_weights
      call pwhg_io_skip_until(iun,'<',nskip)
      call pwhg_io_skip(iun,-1)
      call pwhg_io_read(iun,tag,iret)
      if(iret /= 0) call errr("cannot read!")
      call pwhg_io_backspace(iun)
      if(tag(1:6) == '<rwgt>') then
         call pwhg_io_skip_until(iun,'</rwgt>',nskip)
      elseif(tag == '<weights>') then
         call pwhg_io_skip_until(iun,'</weights>',nskip)
      endif
      call pwhg_io_skip(iun,-nskip)
      if( .not. allocated(buf)) then
         allocate(character(len=nskip) :: buf)
      else
         if(len(buf) < nskip) then
            deallocate(buf)
            allocate(character(len=nskip) :: buf)
         endif
      endif
      call pwhg_io_read_buf(iun,buf(1:nskip))
      if(buf(1:6)=='<rwgt>') then
c     <wgt id='...'> number </wgt> weights
         jpos = 7
         jw = 0
         do
            j = index(buf(jpos:nskip),'<wgt ')
            if(j <= 0) exit
            jpos = jpos + j - 1 + len('<wgt ')
            if(start_equal_strings(buf(1:nskip),'id',jpos)) then
               if(start_equal_strings(buf(1:nskip),'=',jpos)) then
                  call getquotedstringpos(buf(jpos:),i,j)
               else
                  call errr("no id keyword after <wgt")
               endif
            else
               call errr("no = keyword after <wgt id")
            endif
            if(i == -1 .or. j == -1) then
               write(*,*) buf
               write(*,*) buf(jpos:)
               call errr("no id string in <weight  >")
            endif
            jpos = jpos + j + 1
            if( .not. start_equal_strings(buf,'>',jpos)) then
               write(*,*) buf
               write(*,*) buf(jpos:)
               call errr("no id string in <weight >")
            endif
            jw = jw + 1
            if(jw > rwl_num_weights)
     1           call errr("too many weights in event!")
            j=index(buf(jpos:),'</wgt>')
            read(buf(jpos:jpos+j-2),*) rwl_weights(jw)
            jpos = jpos + j - 1 + len('</wgt>')
            if(jpos >= nskip) exit
         enddo
c         if(jw /= rwl_num_weights)
c     1           call errr("not enough weights in event!")
      else
c     buf='<weights> ... </weights>'
         jpos = 10
         jw = 0
         do while (jpos < nskip - len('</weights>'))
            call next_word_in_string(buf,jpos,i,j)
            if(i<0) exit
            jw = jw + 1
            if(jw > rwl_num_weights)
     1           call errr("too many weights in event!")
            read(buf(i:j),*) rwl_weights(jw)
            if(buf(jpos:jpos) == '<') exit
         enddo
c         if(jw /= rwl_num_weights)
c     1           call errr("not enough weights in event!")
      endif
      contains
      subroutine errr(string)
      character(len=*) string
      write(*,'(a)') 'rwl_loadweights: '//string
      call exit(-1)
      end subroutine
      end

      subroutine rwl_readheader(iunin)
      implicit none
      include 'pwhg_rwl.h'
      character(len=:), allocatable :: rwgtbuf
      integer iunin
      INTEGER NCH
      call pwhg_io_skip_until(iunin,'<initrwgt>',nch)
      if(nch <=0 ) call errr("no <initrwgt> in lhe file ")
c     now the file is pointing right after <initrwgt>
      call pwhg_io_skip_until(iunin,'</initrwgt>',nch)
      if(nch <=0 ) call errr('no </initrwgt> in lhe file ')
      call pwhg_io_skip(iunin,-nch)
c     back right after <initrwgt>
      nch = nch - len('</initrwgt>')
      allocate(character(nch) :: rwgtbuf)
      call pwhg_io_read_buf(iunin,rwgtbuf)
      call process_rwgt_info(rwgtbuf)
      deallocate(rwgtbuf)
      contains
      subroutine errr(string)
      character(len=*) string
      write(*,'(a)') 'rwl_readheader: '//string
      call exit(-1)
      end subroutine
      end
      
      subroutine rwl_copyheader(iunin,iunout,num_o_w,num_w)
      implicit none
      include 'pwhg_rwl.h'
      integer iunin,iunout,num_o_w,num_w
      character * (400) string,string0
      character(len=:), allocatable :: rwgtbuf
      integer j,k,iret,nch
      logical written_rwgt
      written_rwgt = .false.
      rwl_num_weights = 0
      rwl_num_groups = 0
      do j=1,1000000
         call pwhg_io_read(iunin,string0,iret)
         if(iret /= 0) goto 999
c         read(unit=iunin,fmt='(a)',end=999,err=999) string0
         string=adjustl(string0)
         if(string(1:10).eq.'<initrwgt>') then
c     we back to the previous line and then position the file right after
c     the <initrwgt>, in case something is written on the same line 
            call pwhg_io_backspace(iunin)
            call rwl_readheader(iunin)
            num_o_w = rwl_num_weights
            call rwl_setup(num_w)
            num_w = rwl_num_weights
            call rwl_write_rwgt_info(iunout)
            written_rwgt = .true. 
            call pwhg_io_skip_until(iunin,'</initrwgt>',nch)
         elseif(string == '<init>' .and. .not. written_rwgt) then
            num_o_w = 0
            call rwl_setup(num_w)
            num_w = rwl_num_weights
            call pwhg_io_write(iunout,'<header>')
            call rwl_write_rwgt_info(iunout)
            call pwhg_io_write(iunout,'</header>')
            written_rwgt = .true.
            call pwhg_io_write(iunout,trim(string0))
         elseif(string == '</init>') then
            call pwhg_io_write(iunout,trim(string0))
            exit
         else
            call pwhg_io_write(iunout,trim(string0))
c     write(iunout,'(a)') trim(string0)
         endif
      enddo
 999  continue
      contains
      subroutine errr(string)
      character(len=*) string
      write(*,'(a)') 'rwl_copyheader: '//string
      call exit(-1)
      end subroutine
      end


      subroutine rwl_write_rwgt_info(iun)
      implicit none
      include 'pwhg_rwl.h'
      include 'pwhg_flg.h'
      integer iun
      integer kg,kw,group
      call pwhg_io_write(iun,"<initrwgt>")
      do kg = 0,rwl_num_groups
         if(kg /= 0) then
            call pwhg_io_write(iun,"<weightgroup name='"//
     1           rwl_groups_array(kg)%name//"' combine='"//
     2           rwl_groups_array(kg)%combine//"' >")
         endif
         do kw = 1,rwl_num_weights
            if(rwl_weights_array(kw)%group == kg) then
               call pwhg_io_write(iun,"<weight id='"//
     1              rwl_weights_array(kw)%id//"' >"//
     2              rwl_weights_array(kw)%desc//"</weight>")
            endif
         enddo
c     If flg_nnlops is set, weights are doubled to include the NNLO corrected ones.
c     We include here a dummy file. NNLOPS generators can override it to call their own routines
c     for filling the NNLOPS weights in the header.
         include 'rwl_write_rwgt_extra_info.f'
         if(kg /= 0) then
            call pwhg_io_write(iun,"</weightgroup>")
         endif
      enddo
      call pwhg_io_write(iun,"</initrwgt>")
      end
            
      subroutine rwl_setup(num_w)
      implicit none
      include 'pwhg_flg.h'
      include 'pwhg_rwl.h'
      integer num_w
      real * 8 powheginput
      character * 300 string
      character(len=:), allocatable :: rwgtbuf
      integer iun,iret,nch,j
c     See where to find reweight information
      call powheginputstring('#rwl_file',string)
      if(string == ' ') then
         num_w= 0
         rwl_num_weights = 0
         flg_rwl = .false.
         return
      endif
      flg_rwl = .true.
      if(string == '-') then
c     in this case the xml weight information is
c     embedded in the powheg.input file.
         call powheginputfile(string)
      endif
      call pwhg_io_open_read(trim(string),iun,iret)
      if(iret /= 0) call errr("cannot open file"//trim(string))
      call pwhg_io_skip_until(iun,'<initrwgt>',nch)
      if(nch <=0 ) call errr("no <initrwgt> in file "
     1        //trim(string))
      call pwhg_io_skip_until(iun,'</initrwgt>',nch)
      if(nch <=0 ) call errr('no </initrwgt> in file '
     1        //trim(string))
      call pwhg_io_skip(iun,-nch)
      nch = nch - len('</initrwgt>')
      allocate(character(nch) :: rwgtbuf)
      call pwhg_io_read_buf(iun,rwgtbuf)
      call pwhg_io_close(iun)
      call process_rwgt_info(rwgtbuf)
      deallocate(rwgtbuf)
      num_w = rwl_num_weights
      contains
      subroutine errr(string)
      character(len=*) string
      write(*,'(a)') 'rwl_setup: '//string
      call exit(-1)
      end subroutine
      end

      subroutine process_rwgt_info(buff)
      implicit none
      character(len=*) buff
      integer l,jpos,j,k,group
      logical start_equal_strings
      l=len(buff)
      group=0
      jpos = 1
      do while(jpos <= l)
         j=index(buff(jpos:),'<')
         if(j>0 .and. j<l-jpos+1) then
            jpos=j+jpos
            if(start_equal_strings(buff,'/weightgroup>',jpos))
     1           then
               group = 0
               cycle
            endif
            if(.not. start_equal_strings(buff,'weight',jpos))
     1           then
               cycle
            endif
            if(start_equal_strings(buff,'group',jpos)) then
c     it is a group
               call rwgt_info_addgroup(buff,jpos,group)
               if(buff(jpos-1:jpos-1)/='>') then
                  call errr('This should not happen ...')
               endif
               cycle
            endif
c     it is a weight
            call rwgt_info_addweight(buff,jpos,group)
            if(buff(jpos-1:jpos-1)/='>') then
               call errr('This should not happen ...')
            endif
         else
            exit
         endif
      enddo
      call rwgt_print_weights
      contains
      subroutine errr(string)
      character(len=*) string
      write(*,'(a)') 'process_rwgt_info: '//string
      call exit(-1)
      end subroutine
      end

      subroutine rwgt_print_weights
      implicit none
      include 'pwhg_rwl.h'
      integer j,k
      character(len=3) group
      do j=1,rwl_num_weights
         write(*,'(a,i4)') 'Weight ',j
         write(*,'(a)') '   id='//rwl_weights_array(j)%id//','
         write(group,'(i3)')  rwl_weights_array(j)%group
         write(*,'(a,i2)') '   group='//trim(adjustl(group))//','
         if(rwl_weights_array(j)%num_keys > 0) then
            write(*,'(a)') '   key value pairs:'
            do k=1,rwl_weights_array(j)%num_keys
               write(*,'(a,d11.5)') '      '//rwl_weights_array(j)%desc(
     1              rwl_weights_array(j)%keys(1,k):
     2              rwl_weights_array(j)%keys(2,k))//'=',
     3              rwl_weights_array(j)%values(k)
            enddo
         endif
      enddo
      end

      subroutine rwl_write_weights(iun)
      implicit none
      integer iun
      include 'pwhg_rwl.h'
      include 'pwhg_flg.h'
      include 'pwhg_rwgtsuda.h'
      character(len=20) string
      character(len=11) tmpstr
      integer jg,jw
      if(rwl_format_rwgt) then
         call pwhg_io_write(iun,'<rwgt>')
         do jg=0,rwl_num_groups
            do jw = 1, rwl_num_weights
               if(rwl_weights_array(jw)%group == jg) then
                  write(tmpstr,'(E11.5)') rwl_weights(jw)*sudawgt
                  call pwhg_io_write(iun,"<wgt id='"
     1                 //rwl_weights_array(jw)%id//"'>"
     2                 //tmpstr//'</wgt>')
               endif
            enddo
         enddo
         call pwhg_io_write(iun,'</rwgt>')
      else
         call pwhg_io_write(iun,'<weights>')
         do jg=0,rwl_num_groups
            do jw = 1, rwl_num_weights
               if(rwl_weights_array(jw)%group == jg) then
                  write(tmpstr,'(E11.5)') rwl_weights(jw)*sudawgt
                  call pwhg_io_write(iun,tmpstr)
               endif
            enddo
c     If flg_nnlops is set, weights are doubled to include the NNLO corrected ones.
c     We include here a dummy file. NNLOPS generators can override it to call their own routines
c     for writing the NNLOPS weights in the event.            
            include 'rwl_write_weights_extra.f'
         enddo
         call pwhg_io_write(iun,'</weights>')
      endif
      end

      subroutine rwgt_info_addgroup(buf,jpos,group)
      implicit none
      include 'pwhg_rwl.h'
      character(len=*) :: buf
      integer jpos,group
      character ch
      integer l,i,j,k,iname,icomb
      logical start_equal_strings
      if(rwl_initialized /= rwl_initialized_const) call rwl_init
      l=len(buf)
      rwl_num_groups = rwl_num_groups + 1
      if(rwl_num_groups > rwl_maxgroups) call errr
     1 ("too many groups")
      group = rwl_num_groups
      iname=0
      icomb=0
c do this twice, see if we find name and combine in any order
      do k=1,2
         if(start_equal_strings(buf,'name',jpos)) then
            if(iname==1) call errr
     1           ("found more than 1 name in weightgroup")
            iname = 1
            if(start_equal_strings(buf,'=',jpos)) then
               call getquotedstringpos(buf(jpos:),i,j)
               if(i==-1) call errr("did not find quote after =")
               if(j==-1) call errr("did not find end quote after =")
               i = i+jpos-1
               j = j+jpos-1
               allocate(character(j-i+1)::
     1              rwl_groups_array(rwl_num_groups)%name)
               rwl_groups_array(rwl_num_groups)%name=buf(i:j)
               jpos = j+2
            else
               call errr
     1           ("did not find = after name in weightgroup")
            endif
         endif
         if(start_equal_strings(buf,'combine',jpos)) then
            if(icomb==1) call errr
     1           ("found more than 1 combine in weightgroup")
            icomb = 1
            if(start_equal_strings(buf,'=',jpos)) then
               call getquotedstringpos(buf(jpos:),i,j)
               if(i==-1) call errr("did not find quote after =")
               if(j==-1) call errr("did not find end quote after =")
               i = i+jpos-1
               j = j+jpos-1
               allocate(character(j-i+1)::
     1              rwl_groups_array(rwl_num_groups)%combine)
               rwl_groups_array(rwl_num_groups)%combine=buf(i:j)
               jpos = j+2
            else
               call errr
     1           ("did not find = after name in weightgroup")
            endif
         endif
      enddo
      if( .not. start_equal_strings(buf,'>',jpos)) then
         call errr
     1        ("found junk in weightgroup tag")
      endif
      if(iname == 0) call errr("Did not find name in weightgroup")
      if(icomb == 0) then
         allocate(character(1)::
     1        rwl_groups_array(rwl_num_groups)%combine)
         rwl_groups_array(rwl_num_groups)%combine=' '
      endif
c     Check that group is not already present
      do k=1,rwl_num_groups-1
         if(rwl_groups_array(rwl_num_groups)%name ==
     1        rwl_groups_array(k)%name) then
            group = k
            deallocate(rwl_groups_array(rwl_num_groups)%name,
     1           rwl_groups_array(rwl_num_groups)%combine)
            rwl_num_groups = rwl_num_groups - 1
            exit
         endif
      enddo
      contains
      subroutine errr(string)
      character(len=*) string
      write(*,'(a)') 'rwgt_info_addgroup: '//string
      call exit(-1)
      end subroutine
      end

      subroutine rwgt_info_addweight(buf,jpos,group)
      implicit none
      include 'pwhg_rwl.h'
      character(len=*) :: buf
      integer jpos,group
      integer l,i,j,k
      real * 8 val
      logical  next_key_value_pair,start_equal_strings
      if(rwl_initialized /= rwl_initialized_const) call rwl_init
      l=len(buf)
      rwl_num_weights = rwl_num_weights + 1
      if(rwl_num_weights>rwl_maxweights) call
     1 errr(' rwgt_info_addweight: too many weights,'//
     2  ' increase maxweights')
      if(start_equal_strings(buf,'id',jpos)) then
         if(start_equal_strings(buf,'=',jpos)) then
            call  getquotedstringpos(buf(jpos:),i,j)
            if(i == -1) call errr
     1 ("can't find beginning of quoted string after id= in a weight")
            if(j == -1) call errr
     1       ("can't find end of quoted string after id= in a weight")
            i=i+jpos-1
            j=j+jpos-1
            allocate(character(j-i+1)::
     1           rwl_weights_array(rwl_num_weights)%id)
            rwl_weights_array(rwl_num_weights)%id=buf(i:j)
            jpos=j+2
            if(start_equal_strings(buf,'>',jpos)) then
               i=jpos
               j = index(buf(jpos:),'</weight>')
               if(j>0) then
                  j=jpos+j-2
                  jpos=j+len('</weight>')+1
                  allocate(character(j-i+1) ::
     1                 rwl_weights_array(rwl_num_weights)%desc)
                  rwl_weights_array(rwl_num_weights)%desc=buf(i:j)
               else
                  call errr ("can't find </weight> after <weight>")
               endif
               rwl_weights_array(rwl_num_weights)%group=group
            else
               call errr ("can't find > after <weight")
            endif
         else
            call errr ("can't find '=' after 'id' in weight")
         endif
      else
         call errr ("can't find 'id' in weight")
      endif
c     Check that weight is not already present
      do k=1,rwl_num_weights-1
         if(rwl_weights_array(rwl_num_weights)%id ==
     1        rwl_weights_array(k)%id) call errr
     2        ('This weight id is already present:'//
     3        rwl_weights_array(k)%id)
      enddo
c     Set up arrays of key-value pairs contained in the desc string
      k=1
      l = 0
      do while(next_key_value_pair(
     1    rwl_weights_array(rwl_num_weights)%desc,k,i,j,val))
         l = l + 1
      enddo
      allocate(rwl_weights_array(rwl_num_weights)%keys(2,l),
     1     rwl_weights_array(rwl_num_weights)%values(l))
      k=1
      l = 0
      do while(next_key_value_pair(
     1     rwl_weights_array(rwl_num_weights)%desc,k,i,j,val))
         l = l + 1
         rwl_weights_array(rwl_num_weights)%keys(1,l) = i
         rwl_weights_array(rwl_num_weights)%keys(2,l) = j
         rwl_weights_array(rwl_num_weights)%values(l) = val
      enddo
      rwl_weights_array(rwl_num_weights)%num_keys = l

      contains
      subroutine errr(string)
      character(len=*) string
      write(*,'(a)') 'rwgt_info_addweight: '//string
      call exit(-1)
      end subroutine
      end

      logical function next_key_value_pair(buf,jpos,iv,jv,val)
c     finds the next key=value pair in buf(jpos:)
c     when returning jpos points at the character right after the found
c     key=value pair, iv and jv are set to the index of the first and last
c     character of the key in buf.
c      If no key=value pair is found it returns false.
      implicit none
      character(len=*) buf
      integer jpos,iv,jv
      real * 8 val
      integer lb,i1,i2,ieq
      lb = len(buf)
 1    continue
      if(jpos <=0 ) then
         write(*,*) 'next_key_value_pair:'
         write(*,*) 'index of first character<=0, exiting ...'
         call exit(-1)
      endif
      if(jpos > lb) goto 998
c this should be the key
      call nextword(jpos,iv,jv)
c No key found
      if(iv>lb) goto 998
      if(buf(iv:jv) == 'default') then
c     In this case we do not expect a value
         val = 0
         next_key_value_pair = .true.
         jpos = jv+1
         return
      endif
c     Next we should find an =
      ieq=jv+1
      do
         if(ieq>lb) exit
         if(buf(ieq:ieq) /= ' ') exit
         ieq=ieq+1
      enddo
c     If there is no =, it is an error
      if(ieq>lb) goto 999
      if(buf(ieq:ieq) /= '=') goto 998
c     Look for number to read
      call nextword(ieq+1,i1,i2)
      read(buf(i1:i2),fmt=*,err=999) val
      jpos=i2+1
      next_key_value_pair = .true.
      return
 998  next_key_value_pair = .false.
      return
 999  write(*,*) 'next_key_value_pair: malformed input buffer'
      write(*,*) '"'//buf//'"'
      call exit(-1)
      contains
      subroutine nextword(jj,ii1,ii2)
      integer jj,ii1,ii2
      ii1=jj
      do
         if(ii1 > lb) exit
         if(buf(ii1:ii1) /= ' ') exit
         ii1=ii1+1
      enddo
      ii2=ii1
      do
         if(ii2 > lb) exit
         if(buf(ii2:ii2) == ' ' .or. buf(ii2:ii2) == '=') exit
         ii2=ii2+1
      enddo
      ii2 = ii2-1
      end subroutine nextword
      end

      subroutine rwl_handle_lhe(task,numevts,count)
      implicit none
      character *(*) task
      integer numevts,count
      include 'nlegborn.h'
      include 'pwhg_rwl.h'
      include 'pwhg_rad.h'
      include 'pwhg_rwgtsuda.h'
      include 'LesHouches.h'
      include 'pwhg_kn.h'
      type(rwl_lhe_block), allocatable, save :: events(:)
      logical,save :: ini=.true.
      integer j
      if(ini) then
         allocate(events(numevts))
         events(:)%nup = 0
         ini = .false.
      endif
      if(task == 'print') then
         do count=1,numevts
            print *, events(count)%nup
         enddo
         return
      endif
      if(task == 'put') then
         if(events(count)%nup == 0) then
            allocate(events(count)%idup(nup),
     1           events(count)%istup(nup),
     2           events(count)%mothup(2,nup),
     3           events(count)%icolup(2,nup),
     4           events(count)%pup(5,nup),
     5           events(count)%vtimup(nup),
     6           events(count)%spinup(nup))
         elseif(nup>events(count)%nup) then
            deallocate(events(count)%idup,
     1           events(count)%istup,
     2           events(count)%mothup,
     3           events(count)%icolup,
     4           events(count)%pup,
     5           events(count)%vtimup,
     6           events(count)%spinup)
            allocate(events(count)%idup(nup),
     1           events(count)%istup(nup),
     2           events(count)%mothup(2,nup),
     3           events(count)%icolup(2,nup),
     4           events(count)%pup(5,nup),
     5           events(count)%vtimup(nup),
     6           events(count)%spinup(nup))
         endif
         events(count)%nup    = nup
         events(count)%idprup = idprup
         events(count)%xwgtup = xwgtup
         events(count)%scalup = scalup
         events(count)%aqedup = aqedup
         events(count)%aqcdup = aqcdup
         events(count)%idup(1:nup) = idup(1:nup)
         events(count)%istup(1:nup) = istup(1:nup)
         events(count)%mothup(:,1:nup) = mothup(:,1:nup)
         events(count)%icolup(:,1:nup) = icolup(:,1:nup)
         events(count)%pup(:,1:nup) = pup(:,1:nup)
         events(count)%vtimup(1:nup) = vtimup(1:nup)
         events(count)%spinup(1:nup) = spinup(1:nup)

         events(count)%rwl_type    = rwl_type
         events(count)%rwl_index   = rwl_index
         events(count)%rwl_seed    = rwl_seed
         events(count)%rwl_n1      = rwl_n1
         events(count)%rwl_n2      = rwl_n2
         events(count)%rwl_weight  = rwl_weight

         events(count)%rad_type  = rad_type
         events(count)%rad_kinreg  = rad_kinreg

         events(count)%sudawgt = sudawgt
         events(count)%uborns = uborns

      elseif(task == 'get') then
         nup = events(count)%nup
         idprup = events(count)%idprup
         xwgtup = events(count)%xwgtup
         scalup = events(count)%scalup
         aqedup = events(count)%aqedup
         aqcdup = events(count)%aqcdup
         idup(1:nup)     = events(count)%idup(1:nup)
         istup(1:nup)    = events(count)%istup(1:nup)
         mothup(:,1:nup) = events(count)%mothup(:,1:nup)
         icolup(:,1:nup) = events(count)%icolup(:,1:nup)
         pup(:,1:nup)    = events(count)%pup(:,1:nup)
         vtimup(1:nup)   = events(count)%vtimup(1:nup)
         spinup(1:nup)   = events(count)%spinup(1:nup)

         rwl_type   = events(count)%rwl_type  
         rwl_index  = events(count)%rwl_index 
         rwl_seed   = events(count)%rwl_seed  
         rwl_n1     = events(count)%rwl_n1    
         rwl_n2     = events(count)%rwl_n2    
         rwl_weight = events(count)%rwl_weight

         rad_type   = events(count)%rad_type
         rad_kinreg = events(count)%rad_kinreg

         sudawgt    = events(count)%sudawgt
         uborns     = events(count)%uborns

      elseif(task == 'getwinfo') then
         xwgtup = events(count)%xwgtup
         rwl_type   = events(count)%rwl_type  
         rwl_index  = events(count)%rwl_index 
         rwl_seed   = events(count)%rwl_seed  
         rwl_n1     = events(count)%rwl_n1    
         rwl_n2     = events(count)%rwl_n2    
         rwl_weight = events(count)%rwl_weight

         rad_type   = events(count)%rad_type
         rad_kinreg = events(count)%rad_kinreg

      endif
      end

      
      subroutine getquotedstringpos(string,i,j)
      implicit none
      character(len=*) string
      character ch
      integer i,j
      integer l
      l=len(string)
      do i=1,l
         ch=string(i:i)
         if(ch=="'" .or. ch=='"') exit
         if(ch /= ' ') then
            write(*,*)
     1    'getquotedstringpos: was expecting a quote, found '//ch
            call exit(-1)
         endif
      enddo
      if(i>l) then
         i=-1
         return
      endif
c skip initial quote
      i = i+1
      j=index(string(i:),ch)
      if(j<=0) then
         j=-1
         return
      endif
c skip final quote
      j = (i-1)+j-1
      contains
      end


c     Finds the next word in string(jpos:), sets i and j to the first
c     and last letter of the word. The symbols < and > are always
c     considered single letter words. Words can be separated by blanks,
c     tabs and newlines.
      subroutine next_word_in_string(string,jpos,i,j)
      use, intrinsic :: ISO_C_BINDING
      implicit none
      character(len=*) string
      integer jpos,i,j
      character ch
      integer l
      l=len(string)
      do i=jpos,l
         select case(string(i:i))
         case(' ',c_new_line,c_horizontal_tab)
            continue
         case default
            exit
         end select
      enddo
      if(i>l) then
         i=-1
         return
      endif
c skip initial quote
      do j=i+1,l
         select case(string(j:j))
         case(' ',c_new_line,c_horizontal_tab,'<','>')
            exit
         end select
      enddo
      jpos = j
      j=j-1
c skip final quote
      end

      logical function start_equal_strings(str1,str2,jpos)
      use, intrinsic :: ISO_C_BINDING
c returns true if str2 is the beginning of str1(jpos:) ignoring leading blanks
c If this is the case jpos is change to point past the end of str2 in str1 
      implicit none
      character(len=*) :: str1,str2
      integer jpos
      integer l1,l2,k
      l1=len(str1)
      l2=len(str2)
      do k=jpos,l1
         select case(str1(k:k))
         case(' ',c_new_line,c_horizontal_tab)
            continue
         case default
            exit
         end select
      enddo
      if( l1-k+1 < l2) then
         start_equal_strings = .false.
      else
         start_equal_strings = str1(k:k+l2-1)==str2
         if(start_equal_strings) jpos = k+l2
      endif
      end
      


      logical function rwl_keypresent(count,key,val)
      character(len=*) key
      integer count
      real * 8 val
      include 'pwhg_rwl.h'
      integer k,l,i,j
      l=rwl_weights_array(count)%num_keys
      do k=1,l
         i=rwl_weights_array(count)%keys(1,k)
         j=rwl_weights_array(count)%keys(2,k)
         if(key == rwl_weights_array(count)%desc(i:j)) then
            val = rwl_weights_array(count)%values(k)
            rwl_keypresent = .true.
            return
         endif
      enddo
      rwl_keypresent = .false.
      end

      subroutine rwl_allocate_weights
      implicit none
      real * 8, allocatable :: tmp(:)
      integer l
      include 'pwhg_rwl.h'
      if(associated(rwl_weights)) then
         l=size(rwl_weights)
         if(l < rwl_num_weights) then
            allocate(tmp(l))
            tmp = rwl_weights
            deallocate(rwl_weights)
            allocate(rwl_weights(rwl_num_weights))
            rwl_weights(1:l) = tmp
            deallocate(tmp)
         endif
      else
         allocate(rwl_weights(rwl_num_weights))
      endif
      end
      
