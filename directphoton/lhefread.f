c...lhefheader(nlf)
c...reads initialization information from a les houches events file on unit nlf. 
      subroutine lhefreadhdr(nlf)
      implicit none
      integer nlf
      character * 200 string
      integer ipr,iret,nch
      include 'LesHouches.h'
      logical ini
      data ini/.true./
      save ini
 1    continue
c     read(nlf,fmt='(a)',err=998,end=998) string
      call pwhg_io_read(nlf,string,iret)
      if(iret == -1) goto 998
      if(adjustl(string(1:10)).eq.'<initrwgt>') then
c     This subroutine is only called by the shower-analysis routines.
c     Here we abandon the old initrwgt handling, and only support the new
c     one. The old handling is only supported by the pwhgreweight.f routines.
         call pwhg_io_backspace(nlf)
         call rwl_readheader(nlf)
         goto 1
      endif
      if(string(1:6).eq.'<init>') then
         call pwhg_io_read(nlf,string,iret)
         if(iret == -1) goto 998
         read(string,*) idbmup(1),idbmup(2),ebmup(1),ebmup(2),
     &        pdfgup(1),pdfgup(2),pdfsup(1),pdfsup(2),idwtup,nprup
         do ipr=1,nprup
            call pwhg_io_read(nlf,string,iret)
            if(iret == -1) goto 998
            read(string,*) xsecup(ipr),xerrup(ipr),xmaxup(ipr),
     &           lprup(ipr)
         enddo
         goto 999
      else
         goto 1
      endif
 998  write(*,*) 'lhefreadhdr: could not find <init> data'
      call exit(1)
 999  end


c...reads event information from a les houches events file on unit nlf. 
      subroutine lhefreadev(nlf)
      implicit none
      integer nlf
      character * 200 string
      include 'LesHouches.h'
      integer i,j,iret
 1    continue
      string=' '
      call pwhg_io_read(nlf,string,iret)
      if(iret /=0 ) goto 666
c      read(nlf,fmt='(a)',err=777,end=666) string
      if(string.eq.'</LesHouchesEvents>') then
         goto 998
      endif
      if(string(1:6).eq.'<event') then
c on error try next event. The error may be caused by merging
c truncated event files. Thus we are tolerant about it.
c Only on EOF return with no event found
         call pwhg_io_read(nlf,string,iret)
         if(iret /=0 ) goto 998
         read(string,*,err=1)nup,idprup,xwgtup,scalup,aqedup,aqcdup
         do i=1,nup
            call pwhg_io_read(nlf,string,iret)
            if(iret /=0 ) goto 998
            read(string,*,err=1) idup(i),istup(i),mothup(1,i),
     &           mothup(2,i),icolup(1,i),icolup(2,i),(pup(j,i),j=1,5),
     &           vtimup(i),spinup(i)
         enddo
         call lhefreadextra(nlf,iret)
         if(iret.ne.0) goto 1
         goto 999
      else
         goto 1
      endif
c no event found:
 777  continue
      print *,"Error in reading"
      print *,string
      nup=0
      return
 666  continue
      print *,"reached EOF"
      print *,string
      nup=0
      return
 998  continue
      print *,"read </LesHouchesEvents>"
      nup=0      
 999  end


      subroutine lhefreadextra(nlf,iret)
      implicit none
      include 'LesHouches.h'
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_rad.h'
      include 'pwhg_st.h'
      include 'pwhg_kn.h'
      include 'pwhg_flg.h'
      include 'pwhg_weights.h'
      include 'pwhg_lhrwgt.h'
      include 'pwhg_rwl.h'
      include 'pwhg_rwgtsuda.h'
      character * 200 string
      integer nlf,iret
      integer iid,iendid
      real * 8 value
      logical lhweights
      logical start_equal_strings
      integer i,j,jpos
      iret = 0
      weights_num = 0
 1    continue
      string=' '
      call pwhg_io_read(nlf,string,iret)
      if(iret /= 0) goto 998
      string = adjustl(string)
      if(string(1:5).eq.'#rwgt') then
         read(string(6:),*) rwl_type,
     $        rwl_index,rwl_weight,rwl_seed,rwl_n1,rwl_n2
      endif
      if(string(1:6).eq.'<rwgt>' .or. string(1:9).eq.'<weights>') then
c     this routine is reached only if flg_rwl_add is true, from the main program,
c     or from the analysis routines. Thus we only enforce the new weight info
c     apparatus. The old apparatus is used only in lhefreadevlhrwgt.
         call pwhg_io_backspace(nlf)
         call rwl_loadweights(nlf,iret)
c on a return with iret != 0 the event will be skipped 
         if(iret.ne.0) return
         goto 1
      endif
      if(string(1:7).eq.'<scales') then
         jpos=9
         if(start_equal_strings(string,'uborns=',jpos)) then
           call getquotedstringpos(string(jpos:),i,j)
           read(string(jpos+i:jpos+j-1),*) uborns
         endif
      endif 
      if(string.eq.'</event>') then
         call pwhg_io_backspace(nlf)
         return
      endif
c Look for old new weight format:
      if(string(1:11).eq.'#new weight') then
         if(weights_num.eq.weights_max) then
            write(*,*) ' too many weights!'
            write(*,*) ' increase weights_max'
            call exit(-1)
         endif
         weights_num = weights_num + 1
         read(string(38:),*) weights_val(weights_num),
     1                       weights_renfac(weights_num),
     2                       weights_facfac(weights_num),
     3                       weights_npdf1(weights_num),
     4                       weights_npdf2(weights_num),
     5                       weights_whichpdf(weights_num)
      endif
      if(string.eq.'# Start extra-info-previous-event') then
         call pwhg_io_read(nlf,string,iret)
         if(iret /= 0) goto 800
         read(string(3:),*) rad_kinreg
         call pwhg_io_read(nlf,string,iret)
 800     if(iret /= 0) then
            write(*,*)
     1           'lhefreadextra:'
            write(*,*)
     1       'found no lines after Start extra-info-previous-event'
            write(*,*) ' exiting ...'
            call exit(-1)
         endif
         read(string(3:),*) rad_type
      endif
      goto 1
      return
 998  continue
      end



c program to handle lh event files with reweight information
c according to the LH standard

      subroutine testlhrwgtstuff
      implicit none
      character * 100 token
      character * 100 values(3)
      integer j,k
      include 'pwhg_lhrwgt.h'

      open(unit=10,file='test.lhe',status='unknown')
      call lhrwgt_readheader(10)
      call lhrwgt_clearheader
      call lhrwgt_headerloader

      write(*,*) lhrwgt_ngroups
      do j=1,lhrwgt_ngroups
         write(*,*)' group',j
         write(*,*)' name=',trim(lhrwgt_group_name_arr(j))
         write(*,*)' combine=',trim(lhrwgt_group_combine_arr(j))
      enddo
      do k=1,lhrwgt_nids
         write(*,*) 'id=',trim(lhrwgt_id_arr(k))
         write(*,*) 'descr=',trim(lhrwgt_descr_arr(k))
         write(*,*) 'group:',
     1        trim(lhrwgt_group_name_arr(lhrwgt_group_ptr(k)))
      enddo

c$$$
c$$$      call next_header_token(token)
c$$$      write(*,'(a)') trim(token)
c$$$      call next_header_token(token)
c$$$      write(*,'(a)') trim(token)
c$$$      call get_values_from_words('combine name',values)
c$$$      write(*,*)'combine:',values(1)
c$$$      write(*,*)'name:',values(2)
c$$$      call next_header_token(token)
c$$$      write(*,'(a)') trim(token)
c$$$      call next_header_token(token)
c$$$      write(*,'(a)') trim(token)
c$$$      call get_values_from_words('id',values)
c$$$      write(*,*)'id:',values(1)
c$$$c      call get_string_up_to('</weight>',string)
c$$$c      write(*,*) string

      end

      subroutine lhrwgt_readheader(iun)
      implicit none
      integer iun
      include 'pwhg_lhrwgt.h'
      character * 200 string
      integer j
      rewind(iun)
      lhrwgt_nheader = 0
 1    read(iun,'(a)') string
      j = index(string,'<initrwgt>')
      if(j.eq.0) then
         if(adjustl(string).eq.'<event>') then
            write(*,*) ' could not find LH reweight info'
            call exit(-1)
         else
            goto 1
         endif
      endif
      string = string(j+len('<initrwgt>'):)
 2    continue
      j = index(string,'</initrwgt>')
      if(j.ne.0) then
         string = string(1:j-1)
         lhrwgt_nheader = lhrwgt_nheader + 1
         lhrwgt_header(lhrwgt_nheader) = string
         return
      else
         lhrwgt_nheader = lhrwgt_nheader + 1
         lhrwgt_header(lhrwgt_nheader) = string         
         read(iun,'(a)') string
         goto 2
      endif
      end


c$$$      subroutine lhrwgt_writeheader(iun)
c$$$      implicit none
c$$$      integer iun
c$$$      include 'pwhg_lhrwgt.h'
c$$$      integer j
c$$$      do j=1,lhrwgt_nheader
c$$$         write(iun,'(a)') trim(lhrwgt_header(j))
c$$$      enddo
c$$$      end
c$$$

      subroutine lhrwgt_clearheader
      implicit none
c should read in the weight information in the Les Houches file header
c should add the group with the NNLOPS reweighting
      include 'pwhg_lhrwgt.h'
c Go through the header, fill the 
c group with no name must be there
      lhrwgt_ngroups = 1
      lhrwgt_nids = 0
      lhrwgt_group_ptr = 0
      lhrwgt_group_name_arr = ' '
      lhrwgt_group_combine_arr = ' '
      lhrwgt_id_arr = ' '
      end

      subroutine lhrwgt_headerloader
      implicit none
c should read in the weight information in the Les Houches file header
c should add the group with the NNLOPS reweighting
      include 'pwhg_lhrwgt.h'
c Go through the header, fill the 
      integer jgroup
      character * 100 values(2),token,string
      logical ini_next_header_token
      common/c_ini_next_header_token/ini_next_header_token
c group with no name must be there
      jgroup = 1
c initialize next_header_token
      ini_next_header_token = .true.
 1    continue
      call next_header_token(token)
      if(token.eq.' ') then
         return
      endif
      if(token(1:1).ne.'<') then
         goto 998
      endif
      call next_header_token(token)
      if(token(1:1).eq.'/') then
         if(token.ne.'/weightgroup') goto 998
         jgroup = 1
         call next_header_token(token)
         if(token.ne.'>') goto 998
         goto 1
      endif
      if(token.eq.'weightgroup') then
         lhrwgt_ngroups = lhrwgt_ngroups + 1
         if(lhrwgt_ngroups.gt.maxgroups) then
            write(*,*)
     1    ' lhrwgt_headerloader: number of groups exceeds maximum'
            write(*,*)
     1    ' increase maxgroups in POWHEG-BOX-V2/include/pwhg_lhrwgt.h'
            write(*,*) 'exiting ... '
            call exit(-1)
         endif
         jgroup = lhrwgt_ngroups
         call get_values_from_words('name combine',values)
         if(values(1).eq.' ') then
            write(*,*) ' lhrwgt_headerloader: group must have a name'
            call exit(-1)
         endif
         lhrwgt_group_name_arr(jgroup)=values(1)
         lhrwgt_group_combine_arr(jgroup)=values(2)
         goto 1
      endif
      if(token.eq.'weight') then
         lhrwgt_nids = lhrwgt_nids + 1
         if(lhrwgt_nids.gt.maxids) then
            write(*,*)
     1    ' lhrwgt_headerloader: number of weights exceeds maximum'
            write(*,*)
     1    ' increase maxids in POWHEG-BOX-V2/include/pwhg_lhrwgt.h'
            write(*,*) 'exiting ... '
            call exit(-1)
         endif
         lhrwgt_group_ptr(lhrwgt_nids) = jgroup
         call get_values_from_words('id',values)
         if(values(1).eq.' ') then
            write(*,*) ' lhrwgt_headerloader: weight must have an id'
            call exit(-1)
         endif
         lhrwgt_id_arr(lhrwgt_nids) = values(1)
c get string up to </weight>
         call get_string_up_to('</weight>',
     1        lhrwgt_descr_arr(lhrwgt_nids))
         goto 1
      endif
      return
 998  write(*,*) '  lhrwgt_headerloader:'
      write(*,*) ' inconsistent header'
      call exit(-1)
      end




     
      subroutine get_values_from_words(words,values)
      implicit none
      character *(*) words, values(*)
      integer maxwords
      parameter (maxwords=10)
      integer ii(2,maxwords+1),k
      character * 100 token
      call break_string_in_words(words,maxwords,ii)
      do k=1,1000
         if(ii(1,k).eq.0) exit
         values(k)=' '
      enddo
 1    continue
      call next_header_token(token)
      if(token(1:1).eq.'>') then
         return
      endif
      do k=1,1000000
         if(ii(1,k).eq.0) exit
         if(token.eq.words(ii(1,k):ii(2,k))) then
            call next_header_token(token)
            if(token.ne.'=') then
               write(*,*) 'get_values_from_words: missing = after '
     1              //words(ii(1,k):ii(2,k))
               call exit(-1)
            endif
            if(values(k).ne.' ') then
               write(*,*) ' get_values_from_words: more than one value'
               write(*,*) ' for ',words(ii(1,k):ii(2,k))
               call exit(-1)
            endif
            call next_header_token(values(k))
            goto 1
         endif
      enddo
      end

      subroutine break_string_in_words(string,maxwords,ii)
c It returns in ii(2,maxwords) the index of the first and last letter
c of each word in string. If there are n words in string, ii(:,n+1) is
c set to zero. Thus
c where outarr is an array of strings. If name="aaaa" surname="bbb"
c dateofbirth="18/03/1952" appear in the following input, in any order,
c it returns outarr(1)="aaaa", outarr(2)="bbb", outarr(3)="18/03/1952".
c If (for example) surname is missing, then outarr(2)=" " (blank).
      implicit none
      character *(*) string
      integer maxwords
      integer ii(2,maxwords+1)
      integer nnames,k,l
c ii(1/2,*) are the first and last caracter of each work in string,
c in order. The last value is zero
      k=1
      l=len(string)
      do nnames=1,1000000
         do while(string(k:k).eq.' ')
            if(k.eq.l) then
               ii(1:2,nnames)=0
               return
            else
               k=k+1
            endif
         enddo
C The last nnames must have ii=0
         if(nnames.eq.maxwords+1) then
            write(*,*) 'get_values: too many names'
            call exit(-1)
         endif
         ii(1,nnames)=k
         do while(string(k:k).ne.' ')
            if(k.eq.l) then
               ii(2,nnames)=l
               ii(:,nnames+1)=0
               return
            else
               k=k+1
            endif
         enddo
         ii(2,nnames)=k-1
      enddo
      end


      subroutine push_back_header_token(token)
      implicit none
      character *(*) token
      include 'pwhg_lhrwgt.h'
      integer jheader
      character * 200 stringtoparse
      common/lhheaderparser/stringtoparse,jheader
      stringtoparse = trim(stringtoparse)//' '//trim(adjustl(token))
      if(stringtoparse.ne.trim(stringtoparse)//' '
     1     //trim(adjustl(token))) then
         write(*,*)' push_back_header_token: no room to push back'
         call exit(-1)
      endif
      end

      subroutine next_header_token(token)
      implicit none
      character *(*) token
      include 'pwhg_lhrwgt.h'
      integer jheader
      character * (lhrwgt_max_header_columns) stringtoparse
      common/lhheaderparser/stringtoparse,jheader
      logical ini_next_header_token
      common/c_ini_next_header_token/ini_next_header_token
      if(ini_next_header_token) then
         jheader = 1
         stringtoparse = lhrwgt_header(jheader)
         ini_next_header_token=.false.
      endif
 1    continue
      call next_token_from_string(stringtoparse,token)
      if(token.eq.' ') then
         if(jheader .lt. lhrwgt_nheader) then
            jheader = jheader + 1
            stringtoparse = lhrwgt_header(jheader)
            goto 1
         else
c nothing else to parse; return the empty token
            return
         endif
      endif
      end


      subroutine get_string_up_to(tag,string)
      implicit none
      character *(*) string,tag
      include 'pwhg_lhrwgt.h'
      integer jheader
      character * 200 stringtoparse
      common/lhheaderparser/stringtoparse,jheader
      integer lll,left,j
      logical more
      string = ' '
 1    j = index(stringtoparse,tag)
c room left in string
      left = len(string) - len(trim(string))
      if(j.ne.0) then
c length of stuff to add
         lll = len(trim(stringtoparse(1:j-1)))
         more = .false.
      else
         lll = len(trim(stringtoparse))
         more = .true.
      endif
      if(lll.le.left-1) then
         string = trim(string)//' '//stringtoparse(1:lll)
         if(more) then
            jheader = jheader+1
            stringtoparse = lhrwgt_header(jheader)
            goto 1
         else
            stringtoparse = stringtoparse(j+len(tag):)
         endif
      else
         write(*,*) ' get_string_up_to: string too short'
         call exit(-1)
      endif
      end


      subroutine next_token_from_string(string,token)
c Returns the next token in string, taking it away from string
      implicit none
      character * (*) string, token
      character * 1 sep
      integer iend,itmp
      if(string.eq.' ') then
         token=' '
         return
      endif
c     get rid of leading spaces
      string=adjustl(string)
      select case (string(1:1))
      case ('=','>','<')
         token=string(1:1)
         string = string(2:)
         return
      case ("'",'"')
         sep = string(1:1)
         iend = index(string(2:),sep)+1
         if(iend.eq.1) then
            write(*,*) ' next_token: unterminated string'
            call exit(-1)
         endif
         token = string(1:iend)
         if(token.ne.string(1:iend)) then
            write(*,*) ' next_token: token string too short'
            call exit(-1)
         endif
         string=string(iend+1:)
         return
      case default
         iend = 1000000
         itmp = index(string,' ')
         if(itmp.gt.0) iend = min(itmp,iend)
         itmp = index(string,'=')
         if(itmp.gt.0) iend = min(itmp,iend)
         itmp = index(string,'<')
         if(itmp.gt.0) iend = min(itmp,iend)
         itmp = index(string,'>')
         if(itmp.gt.0) iend = min(itmp,iend)
         iend = iend - 1
         if(iend.eq.1000000-1) then
            iend = len(string)
         endif
         token=string(1:iend)
         if(token.ne.string(1:iend)) then
            write(*,*) ' next_token: token string too short'
            call exit(-1)
         endif
         string=adjustl(string(iend+1:))
      end select
      end


      subroutine lhrwgt_loadweights(iun,iret)
      implicit none
      integer iun,iret
c assumes that the '<rwgt>' line has already been red
      include 'pwhg_lhrwgt.h'
      integer nweights,iid,iendid
      real * 8 value
      character * 100 string
      nweights = 0
 1    read(unit=iun,fmt='(a)',end=998) string
      call lhrwgt_id_value_ind(string,iid,iendid,value)
      if(iid.lt.0) goto 996
      nweights = nweights+1
      if(nweights.gt.lhrwgt_nids) then
         write(*,*) ' lhrwgt_loadweights: more weights in event'
         write(*,*) ' than in the header'
         write(*,*) ' exiting ...'
         call exit(-1)
      endif
      if(string(iid:iendid).eq.lhrwgt_id_arr(nweights)) then
         lhrwgt_weights(nweights) = value
      else
         write(*,*) ' lhrwgt_loadweights: id in event'
         write(*,*) ' does not match id in declaration'
         write(*,*) ' exiting ...'
         call exit(-1)
      endif
      goto 1
 996  continue
c did we find </rwgt>?
      if(adjustl(string).eq.'</rwgt>') then
         if(nweights.eq.lhrwgt_nids) then
c     found all weights;
            iret = 0
         else
            write(*,*) ' lhrwgt_loadweights: did not find all weights!'
            write(*,*) ' exiting ...'
            call exit(-1)
         endif
      else
         iret = -1
      endif
      return
 998  continue
c end of file
      iret = -1
      end


      subroutine lhrwgt_id_value_ind(string,iid,iendid,value)
      implicit none
      character *(*) string
      integer iid,iendid
      real * 8 value
      integer j,k,l,ios
      j=1
      l=len(string)
      do while(string(j:j).eq.' '.and.j.lt.l)
         j = j + 1
      enddo
      if(string(j:j+4).ne.'<wgt ') goto 777
      j=j+5
      do while(string(j:j).eq.' '.and.j.lt.l)
         j = j + 1
      enddo
      if(string(j:j+1).ne.'id') goto 777
      j=j+2
      do while(string(j:j).eq.' '.and.j.lt.l)
         j = j + 1
      enddo
      if(string(j:j).ne.'=') goto 777
      j = j + 1
      do while(string(j:j).eq.' '.and.j.lt.l)
         j = j + 1
      enddo
      select case (string(j:j))
      case ('"',"'")
         iid = j
         iendid = index(string(iid+1:),string(iid:iid))+iid
         j=iendid+1
      case default
         goto 777
      end select
      do while(string(j:j).eq.' '.and.j.lt.l)
         j = j + 1
      enddo
      if(string(j:j).ne.'>') goto 777
      j = j+1
      do while(string(j:j).eq.' '.and.j.lt.l)
         j = j + 1
      enddo
      k = index(string,'</wgt>')
      if(k.eq.0) goto 777
      read(unit=string(j:k-1),fmt=*,iostat=ios) value
      if(ios.ne.0) goto 777
      return
c Something went wrong:
 777  iid = -1
      end

      subroutine copy_lhrw_to_weights(num)
      implicit none
      include 'pwhg_weights.h'
      include 'pwhg_lhrwgt.h'
      integer num
      if(num.gt.lhrwgt_nids) then
         write(*,*) ' copy_lhrw_to_weights: error, num>lhrwgt_nids'
         call exit(-1)
      endif
      weights_num = num
      if(num .gt. weights_max) then
         write(*,*) ' copy_lhrw_to_weights: too many weights,'
         write(*,*) ' increase weights_max in pwhg_weights.h'
         write(*,*) ' exiting ...'
      endif
      weights_val(:num) = lhrwgt_weights(:num)
      end
