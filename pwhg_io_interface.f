      subroutine pwhg_io_open_read(path,iun,iret)
      use, intrinsic :: ISO_C_BINDING
      implicit none
      include 'pwhg_io_interface.h'
      character(len=*) :: path
      integer iun,iret
      call pwhg_io_open('rb',path,iun,iret)
      allocate(character(len=szchunk*2)::pwhg_io_buffer(iun)%buffer)
c     current is the pointer to the file position. Reading will start at
c     the next character
      pwhg_io_buffer(iun)%current=0
c     upper is a pointer to the last character already in the buffer
      pwhg_io_buffer(iun)%upper=0
c     we must always have: 0 <= current <= upper.
c     length buffer
      pwhg_io_buffer(iun)%length=szchunk*2
      pwhg_io_buffer(iun)%mode='rb'
      end

      subroutine pwhg_io_rewind(iun)
      use, intrinsic :: ISO_C_BINDING
      implicit none
      include 'pwhg_io_interface.h'
      integer iun,iret
      if(.not. pwhg_io_buffer(iun)%opened) then
         write(*,*) ' pwhg_io_rewind: unit ',iun,
     1        ' is not opened, exiting ...'
         call exit(-1)
      endif
      iret = gzrewind(files_handle(iun))
      if(pwhg_io_buffer(iun)%mode == 'rb') then
c     empty buffer
         pwhg_io_buffer(iun)%current = 0
         pwhg_io_buffer(iun)%upper = 0
      endif
      end

      subroutine pwhg_io_open_write(path,iun,compress,iret)
      use, intrinsic :: ISO_C_BINDING
      implicit none
      include 'pwhg_io_interface.h'
      character(len=*) :: path
      integer iun,iret
      logical compress
      if(compress) then
         call pwhg_io_open('wb',path,iun,iret)
         pwhg_io_buffer(iun)%mode='wb'
      else
         call pwhg_io_open('wT',path,iun,iret)
         pwhg_io_buffer(iun)%mode='wT'
      endif
      end

      subroutine pwhg_io_open(mode,path,iun,iret)
      use, intrinsic :: ISO_C_BINDING
      implicit none
      include 'pwhg_io_interface.h'
      character(*) path,mode
      integer iun,iret
      integer j
      logical ok
      call init
      iun=0
      do j=1,pwhg_io_maxfiles
c unit 6 is always opened if init was called
         if(.not. pwhg_io_buffer(j)%opened ) then
            iun=j
            exit
         endif
      enddo
      if(iun==0) then
         write(*,*)
     1        'pwhg_io_open : too many opened files'
         write(*,*) 'increase pwhg_io_maxfiles'
         write(*,*) 'exiting ...'
         call exit(-1)
      endif
      if( mode == 'wT') then
c     'wT' yields to uncompressed files. On some zlib implementation
c     this does not work. Use standard fortran style io
         do j = 7,99
            inquire(unit=j,opened=ok)
            if(.not.ok) then
               exit
            endif
         enddo
         pwhg_io_buffer(iun)%fortran_unit=j
         open(file=trim(path),unit=j,status='unknown')
      else
         pwhg_io_buffer(iun)%fortran_unit=0
         files_handle(iun) = gzOpen(trim(path)//c_null_char,
     1        trim(mode)//c_null_char)
         if(.not. c_associated(files_handle(iun))) then
            write(*,*) ' pwhg_io_open: cannot open ',trim(path)
            iret = -1
            return
         endif
      endif
      pwhg_io_buffer(iun)%opened = .true.
      iret = 0
      contains
      subroutine init
      if(pwhg_io_initialized /= pwhg_io_weirdnum) then
         pwhg_io_initialized = pwhg_io_weirdnum
         pwhg_io_buffer(:)%opened = .false.
c Unit 6 is always open (standard output)
         pwhg_io_buffer(6)%opened = .true.
      endif
      end subroutine init
      end

      subroutine pwhg_io_close(iun)
      use, intrinsic :: ISO_C_BINDING
      implicit none
      include 'pwhg_io_interface.h'
      integer iun,iret
      if(.not. pwhg_io_buffer(iun)%opened) then
         write(*,*) ' pwhg_io_close: unit ',iun,
     1        ' is not opened, exiting ...'
         call exit(-1)
      endif
      if(pwhg_io_buffer(iun)%fortran_unit>0) then
         close(pwhg_io_buffer(iun)%fortran_unit)
         pwhg_io_buffer(iun)%fortran_unit = 0
      else
         iret = gzClose(files_handle(iun))
         if(iret /= 0) then
            write(*,*)
     1           ' pwhg_io_close: some error while closing unit ',iun
            write(*,*) ' iret=',iret
         endif
         files_handle(iun) = c_null_ptr
         if(pwhg_io_buffer(iun)%buffer == 'rb') then
            deallocate(pwhg_io_buffer(iun)%buffer)
         endif
      endif
      pwhg_io_buffer(iun)%opened = .false.
      end

      subroutine pwhg_io_read_buf(iun,buffer)
c     reads len(buffer) characters from unit iun. No error checking (must make sure
c     that the read operation can be performed
      use, intrinsic :: ISO_C_BINDING
      implicit none
      include 'pwhg_io_interface.h'
      integer iun,iret
      character(len=*) :: buffer
      if(.not. pwhg_io_buffer(iun)%opened) then
         write(*,*) ' pwhg_io_read_buf: unit ',iun,
     1        ' is not opened, exiting ...'
         call exit(-1)
      endif
      buffer=pwhg_io_buffer(iun)%buffer(pwhg_io_buffer(iun)%current+1:
     1     pwhg_io_buffer(iun)%current + len(buffer))
      pwhg_io_buffer(iun)%current = pwhg_io_buffer(iun)%current +
     1     len(buffer)
      end

      subroutine pwhg_io_read(iun,buffer,iret)
      use, intrinsic :: ISO_C_BINDING
      implicit none
      include 'pwhg_io_interface.h'
      integer iun,iret
      character(len=*) :: buffer
      integer nchar
      if(.not. pwhg_io_buffer(iun)%opened) then
         write(*,*) ' pwhg_io_read: unit ',iun,
     1        ' is not opened, exiting ...'
         call exit(-1)
      endif
      call pwhg_io_skip_until(iun,c_new_line,nchar)
      if(nchar < 0) then
         if(pwhg_io_buffer(iun)%current == 0) then
            iret=-1
         else
            buffer =
     1   pwhg_io_buffer(iun)%buffer(1:pwhg_io_buffer(iun)%current-2)
            iret = 0
         endif
      else
         buffer=pwhg_io_buffer(iun)%buffer(1:nchar-1)
         iret = 0
      endif
      end


      subroutine pwhg_io_write(iun,buffer)
      use, intrinsic :: ISO_C_BINDING
      implicit none
      include 'pwhg_io_interface.h'
      integer iun
      character(*) buffer
      integer j,l
      if(.not. pwhg_io_buffer(iun)%opened) then
         if(iun /= 6) then
c     unit 6 may be used before the pwhg_io interface is ever accessed
            write(*,*) ' pwhg_io_write: unit ',iun,
     1           ' is not opened, exiting ...'
            call exit(-1)
         endif
      endif
      l=len(buffer)
      if(iun == 6) then
         write(*,'(a)') trim(buffer)
      elseif(pwhg_io_buffer(iun)%fortran_unit > 0) then
         write(pwhg_io_buffer(iun)%fortran_unit,'(a)') trim(buffer)
      else
         j = gzwrite(files_handle(iun),buffer,l)
         if(j /= l ) then
            write(*,*) ' pwhg_io_write: error writing'
            write(*,*) ' exiting ...'
            call exit(-1)
         endif
         j=gzPutc(files_handle(iun),c_new_line)
         if(j == -1 ) then
            write(*,*) ' pwhg_io_write: error writing new line'
            write(*,*) ' exiting ...'
            call exit(-1)         
         endif
      endif
      end

      subroutine pwhg_io_skip(iun,nskip)
c     moves file position by nskip bytes
      use, intrinsic :: ISO_C_BINDING
      implicit none
      include 'pwhg_io_interface.h'
      integer iun,nskip
      integer k
      if(.not. pwhg_io_buffer(iun)%opened) then
         write(*,*) ' pwhg_io_skip: unit ',iun,
     1        ' is not opened, exiting ...'
         call exit(-1)
      endif
      k = pwhg_io_buffer(iun)%current + nskip
      if(k<0 .or. k>pwhg_io_buffer(iun)%upper) then
         write(*,*) ' pwhg_io_skip: cannot skip outside of buffer'
         write(*,*) ' exiting ...'
         call exit(-1)
      endif
      pwhg_io_buffer(iun)%current = k
      end
      
      subroutine pwhg_io_backspace(iun)
c     moves the file position back before the previous operation
      use, intrinsic :: ISO_C_BINDING
      implicit none
      include 'pwhg_io_interface.h'
      integer iun,nskip
      integer(z_off_t) offset,longskip
      if(.not. pwhg_io_buffer(iun)%opened) then
         write(*,*) ' pwhg_io_backspace: unit ',iun,
     1        ' is not opened, exiting ...'
         call exit(-1)
      endif
      if(pwhg_io_buffer(iun)%current == 0) then
         write(*,*) ' pwhg_io_backspace: cannot backspace!'
         write(*,*) ' curremt:',pwhg_io_buffer(iun)%current
         write(*,*) ' exiting ...'
         call exit(-1)
      endif
      pwhg_io_buffer(iun)%current=0
      end
      
      subroutine pwhg_io_skip_until(iun,mark,ntot)
c     Search for the string <mark> in the file, and moves the file position
c     at the first character after mark. Returns in ntot the number of bytes
c     up to the last character in mark. If the mark is not found,
c     typically because an end of file is found before the mark, ntot is set to -1,
c     the buffer contains all characters that was able to get, and the current pointer
c     is set to right after the last valid character in the buffer
      use, intrinsic :: ISO_C_BINDING
      implicit none
      include 'pwhg_io_interface.h'
      integer iun,ntot
      character(len=*) :: mark
      character(len=szchunk) buffer
      integer upper,inum,imark
      if(.not. pwhg_io_buffer(iun)%opened) then
         write(*,*) ' pwhg_io_skip_until: unit ',iun,
     1        ' is not opened, exiting ...'
         call exit(-1)
      endif
c     discarde saved result from previous read
      if(pwhg_io_buffer(iun)%current > 0) then
         pwhg_io_buffer(iun)%buffer=pwhg_io_buffer(iun)%buffer(
     1        pwhg_io_buffer(iun)%current+1:pwhg_io_buffer(iun)%upper)
         pwhg_io_buffer(iun)%upper=pwhg_io_buffer(iun)%upper -
     1        pwhg_io_buffer(iun)%current
         pwhg_io_buffer(iun)%current=0
      endif
c
      do
         upper = pwhg_io_buffer(iun)%upper
         imark=index(pwhg_io_buffer(iun)%buffer(1:upper),mark)
         if(imark > 0) then
            pwhg_io_buffer(iun)%current = imark + len(mark)-1
            ntot = pwhg_io_buffer(iun)%current
            exit
         else
            if(upper+szchunk >  pwhg_io_buffer(iun)%length) then
               call pwhg_io_enlarge_buffer
            endif
c read one more chunk into buffer. If it doesn't fit, increase buffer size.
            inum = gzRead(files_handle(iun),buffer(1:szchunk),szchunk)
            if(inum > 0) then
               pwhg_io_buffer(iun)%buffer(upper+1:upper + inum) = buffer(1:inum)
               pwhg_io_buffer(iun)%upper = upper + inum
               cycle
            else
               ntot = -1
               exit
            endif
         endif
      enddo
      contains
      subroutine pwhg_io_enlarge_buffer
      character (len=:), pointer :: buf
      integer i1,h1,i2,h2,length,nu,nl,lg,olg,ibuf,cur,icur
      logical, save :: ini=.true.
      real(kind=8) :: powheginput
      integer, save :: max_io_bufsize
      if(ini) then
         max_io_bufsize = powheginput("#max_io_bufsize")
         if(max_io_bufsize < 0) max_io_bufsize=100000
         ini = .false.
      endif
c     Put safety limit on buffer size
      olg=pwhg_io_buffer(iun)%length
      lg = olg + 1000
      if(lg>max_io_bufsize) then
         write(*,*) " pwhg_io_skip_until: we can't find mark '"
     1        //mark//"'"
         write(*,*) " and are above max_io_bufsize limit."
         write(*,*) " if the mark is really there, you must."
         write(*,*) " set max_io_bufsize to a value larger than ",
     1        max_io_bufsize," in powheg.input"
         write(*,*) " exting ..."
         call exit(-1)
      endif
      allocate(character(len=lg)::buf)
      buf=pwhg_io_buffer(iun)%buffer(1:olg)
      deallocate(pwhg_io_buffer(iun)%buffer)
      pwhg_io_buffer(iun)%buffer => buf
      pwhg_io_buffer(iun)%length = lg
      end subroutine pwhg_io_enlarge_buffer
     
      end


