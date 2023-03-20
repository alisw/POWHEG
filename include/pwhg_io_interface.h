c     -*- Fortran -*-
      interface
         function gzOpen(path,mode) bind(C, NAME='gzopen')
         use, intrinsic :: ISO_C_BINDING
         character(c_char) path,mode
         type(c_ptr) :: gzOpen
         end function
      end interface
      interface
         function gzClose(file) bind(C, NAME='gzclose')
         use, intrinsic :: ISO_C_BINDING
         integer(c_int) :: gzClose
         type(c_ptr), value :: file
         end function
      end interface
      interface
         function gzGets(file,buf,len) bind(C, NAME='gzgets')
         use, intrinsic :: ISO_C_BINDING
         character(c_char) buf(*)
         type(c_ptr) :: gzGets
         type(c_ptr), value :: file
         integer(c_int),value :: len
         end function
      end interface
      interface
         function gzRead(file,buf,len) bind(C, NAME='gzread')
         use, intrinsic :: ISO_C_BINDING
         character(c_char) buf(*)
         integer :: gzRead
         type(c_ptr), value :: file
         integer(c_int),value :: len
         end function
      end interface
      interface
         function gzSeek(file,offset,whence) bind(C, NAME='gzseek')
         use, intrinsic :: ISO_C_BINDING
         integer, parameter :: z_off_t=8
         integer(z_off_t), value :: offset
         integer(z_off_t) gzSeek
         type(c_ptr), value :: file
         integer(c_int),value :: whence
         end function
      end interface
      interface
         function gzRewind(file) bind(C, NAME='gzrewind')
         use, intrinsic :: ISO_C_BINDING
         type(c_ptr), value :: file
         integer(c_int) :: gzRewind
         end function
      end interface
      interface
         function gzWrite(file,buf,len) bind(C, NAME='gzwrite')
         use, intrinsic :: ISO_C_BINDING
         character(c_char) buf(*)
         integer(c_int) :: gzWrite
         type(c_ptr), value :: file
         integer(c_int),value :: len
         end function
c     Writes ch, converted to an unsigned char, into the compressed file.
c     gzputc returns the value that was written, or –1 in case of error.
         function gzPutc(file,ch) bind(C, NAME='gzputc')
         use, intrinsic :: ISO_C_BINDING
         character(c_char), value :: ch
         integer(c_int) :: gzPutc
         type(c_ptr), value :: file
         end function
      end interface
      integer, parameter :: SEEK_CUR=1
      integer, parameter :: z_off_t=8
      integer, parameter :: z_ok=0
      integer, parameter:: pwhg_io_maxfiles=10
      integer, parameter:: szchunk=100
      type(c_ptr) :: files_handle(pwhg_io_maxfiles)=c_null_ptr

      type pwhg_io_buffer_type
      sequence
      character(len=:), pointer :: buffer
c     length: length of buffer
c     current: current position in buffer (next read starts at next byte)
c     upper: position of last byte already red
c     current is either 0 or equal to upper
      integer :: length,current,upper
      character(len=4) :: mode
      logical :: opened
      integer :: fortran_unit
      end type pwhg_io_buffer_type
      integer pwhg_io_initialized
      integer, parameter :: pwhg_io_weirdnum=76391261
      
      type(pwhg_io_buffer_type) :: pwhg_io_buffer(pwhg_io_maxfiles)
      common/pwhg_io_interface/files_handle,pwhg_io_buffer,
     1     pwhg_io_initialized

