      program versioninclude
      implicit none
      character * 400 string
      integer lchunks
      parameter (lchunks=50)
      integer ios,j,lll
      open(unit=11,file='svnversion.txt',iostat=ios,status='old')
      if(ios.ne.0) then
         write(*,*) ' cannot find svnversion.txt'
         call exit(-1)
      endif
      open(unit=13,file='svnversion.tmp',status='unknown')
      write(13,'(a)') 'c -*- Fortran -*-'
      write(13,'(a)')
     1 "      call pwhg_io_write(iun, 'SVN status Begin')"
      do j=1,1000000
         read(unit=11,fmt='(a)',iostat=ios,end=99) string
         lll = len(trim(string))
         write(13,'(a)') "      call pwhg_io_write(iun,"
         do while(lll.gt.lchunks)
            write(13,'(a)') "     1'"//string(1:lchunks)//"'//"
            string = string(lchunks+1:)
            lll = lll-lchunks
         enddo
         write(13,'(a)') "     1'"//trim(string)//"')"
      enddo
 99   continue
      write(13,'(a)') "      call pwhg_io_write(iun, 'SVN status End')"
      close(11)
      close(13)
      end
