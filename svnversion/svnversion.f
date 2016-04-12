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
      open(unit=12,file='svnversion.tmp',status='unknown')
      write(12,'(a)') 'c -*- Fortran -*-'
      write(12,'(a)') "      write(iun,'(a)') 'SVN status Begin'"
      do j=1,1000000
         read(unit=11,fmt='(a)',iostat=ios,end=99) string
         lll = len(trim(string))
         write(12,'(a)') "      write(iun,'(a)')"
         do while(lll.gt.lchunks)
            write(12,'(a)') "     1'"//string(1:lchunks)//"'//"
            string = string(lchunks+1:)
            lll = lll-lchunks
         enddo
         write(12,'(a)') "     1'"//trim(string)//"'"
      enddo
 99   continue
      write(12,'(a)') "      write(iun,'(a)') 'SVN status End'"
      end
