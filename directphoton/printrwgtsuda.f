      subroutine printrwgtsuda(nlf,weight)
      implicit none
      integer nlf
      double precision weight
      character(len=300) string
      write(string,'(a,e16.9,a)')"<wgt id='sudakovwgt'> ",
     1     weight,' </wgt>' 
      call pwhg_io_write(nlf,trim(string))
      end
