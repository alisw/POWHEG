      subroutine multi_plot_setup(dsig0,dsig,ndim)
      implicit none
      integer ndim
      real * 8 dsig0,dsig(ndim)
      logical, save :: ini = .true.
      character * 6 whcprg
      common/cwhcprg/whcprg
      integer, save :: nweights = 0
      integer maxweights, numweights,radtype
      parameter(maxweights=50)
      real *8 weight
      common/weights/weight(maxweights),numweights,radtype
      integer j
      include 'herwigsettings.h'
      integer i
      if(ini) then
         if(whcprg /= 'NLO') then
            if(numweights > ndim + 1) then
               write(*,*) '  multi_plot_setup: there are ',
     $              numweights,' weights, but ndim=',ndim
               write(*,*) ' cannot handle, exiting ...'
               call exit(-1)
            endif
            call setupmulti(numweights+1)
         else
            ini = .false.
         endif
      endif

c     rescale the weight of the event depending on the rad_type
c     (1..btilde, 2..remn)
c     using the ub_..._corr factors
      if (radtype == 1) then
         ub_corr = ub_btilde_corr
      else if (radtype == 2) then
         ub_corr = ub_remn_corr
      else 
         ub_corr = 1d0
      endif
      dsig=0
      dsig(1) = dsig0
      if(numweights > 0) then
         dsig(2:numweights+1)=weight(1:numweights)
      endif

c      write(*,*) ' weights:'
c      do j=1,numweights+1
c         write(*,'(e10.5)') dsig(j)
c      enddo
c      write(*,*) ' end weights:'
      
      dsig=ub_corr * dsig


      end
