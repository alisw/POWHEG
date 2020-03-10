      subroutine multi_plot_setup(dsig0,dsig,ndim)
      implicit none
      integer ndim
      real * 8 dsig0,dsig(ndim)
      logical, save :: ini = .true.
      include 'pwhg_weights.h'
      include 'pwhg_rwl.h'
      character * 6 whcprg
      common/cwhcprg/whcprg
      integer, save :: nweights = 0
      if(ini) then
         if(whcprg /= 'NLO') then
            if(rwl_initialized == rwl_initialized_const) then
               nweights = rwl_num_weights
               if(nweights + 1 > ndim) then
                  write(*,*) '  multi_plot_setup: there are ',
     1                 rwl_num_weights,' weights, but ndim=',ndim
                  write(*,*) ' cannot handle, exiting ...'
                  call exit(-1)
               endif
               call setupmulti(nweights+1)
            endif
         else
            call setupmulti(1)
         endif
         ini = .false.
      endif

      dsig=0
      dsig(1) = dsig0
      if(nweights.gt.0) then
         dsig(2:nweights+1)=rwl_weights(1:nweights)
      endif

      end
