      subroutine multi_plot_setup(dsig0,dsig,ndim)
      implicit none
      integer ndim
      real * 8 dsig0,dsig(ndim)
      logical, save :: ini = .true.
      include 'pwhg_weights.h'
      include 'pwhg_lhrwgt.h'
      character * 6 whcprg
      common/cwhcprg/whcprg

      if(ini) then
         if(whcprg /= 'NLO') then
            if(lhrwgt_nids.gt.0) then
               if(lhrwgt_nids > ndim) then
                  write(*,*) '  multi_plot_setup: there are ',
     1                 lhrwgt_nids,' weights, but ndim=',ndim
                  write(*,*) ' cannot handle, exiting ...'
                  call exit(-1)
               endif
               call setupmulti(lhrwgt_nids)
            elseif(weights_num.gt.0) then
               if(weights_num + 1 > ndim) then
                  write(*,*) '  multi_plot_setup: there are ',
     1                 weights_num + 1,' weights, but ndim=',ndim
                  write(*,*) ' cannot handle, exiting ...'
                  call exit(-1)
               endif
               call setupmulti(weights_num+1)
            else
               if(ndim <= 0) then
                  write(*,*) '  multi_plot_setup: ndim=',ndim
                  write(*,*) ' cannot handle, exiting ...'
                  call exit(-1)
               endif
               call setupmulti(1)
            endif
         else
            call setupmulti(1)
         endif
      endif

      dsig=0
      if(lhrwgt_nids.gt.0) then
         dsig(1:lhrwgt_nids)=lhrwgt_weights(1:lhrwgt_nids)
      elseif(weights_num.gt.0) then
         dsig(1)=dsig0
         dsig(2:weights_num+1)=weights_val(1:weights_num)
      else
         dsig(1)=dsig0
      endif

      end
