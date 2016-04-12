c -*- Fortran -*-
      type histptr
         sequence
         character * 20 id
         real * 8, pointer :: xhistarr(:),yhistarr(:,:),yhistarr1(:,:),
     1        errhistarr1(:,:),yhistarr2(:,:),errhistarr2(:,:)
         integer, pointer :: nhits(:)
         integer nbins,ient1,nmulti
      end type histptr
      type(histptr), pointer :: hist_ptr(:)
      integer nhist,nmulti,jhist
      common /histnew/hist_ptr,nhist,nmulti,jhist

