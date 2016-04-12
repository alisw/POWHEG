c -*-Fortran-*-
      integer nintervals
      parameter (nintervals=50)
      real * 8 xgrid(0:nintervals,ndiminteg),
     1     xgrid0(0:nintervals,ndiminteg),ymax(nintervals,ndiminteg),
     2     ymaxrat(nintervals,ndiminteg),xacc(0:nintervals,ndiminteg),
     3     xmmm(0:nintervals,ndiminteg),
     4     xint,gen_sigma,gen_sigma2,gen_totev,
     5     xgridrm(0:nintervals,ndiminteg),
     6     xgrid0rm(0:nintervals,ndiminteg),
     7     ymaxrm(nintervals,ndiminteg),
     8     ymaxratrm(nintervals,ndiminteg),
     9     xaccrm(0:nintervals,ndiminteg),
     $     xmmmrm(0:nintervals,ndiminteg),
     $     xintrm,gen_sigmarm,gen_sigma2rm,gen_totevrm
      integer nhits(1:nintervals,ndiminteg),gen_isigma,gen_mcalls,
     1     nhitsrm(1:nintervals,ndiminteg), gen_isigmarm,gen_mcallsrm
      integer ifold(ndiminteg),ifoldrm(ndiminteg)
      common/cgengrids/xgrid,xgrid0,ymax,ymaxrat,xacc,xmmm,
     1     xint,gen_sigma,gen_sigma2,gen_totev,
     1     xgridrm,xgrid0rm,ymaxrm,ymaxratrm,xaccrm,xmmmrm,
     1     xintrm,gen_sigmarm,gen_sigma2rm,gen_totevrm,
     1     ifold,nhits,gen_isigma,gen_mcalls,
     1     ifoldrm,nhitsrm,gen_isigmarm,gen_mcallsrm
      save /cgengrids/
