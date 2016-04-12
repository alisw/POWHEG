c -*- Fortran -*-
      integer pdf_ih1,pdf_ih2,pdf_ndns1,pdf_ndns2,pdf_nparton
      integer pdf_ndns1lhe,pdf_ndns2lhe
      logical pdf_dis_photon
      real * 8 pdf_q2min
      common/pwhg_pdf/pdf_q2min,pdf_ih1,pdf_ih2,pdf_ndns1,pdf_ndns2,
     1     pdf_nparton,pdf_dis_photon,pdf_ndns1lhe,pdf_ndns2lhe
      save /pwhg_pdf/
