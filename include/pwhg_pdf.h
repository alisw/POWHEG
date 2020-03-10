c -*- Fortran -*-
      integer pdf_ih1,pdf_ih2,pdf_ndns1,pdf_ndns2,pdf_nparton
      integer pdf_ndns1lhe,pdf_ndns2lhe
      logical pdf_dis_photon,pdf_alphas_from_PDF
      real * 8 pdf_q2min,pdf_cutoff_fact
      common/pwhg_pdf/pdf_q2min,pdf_cutoff_fact,
     1     pdf_ih1,pdf_ih2,pdf_ndns1,pdf_ndns2,
     2     pdf_nparton,pdf_dis_photon,pdf_alphas_from_PDF,
     3     pdf_ndns1lhe,pdf_ndns2lhe
      save /pwhg_pdf/
