c -*-Fortran-*-
      double precision modlog_p,resscfact
      logical flg_modlog
     $     ,flg_include_delta_terms,flg_distribute_by_ub
     $     ,flg_distribute_by_ub_AP,flg_uubornonly
     $     ,flg_hoppet_initialized,flg_use_NNLOPS_pdfs
     $     ,flg_dtermsallorders
      character * 4 flg_minnloproc
      common/minnlo_flg/modlog_p,resscfact,flg_modlog
     $     ,flg_include_delta_terms
     $     ,flg_distribute_by_ub,flg_distribute_by_ub_AP
     $     ,flg_uubornonly,flg_minnloproc,flg_dtermsallorders
      common/hoppet_pdf_flg/ flg_hoppet_initialized
     $     ,flg_use_NNLOPS_pdfs
