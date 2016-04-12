c -*- Fortran -*-

c Maximum number of weights
      integer lhrwgt_maxnheader
      parameter (lhrwgt_maxnheader=200)
      integer  lhrwgt_max_header_columns
      parameter ( lhrwgt_max_header_columns=200)
      character * (lhrwgt_max_header_columns)
     1     lhrwgt_header(lhrwgt_maxnheader)
      character * 100 lhrwgt_id,lhrwgt_group_name,
     1     lhrwgt_group_combine
      integer lhrwgt_nheader
      character * 400 lhrwgt_descr

      integer maxgroups
      parameter (maxgroups=20)
      integer maxids
      parameter (maxids=200)
      integer lhrwgt_ngroups,lhrwgt_nids
      character * 100 lhrwgt_group_name_arr(maxgroups),
     1                lhrwgt_group_combine_arr(maxgroups),
     2                lhrwgt_id_arr(maxids)
      character * 400 lhrwgt_descr_arr(maxids)
      real * 8 lhrwgt_weights(maxids)
      integer lhrwgt_group_ptr(maxids)


      common/pwhg_lhrwgt/lhrwgt_nheader,lhrwgt_header,
     1     lhrwgt_id,lhrwgt_descr,lhrwgt_group_name,
     2     lhrwgt_group_combine,
     3     lhrwgt_weights,lhrwgt_ngroups,lhrwgt_nids,
     4     lhrwgt_group_ptr,lhrwgt_group_name_arr,
     5     lhrwgt_group_combine_arr,lhrwgt_id_arr,lhrwgt_descr_arr




      
