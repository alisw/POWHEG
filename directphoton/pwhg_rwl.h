c -*- Fortran -*-

      type rwl_lhe_block
         sequence
         integer nup,idprup,rad_type,rad_kinreg
         real * 8 xwgtup,scalup,aqedup,aqcdup
         integer, pointer :: idup(:),istup(:),mothup(:,:),icolup(:,:)
         real * 8, pointer :: pup(:,:),vtimup(:),spinup(:)
         integer rwl_type,rwl_index,rwl_seed,rwl_n1,rwl_n2
         real * 8 rwl_weight
	 real * 8 sudawgt
	 real * 8 uborns
      end type rwl_lhe_block

      type rwl_weight_info
         sequence
         character(len=:), pointer :: id, desc
         integer :: group, num_keys
         integer, pointer :: keys(:,:)
         real * 8 , pointer :: values(:)
         logical :: filled
      end type rwl_weight_info

      type rwl_group_info
         sequence
         character(len=:), pointer :: name, combine
      end type rwl_group_info

      real * 4, pointer :: rwl_weights(:)

      integer, parameter :: rwl_maxweights=500,rwl_maxgroups=50
      integer :: rwl_num_weights, rwl_num_groups
      type(rwl_weight_info) rwl_weights_array(rwl_maxweights)
      type(rwl_group_info) rwl_groups_array(rwl_maxgroups)
      integer rwl_type,rwl_index,rwl_seed,rwl_n1,rwl_n2
c     if rwl_initialized /= rwl_initialized_const, it means
c     that 
      integer rwl_initialized
      integer, parameter :: rwl_initialized_const = 1729178316
      real * 8  rwl_weight
      logical rwl_format_rwgt
      common/pwhg_rwl_common/rwl_weight,rwl_type,rwl_index,rwl_seed,
     1     rwl_n1,rwl_n2,rwl_num_weights,rwl_weights_array,
     2     rwl_groups_array,rwl_weights,rwl_num_groups,
     3     rwl_format_rwgt,rwl_initialized
      save /pwhg_rwl_common/
