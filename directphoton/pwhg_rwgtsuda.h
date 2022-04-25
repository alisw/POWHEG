c     Common block to save the weight of an event generated with a
c     modified Sudakov.
c
c     flg_sudakov_rwgt activates reweighting with modified Sudakov.
c     flg_sudarw_tmp is used to decide if a certain underlying Born
c     allows for a radiation that is condition to the reweighting
c     (i.e. in our case if the uborn is pure QCD)
c
c     sudawgt contains the event weight associated with
c     the modified Sudakov
c     sudarwgtfact contains the constant with which the Sudakov exponent has
c     been scaled
      logical flg_sudakov_rwgt, flg_sudarw_tmp
      common/sudarwgt_flg/ flg_sudakov_rwgt, flg_sudarw_tmp
      double precision sudawgt, sudarwgtfac, uborns
      common/sudarwgt/ sudawgt, sudarwgtfac, uborns
