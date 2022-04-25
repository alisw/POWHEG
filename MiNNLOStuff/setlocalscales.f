      subroutine setlocalscales(iuborn,imode,rescfac)
c     returns the rescaling factor including sudakov form factors and
c     coupling rescaling, for Born (imode=1) and NLO corrections (imode=2)
c     (imode=3) is as imode=1, but with the computation of the d3 term. 
c      use libnnlops
      use rad_tools
      use pdfs_tools
      use sudakov_radiators
      use internal_parameters, pi_hoppet => pi, cf_hoppet => cf, ca_hoppet => ca, tf_hoppet => tf
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_st.h'
      include 'pwhg_pdf.h'
      include 'pwhg_flg.h'
      include 'pwhg_math.h'
      include 'minnlo_flg.h'
      integer iuborn,imode
      real * 8 rescfac,expsud,sudakov,pwhg_alphas,resc_coupl
      real * 8 ptb2,mb2,mu2,alphas,alphas_ptb2,b0,optb2,omb2,orescfac
     $     ,omuf2,alphas_dterms,alphas_sudakov,mu2dterms,mu2fact
      real * 8 pb(0:3), ptmax, yh, sudakov_form_factor
      integer oimode,i,flav
      common/flav_initial_state/ flav
      save optb2,omb2,orescfac,oimode,omuf2
      data optb2/-1d0/
      logical ini,bmass_in_minlo_flg
      common/c_bmass_in_minlo_flg/bmass_in_minlo_flg
      data ini/.true./
      save ini
      real * 8 powheginput,factsc2min,frensc2min,as,y,b1,tmp,bfact
      save factsc2min,frensc2min,b0,b1
      integer imax
      double precision D1, D2, D3, D1Q, D2Q, D3Q, xx1, xx2
      integer nflav, pdf_set
      parameter(nflav=5) ! BEWARE OF THE FLAVOUR: nflav must match LHAPDF's set
      double precision msqB(-nflav:nflav,-nflav:nflav), msqV1(-nflav:nflav,-nflav:nflav), msqV2(-nflav:nflav,-nflav:nflav)
      double precision pborn(1:nlegborn-1,1:4), L, pt, jacob
      character*100 pdf_name
      double precision d3terms
      common/d3terms/d3terms
      double precision pborn_UUB(0:3,nlegborn-1),e2,e,ebeam
      double precision Delta1,Delta2,modlog
      double precision s3, kappar, kappaf, kappaq, Q0_ini
      double precision mb2_sav, kappar_sav, kappaf_sav, kappaq_sav
      double precision kappaq_ini      
      save mb2_sav, kappar_sav, kappaf_sav, kappaq_sav, kappaq, kappaq_ini, Q0_ini
      double precision, parameter :: accuracy=1d-2
      double precision alphas_cutoff_fact, DSmear
      save alphas_cutoff_fact
      double precision, save :: Qsmear
      logical flg_profiledscales,flg_largeptscales,flg_rescaleQ0
      save flg_profiledscales,flg_largeptscales,flg_rescaleQ0
      real *8 ckappaq
      common/common_kappaQ/ckappaq
      real *8 ckill_D3
      common/common_kill_D3/ckill_D3      
      logical kill_D3,kill_H2,kill_D3_ini
      data kill_D3/.false./,kill_H2/.false./
      save kill_D3,kill_H2,kill_D3_ini

      logical flg_FOatQ         !: temporary placeholder for an option not yet fully implemented (eventually it'll be coded in common block)
      save flg_FOatQ
      
      
c     the only purpose of this is to invalidate the currently stored results
c     of setlocalscales (see the line  if(imode.eq.oimode ...))
      if(imode == -1) then
         oimode = -1
         return
      endif      

c     All the process dependence in this subroutine is and must be encoded through 
c     flg_minnloproc

      if (flg_minnloproc == 'H') then 
c     gg -> H production
c     Sudakov for a gluon 
         flav=0
      elseif (flg_minnloproc == 'Z') then 
c     Sudakov for a quark
         flav=1                 ! any value different from zero
      else
         write(*,*) ' setlocalscales: flg_minnloproc is ',flg_minnloproc
         write(*,*) ' not handled now, exiting ...'
         call exit(-1)
      endif

      if(ini) then
        
         ckappaq = -1d0
         ckill_D3 = -1d0
         if(powheginput("#kill_D3").eq.1) kill_D3 = .true.
         kill_D3_ini = kill_D3
         if(powheginput("#kill_H2").eq.1) kill_H2 = .true.         
c     set the number of flavour used in the Dterms
         call set_nflav(nflav)
         
c     this to make sure that at the line   if (.. abs(kappar-kappar_sav) > 1d-7 )
c     the if block is entered on the first call
         mb2_sav=-1d0
         kappar_sav=-1d0
         kappaf_sav=-1d0
         kappaq_sav=-1d0

         bmass_in_minlo_flg = powheginput("#bmass_in_minlo").eq.1

         b0=(33d0-2d0*st_nlight)/(12*pi)
         b1=(153d0-19d0*st_nlight)/(24*pi**2)
         ini = .false.

c     read the modified logarithms parameter p, if it's not found, set it to its default (p=6)
         modlog_p = powheginput("#modlog_p")
         if(modlog_p < -100d0) modlog_p = 6d0 !>> use standard modified logs by default
         write(*,*) '============================================='
         write(*,*) "using modlog_p = ", modlog_p
         write(*,*) '============================================='  

c     alphas_cutoff_fact**2 * pdf_q2min is the scale (in GeV^2) at which we freeze the running of alphas
c     pdf_q2min is the LHAPDF cutoff (e.g. for NNPDF30 it's 1GeV^2)
         alphas_cutoff_fact = powheginput('#alphas_cutoff_fact')
         
         if (.not. flg_hoppet_initialized) then            
            pdf_cutoff_fact = powheginput('#pdf_cutoff_fact')
c     If alphas_cutoff_fact is not set in the input card, set it to pdf_cutoff_fact (deprecated)
            if(alphas_cutoff_fact.lt.0d0) alphas_cutoff_fact=pdf_cutoff_fact
c     The conditions below should ultimately be replaced by a condition on pt/M
            if (flg_minnloproc == 'Z') then 
               if(pdf_cutoff_fact.lt.0d0) pdf_cutoff_fact=1.8d0
            elseif (flg_minnloproc == 'H') then 
               if(pdf_cutoff_fact.lt.0d0) pdf_cutoff_fact=2.5d0
            else
               write(*,*) ' setlocalscales: flg_minnloproc is ',flg_minnloproc
               call exit(-1)
            endif
         endif

c     If alphas_cutoff_fact is not set in the input card, use zero (Q0 will protect from very small pt values)
         if(alphas_cutoff_fact.lt.0d0) alphas_cutoff_fact=0d0

         write(*,*) 'freezing of alphas in setlocalscales at [GeV] '
     $        ,sqrt(alphas_cutoff_fact**2 *pdf_q2min)

         write(*,*)'freezing of PDF evolution in '
     $        ,'setlocalscales at [GeV]'
     $        ,sqrt(pdf_cutoff_fact**2 *pdf_q2min)

         if (flg_minnlo) then
            write(*,*) '========================='
            write(*,*) '*MINNLO* ACTIVATED'
            write(*,*) '========================='
         elseif (flg_minlo) then
            write(*,*) '====================================='
            write(*,*) 'MINLO ACTIVATED (BUT *MINNLO* IS OFF)'
            write(*,*) '====================================='
         endif
            
         if (flg_minlo.or.flg_minnlo) then
!     set the MiNNLO flags
            flg_include_delta_terms = .false.
            if(powheginput('#inc_delta_terms').eq.1)
     $           flg_include_delta_terms=.true.
            write(*,*) "include_delta_terms is ",flg_include_delta_terms

            flg_distribute_by_ub = .true.
            if(powheginput('#distribute_by_ub').eq.0)
     $           flg_distribute_by_ub = .false.
            write(*,*) "distribute_by_ub is ",flg_distribute_by_ub

            flg_distribute_by_ub_AP = .true.
            if(powheginput('#distribute_by_ub_AP').eq.0)
     $           flg_distribute_by_ub_AP = .false.
            write(*,*) "distribute_by_ub_AP is ",flg_distribute_by_ub_AP

            
            flg_dtermsallorders = .true.
            if(powheginput('#d3allorders').eq.0)
     $           flg_dtermsallorders = .false.
            write(*,*) '============================================='
            write(*,*) "flg_dtermsallorders = ",flg_dtermsallorders
            write(*,*) '============================================='

            kappaq = powheginput('#kappaQ')
            if(kappaq.lt.0d0) kappaq = 1d0
            write(*,*) '============================================='
            write(*,*) "using kappaQ = ",kappaq
            write(*,*) '============================================='
            kappaq_ini = kappaq
            
c This also calls init_anom_dim
            call init_Dterms(flav,flg_minnlo,kappaq)

c Setup profiled scales (main switch and parameters)    
            flg_profiledscales=.true.
            if(powheginput('#profiledscales').eq.0) flg_profiledscales=.false.
            Q0_ini = powheginput('#Q0')
            if(Q0_ini.lt.0d0) Q0_ini = 0d0
            if(Q0_ini.eq.0d0) Q0_ini = 1d-10
            npow = powheginput('#npow')
            if(npow.lt.0d0) npow = 1d0              
            call init_profiled_scales_parameters(flg_profiledscales,Q0_ini,npow)
            !>> flag that decides whether Q0 is rescaled by kappaf
            flg_rescaleQ0 = powheginput('#rescaleQ0').eq.1

            if(flg_profiledscales) then
              write(*,*) '============================================='
              write(*,*) 'Using profiled scales: (Q0,npow) = ',Q0_ini,npow
              write(*,*) '============================================='
            endif

c     Option to use mur and muf = pt in fixed order (Bbar) at large pt
            flg_largeptscales=powheginput('#largeptscales').eq.1
            if(flg_largeptscales) then
              write(*,*) '============================================='
               write(*,*) 'Using mu=pt in fixed order'
              write(*,*) '============================================='
            endif

            flg_FOatQ=powheginput('#FOatQ').eq.1
            if(flg_FOatQ) then
              write(*,*) '============================================='
            write(*,*)'Using hard scale in finite terms of fixed order'
              write(*,*) '============================================='
               write(*,*) 'this option is not functional at the moment'
               call exit(-1)
               if(flg_include_delta_terms) then
                  write(*,*) "ERROR: include_delta_terms not supported ",
     c                 "when flg_FOatQ is activated"
                  call exit(-1)
               endif
            endif

            Qsmear = powheginput("#Qsmear")
            if(Qsmear .gt. 0d0) then
              write(*,*) '============================================='
              write(*,*) 'Using smearing funtion with Qsmear = ', Qsmear
              write(*,*) '============================================='
              if(.not.flg_minnlo) then
                 write(*,*) 'Smearing is not functional with MiNLO'
                 call exit(-1)
              endif
           endif
        endif
      endif

      rescfac = 1

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC   These lines are process dependent!!!!!!
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c     pb(0:3) is the colourless "boson" momentum
      pb(:)=0d0
c     sum over colourless particles (they must all come from a single boson decay)
      do i=3,nlegborn
c     the sequence of colourless particles is unchanged in the Born and in the real flavour list
         if (abs(flst_born(i,iuborn)).gt.6) then
            if(flg_minlo_real) then
               pb(:)=pb(:) + kn_cmpreal(:,i)
!               pb(:)=pb(:) + kn_cmpborn(:,i) ! alternative (consistent) option to define pT in reals
            else
               pb(:)=pb(:) + kn_cmpborn(:,i)
            endif
         endif
      enddo

      ptb2 = pb(1)**2 + pb(2)**2 ! transverse momentum squared
      mb2  = pb(0)**2 - pb(3)**2 - ptb2 ! invariant mass squared

      if(imode.eq.oimode.and.ptb2.eq.optb2.and.mb2.eq.omb2) then
         rescfac=orescfac
         st_mufact2=omuf2
         return
      else
         optb2=ptb2
         omb2=mb2
         oimode=imode
      endif



c     first define the logarithms
      pt = sqrt(ptb2)


c     Here we assume that B(2) is fixed (not dependent upon the kinematics).
c     In this case, all the following is needed only when the D3 contribution
c     is evaluated.
      if(imode == 3) then
         call invISRmap(kn_pborn,nlegborn,pborn_UUB)
         call get_B_V1_V2(pborn_UUB,msqB,msqV1,msqV2)
         e2 =  kn_sbeams
         e = sqrt(e2)
         ebeam = e/2d0        
         xx1 = pborn_UUB(0,1)/ebeam
         xx2 = pborn_UUB(0,2)/ebeam   
      endif
C     this call arranges for the A and B resummation coefficients to be set
C     as well as the virtual corrections
      kappar = st_renfact
      kappaf = st_facfact

      !>> update Q0 scale according to KF (comment to avoid rescaling)
      if (flg_rescaleQ0) call reset_profiled_scales_parameters(Q0_ini/kappaf)
      
      if(ckappaq.gt.0) then
        kappaq = ckappaq
      else
        kappaq = kappaq_ini
      endif
      if (abs(mb2-mb2_sav) > 0d0 .or. abs(kappar-kappar_sav) > 0d0
     1     .or. abs(kappaf-kappaf_sav) > 0d0
     2     .or. abs(kappaq-kappaq_sav) > 0d0) then
         call reset_dynamical_parameters(dsqrt(mb2), kappar, kappaf, kappaq)
         call init_anom_dim
         mb2_sav = mb2
         kappar_sav = kappar
         kappaf_sav = kappaf
         kappaq_sav = kappaq
      endif

      if(imode == 3) then
         call virtual_scale_dep(msqB, msqV1, msqV2)
      endif

c     Here we compute the modified log with KR=1 as the KR dependence is
c     explicitly accounted for in the computation of the D terms
      call log_and_jacob(pt,cs%Q,L,jacob)
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     alpha_s reweighting
      if(flg_profiledscales) then
         ! scale entering the alphas passed to the dterms
         mu2dterms = st_renfact**2/kappaq**2 * (kappaq*sqrt(mb2) * exp(-L) + Q0 / (1d0 +
     $        (kappaq*sqrt(mb2)/Q0*exp(-L))**npow))**2
         ! scale entering pdfs of the fixed order
         mu2fact   = st_facfact**2/kappaq**2 * (kappaq*sqrt(mb2) * exp(-L) + Q0 / (1d0 +
     $        (kappaq*sqrt(mb2)/Q0*exp(-L))**npow))**2
      else
         mu2dterms = st_renfact**2*mb2*exp(-2d0*L)
         mu2fact   = st_facfact**2*mb2*exp(-2d0*L)
      endif

      if(flg_largeptscales) then
         ! those scales only enter the fixed order
         if(flg_profiledscales) then
            ! this now has a kappaq dependence at large pT --> not desirable
            mu2        = st_renfact**2/kappaq**2 * (sqrt(ptb2) + Q0 / (1d0 +
     $           (sqrt(ptb2)/Q0 )**npow))**2
            mu2fact = st_facfact**2/kappaq**2 * (sqrt(ptb2) + Q0 / (1d0 +
     $           (sqrt(ptb2)/Q0 )**npow))**2
         else
            ! this now has a kappaq dependence at large pT --> not desirable
            mu2     = st_renfact**2/kappaq**2 * ptb2
            mu2fact = st_facfact**2/kappaq**2 * ptb2
         endif
      else
c     This is needed when running with modified logs in the fixed-order
c     (Bbar) part
         mu2 = mu2dterms
      endif

c     If FOatQ is true, then the PDFs in the Bbar function have to be
c     evaluated at the hard scale, so we don't change st_mufact2.
      if(.not.flg_FOatQ .or.imode == 3) st_mufact2 = mu2fact
      omuf2=st_mufact2
      if(pdf_alphas_from_pdf) then
         alphas         = pwhg_alphas(mu2,st_lambda5MSB,st_nlight)
         alphas_dterms  = pwhg_alphas(mu2dterms,st_lambda5MSB,st_nlight)
         if((mu2.lt.alphas_cutoff_fact**2 *pdf_q2min) .and. (.not.flg_profiledscales)) then
            alphas         = pwhg_alphas(alphas_cutoff_fact**2 *pdf_q2min,st_lambda5MSB,st_nlight)
         endif
         if((mu2dterms.lt.alphas_cutoff_fact**2*pdf_q2min) .and. (.not.flg_profiledscales)) then
            alphas_dterms  = pwhg_alphas(alphas_cutoff_fact**2 *pdf_q2min,st_lambda5MSB,st_nlight)
         endif
      else
         write(*,*) "ERROR: pdf_alphas_from_pdf=false ",
     c        "not supported for minnlops"
         call exit(-1)
      endif

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      !>> ER+PM: use numerical integration for Sudakov
      sudakov_form_factor = Sudakov_pt_exact(L, alphas_cutoff_fact*sqrt(pdf_q2min), analytic_alphas=.false.)
      if(Qsmear .gt. 0d0) then
         sudakov_form_factor = sudakov_form_factor * Smearing(exp(-L),Qsmear,modlog_p)
      endif

      if(flg_FOatQ) then
c     It seems that in this case also as(pt) would be correct (??)
         Delta1  = expsudakov_pt(L)  * st_alpha
      else
         Delta1  = expsudakov_pt(L)  * alphas_dterms
      endif
      Delta2  = exp2sudakov_pt(L) * alphas_dterms**2
      rescfac = sudakov_form_factor

      if(imode.eq.2) then
         if(.not.flg_FOatQ) rescfac = rescfac * (alphas/st_alpha)**2
         if (flg_include_delta_terms.and.flg_minnlo) rescfac = rescfac * (1d0 + Delta1)
      elseif(imode == 1) then
         if(.not.flg_FOatQ) rescfac = rescfac * (alphas/st_alpha)
         if(.not.flg_bornonly) then
            if (.not. flg_minnlo) then 
               rescfac = rescfac * 
     1              (1 + Delta1 + alphas*b0*log(mu2/st_muren2))
            else
               if (flg_include_delta_terms) then
                  rescfac = rescfac *
     1                 (1d0 + (Delta1+alphas*b0*log(mu2/st_muren2) * (1d0 + Delta1)) + (Delta1**2/2d0+Delta2))
               else
                  if(flg_FOatQ) then
                     rescfac = rescfac * (1d0 + Delta1)
                  else
                     rescfac = rescfac * (1d0 + (Delta1+alphas*b0*log(mu2/st_muren2)))
                  endif
               end if
            endif
         endif
      endif

      if(ckill_D3.eq.0) then
        kill_D3 = .false.
      elseif(ckill_D3.eq.1) then
        kill_D3 = .true.
      else
        kill_D3 = kill_D3_ini
      endif
      
      if(imode.eq.3 .and. kill_D3) then
        d3terms = 0d0
        return
      endif
      
      if (flg_minnlo .and. imode .eq. 3) then
         
         call reset_dynamical_C_matxs()

         if(kill_H2) msqV2 = 0d0

         if(flg_dtermsallorders) then
            call DtermsAllOrders(D1, D2, D3, exp(-L), xx1, xx2, msqB, msqV1, msqV2, alphas_dterms)!, alphas_sudakov)
            if(flg_FOatQ) call D1D2atQ(D1Q, D2Q, exp(-L), xx1, xx2, msqB, msqV1, msqV2, st_alpha)
         else
            call Dterms(D1, D2, D3, exp(-L), xx1, xx2, msqB, msqV1, msqV2, alphas_dterms)!, alphas_sudakov)
            if(flg_FOatQ) call D1D2atQ(D1Q, D2Q, exp(-L), xx1, xx2, msqB, msqV1, msqV2, st_alpha)
         endif
        
         
         if (flg_include_delta_terms) then
            d3terms = D3 - D2 * Delta1 + D1 * Delta1**2/2d0 - D1 * Delta2
         else if (flg_uubornonly) then
            d3terms = D1+D2+D3
         else
            if(flg_FOatQ) then
c$$$               print*, D1*(1d0-st_alpha/alphas_dterms),D2*(1d0-(st_alpha/alphas_dterms)**2),D3*(1d0-(st_alpha/alphas_dterms)**3)
c$$$               print*, D1-D1Q,D2-D2Q,D3-D3Q
c$$$               print*, '------------------------------'
c$$$               d3terms = D3 + D1*(1d0-st_alpha/alphas_dterms)
c$$$     c                      + D2*(1d0-(st_alpha/alphas_dterms)**2)
c$$$     c                      + b0*D1*log(mu2/st_muren2)*(st_alpha/alphas_dterms) * st_alpha

               d3terms = D3 + D1 - D1Q
     c                      + D2 - D2Q
!     c                      + b0*D1Q*log(mu2/st_muren2)*st_alpha
!               print*, 2d0*log(exp(-L)),log(mu2/st_muren2)
c$$$               print*, b0, D1Q, log(mu2/st_muren2), st_alpha
c$$$               print*, b0*D1Q*log(mu2/st_muren2)*st_alpha
c$$$               print*, "======================================"               
            else
               d3terms = D3
            endif
         end if
         d3terms = d3terms * jacob * sudakov_form_factor

         if(Qsmear .gt. 0d0) then
            call DSmearing(DSmear, exp(-L), xx1, xx2, msqB, msqV1, msqV2, Qsmear, modlog_p)
            d3terms = d3terms + sudakov_form_factor * jacob * DSmear
         endif
         if(kill_D3) then
            print*, "ERROR: SHOULD NEVER GET HERE !!!!"
            d3terms = 0d0
         endif 
         return
      endif

      if(bmass_in_minlo_flg) then
         call bmass_in_minlo(bfact,alphas)
         rescfac = rescfac * bfact
      endif
      orescfac=rescfac

      end


      subroutine init_Dterms(flav,minnlo,kappaq)
      use rad_tools
      use pdfs_tools
      use internal_parameters, pi_hoppet => pi, cf_hoppet => cf, ca_hoppet => ca, tf_hoppet => tf
      implicit none
      include 'nlegborn.h'
      include 'PhysPars.h'
      include 'pwhg_pdf.h'
      include 'pwhg_kn.h'
      include 'pwhg_st.h'
      include 'minnlo_flg.h'
      include 'pwhg_flg.h'
      double precision D1, D2, D3, xx1, xx2, get_M_for_init_Dterms
      integer pdf_set, flav
      character*100 pdf_name
      character*100 string
      integer stringlength
      real *8 kappar,kappaf,kappaq
      logical ini,minnlo
      data ini/.true./
      save ini
      real *8 powheginput
      logical check_pdf
      parameter (check_pdf = .true.)
      

      !************************************************************************
      ! initialise pdf (already done in hoppetif.f, see there)
      string = " "
      call lhapdfname(pdf_ndns1,string,pdf_set)
      pdf_name = trim(string(1:stringlength(string)-1))//trim(".LHgrid")
      ! set the global mass of the colour singlet (this has to happen before the pdf initialisation)
      cs%M = get_M_for_init_Dterms()
      cs%Q = kappaq*get_M_for_init_Dterms()
      call init_masses_Dterms() ! initializes hard-coded mass values in NNLOPS_plugin
      if(.not. flg_hoppet_initialized) then
         if(flg_use_NNLOPS_pdfs) then
            call init_pdfs_NNLOPS(pdf_name, pdf_set, cutoff_high_fact=pdf_cutoff_fact)
         else
            call init_pdfs_from_LHAPDF(pdf_name, pdf_set, pdf_cutoff_fact)
         endif         
      endif
      
      if(check_pdf) call checkpdf


      !************************************************************************

      kappar=st_renfact
      kappaf=st_facfact

      if (flav.eq.0) then ! gg-initiated
         call set_process_and_parameters('pp', 'gg', sqrt(kn_sbeams), get_M_for_init_Dterms(), kappar, kappaf, kappaq) ! M, KR, KF, KQ
      else ! qqbar-initiated
         call set_process_and_parameters('pp', 'qq', sqrt(kn_sbeams), get_M_for_init_Dterms(), kappar, kappaf, kappaq) ! M, KR, KF, KQ
      endif

c     initialise anomalous dimensions
      call init_anom_dim
c     initialise coefficient functions
      if(flg_minnlo) call init_C_matxs()
      end

      subroutine log_and_jacob(pt,M,L,jacob)
      include 'minnlo_flg.h'
      double precision pt, M, L, jacob
      double precision a, b, c, d, g

      if (modlog_p > 0d0) then
         !>> standard modified log
         L     = 1d0/modlog_p*log((M/pt)**modlog_p + 1d0)
         jacob = 1d0/pt*(M/pt)**modlog_p/(1d0 + (M/pt)**modlog_p) ! dL/dpt
      else if (modlog_p .eq. -1d0) then
         !>> piecewise modified log
         if ((M .ge. pt).and.(pt .ge. M/2d0)) then
            a = 5d0
            b = -8d0/M
            c = 4d0/M**2
            g  = a + b*pt + c*pt**2
            L     = log(g)
            jacob = abs((-8*(M - pt))/(5*M**2 - 8*M*pt + 4*pt**2))
         else if (pt .ge. M) then
            L     = 0d0
            jacob = 0d0
         else
            L     = log(M/pt)
            jacob = 1d0/pt
         end if
      else if (modlog_p .eq. -2d0) then
         !>> piecewise modified log
         if ((M .ge. pt).and.(pt .ge. M/2d0)) then
            a = 4d0*(1d0-log(2d0))
            b = 8d0*(-2d0+3d0*log(2d0))/M
            c = 4d0*(5d0-9d0*log(2d0))/M**2
            d = 8d0*(-1d0+2d0*log(2d0))/M**3
            g  = a + b*pt + c*pt**2 + d*pt**3
            L     = g
            jacob = b + 2d0*c*pt + 3d0*d*pt**2
         else if (pt .ge. M) then
            L     = 0d0
            jacob = 0d0
         else
            L     = log(M/pt)
            jacob = 1d0/pt
         end if
      else if (modlog_p .eq. 0d0) then
         !>> theta function with unmodified log
         L     = 0d0
         jacob = 0d0
         if (pt <= M) then
            L     = log(M/pt)
            jacob = 1d0/pt
         end if
      else
         write(*,*) "ERROR: log_and_jacob, unrecognized modlog_p"
         call exit(-1)
      end if
      end

      subroutine sudakov_exponent
      write(*,*) ' sudakov_exponent should not be called here'
      call exit(-1)
      end
      

      subroutine checkpdf
      include 'pwhg_st.h'
      include 'pwhg_pdf.h'
      real *8 xxx(1:4),muF,pdf1(-pdf_nparton:pdf_nparton)
      integer i

      write(*,*) '==============================================='
      write(*,*) 'LHAPDF'
      write(*,*) '==============================================='

      xxx(1)=0.001d0
      xxx(2)=0.01d0
      xxx(3)=0.05d0
      xxx(4)=0.1d0

      muF=0.75d0
      st_mufact2 = muF**2
      write(*,*) 'muF = ',muF
      do i=1,4
         call pdfcall(1,xxx(i),pdf1)
         write(*,*) 'x,g,d,u,c,b -> ',xxx(i),pdf1(0),pdf1(1),pdf1(2),pdf1(4),pdf1(5)
      enddo

      muF=1.d0
      st_mufact2 = muF**2
      write(*,*) 'muF = ',muF
      do i=1,4
         call pdfcall(1,xxx(i),pdf1)
         write(*,*) 'x,g,d,u,c,b -> ',xxx(i),pdf1(0),pdf1(1),pdf1(2),pdf1(4),pdf1(5)
      enddo

      muF=2.d0
      st_mufact2 = muF**2
      write(*,*) 'muF = ',muF
      do i=1,4
         call pdfcall(1,xxx(i),pdf1)
         write(*,*) 'x,g,d,u,c,b -> ',xxx(i),pdf1(0),pdf1(1),pdf1(2),pdf1(4),pdf1(5)
      enddo
      
      muF=5.d0
      st_mufact2 = muF**2
      write(*,*) 'muF = ',muF
      do i=1,4
         call pdfcall(1,xxx(i),pdf1)
         write(*,*) 'x,g,d,u,c,b -> ',xxx(i),pdf1(0),pdf1(1),pdf1(2),pdf1(4),pdf1(5)
      enddo
      
      muF=91.d0
      st_mufact2 = muF**2
      write(*,*) 'muF = ',muF
      do i=1,4
         call pdfcall(1,xxx(i),pdf1)
         write(*,*) 'x,g,d,u,c,b -> ',xxx(i),pdf1(0),pdf1(1),pdf1(2),pdf1(4),pdf1(5)
      enddo
      
      write(*,*) '==============================='
      end
      
      

