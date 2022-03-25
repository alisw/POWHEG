      subroutine btildereal(xrad,resreal,www)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_flg.h'
      include 'pwhg_par.h'
      real * 8 xrad(3),resreal(maxprocborn),www
      real * 8 rr(maxalr),rc(maxalr),rp(maxalr),rm(maxalr),
     # rs(maxalr),rcs(maxalr),rps(maxalr),rms(maxalr),xl,xlp,xlm,
     # jac_over_csi,jac_over_csi_coll,jac_over_csi_soft,
     # jac_over_csi_p,jac_over_csi_m,rrr,rrrc,rrrs,rrrcs,
     # rrrp,rrrps,rrrm,rrrms,remnant,out0,out1
      real * 8, allocatable :: out1arr(:),out0arr(:)
      integer j,iuborn
      logical valid_emitter
      external valid_emitter
      logical pwhg_isfinite
      external pwhg_isfinite
      logical, save :: ini=.true.

      if(ini) then
         if(flg_analysisextrainfo) then
            allocate(out1arr(maxalr),out0arr(maxalr))
         endif
         ini = .false.
      endif

      do j=1,flst_nborn
         resreal(j)=0
      enddo
      do kn_emitter=0,nlegborn
c output values for analysis_driver
         out0=0
         out1=0
         if(flg_analysisextrainfo) then
            out0arr = 0
            out1arr = 0
         endif
c check that emitter is valid
         if(valid_emitter(kn_emitter)) then
            if(kn_emitter.gt.2) then
c     final state radiation
               call gen_real_phsp_fsr(xrad,jac_over_csi,
     #              jac_over_csi_coll,jac_over_csi_soft)
c This subroutine may set the scales with values depending
c upon the real emission kinematics
               call setscalesbtlreal
c sigreal fills the array rr with the value of the R_alpha contribution
c that have emitter equal to kn_emitter. All other contributions are set
c to zero. 
               call sigreal_btl(rr)
               if(flg_withsubtr) then
c We may prefer to set the counterterms scales different from the real scales
                  call setscalesbtlct
                  call collfsr(rc)
c     soft subtraction
                  call soft(rs)
                  call softcollfsr(rcs)
c     in final state radiation csimax is independent of y
                  xl=log(kn_csimax/par_csicut)
               endif
               do j=1,flst_nalr
                  iuborn=flst_alr2born(j)
                  rrr=rr(j)*kn_jacborn
     #                 *jac_over_csi/(1-kn_y)/kn_csitilde
                  if(flg_withsubtr) then
                     rrrc=rc(j)*kn_jacborn
     #                 *jac_over_csi_coll/(1-kn_y)/kn_csitilde
                     rrrs=rs(j)*kn_jacborn
     #                 *jac_over_csi_soft/(1-kn_y)/kn_csitilde
                     rrrcs=rcs(j)*kn_jacborn
     #                 *jac_over_csi_soft/(1-kn_y)/kn_csitilde
                     remnant=(rrrs-rrrcs)*xl*kn_csitilde
                  endif
                  if(flg_withsubtr) then
                     resreal(iuborn)= resreal(iuborn)+rrr-rrrc
     #-rrrs+rrrcs+remnant
                  else
c     provide a damping factor for the singular region,
c     to avoid divergent integral (25 is an ad hoc value
                     resreal(iuborn)= resreal(iuborn)
     #               +rrr*(1-kn_y**2)*kn_csi/
     #                  (25/kn_sbeams+(1-kn_y**2)*kn_csi)
                  endif
                  if(flg_nlotest) then
                     out1=out1+rrr
                     if(flg_analysisextrainfo) then
                        out1arr(j) = rrr
                     endif
                     if(flg_withsubtr) then
                        out0=out0-rrrc-rrrs+rrrcs+remnant
                        if(flg_analysisextrainfo) then
                           out0arr(j) = -rrrc-rrrs+rrrcs+remnant
                        endif
                     endif
                  endif
               enddo
            else
c     initial state singularities.
c     Regions that have only + (-) collinear singularity should return
c     zero rm (rp).
               call gen_real_phsp_isr
     #(xrad,jac_over_csi,jac_over_csi_p,jac_over_csi_m,
     #jac_over_csi_soft)
               call setscalesbtlreal
               call sigreal_btl(rr)
               if(flg_withsubtr) then
                  call setscalesbtlct
                  call soft(rs)
                  if(kn_emitter.ne.2) then
                     call collisrp(rp)
                     call softcollisrp(rps)
                  endif
                  if(kn_emitter.ne.1) then
                     call collisrm(rm)
                     call softcollisrm(rms)
                  endif
c     remnants (see xscaled.pdf in docs directory)
                  xl =log(kn_csimax/par_csicut)
                  xlp=log(kn_csimaxp/par_csicut)
                  xlm=log(kn_csimaxm/par_csicut)
               endif
               do j=1,flst_nalr
                  rrr=rr(j)*kn_jacborn
     #                 *jac_over_csi/(1-kn_y**2)/kn_csitilde
                  if(flg_withsubtr) then
                     rrrs=rs(j)*kn_jacborn
     #                 *jac_over_csi_soft/(1-kn_y**2)/kn_csitilde
                     remnant=rrrs*xl*kn_csitilde
                     if(kn_emitter.ne.2) then
                        rrrp=rp(j)*kn_jacborn
     #                 *jac_over_csi_p/(1-kn_y)/kn_csitilde/2
                        rrrps=rps(j)*kn_jacborn
     #                 *jac_over_csi_soft/(1-kn_y)/kn_csitilde/2
                        remnant=remnant-rrrps*xlp*kn_csitilde
                     else
                        rrrp=0
                        rrrps=0
                     endif
                     if(kn_emitter.ne.1) then
                        rrrm=rm(j)*kn_jacborn
     #                 *jac_over_csi_m/(1+kn_y)/kn_csitilde/2
                        rrrms=rms(j)*kn_jacborn
     #                 *jac_over_csi_soft/(1+kn_y)/kn_csitilde/2
                        remnant=remnant-rrrms*xlm*kn_csitilde
                     else
                        rrrm=0
                        rrrms=0
                     endif
                  endif
                  iuborn=flst_alr2born(j)
                  if(flg_withsubtr) then
                     resreal(iuborn)= resreal(iuborn)+rrr
     #              -rrrs-rrrp-rrrm+rrrps+rrrms+remnant
                  else
c     provide a damping factor for the singular region,
c     to avoid divergent integral (25 is an ad hoc value)
                     resreal(iuborn)= resreal(iuborn)
     #               +rrr*(1-kn_y**2)*kn_csi/
     #                  (25/kn_sbeams+(1-kn_y**2)*kn_csi)
                  endif
                  if(flg_nlotest) then
                     out1=out1+rrr
                     if(flg_analysisextrainfo) then
                        out1arr(j) = rrr
                     endif
                     if(flg_withsubtr) then
                        out0=out0-rrrs-rrrp-rrrm+rrrps+rrrms+remnant
                        if(flg_analysisextrainfo) then
                           out0arr(j) = -rrrs-rrrp-rrrm+rrrps+rrrms
     1                          +remnant
                        endif
                     endif
                  endif
               enddo
            endif
         endif
         if (.not.pwhg_isfinite(out0).or..not.pwhg_isfinite(out1)) then
            out0 = 0d0 
            out1 = 0d0 
            resreal = 0d0 
            if(flg_analysisextrainfo) then
               out1arr = 0
               out0arr = 0
            endif
         endif
         if(flg_nlotest) then
            out0=out0*www
            out1=out1*www
            call analysis_extrainfo('realct',flst_nalr,out0arr,www)
            if(out0.ne.0d0) call analysis_driver(out0,0)
            call analysis_extrainfo('real',flst_nalr,out1arr,www)
            if(out1.ne.0d0) call analysis_driver(out1,1)
         endif
      enddo
      end

      subroutine checklims(iun)
      implicit none
      integer iun
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_dbg.h'
       include 'pwhg_rnd.h'
       include 'pwhg_flg.h'
       include 'pwhg_par.h'
      real * 8 powheginput
      external sigreal_btl,soft,collfsr,softcollfsr, collisrp,
     $     softcollisrp,collisrm,softcollisrm,powheginput
      real * 8 tinyforlims, save_isrtinycsi,save_isrtinyy,
     1     save_fsrtinycsi,save_fsrtinyy
      logical savewithdamp
      savewithdamp = flg_withdamp
      flg_withdamp = .false.
      call randomsave
      if(rnd_cwhichseed /= 'none') then
         call setrandom(rnd_initialseed,rnd_i1,rnd_i2)
      endif
      if(powheginput("#chklimseed") .gt. 0) then
         call setrandom(nint(powheginput("#chklimseed")),0,0)
      endif
      tinyForLims=powheginput("#tinyForLims")
      if(tinyForLims .ge. 0) then
         save_isrtinycsi = par_isrtinycsi
         save_isrtinyy   = par_isrtinyy  
         save_fsrtinycsi = par_fsrtinycsi
         save_fsrtinyy   = par_fsrtinyy  
         
         par_isrtinycsi = tinyforlims 
         par_isrtinyy   = tinyforlims 
         par_fsrtinycsi = tinyforlims 
         par_fsrtinyy   = tinyforlims 
      endif

      if(dbg_softtest) then
         write(iun,*) '******************************************'   
         write(iun,*) '           CHECK  SOFT LIMITS             '     
         write(iun,*)
         do kn_emitter=0,nlegborn
            call checksoft(sigreal_btl,soft,' soft',iun)
         enddo
         write(iun,*) '******************************************'
         write(iun,*)
      endif
      
      if(dbg_colltest) then
         write(iun,*) '******************************************'   
         write(iun,*) '      CHECK  COLL. LIMITS FOR FSR       '         
         write(iun,*)
         do kn_emitter=3,nlegborn
            call checkcoll(sigreal_btl,collfsr,1,' coll',iun)
         enddo
         write(iun,*) '******************************************'   
         write(iun,*)
         write(iun,*) '******************************************'   
         write(iun,*) '      CHECK  COLL. LIMITS FOR ISR       '         
         write(iun,*)
         do kn_emitter=0,2
c            call randomsave
            if(kn_emitter.ne.2) call checkcoll(sigreal_btl,collisrp,1
     $           ,' coll-plus',iun)
c            call randomrestore
            if(kn_emitter.ne.1) call checkcoll(sigreal_btl,collisrm,-1
     $           ,' coll-minus',iun)
         enddo
         write(iun,*) '******************************************'   
         write(iun,*)
      endif
      
      if(dbg_softtest.and.dbg_colltest) then  
         write(iun,*) '******************************************'   
         write(iun,*) '   CHECK  SOFT-COLL. LIMITS FOR FSR     '         
         write(iun,*)
         do kn_emitter=3,nlegborn
            call checksoft(collfsr,softcollfsr,' soft-coll',iun)
         enddo
         write(iun,*) '******************************************'
         write(iun,*)
         write(iun,*) '******************************************'   
         write(iun,*) '   CHECK  SOFT-COLL. LIMITS FOR ISR +   '         
         write(iun,*)
         do kn_emitter=0,2
            if(kn_emitter.ne.2)call checksoft(collisrp,softcollisrp,
     $           ' soft-coll-plus',iun)
         enddo
         write(iun,*) '******************************************'
         write(iun,*)
         write(iun,*) '******************************************'   
         write(iun,*) '   CHECK  SOFT-COLL. LIMITS FOR ISR -   '         
         write(iun,*)
         do kn_emitter=0,2
            if(kn_emitter.ne.1)call checksoft(collisrm,softcollisrm,
     $           ' soft-coll-minus',iun)
         enddo
         write(iun,*) '******************************************'
         write(iun,*)

         
         write(iun,*) '******************************************'   
         write(iun,*) '   CHECK  COLL.-SOFT LIMITS FOR FSR     '         
         write(iun,*)
         do kn_emitter=3,nlegborn
            call checkcoll(soft,softcollfsr,1,' coll-soft',iun)
         enddo
         write(iun,*) '******************************************'
         write(iun,*)
         write(iun,*) '******************************************'   
         write(iun,*) '   CHECK  COLL.-SOFT LIMITS FOR ISR     '         
         write(iun,*)
         do kn_emitter=0,2
            if(kn_emitter.ne.2) call checkcoll(soft,softcollisrp,1
     $           ,' coll-plus-soft',iun)
            if(kn_emitter.ne.1) call checkcoll(soft,softcollisrm,-1
     $           ,' coll-minus-soft',iun)
         enddo
      endif
      call randomrestore
      flg_withdamp = savewithdamp
      if(tinyForLims .ge. 0) then
         par_isrtinycsi = save_isrtinycsi
         par_isrtinyy   = save_isrtinyy  
         par_fsrtinycsi = save_fsrtinycsi
         par_fsrtinyy   = save_fsrtinyy  
      endif
      end

      subroutine checkborn(iun)
c Check if Born, colour correlated born and spin correlated Born
c are consistent with total Born
      implicit none
      integer iun
      include 'nlegborn.h'
      include 'pwhg_math.h'
      include 'pwhg_flst.h'
      include 'pwhg_flg.h'
      include 'pwhg_br.h'
      include 'pwhg_st.h'
      integer  iborn,j,k,mu,kres,ires
      real * 8 tot
      real * 8 gtens(0:3,0:3),ap
      data gtens/1d0, 0d0, 0d0, 0d0,
     #           0d0,-1d0, 0d0, 0d0,
     #           0d0, 0d0,-1d0, 0d0,
     #           0d0, 0d0, 0d0,-1d0/
      save gtens
      logical colcorr
      external colcorr
      do ires=1,flst_nreson         
         do iborn=1,flst_nborn
            kres=flst_reslist(ires)
            do j=1,nlegborn
               if(colcorr(j,iborn,kres)) then
                  tot=0
                  do k=1,nlegborn
                     if(colcorr(k,iborn,kres)) then
                        if(k.ne.j) then
                           tot=tot+br_bornjk(j,k,iborn)
                        endif
                     endif
                  enddo
                  if(flst_born(j,iborn).eq.0) then
                     tot=tot/(ca*br_born(iborn))
                  else
                     tot=tot/(cf*br_born(iborn))
                  endif
                  if(abs((tot-1)/tot).gt.1d-8) then
                     write(iun,'(f6.3,a,20(i3,1x))') tot,
     1                    ' colour check fails for flav. struct:',kres,
     2                    (flst_born(k,iborn),k=1,nlegborn)
                  endif
               endif
            enddo
         enddo
      enddo
      do iborn=1,flst_nborn
         do j=1,nlegborn
            if(flst_born(j,iborn).eq.0) then
               tot=0
               do mu=0,3
                  tot=tot-gtens(mu,mu)*br_bmunu(mu,mu,j,iborn)
               enddo
               tot=tot/br_born(iborn)
               if(abs((tot-1)/tot).gt.1d-8) then
                  write(iun,'(f6.3,a,i2,a,20(i3,1x))')
     1    tot, ' spin correlated amplitude'//
     2 ' wrong for leg', j, ' flavour struct:',
     3                 (flst_born(k,iborn),k=1,nlegborn)
               endif
            endif
         enddo
      enddo
      end



      subroutine checksoft(sig,sigs,label,iun)
      implicit none
      include 'pwhg_dbg.h'
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_st.h'
      character *(*) label
      integer iun
      real * 8 xborn(ndiminteg-3),xrad(3)
      integer nexp
      parameter (nexp=5)
      real * 8 jac_over_csi,
     #jac_over_csi_coll,jac_over_csi_soft,rr(maxalr,nexp),
     #rs(maxalr,nexp),jac_over_csi_p,jac_over_csi_m
      integer j,jexp,alr,alrp
      character * 15 flag
      character * 32 fff
      logical ident(maxalr)
      real * 8 random,dotp
      external random,dotp
      logical valid_emitter,iszero,isnonzero,isequal,bad_born_kin,
     1     bad_real_kin
      external valid_emitter,bad_real_kin
      real * 8 tolparb,tolpar
c     The tolpar parameter sets the minimum in the ratio (transverse momentum)/(total energy)
c     for a kinematic configuration to be accepted for collinear tests. It is started
c     at 1, and for each failed configuration it is reduced by a factor 0.95 until an
c     acceptable configuration is found
      tolparb = 1
 1    continue
      tolparb = tolparb * 0.95d0
      if(tolparb == 0) then
c     We must abort this test
         write(*,*) 'checksoft: cannot complete the test; returning'
         return
      endif
c
      do
         do j=1,ndiminteg-3
            xborn(j)=random()
         enddo
         call gen_born_phsp(xborn)
         if(kn_jacborn /= 0d0) exit
      enddo
      
      if(bad_born_kin(tolparb)) goto 1
      call setscalesbtilde
      call allborn
      call checkborn(iun)
c      write(iun,*)' mass',sqrt(2*dotp(kn_pborn(0,3),kn_pborn(0,4)))
      tolpar = 1
 20   do j=1,3
         xrad(j)=random()
      enddo
      tolpar = tolpar * 0.95d0
c     If we cannot pass the soft limit tolerance test, regenerate
c     born configuration
      if(tolpar == 0) then
         goto 1
      endif
c Check soft limits
      if(valid_emitter(kn_emitter)) then
          do jexp=1,nexp
            xrad(1)=10d0**(-jexp)
            if(kn_emitter.gt.2) then
               call gen_real_phsp_fsr(xrad,jac_over_csi,
     $              jac_over_csi_coll,jac_over_csi_soft)
            else
               call gen_real_phsp_isr (xrad,jac_over_csi,jac_over_csi_p,
     $              jac_over_csi_m,jac_over_csi_soft)
            endif
            if(jexp == 1) then
               if(bad_real_kin(tolpar)) goto 20
               write(iun,*) ' Random Born variables ====> ',xborn
               write(iun,*) ' Random radiation variables ====> ',xrad
            endif
            write(iun,*) '### Check soft',xrad(1)
            call sig(rr(1,jexp))
            call sigs(rs(1,jexp))
         enddo
         do alr=1,flst_nalr
            ident(alr)=.false.
         enddo
         do alr=1,flst_nalr
            if(kn_emitter /= flst_emitter(alr)) cycle
c only radiated gluons or photons
            if(flst_alr(nlegreal,alr).ne.0 .and.
     1         flst_alr(nlegreal,alr).ne.22  ) cycle
            if(ident(alr)) cycle
            fff = '(a,1x,i3,1x,a, 20(1x,i3),a,a,a)'
            write(fff(15:17),'(i3)') nlegreal
            write(iun,fff)
     $           ' emitter ',kn_emitter, ', process ',
     $           (flst_alr(j,alr),j=1,nlegreal) !,', ',label,':'
            do alrp=alr+1,flst_nalr
               if(kn_emitter /= flst_emitter(alrp)) cycle
               isequal=.true.
               do jexp=1,nexp
                  if(rr(alr,jexp).ne.rr(alrp,jexp).or. rs(alr,jexp)
     $                 .ne.rs(alrp,jexp)) isequal=.false.
               enddo
               if(isequal) then
                  write(iun,'(2a,i2,a,20(1x,i3))') label,' emitter ',
     1                 flst_emitter(alrp),', process ',
     2             (flst_alr(j,alrp),j=1,nlegreal)
                  write(iun,'(2a,i2,a,20(1x,i3))') ' ==  ',' emitter ',
     1                 flst_emitter(alr),', process ',
     2             (flst_alr(j,alr),j=1,nlegreal)
                  ident(alrp)=.true.
               endif
            enddo
            do jexp=1,nexp
               call setwarnflag(
     1              rs(alr,jexp),rr(alr,jexp),jexp,2,flag)
               if(rs(alr,jexp) == 0 .or. rr(alr,jexp) /= rr(alr,jexp)
     1              .or. rs(alr,jexp) /= rs(alr,jexp)) then
                  write(iun,*) 'alr=',alr,rr(alr,jexp),
     2             '/',rs(alr,jexp),flag
               else
                  write(iun,*) 'alr=',alr,
     1             rr(alr,jexp)/rs(alr,jexp),flag
               endif
            enddo
         enddo
      endif
      end


      subroutine checkcoll(sig,sigc,idir,label,iun)
      implicit none
      integer iun
      include 'pwhg_dbg.h'
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      character *(*) label
      integer idir
      real * 8 xborn(ndiminteg-3),xrad(3)
      integer nexp
      parameter (nexp=8)
      real * 8 jac_over_csi,
     #jac_over_csi_coll,jac_over_csi_soft,rr(maxalr,nexp),
     #rc(maxalr,nexp),jac_over_csi_p,jac_over_csi_m
      integer j,jexp,alr,alrp
      real * 8 random
      external random
      logical ident(maxalr)
      character * 15 flag
      logical valid_emitter,iszero,isnonzero,isequal,bad_born_kin,
     1     bad_real_kin
      external valid_emitter,bad_born_kin,bad_real_kin
c     The tolpar parameter sets the minimum in the ratio (transverse momentum)/(total energy)
c     for a kinematic configuration to be accepted for collinear tests. It is started
c     at 1, and for each failed configuration it is reduced by a factor 0.95 until an
c     acceptable configuration is found
      real * 8 tolparb,tolpar
      tolparb = 1
 19   continue
      tolparb = tolparb * 0.95d0
      if(tolparb == 0) then
c     We must abort this test
         write(*,*) 'checkcoll: cannot complete the test; returning'
         return
      endif
      do
         do j=1,ndiminteg-3
            xborn(j)=random()
         enddo
         call gen_born_phsp(xborn)         
         if (kn_jacborn /= 0d0) exit
      enddo

      if(bad_born_kin(tolparb)) goto 19
      call setscalesbtilde
      call allborn
      tolpar = 1
 20   continue
      tolpar = tolpar * 0.95d0
c     If we cannot pass the soft limit tolerance test, regenerate
c     born configuration
      if(tolpar == 0) then
         goto 19
      endif
      
      do j=1,3
         xrad(j)=random()
      enddo
      if(valid_emitter(kn_emitter)) then
         do jexp=1,nexp
            if(idir.ne.-1) then
               xrad(2)=10d0**(-jexp)
            else
               xrad(2)=1-10d0**(-jexp)
            endif
            if(kn_emitter.gt.2) then
               call gen_real_phsp_fsr(xrad,jac_over_csi,
     #jac_over_csi_coll,jac_over_csi_soft)
            else
               call gen_real_phsp_isr
     #(xrad,jac_over_csi,jac_over_csi_p,jac_over_csi_m,
     #jac_over_csi_soft)
            endif
            if(jexp == 1) then
               if(bad_real_kin(tolpar)) goto 20
               write(iun,*) ' Random Born variables ====> ',xborn
               write(iun,*) ' Random radiation variables ====> ',xrad
            endif
            write(iun,*) '######### Check coll',xrad(2)
            call sig(rr(1,jexp))
            call sigc(rc(1,jexp))
         enddo
         do alr=1,flst_nalr
            ident(alr)=.false.
         enddo
         do alr=1,flst_nalr
            do jexp=1,nexp
               if(rc(alr,jexp).ne.rc(alr,1)) then
                  write(iun,*)
     #' checklims error : coll lim depends upon coll variable'
               endif
            enddo
         enddo
         do alr=1,flst_nalr
            if(kn_emitter /= flst_emitter(alr)) cycle
            if(ident(alr)) cycle
            write(iun,'(2a,i2,a,20(1x,i3))') label,'  emitter ',
     1           kn_emitter,', process ',(flst_alr(j,alr),j=1,nlegreal)
            do alrp=alr+1,flst_nalr
               if(kn_emitter /= flst_emitter(alrp)) cycle
               isequal=.true.
               do jexp=1,nexp
                  if(rr(alr,jexp).ne.rr(alrp,jexp).or.
     1                 rc(alr,jexp).ne.rc(alrp,jexp)) isequal=.false.
               enddo
               if(isequal) then
                  write(iun,'(2a,i2,a,20(1x,i3))') label,' emitter ',
     1                 flst_emitter(alrp),', process ',
     2             (flst_alr(j,alrp),j=1,nlegreal)
                  write(iun,'(2a,i2,a,20(1x,i3))') ' ==  ',' emitter ',
     1                 flst_emitter(alr),', process ',
     2             (flst_alr(j,alr),j=1,nlegreal)
                  ident(alrp)=.true.                     
               endif
            enddo
            do jexp=1,nexp
               call setwarnflag(rr(alr,jexp),rc(alr,jexp),
     1              jexp,1,flag)
c     Added this 'if' to be sure that no division by zero occurs
c     if((rr(alr,jexp-1)-rc(alr,jexp-1)).ne.0d0) then
c     write(iun,*) (rr(alr,jexp)-rc(alr,jexp))/
c     #(rr(alr,jexp-1)-rc(alr,jexp-1)),rc(alr,jexp)/rr(alr,jexp),flag
c     endif
               if(rc(alr,jexp) == 0 .or. rr(alr,jexp) /= rr(alr,jexp)
     1              .or. rc(alr,jexp) /= rc(alr,jexp)) then
                  write(iun,*) 'alr=',alr,rr(alr,jexp),'/',
     1             rc(alr,jexp),flag
               else
                  write(iun,*) 'alr=',alr,rr(alr,jexp)/rc(alr,jexp),flag
               endif
 1          enddo
         enddo
      endif
      end

      subroutine setwarnflag(r1,r2,jexp,jexpfirst,flag)
      implicit none
      character * 15 flag
      real * 8 r1,r2,dist
      integer jexp,jexpfirst
      if(r1 == 0 .and. r2 /= 0 .or. r1 /= 0 .and. r2 == 0) then
         dist = 1
      elseif(r1 == 0 .and. r2 == 0) then
         dist = 0
      else
         dist = abs(r1/r2-1)
      endif
      flag=' '
      if(dist.gt.0.01) then
         if(jexp.eq.jexpfirst) then
            if(dist.lt.0.1) then
               flag='*-WARN-*'
            else
               flag='*-WWARN-*'
            endif
         elseif(jexp.eq.3) then
            if(dist.lt.0.1) then
               flag='*-WWARN-*'
            else
               flag='*-WWWARN-*'
            endif
         elseif(jexp.ge.4) then
            if(dist.lt.0.1) then
               flag='*-WWWARN-*'
            elseif(dist.lt.0.3) then
               flag='*-WWWWARN-*'
            else
               flag='*-WWWWWARN-*'
            endif
         endif
      endif
      end

      logical function bad_born_kin(tolpar)
      implicit none
      real * 8 tolpar
      include 'pwhg_dbg.h'
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_st.h'
      include 'pwhg_pdf.h'
      logical isnotres(nlegborn)
      real * 8 p(0:3,nlegborn)
      real * 8 pdf(-pdf_nparton:pdf_nparton), save_stmufact2, e
      integer m,l,k
      bad_born_kin = .true.
c     In V2 the resonance structure is the same for all.
      isnotres = .true.
      do k=3,nlegborn
         if(flst_bornres(k,1)/=0) then
            isnotres(flst_bornres(k,1)) = .false.
         endif
      enddo
      l=0    
      do k=1,nlegborn
         if(isnotres(k)) then
            l=l+1
            p(:,l) = kn_cmpborn(:,k)
         endif
      enddo
c     exclude zero pdf's
      e = p(0,1)+p(0,2)
      save_stmufact2=st_mufact2
      st_mufact2 = e**2
      call pdfcall(1,kn_xb1,pdf)
      do k=-3,3
         if(pdf(k) == 0) return
      enddo
      call pdfcall(2,kn_xb2,pdf)
      st_mufact2 = save_stmufact2
      do k=-3,3
         if(pdf(k) == 0) return
      enddo

c     if any transverse momentum is too small, return
      if(l > 3) then ! exclude single resonance production!
         do k=3,l
            if(p(1,k)**2+p(2,k)**2 < (e*tolpar)**2) return
            if(l > 2) then
c     any transverse momentum of a pair is too small
               do m=k+1,l
                  if((p(1,k)+p(1,m))**2+(p(2,k)+p(2,m))**2
     1                 < (e*tolpar)**2) return
               enddo
            endif
         enddo
      endif
      bad_born_kin = .false.
      end

      logical function bad_real_kin(tolpar)
      implicit none
      real * 8 tolpar
      include 'pwhg_dbg.h'
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_st.h'
      include 'pwhg_pdf.h'
      logical isnotres(nlegreal)
      real * 8 p(0:3,nlegreal)
      real * 8 pdf(-pdf_nparton:pdf_nparton), save_stmufact2, e
      integer m,l,k
      
      bad_real_kin = .true.
c     In V2 the resonance structure is the same for all.
      isnotres = .true.
      do k=3,nlegreal
         if(flst_realres(k,1)/=0) then
            isnotres(flst_realres(k,1)) = .false.
         endif
      enddo
      l=0    
      do k=1,nlegreal
         if(isnotres(k)) then
            l=l+1
            p(:,l) = kn_cmpreal(:,k)
         endif
      enddo
c     exclude zero pdf's
      e = p(0,1)+p(0,2)
      save_stmufact2=st_mufact2
      st_mufact2 = e**2
      call pdfcall(1,kn_x1,pdf)
c tolerate zero pdf only for heavy quarks
      do k=-3,3
         if(pdf(k) == 0) return
      enddo
      call pdfcall(2,kn_x2,pdf)
      st_mufact2 = save_stmufact2
      do k=-3,3
         if(pdf(k) == 0) return
      enddo

c     if any transverse momentum is too small, return
      if(l > 3) then ! exclude single resonance production!
         do k=3,l
            if(p(1,k)**2+p(2,k)**2 < (e*tolpar)**2) return
            if(l > 2) then
c     any transverse momentum of a pair is too small
               do m=k+1,l
                  if((p(1,k)+p(1,m))**2+(p(2,k)+p(2,m))**2
     1                 < (e*tolpar)**2) return
               enddo
            endif
         enddo
      endif
      bad_real_kin = .false.
      end

      subroutine sigreal_rad(sig)
      implicit none
      real * 8 sig
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_rad.h'
      include 'pwhg_flg.h'
      include 'pwhg_par.h'
      include 'pwhg_pdf.h'
      real * 8 rr(maxalr),rc(maxalr),rs(maxalr),rcs(maxalr)
      integer alr,alrpr,iret,em
      integer nmomset,emitter
      parameter (nmomset=10)
      real * 8 res(nmomset,maxalr),preal(0:3,nlegreal,nmomset),cprop
      integer equivto(maxalr)
      real * 8 equivcoef(maxalr)
      integer j,k
      real * 8 sumdijinv,dampfac,r
      real * 8 pdf1(-pdf_nparton:pdf_nparton),
     1         pdf2(-pdf_nparton:pdf_nparton)
      real * 8 ptsq,pwhg_pt2,dijterm
      logical computed(maxalr)
      logical condition
      logical ini
      data ini/.true./
      save ini,equivto,equivcoef
      real * 8 powheginput
      external pwhg_pt2,dijterm,powheginput
      if(ini) then
         do alr=1,flst_nalr
            equivto(alr)=-1
         enddo
         if(flg_smartsig) then
            flg_in_smartsig = .true.
            call fillequivarrayrad(flst_nalr,equivto,equivcoef,iret)
            if(iret<0) then
               call randomsave
c     generate "nmomset" random real-phase space configurations
               call fillmomenta(nlegreal,nmomset,kn_masses,preal)
               do alr=1,flst_nalr
                  flst_cur_alr = alr
                  if(kn_emitter.eq.0) then
                     kn_resemitter=0
                  else
                     kn_resemitter=flst_alrres(nlegreal,alr)
                  endif
                  do j=1,nmomset
                     call realgr(
     1                    flst_alr(1,alr),preal(0,1,j),res(j,alr))
                  enddo
                  call compare_vecs(nmomset,alr,res,1,alrpr,cprop,iret)
                  if(iret.eq.0) then
c     they are equal
                     equivto(alr)=alrpr
                     equivcoef(alr)=1
                  elseif(iret.eq.1) then
c     they are proportional
                     equivto(alr)=alrpr
                     equivcoef(alr)=cprop
                  else
c     < 0 for unequal:
                     equivto(alr)=-1
                  endif
               enddo
               call randomrestore
               call printrealequivregions("rad",equivto,equivcoef)
c     Write equiv file, if required
               if(powheginput('#writeequivfile') == 1) then
                  call writeequivfile('rad',flst_nalr,equivto,equivcoef)
               endif
            endif
         endif
         flg_in_smartsig = .false.
         ini=.false.
      endif
c End initialization phase; compute graphs
      do alr=1,flst_nalr
         rr(alr)=0
      enddo
      call pdfcall(1,kn_x1,pdf1)
      call pdfcall(2,kn_x2,pdf2)
      if(flg_withdamp) then
         call collrad(rc)
         call collsoftrad(rcs)
         call softrad(rs)
      endif
      do alr=1,flst_nalr
         computed(alr)=.false.
      enddo
      do j=1,rad_alr_nlist
         alr=rad_alr_list(j)
         em=flst_emitter(alr)
c check if emitter corresponds to current radiation region (i.e. rad_kinreg):
         if((rad_kinreg.eq.1.and.em.le.2).or.(em.gt.2.and.
     #       flst_lightpart+rad_kinreg-2.eq.em))then
c check if we have a g -> Q Qbar splitting below threshold:
            if(em.gt.0) then
               if(flst_alr(em,alr)+flst_alr(nlegreal,alr).eq.0.and.
     #abs(flst_alr(em,alr)).ge.4) then
                  ptsq=pwhg_pt2()
                  if(abs(flst_alr(em,alr)).eq.4
     #  .and.ptsq.lt.rad_charmthr2.or.
     # abs(flst_alr(em,alr)).eq.5.and.ptsq.lt.rad_bottomthr2) then
                     rr(alr)=0
                     goto 995
                  endif
               endif
            endif
c ----------------
c Gimnastic to avoid problem with non-lazy evaluation of logical
c expressions in gfortran; replaces the line
c            if(equivto(alr).lt.0.or..not.computed(equivto(alr))) then
            if(equivto(alr).lt.0) then
               condition=.true.
            elseif(.not.computed(equivto(alr))) then
               condition=.true.
            else
               condition=.false.
            endif
            if(condition) then
               if(kn_emitter.eq.0) then
                  kn_resemitter=0
               else
                  kn_resemitter=flst_alrres(nlegreal,alr)
               endif
               flst_cur_alr = alr
               call realgr(flst_alr(1,alr),kn_cmpreal,rr(alr))
               sumdijinv=0
               do k=1,flst_allreg(1,0,alr)
                  sumdijinv=sumdijinv
     #+1/dijterm(flst_allreg(1,k,alr),flst_allreg(2,k,alr),alr)
               enddo
               rr(alr)=rr(alr)/dijterm(em,nlegreal,alr)/sumdijinv
               if(em.gt.2) then
                  if(flg_doublefsr) then
c    supply a factor E_em/(E_em+E_rad), times 2 if both gluons
                     rr(alr)=rr(alr)
     1                    *kn_cmpreal(0,kn_emitter)**par_2gsupp/
     2                    (kn_cmpreal(0,kn_emitter)**par_2gsupp
     3                    +kn_cmpreal(0,nlegreal)**par_2gsupp)
                     if(flst_alr(kn_emitter,alr).eq.0.and.
     1                    flst_alr(nlegreal,alr).eq.0) then
                        rr(alr)=rr(alr)*2
                     endif
                  else
c     If the emitter is in the final state, and if the emitted and emitter
c     are both gluons, supply a factor E_em/(E_em+E_rad) * 2
                     if(flst_alr(kn_emitter,alr).eq.0.and.
     1                    flst_alr(nlegreal,alr).eq.0) then
                        rr(alr)=rr(alr)*2
     1                       *kn_cmpreal(0,kn_emitter)**par_2gsupp/
     2                       (kn_cmpreal(0,kn_emitter)**par_2gsupp
     3                       +kn_cmpreal(0,nlegreal)**par_2gsupp)
                     endif
                  endif
               endif
               rr(alr)=rr(alr)*flst_mult(alr)
c supply Born zero damping factor, if required
               if(flg_withdamp) then
                  r=rr(alr)
                  call bornzerodamp(alr,r,rc(alr),rs(alr),rcs(alr),
     1                 dampfac)
                  rr(alr)=rr(alr) * dampfac
               endif
               computed(alr)=.true.
               if(equivto(alr).gt.0) then
                  rr(equivto(alr))=rr(alr)/equivcoef(alr)
                  computed(equivto(alr))=.true.
               endif
            else
               rr(alr)=rr(equivto(alr))*equivcoef(alr)
            endif
         else
            rr(alr)=0
         endif
 995     continue
      enddo
      sig=0
      do j=1,rad_alr_nlist
         alr=rad_alr_list(j)
         if(rr(alr).ne.0) then
            rr(alr)=rr(alr)*pdf1(flst_alr(1,alr))*pdf2(flst_alr(2,alr))
            sig=sig+rr(alr)
            rad_real_arr(j)=rr(alr)
         else
            rad_real_arr(j)=0
         endif
      enddo
      end


      subroutine sigreal_btl(rr)
      implicit none
      real * 8 rr(*)
      call sigreal_btl0(rr,0)
      end

c Real cross section, required by btilde;
c fills the array rr(alr) with the invariant cross section, multiplied
c by csi^2 (1-y^2) for ISR regions
c    csi^2 (1-y)   for FSR regions
      subroutine sigreal_btl0(rr,imode)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_flg.h'
      include 'pwhg_par.h'
      include 'pwhg_pdf.h'
      integer imode
      real * 8 rr(maxalr)
      real * 8 rc(maxalr),rs(maxalr),rcs(maxalr),r
      integer alr,alrpr,iret
      integer nmomset
      parameter (nmomset=10)
      real * 8 res(nmomset,maxalr),preal(0:3,nlegreal,nmomset),cprop
      integer equivto(maxalr),markused(maxalr)
      real * 8 equivcoef(maxalr)
      common/cequivtoreal/equivto
      common/cequivcoefreal/equivcoef
      integer j,k
      real * 8 sumdijinv,dampfac
      real * 8 pdf1(-pdf_nparton:pdf_nparton),
     1         pdf2(-pdf_nparton:pdf_nparton)
      real * 8 rescfac
      logical ini
      data ini/.true./
      save ini,/cequivtoreal/,/cequivcoefreal/
      real * 8 dijterm,powheginput
      external dijterm,powheginput
      if(ini) then
         do alr=1,flst_nalr
            equivto(alr)=-1
         enddo
         if(flg_smartsig) then
            flg_in_smartsig = .true.
            call fillequivarraybtl(flst_nalr,equivto,equivcoef,iret)
            if(iret<0) then
               call randomsave
c     generate "nmomset" random real-phase space configurations
               call fillmomenta(nlegreal,nmomset,kn_masses,preal)
               do alr=1,flst_nalr
                  flst_cur_alr = alr
                  do j=1,nmomset
                     if(kn_emitter.eq.0) then
                        kn_resemitter=0
                     else
                        kn_resemitter=flst_alrres(nlegreal,alr)
                     endif
                     call realgr(
     1                    flst_alr(1,alr),preal(0,1,j),res(j,alr))
                  enddo
                  call compare_vecs(nmomset,alr,res,0,alrpr,cprop,iret)
                  if(iret.eq.0) then
c     they are equal:
                     equivto(alr)=alrpr
                     equivcoef(alr)=1
                  elseif(iret.eq.1) then
c     they are proportional:
                     equivto(alr)=alrpr
                     equivcoef(alr)=cprop
                  else
c     < 0 for unequal:
                     equivto(alr)=-1
                  endif
               enddo
               call randomrestore
               call printrealequivregions('btl',equivto,equivcoef)
c     Write equiv file, if required
               if(powheginput('#writeequivfile') == 1) then
                  call writeequivfile('btl',flst_nalr,equivto,equivcoef)
               endif
            endif
         flg_in_smartsig = .false.
         endif
         ini=.false.
      endif
c End initialization phase; compute graphs
      do alr=1,flst_nalr
         rr(alr)=0
         markused(alr)=0
      enddo
      if(flg_withdamp) then
         call collbtl(rc)
         call collsoftbtl(rcs)
         call softbtl(rs)
      endif
      do alr=1,flst_nalr
c Only R_alpha (namely alr) with the current emitter: 
         if(flst_emitter(alr).eq.kn_emitter) then
            if(equivto(alr).lt.0) then
c Not equal to any previous one, compute explicitly.
c First mark as being computed
               markused(alr)=1
               if(kn_emitter.eq.0) then
                  kn_resemitter=0
               else
                  kn_resemitter=flst_alrres(nlegreal,alr)
               endif
               flst_cur_alr = alr
               call realgr(flst_alr(1,alr),kn_preal,rr(alr))
c Supply FKS factor to separate singular region:
               sumdijinv=0
c Loop over all singular regions of the given contribution
               do k=1,flst_allreg(1,0,alr)
c flst_allreg({1,2},...) are the two legs that identify the k'th region
                  sumdijinv=sumdijinv
     1      +1/dijterm(flst_allreg(1,k,alr),flst_allreg(2,k,alr),alr)
               enddo
               rr(alr)=rr(alr)/dijterm(kn_emitter,nlegreal,alr)
     1              /sumdijinv
c If the emitter is in the final state, and if the emitted and emitter
c are both gluons, supply a factor E_em/(E_em+E_rad) * 2
               if(kn_emitter.gt.2) then
                  if(flg_doublefsr) then
c    supply a factor E_em/(E_em+E_rad), times 2 if both gluons
                     rr(alr)=rr(alr)
     1                    *kn_cmpreal(0,kn_emitter)**par_2gsupp/
     2                    (kn_cmpreal(0,kn_emitter)**par_2gsupp
     3                    +kn_cmpreal(0,nlegreal)**par_2gsupp)
                     if(flst_alr(kn_emitter,alr).eq.0.and.
     1                    flst_alr(nlegreal,alr).eq.0) then
                        rr(alr)=rr(alr)*2
                     endif
                  else
c     If the emitter is in the final state, and if the emitted and emitter
c     are both gluons, supply a factor E_em/(E_em+E_rad) * 2
                     if(flst_alr(kn_emitter,alr).eq.0.and.
     1                    flst_alr(nlegreal,alr).eq.0) then
                        rr(alr)=rr(alr)*2
     1                       *kn_cmpreal(0,kn_emitter)**par_2gsupp/
     2                       (kn_cmpreal(0,kn_emitter)**par_2gsupp
     3                       +kn_cmpreal(0,nlegreal)**par_2gsupp)
                     endif
                  endif
               endif
c supply Born zero damping factor, if required
               if(flg_withdamp) then
                  if(kn_emitter.gt.2) then
                     r=rr(alr)*(1-kn_y)*kn_csi**2
                  else
                     r=rr(alr)*(1-kn_y**2)*kn_csi**2
                  endif
                  r=r*flst_mult(alr)
                  call bornzerodamp(alr,r,rc(alr),rs(alr),rcs(alr),
     1                 dampfac)
                  if(imode.eq.0) then
                     rr(alr) =rr(alr) * dampfac
                  elseif(imode.eq.1) then
                     rr(alr) =rr(alr) * (1-dampfac)
                  else
                     write(*,*) ' sigreal_btl0: improper call'
                  endif
               endif
            else
               if(markused(equivto(alr)).ne.1) then
                  write(*,*) ' error: sigreal_btl flg_smartsig bug'
                  call exit(1)
               endif
               rr(alr)=rr(equivto(alr))*equivcoef(alr)
               markused(alr)=1
            endif
         endif
      enddo
      if(.not.flg_minlo) then
         rescfac = 1
         call pdfcall(1,kn_x1,pdf1)
         call pdfcall(2,kn_x2,pdf2)
      endif
      do alr=1,flst_nalr
         if(flg_minlo) then
            flg_minlo_real=.true.
            call setlocalscales(flst_alr2born(alr),2,rescfac)
            flg_minlo_real=.false.
            call pdfcall(1,kn_x1,pdf1)
            call pdfcall(2,kn_x2,pdf2)
         endif
         rr(alr)=rr(alr)*flst_mult(alr)
         if(kn_emitter.gt.2) then
            rr(alr)=rr(alr)*(1-kn_y)*kn_csi**2
         else
            rr(alr)=rr(alr)*(1-kn_y**2)*kn_csi**2
         endif
c include pdf's
         rr(alr)=rr(alr)*pdf1(flst_alr(1,alr))*pdf2(flst_alr(2,alr))
     1        *rescfac
      enddo
      end





      subroutine fillmomenta(nparticles,nmomset,masses,p)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_kn.h'
      integer nparticles,nmomset
      real * 8 p(0:3,nparticles,nmomset),masses(nparticles)
      integer mu,j,k,i
      real * 8 pl,en,scale,ptki
      real * 8 random
      external random
      real * 8 ptmin,pt,ptcut,beta, vec(3),plab(0:3,nparticles),
     $     pcm(0:3,nparticles)!,costh,modk,modi
      integer last,sign
      logical debug
      parameter (debug=.false.)
      real * 8 dotp
      external dotp

c     produce ONLY events with minimum pt larger than ptcut
      ptcut=sqrt(kn_sbeams)/((nparticles+2)*10)
      if(nparticles.eq.3) ptcut = 0d0
      if (debug)  write(*,*) 'ptcut min in fillmomenta = ',ptcut

c     average value of the generated momentum.
      scale = sqrt(kn_sbeams)/(nparticles+2)
      do j=1,nmomset
 1       ptmin=1d30
         last = (nparticles-2)*random()+1
         last=last+2
         do k=3,nparticles
            if (k.ne.last) then
               do mu=1,3
                  plab(mu,k)=scale*(random()-0.5d0)
               enddo
            endif
         enddo
         plab(3,last)=scale*(random()-0.5d0)
         plab(1,last)=0d0
         plab(2,last)=0d0
         do k=3,nparticles
            if (k.ne.last) then
               plab(1,last)=plab(1,last)-plab(1,k)
               plab(2,last)=plab(2,last)-plab(2,k)
            endif
         enddo
         pl=0d0
         en=0d0
         do k=3,nparticles
            plab(0,k)=sqrt(masses(k)**2+
     #         plab(1,k)**2+plab(2,k)**2+plab(3,k)**2)
            pl=pl+plab(3,k)
            en=en+plab(0,k)
         enddo
         plab(0,1)=(en+pl)/2
         plab(0,2)=(en-pl)/2
         plab(3,1)=plab(0,1)
         plab(3,2)=-plab(0,2)
         plab(1,1)=0d0
         plab(1,2)=0d0
         plab(2,1)=0d0
         plab(2,2)=0d0

c         do k=1,nparticles
c            do mu=0,3
c               plab(mu,k)=plab(mu,k)*scale
c            enddo
c         enddo
c     boost momenta in the center-of-mass frame

         beta=-(plab(0,1)-plab(0,2))/(plab(0,1)+plab(0,2))
         vec(1)=0d0
         vec(2)=0d0
         vec(3)=1d0
         call mboost(nparticles,vec,beta,plab(0,1),pcm(0,1))
         
c     compute the minimum pt wrt the beam axis
         do k=3,nparticles
            pt=sqrt(pcm(1,k)**2+pcm(2,k)**2)
            ptmin=min(ptmin,pt)
         enddo
         if (ptmin.lt.ptcut) then
c            write(*,*) 'ptmin = ',ptmin            
            goto 1
         endif
c     compute the minimum pt among each pair of final-state momenta
         do k=3,nparticles-1
            do i=k+1,nparticles
c               modk=sqrt(pcm(1,k)**2+pcm(2,k)**2+pcm(3,k)**2)
c               modi=sqrt(pcm(1,i)**2+pcm(2,i)**2+pcm(3,i)**2)
c               costh=
c     $          (pcm(1,k)*pcm(1,i)+pcm(2,k)*pcm(2,i)+pcm(3,k)*pcm(3,i))/
c     $              modk/modi
c               ptmin=min(ptmin,modk*sqrt(abs(1-costh**2)),
c     $                         modi*sqrt(abs(1-costh**2)))
               ptki = sqrt(abs(2*dotp(pcm(0,k),pcm(0,i))*
     $              pcm(0,k)*pcm(0,i)/(pcm(0,k)+pcm(0,i))**2))               
               ptmin=min(ptmin,ptki)
            enddo
         enddo
         if (ptmin.lt.ptcut) then
c            write(*,*) 'ptmin = ',ptmin            
            goto 1
         endif

         if (debug) then
            write(*,*) 'set of momenta # ',j
            do k=1,nparticles
               write(*,'(10d12.4)') pcm(:,k), 
     $              sqrt(abs(pcm(0,k)**2-pcm(1,k)**2-pcm(2,k)**2-
     $              pcm(3,k)**2))
            enddo
            write(*,*) 'last particle ',last
         endif         
         p(:,:,j)=pcm(:,:)
      enddo
      end
 

      subroutine compare_vecs(nmomset,alr,res,imode,alrpr,cprop,iret)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      real * 8 ep
      integer nmomset,alr,alrpr,imode,iret,j,k
      real * 8 res(nmomset,*),cprop,rat
      real * 8 powheginput,tmp
c imode=0 when called from btilde,
c imode=1 when called for radiation. In the latter
c case, graphs that do not have the same underlying Born
c     are not considered.
      tmp = powheginput("#compare_vecsr_ep")
      if(tmp>0) then
         ep = tmp
      else
         ep = 1d-8
      endif
      do j=1,alr-1
         if(flst_emitter(j).ne.flst_emitter(alr)) goto 10
         if(imode.eq.1.and.flst_alr2born(j).ne.flst_alr2born(alr))
     1        goto 10
         rat=res(1,alr)/res(1,j)
         do k=1,nmomset
            if(abs(1-res(k,alr)/res(k,j)/rat).gt.ep) goto 10
         enddo
         if(abs(1-rat).lt.ep) then
            iret=0
            cprop=1
         else
            iret=1
            cprop=rat
         endif
         alrpr=j
         return
 10      continue
      enddo
      iret=-1
      end

      subroutine realgr(rflav,p,res)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      integer rflav(nlegreal)
      real * 8 p(0:3,nlegreal),res
      call real_ampsq(p,rflav,res)
c flux factor
      res=res/(8*p(0,1)*p(0,2))
      if(res.eq.0) then
c         write(*,*) 'realgr:', rflav
         continue
      endif
      end



      subroutine real_ampsq(p,rflav,amp2)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_math.h'
      include 'pwhg_st.h'
      real * 8 p(0:3,nlegreal)
      integer rflav(nlegreal)
      real * 8 amp2 
      logical pwhg_isfinite
      external pwhg_isfinite
      call setreal(p,rflav,amp2)
c     check if amp2 is finite
      if (.not.pwhg_isfinite(amp2)) amp2=0d0
      amp2 = amp2*st_alpha/(2*pi)
      end


      subroutine printrealequivregions(flag,equivto,equivcoef)
c it prints the set of equivalent alr regions
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_rnd.h'
      character * 3 flag
      integer equivto(maxalr)
c      common/cequivtoreal/equivto
      real * 8 equivcoef(maxalr)
c      common/cequivcoefreal/equivcoef
      integer j,k,iun,count
      save count
      data count/0/
      call newunit(iun)
      write(*,*) 'Writing realequivregions-'//flag//' file...'
      if(rnd_cwhichseed == 'none') then
         open(unit=iun,file='realequivregions-'//flag,status='unknown')
      else
         open(unit=iun,file='realequivregions-'//flag//"-"
     1        //trim(rnd_cwhichseed),status='unknown')
      endif
      do j=1,flst_nalr
         if(equivto(j).eq.-1) then
            write(iun,'(a)')
     1           'Beginning sequence of equivalent amplitudes'
c            write(iun,100) 1d0,j, flst_alr(:,j)
            write(iun,101) j, 1d0, flst_alr(:,j)
            do k=1,flst_nalr
               if(equivto(k).eq.j) then
c                  write(iun,100) equivcoef(k),k,flst_alr(:,k)
                  write(iun,101) k,equivcoef(k),flst_alr(:,k)
               endif
            enddo
            count=count+1
         endif
      enddo
      write(iun,*) ''
      write(iun,'(a,i4,a)') 'Found ',count, ' equivalent groups'
      close(iun)
      write(*,*) 'Done'
c 100  format(d11.4,5x,i4,5x,100(i4,1x))
 101  format(i4,5x,d11.4,5x,100(i4,1x))
      end


