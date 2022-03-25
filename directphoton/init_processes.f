      subroutine init_processes
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_flg.h'
      include 'pwhg_rwgtsuda.h'
      integer i, j, k
      double precision powheginput
c     number of quark flavours
      integer NF
      parameter (NF = 5)

c     position of first light parton in flavour structure
c     the photon is considered a light parton in case of QED radiation
      flst_lightpart = 3        ! all quarks massless

c     number of inequivalent flavour structures
      flst_nborn = (6*NF+11)*NF ! QED = 6*NF, QED+QCD = (6*NF+11)*NF
      flst_nreal = (6*NF+5)*NF  ! (6*NF+5)*NF

      if (flst_nborn.GT.maxprocborn) then
         write(*,*) 'init_processes: increase maxprocborn'
         call exit(1)
      endif
      if (flst_nreal.GT.maxprocreal) then
         write(*,*) 'init_processes: increase maxprocreal'
         call exit(1)
      endif

      k = 0
      do i = 1, NF
c     Born QED
c     --------- q qbar -> gamma g
         flst_born(1, i) = i    ! q
         flst_born(2, i) = -i   ! qbar
         flst_born(3, i) = 22   ! gluon
         flst_born(4, i) = 0    ! photon
c     --------- qbar q -> gamma g
         flst_born(1, i+NF) = -i
         flst_born(2, i+NF) = i
         flst_born(3, i+NF) = 22
         flst_born(4, i+NF) = 0
c     --------- q g -> gamma q
         flst_born(1, i+2*NF) = i
         flst_born(2, i+2*NF) = 0
         flst_born(3, i+2*NF) = 22
         flst_born(4, i+2*NF) = i
c     --------- qbar g -> gamma qbar
         flst_born(1, i+3*NF) = -i
         flst_born(2, i+3*NF) = 0
         flst_born(3, i+3*NF) = 22
         flst_born(4, i+3*NF) = -i
c     --------- g q -> gamma q
         flst_born(1, i+4*NF) = 0
         flst_born(2, i+4*NF) = i
         flst_born(3, i+4*NF) = 22
         flst_born(4, i+4*NF) = i
c     --------- g qbar -> gamma qbar
         flst_born(1, i+5*NF) = 0
         flst_born(2, i+5*NF) = -i
         flst_born(3, i+5*NF) = 22
         flst_born(4, i+5*NF) = -i

c     Born QCD
c     q qbar annihilation
c     --------- q qbar -> g g
         flst_born(1, i+6*NF) = i 
         flst_born(2, i+6*NF) = -i
         flst_born(3, i+6*NF) = 0
         flst_born(4, i+6*NF) = 0
c     --------- qbar q -> g g
         flst_born(1, i+7*NF) = -i
         flst_born(2, i+7*NF) = i
         flst_born(3, i+7*NF) = 0
         flst_born(4, i+7*NF) = 0
         do j = 1, NF
c     --------- q qbar -> q' qbar' (incl. q qbar -> q qbar)
            flst_born(1, i+7*NF+j*NF) = i
            flst_born(2, i+7*NF+j*NF) = -i
            flst_born(3, i+7*NF+j*NF) = j
            flst_born(4, i+7*NF+j*NF) = -j
c     --------- qbar q -> q' qbar' (incl. qbar q -> q qbar)
            flst_born(1, i+(NF+7)*NF+j*NF) = -i
            flst_born(2, i+(NF+7)*NF+j*NF) = i
            flst_born(3, i+(NF+7)*NF+j*NF) = j
            flst_born(4, i+(NF+7)*NF+j*NF) = -j
c     q q' scattering
c     --------- q q' -> q q' (incl. q q -> q q)
            flst_born(1, i+(2*NF+7)*NF+j*NF) = i
            flst_born(2, i+(2*NF+7)*NF+j*NF) = j
            flst_born(3, i+(2*NF+7)*NF+j*NF) = i
            flst_born(4, i+(2*NF+7)*NF+j*NF) = j
c     --------- qbar qbar' -> qbar qbar' (incl. qb qb -> qb qb)
            flst_born(1, i+(3*NF+7)*NF+j*NF) = -i
            flst_born(2, i+(3*NF+7)*NF+j*NF) = -j
            flst_born(3, i+(3*NF+7)*NF+j*NF) = -i
            flst_born(4, i+(3*NF+7)*NF+j*NF) = -j

            if (i.NE.j) then
               k = k + 1
c     --------- q qbar' -> q qbar'
               flst_born(1, (4*NF+8)*NF+k) = i
               flst_born(2, (4*NF+8)*NF+k) = -j
               flst_born(3, (4*NF+8)*NF+k) = i
               flst_born(4, (4*NF+8)*NF+k) = -j
c     --------- qbar q' -> q' qbar
               flst_born(1, (5*NF+7)*NF+k) = -i
               flst_born(2, (5*NF+7)*NF+k) = j
               flst_born(3, (5*NF+7)*NF+k) = j
               flst_born(4, (5*NF+7)*NF+k) = -i
            endif
         enddo

c     q g scattering
c     --------- q g -> q g
         flst_born(1, i+(6*NF+6)*NF) = i
         flst_born(2, i+(6*NF+6)*NF) = 0
         flst_born(3, i+(6*NF+6)*NF) = i
         flst_born(4, i+(6*NF+6)*NF) = 0
c     --------- g q -> q g
         flst_born(1, i+(6*NF+7)*NF) = 0
         flst_born(2, i+(6*NF+7)*NF) = i
         flst_born(3, i+(6*NF+7)*NF) = i
         flst_born(4, i+(6*NF+7)*NF) = 0
c     --------- qbar g -> qbar g
         flst_born(1, i+(6*NF+8)*NF) = -i
         flst_born(2, i+(6*NF+8)*NF) = 0
         flst_born(3, i+(6*NF+8)*NF) = -i
         flst_born(4, i+(6*NF+8)*NF) = 0
c     --------- g qbar -> qbar g
         flst_born(1, i+(6*NF+9)*NF) = 0
         flst_born(2, i+(6*NF+9)*NF) = -i
         flst_born(3, i+(6*NF+9)*NF) = -i
         flst_born(4, i+(6*NF+9)*NF) = 0

c     gluon fusion
c     --------- g g -> q qbar
         flst_born(1, i+(6*NF+10)*NF) = 0
         flst_born(2, i+(6*NF+10)*NF) = 0
         flst_born(3, i+(6*NF+10)*NF) = i
         flst_born(4, i+(6*NF+10)*NF) = -i
      enddo
c$$$c     --------- g g -> g g
c$$$      flst_born(1, (6*NF+11)*NF+1) = 0
c$$$      flst_born(2, (6*NF+11)*NF+1) = 0
c$$$      flst_born(3, (6*NF+11)*NF+1) = 0
c$$$      flst_born(4, (6*NF+11)*NF+1) = 0

c     Real
      k = 0
      do i = 1, NF
c     q qbar annihilation
c     --------- q qbar -> gamma g g
         flst_real(1, i) = i 
         flst_real(2, i) = -i
         flst_real(3, i) = 22
         flst_real(4, i) = 0
         flst_real(5, i) = 0
c     --------- qbar q -> gamma g g
         flst_real(1, i+NF) = -i
         flst_real(2, i+NF) = i
         flst_real(3, i+NF) = 22
         flst_real(4, i+NF) = 0
         flst_real(5, i+NF) = 0
         do j = 1, NF
c     --------- q qbar -> gamma q' qbar' (incl. q qbar -> gamma q qbar)
            flst_real(1, i+NF+j*NF) = i
            flst_real(2, i+NF+j*NF) = -i
            flst_real(3, i+NF+j*NF) = 22
            flst_real(4, i+NF+j*NF) = j
            flst_real(5, i+NF+j*NF) = -j
c     --------- qbar q -> gamma q' qbar' (incl. qbar q -> gamma q qbar)
            flst_real(1, i+(NF+1)*NF+j*NF) = -i
            flst_real(2, i+(NF+1)*NF+j*NF) = i
            flst_real(3, i+(NF+1)*NF+j*NF) = 22
            flst_real(4, i+(NF+1)*NF+j*NF) = j
            flst_real(5, i+(NF+1)*NF+j*NF) = -j
c     q q' scattering
c     --------- q q' -> gamma q q' (incl. q q -> gamma q q)
            flst_real(1, i+(2*NF+1)*NF+j*NF) = i
            flst_real(2, i+(2*NF+1)*NF+j*NF) = j
            flst_real(3, i+(2*NF+1)*NF+j*NF) = 22
            flst_real(4, i+(2*NF+1)*NF+j*NF) = i
            flst_real(5, i+(2*NF+1)*NF+j*NF) = j
c     --------- qbar qbar' -> gamma qbar qbar' (incl. qb qb -> gamma qb qb)
            flst_real(1, i+(3*NF+1)*NF+j*NF) = -i
            flst_real(2, i+(3*NF+1)*NF+j*NF) = -j
            flst_real(3, i+(3*NF+1)*NF+j*NF) = 22
            flst_real(4, i+(3*NF+1)*NF+j*NF) = -i
            flst_real(5, i+(3*NF+1)*NF+j*NF) = -j
            
            if (i.NE.j) then
               k = k + 1
c     --------- q qbar' -> gamma q qbar'
               flst_real(1, (4*NF+2)*NF+k) = i
               flst_real(2, (4*NF+2)*NF+k) = -j
               flst_real(3, (4*NF+2)*NF+k) = 22
               flst_real(4, (4*NF+2)*NF+k) = i
               flst_real(5, (4*NF+2)*NF+k) = -j
c     --------- qbar q' -> gamma q' qbar
               flst_real(1, (5*NF+1)*NF+k) = -i
               flst_real(2, (5*NF+1)*NF+k) = j
               flst_real(3, (5*NF+1)*NF+k) = 22
               flst_real(4, (5*NF+1)*NF+k) = j
               flst_real(5, (5*NF+1)*NF+k) = -i
            endif
         enddo

c     q g scattering
c     --------- q g -> gamma q g
               flst_real(1, 6*NF**2+i) = i
               flst_real(2, 6*NF**2+i) = 0
               flst_real(3, 6*NF**2+i) = 22
               flst_real(4, 6*NF**2+i) = i
               flst_real(5, 6*NF**2+i) = 0
c     --------- g q -> gamma q g
               flst_real(1, (6*NF+1)*NF+i) = 0
               flst_real(2, (6*NF+1)*NF+i) = i
               flst_real(3, (6*NF+1)*NF+i) = 22
               flst_real(4, (6*NF+1)*NF+i) = i
               flst_real(5, (6*NF+1)*NF+i) = 0
c     --------- qbar g -> gamma qbar g
               flst_real(1, (6*NF+2)*NF+i) = -i
               flst_real(2, (6*NF+2)*NF+i) = 0
               flst_real(3, (6*NF+2)*NF+i) = 22
               flst_real(4, (6*NF+2)*NF+i) = -i
               flst_real(5, (6*NF+2)*NF+i) = 0
c     --------- g qbar -> gamma qbar g
               flst_real(1, (6*NF+3)*NF+i) = 0
               flst_real(2, (6*NF+3)*NF+i) = -i
               flst_real(3, (6*NF+3)*NF+i) = 22
               flst_real(4, (6*NF+3)*NF+i) = -i
               flst_real(5, (6*NF+3)*NF+i) = 0

c     gluon fusion
c     --------- g g -> gamma q qbar
               flst_real(1, (6*NF+4)*NF+i) = 0
               flst_real(2, (6*NF+4)*NF+i) = 0
               flst_real(3, (6*NF+4)*NF+i) = 22
               flst_real(4, (6*NF+4)*NF+i) = i
               flst_real(5, (6*NF+4)*NF+i) = -i
      enddo


c     Do we want to produce NLO distributions prior to event generation?
      if (powheginput('#testplots').eq.1.d0) then
         flg_analysisextrainfo = .TRUE.
         write(*,*) 'testplots: Using analysisextrainfo.'
      endif

c     Do we want to include electromagnetic soft-virtual contributions?
      if (powheginput('#emvirtual').eq.1.d0) then
         flg_with_em = .TRUE.
         write(*,*) 'emvirtual: Computing virtual QED corrections.'
      endif

c     Are we doing a run with enhanced QED radiation?
      sudarwgtfac = powheginput('#enhancedradfac')
      if (sudarwgtfac.gt.0.d0) then
         flg_sudakov_rwgt = .TRUE.
      else
         flg_sudakov_rwgt = .FALSE.
         sudarwgtfac = 1.d0 ! just to be safe...
      endif

      return
      end
