c     Initialize the histograms
      subroutine init_hist
      implicit none
      double precision powheginput
      double precision ptmin, ptmax, ptwidth, ymin, ymax, ywidth
      common /histpty/ ptmin, ptmax, ptwidth, ymin, ymax, ywidth

      call inihists

c     Get parameters for pt and y histograms
      ptmin = powheginput("#histptmin")
      ptmax = powheginput("#histptmax")
      ptwidth = powheginput("#histptwidth")
      ymin = powheginput("#histymin")
      ymax = powheginput("#histymax")
      ywidth = powheginput("#histywidth")

      if ((ptmin.LE.-1d6).OR.(ptmax.LE.-1d6).OR.(ptwidth.LE.-1d6)) then
         write(*,*) 'histptmin, -max or -width missing.'
         write(*,*) 'Assuming standard values 0, 7000, 100.'
         call bookupeqbins('Photon pT',100.d0,0.d0,7000.d0)
      else
         call bookupeqbins('Photon pT',ptwidth,ptmin,ptmax)
      endif

      if ((ymin.LE.-1d6).OR.(ymax.LE.-1d6).OR.(ywidth.LE.-1d6)) then
         write(*,*) 'histymin, -max or -width missing.'
         write(*,*) 'Assuming standard values -6, 6, 1.'
         call bookupeqbins('Photon y',1.d0,-6.d0,6.d0)
      else
         call bookupeqbins('Photon y',ywidth,ymin,ymax)
      endif

      end


c     Analysis routine using analysisextrainfo to search for photon
c     final states and produce pt and y histograms for them
      subroutine analysis(dsig0)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_anexinf.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      character * 12 type
      integer i, ypos
      logical isqed, isuqcdborn, fs
      double precision dsig0
      double precision y, pt
      integer qedbornpos(2,maxprocborn), qedalrpos(2,maxalr),
     &     qedemitter(maxalr), uqcdborn(maxprocborn)
      save qedbornpos, qedalrpos, qedemitter, uqcdborn
      logical ini
      data ini/.TRUE./
      save ini
      integer nqedborn, nqedalr, nuqcdborn
      data nqedborn, nqedalr, nuqcdborn/0,0,0/
      save nqedborn, nqedalr, nuqcdborn

c     produce arrays that give index of QED process and position
c     of photon and a list of underlying QCD Borns
      if (ini) then
         do i = 1, flst_nborn
            if (isqed(flst_born(1,i), nlegborn, ypos)) then
               nqedborn = nqedborn + 1 ! count no of QED Borns
               qedbornpos(1,nqedborn) = i ! pos of QED Born in flst_born
               qedbornpos(2,nqedborn) = ypos ! position of photon
            elseif (isuqcdborn(i,fs)) then
c$$$               if (fs) then     ! for now only accept uborn for FS emission
                  nuqcdborn = nuqcdborn + 1
                  uqcdborn(nuqcdborn) = i
c$$$               endif
            endif
         enddo

         qedemitter = -1        ! initialize nonsense
         do i = 1, flst_nalr
            if (isqed(flst_alr(1,i), nlegreal, ypos)) then
               nqedalr = nqedalr + 1 ! count no of QED alrs
               qedalrpos(1,nqedalr) = i ! pos of QED alr in flst_alr
               qedalrpos(2,nqedalr) = ypos ! position of photon
               if (ypos.EQ.5) then ! if photon emitted save position of emitter
                  qedemitter(nqedalr) = flst_emitter(i)
               endif
            endif
         enddo

         if (.NOT.nqedalr.EQ.flst_nalr) then
            write(*,*) 'analysis: nqedalr not equal to flst_nalr.'
            write(*,*) 'exit.'
            call exit(-1)
         endif

         ini = .FALSE.
      endif

c     get photon information depending on type of contribution
      type = anexinf_type
      if (type.eq.'born') then
         do i = 1, nqedborn
            call getypt(kn_pborn(0,qedbornpos(2,i)),y,pt)
            call checkandfill(pt, y, anexinf_sigarr(qedbornpos(1,i)))
         enddo

      elseif ((type.eq.'virt').OR.(type.eq.'colr')) then
         do i = 1, nqedborn     ! virtuals with FS photon
            call getypt(kn_pborn(0,qedbornpos(2,i)),y,pt)
            call checkandfill(pt, y, anexinf_sigarr(qedbornpos(1,i)))
         enddo
         do i = 1, nuqcdborn    ! integrated CTs for QED emission (?)
            call getypt(kn_pborn(0,3),y,pt) ! leg unimportant for pt since 2->2?
            call checkandfill(pt, y, anexinf_sigarr(uqcdborn(i)))
         enddo

      elseif (type.eq.'realct') then
         do i = 1, nqedalr
            if (qedemitter(i).EQ.-1) then ! photon not emitted
               call getypt(kn_preal(0,qedalrpos(2,i)),y,pt)
               call checkandfill(pt, y, anexinf_sigarr(qedalrpos(1,i)))
            else if (qedemitter(i).GE.0) then ! photon emission
               call getypt(kn_preal(0,qedemitter(i)),y,pt)
               call checkandfill(pt, y, anexinf_sigarr(qedalrpos(1,i)))
            endif
         enddo

      elseif (type.eq.'real') then
         do i = 1, nqedalr
            call getypt(kn_preal(0,qedalrpos(2,i)),y,pt)
            call checkandfill(pt, y, anexinf_sigarr(qedalrpos(1,i)))
         enddo

      elseif ((type.eq.'reg').OR.(type.eq.'remn')) 
     &        then
c$$$  no idea what this is... turns up only for event generation
c$$$  with this extra if-case one does not run in the case below

      else
         write(*,*) 'analysis: anexinf_type unknown:', type
         write(*,*) 'exit!'
         call exit(-1)
      endif

      return
      end


      subroutine getypt(p,y,pt)
      implicit none
      double precision p(0:3),y,pt

      y=0.5d0*log((p(0)+p(3))/(p(0)-p(3)))
      pt=sqrt(p(1)**2+p(2)**2)

      end


c     Check if a photon is a final state of the flavour structure
c     flst(nfl) and return its position
      logical function isqed(flst, nfl, ypos)

      implicit none
      include 'nlegborn.h'
      integer flst(nlegreal), nfl, ypos
      integer i

      isqed = .FALSE.
      do i = 3, nfl
         if (flst(i).eq.22) then
            isqed = .TRUE.
            ypos = i
            return
         endif
      enddo

      return
      end


c     Check if Born flst_born(1,iborn) is underlying Born for QED
c     radiation.
c     Makes no sense to save emitter since one uborn corresponds to many
c     alrs with different emitters...
      logical function isuqcdborn(iborn,fs)

      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      integer iborn
      logical fs
      integer i, jalr

      isuqcdborn = .FALSE.
      fs = .FALSE.
      do i = 1, flst_born2alr(0,iborn)
         jalr = flst_born2alr(i,iborn) ! position of alr with underlying iborn
         if (flst_alr(nlegreal,jalr).EQ.22) then ! QED radiation
            isuqcdborn = .TRUE.
            if (flst_emitter(jalr).GE.3) then ! check FS emission
               fs = .TRUE.
            endif
            return
         endif
      enddo

      return
      end

      


c     Check limits and fill histogram accordingly
      subroutine checkandfill(pt,y,dsig)

      implicit none
      double precision pt, y,dsig
      double precision ptmin, ptmax, ptwidth, ymin, ymax, ywidth
      common /histpty/ ptmin, ptmax, ptwidth, ymin, ymax, ywidth

      if ((pt.GE.ptmin).AND.(pt.LE.ptmax)
     &     .AND.(y.GE.ymin).AND.(y.LE.ymax)) then
         call filld('Photon pT', pt, dsig)
         call filld('Photon y', y, dsig)
      endif

      return
      end
