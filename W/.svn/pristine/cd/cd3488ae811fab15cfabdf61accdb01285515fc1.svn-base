      subroutine bornzerodamp(alr,r0,rc,rs,rcs,dampfac)
c given the R_alpha region (i.e. the alr) and the associated
c real contribution r (without pdf factor),
c returns in dampfac the damping factor to be applied to
c the real contribution to implement Born zero suppression
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_flg.h'
      include 'pwhg_math.h'
      integer alr
      real * 8 r0,rc,rs,rcs,rapp,dampfac,h,oh,m2,pt2,
     1     powheginput,dotp,amp,ratio,hnew_damp,omcth
      logical ini
      data ini/.true./
      logical angcorr_damp,new_damp,theta_damp
      integer numopt
      save ini,h,angcorr_damp,new_damp,hnew_damp,theta_damp
      external powheginput,dotp
      if(ini) then
         angcorr_damp = powheginput("#angcorr_damp") .eq. 1
         new_damp = powheginput("#new_damp") .eq. 1
         theta_damp = powheginput("#theta_damp") .eq. 1
         numopt = 0         
         if(angcorr_damp) numopt = numopt + 1
         if(new_damp) numopt = numopt + 1
         if(theta_damp) numopt = numopt + 1
         if(numopt.gt.1) then
            write(*,*) ' bornzerodamp:'
            write(*,*) ' you should specify only one of'//
     1           'angcorr_damp, theta_damp, new_damp'
            write(*,*) ' exiting ...'
            call exit(-1)
         endif
         if(angcorr_damp) then
            write(*,*) ' using angular correlations aware damp function'
         elseif(theta_damp) then
            write(*,*) ' using theta dependent damp function'
         elseif(new_damp) then
            write(*,*) ' using new, better default damp function'            
c if less than 0 will not be used
            hnew_damp =  powheginput("#hnew_damp")
            if(hnew_damp.gt.0) then
               write(*,*) ' Using hnew_damp=',hnew_damp
            endif
         endif
         h=powheginput("#hdamp")
         if(h.lt.0) then
            h=powheginput("#hfact")
            if(h.gt.0) then
               write(*,*)'***************************************'
               write(*,*)' Warning: hfact is here for backward'
               write(*,*)' compatibility with older inplementations'
               write(*,*)' New implementations will use hdamp and'
               write(*,*)' bornzerodamp independently.'
               write(*,*)'***************************************'
            endif
         endif
         if(h.gt.0) then
            write(*,*)'***************************************'
            write(*,*)' Using a damping factor h**2/(pt2+h**2)'
            write(*,*)' to separate real contributions between'
            write(*,*)' Sudakov and remnants    h=',h,' GeV   ' 
            write(*,*)'***************************************'
         endif
          ini=.false.
      endif

      if(flg_bornzerodamp) then
         if(angcorr_damp) then
            call ampwj(kn_cmpreal,flst_alr(:,alr),amp,dampfac)
            if(dampfac.lt.0) then
               dampfac = 0
            endif
         elseif(new_damp) then
            rapp = rc+rs-rcs
            if(hnew_damp.gt.0) then
               pt2 = kn_cmpreal(1,5)**2+kn_cmpreal(2,5)**2
               m2 = 2*dotp(kn_cmpreal(:,3),kn_cmpreal(:,4))
               rapp = rapp*hnew_damp**2*m2/(pt2+m2*hnew_damp**2)
            endif
            dampfac= min(1d0,rapp/r0)
            dampfac = max(dampfac,0d0)
         elseif(theta_damp) then
            if(flst_uborn(1,alr)*flst_uborn(3,alr).gt.0) then
c The incoming parton has the same sign as the outgoing lepton
               omcth = dotp(kn_cmpborn(:,1),kn_cmpborn(:,4))
            else
               omcth = dotp(kn_cmpborn(:,1),kn_cmpborn(:,3))
            endif
            dampfac=omcth/(pt2+omcth)
         elseif(r0.gt.5*rc.and.r0.gt.5*rs) then
            dampfac=0
c we might as well return here, dampfac already zero ...
            return
         else
            dampfac = 1
         endif
      else
         dampfac = 1
      endif

c local variables
      if(h.gt.0) then
         pt2 = kn_cmpreal(1,5)**2+kn_cmpreal(2,5)**2
         m2 = 2*dotp(kn_cmpreal(:,3),kn_cmpreal(:,4))
         dampfac = dampfac*h**2*m2/(pt2+m2*h**2)
      endif

      if(dampfac.gt.1) dampfac = 1

      end




      subroutine ampWj(p,flav,amp,ratio)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_math.h'
      include 'PhysPars.h'
      include 'pwhg_st.h'
      real * 8 p(0:3,5),amp,ratio
      integer flav(5)
      integer i,j     
      real * 8 pCS(0:3,5),q(0:3),qCS(0:3),pp(0:3,5)
c      real * 8 p14,p15,p24,p25,p45
      real * 8 p12,q2,E,qt,qz,q0,gw,colfac
      real * 8 phi,sth,cth,sph,cph,c2ph,opcth2,a0th,propW
      real * 8 s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,t0
      real * 8 sphi,cphi,norm,globfac,num,den,ampij
      real * 8 a00,a0,a1,a2,a3,a4,a4born
      real * 8 b00,b0,b1,b2,b3,b4
c      real * 8 SGN
      real * 8 dotp
      external dotp

      gw=ph_unit_e/ph_sthw
c     the gluon must have transverse momentum aligned along the -x axis
      pp=p
      phi = pi - atan2(p(2,5),p(1,5))
      cphi = cos(phi)
      sphi = sin(phi)
      do i=3,5
         pp(1,i) = cphi*p(1,i)-sphi*p(2,i)
         pp(2,i) = sphi*p(1,i)+cphi*p(2,i)
      enddo

      q = pp(:,3) + pp(:,4)    
      qz = q(3)
c      SGN = sign(1d0,qz)
      qt = sqrt(q(1)**2+q(2)**2)
      q0 = q(0)

      call CollinsSoper_frame(pp,pCS)

      qCS = pCS(:,3) + pCS(:,4)
      p12=dotp(pCS(:,1),pCS(:,2))
c      p14=dotp(pCS(:,1),pCS(:,3))
c      p15=dotp(pCS(:,1),pCS(:,4))
c      p24=dotp(pCS(:,2),pCS(:,3))
c      p25=dotp(pCS(:,2),pCS(:,4))
c      p45=dotp(pCS(:,3),pCS(:,4))
      q2=dotp(qCS,qCS)
      E = sqrt(2*p12)

      cth = pCS(3,3)/pCS(0,3)
      sth = sqrt(1-cth**2)
      phi = atan2(pCS(2,3),pCS(1,3))
      sph = sin(phi)
      cph = cos(phi)

      c2ph   = cph**2-sph**2
      opcth2 = 1+cth**2
      a0th   = (1-3*cth**2)/2d0

      propW = 1d0/((q2-ph_Wmass**2)**2 + ph_WmWw**2)

      if(flav(3).gt.0) then
      if(flav(1).gt.0.and.flav(2).lt.0) then
      a4born = 1
      a00 = 1
      a0 = qt**2/(q2+qt**2)
      a1 = 0.2D1*qt*sqrt(q2)*E*(q2+E**2)*qz/(q2+qt**2)/(q2**2+E**4-2*E**
     #2*qt**2)
      a2 = qt**2/(q2+qt**2)
      a3 = 4*qt*E*(q2+E**2)*qz/sqrt(q2+qt**2)/(q2**2+E**4-2*E**2*qt**2)
      a4 = 2*sqrt(q2)/sqrt(q2+qt**2)
      norm = 8*(q2**2+E**4-2*E**2*qt**2)*q2*propW
      globfac = -4/(q2-E*q0-E*qz)/E/(-q0-qz+E)
      elseif(flav(1).lt.0.and.flav(2).gt.0) then
      a4born = -1
      a00 = 1
      a0 = qt**2/(q2+qt**2)
      a1 = 0.2D1*qt*sqrt(q2)*E*(q2+E**2)*qz/(q2+qt**2)/(q2**2+E**4-2*E**
     #2*qt**2)
      a2 = qt**2/(q2+qt**2)
      a3 = -4*qt*E*(q2+E**2)*qz/sqrt(q2+qt**2)/(q2**2+E**4-2*E**2*qt**2)
      a4 = -2*sqrt(q2)/sqrt(q2+qt**2)
      norm = 8*(q2**2+E**4-2*E**2*qt**2)*q2*propW
      globfac = -4/(q2-E*q0-E*qz)/E/(-q0-qz+E)
      elseif(flav(1).eq.0.and.flav(2).lt.0) then
      a4born = 1
      a00 = 1
      a0 = qt**2*(3*E**4+3*q2**2+4*E**2*q2+2*q2*E*qz+2*E**3*qz-2*E**2*qt
     #**2)/(q2+qt**2)/(-4*E**2*q2+3*q2**2+3*E**4+2*q2*E*qz+2*E**3*qz-2*E
     #**2*qt**2)
      a1 = 0.1D1*qt*sqrt(q2)*(q2**2+E**4+6*q2*E*qz+6*E**3*qz-2*E**2*qt**
     #2)/(q2+qt**2)/(-4*E**2*q2+3*q2**2+3*E**4+2*q2*E*qz+2*E**3*qz-2*E**
     #2*qt**2)
      a2 = qt**2*(3*E**4+3*q2**2+4*E**2*q2+2*q2*E*qz+2*E**3*qz-2*E**2*qt
     #**2)/(q2+qt**2)/(-4*E**2*q2+3*q2**2+3*E**4+2*q2*E*qz+2*E**3*qz-2*E
     #**2*qt**2)
      a3 = 2*qt*(3*q2**2-E**4+2*q2*E*qz+2*E**3*qz-2*E**2*qt**2)/sqrt(q2+
     #qt**2)/(-4*E**2*q2+3*q2**2+3*E**4+2*q2*E*qz+2*E**3*qz-2*E**2*qt**2
     #)
      a4 = 2*sqrt(q2)*(q2**2+E**4+6*q2*E*qz-2*E**3*qz-2*E**2*qt**2)/sqrt
     #(q2+qt**2)/(-4*E**2*q2+3*q2**2+3*E**4+2*q2*E*qz+2*E**3*qz-2*E**2*q
     #t**2)
      norm = -8*(-4*E**2*q2+3*q2**2+3*E**4+2*q2*E*qz+2*E**3*qz-2*E**2*qt
     #**2)*propW*q2/E**2
      globfac = -2/(q2-E*q0-E*qz)
      elseif(flav(1).gt.0.and.flav(2).eq.0) then
      a4born = 1
      a00 = 1
      a0 = qt**2*(-2*q2*E*qz+4*E**2*q2+3*q2**2+3*E**4-2*E**3*qz-2*E**2*q
     #t**2)/(q2+qt**2)/(3*q2**2-4*E**2*q2-2*q2*E*qz-2*E**3*qz-2*E**2*qt*
     #*2+3*E**4)
      a1 = -0.1D1*qt*sqrt(q2)*(q2**2+E**4-6*q2*E*qz-6*E**3*qz-2*E**2*qt*
     #*2)/(q2+qt**2)/(3*q2**2-4*E**2*q2-2*q2*E*qz-2*E**3*qz-2*E**2*qt**2
     #+3*E**4)
      a2 = qt**2*(-2*q2*E*qz+4*E**2*q2+3*q2**2+3*E**4-2*E**3*qz-2*E**2*q
     #t**2)/(q2+qt**2)/(3*q2**2-4*E**2*q2-2*q2*E*qz-2*E**3*qz-2*E**2*qt*
     #*2+3*E**4)
      a3 = -2*qt*(3*q2**2-2*q2*E*qz-E**4-2*E**3*qz-2*E**2*qt**2)/sqrt(q2
     #+qt**2)/(3*q2**2-4*E**2*q2-2*q2*E*qz-2*E**3*qz-2*E**2*qt**2+3*E**4
     #)
      a4 = 2*sqrt(q2)*(q2**2-6*q2*E*qz+E**4+2*E**3*qz-2*E**2*qt**2)/sqrt
     #(q2+qt**2)/(3*q2**2-4*E**2*q2-2*q2*E*qz-2*E**3*qz-2*E**2*qt**2+3*E
     #**4)
      norm = -8*(3*q2**2-4*E**2*q2-2*q2*E*qz-2*E**3*qz-2*E**2*qt**2+3*E*
     #*4)*q2*propW/E**2
      globfac = 2/E/(-q0-qz+E)
      elseif(flav(1).eq.0.and.flav(2).gt.0) then
      a4born = -1
      a00 = 1
      a0 = qt**2*(3*E**4+3*q2**2+4*E**2*q2+2*q2*E*qz+2*E**3*qz-2*E**2*qt
     #**2)/(q2+qt**2)/(-4*E**2*q2+3*q2**2+3*E**4+2*q2*E*qz+2*E**3*qz-2*E
     #**2*qt**2)
      a1 = 0.1D1*qt*sqrt(q2)*(q2**2+E**4+6*q2*E*qz+6*E**3*qz-2*E**2*qt**
     #2)/(q2+qt**2)/(-4*E**2*q2+3*q2**2+3*E**4+2*q2*E*qz+2*E**3*qz-2*E**
     #2*qt**2)
      a2 = qt**2*(3*E**4+3*q2**2+4*E**2*q2+2*q2*E*qz+2*E**3*qz-2*E**2*qt
     #**2)/(q2+qt**2)/(-4*E**2*q2+3*q2**2+3*E**4+2*q2*E*qz+2*E**3*qz-2*E
     #**2*qt**2)
      a3 = -2*qt*(3*q2**2-E**4+2*q2*E*qz+2*E**3*qz-2*E**2*qt**2)/sqrt(q2
     #+qt**2)/(-4*E**2*q2+3*q2**2+3*E**4+2*q2*E*qz+2*E**3*qz-2*E**2*qt**
     #2)
      a4 = -2*sqrt(q2)*(q2**2+E**4+6*q2*E*qz-2*E**3*qz-2*E**2*qt**2)/sqr
     #t(q2+qt**2)/(-4*E**2*q2+3*q2**2+3*E**4+2*q2*E*qz+2*E**3*qz-2*E**2*
     #qt**2)
      norm = -8*(-4*E**2*q2+3*q2**2+3*E**4+2*q2*E*qz+2*E**3*qz-2*E**2*qt
     #**2)*propW*q2/E**2
      globfac = -2/(q2-E*q0-E*qz)
      elseif(flav(1).lt.0.and.flav(2).eq.0) then
      a4born = -1
      a00 = 1
      a0 = qt**2*(-2*q2*E*qz+4*E**2*q2+3*q2**2+3*E**4-2*E**3*qz-2*E**2*q
     #t**2)/(q2+qt**2)/(3*q2**2-4*E**2*q2-2*q2*E*qz-2*E**3*qz-2*E**2*qt*
     #*2+3*E**4)
      a1 = -0.1D1*qt*sqrt(q2)*(q2**2+E**4-6*q2*E*qz-6*E**3*qz-2*E**2*qt*
     #*2)/(q2+qt**2)/(3*q2**2-4*E**2*q2-2*q2*E*qz-2*E**3*qz-2*E**2*qt**2
     #+3*E**4)
      a2 = qt**2*(-2*q2*E*qz+4*E**2*q2+3*q2**2+3*E**4-2*E**3*qz-2*E**2*q
     #t**2)/(q2+qt**2)/(3*q2**2-4*E**2*q2-2*q2*E*qz-2*E**3*qz-2*E**2*qt*
     #*2+3*E**4)
      a3 = 2*qt*(3*q2**2-2*q2*E*qz-E**4-2*E**3*qz-2*E**2*qt**2)/sqrt(q2+
     #qt**2)/(3*q2**2-4*E**2*q2-2*q2*E*qz-2*E**3*qz-2*E**2*qt**2+3*E**4)
      a4 = -2*sqrt(q2)*(q2**2-6*q2*E*qz+E**4+2*E**3*qz-2*E**2*qt**2)/sqr
     #t(q2+qt**2)/(3*q2**2-4*E**2*q2-2*q2*E*qz-2*E**3*qz-2*E**2*qt**2+3*
     #E**4)
      norm = -8*(3*q2**2-4*E**2*q2-2*q2*E*qz-2*E**3*qz-2*E**2*qt**2+3*E*
     #*4)*q2*propW/E**2
      globfac = 2/E/(-q0-qz+E)
      endif
      else
      if(flav(1).gt.0.and.flav(2).lt.0) then
      a4born = -1
      a00 = 1
      a0 = qt**2/(q2+qt**2)
      a1 = 0.2D1*qt*sqrt(q2)*E*(q2+E**2)*qz/(q2+qt**2)/(q2**2+E**4-2*E**
     #2*qt**2)
      a2 = qt**2/(q2+qt**2)
      a3 = -4*qt*E*(q2+E**2)*qz/sqrt(q2+qt**2)/(q2**2+E**4-2*E**2*qt**2)
      a4 = -2*sqrt(q2)/sqrt(q2+qt**2)
      norm = 8*(q2**2+E**4-2*E**2*qt**2)*q2*propW
      globfac = -4/(q2-E*q0-E*qz)/E/(-q0-qz+E)
      elseif(flav(1).lt.0.and.flav(2).gt.0) then
      a4born = 1
      a00 = 1
      a0 = qt**2/(q2+qt**2)
      a1 = 0.2D1*qt*sqrt(q2)*E*(q2+E**2)*qz/(q2+qt**2)/(q2**2+E**4-2*E**
     #2*qt**2)
      a2 = qt**2/(q2+qt**2)
      a3 = 4*qt*E*(q2+E**2)*qz/sqrt(q2+qt**2)/(q2**2+E**4-2*E**2*qt**2)
      a4 = 2*sqrt(q2)/sqrt(q2+qt**2)
      norm = 8*(q2**2+E**4-2*E**2*qt**2)*q2*propW
      globfac = -4/(q2-E*q0-E*qz)/E/(-q0-qz+E)
      elseif(flav(1).eq.0.and.flav(2).lt.0) then
      a4born = -1
      a00 = 1
      a0 = qt**2*(3*E**4+3*q2**2+4*E**2*q2+2*q2*E*qz+2*E**3*qz-2*E**2*qt
     #**2)/(q2+qt**2)/(-4*E**2*q2+3*q2**2+3*E**4+2*q2*E*qz+2*E**3*qz-2*E
     #**2*qt**2)
      a1 = 0.1D1*qt*sqrt(q2)*(q2**2+E**4+6*q2*E*qz+6*E**3*qz-2*E**2*qt**
     #2)/(q2+qt**2)/(-4*E**2*q2+3*q2**2+3*E**4+2*q2*E*qz+2*E**3*qz-2*E**
     #2*qt**2)
      a2 = qt**2*(3*E**4+3*q2**2+4*E**2*q2+2*q2*E*qz+2*E**3*qz-2*E**2*qt
     #**2)/(q2+qt**2)/(-4*E**2*q2+3*q2**2+3*E**4+2*q2*E*qz+2*E**3*qz-2*E
     #**2*qt**2)
      a3 = -2*qt*(3*q2**2-E**4+2*q2*E*qz+2*E**3*qz-2*E**2*qt**2)/sqrt(q2
     #+qt**2)/(-4*E**2*q2+3*q2**2+3*E**4+2*q2*E*qz+2*E**3*qz-2*E**2*qt**
     #2)
      a4 = -2*sqrt(q2)*(q2**2+E**4+6*q2*E*qz-2*E**3*qz-2*E**2*qt**2)/sqr
     #t(q2+qt**2)/(-4*E**2*q2+3*q2**2+3*E**4+2*q2*E*qz+2*E**3*qz-2*E**2*
     #qt**2)
      norm = -8*(-4*E**2*q2+3*q2**2+3*E**4+2*q2*E*qz+2*E**3*qz-2*E**2*qt
     #**2)*propW*q2/E**2
      globfac = -2/(q2-E*q0-E*qz)
      elseif(flav(1).gt.0.and.flav(2).eq.0) then
      a4born = -1
      a00 = 1
      a0 = qt**2*(-2*q2*E*qz+4*E**2*q2+3*q2**2+3*E**4-2*E**3*qz-2*E**2*q
     #t**2)/(q2+qt**2)/(3*q2**2-4*E**2*q2-2*q2*E*qz-2*E**3*qz-2*E**2*qt*
     #*2+3*E**4)
      a1 = -0.1D1*qt*sqrt(q2)*(q2**2+E**4-6*q2*E*qz-6*E**3*qz-2*E**2*qt*
     #*2)/(q2+qt**2)/(3*q2**2-4*E**2*q2-2*q2*E*qz-2*E**3*qz-2*E**2*qt**2
     #+3*E**4)
      a2 = qt**2*(-2*q2*E*qz+4*E**2*q2+3*q2**2+3*E**4-2*E**3*qz-2*E**2*q
     #t**2)/(q2+qt**2)/(3*q2**2-4*E**2*q2-2*q2*E*qz-2*E**3*qz-2*E**2*qt*
     #*2+3*E**4)
      a3 = 2*qt*(3*q2**2-2*q2*E*qz-E**4-2*E**3*qz-2*E**2*qt**2)/sqrt(q2+
     #qt**2)/(3*q2**2-4*E**2*q2-2*q2*E*qz-2*E**3*qz-2*E**2*qt**2+3*E**4)
      a4 = -2*sqrt(q2)*(q2**2-6*q2*E*qz+E**4+2*E**3*qz-2*E**2*qt**2)/sqr
     #t(q2+qt**2)/(3*q2**2-4*E**2*q2-2*q2*E*qz-2*E**3*qz-2*E**2*qt**2+3*
     #E**4)
      norm = -8*(3*q2**2-4*E**2*q2-2*q2*E*qz-2*E**3*qz-2*E**2*qt**2+3*E*
     #*4)*q2*propW/E**2
      globfac = 2/E/(-q0-qz+E)
      elseif(flav(1).eq.0.and.flav(2).gt.0) then
      a4born = 1
      a00 = 1
      a0 = qt**2*(3*E**4+3*q2**2+4*E**2*q2+2*q2*E*qz+2*E**3*qz-2*E**2*qt
     #**2)/(q2+qt**2)/(-4*E**2*q2+3*q2**2+3*E**4+2*q2*E*qz+2*E**3*qz-2*E
     #**2*qt**2)
      a1 = 0.1D1*qt*sqrt(q2)*(q2**2+E**4+6*q2*E*qz+6*E**3*qz-2*E**2*qt**
     #2)/(q2+qt**2)/(-4*E**2*q2+3*q2**2+3*E**4+2*q2*E*qz+2*E**3*qz-2*E**
     #2*qt**2)
      a2 = qt**2*(3*E**4+3*q2**2+4*E**2*q2+2*q2*E*qz+2*E**3*qz-2*E**2*qt
     #**2)/(q2+qt**2)/(-4*E**2*q2+3*q2**2+3*E**4+2*q2*E*qz+2*E**3*qz-2*E
     #**2*qt**2)
      a3 = 2*qt*(3*q2**2-E**4+2*q2*E*qz+2*E**3*qz-2*E**2*qt**2)/sqrt(q2+
     #qt**2)/(-4*E**2*q2+3*q2**2+3*E**4+2*q2*E*qz+2*E**3*qz-2*E**2*qt**2
     #)
      a4 = 2*sqrt(q2)*(q2**2+E**4+6*q2*E*qz-2*E**3*qz-2*E**2*qt**2)/sqrt
     #(q2+qt**2)/(-4*E**2*q2+3*q2**2+3*E**4+2*q2*E*qz+2*E**3*qz-2*E**2*q
     #t**2)
      norm = -8*(-4*E**2*q2+3*q2**2+3*E**4+2*q2*E*qz+2*E**3*qz-2*E**2*qt
     #**2)*propW*q2/E**2
      globfac = -2/(q2-E*q0-E*qz)
      elseif(flav(1).lt.0.and.flav(2).eq.0) then
      a4born = 1
      a00 = 1
      a0 = qt**2*(-2*q2*E*qz+4*E**2*q2+3*q2**2+3*E**4-2*E**3*qz-2*E**2*q
     #t**2)/(q2+qt**2)/(3*q2**2-4*E**2*q2-2*q2*E*qz-2*E**3*qz-2*E**2*qt*
     #*2+3*E**4)
      a1 = -0.1D1*qt*sqrt(q2)*(q2**2+E**4-6*q2*E*qz-6*E**3*qz-2*E**2*qt*
     #*2)/(q2+qt**2)/(3*q2**2-4*E**2*q2-2*q2*E*qz-2*E**3*qz-2*E**2*qt**2
     #+3*E**4)
      a2 = qt**2*(-2*q2*E*qz+4*E**2*q2+3*q2**2+3*E**4-2*E**3*qz-2*E**2*q
     #t**2)/(q2+qt**2)/(3*q2**2-4*E**2*q2-2*q2*E*qz-2*E**3*qz-2*E**2*qt*
     #*2+3*E**4)
      a3 = -2*qt*(3*q2**2-2*q2*E*qz-E**4-2*E**3*qz-2*E**2*qt**2)/sqrt(q2
     #+qt**2)/(3*q2**2-4*E**2*q2-2*q2*E*qz-2*E**3*qz-2*E**2*qt**2+3*E**4
     #)
      a4 = 2*sqrt(q2)*(q2**2-6*q2*E*qz+E**4+2*E**3*qz-2*E**2*qt**2)/sqrt
     #(q2+qt**2)/(3*q2**2-4*E**2*q2-2*q2*E*qz-2*E**3*qz-2*E**2*qt**2+3*E
     #**4)
      norm = -8*(3*q2**2-4*E**2*q2-2*q2*E*qz-2*E**3*qz-2*E**2*qt**2+3*E*
     #*4)*q2*propW/E**2
      globfac = 2/E/(-q0-qz+E)
      endif
      endif

      b00 = opcth2      
      b0 = a0th
      b1 = 2*sth*cph*cth
      b2 = 0.5D0*c2ph*sth**2
      b3 = sth*cph
      b4 = cth

      a4born = a4born*2

      num = a00*b00 + a4born*b4 + a1*b1 + a2*b2 + a3*b3
      den = a00*b00 + a0*b0 + a1*b1 + a2*b2 + a3*b3 + a4*b4
      ampij = den*norm*globfac
      ratio = num/den

c     COLOUR FACTORS, QUARK-GLUON SWITCH AND CKM MATRIX
c     colour factors: CF*n from sum over initial colours, 1/4 from
c     average over initial spins, 1/n from average over quark colours
c     and 1/(n^2-1) from average over gluon colours
      if(flav(5).ne.0) then
         amp=-1d0   !quark-gluon switch
         colfac=CF*nc/4d0/nc/(nc**2-1d0)
         j=flav(5)
         if (flav(1).eq.0) then
            i=flav(2)
         else
            i=flav(1)
         endif
      else
         amp=1d0   !no quark-gluon switch
         colfac=CF*nc/4d0/nc**2
         i=flav(1)
         j=flav(2)
      endif
      if (mod(i,2).eq.0) then
         amp=amp*ph_CKM(abs(i)/2,(abs(j)+1)/2)**2
      else
         amp=amp*ph_CKM(abs(j)/2,(abs(i)+1)/2)**2
      endif

      amp=amp*colfac*(gw**2/8d0)**2*(4*pi*st_alpha)*ampij
      amp=amp/(st_alpha/(2*pi))
      end


      subroutine CollinsSoper_frame(pin,pout)
      implicit none
      real* 8 dotp
      external dotp
      real * 8 pin(0:3,5),ptemp(0:3,5),pout(0:3,5),q(0:3),vec(3),beta,
     $     qt,qmod      
      q = pin(:,3) + pin(:,4)
c     first transverse boost to obtain qz=0
      beta = -q(3)/q(0)
      vec(1) =  0
      vec(2) =  0
      vec(3) = 1
      call mboost(5,vec,beta,pin(0,1),ptemp(0,1))
      q = ptemp(:,3) + ptemp(:,4)
c     second longitudinal boost to obtain qt=0
      qt = sqrt(q(1)**2+q(2)**2)
      qmod = sqrt(q(1)**2+q(2)**2+q(3)**2)
      beta = -qt/q(0)
      if (qt.ne.0d0) then         
         vec(1) = q(1)/qt
         vec(2) = q(2)/qt
         vec(3) = q(3)/qt
         call mboost(5,vec,beta,ptemp(0,1),pout(0,1))
      else
         pout=ptemp
      endif    
      end


