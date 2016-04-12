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
     1     powheginput,dotp,amp,ratio,hnew_damp
      logical ini
      data ini/.true./
      logical angcorr_damp,new_damp
      integer numopt
      save ini,h,angcorr_damp,new_damp,hnew_damp
      external powheginput,dotp
      if(ini) then
         angcorr_damp = powheginput("#angcorr_damp") .eq. 1
         new_damp = powheginput("#new_damp") .eq. 1
         numopt = 0
         if(angcorr_damp) numopt = numopt + 1
         if(new_damp) numopt = numopt + 1
         if(numopt .gt. 1) then
            write(*,*) ' bornzerodamp:'
            write(*,*) ' you should specify only one of'//
     1           'angcorr_damp, new_damp'
            write(*,*) ' exiting ...'
            call exit(-1)
         endif
         if(angcorr_damp) then
            write(*,*) ' using angular correlations aware damp function'
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
            call ampzj(kn_cmpreal,flst_alr(:,alr),amp,dampfac)
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




      subroutine ampZj(p,flav,amp,ratio)
      implicit none
      include 'pwhg_math.h'
      include 'PhysPars.h'
      include 'pwhg_st.h'
      real * 8 p(0:3,5),amp,ratio
      integer flav(5)
      real *8 pCS(0:3,5),q(0:3),qCS(0:3),pp(0:3,5)
      real *8 p12,q2,qz,qt,q0,E,t0,colfac,couplz
c      real * 8 p14,p15,p24,p25,p45
      real *8 phi,sth,cth,sph,cph,c2ph,opcth2,a0th,prop_int
      real *8 s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15
      real *8 T3L,T3Q,chargeL,chargeQ,VL,AL,VH,AH,chhad,chlep,propZ
      real * 8 sphi,cphi,norm,globfac,num,den,ampij
      real * 8 a00,a0,a1,a2,a3,a4,a4born
      real * 8 b00,b0,b1,b2,b3,b4
      integer i
      real*8 dotp
      external dotp
c      real * 8 SGN

      ph_cthw = sqrt(1-ph_sthw**2)
      couplz = 1/(2*ph_sthw*ph_cthw)
c     the fifth particle must have transverse momentum aligned along the -x axis
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

      if (mod(abs(flav(3)),2).eq.1) then
         chargeL = -1
         T3L = -1d0/2d0
      elseif (mod(abs(flav(3)),2).eq.0) then
         chargeL = 0
         T3L = 1d0/2d0
      endif

      if (flav(5).eq.0) then
         if (mod(abs(flav(1)),2).eq.0) then
            chargeQ = 2d0/3d0
            T3Q = 1d0/2d0
         elseif (mod(abs(flav(1)),2).eq.1) then
            chargeQ = -1d0/3d0
            T3Q = -1d0/2d0
         endif
      else
         if (mod(abs(flav(5)),2).eq.0) then
            chargeQ = 2d0/3d0
            T3Q = 1d0/2d0
         elseif (mod(abs(flav(5)),2).eq.1) then
            chargeQ = -1d0/3d0
            T3Q = -1d0/2d0
         endif
      endif
      VL = T3L - 2*chargeL*ph_sthw**2
      AL = -T3L
      VH = T3Q - 2*chargeQ*ph_sthw**2
      AH = -T3Q
      chhad = -chargeQ
      chlep = -chargeL

      propZ = 1/((q2-ph_Zmass**2)**2 + ph_ZmZw**2)*couplz**4
      prop_int = 1/q2 * (q2-ph_Zmass2)/((q2-ph_Zmass2)**2+ph_ZmZw**2)*
     $     couplz**2

      if(flav(1).gt.0.and.flav(2).lt.0) then
      a4born = 1
      a0 = qt**2/(q2+qt**2)
      a1 = -0.2D1*qt*E*(E**2+q2)*qz*sqrt(q2)/(q2+qt**2)/(-q2**2+E**2*(-E
     #**2+2*qt**2))
      a2 = qt**2/(q2+qt**2)
      a3 = -8*qt*E*q2**2*AL*AH*qz*(E**2+q2)*(2*propZ*VL*VH+prop_int*chle
     #p*chhad)/sqrt(q2+qt**2)/(-q2**2+E**2*(-E**2+2*qt**2))/((VL**2+AL**
     #2)*(VH**2+AH**2)*q2**2*propZ+chlep**2*chhad**2+2*q2**2*chhad*chlep
     #*prop_int*VL*VH)
      a4 = 4*sqrt(q2)**5*AL*AH*(2*propZ*VL*VH+prop_int*chlep*chhad)/sqrt
     #(q2+qt**2)/((VL**2+AL**2)*(VH**2+AH**2)*q2**2*propZ+chlep**2*chhad
     #**2+2*q2**2*chhad*chlep*prop_int*VL*VH)
      a00 = 1
      norm = (2*(VL**2+AL**2)*(VH**2+AH**2)*q2**3-2*E**2*(-E**2+2*qt**2)
     #*(VL**2+AL**2)*(VH**2+AH**2)*q2)*propZ+4*prop_int*chlep*chhad*q2**
     #3*VL*VH+(-4*E**2*(-E**2+2*qt**2)*chhad*chlep*VL*VH*prop_int+2*chle
     #p**2*chhad**2)*q2-2*E**2*(-E**2+2*qt**2)*chlep**2*chhad**2/q2
      globfac = -4/(E*q0+E*qz-q2)/E/(-E+q0+qz)
      elseif(flav(1).lt.0.and.flav(2).gt.0) then
      a4born = -1
      a0 = qt**2/(q2+qt**2)
      a1 = -0.2D1*qt*E*(E**2+q2)*qz*sqrt(q2)/(q2+qt**2)/(-q2**2+E**2*(-E
     #**2+2*qt**2))
      a2 = qt**2/(q2+qt**2)
      a3 = 8*qt*E*q2**2*AL*AH*qz*(E**2+q2)*(2*propZ*VL*VH+prop_int*chlep
     #*chhad)/sqrt(q2+qt**2)/(-q2**2+E**2*(-E**2+2*qt**2))/((VL**2+AL**2
     #)*(VH**2+AH**2)*q2**2*propZ+chlep**2*chhad**2+2*q2**2*chhad*chlep*
     #prop_int*VL*VH)
      a4 = -4*sqrt(q2)**5*AL*AH*(2*propZ*VL*VH+prop_int*chlep*chhad)/sqr
     #t(q2+qt**2)/((VL**2+AL**2)*(VH**2+AH**2)*q2**2*propZ+chlep**2*chha
     #d**2+2*q2**2*chhad*chlep*prop_int*VL*VH)
      a00 = 1
      norm = (2*(VL**2+AL**2)*(VH**2+AH**2)*q2**3-2*E**2*(-E**2+2*qt**2)
     #*(VL**2+AL**2)*(VH**2+AH**2)*q2)*propZ+4*prop_int*chlep*chhad*q2**
     #3*VL*VH+(-4*E**2*(-E**2+2*qt**2)*chhad*chlep*VL*VH*prop_int+2*chle
     #p**2*chhad**2)*q2-2*E**2*(-E**2+2*qt**2)*chlep**2*chhad**2/q2
      globfac = -4/(E*q0+E*qz-q2)/E/(-E+q0+qz)
      elseif(flav(1).eq.0.and.flav(2).lt.0) then
      a4born = 1
      a0 = qt**2*(3*q2**2+2*E*(2*E+qz)*q2+E**2*(3*E**2+2*E*qz-2*qt**2))/
     #(q2+qt**2)/(3*q2**2+2*E*(-2*E+qz)*q2+E**2*(3*E**2+2*E*qz-2*qt**2))
      a1 = 0.1D1*qt*(q2**2+6*q2*E*qz+E**2*(E**2+6*E*qz-2*qt**2))*sqrt(q2
     #)/(q2+qt**2)/(3*q2**2+2*E*(-2*E+qz)*q2+E**2*(3*E**2+2*E*qz-2*qt**2
     #))
      a2 = qt**2*(3*q2**2+2*E*(2*E+qz)*q2+E**2*(3*E**2+2*E*qz-2*qt**2))/
     #(q2+qt**2)/(3*q2**2+2*E*(-2*E+qz)*q2+E**2*(3*E**2+2*E*qz-2*qt**2))
      a3 = 4*qt*q2**2*AL*AH*(3*q2**2+2*q2*E*qz+E**2*(-E**2+2*E*qz-2*qt**
     #2))*(2*propZ*VL*VH+prop_int*chlep*chhad)/sqrt(q2+qt**2)/(3*q2**2+2
     #*E*(-2*E+qz)*q2+E**2*(3*E**2+2*E*qz-2*qt**2))/((VL**2+AL**2)*(VH**
     #2+AH**2)*q2**2*propZ+chlep**2*chhad**2+2*q2**2*chhad*chlep*prop_in
     #t*VL*VH)
      a4 = 4*sqrt(q2)**5*AL*AH*(q2**2+6*q2*E*qz-E**2*(-E**2+2*E*qz+2*qt*
     #*2))*(2*propZ*VL*VH+prop_int*chlep*chhad)/sqrt(q2+qt**2)/(3*q2**2+
     #2*E*(-2*E+qz)*q2+E**2*(3*E**2+2*E*qz-2*qt**2))/((VL**2+AL**2)*(VH*
     #*2+AH**2)*q2**2*propZ+chlep**2*chhad**2+2*q2**2*chhad*chlep*prop_i
     #nt*VL*VH)
      a00 = 1
      norm = (-6*(VL**2+AL**2)*(VH**2+AH**2)/E**2*q2**3-4*(-2*E+qz)*(VL*
     #*2+AL**2)*(VH**2+AH**2)/E*q2**2-2*(3*E**2+2*E*qz-2*qt**2)*(VL**2+A
     #L**2)*(VH**2+AH**2)*q2)*propZ-12*chhad*chlep*prop_int*VL*VH/E**2*q
     #2**3-8/E*(-2*E+qz)*chhad*chlep*VL*VH*prop_int*q2**2+(-4*(3*E**2+2*
     #E*qz-2*qt**2)*chhad*chlep*VL*VH*prop_int-6*chlep**2*chhad**2/E**2)
     #*q2-4/E*(-2*E+qz)*chlep**2*chhad**2-2*(3*E**2+2*E*qz-2*qt**2)*chle
     #p**2*chhad**2/q2
      globfac = 2/(E*q0+E*qz-q2)
      elseif(flav(1).gt.0.and.flav(2).eq.0) then
      a4born = 1
      a0 = qt**2*(-3*q2**2+2*E*(-2*E+qz)*q2+E**2*(-3*E**2+2*E*qz+2*qt**2
     #))/(q2+qt**2)/(-3*q2**2+2*E*(2*E+qz)*q2+E**2*(-3*E**2+2*E*qz+2*qt*
     #*2))
      a1 = -0.1D1*qt*(-q2**2+6*q2*E*qz+E**2*(-E**2+6*E*qz+2*qt**2))*sqrt
     #(q2)/(q2+qt**2)/(-3*q2**2+2*E*(2*E+qz)*q2+E**2*(-3*E**2+2*E*qz+2*q
     #t**2))
      a2 = qt**2*(-3*q2**2+2*E*(-2*E+qz)*q2+E**2*(-3*E**2+2*E*qz+2*qt**2
     #))/(q2+qt**2)/(-3*q2**2+2*E*(2*E+qz)*q2+E**2*(-3*E**2+2*E*qz+2*qt*
     #*2))
      a3 = -4*qt*q2**2*AL*AH*(-3*q2**2+2*q2*E*qz+E**2*(E**2+2*E*qz+2*qt*
     #*2))*(2*propZ*VL*VH+prop_int*chlep*chhad)/sqrt(q2+qt**2)/(-3*q2**2
     #+2*E*(2*E+qz)*q2+E**2*(-3*E**2+2*E*qz+2*qt**2))/((VL**2+AL**2)*(VH
     #**2+AH**2)*q2**2*propZ+chlep**2*chhad**2+2*q2**2*chhad*chlep*prop_
     #int*VL*VH)
      a4 = 4*sqrt(q2)**5*AL*AH*(-q2**2+6*q2*E*qz-E**2*(E**2+2*E*qz-2*qt*
     #*2))*(2*propZ*VL*VH+prop_int*chlep*chhad)/sqrt(q2+qt**2)/(-3*q2**2
     #+2*E*(2*E+qz)*q2+E**2*(-3*E**2+2*E*qz+2*qt**2))/((VL**2+AL**2)*(VH
     #**2+AH**2)*q2**2*propZ+chlep**2*chhad**2+2*q2**2*chhad*chlep*prop_
     #int*VL*VH)
      a00 = 1
      norm = (-6*(VL**2+AL**2)*(VH**2+AH**2)/E**2*q2**3+4*(2*E+qz)*(VL**
     #2+AL**2)*(VH**2+AH**2)/E*q2**2+2*(-3*E**2+2*E*qz+2*qt**2)*(VL**2+A
     #L**2)*(VH**2+AH**2)*q2)*propZ-12*chhad*chlep*prop_int*VL*VH/E**2*q
     #2**3+8/E*(2*E+qz)*chhad*chlep*VL*VH*prop_int*q2**2+(4*(-3*E**2+2*E
     #*qz+2*qt**2)*chhad*chlep*VL*VH*prop_int-6*chlep**2*chhad**2/E**2)*
     #q2+4/E*(2*E+qz)*chlep**2*chhad**2+2*(-3*E**2+2*E*qz+2*qt**2)*chlep
     #**2*chhad**2/q2
      globfac = -2/E/(-E+q0+qz)
      elseif(flav(1).eq.0.and.flav(2).gt.0) then
      a4born = -1
      a0 = qt**2*(3*q2**2+2*E*(2*E+qz)*q2+E**2*(3*E**2+2*E*qz-2*qt**2))/
     #(q2+qt**2)/(3*q2**2+2*E*(-2*E+qz)*q2+E**2*(3*E**2+2*E*qz-2*qt**2))
      a1 = 0.1D1*qt*(q2**2+6*q2*E*qz+E**2*(E**2+6*E*qz-2*qt**2))*sqrt(q2
     #)/(q2+qt**2)/(3*q2**2+2*E*(-2*E+qz)*q2+E**2*(3*E**2+2*E*qz-2*qt**2
     #))
      a2 = qt**2*(3*q2**2+2*E*(2*E+qz)*q2+E**2*(3*E**2+2*E*qz-2*qt**2))/
     #(q2+qt**2)/(3*q2**2+2*E*(-2*E+qz)*q2+E**2*(3*E**2+2*E*qz-2*qt**2))
      a3 = -4*qt*q2**2*AL*AH*(3*q2**2+2*q2*E*qz+E**2*(-E**2+2*E*qz-2*qt*
     #*2))*(2*propZ*VL*VH+prop_int*chlep*chhad)/sqrt(q2+qt**2)/(3*q2**2+
     #2*E*(-2*E+qz)*q2+E**2*(3*E**2+2*E*qz-2*qt**2))/((VL**2+AL**2)*(VH*
     #*2+AH**2)*q2**2*propZ+chlep**2*chhad**2+2*q2**2*chhad*chlep*prop_i
     #nt*VL*VH)
      a4 = -4*sqrt(q2)**5*AL*AH*(q2**2+6*q2*E*qz-E**2*(-E**2+2*E*qz+2*qt
     #**2))*(2*propZ*VL*VH+prop_int*chlep*chhad)/sqrt(q2+qt**2)/(3*q2**2
     #+2*E*(-2*E+qz)*q2+E**2*(3*E**2+2*E*qz-2*qt**2))/((VL**2+AL**2)*(VH
     #**2+AH**2)*q2**2*propZ+chlep**2*chhad**2+2*q2**2*chhad*chlep*prop_
     #int*VL*VH)
      a00 = 1
      norm = (-6*(VL**2+AL**2)*(VH**2+AH**2)/E**2*q2**3-4*(-2*E+qz)*(VL*
     #*2+AL**2)*(VH**2+AH**2)/E*q2**2-2*(3*E**2+2*E*qz-2*qt**2)*(VL**2+A
     #L**2)*(VH**2+AH**2)*q2)*propZ-12*chhad*chlep*prop_int*VL*VH/E**2*q
     #2**3-8/E*(-2*E+qz)*chhad*chlep*VL*VH*prop_int*q2**2+(-4*(3*E**2+2*
     #E*qz-2*qt**2)*chhad*chlep*VL*VH*prop_int-6*chlep**2*chhad**2/E**2)
     #*q2-4/E*(-2*E+qz)*chlep**2*chhad**2-2*(3*E**2+2*E*qz-2*qt**2)*chle
     #p**2*chhad**2/q2
      globfac = 2/(E*q0+E*qz-q2)
      elseif(flav(1).lt.0.and.flav(2).eq.0) then
      a4born = -1
      a0 = qt**2*(-3*q2**2+2*E*(-2*E+qz)*q2+E**2*(-3*E**2+2*E*qz+2*qt**2
     #))/(q2+qt**2)/(-3*q2**2+2*E*(2*E+qz)*q2+E**2*(-3*E**2+2*E*qz+2*qt*
     #*2))
      a1 = -0.1D1*qt*(-q2**2+6*q2*E*qz+E**2*(-E**2+6*E*qz+2*qt**2))*sqrt
     #(q2)/(q2+qt**2)/(-3*q2**2+2*E*(2*E+qz)*q2+E**2*(-3*E**2+2*E*qz+2*q
     #t**2))
      a2 = qt**2*(-3*q2**2+2*E*(-2*E+qz)*q2+E**2*(-3*E**2+2*E*qz+2*qt**2
     #))/(q2+qt**2)/(-3*q2**2+2*E*(2*E+qz)*q2+E**2*(-3*E**2+2*E*qz+2*qt*
     #*2))
      a3 = 4*qt*q2**2*AL*AH*(-3*q2**2+2*q2*E*qz+E**2*(E**2+2*E*qz+2*qt**
     #2))*(2*propZ*VL*VH+prop_int*chlep*chhad)/sqrt(q2+qt**2)/(-3*q2**2+
     #2*E*(2*E+qz)*q2+E**2*(-3*E**2+2*E*qz+2*qt**2))/((VL**2+AL**2)*(VH*
     #*2+AH**2)*q2**2*propZ+chlep**2*chhad**2+2*q2**2*chhad*chlep*prop_i
     #nt*VL*VH)
      a4 = -4*sqrt(q2)**5*AL*AH*(-q2**2+6*q2*E*qz-E**2*(E**2+2*E*qz-2*qt
     #**2))*(2*propZ*VL*VH+prop_int*chlep*chhad)/sqrt(q2+qt**2)/(-3*q2**
     #2+2*E*(2*E+qz)*q2+E**2*(-3*E**2+2*E*qz+2*qt**2))/((VL**2+AL**2)*(V
     #H**2+AH**2)*q2**2*propZ+chlep**2*chhad**2+2*q2**2*chhad*chlep*prop
     #_int*VL*VH)
      a00 = 1
      norm = (-6*(VL**2+AL**2)*(VH**2+AH**2)/E**2*q2**3+4*(2*E+qz)*(VL**
     #2+AL**2)*(VH**2+AH**2)/E*q2**2+2*(-3*E**2+2*E*qz+2*qt**2)*(VL**2+A
     #L**2)*(VH**2+AH**2)*q2)*propZ-12*chhad*chlep*prop_int*VL*VH/E**2*q
     #2**3+8/E*(2*E+qz)*chhad*chlep*VL*VH*prop_int*q2**2+(4*(-3*E**2+2*E
     #*qz+2*qt**2)*chhad*chlep*VL*VH*prop_int-6*chlep**2*chhad**2/E**2)*
     #q2+4/E*(2*E+qz)*chlep**2*chhad**2+2*(-3*E**2+2*E*qz+2*qt**2)*chlep
     #**2*chhad**2/q2
      globfac = -2/E/(-E+q0+qz)
      endif

c Born value of a4 equal to qqb case at qt=0
      a4born = a4born*
     1 4*q2**2*AL*AH*(2*propZ*VL*VH+prop_int*chlep*chhad)
     2 /((VL**2+AL**2)*(VH**2+AH**2)*q2**2*propZ+chlep**2*chhad
     3 **2+2*q2**2*chhad*chlep*prop_int*VL*VH)

      b00 = opcth2      
      b0 = a0th
      b1 = 2*sth*cph*cth
      b2 = 0.5D0*c2ph*sth**2
      b3 = sth*cph
      b4 = cth     
      num = a00*b00 + a4born*b4 + a1*b1 + a3*b3 + a2*b2
      den = a00*b00 + a4*b4 + a0*b0 + a1*b1 + a2*b2 + a3*b3
c      ampij = den*norm*globfac
      ampij = den*norm*globfac
      ratio = num/den

c     COLOUR FACTORS, QUARK-GLUON SWITCH
c     colour factors: CF*n from sum over initial colours, 1/4 from
c     average over initial spins, 1/n from average over quark colours
c     and 1/(n^2-1) from average over gluon colours
      if(flav(5).ne.0) then
         amp=-1d0   !quark-gluon switch
         colfac=CF*nc/4d0/nc/(nc**2-1d0)
      else
         amp=1d0    !no quark-gluon switch
         colfac=CF*nc/4d0/nc**2
      endif

      amp=amp*colfac*(4*pi*st_alpha)*ph_unit_e**4*ampij
      amp=amp/(st_alpha/(2*pi))
      end


      subroutine CollinsSoper_frame(pin,pout)
      implicit none
      real*8 dotp
      external dotp
      real *8 pin(0:3,5),ptemp(0:3,5),pout(0:3,5),q(0:3),vec(3),beta,qt
      
      q = pin(:,3) + pin(:,4)
c     first transverse boost to obtain qz=0
      beta = q(3)/q(0)
      vec(1) =  0
      vec(2) =  0
      vec(3) = -1
      call mboost(5,vec,beta,pin(0,1),ptemp(0,1))
      q = ptemp(:,3) + ptemp(:,4)
c     second longitudinal boost to obtain qt=0
      qt = sqrt(q(1)**2+q(2)**2)
      beta = qt/q(0)
      vec(1) = -1
      vec(2) =  0
      vec(3) =  0
      call mboost(5,vec,beta,ptemp(0,1),pout(0,1))
      return
      end

