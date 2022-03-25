      subroutine init_couplings
      implicit none
      include 'pwhg_em.h'
      include 'PhysPars.h'
      include 'pwhg_st.h'
      include 'pwhg_math.h'
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_rad.h'
      logical verbose
      parameter(verbose=.true.)
      integer aemrun
      real *8 powheginput
      external powheginput
      real *8 alfaem,pwhg_alphas
      external alfaem,pwhg_alphas
      real *8 alphaem_inv
      common/calphaem_inv/alphaem_inv

c     number of light flavors
      st_nlight = 5

c     alphaem
c     typical values for alphaem:
c     Thompson value:    1/137.0359895d0
c     at z mass (91.188) 1/127.934 (?)
c     at top mass (175)  1/127.011989

      aemrun=0
c     definition of alphaem_pow value, according to aemrun
      if(aemrun.eq.0) then
         alphaem_inv=powheginput('#alphaem_inv')
         if(alphaem_inv.lt.0) alphaem_inv=127.011989
         em_alpha=1d0/alphaem_inv
         zmass_pow=91.188d0     !Not relevant in POWHEG; needed only by set_madgraph_parameters
      elseif(aemrun.eq.1) then
         write(*,*) 'Invalid option for aemrun: program stops'
         call exit(1)
      else
         write(*,*) 'Error while setting aemrun'
         call exit(1)
      endif


      if(verbose) then
         write(*,*) '--------------------------------------'
         write(*,*) 'POWHEG: RELEVANT PARAMETERS'
         write(*,*) '1/em_alpha      ',1.d0/em_alpha
         write(*,*) 'lambda_QCD     ',st_lambda5MSB
         write(*,'(1X,A,f7.3,A,f15.7)') 'alpha_s(',91.2d0,')'
     $,pwhg_alphas(91.2d0**2,st_lambda5MSB,st_nlight)
         write(*,*) '--------------------------------------'
      endif

      end

c-------------------------------------------------------------------------
      function alfaem(q2)
c Alpha_em(MSbar) at the scale q2 = q^2. 
c Uses alpha_Thomson below the electron mass, alpha(mass) below
c mu_mass and m_tau, and the evolution equation above m_tau, comnsidering the b threshold
c This function is taken from the MC@NLO and modified by SA&ER
c-------------------------------------------------------------------------
      implicit none
      include 'pwhg_math.h'
      include 'PhysPars.h'
      integer npoints,ideg
      parameter (npoints=3,ideg=3)
      real*8 ooa(npoints),xlogmu(npoints)
c 1/alpha_em at m_e=0.000511,m_mu=0.1056,m_tau=1.777      
      data ooa     / 137.036, 135.95, 133.513 /
c logs of sqrt(q2) at m_e=0.000511,m_mu=0.1056,m_tau=1.777      
      data xlogmu  / -7.57914, -2.2481, 0.574927 /
      real *8 zm
      real*8 ooaz,xlq,b,q2
      real *8 alfaem

      real *8 alphaem_inv
      common/calphaem_inv/alphaem_inv

      zm=zmass_pow
      ooaz=alphaem_inv


      if(q2.lt.exp(2.*xlogmu(1))) then
         alfaem = 1.d0/ooa(1)	 
      elseif(q2.lt.exp(2.*xlogmu(2))) then
         xlq = log(q2)/2.d0
         alfaem = 1.d0/ooa(2)
      elseif(q2.lt.exp(2.*xlogmu(3))) then
         xlq = log(q2)/2.d0
         alfaem = 1.d0/ooa(3)
      elseif(q2.lt.5.**2) then
         b = 3 + 2*nc*(1d0/3d0)**2 + 2*nc*(2d0/3d0)**2
         xlq = log(q2) - 2.*xlogmu(3)
         alfaem = 1d0/ooa(3)/(1.d0 - 1.d0/3.d0/pi/ooa(3)*b*xlq)
      else
         b = 3 + 3*nc*(1d0/3d0)**2 + 2*nc*(2d0/3d0)**2
         xlq = log(q2/zm**2)
         alfaem = 1d0/ooaz/(1.d0 - 1.d0/3.d0/pi/ooaz*b*xlq)
      endif
      return
      end

