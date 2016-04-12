c Find the underlying Born momenta from the real momenta and the
c emitter-readiated pair
      subroutine compuborn(em,rad,reslist,cmppborn)
      implicit none
      include 'nlegborn.h'
      integer em,rad,reslist(nlegreal)
      include 'pwhg_flst.h'
      real * 8 cmppborn(0:3,nlegreal)
      if(em.lt.3) then
         call findubisr(rad,cmppborn)
      else
         call findubfsr(em,rad,reslist,cmppborn)
      endif
      end
      
      subroutine findubfsr(em,rad,reslist,cmppborn)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      integer em,rad,reslist(nlegreal)
      real * 8 cmppborn(0:3,nlegreal)
      include 'pwhg_kn.h'
      real * 8 krecv(3),q0,q2,krec,beta,      
     1     k0rec,k,vec(3)
      integer res,j
      res=reslist(em)
      if(res.ne.reslist(rad)) then
         write(*,*) ' findubfsr: error'
         call pwhg_exit(-1)
      endif
      cmppborn=kn_cmpreal
      if(res.ne.0) then
         call boost2reson(kn_cmpreal(:,res),nlegreal,
     1        kn_cmpreal,cmppborn)
         q0=cmppborn(0,res)
      else
         q0=2*cmppborn(0,1)
      endif
      q2=q0**2
c recoil system momentum 
      k0rec=q0-cmppborn(0,em)-cmppborn(0,rad)
      krecv=-cmppborn(1:3,em)-cmppborn(1:3,rad)
      krec=sqrt(krecv(1)**2+krecv(2)**2+krecv(3)**2)
      beta=(q2-(k0rec+krec)**2)/(q2+(k0rec+krec)**2)
      vec=krecv/krec
      if(res.eq.0) then
         call mboost(nlegreal-2,vec,beta,
     1        cmppborn(:,3:nlegreal),cmppborn(:,3:nlegreal))
      else
         do j=3,nlegreal
            if(reslist(j).eq.res) then
               call mboost(1,vec,beta,
     1              cmppborn(:,j:j),cmppborn(:,j:j))
            endif
         enddo
      endif
      k=(q2-(k0rec**2-krec**2))/(2*q0)
      cmppborn(0,em)=k
      cmppborn(1:3,em)=-vec*k
      cmppborn(:,rad)=0
      end

      subroutine findubisr(j,cmppborn)
      implicit none
      integer j
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      real * 8 cmppborn(0:3,nlegreal)
      include 'pwhg_kn.h'
      real * 8 krecv(3),q0,q2,krec,k0rec,
     1     krecperp,mrec2,beta,vec(3)
      cmppborn(0:3,1)=kn_cmpreal(0:3,1)
      cmppborn(0:3,2)=kn_cmpreal(0:3,2)
      q0=2*cmppborn(0,1)
      q2=q0**2
c recoil system momentum 
      k0rec=q0-kn_cmpreal(0,j)
      krecv=-kn_cmpreal(1:3,j)
      krec=sqrt(krecv(1)**2+krecv(2)**2+krecv(3)**2)
      mrec2=(k0rec**2-krec**2)
      beta=-krecv(3)/k0rec
      vec(1)=0
      vec(2)=0
      vec(3)=1
      call mboost(nlegreal-2,vec,beta,
     1     kn_cmpreal(:,3:nlegreal),cmppborn(:,3:nlegreal))
c Now the transverse boost
      krecperp=sqrt(krecv(1)**2+krecv(2)**2)
      vec(3)=0
      vec(1:2)=krecv(1:2)/krecperp
      beta=-krecperp/sqrt(mrec2+krecperp**2)
      call mboost(nlegreal-2,vec,beta,
     1     cmppborn(:,3:nlegreal),cmppborn(:,3:nlegreal))
      cmppborn(:,j)=0
      end

      function dijterm(em,rad,alr)
      implicit none
      real * 8 dijterm
      integer em,rad,alr
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_par.h'
      integer rflav(nlegreal),reslist(nlegreal),res
      real * 8 cmppborn(0:3,nlegreal)
      real * 8 avub,getdistance,dalr,tmp
      integer nub,mergeisr,mergefsr,i,j,k,ifl1,ifl2,onem
      parameter (onem=1000000)
      logical ini,olddij
      data ini/.true./
      save ini,olddij
      real * 8 powheginput
      external powheginput
      if(ini) then
         if(powheginput("#olddij").eq.1) then
            olddij=.true.
         else
            olddij=.false.
         endif
         ini=.false.
      endif
      if(olddij) then
         dijterm=kn_dijterm(em,rad)
         return
      endif
      rflav(:)=flst_alr(:,alr)
c we only consider singularities within the same
c resonance decay (or in production)
      reslist(:)=flst_alrres(:,alr)
      res=reslist(nlegreal)
c find UB flavour
      if(em.eq.0) then
         continue
      elseif(em.eq.1) then
         rflav(1)=mergeisr(rflav(1),rflav(rad))
      elseif(em.eq.2) then
         rflav(2)=mergeisr(rflav(2),rflav(rad))
      else
         rflav(em)=mergefsr(rflav(em),rflav(rad))
      endif
c invalidater rad parton with impossible pdg code
      rflav(rad)=onem
      call compuborn(em,rad,reslist,cmppborn)
c looop over all possible singularities
      dalr=getdistance(em,rad,res,kn_cmpreal)
      if(abs(dalr/kn_dijterm(em,rad)-1).gt.1d-6) then
         write(*,*) 'dalr', dalr/kn_dijterm(em,rad)
         write(*,*) 'This may imply that running the program ',
     1              'with olddij is wrong; it may also imply ',
     2              'that kn_dijterm_soft() are incorrect.',
     3              'A dijterm_soft() function should be built instead'
         call pwhg_exit(-1)
      endif
c get average singularity from underlying Born
      avub=1
      nub=0
      do j=3,nlegreal
         if(reslist(j).ne.res) cycle
         if(res.eq.0) then
            ifl1=mergeisr(rflav(1),rflav(j))
            ifl2=mergeisr(rflav(2),rflav(j))
            if(ifl1.eq.rflav(1).and.ifl2.eq.rflav(2)) then
               tmp = getdistance(0,j,res,cmppborn)
               if(tmp.eq.0) goto 999
               avub=avub*(1+dalr/tmp)
               nub=nub+1
            else
               if(ifl1.ne.onem) then
               tmp = getdistance(1,j,res,cmppborn)
               if(tmp.eq.0) goto 999
                  avub=avub*(1+dalr/tmp)
                  nub=nub+1
               endif
               if(ifl2.ne.onem) then
               tmp = getdistance(2,j,res,cmppborn)
               if(tmp.eq.0) goto 999
                  avub=avub*(1+dalr/tmp)
                  nub=nub+1
               endif
            endif
         endif
         do k=j+1,nlegreal
            if(reslist(k).ne.res) cycle
            if(mergefsr(rflav(j),rflav(k)).ne.onem) then
               tmp = getdistance(j,k,res,cmppborn)
               if(tmp.eq.0) goto 999
               avub=avub*(1+dalr/tmp)
               nub=nub+1
            endif
         enddo
      enddo
      dijterm=dalr*avub
      return
 999  continue
      dijterm = 1d200
      if(dalr.eq.0) then
         call increasecnt("dalr and ub distance = 0 in ubprojections")
      else
         call increasecnt("ub distance = 0 in ubprojections")
      endif
      end

      function getdistance(em,rad,res,cmp)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_par.h'
      include 'pwhg_kn.h'
      real * 8 getdistance
      integer em,rad,res
      real * 8 cmp(0:3,nlegreal),y,e_em,e_rad
      real * 8 dotp
      external dotp
      if(em.lt.3) then
         y=1-dotp(cmp(:,1),cmp(:,rad))/(cmp(0,1)*cmp(0,rad))
         if(em.eq.0) then
            getdistance=(cmp(0,rad)**2*(1-y**2))**par_diexp
         elseif(em.eq.1) then
            getdistance=(cmp(0,rad)**2*2*(1-y))**par_diexp
         elseif(em.eq.2) then
            getdistance=(cmp(0,rad)**2*2*(1+y))**par_diexp
         endif
      else
         if(res.eq.0) then
            e_em=cmp(0,em)
            e_rad=cmp(0,rad)
         else
            e_em=dotp(cmp(:,res),cmp(:,em))
            e_rad=dotp(cmp(:,res),cmp(:,rad))
         endif
         if(kn_masses(em).gt.0.and.kn_masses(rad).eq.0) then
            getdistance=(2*dotp(cmp(:,em),cmp(:,rad))*
     1        e_rad/e_em
     2        )**par_dijexp
         else
            getdistance=(2*dotp(cmp(:,em),cmp(:,rad))*
     1           e_em*e_rad/(e_em+e_rad)**2
     2           )**par_dijexp
         endif
      endif
      end

      function mergeisr(i,j)
      implicit none
      integer mergeisr,i,j,onem
      parameter (onem = 1000000)
      if(abs(i).gt.5.or.abs(j).gt.5) then
         mergeisr = onem
      elseif(j.eq.0) then
         mergeisr = i
      elseif(i.eq.0) then
         mergeisr = -j
      elseif(i.eq.j) then
         mergeisr = 0
      else
         mergeisr = onem
      endif
      end

      function mergefsr(i,j)
      implicit none
      integer mergefsr,i,j,onem
      parameter (onem = 1000000)
      if(abs(i).gt.5.or.abs(j).gt.5) then
         mergefsr = onem
      elseif(j.eq.0) then
         mergefsr = i
      elseif(i.eq.0) then
         mergefsr = j
      elseif(i.eq.-j) then
         mergefsr = 0
      else
         mergefsr = onem
      endif
      end
