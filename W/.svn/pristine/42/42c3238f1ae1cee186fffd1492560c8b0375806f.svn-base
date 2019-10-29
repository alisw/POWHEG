c  The next subroutines, open some histograms and prepare them 
c      to receive data 
c  You can substitute these  with your favourite ones
c  init   :  opens the histograms
c  topout :  closes them
c  pwhgfill  :  fills the histograms with data

      subroutine init_hist
      implicit none
      include  'LesHouches.h'
      include 'pwhg_math.h'
      integer j,k
      character * 1 cnum(9)
      data cnum/'1','2','3','4','5','6','7','8','9'/

      call inihists

      call bookupeqbins('W_pt',0.2d0,0d0,100d0)
      call bookupeqbins('W_pt-y-gt-0',0.2d0,0d0,100d0)
      call bookupeqbins('W_pt-y-lt-0',0.2d0,0d0,100d0)
      call bookupeqbins('Y_W',1d0,-5d0,5d0)

      end
     
      subroutine analysis(dsig0)
      implicit none
      real * 8 dsig0
      include 'pwhg_weights.h'
      include 'pwhg_lhrwgt.h'
      real * 8 dsig(10)
      include 'hepevt.h'
      include 'pwhg_math.h' 
      include  'LesHouches.h'
      character * 6 WHCPRG
      common/cWHCPRG/WHCPRG
      data WHCPRG/'NLO   '/
      real * 8 pw(4),pl(4),pnu(4)
      real * 8 y,eta,pt,m
      integer ihep,absid
      real * 8 powheginput,dotp
      external powheginput,dotp
      logical ini
      data ini/.true./
      save ini

      if(ini) then
         if(lhrwgt_nids.gt.0) then
            call setupmulti(lhrwgt_nids)
         elseif(weights_num.gt.0) then
            call setupmulti(weights_num+1)
         else
            call setupmulti(1)
         endif
      endif

      dsig=0
      if(lhrwgt_nids.gt.0) then
         dsig(1:lhrwgt_nids)=lhrwgt_weights(1:lhrwgt_nids)
      elseif(weights_num.gt.0) then
         dsig(1)=dsig0
         dsig(2:weights_num+1)=weights_val(1:weights_num)
      else
         dsig(1)=dsig0
      endif

      if(sum(abs(dsig)).eq.0) return

      do ihep=1,nhep
c p_W = p_l + p_nu
         absid = abs(idhep(ihep))
         if(isthep(ihep) .eq. 1) then
            if( absid .ge. 11 .and. absid .le. 16) then
               if(mod(absid,2) .eq. 0) then
                  pnu  = phep(1:4,ihep)
               else
                  pl = phep(1:4,ihep)
               endif
            endif
         endif
      enddo

      pw = pl + pnu

      call getyetaptmass(pw,y,eta,pt,m)

      call filld('W_pt',pt,dsig)
      if(y .gt. 0) then
         call filld('W_pt-y-gt-0',pt,dsig)
      else
         call filld('W_pt-y-lt-0',pt,dsig)
      endif
      call filld('Y_W',y,dsig)
      end

      subroutine yetaptmassplot(p,dsig,prefix)
      implicit none
      real * 8 p(4),dsig
      character *(*) prefix
      real * 8 y,eta,pt,m
      call getyetaptmass(p,y,eta,pt,m)
      call filld(prefix//'_y',y,dsig)
      call filld(prefix//'_eta',eta,dsig)
      call filld(prefix//'_pt',pt,dsig)
      call filld(prefix//'_m',m,dsig)
      end

      subroutine deltaplot(p1,p2,dsig,prefix)
      implicit none
      real * 8 p1(4),p2(4),dsig
      character *(*) prefix
      real * 8 dy,deta,delphi,dr
      call getdydetadphidr(p1,p2,dy,deta,delphi,dr)
      call filld(prefix//'_dy',dy,dsig)
      call filld(prefix//'_deta',deta,dsig)
      call filld(prefix//'_delphi',delphi,dsig)
      call filld(prefix//'_dr',dr,dsig)
      end


      subroutine getyetaptmass(p,y,eta,pt,mass)
      implicit none
      real * 8 p(4),y,eta,pt,mass,pv
      real * 8 gety,getpt,geteta,getmass
      external gety,getpt,geteta,getmass
      y  = gety(p)
      pt = getpt(p)
      eta = geteta(p)
      mass = getmass(p)
      end


      function gety(p)
      implicit none
      real * 8 gety,p(4)
      gety=0.5d0*log((p(4)+p(3))/(p(4)-p(3)))
      end

      function getpt(p)
      implicit none
      real * 8 getpt,p(4)
      getpt = sqrt(p(1)**2+p(2)**2)
      end

      function getmass(p)
      implicit none
      real * 8 getmass,p(4)
      getmass=sqrt(abs(p(4)**2-p(3)**2-p(2)**2-p(1)**2))
      end

      function geteta(p)
      implicit none
      real * 8 geteta,p(4),pv
      real * 8 tiny
      parameter (tiny=1.d-5)
      pv = sqrt(p(1)**2+p(2)**2+p(3)**2)
      if(pv.lt.tiny)then
         geteta=sign(1.d0,p(3))*1.d8
      else
         geteta=0.5d0*log((pv+p(3))/(pv-p(3)))
      endif
      end



      subroutine getdydetadphidr(p1,p2,dy,deta,dphi,dr)
      implicit none
      real * 8 p1(*),p2(*),dy,deta,dphi,dr
      real * 8 getdy,getdeta,getdphi,getdr
      external getdy,getdeta,getdphi,getdr
      dy=getdy(p1,p2)
      deta=getdeta(p1,p2)
      dphi=getdphi(p1,p2)
      dr=getdr(deta,dphi)
      end

      function getdy(p1,p2)
      implicit none
      real*8 p1(*),p2(*),getdy
      real*8 y1,y2
      real*8 gety
      external gety
      y1 = gety(p1)
      y2 = gety(p2)
      getdy = y1-y2
      end

      function getdeta(p1,p2)
      implicit none
      real*8 p1(*),p2(*),getdeta
      real*8 eta1,eta2
      real*8 geteta
      external geteta
      eta1 = geteta(p1)
      eta2 = geteta(p2)
      getdeta = eta1-eta2
      end

      function getdphi(p1,p2)
      implicit none
      include 'pwhg_math.h' 
      real*8 p1(*),p2(*),getdphi
      real*8 phi1,phi2
      real*8 geteta
      external geteta
      phi1=atan2(p1(2),p1(1))
      phi2=atan2(p2(2),p2(1))
      getdphi=abs(phi1-phi2)
      getdphi=min(getdphi,2d0*pi-getdphi)
      end

      function getdr(deta,dphi)
      implicit none
      real*8 getdr,deta,dphi 
      getdr=sqrt(deta**2+dphi**2)
      end

      function islept(j)
      implicit none
      logical islept
      integer j
      if(abs(j).ge.11.and.abs(j).le.15) then
         islept = .true.
      else
         islept = .false.
      endif
      end


      function dotp3(p1,p2)
      implicit none
      real * 8 dotp3,p1(3),p2(3)
      dotp3 = p1(1)*p2(1) + p1(2)*p2(2) + p1(3)*p2(3)
      end

      subroutine cross3(p1,p2,p3)
      implicit none
      real * 8 p3(3),p1(3),p2(3)
      p3(1) = p1(2)*p2(3)-p1(3)*p2(2)
      p3(2) = p1(3)*p2(1)-p1(1)*p2(3)
      p3(3) = p1(1)*p2(2)-p1(2)*p2(1)
      end
