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
      integer j,k,i
      real * 8 dy,dpt,dr
      character * 1 cnum(9)
      data cnum/'1','2','3','4','5','6','7','8','9'/
      integer maxjet
      parameter (maxjet=4)
      integer nptmin
      parameter (nptmin=6)
      character * 4 cptmin(nptmin)
      real * 8 ptminarr(nptmin)
      data cptmin/  '-000','-020','-040','-060','-080','-100'/
      data ptminarr/ 0d0, 20d0,  40d0,  60d0,  80d0, 100d0/
      common/infohist/ptminarr,cnum,cptmin
      save /infohist/
      real * 8 powheginput
      external powheginput
      integer nptbins1
      parameter (nptbins1 = 20)
      integer nptbins2
      parameter (nptbins2 = 20)
      integer netabins
      parameter (netabins = 20)
      real*8 ptbins1(nptbins1 + 1)
      real*8 ptbins2(nptbins2 + 1)
      real*8 etabins(netabins + 1)
      data ptbins1/60d0, 86.8d0, 117.8d0,153.1d0,192.6d0,236.3d0,
     1             284.3d0,336.5d0,392.9d0,453.5d0,518.4d0,587.5d0,
     2             660.9d0,738.5d0,820.3d0,906.3d0,996.6d0,1091.1d0,
     3             1189.8d0,1292.8d0,1400d0/
      data ptbins2/60d0,72.8d0,87.6d0,104.5d0,123.3d0,144.2d0,
     1             167.1d0,192d0,219d0,248d0,278.9d0,312d0,347d0,
     2             384d0,423.1d0,464.2d0,507.3d0,552.5d0,599.6d0,
     3             648.8d0,700d0/
      data etabins/-2.8d0,-2.52d0,-2.24d0,-1.96d0,-1.68d0,-1.4d0,
     1     -1.12d0,-0.84d0,-0.56d0,-0.28d0,0.0d0,0.28d0,0.56d0,
     2     0.84d0,1.12d0,1.4d0,1.68d0,1.96d0,2.24d0,2.52d0,2.8d0/
      call inihists

      dy=0.5d0
      dpt=10d0
      dr=0.2d0
      
      do i=1,nptmin
      call bookupeqbins('sigtot'//cptmin(i),1d0,0.5d0,1.5d0)      
      call bookupeqbins('Njet'//cptmin(i),1d0,-0.5d0,5.5d0)

      do j=1,maxjet
         call bookupeqbins('j'//cnum(j)//'-y'//cptmin(i),dy,-5d0,5d0)
         call bookupeqbins('j'//cnum(j)//'-eta'//cptmin(i),dy,-5d0,5d0)
         call bookupeqbins('j'//cnum(j)//'-pt'//cptmin(i),dpt,0d0,400d0)
         call bookupeqbins('j'//cnum(j)//'-ptzoom'//cptmin(i),
     $        2d0,1d0,151d0)
         call bookupeqbins('j'//cnum(j)//'-m'//cptmin(i),dpt,0d0,400d0) 
         call bookupeqbins('j'//cnum(j)//'-ptzoom2'//cptmin(i),
     $        0.5d0,0d0,20d0)
      enddo

c      goto 300

      do j=1,maxjet-1
         do k=j+1,maxjet
            call bookupeqbins('j'//cnum(j)//'j'//cnum(k)//
     1           '-y'//cptmin(i),dy,-5d0,5d0)  
            call bookupeqbins('j'//cnum(j)//'j'//cnum(k)//
     1           '-eta'//cptmin(i),dy,-5d0,5d0)
            call bookupeqbins('j'//cnum(j)//'j'//cnum(k)//
     1           '-pt'//cptmin(i),dpt,0d0,400d0)
            call bookupeqbins('j'//cnum(j)//'j'//cnum(k)//
     1           '-m'//cptmin(i),dpt,0d0,400d0)  
         enddo
      enddo
   
      do j=1,maxjet-1
         do k=j+1,maxjet
            call bookupeqbins('j'//cnum(j)//'j'//cnum(k)//
     1           '-dy'//cptmin(i),dy,-5d0,5d0)  
            call bookupeqbins('j'//cnum(j)//'j'//cnum(k)//
     1           '-deta'//cptmin(i),dy,-5d0,5d0)
            call bookupeqbins('j'//cnum(j)//'j'//cnum(k)//
     1           '-delphi'//cptmin(i),pi/20,0d0,pi)
            call bookupeqbins('j'//cnum(j)//'j'//cnum(k)//
     1           '-dr'//cptmin(i),dr,0d0,10d0)  
         enddo
      enddo
  
      do j=1,maxjet
         call bookupeqbins('ptrel'//cnum(j)//cptmin(i),0.5d0,0d0,20d0)
      enddo   

 300  continue

      enddo
      
      call bookupeqbins('sigma',1.0d0,0.0d0,3.0d0)
      call bookupeqbins('njet',1d0,-0.5d0,5.5d0)
      call bookup('J pt 1st',nptbins1,ptbins1)
      call bookup('J pt 2nd',nptbins1,ptbins1)
      call bookup('J pt 3rd',nptbins2,ptbins2)
      call bookup('y 1st',netabins,etabins)
      call bookup('y 2nd',netabins,etabins)
      call bookup('y 3rd',netabins,etabins)
      call bookup('eta 1st',netabins,etabins)
      call bookup('eta 2nd',netabins,etabins)
      call bookup('eta 3rd',netabins,etabins)

      end
     

      subroutine analysis(dsig0)
      implicit none
      real * 8 dsig0
      include 'hepevt.h'
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_math.h' 
      include 'pwhg_rad.h' 
      include 'pwhg_weights.h'
c      include 'pwhg_kn.h' 
c      include 'pwhg_flg.h'
c      include 'LesHouches.h'
      integer   maxjet,mjets,njets,numjets,ntracks
      parameter (maxjet=2048)
      real * 8  ktj(maxjet),etaj(maxjet),rapj(maxjet),
     1    phij(maxjet),pj(4,maxjet),rr,ptrel(4),ptot(4)
      integer maxtrack
      parameter (maxtrack=2048)
      real * 8  ptrack(4,maxtrack)
      integer   jetvec(maxtrack),itrackhep(maxtrack)
      character * 1 cnum(9)
      integer nptmin
      parameter (nptmin=6)
      character * 4 cptmin(nptmin)
      real * 8 ptminarr(nptmin)      
      common/infohist/ptminarr,cnum,cptmin
      save /infohist/
      integer j,k,i,jj,ii
c     we need to tell to this analysis file which program is running it
      character * 6 WHCPRG
      common/cWHCPRG/WHCPRG
      data WHCPRG/'NLO   '/
      integer ih,il,inu
      real * 8 psum(4)
      real * 8 httot,y,eta,pt,m
      real * 8 dy,deta,delphi,dr
      integer ihep
      real * 8 powheginput,dotp
      external powheginput,dotp
      real * 8 ptmin,mininvmass
      integer maxnumlep
      parameter (maxnumlep=10)
      real * 8 pvl(4,maxnumlep),plep(4,maxnumlep)
      integer mu,ilep,ivl,nlep,nvl
      logical is_W
      real * 8 mV2,ptvb,mvb,ptlep,ptminfastjet,ptvl,R,ylep,yvb,yvl
      real * 8 Wmass,Wwidth,Wmasslow,Wmasshigh
      integer jpart, jjet
      real * 8 palg
c      real * 8 rescfac1,rescfac2
c      common /crescfac/rescfac1,rescfac2
      real * 8 dsig(7)
      integer nweights
      logical ini
      data ini/.true./
      logical inimulti
      data inimulti/.true./
      integer  minlo
      data minlo/0/
      save inimulti,minlo,ini
      integer minnumjets,maxnumjets
      logical pwhg_isfinite
      external pwhg_isfinite
      character * 4 prefix
      real*8 MinJetPt
      parameter (MinJetPt = 60d0)
      real*8 Min1stJetPt
      parameter (Min1stJetPt = 80d0)
      real*8 MaxEtaJet
      parameter (MaxEtaJet = 2.8d0)
      real*8 MaxRapJet
      parameter (MaxRapJet = 2.8d0)
      real*8 Rpar
      parameter (Rpar = 0.4d0)
      real*8 pj1(4),pj2(4),pj3(4)
      real*8 ptj1,ptj2,ptj3
      real*8 etaj1,etaj2,etaj3
      real*8 yj1,yj2,yj3
      real*8 mj1,mj2,mj3

      if (.not.pwhg_isfinite(dsig0)) then
         write(*,*) '*** PROBLEMS in subroutine analysis ***'
c     in general the following quantities are NOT available. They should be at the NLO stage
         return
      endif

c      write(*,*) 'nhep ',nhep

      maxnumjets=4

c      call reweightifneeded(dsig0,dsig)

      if(inimulti) then
         if(weights_num.eq.0) then
            call setupmulti(1)
         else
            call setupmulti(weights_num)
         endif
         inimulti=.false.
      endif

      dsig=0
      if(weights_num.eq.0) then
         dsig(1)=dsig0
         nweights=1
      else
         dsig(1:weights_num)=weights_val(1:weights_num)
          nweights=weights_num
      endif

      if(sum(abs(dsig)).eq.0) return

      if (ini) then
         minlo=powheginput('#minlo')
         if (minlo.eq.1) then
            minnumjets = 2
         else
            minnumjets = 3
         endif
         ini=.false.
      endif

c     set up arrays for jet finding
c      do jpart=1,maxtrack
c         do mu=1,4
c            ptrack(mu,jpart)=0d0
c         enddo
c         jetvec(jpart)=0
c      enddo      
c      do jjet=1,maxjet
c         do mu=1,4
c            pj(mu,jjet)=0d0
c         enddo
c      enddo

      ntracks=0
      mjets=0
c     Loop over final state particles to find jets 
      do ihep=1,nhep
         if (isthep(ihep).eq.1) then
           if (ntracks.eq.maxtrack) then
              write(*,*) 'Too many particles. Increase maxtrack.'//
     #             ' PROGRAM ABORTS'
              call pwhg_exit(-1)
           endif
c     copy momenta to construct jets 
           ntracks=ntracks+1
           do mu=1,4
              ptrack(mu,ntracks)=phep(mu,ihep)
           enddo
        endif
      enddo



c      goto 111


      if (ntracks.eq.0) then
         numjets=0
      else
c     palg=1 is standard kt, -1 is antikt
         palg = -1d0
         R = 0.5d0              ! radius parameter
c         ptminfastjet = 1d0
         ptminfastjet = ptminarr(1)
         call fastjetppgenkt(ptrack,ntracks,R,palg,ptminfastjet,
     $        pj,numjets,jetvec)
c         call fastjetktwhich(ptrack,ntracks,ptminfastjet,R,
c     $        pj,mjets,jetvec) 
      endif

c      write(*,*) 'numjets ',numjets


      do i=1,nptmin        
         njets=0    
         do j=1,min(maxnumjets,numjets)
            ktj(j) = sqrt(pj(1,j)**2 + pj(2,j)**2 )            
            if (ktj(j).gt.ptminarr(i)) then
               njets=njets+1
            endif
         enddo


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         if (njets.lt.minnumjets) exit
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


         mininvmass=1d30
         do ii=1,njets
            do j=ii+1,njets
               ptot(:)=pj(:,ii)+pj(:,j)
               call pwhg_getinvmass(ptot(:),m)
               mininvmass=min(mininvmass,m)
            enddo
         enddo
         

c********************************************************
c         if (mininvmass.lt.100d0) exit
c********************************************************



c         write(*,*) 'njets ',njets


         call filld('sigtot'//cptmin(i),1d0,dsig)
                  
         if(njets.eq.0) then
            call filld('Njet'//cptmin(i),0d0,dsig)
         elseif(njets.eq.1) then
            call filld('Njet'//cptmin(i),1d0,dsig)
         elseif(njets.eq.2) then
            call filld('Njet'//cptmin(i),2d0,dsig)
         elseif(njets.eq.3) then
            call filld('Njet'//cptmin(i),3d0,dsig)
         elseif(njets.eq.4) then
            call filld('Njet'//cptmin(i),4d0,dsig)
         elseif(njets.eq.5) then
            call filld('Njet'//cptmin(i),5d0,dsig)
         else
c     write(*,*) ' Njet?',mjets
         endif
         
c     jets
         mjets=min(njets,maxnumjets)
         
         do j=1,mjets
            call getyetaptmass(pj(:,j),y,eta,pt,m)
            call filld('j'//cnum(j)//'-y'//cptmin(i),     y, dsig)
            call filld('j'//cnum(j)//'-eta'//cptmin(i), eta, dsig)
            call filld('j'//cnum(j)//'-pt'//cptmin(i),   pt, dsig)
            call filld('j'//cnum(j)//'-ptzoom'//cptmin(i),   pt, dsig)
            call filld('j'//cnum(j)//'-ptzoom2'//cptmin(i),   pt, dsig)
            call filld('j'//cnum(j)//'-m'//cptmin(i),     m, dsig)
c     call filld('ptrel'//cnum(j)//cptmin(i),ptrel(j), dsig)         
         enddo
         


         do j=1,mjets
            do k=j+1,mjets
               call getyetaptmass(pj(:,j)+pj(:,k),y,eta,pt,m)
               call filld('j'//cnum(j)//'j'//cnum(k)//'-y'//cptmin(i),
     $              y, dsig)
               call filld('j'//cnum(j)//'j'//cnum(k)//'-eta'//cptmin(i),
     $              eta, dsig)
               call filld('j'//cnum(j)//'j'//cnum(k)//'-pt'//cptmin(i),
     $              pt, dsig)
               call filld('j'//cnum(j)//'j'//cnum(k)//'-m'//cptmin(i), 
     $              m, dsig)
            enddo
         enddo
         
         do j=1,mjets
            do k=j+1,mjets
               prefix = 'j'//cnum(j)//'j'//cnum(k)               
               call getdydetadphidr(pj(:,j),pj(:,k),dy,deta,delphi,dr)
               call filld(prefix//'-dy'//cptmin(i),dy,dsig)
               call filld(prefix//'-deta'//cptmin(i),deta,dsig)
               call filld(prefix//'-delphi'//cptmin(i),delphi,dsig)
               call filld(prefix//'-dr'//cptmin(i),dr,dsig)
            enddo
         enddo

         
      enddo

c      call fastjetppgenkt_BBUY(ptrack,ntracks,Rpar,-1d0,MinJetPt,
c     $     MaxRapJet,pj,njets,jetvec)




 111  continue

      palg = -1d0
      R = 0.4d0                 ! radius parameter
c     ptminfastjet = 1d0
      ptminfastjet = 60d0
      call fastjetppgenkt(ptrack,ntracks,R,palg,ptminfastjet,
     $     pj,njets,jetvec)
      
c We need at least 3 jets:
      if (njets.lt.3) return
c
c We calculate the pt's, eta's and rapidities for the three leading jets:
      pj1 = pj(:,1)
      pj2 = pj(:,2)
      pj3 = pj(:,3)
c
      call getyetaptmass(pj1,yj1,etaj1,ptj1,mj1)
      call getyetaptmass(pj2,yj2,etaj2,ptj2,mj2)
      call getyetaptmass(pj3,yj3,etaj3,ptj3,mj3)

c     we require the first to be harder
      if (ptj1.lt.80d0) return
      if (abs(yj1).gt.2.8.or.abs(yj2).gt.2.8.or.abs(yj3).gt.2.8) return

      call filld('sigma',1.5d0,dsig)

      if(njets.eq.0) then
         call filld('njet',0d0,dsig)
      elseif(njets.eq.1) then
         call filld('njet',1d0,dsig)
      elseif(njets.eq.2) then
         call filld('njet',2d0,dsig)
      elseif(njets.eq.3) then
         call filld('njet',3d0,dsig)
      elseif(njets.eq.4) then
         call filld('njet',4d0,dsig)
      elseif(njets.eq.5) then
         call filld('njet',5d0,dsig)
      endif
      call filld('J pt 1st',ptj1,dsig)
      call filld('J pt 2nd',ptj2,dsig)
      call filld('J pt 3rd',ptj3,dsig)
      call filld('y 1st',yj1,dsig)
      call filld('y 2nd',yj2,dsig)
      call filld('y 3rd',yj3,dsig)
      call filld('eta 1st',etaj1,dsig)
      call filld('eta 2nd',etaj2,dsig)
      call filld('eta 3rd',etaj3,dsig)




      end




c      subroutine yetaptmassplot(p,dsig,prefix)
c      implicit none
c      real * 8 p(4),dsig
c      character *(*) prefix
c      real * 8 y,eta,pt,m
c      call getyetaptmass(p,y,eta,pt,m)
c      call filld(prefix//'-y',y,dsig)
c      call filld(prefix//'-eta',eta,dsig)
c      call filld(prefix//'-pt',pt,dsig)
c      call filld(prefix//'-m',m,dsig)
c      end

      subroutine deltaplot(p1,p2,dsig,prefix,postfix)
      implicit none
      real * 8 p1(4),p2(4),dsig(7)
      character *(*) prefix,postfix
      real * 8 dy,deta,delphi,dr
      call getdydetadphidr(p1,p2,dy,deta,delphi,dr)
      call filld(prefix//'-dy'//postfix,dy,dsig)
      call filld(prefix//'-deta'//postfix,deta,dsig)
      call filld(prefix//'-delphi'//postfix,delphi,dsig)
      call filld(prefix//'-dr'//postfix,dr,dsig)
      end

      subroutine getyetaptmass(p,y,eta,pt,mass)
      implicit none
      real * 8 p(4),y,eta,pt,mass
      call pwhg_getrapidity(p,y)      
      pt=sqrt(p(1)**2+p(2)**2)
      call pwhg_getpseudorapidity(p,eta)
      call pwhg_getinvmass(p,mass)
      end

      subroutine getdydetadphidr(p1,p2,dy,deta,dphi,dr)
      implicit none
      include 'pwhg_math.h' 
      real * 8 p1(*),p2(*),dy,deta,dphi,dr
      real * 8 y1,eta1,pt1,mass1,phi1
      real * 8 y2,eta2,pt2,mass2,phi2
      call getyetaptmass(p1,y1,eta1,pt1,mass1)
      call getyetaptmass(p2,y2,eta2,pt2,mass2)
      dy=y1-y2
      deta=eta1-eta2
      phi1=atan2(p1(1),p1(2))
      phi2=atan2(p2(1),p2(2))
      dphi=abs(phi1-phi2)
      dphi=min(dphi,2d0*pi-dphi)
      dr=sqrt(deta**2+dphi**2)
      end


      subroutine sortbypt(n,iarr)
      implicit none
      integer n,iarr(n)
      include 'hepevt.h'
      integer j,k
      real * 8 tmp,pt(nmxhep)
      logical touched
      do j=1,n
         pt(j)=sqrt(phep(1,iarr(j))**2+phep(2,iarr(j))**2)
      enddo
c bubble sort
      touched=.true.
      do while(touched)
         touched=.false.
         do j=1,n-1
            if(pt(j).lt.pt(j+1)) then
               k=iarr(j)
               iarr(j)=iarr(j+1)
               iarr(j+1)=k
               tmp=pt(j)
               pt(j)=pt(j+1)
               pt(j+1)=tmp
               touched=.true.
            endif
         enddo
      enddo
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

      function get_ptrel(pin,pjet)
      implicit none
      real * 8 get_ptrel,pin(0:3),pjet(0:3)
      real * 8 pin2,pjet2,cth2,scalprod
      pin2  = pin(1)**2 + pin(2)**2 + pin(3)**2
      pjet2 = pjet(1)**2 + pjet(2)**2 + pjet(3)**2
      scalprod = pin(1)*pjet(1) + pin(2)*pjet(2) + pin(3)*pjet(3)
      cth2 = scalprod**2/pin2/pjet2
      get_ptrel = sqrt(pin2*abs(1d0 - cth2))
      end

      subroutine pwhgfinalopshist
      end
