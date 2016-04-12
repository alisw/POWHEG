      subroutine analysis_driver(dsig0,ikin)
      implicit none
      real * 8 dsig0
      integer ikin
      integer jpart,mu
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'hepevt.h'
      real * 8 powheginput
      logical testplots,ini
      data testplots/.true./
      data ini/.true./
      save testplots,ini
      if(ini) then
         if (powheginput('#testplots').eq.1d0) then
            testplots=.true.
         else
            testplots=.false.
         endif
         ini=.false.
      endif
      if(.not.testplots) return
      if(ikin.eq.0) then
         nhep=nlegborn
         do jpart=1,nhep
            do mu=1,3
               phep(mu,jpart)=kn_pborn(mu,jpart)
            enddo
            phep(4,jpart)=kn_pborn(0,jpart)
            phep(5,jpart)=sqrt(abs(phep(4,jpart)**2-
     #           phep(1,jpart)**2-phep(2,jpart)**2-phep(3,jpart)**2))
            if(jpart.le.2) then
               isthep(jpart)=-1
            else
               isthep(jpart)=1
            endif
c we may not know the flavour af all partons
            idhep(jpart)=0
         enddo
      else
         nhep=nlegreal
         do jpart=1,nhep
            do mu=1,3
               phep(mu,jpart)=kn_preal(mu,jpart)
            enddo
            phep(4,jpart)=kn_preal(0,jpart)
            phep(5,jpart)=sqrt(abs(phep(4,jpart)**2-
     #           phep(1,jpart)**2-phep(2,jpart)**2-phep(3,jpart)**2))
            if(jpart.le.2) then
               isthep(jpart)=-1
            else
               isthep(jpart)=1
            endif
c we may not know the flavour af all partons
            idhep(jpart)=0
         enddo
      endif
c     assign colorless particles'ID and massive partons ID
c     NB: all regions/flavstruct have the same white and colored particles
c     We can then use the flst_born(jpart,1) of the FIRST Born
      do jpart=3,flst_lightpart-1
         idhep(jpart)=flst_born(jpart,1)
      enddo
c     If a massive parton is treated as massless, in find_regions, we add it to idhep
       do jpart=flst_lightpart,nlegborn
          if (kn_masses(jpart).ne.0d0) then
             idhep(jpart)=flst_born(jpart,1)
          endif
       enddo


c     call analysis routine
      call analysis(dsig0)
      end


     
      subroutine lhtohep
      implicit none
      include 'hepevt.h'
      include 'LesHouches.h'
      integer j
      nhep=nup
      do j=1,nup
         idhep(j)=idup(j)
         isthep(j)=istup(j)
         phep(:,j)=pup(:,j)
         jmohep(:,j)=mothup(:,j)
      enddo
      end

c     This is to store in common block extra information on the
c     partonic NLO event
C     type is either: 'born', 'virt', 'real', 'realct'
C     narr number of entries  in the array arr 
C     weight is the integration weight 
      subroutine analysis_extrainfo(type,narr,arr,weight)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flg.h'
      include 'pwhg_anexinf.h'
      character *(*) type
      integer narr
      real * 8 arr(narr),weight
      integer marr,dummy
      
      if(.not.flg_analysisextrainfo) return
      if(narr.gt.maxalr) then
         write(*,*) ' analysis_extrainfo: arr dimension too big '
         write(*,*) ' exiting ...'
         call exit(-1) 
      endif
      anexinf_sigarr(1:narr) = arr * weight
      anexinf_narr = narr
      anexinf_type = type
      if(anexinf_type.ne.type) then
         write(*,*) 'anexinf_type',anexinf_type
         write(*,*) 'type', type
         write(*,*) ' analysis_extrainfo: type string too long '
         write(*,*) ' exiting ...'
         call exit(-1) 
      endif
      end
