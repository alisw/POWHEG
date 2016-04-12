c determine if entry k is a descendant of entry
c m in the hep event structure. It looks back
c up to ngenerations, default 4
      function comesfrom(m,k)
      implicit none
      logical comesfrom
      integer m,k
      include  'hepevt.h'
      integer j,kcurr
      integer ngenerations
      parameter (ngenerations=4)
      kcurr = k
      do j=1,ngenerations
         if(kcurr.eq.0) then
            comesfrom = .false.
            return
         endif
         if(kcurr.eq.m) then
            comesfrom = .true.
            return
         endif
         kcurr=jmohep(1,kcurr)
      enddo
      comesfrom=.false.
      end

c determine if entry k is a descendant of a
c particle with PDG idntifier=id in the hep event structure.
c  It looks backup to ngenerations, default 4
      function comesfromid(id,k)
      implicit none
      logical comesfromid
      integer id,k
      include  'hepevt.h'
      integer j,kcurr
      integer ngenerations
      parameter (ngenerations=4)
      kcurr = k
      do j=1,ngenerations
         if(kcurr.eq.0) then
            comesfromid = .false.
            return
         endif
         if(idhep(kcurr).eq.id) then
            comesfromid = .true.
            return
         endif
         kcurr = jmohep(1,kcurr)
      enddo
      comesfromid=.false.
      end


c logical functions to identify quarks, charged leptons and neutrinos

      function isquark(id)
      implicit none
      logical isquark
      integer id,aid
      aid=abs(id)
      if(aid.ge.1.and.aid.le.6) then
         isquark=.true.
      else
         isquark=.false.
      endif
      end

      function isutype(id)
      implicit none
      logical isutype
      integer id,aid
      aid=abs(id)
      if(aid.eq.2.or.aid.eq.4.or.aid.eq.6) then
         isutype=.true.
      else
         isutype=.false.
      endif
      end

      function isdtype(id)
      implicit none
      logical isdtype
      integer id,aid
      aid=abs(id)
      if(aid.eq.1.or.aid.eq.3.or.aid.eq.5) then
         isdtype=.true.
      else
         isdtype=.false.
      endif
      end

      function islepton(id)
      implicit none
      logical islepton
      integer id,aid
      aid=abs(id)
      if(aid.eq.11.or.aid.eq.13.or.aid.eq.15) then
         islepton=.true.
      else
         islepton=.false.
      endif
      end

      function isnu(id)
      implicit none
      logical isnu
      integer id,aid
      aid=abs(id)
      if(aid.eq.12.or.aid.eq.14.or.aid.eq.16) then
         isnu=.true.
      else
         isnu=.false.
      endif
      end

      function isewup(id)
c id is in up position in ew doublet
      implicit none
      logical isewup
      integer id
      isewup=2*(id/2).eq.id
      end


      function isewdo(id)
c id is in up position in ew doublet
      implicit none
      logical isewdo
      integer id
      isewdo=.not.2*(id/2).eq.id
      end

      integer function ewgeneration(id)
      integer id,aid
      aid=abs(id)
      if(aid.eq.1.or.aid.eq.2.or.aid.eq.11.or.aid.eq.12) then
         ewgeneration = 1
      elseif(aid.eq.3.or.aid.eq.4.or.aid.eq.13.or.aid.eq.14) then
         ewgeneration = 2
      elseif(aid.eq.5.or.aid.eq.6.or.aid.eq.15.or.aid.eq.16) then
         ewgeneration = 3
      else
         write(*,*)"ewgeneration requested with unknown flavour",id
         call pwhg_exit(-1)
      endif
      end


      function chargeofid(id)
      implicit none
      real * 8 chargeofid
      integer id
      logical isutype,isdtype,islepton,isnu
      if(abs(id).eq.24) then
         chargeofid = sign(1,id)
      elseif(id.eq.21.or.id.eq.22.or.id.eq.23.or.id.eq.25) then
         chargeofid = 0
      elseif(isutype(id)) then
         chargeofid = 2d0/3*sign(1,id)
      elseif(isdtype(id)) then
         chargeofid = -1d0/3*sign(1,id)
      elseif(islepton(id)) then
         chargeofid = -sign(1,id)
      elseif(isnu(id)) then
         chargeofid = 0
      else
         write(*,*) ' chargeofid: warning, charge not known for id=',
     1        id
         write(*,*) ' returning absurd value 1000000'
         chargeofid = 1000000d0
      endif
      end

      subroutine flavourconj(id)
      implicit none
      integer id
      if(abs(id).le.18.or.abs(id).eq.24) then
         id=-id
      endif
      end

      function i3chargeofid(id)
      implicit none
      integer i3chargeofid
      integer id
      logical isutype,isdtype,islepton
      if(abs(id).gt.16) then
         write(*,*) ' this only works for light fermions'
         write(*,*) ' if you need to extend it, go ahead'
         call pwhg_exit(-1)
      endif
      if(isutype(id)) then
         i3chargeofid = 2*sign(1,id)
      elseif(isdtype(id)) then
         i3chargeofid = -1*sign(1,id)
      elseif(islepton(id)) then
         i3chargeofid = -3*sign(1,id)
      else
         i3chargeofid = 0
      endif
      end


      subroutine swaplh(i,j)
c swap elements i and j in the LH record
      implicit none
      include 'LesHouches.h'
      integer i,j,k,l
      real * 8 v(5),v1
      integer itmp,itmp2(2)      
      if(.not. (i.ge.1.and.i.le.nup.and.
     1     j.ge.1.and.j.le.nup)) then
         write(*,*) ' swaplh: invalid entries ii,j=',i,j
         write(*,*) '       : nup=',nup
         call pwhg_exit(-1)
      endif
      if(i.eq.j) return

c Swap:

      itmp=idup(i)
      idup(i)=idup(j)
      idup(j)=itmp

      v=pup(:,i)
      pup(:,i)=pup(:,j)
      pup(:,j)=v
      
      itmp2=mothup(:,i)
      mothup(:,i)=mothup(:,j)
      mothup(:,j)=itmp2
      
      itmp2=icolup(:,i)
      icolup(:,i)=icolup(:,j)
      icolup(:,j)=itmp2
      
      v1=vtimup(i)
      vtimup(i)=vtimup(j)
      vtimup(j)=v1
      
      v1=spinup(i)
      spinup(i)=spinup(j)
      spinup(j)=v1

c update mothup links:
      do k=1,nup
         do l=1,2
            if(k.ne.i.and.k.ne.j) then
               if(mothup(l,k).eq.i) then
                  mothup(l,k)=j
               elseif(mothup(l,k).eq.j) then
                  mothup(l,k)=i
               endif
            endif
         enddo
      enddo

      end

