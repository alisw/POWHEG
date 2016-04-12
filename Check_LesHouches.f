
      subroutine check_leshouches
      implicit none
      include 'LesHouches.h'
      integer chain(100)
      integer j,m

c First check that colour and flavours are consistent

      do j=1,nup
         if( (idup(j).eq.21 .and. icolup(1,j)*icolup(2,j).eq.0) .or.
     1       (idup(j).le.6.and.idup(j).ge.1 .and.
     2                 (icolup(1,j).eq.0.or.icolup(2,j).ne.0) ) .or.
     3       (idup(j).le.-1.and.idup(j).ge.-6 .and.
     4                 (icolup(2,j).eq.0.or.icolup(1,j).ne.0) ) .or.
     5       ((abs(idup(j)).gt.6.and.idup(j).ne.21) .and.
     6          (icolup(2,j).ne.0.or.icolup(1,j).ne.0)) ) then
            write(*,*) ' check_leshouches:'
            write(*,*) ' parton ',j,'had idup=',idup(j),
     1           'and icolup=',icolup(:,j)
         endif
      enddo 

c chain is an integer array for marking a group of particles involved in a
c certain stage of the process. For example, for the production stage
c chain(1) and chain(2) will be set equal to -1 (incoming), for all particles j
c with mothup(1,j)=1 we will have chain(j)=1. and chain(k)=0 for all the others.
c For the decay of the resonance m, chain(m)=-1, chain(j)=1 for all j with mothup(j)=m,
c and zero for all the others.

      chain = 0
      chain(1) = -1
      chain(2) = -1
      do j=3,nup
         if(mothup(1,j).eq.1) then
            chain(j) = 1
         endif
      enddo

      call check_leshouches_chain(chain)

c Now do it for all resonances
      do m=1,nup
         if(istup(m).eq.2) then
            chain = 0
            chain(m)=-1
            do j=3,nup
               if(mothup(1,j).eq.m) then
                  chain(j) = 1
               endif
            enddo
         endif

         call check_leshouches_chain(chain)

      enddo

      end


      subroutine check_leshouches_chain(chain)
      implicit none
      include 'LesHouches.h'
      real * 8 ptot(1:4),charge
      integer chain(100)
      integer j,k,itmp,col,acol,ncol,nacol
      real * 8 chargeofid
      external chargeofid
c invert colour, flavour and momentum for incoming particles
      do j=1,nup
         if(chain(j).eq.-1) then
c invert colour
            itmp = icolup(1,j)
            icolup(1,j) = icolup(2,j)
            icolup(2,j) = itmp
c invert flavour
            call flavourconj(idup(j))
c invert momenta
            pup(:,j)=-pup(:,j)
         endif
      enddo
c check momentum conservation
      ptot = 0
      do j = 1,nup
         if(chain(j).ne.0) then
            ptot = ptot + pup(1:4,j)
         endif
      enddo
      if(sum(abs(ptot)).gt.1d-6) then
         write(*,*) ' check_leshouches_chain:'
         write(*,*) ' momentum violation for system '
         write(*,'(500(1x,i3))') chain(1:nup)
      endif
c colour conservation: a given colour
c can match only a single anticolour
      do j = 1,nup
         if(chain(j).ne.0) then
            col = icolup(1,j)
            acol = icolup(2,j)
            nacol = 0
            ncol = 0
            do k=1,nup
               if(chain(k).ne.0) then
                  if(acol.ne.0.and.icolup(1,k).eq.acol) nacol = nacol+1
                  if(col.ne.0.and.icolup(2,k).eq.col) ncol = ncol+1
               endif
            enddo
            if(       (col .ne.0 .and. ncol.ne.1)
     1           .or. (acol.ne.0 .and. nacol.ne.1) ) then
               write(*,*) ' check_leshouches_chain:'
               write(*,*) ' inconsistent colour assignment for system '
               write(*,'(500(1x,i3))') chain(1:nup)
            endif
         endif
      enddo
c Check electric charge conservation
      charge = 0
      do j = 1,nup
         if(chain(j).ne.0) then
            charge = charge + chargeofid(idup(j))
         endif
      enddo
      if(abs(charge).gt.1d-6) then
         write(*,*) ' check_leshouches_chain:'
         write(*,*) ' charge not conserved for system '
         write(*,'(500(1x,i3))') chain(1:nup)
      endif

c invert back colour, flavour and momentum for incoming particles
      do j=1,nup
         if(chain(j).eq.-1) then
c invert colour
            itmp = icolup(1,j)
            icolup(1,j) = icolup(2,j)
            icolup(2,j) = itmp
c invert flavour
            call flavourconj(idup(j))
c invert momenta
            pup(:,j)=-pup(:,j)
         endif
      enddo

      end




      function chargeofid(id)
      implicit none
      real * 8 chargeofid
      integer id,aid
      logical isutype,isdtype,islepton,isnu
      aid = abs(id)
      if(aid.eq.24) then
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
         write(*,*)
     1        ' this only works for fermions, W, Z, photons and guons'
         write(*,*) ' if you need to extend it, go ahead'
         call pwhg_exit(-1)
      endif
      end

      subroutine flavourconj(id)
      implicit none
      integer id
      if(abs(id).le.18.or.abs(id).eq.24) then
         id=-id
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

