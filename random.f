
      function random()
      real * 8 random
      real * 8 saverandom
      logical fixed
      COMMON/tmpfixed/fixed
      data fixed/.false./
      save saverandom
      if(fixed) then
         random=saverandom
         return
      endif
      call rm48(random,1)
      saverandom=random
      end


      subroutine resetrandom
      call RM48IN(54217137,0,0)
      end

      subroutine randomsave
      implicit none
      integer j,ipar(3,10)
      data j/0/
      save j,ipar
      j=j+1
      if(j.gt.10) then
         write(*,*) ' Too many recursive calls to randomsave'
         stop
      endif
      call rm48ut(ipar(1,j),ipar(2,j),ipar(3,j))
      return
      entry randomrestore
      if(j.le.0) then
         write(*,*) ' Too many calls to randomrestore'
         stop
      endif
      call setrandom(ipar(1,j),ipar(2,j),ipar(3,j))
      j=j-1
      return
      end

      subroutine readcurrentrandom(i1,n1,n2)
      implicit none
      integer i1,n1,n2
      call rm48ut(i1,n1,n2)
      end

      subroutine setrandom(i1,n1,n2)
      implicit none
      integer, parameter :: modcns=1000000000
      integer i1,n1,n2
      integer i1cur,n1cur,n2cur,j,k,seconds,iret
      real * 8 tmp
      if(i1.eq.0) then
c This is used for complete initialization
         call resetrandom
         return
      endif
      call reset_timer
c Reinitializing the random number may be expensive;
c If we need a sequence number greater than the current status
c just call the generator enough times to get there. 
      call rm48ut(i1cur,n1cur,n2cur)
      if( .not. (i1.eq.i1cur.and.(n2.gt.n2cur .or.
     1     (n2.eq.n2cur .and.n1.ge.n1cur)))) then
         call cachesetfullrandomstatus(i1cur,n1cur,n2cur,iret)
         if(iret == -1) goto 1
         call rm48ut(i1cur,n1cur,n2cur)
      endif
      do k=n2cur,n2-1
         write(*,*) 'setrandom: skipping ', modcns-n1cur, ' ...'
         do j=n1cur,modcns-1
            call rm48(tmp,1)
         enddo
         n1cur = 0
      enddo
      do j=n1cur,n1-1
         call rm48(tmp,1)
      enddo
      call rm48ut(i1cur,n1cur,n2cur)
      if(i1.eq.i1cur.and.n2.eq.n2cur.and.n1.eq.n1cur) then
c     Succeded
         goto 999
      else
c Failed!
         write(*,*) ' setrandom: failed ...'
         call exit(-1)
      endif
c reinitialize from scratch
 1    if (I1.gt.0.and.n1.ge.0.and.n2.ge.0) then
c     restart a previous run
            call rm48in(I1,N1,N2)
      else
         write(*,*) 'ERROR: setrandom called with',i1,n1,n2
         call exit(-1)
      endif
 999  call get_timer(seconds)
      call addtocnt('setrandom time (sec)',seconds)
      end


      subroutine savecurrentrandom
      implicit none
      integer ipar(3)
      common/crandom/ipar
      call rm48ut(ipar(1),ipar(2),ipar(3))
      end


      subroutine getcurrentrandom(i1,n1,n2)
      implicit none
      integer i1,n1,n2
      integer ipar(3)
      common/crandom/ipar
      i1 = ipar(1)
      n1 = ipar(2)
      n2 = ipar(3)
      end

      subroutine printcurrentrandom
      implicit none
      integer ipar(3)
      common/crandom/ipar
      write(*,*) 'Random number seeds: ',ipar(1),ipar(2), ipar(3)
      end

      subroutine cachefullrandomstatus
      implicit none
      integer,parameter :: ncache=10
      real * 8 U(104)
      integer I97,J97,NTOT, NTOT2, IJKL
      COMMON/R48ST1/U,I97,J97,NTOT,NTOT2,IJKL
      real * 8 cache(104,ncache)
      integer cache_I97(ncache),cache_J97(ncache),
     1     cache_NTOT(ncache), cache_NTOT2(ncache),
     1     cache_IJKL(ncache),ptrcache,filled
      common/randomcache/cache,cache_NTOT,cache_NTOT2,
     1     cache_I97,cache_J97,cache_IJKL,ptrcache,filled
      integer j
      logical, save :: ini=.true.
      if(ini) then
         ptrcache = 0
         filled = 0
         ini = .false.
      endif
      do j=1,filled
         if(ntot == cache_NTOT(j) .and.
     1        ntot2 == cache_NTOT2(j) .and. cache_IJKL(j) == IJKL) then
            return
         endif
      enddo
      if(ptrcache == ncache) then
         ptrcache = 1
      else
         ptrcache = ptrcache+1
         filled = max(filled,ptrcache)
      endif
      cache(:,ptrcache)=U
      cache_J97(ptrcache)=J97
      cache_I97(ptrcache)=I97
      cache_NTOT(ptrcache)=NTOT
      cache_NTOT2(ptrcache)=NTOT2
      cache_IJKL(ptrcache)=IJKL
      end

      subroutine cachesetfullrandomstatus(iseed,n1,n2,iret)
      implicit none
      integer iseed,n1,n2,iret
      integer,parameter :: ncache=10
      real * 8 U(104)
      integer I97,J97,NTOT, NTOT2, IJKL
      COMMON/R48ST1/U,I97,J97,NTOT,NTOT2,IJKL
      real * 8 cache(104,ncache)
      integer cache_I97(ncache),cache_J97(ncache),
     1     cache_NTOT(ncache), cache_NTOT2(ncache),
     1     cache_IJKL(ncache),ptrcache,filled
      common/randomcache/cache,cache_NTOT,cache_NTOT2,
     1     cache_I97,cache_J97,cache_IJKL,ptrcache,filled
      integer j,jfound,dist1,dist2
      jfound = -1
      dist1 = -1
      dist2 = -1
      do j=1,filled
         if(cache_IJKL(j) == iseed .and. cache_NTOT2(j) <= n2
     1    .and. cache_NTOT(j) <= n1) then
            if(dist1 < 0) then
               jfound = j
               dist1 = n1 - cache_NTOT(j)
               dist2 = n2 - cache_NTOT2(j)
            else
               if(n2 - cache_NTOT2(j) > dist2) then
                  jfound = j
                  dist1 = n1 - cache_NTOT(j)
                  dist2 = n2 - cache_NTOT2(j)
               elseif(n1 - cache_NTOT(j) > dist1) then
                  jfound = j
                  dist1 = n1 - cache_NTOT(j)
                  dist2 = n2 - cache_NTOT2(j)
               endif
            endif
         endif
      enddo
      if(jfound /= -1) then
         U      =  cache(:,jfound)
         I97    =  cache_I97(jfound)
         J97    =  cache_J97(jfound)
         NTOT   =  cache_NTOT(jfound)
         NTOT2  =  cache_NTOT2(jfound)
         IJKL   =  cache_IJKL(jfound)
         iret   = 0
      else
         iret = -1
      endif
      end

