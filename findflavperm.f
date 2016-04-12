      program test
      implicit none
      integer nlegs
      parameter (nlegs = 5)
      integer flav(nlegs),arrflav(nlegs,2),perm(nlegs)
     1     ,j
      arrflav(:,1) = (/ -2,3,4,7,8 /)
      arrflav(:,2) = (/ -2,7,5,6,8 /)
 1    write(*,*) arrflav(:,1)
      write(*,*) arrflav(:,2)
      read(*,*) flav
      call permflav(nlegs,flav,2,arrflav,perm)
      if(perm(1).lt.0) then
         write(*,*) ' not found'
      else
         write(*,*) (flav(perm(j)),j=1,nlegs)
         write(*,*)
      endif
      goto 1
      end



      subroutine permflav(nlegs,flav,nlist,arrflav,perm)
c Given a flavour array flav(nlegs), finds the permutatino perm such that
c flav(perm(j)) is present in arrflav.
      implicit none
      integer nlegs,flav(nlegs),nlist,arrflav(nlegs,nlist),perm(nlegs)
      integer jlist,flavtmp(nlegs),j,k,l,itmp
      do j=1,nlegs
         perm(j)=j
      enddo
      do jlist=1,nlist
         if(all(flav(1:2) == arrflav(1:2,jlist))) then 
            flavtmp = flav
            do k=3,nlegs
               do l=k,nlegs
                  if (flavtmp(l).eq.arrflav(k,jlist)) then
                     if(l.ne.k) then
                        itmp = flavtmp(k)
                        flavtmp(k) = flavtmp(l)
                        flavtmp(l) = itmp
                        itmp=perm(k)
                        perm(k)=perm(l)
                        perm(l)=itmp
                     endif
                     exit
                  endif
               enddo
               if(l.eq.nlegs+1) goto 10
            enddo
         else
            goto 10
         endif
c check:
         do j=1,nlegs
            if(flav(perm(j)).ne.arrflav(j,jlist)) then
               write(*,*) ' not working ...'
               call exit(-1)
            endif
         enddo
         return
 10      continue
      enddo
c fails:
      perm(1)=-1
      end
