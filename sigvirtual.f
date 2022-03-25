      subroutine sigvirtual(virt_arr)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_br.h'
      include 'pwhg_flg.h'
      real * 8 virt_arr(maxprocborn)
      integer equivto(maxprocborn)
      common/cequivtovirt/equivto
      real * 8 equivcoef(maxprocborn)
      common/cequivcoefvirt/equivcoef
      integer nmomset
      parameter (nmomset=10)
      real * 8 pborn(0:3,nlegborn,nmomset),cprop
      real * 8 virtual(nmomset,maxprocborn)
      integer iborn,ibornpr,mu,nu,k,j,iret
      logical ini
      data ini/.true./
      save ini,/cequivtovirt/,/cequivcoefvirt/
      logical pwhg_isfinite
      real * 8 save_kn_jacborn
      real * 8 powheginput
      external pwhg_isfinite,powheginput
      if(ini) then
         do iborn=1,flst_nborn
            equivto(iborn)=-1
         enddo
         if(flg_smartsig.and..not.flg_novirtual) then
            flg_in_smartsig = .true.
            call fillequivarrayvirt(flst_nborn,equivto,equivcoef,iret)
            if(iret<0) then
               call randomsave
               save_kn_jacborn = kn_jacborn
c     now set kn_jacborn to a value different from zero, in order for the virtual
c     to be evaluated, when called wth the momenta returned by fillmomenta
               kn_jacborn = 1d0               
               call fillmomenta(nlegborn,nmomset,kn_masses,pborn)
               do iborn=1,flst_nborn
                  do j=1,nmomset
                     flst_cur_iborn = iborn 
                     call setvirtual(pborn(0,1,j),flst_born(1,iborn),
     1                    virtual(j,iborn))
c     check if virtual(j,iborn) is finite
                     if (.not.pwhg_isfinite(virtual(j,iborn))) 
     1                    virtual(j,iborn)=0d0
                  enddo
                  call compare_vecsv(nmomset,iborn,virtual,ibornpr,
     1                 cprop,iret)
                  if(iret.eq.0) then
                     equivto(iborn)=ibornpr
                     equivcoef(iborn)=1
                  elseif(iret.eq.1) then
                     equivto(iborn)=ibornpr
                     equivcoef(iborn)=cprop
                  endif
               enddo
               call randomrestore
               kn_jacborn = save_kn_jacborn 
               call printvirtequiv
c     Write equiv file, if required
               if(powheginput('#writeequivfile') == 1) then
                  call writeequivfile('virt',
     1                 flst_nborn,equivto,equivcoef)
               endif
            endif
         endif
         flg_in_smartsig = .false.
         ini=.false.
      endif
      do iborn=1,flst_nborn
         if(equivto(iborn).lt.0) then
            flst_cur_iborn = iborn
            call setvirtual(kn_cmpborn,flst_born(1,iborn),
     #           virt_arr(iborn))
c     check if virt_arr(iborn) is finite
                  if (.not.pwhg_isfinite(virt_arr(iborn))) 
     #                 virt_arr(iborn)=0d0
            virt_arr(iborn)=virt_arr(iborn)/(2*kn_sborn)
         else
            virt_arr(iborn)=virt_arr(equivto(iborn))*equivcoef(iborn)
         endif
      enddo
      end

      subroutine compare_vecsv(nmomset,iborn,res,ibornpr,cprop,iret)
      implicit none
      real * 8 ep
      integer nmomset,iborn,ibornpr,iret,j,k
      real * 8 res(nmomset,iborn),cprop,rat
      real * 8 powheginput,tmp
      tmp = powheginput("#compare_vecsv_ep")
      if(tmp>0) then
         ep = tmp
      else
         ep = 1d-8
      endif
      do j=1,iborn-1
         rat=res(1,iborn)/res(1,j)
         do k=1,nmomset
            if(abs(1-res(k,iborn)/res(k,j)/rat).gt.ep) goto 10
         enddo
         if(abs(1-rat).lt.ep) then
            iret=0
            cprop=1
         else
            iret=1
            cprop=rat
         endif
         ibornpr=j
         return
 10      continue
      enddo
      iret=-1
      end


      subroutine printvirtequiv
c When invoked after the first call to virtuals,
c it prints the set of equivalent virtual configurations
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_rnd.h'
      integer equivto(maxprocborn)
      common/cequivtovirt/equivto
      real * 8 equivcoef(maxprocborn)
      common/cequivcoefvirt/equivcoef
      integer j,k,iun,count
      save count
      data count/0/
      call newunit(iun)
      if(rnd_cwhichseed == 'none') then
         open(unit=iun,file='virtequiv',status='unknown')
      else
         open(unit=iun,file='virtequiv-'
     1        //trim(rnd_cwhichseed),status='unknown')
      endif
      write(*,*) 'Writing virtequiv file...'
      do j=1,flst_nborn
         if(equivto(j).eq.-1) then
            write(iun,'(a)')
     1           'Beginning sequence of equivalent amplitudes'
c            write(iun,100) 1d0,j, flst_born(:,j)
            write(iun,101) j, 1d0, flst_born(:,j)
            do k=1,flst_nborn
               if(equivto(k).eq.j) then
c                  write(iun,100) equivcoef(k),k,flst_born(:,k)
                  write(iun,101) k,equivcoef(k),flst_born(:,k)
               endif
            enddo
            count=count+1
         endif
      enddo
      write(iun,*) ''
      write(iun,'(a,i4,a)') 'Found ',count, ' equivalent groups'
      close(iun)
      write(*,*) 'Done'
c 100  format(d11.4,5x,i4,5x,100(i4,1x))
 101  format(i4,5x,d11.4,5x,100(i4,1x))
      end
