c These subroutines are to be used in conjunction with the
c pwhg_bookhist-multi.f package.
c They allow to read a histogram file, reloading the histograms in the
c common blocks, and to manipulate them.
c
c pwhgloadhistos(fname): loads the histograms in previously saved file.
c                        The histograms are added to those already
c                        present in the common block.
c pwhgemptyduphisto(str1,str2): makes an empty duplicate of the histogram
c                               tagged as str1 into str2, creating the latter.
c pwhgoperatehisto(op,str1,str2,str3): str3 = str1 op str2, where op
c                                      is an arithmetic operation (at the
c                                      moment op = '/' only; feel free to add
c                                      more
c pwhgrebin(str1,str2): copy histogram str1 into histogram str2, where str2
c                       can have different binning. Bin edges must match
c                       exactly (i.e. bin edges in str2 must be also edges
c                       of str1). The range of str2 can be smaller than the
c                       range in str1, in which case str1 is clipped.
c
c One example, to be compiled together with pwhg_bookhist-multi.f and
c newunit.f:
c      implicit none
c      integer j
c      call inihists
c      call pwhgloadhisos('pwg-NLO.top')
c      call bookupeqbins('a0-rebin',4d0,0d0,100d0)
c      call bookupeqbins('pt2_rebin',4d0,0d0,100d0)
c      call pwhgrebin('a0','a0-rebin')
c      call pwhgrebin('pt2','pt2_rebin')
c      call pwhgemptyduphisto('a0-rebin','a0-norm')
c      call pwhgoperatehisto('/','a0-rebin','V_pt2_rebin','a0-norm')
c      call pwhgtopout("pwg-NLO-1")
c      end
cc This is to avoid having to link to much stuff
c      subroutine pwhg_exit(int)
c      integer int
c      call exit(int)
c      end
c



      subroutine pwhgloadhistos(fname)
c reads a file with output histograms and reloads them
      implicit none
      character *(*) fname
      include 'pwhg_bookhist-multi-new.h'
      integer maxbins,nchunk
      parameter (nchunk=100)
      real * 8, allocatable :: x1(:),x2(:),v(:),e(:),tmp(:)
      character * 100 string,line
      integer j,k,l,jh,iun
      integer indexhist
      procedure() :: indexhist
      call newunit(iun)
      open(unit=iun,file=fname,status='old',err=700)
      maxbins = nchunk
      allocate(x1(maxbins+1),x2(maxbins),v(maxbins),e(maxbins))
      do j=1,1000000
         read(10,'(a)',end=800) string
         if(string(1:1).eq.'#') then
c find name of histo
            string=adjustl(string(2:))
c     with the (optional) third argument set to .true., the last rather than the first occurrence
c     of the string is returned.
            k=index(string,' index ',.true.)
            if(k>0) then
               string = string(1:k-1)
            endif
c count the lines
            k=0
            read(10,'(a)') line
            do while(line.ne.' ')
               k=k+1
               if(k.gt.maxbins) call resizearrs
               read(line,*) x1(k),x2(k),v(k),e(k)
               read(10,'(a)') line
            enddo
            x1(k+1)=x2(k)
            call bookup(trim(string),k,x1)
            jh = indexhist(string)
            do l=1,k
               hist_ptr(jh)%yhistarr2(1,l) = v(l)
               hist_ptr(jh)%errhistarr2(1,l) = e(l)
            enddo
         endif
      enddo
      write(*,*) ' pwhgloadhistos: more than 10^6 lines in '//fname
      write(*,*) ' something wrong ... exiting ...'
      call exit(-1)
 700  continue
      write(*,*) ' pwhgloadhistos: cannot open '//fname
      write(*,*) ' exiting ...'
      call exit(-1)
 800  continue
      close(iun)
      return
      contains
      subroutine resizearrs
         allocate(tmp(maxbins+1))
         tmp(1:maxbins) = x1(1:maxbins)
         deallocate(x1)
         allocate(x1(maxbins+1+nchunk))
         x1(1:maxbins) = tmp(1:maxbins)
         tmp(1:maxbins) = x2(1:maxbins)
         deallocate(x2)
         allocate(x2(maxbins+nchunk))
         x2(1:maxbins) = tmp(1:maxbins)
         tmp(1:maxbins) = v(1:maxbins)
         deallocate(v)
         allocate(v(maxbins+nchunk))
         v(1:maxbins) = tmp(1:maxbins)
         tmp(1:maxbins) = e(1:maxbins)
         deallocate(e)
         allocate(e(maxbins+nchunk))
         e(1:maxbins) = tmp(1:maxbins)
         deallocate(tmp)
         maxbins = maxbins + nchunk
      end subroutine resizearrs

      end

      subroutine pwhgemptyduphisto(str1,str2)
c makes an empty duplicate of histogram str1 into histogram str2
      implicit none
      include 'pwhg_bookhist-multi-new.h'
      character * (*) str1,str2
      integer j,k,n,ind1,ind2
      do j=1,jhist
         if(hist_ptr(j)%id.eq.str1) then
            ind1=j
            goto 11
         endif
      enddo
      write(*,*) ' pwhgemptyduphisto: hist. ',trim(str1),
     1           'not found! exiting ...'      
 11   continue
      do j=1,jhist
         if(hist_ptr(j)%id.eq.str2) then
            write(*,*) ' pwhgemptyduphisto: hist. ',trim(str2),
     1           'already booked! exiting ...'
            call exit(-1)
         endif
      enddo
      call bookup(str2,hist_ptr(ind1)%nbins,hist_ptr(ind1)%xhistarr)
      end

      
      subroutine pwhgoperatehisto(op,str1,str2,str3)
c makes an empty duplicate of histogram str1 into histogram str2
      implicit none
      include 'pwhg_bookhist-multi-new.h'
      character * (*) op,str1,str2,str3
      integer j,k,n,ind1,ind2,ind3
      ind1=0
      ind2=0
      ind3=0
      do j=1,jhist
         if(hist_ptr(j)%id.eq.str1) then
            ind1=j
         elseif(hist_ptr(j)%id.eq.str2) then
            ind2=j
         elseif(hist_ptr(j)%id.eq.str3) then
            ind3=j
         endif
      enddo
      if(ind1.eq.0) then
         write(*,*) ' pwhgoperatehisto: histogram ',trim(str1),
     1        ' not found, exiting ...'
         call exit(-1)
      endif
      if(ind2.eq.0) then
         write(*,*) ' pwhgoperatehisto: histogram ',trim(str2),
     1        ' not found, exiting ...'
         call exit(-1)
      endif
      if(ind3.eq.0) then
         write(*,*) ' pwhgoperatehisto: histogram ',trim(str3),
     1        ' not found, exiting ...'
         call exit(-1)
      endif
      n = min(hist_ptr(ind1)%nbins,hist_ptr(ind2)%nbins,hist_ptr(ind3)%nbins)
c check that x values are compatible
      do k=1,hist_ptr(ind1)%nbins
         if(hist_ptr(ind1)%xhistarr(k).ne.hist_ptr(ind2)%xhistarr(k).or.
     1        hist_ptr(ind1)%xhistarr(k).ne.hist_ptr(ind3)%xhistarr(k)) then
            write(*,*)' pwhgoperatehisto: x array no compatible among ',
     1           trim(str1),',',trim(str2),',',trim(str3)
            write(*,*) ' exiting'
            call exit(-1)
         endif
      enddo

      if(op.eq.'/') then
         do k=1,hist_ptr(ind1)%nbins
            hist_ptr(ind3)%yhistarr2(1,k) =  
     1           hist_ptr(ind1)%yhistarr2(1,k) / hist_ptr(ind2)%yhistarr2(1,k) 
            hist_ptr(ind3)%errhistarr2(1,k) = hist_ptr(ind3)%yhistarr2(1,k) *
     1           sqrt((hist_ptr(ind1)%errhistarr2(1,k)/hist_ptr(ind1)%yhistarr2(1,k))**2
     2               +(hist_ptr(ind2)%errhistarr2(1,k)/hist_ptr(ind2)%yhistarr2(1,k))**2)
         enddo
      else
         write(*,*) ' pwhgoperatehisto: unknown operator ',op
         write(*,*) ' exiting'
         call exit(-1)
      endif
      end

      subroutine pwhgoperatehisto1(op,val,str)
c applies op val to the histogram in str.
      implicit none
      include 'pwhg_bookhist-multi-new.h'
      character * (*) op,str
      real * 8 val
      integer j,k,n,ind
      ind=0
      do j=1,jhist
         if(hist_ptr(j)%id.eq.str) then
            ind=j
         endif
      enddo
      if(ind.eq.0) then
         write(*,*) ' pwhgoperatehisto: histogram ',trim(str),
     1        ' not found, exiting ...'
         call exit(-1)
      endif
      if(op.eq.'/') then
         do k=1,hist_ptr(ind)%nbins
            hist_ptr(ind)%yhistarr2(1,k)   = hist_ptr(ind)%yhistarr2(1,k) / val
            hist_ptr(ind)%errhistarr2(1,k) = hist_ptr(ind)%errhistarr2(1,k) / val
         enddo
      else
         write(*,*) ' pwhgoperatehisto: unknown operator ',op
         write(*,*) ' exiting'
         call exit(-1)
      endif
      end      
      
      subroutine pwhgrebin(str1,str2)
c makes an empty duplicate of histogram str1 into histogram str2
      implicit none
      include 'pwhg_bookhist-multi-new.h'
      character * (*) str1,str2
      integer j,k,n,ind1,ind2,istage
      real * 8 delk,delj
      logical pwhg_app
      ind1=0
      ind2=0
      do j=1,jhist
         if(hist_ptr(j)%id.eq.str1) then
            ind1=j
         elseif(hist_ptr(j)%id.eq.str2) then
            ind2=j
         endif
      enddo
      if(ind1.eq.0) then
         write(*,*) ' pwhgoperatehisto: histogram ',trim(str1),
     1        ' not found, exiting ...'
         call exit(-1)
      endif
      if(ind2.eq.0) then
         write(*,*) ' pwhgoperatehisto: histogram ',trim(str2),
     1        ' not found, exiting ...'
         call exit(-1)
      endif
c check that x values are compatible
      j = 1
      istage = 0
      do k=1,hist_ptr(ind1)%nbins
         if(istage.eq.0) then
            if(pwhg_app('le',hist_ptr(ind1)%xhistarr(k+1),hist_ptr(ind2)%xhistarr(j))) cycle
            if(pwhg_app('ne',hist_ptr(ind1)%xhistarr(k),hist_ptr(ind2)%xhistarr(j)))
     1           goto 888
            istage = 1
         endif
         if(.not.(pwhg_app('ge',hist_ptr(ind1)%xhistarr(k),hist_ptr(ind2)%xhistarr(j)).and.
     1      pwhg_app('le',hist_ptr(ind1)%xhistarr(k+1),hist_ptr(ind2)%xhistarr(j+1)))) then
            j=j+1
            if(j.gt.hist_ptr(ind2)%nbins) then
               exit
            else
               if(.not.(pwhg_app('ge',hist_ptr(ind1)%xhistarr(k),hist_ptr(ind2)%xhistarr(j))
     1   .and. pwhg_app('le',hist_ptr(ind1)%xhistarr(k+1),hist_ptr(ind2)%xhistarr(j+1))))
     2              then
                  goto 888
               endif
            endif
         endif
         delk = hist_ptr(ind1)%xhistarr(k+1)-hist_ptr(ind1)%xhistarr(k)
         delj = hist_ptr(ind2)%xhistarr(j+1)-hist_ptr(ind2)%xhistarr(j)
         hist_ptr(ind2)%yhistarr2(1,j) = hist_ptr(ind2)%yhistarr2(1,j) +
     1        hist_ptr(ind1)%yhistarr2(1,k)*delk/delj
         hist_ptr(ind2)%errhistarr2(1,j) = sqrt(hist_ptr(ind2)%errhistarr2(1,j)**2+
     1        (hist_ptr(ind1)%errhistarr2(1,k)*delk/delj)**2)
      enddo
      return
 888  continue
      write(*,*) ' pwhgrebin: incompatible binning between ',
     1     trim(str1),' and ',trim(str2),', exiting...'
      call exit(-1)
      end


      logical function pwhg_app(op,a,b)
      implicit none
      character * 2 op
      real * 8 a,b
      real * 8 r
      if(a.eq.0.and.b.eq.0) then
         r=0
      else
         r = (a-b)/(abs(a)+abs(b))
         if(abs(r).lt.1d-6) r = 0
      endif
      if(op.eq.'ge') then
         pwhg_app = r .ge. 0
      elseif(op.eq.'le') then
         pwhg_app = r .le. 0
      elseif(op.eq.'eq') then
         pwhg_app = r .eq. 0
      elseif(op.eq.'ne') then
         pwhg_app = r .ne. 0
      elseif(op.eq.'lt') then
         pwhg_app = r .lt. 0
      elseif(op.eq.'gt') then
         pwhg_app = r .gt. 0
      endif
      end


      subroutine pwhginteghisto(dir,str1,str2)
c Integrates histogram str1 into histogram str2 (that is created).
c If dir=1, integrates from the low bin forward,
c if dir=-1, it integrates from the high bin backward
      implicit none
      include 'pwhg_bookhist-multi-new.h'
      character * (*) str1,str2
      integer dir,j,n,ind1,ind2
      real * 8 del
      ind1=0
      do j=1,jhist
         if(hist_ptr(j)%id.eq.str1) then
            ind1=j
            exit
         endif
      enddo
      if(ind1.eq.0) then
         write(*,*) ' pwhginteghisto: histogram ',trim(str1),
     1        ' not found, exiting ...'
         call exit(-1)
      endif
      call pwhgemptyduphisto(str1,str2)
      ind2=0
      do j=1,jhist
         if(hist_ptr(j)%id.eq.str2) then
            ind2=j
            exit
         endif
      enddo
      n = hist_ptr(ind1)%nbins
      if(dir.eq.-1) then
c each entry contains the integral from xhistarr(j) up to infinity (including
c the overflow). When printed with pwhgtopout, the first column contains
c the x, the third one the corresponding value of the integral.
         hist_ptr(ind2)%yhistarr2(1,n+1) = hist_ptr(ind1)%yhistarr2(1,n+1) 
         hist_ptr(ind2)%errhistarr2(1,n+1) = hist_ptr(ind1)%errhistarr2(1,n+1) 
         do j=n,1,-1
            del = hist_ptr(ind1)%xhistarr(j+1)-hist_ptr(ind1)%xhistarr(j)
            hist_ptr(ind2)%yhistarr2(1,j) =hist_ptr(ind2)% yhistarr2(1,j+1)
     1           +hist_ptr(ind1)%yhistarr2(1,j)*del
            hist_ptr(ind2)%errhistarr2(1,j) = sqrt(
     1           hist_ptr(ind2)%errhistarr2(1,j+1)**2
     2           +hist_ptr(ind1)%errhistarr2(1,j)**2*del**2)
         enddo
      elseif(dir.eq.1) then
c each entry j contains the integral up to xhistarr(j);
c thus j=1 contains only the underflow. When printed with pwhgtopout,
c the first column contains the x, the third one the corresponding
c value of the integral.
         hist_ptr(ind2)%yhistarr2(1,1) = hist_ptr(ind1)%yhistarr2(1,0) 
         hist_ptr(ind2)%errhistarr2(1,1) = hist_ptr(ind1)%errhistarr2(1,0) 
         do j=2,n
            del = hist_ptr(ind1)%xhistarr(j)-hist_ptr(ind1)%xhistarr(j-1)
            hist_ptr(ind2)%yhistarr2(1,j) = hist_ptr(ind2)%yhistarr2(1,j-1)
     1           +hist_ptr(ind1)%yhistarr2(1,j-1)*del
            hist_ptr(ind2)%errhistarr2(1,j) = sqrt(
     1           hist_ptr(ind2)%errhistarr2(1,j-1)**2
     2           +hist_ptr(ind1)%errhistarr2(1,j-1)**2*del**2)
         enddo
      else
         write(*,*) ' pwhginteghisto: error: first argument should be'
         write(*,*) ' either 1 or -1; got ',dir
         call exit(-1)
      endif
      end

