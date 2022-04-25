C$$      implicit none
C$$      real * 8 x1b,x2b,fourPtsqOsh,res
C$$      x1b = -1
C$$      call  borndenomint(x1b,x2b,fourPtsqOsh,res)
C$$      end
C$$
c     Here we find the region in y where
c     csimax - 2 ptsq/(shat (1-y^2)) > 0.
c     This leads to
c      1 - max(g(y,x1sq),g(-y,x2sq)) - 2 ptsq/(shat (1-y^2)) > 0
c     where
c     g(y,xsq)=2*(1+y)*xsq/sqrt((1+xsq)^2*(1-y)^2+16*y*xsq+(1-y)*(1-xsq)).
c     We can rewrite it as
c     min(fp(y),fm(y)) > 2 ptsq/shat
c     where
c     fp(y)=(1 - g(y,x1sq) )*(1-y^2)
c     fm(y)=(1 - g(-y,x2sq))*(1-y^2)
c     We also introduce the derivatives difffp and difffm
c     The main program tests for various assumption on the behaviour of the functions
c     for random choice of the parameters.
c     The contained function findyrange(y1,y2) find the range of integration in y
c     at fixed x1sq,x2sq (the square of the barred x1,x2) and fourPtsqOsh (twice the pt^2 over s_hat).
c     In order to insert it in the code, turn the main program into a subroutine that calls findyrange(y1,y2)
c
c     At the moment the main program checks itself for consistency
      subroutine borndenomint(x1b,x2b,fourPtsqOsh,res)
      implicit none
      real * 8 x1b,x2b,fourPtsqOsh,res,res_bounds
      real * 8 x1sq,x2sq,yp,yd,ym,y,y0,y1,y2,tmp,tmp1,tmp2
      integer j,k
      character * 100 string
      real * 8 dgauss
      include 'brinclude.h'
      include 'minnlo_flg.h'
      logical theta_boundaries
      common/boundaries/theta_boundaries
      real *8 deltaymaxvalue
      common/cdeltaymax/deltaymaxvalue
      real *8 ymax,ymin,a,c,dgauss_accuracy,powheginput
      real * 8 accuracy_factor
      common/accuracy_factor/ accuracy_factor
      logical ini
      data ini/.true./
      save ini,dgauss_accuracy

      if(ini) then
         dgauss_accuracy = powheginput("#dgauss_accuracy")
         if(dgauss_accuracy.lt.0d0) dgauss_accuracy = 1d-3
         theta_boundaries=.false.
         ini=.false.
         write(*,*) '============================================='
         write(*,*) 'dgauss accuracy in borndenomint (def=0.001) = ',dgauss_accuracy
         write(*,*) 'theta_boundaries (def=F) = ',theta_boundaries
         write(*,*) '============================================='
      endif

c      dgauss_accuracy = 1d-3
c      theta_boundaries=.false.








c     First show that fp(y) is monotonically decreasing
c     (that also implies that fm(y) is monotonically increasing)
      if(x1b < 0) then
         do j=1,10000
            call random_number(x1sq)
            call  random_number(y1)
            y1 = -1 + 2*y1
            tmp1 = fp(y1)
            call  random_number(y2)
            y2 = -1 + 2*y2
            tmp2 = fp(y2)
            if((tmp2-tmp1)*(y2-y1) > 0) then
               write(*,*) 'fp is not decreasing!'
               call plot(fp,'dots')
            endif
         enddo

         write(*,*) 'Shown that fp is monotonically decreasing for ',
     1        j-1,' random calls'
c     In the following we plot fp, fm, fpm, their minimum, the
c     minimum we calculated, etc.
         do j=1,100
            write(11,*) ' set limits x -1 1 y 0 1'
            call random_number(x1sq)
            call random_number(x2sq)
            call random_number(fourPtsqOsh)
c     favour smaller values
            fourPtsqOsh = fourPtsqOsh**2/(1-fourPtsqOsh)
            
            call findzero(fpmfm,-1d0,1d0,yd,1d-8)
            if(yd == -2) then
               write(*,*) ' Error! yd=-2'
            endif
            call findyrange(y1,y2)
            call plot(fmp,'')
            call plot(fm,'dashes')
            call plot(fp,'dots')
            call plotv(yd,fmp(yd),'')
            call plotv(y1,fourPtsqOsh,'dotdash')
            call plotv(y2,fourPtsqOsh,'dotdash')
            call ploto(-1d0,1d0,fourPtsqOsh,'')
            if(.not. y1 == -2) then
               write(string,'(a,3(d10.4,1x),a)') 'title "Integral:',
     1              y1,y2,integral(y2)-integral(y1),'"'
               call plotstring(string)
            endif
            if(y1 /= -2) then
               tmp = dgauss(integrand,y1,y2,1d-6)
               write(*,*) (integral(y2)-integral(y1))/tmp
            endif
            write(11,*) 'newplot'
c     endif
         enddo
      else
         x1sq = x1b**2
         x2sq = x2b**2
         if ((.not.theta_boundaries).or.(.not.flg_distribute_by_ub)) then
            call findyrange(y1,y2)
            if(y1 == -2 .or. y2 == -2 .or. y1 > y2) then
               !write(*,*) ' borndenomint: no range for y!'
               !write(*,*) ' y1, y2 ',y1,y2

               call findzero(fpmfm,-1d0,1d0,yd,1d-8)
               if(yd == -2) then
                  write(*,*) ' Error! yd=-2'
               endif
               call findyrange(y1,y2)
               call plot(fmp,'')
               call plot(fm,'dashes')
               call plot(fp,'dots')
               call plotv(yd,fmp(yd),'')
               call plotv(y1,fourPtsqOsh,'dotdash')
               call plotv(y2,fourPtsqOsh,'dotdash')
               call ploto(-1d0,1d0,fourPtsqOsh,'')
               if(.not. y1 == -2) then
                  write(string,'(a,3(d10.4,1x),a)') 'title "Integral:',
     1                 y1,y2,integral(y2)-integral(y1),'"'
                  call plotstring(string)
               endif
               if(y1 /= -2) then
                  tmp = dgauss(integrand,y1,y2,1d-6)
                  write(*,*) (integral(y2)-integral(y1))/tmp
               endif
               
               write(*,*) ' exiting ...'
               call exit(-1) 
               return
            endif
         endif
         if (flg_distribute_by_ub) then
            if(theta_boundaries) then
c               res = dgauss(integrand_with_HJ,-1d0,1d0,1d-4)/1d8 !ER: using 1d-4 seems to make code a bit faster...
!               res = dgauss(integrand_with_HJ,-1d0,1d0,1d-2)/1d8 !ER: 1d-2 is the only way to have the DY code running reasonably quickly
               res = dgauss(integrand_with_HJ,-1d0,1d0,dgauss_accuracy)/1d8 !MW: with the AP approximation in the spreading we could try higher values again, eg 1d-3
            else
c$$$               fourPtsqOsh = 3.1891834005683616d-5 
c$$$               y1 = -0.99561386051320455d0
c$$$               y2 = 0.99999999916711657d0
c$$$               x1sq = 2.8394204673965019d0-9
c$$$               x2sq = 0.88653038768397741d0
c$$$               brkn_xb1 = 5.3286212732718225d-5
c$$$               brkn_xb2 = 0.94155742665223419
c$$$               fourPtsqOsh = 3.1891834005683616d-5 
c$$$               y1 = -0.99996148031407717d0
c$$$               y2 = 0.99999999832494535d0
c$$$               x1sq = 1.5986849383634349d-8
c$$$               x2sq =  0.31113270708920993d0 
c$$$               brkn_xb1 = 1.2643911334565087d-4
c$$$               brkn_xb2 = 0.55779270978492534d0
               accuracy_factor = 1d8
               res = dgauss(integrand_with_HJ,y1,y2,dgauss_accuracy)/accuracy_factor
               if(res .eq. 0d0) then
                  accuracy_factor = 1d0
                  res = dgauss(integrand_with_HJ,y1,y2,dgauss_accuracy)/accuracy_factor
               endif
               if(res .eq. 0d0) print*, fourPtsqOsh,y1,y2,x1sq,x2sq,brkn_xb1,brkn_xb2
               res_bounds=res
            endif
c$$$            if(abs(res/res_bounds-1d0).gt. 1d-3) then
c$$$               print*,"-----------------"
c$$$               print*, "res bounds = ",res_bounds
c$$$               print*, "res theta  = ",res
c$$$               print*, "ratio = ",res/res_bounds
c$$$            endif
         else

c$$$            A=sqrt(fourPtsqOsh/4d0)
c$$$
c$$$            c= exp(2*deltaymaxvalue)
c$$$
c$$$            ymin = ( (1+a**2)*(1-c**2) + 2* sqrt(a**2*(1+a**2)*(c-1)**2*c))
c$$$     1           /(a**2*(c-1)**2+(c+1)**2)
c$$$
c$$$            ymax = - ymin
c$$$
c$$$            if(ymax.lt.ymin) then
c$$$               print*, 'ymax<ymin! ',ymax,ymin
c$$$            endif
c$$$
c$$$            res = integral(min(y2,ymax)) - integral(max(y1,ymin))
!            res = dgauss(integrand,y1,y2,1d-6) // MW: checked to work equally well to compute d3terms (without weighting to HJ)

            res = integral(y2) - integral(y1)
!            res = dgauss(integrand,y1,y2,1d-6) // MW: checked to worke equally to compute d3terms (without weighting to HJ)
         endif
      endif

      contains

      subroutine plotstring(string)
      character *(*) string
      write(11,*) trim(string)
      end subroutine plotstring
      
      subroutine plot(f,string)
      real * 8 f,y
      character *(*) string
      integer k
      write(11,*) '# '//trim(string)
      do k=1,200
         y = 1-k/101d0
         write(11,*) y,f(y)
      enddo
      write(11,*)
      write(11,*)
c      write(11,*) ' join '//string 
      end

      subroutine plotv(y,val,string)
      real * 8 y,val
      character *(*) string
      write(11,*) '# '//trim(string)
      write(11,*) y,0
      write(11,*) y,val
      write(11,*)
      write(11,*)
c      write(11,*) ' join '//string
      end subroutine plotv

      subroutine ploto(y1,y2,val,string)
      real * 8 y1,y2,val
      character *(*) string
      write(11,*) '# '//trim(string)
      write(11,*) y1,val
      write(11,*) y2,val
      write(11,*)
      write(11,*)
c      write(11,*) ' join '//string
      end subroutine ploto


      subroutine findyrange(y1,y2)
      real * 8 y1,y2
      real * 8  yd,accuracy

      accuracy = 1d-10

      call findzero(fpmfm,-1d0,1d0,yd,accuracy)
      if(fmpd(yd)>0) then
         call findzero(fmpd,yd,-1d0,y1,accuracy)
         call findzero(fmpd,yd,1d0,y2,accuracy)
         if(y1 == -2 .or. y2 == -2) then
            write(*,*) ' findyrange: error!'
            call exit(-1)
         endif
!         y1 = y1 + accuracy
!         y2 = y2 - accuracy
      else
         y1 = -2
         y2 = -2
      endif

      end subroutine findyrange
      
      function fpmfm(y)
      real * 8 fpmfm,y
      fpmfm = fp(y) - fm(y)
      end function fpmfm
      
      function fp(y)
      real * 8 fp,y
      real * 8 xip,omxipoopy
      omxipoopy = 2*x1sq/
     1     (sqrt((1+x1sq)**2*(1-y)**2+16*y*x1sq)+(1-y)*(1-x1sq))
      xip = 1 - (1+y)*omxipoopy
      fp = xip**2/omxipoopy * (1-y)
      end function fp

      function fm(y)
      real * 8 fm,y
      real * 8 omximoomy,xim
      omximoomy = 2*x2sq/
     1     (sqrt((1+x2sq)**2*(1+y)**2-16*y*x2sq)+(1+y)*(1-x2sq))
      xim = 1 - (1-y) * omximoomy
      fm = xim**2/omximoomy * (1+y)
      end function fm

      function fmp(y)
      real * 8 fmp,y
      fmp = min(fp(y),fm(y))
      end function fmp

      function fmpd(y)
      real * 8 fmpd,y
      fmpd = min(fp(y),fm(y)) - fourPtsqOsh
      end function fmpd

      subroutine findzero(f,yli,yri,y0,accuracy)
c finds the zero of a decreasing function between yl,yr
      implicit none
      real * 8 f,yli,yri,y0
      real * 8 yl,yr,ytmp,rtmp,accuracy
      integer k

      if(f(yli) < 0 .or. f(yri) > 0) then
c     no solution
         y0=-2
      else
         yl = yli
         yr = yri
         do k=1,1000
            if(abs(yr-yl)<accuracy) exit
            ytmp = (yl+yr)/2
            rtmp = f(ytmp)
            if(rtmp > 0) then
               yl = ytmp
            elseif(rtmp < 0) then
               yr = ytmp
            else
               y0 = ytmp
               return
            endif
         enddo
         y0 = (yr+yl)/2
      endif
      end subroutine findzero


      function integrand(y)
      implicit none
      real * 8 integrand,y,k
      include 'pwhg_math.h'
      
      k = fourPtsqOsh
      integrand = 1/(2*(1-y**2)+k/2-sqrt(k*(1-y**2)+k**2/4))
      end function integrand

cccccccccccccccccccccccccccccccc
c     original integrand, just look below for the latest version of the
c     function with the same name 

c$$$      function integrand_with_HJ(y)
c$$$      implicit none
c$$$      real * 8 integrand_with_HJ,y
c$$$      real * 8 k,t,csi
c$$$      real * 8 random
c$$$      include 'brinclude.h'
c$$$      include 'pwhg_math.h'
c$$$      include 'pwhg_st.h'
c$$$      include 'pwhg_pdf.h'
c$$$      include 'pwhg_kn.h'
c$$$      integer bflav(br_nlegreal)
c$$$      real * 8 bornjk(br_nlegreal,br_nlegreal)
c$$$      real * 8 bmunu(0:3,0:3,br_nlegreal),born
c$$$      real * 8 pdf1(-pdf_nparton:pdf_nparton),
c$$$     1         pdf2(-pdf_nparton:pdf_nparton)
c$$$      integer iborn
c$$$      real * 8 my_mufact2
c$$$      common/my_iborn/my_mufact2,iborn
c$$$      real * 8 savemu2
c$$$
c$$$      k = fourPtsqOsh
c$$$      t = k/(1-y**2)
c$$$      csi = -t/2+sqrt(t**2/4+t)
c$$$      brkn_csi = csi
c$$$      brkn_y = y
c$$$      brkn_azi = 2*pi*random()
c$$$      call br_real_phsp_isr_rad
c$$$
c$$$      bflav(1:br_nlegreal) = flst_born(1:br_nlegreal,iborn)
c$$$      call setborn(brkn_preal,bflav,born,bornjk,bmunu)
c$$$      born=born/(2*kn_sborn)
c$$$
c$$$      savemu2 = st_mufact2
c$$$      st_mufact2 = my_mufact2
c$$$      !st_alpha
c$$$
c$$$      if (brkn_x1.gt.1 .or. brkn_x2.gt.1) then
c$$$         integrand_with_HJ = 0d0
c$$$         return
c$$$      endif
c$$$      call pdfcall(1,brkn_x1,pdf1)
c$$$      call pdfcall(2,brkn_x2,pdf2)
c$$$      st_mufact2 = savemu2
c$$$
c$$$      integrand_with_HJ = 1/(2*(1-y**2)+k/2-sqrt(k*(1-y**2)+k**2/4)) * born
c$$$     1     * pdf1(flst_born(1,iborn)) * pdf2(flst_born(2,iborn)) * 1d8
c$$$
c$$$      end function integrand_with_HJ

      function integrand_with_HJ(y)
      implicit none
      real * 8 integrand_with_HJ,y
      real * 8 k,t,csi
      real * 8 random
      include 'brinclude.h'
      include 'pwhg_math.h'
      include 'pwhg_st.h'
      include 'pwhg_pdf.h'
      include 'pwhg_kn.h'
      logical theta_boundaries
      common/boundaries/theta_boundaries

      integer bflav(br_nlegreal)
      real * 8 bornjk(br_nlegreal,br_nlegreal)
      real * 8 bmunu(0:3,0:3,br_nlegreal),born
      real * 8 pdf1(-pdf_nparton:pdf_nparton),
     1         pdf2(-pdf_nparton:pdf_nparton)
      integer iborn
      real * 8 my_mufact2
      common/my_iborn/my_mufact2,iborn
      real * 8 accuracy_factor
      common/accuracy_factor/ accuracy_factor
      real * 8 savemu2, xim, xip, tmax, rndazi
      integer fb
      real *8 born_times_pdf

      k = fourPtsqOsh
      t = k/(1-y**2)
      csi = -t/2+sqrt(t**2/4+t)

      if(theta_boundaries) then
         xim = 1d0 - 2*(1-y)*x2sq/(sqrt((1+x2sq)**2*(1+y)**2-16*y*x2sq)+(1+y)*(1-x2sq))
         xip = 1d0 - 2*(1+y)*x1sq/(sqrt((1+x1sq)**2*(1-y)**2+16*y*x1sq)+(1-y)*(1-x1sq))
         tmax = min(xim**2/(1-xim),xip**2/(1-xip))
         
         if ((tmax*(1-y**2)-k).lt.0) then 
            integrand_with_HJ = 0d0!1d-5
            return
         endif
      endif

      brkn_csi = csi
      brkn_y = y
      call random_number(rndazi)
      brkn_azi = 2*pi*rndazi
      call br_real_phsp_isr_rad


      if (brkn_x1.gt.1 .or. brkn_x2.gt.1) then
c     check the condition always, although it should never enter here if
c     theta_boundaries is true
         integrand_with_HJ = 0d0
         return
      endif

      savemu2 = st_mufact2
      st_mufact2 = my_mufact2
      !st_alpha

      call pdfcall(1,brkn_x1,pdf1)
      call pdfcall(2,brkn_x2,pdf2)
      st_mufact2 = savemu2

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ONLY COMPUTE BORN FOR FLAVOUR STRUCTURE ASSOCIATED TO iborn
c$$$      bflav(1:br_nlegreal) = flst_born(1:br_nlegreal,iborn)
c$$$      call setborn(brkn_preal,bflav,born,bornjk,bmunu)
c$$$      born=born/(2*kn_sborn)
c$$$      integrand_with_HJ = 1/(2*(1-y**2)+k/2-sqrt(k*(1-y**2)+k**2/4)) * born
c$$$     1     * pdf1(flst_born(1,iborn)) * pdf2(flst_born(2,iborn)) * 1d8
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     COMPUTE THE SUM OVER ALL FLAVOUR STRUCTURE
      born_times_pdf=0d0
      do fb=1,flst_nborn
         bflav(1:br_nlegreal) = flst_born(1:br_nlegreal,fb)
c     print*, fb,bflav(1:br_nlegreal)
         if(flg_distribute_by_ub_AP) then
c     Approximate the ZJ using the AP approximation, see notes BJapprox.tm
            call setbornAP(csi,y,flst_born(1,fb),flst_born(2,fb),born)
         else
            call setborn(brkn_preal,bflav,born,bornjk,bmunu)
            born = born/(2*brkn_sreal) ! include flux here as it is also included in numerator
         endif
         born_times_pdf = born_times_pdf + pdf1(flst_born(1,fb)) * pdf2(flst_born(2,fb)) * born
      enddo

      integrand_with_HJ = 1/(2*(1-y**2)+k/2-sqrt(k*(1-y**2)+k**2/4)) * born_times_pdf * accuracy_factor
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      end function integrand_with_HJ

      function integral(y)
      real * 8 integral,y
      real * 8 ik
      ik = 1/fourPtsqOsh
      integral = log((1+y)/(1-y))/2
     1     +log((sqrt(ik*(1-y**2)+1d0/4)+1d0/2+2*ik*(1-y))
     2 / (sqrt(ik*(1-y**2)+1d0/4)+1d0/2+2*ik*(1+y)))/4
      end function integral
      end


      subroutine setbornAP(csi,y,fl1,fl2,born)
      implicit none
      real * 8 csi,y,born
      integer fl1,fl2
      real * 8, parameter :: ca=3d0,cf=4d0/3d0,tf=0.5d0
      integer flav
      common/flav_initial_state/ flav

      if(flav .eq. 0) then ! gg case (Higgs)
         if(fl1 .ne. 0) then ! qg initial state
            born = cf * (1d0+csi**2)/(1d0-csi)/(1d0-y)/csi
         elseif(fl2 .ne. 0) then ! gq initial state
            born = cf * (1d0+csi**2)/(1d0-csi)/(1d0+y)/csi
         else ! gg initial state
c           q qbar 1/(1-y)+1/(1+y)=2/(1-y**2)
            born = 2d0 * ca * ((1d0-csi)/csi + csi/(1d0-csi) + (1d0-csi)*csi) * 2d0/(1d0-y**2)/csi
            ! is there this stupid 2 or not in the regularized splitting function?!?
         endif
      elseif(flav .eq. 1) then ! qq case (Drell-Yan)
         if(fl1 .eq. 0) then ! gq initial state
            born = tf * (csi**2+(1d0-csi)**2)/(1d0-y)/csi
         elseif(fl2 .eq. 0) then ! qg initial state
            born = tf * (csi**2+(1d0-csi)**2)/(1d0+y)/csi
         else ! qq initial state
c           q qbar 1/(1-y)+1/(1+y)=2/(1-y**2)
            born = cf * (1d0+(1d0-csi)**2)/csi**2 * 2d0/(1d0-y**2)
         endif
      else
         write(*,*) 'ERROR: value of flavour in setbornAP wrong.'
         stop
      endif
c from the flux factor 1/(2*Ecm^2)=z/(2*Q^2)=(1-csi)/(2*Q^2)      
      born = born * (1-csi)
      end

      subroutine setbornAPfrombrphsp(fl1,fl2,born)
      implicit none
      integer fl1,fl2
      real * 8 born
      include 'nlegborn.h'
      include 'pwhg_kn.h'
      real * 8 sqrts,k0,csi,y
      sqrts = 2*kn_cmpborn(0,1)
      k0 = kn_cmpborn(0,nlegborn)
      csi = 2*k0/sqrts
      y = kn_cmpborn(3,nlegborn)/k0
      call setbornAP(csi,y,fl1,fl2,born)
      end
