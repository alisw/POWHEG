      subroutine sigd3term(res)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_br.h'
      include 'pwhg_st.h'
      include 'pwhg_flg.h'
      include 'pwhg_pdf.h'
      include 'minnlo_flg.h'
c Cross section in the underlying Born of the Born
      real * 8 uubsigma
      real * 8 pdf1(-pdf_nparton:pdf_nparton),
     1         pdf2(-pdf_nparton:pdf_nparton)
      real * 8 rescfac,tot0
      real * 8 res0(flst_nborn), res(flst_nborn)
      real * 8 d3terms,uubjakob
      common/d3terms/ d3terms
      integer j
      real * 8 mufact2(flst_nborn)
      integer iborn,counter
      real * 8 my_mufact2
      common/my_iborn/my_mufact2,iborn
      real *8 tmp, bornAP
      logical test_rapidity
      real *8 powheginput
      
      test_rapidity=powheginput("#test_rapidity").eq.1


      tot0=0d0
      d3terms = 0d0
      do j=1,flst_nborn
c     Here in fact setlocalscales works only on its first call;
c     in all remaining ones it returns the previous values, and
c     the d3term is set on the first call
         call setlocalscales(j,3,rescfac)
c     remember the scales used for each subprocess;
c     they are needed in the flg_distribute_by_ub case
         mufact2(j)=st_mufact2
         call pdfcall(1,kn_xb1,pdf1)
         call pdfcall(2,kn_xb2,pdf2)

         if(flg_distribute_by_ub) then
            res0(j)=br_born(j) *
     1           pdf1(flst_born(1,j))*pdf2(flst_born(2,j))*kn_jacborn
            tot0=tot0+res0(j)
         endif
      enddo
      if(flg_minnlo) then
         if(flg_distribute_by_ub) then
c     arbitrarily choose mufact2(1) to fill my_mufact2 (needed by
c     evaluubjakob), as in any case, for now, all the entries of mufact2
c     are the same. For sure this can be coded better.
            my_mufact2 = mufact2(1)
c     evaluate uubjakob.  
c     If one wants to include d3 in only some
c     partonic channel, the "iborn" variable (see common block above) can ne used to pass this information to the subroutines in ProjKinematics.f
            if(test_rapidity) then
               uubjakob=1d0
            else
               call evaluubjakob(uubjakob)
            endif
            if (uubjakob .eq. -1d10) then
               counter=counter+1
               write(*,*) 'ERROR: uubjakob not correct, counter:', counter
            endif
c     uubjakob includes, in the denominator, the sum over all underlying
c     Born (pp->Hj) channels.  Now supply the d3 terms to each channel,
c     accordingly to its underlying Born matrix element.
            do j=1,flst_nborn
               if (res0(j) .ne. 0) then
                  if (uubjakob .ne. -1d10) then
                     if(flg_distribute_by_ub_AP) then
                        call setbornAPfrombrphsp(flst_born(1,j),flst_born(2,j),bornAP)
                        res(j) = res(j) + d3terms * uubjakob * res0(j)/br_born(j)*bornAP
                     else
                        res(j) = res(j) + d3terms * uubjakob * res0(j)
                     endif
                  endif
               endif
               if (res(j) > 1d-4) then
                  print*, "WARNING:"
                  print*, "res(j) too large; usually < ~10^-4:", res(j)
                  print*, "setting res(j)=0 for this event/flavour"
                  print*, d3terms, uubjakob, res0(j)
                  res(j) = 0d0
               endif
               if ((res(j)+1 .eq. res(j)) .or. (res(j)-1 .eq. res(j))) then
                  print*, "WARNING:"
                  print*, "res(j) NaN or Infinity", res(j)
                  print*, "setting res(j)=0 for this event/flavour"
                  print*, d3terms, uubjakob, res0(j)
                  res(j) = 0d0
               endif
            enddo
         else
c$$$c     start dbugging instruction to limit the cos theta=y range
c$$$c     first find vector boson rapidity
c$$$            call find_j_lims(tmp)
c$$$c     end debugging
            call evaluubjakob(uubjakob)
            res(:) = res(:) + res0(:)*(d3terms*tmp/tot0*kn_jacborn*uubjakob)
         endif
      endif
      end

c$$$c only for debugging (REMOVE?) >>>
c$$$c     start debugging
c$$$      subroutine find_j_lims(tmp)
c$$$      implicit none
c$$$      real * 8 tmp
c$$$      include 'nlegborn.h'
c$$$      include 'pwhg_kn.h'
c$$$      real * 8 pv(0:3),pp(0:3),yp,yv
c$$$      real *8, save :: deltaymax=-1d0
c$$$      real *8 deltaymaxvalue
c$$$      common/cdeltaymax/deltaymaxvalue
c$$$      real *8 powheginput
c$$$      
c$$$      if(deltaymax.eq.-1) then
c$$$         deltaymax=powheginput('#deltaymax')
c$$$      endif
c$$$      if(deltaymax == -1000000) then
c$$$         tmp=1d0
c$$$         deltaymaxvalue=log(1d10)
c$$$         return
c$$$      else
c$$$         deltaymaxvalue=deltaymax
c$$$      endif
c$$$
c$$$c     compute the V rapidity
c$$$      pv = kn_cmpborn(:,3) + kn_cmpborn(:,4)
c$$$      yv = 0.5d0 * log( (pv(0) + pv(3)) / (pv(0) - pv(3)))
c$$$      pp =  kn_cmpborn(:,5)
c$$$      yp = 0.5d0 * log( (pp(0) + pp(3)) / (pp(0) - pp(3)))
c$$$      if( abs(yv-yp) > deltaymax ) then
c$$$         tmp = 0
c$$$      else
c$$$         tmp = 1
c$$$      endif
c$$$      end
c$$$c <<< only for debugging (REMOVE?)
