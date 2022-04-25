      subroutine mboost(m,vec,beta,vin,vout)
c     boosts the m vectors vin(0:3,m) into the vectors vout(0:3,m) (that can
c     be the same) in the direction of vec(3) (|vec|=1) with velocity
c     beta.  Lorents convention: (t,x,y,z).
      implicit none
      integer m
      real * 8 vec(3),beta,vin(0:3,m),vout(0:3,m)
      real * 8 betav,gamma
      real * 8 vdotb
      integer ipart,idim
      real * 8 tiny
      parameter (tiny=1d-14)
c      if (abs(beta).ge.1d0) then
c         write(*,*) '********** WARNING ***************'
c         write(*,*) 'mboost called with beta=',beta
c         write(*,*) '**********************************'
c     endif
      if(beta == 0) then
         vout=vin
         return
      endif
      if (beta.ge.1d0) then
         beta = 1-tiny
      elseif (beta.le.-1d0) then
         beta = -1+tiny
      endif
      gamma=1/sqrt(1-beta**2)
      do ipart=1,m
         vdotb=vin(1,ipart)*vec(1)
     #         +vin(2,ipart)*vec(2)+vin(3,ipart)*vec(3)
         do idim=1,3
            vout(idim,ipart)=vin(idim,ipart)
     #           +vec(idim)*((gamma-1)*vdotb
     #           +gamma*beta*vin(0,ipart))
         enddo
         vout(0,ipart)=gamma*(vin(0,ipart)+vdotb*beta)
      enddo
      end

      subroutine mrotate(dir,sinphi,cosphi,vec)
c Rotates vector vec counterclockwise around the direction
c dir (|dir|=1) with angle phi, given sin phi and cos phi.
      implicit none
      real * 8 sinphi,cosphi,dir(3),vec(3)
      real * 8 dircrossvec(3),dirdotvec
      integer i
      dircrossvec(1)=dir(2)*vec(3)-dir(3)*vec(2)
      dircrossvec(2)=dir(3)*vec(1)-dir(1)*vec(3)
      dircrossvec(3)=dir(1)*vec(2)-dir(2)*vec(1)
      dirdotvec=dir(1)*vec(1)+dir(2)*vec(2)+dir(3)*vec(3)
      do i=1,3
         vec(i)=vec(i)+sinphi*dircrossvec(i)
     #        -(1-cosphi)*(vec(i)-dir(i)*dirdotvec)
      enddo
      end

      subroutine rotate3tovec(v,vec)
c Rotate vec with the rotation that brings the positive
c third direction along v.
      implicit none
      real * 8 v(3),vec(3)
      real * 8 d(3),dir(3),cosphi,sinphi
      d = v/sqrt(v(1)**2+v(2)**2+v(3)**2)
      cosphi = d(3)
      dir(1) = d(2)
      dir(2) = -d(1)
      dir(3) = 0
      sinphi = sqrt(d(1)**2+d(2)**2)
      if(sinphi.lt.1d-9) return
      dir = - dir/sinphi
      call mrotate(dir,sinphi,cosphi,vec)
      end

      subroutine boost2reson(pres,nm,pin,pout)
      implicit none
      integer nm
      real * 8 pres(0:3),pin(0:3,nm),pout(0:3,nm)
      real * 8 vec(3),beta
      beta=sqrt(pres(1)**2+pres(2)**2+pres(3)**2)/pres(0)
      if(beta < 1d-30) then
         pout = pin
      else
         vec(1)=pres(1)/(beta*pres(0))
         vec(2)=pres(2)/(beta*pres(0))
         vec(3)=pres(3)/(beta*pres(0))
         call mboost(nm,vec,-beta,pin,pout)
      endif
      end

      subroutine boost2resoninv(pres,nm,pin,pout)
      implicit none
      integer nm
      real * 8 pres(0:3),pin(0:3,nm),pout(0:3,nm)
      real * 8 vec(3),beta
      beta=sqrt(pres(1)**2+pres(2)**2+pres(3)**2)/pres(0)
      vec(1)=pres(1)/(beta*pres(0))
      vec(2)=pres(2)/(beta*pres(0))
      vec(3)=pres(3)/(beta*pres(0))
      call mboost(nm,vec,beta,pin,pout)
      end
