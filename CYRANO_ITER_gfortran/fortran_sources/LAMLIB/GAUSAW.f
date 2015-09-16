      subroutine gausaw(ngauss, aga, wga)

      implicit none
      integer ngauss
      double precision aga(ngauss), wga(ngauss)

c     Given the order ngauss, returns the gaussian abscissae and weights
c     on (0,1).
c     Implemented for 1 <= ngauss <= 5.

      integer i, jmax
c     sizes to change if higher degree:
      double precision gabs(5,3), gwe(5,3)

c     Abscissae and weights on (-1,1) from Abramowitz & Stegun:
      data gabs(1,1)/0.0d0/
     ;    ,gabs(2,1)/0.577350269189626d0/,gabs(2,2)/0.0d0/
     ;    ,gabs(3,1)/0.0d0/              ,gabs(3,2)/0.774596669241483d0/
     ;    ,gabs(4,1)/0.339981043584856d0/,gabs(4,2)/0.861136311594053d0/
     ;    ,gabs(5,1)/0.0d0/              ,gabs(5,2)/0.538469310105683d0/
     ;    ,gabs(5,3)/0.906179845938664d0/
      data gwe(1,1)/2.d0/
     ;    ,gwe(2,1)/1.d0/
     ;    ,gwe(3,1)/0.888888888888889d0/,gwe(3,2)/0.555555555555556d0/
     ;    ,gwe(4,1)/0.652145154862546d0/,gwe(4,2)/0.347854845137454d0/
     ;    ,gwe(5,1)/0.568888888888889d0/,gwe(5,2)/0.478628670499366d0/
     ;    ,gwe(5,3)/0.236926885056189d0/


      write(603,*)'gausaw: ngauss= ', ngauss

      if(ngauss.lt.1)stop 'error: ngauss must be >= 1'
      if(ngauss.gt.5)stop 'error: ngauss must be <= 5'

c     Abscissae and weights on interval (0,1):

      if(mod(ngauss,2).eq.0)then
      jmax = ngauss / 2
        do i = 1, jmax
        aga(i) = (1.d0 - gabs(ngauss,jmax+1-i)) * 0.5d0
        wga(i) = 0.5d0 * gwe(ngauss,jmax+1-i)
        aga(jmax+i) = (gabs(ngauss,i) + 1.d0) * 0.5d0
        wga(jmax+i) = 0.5d0 * gwe(ngauss,i)
        end do

      else
      jmax = (ngauss+1) / 2
        do i = 1, jmax
        aga(i) = (1.d0 - gabs(ngauss,jmax+1-i)) * 0.5d0
        wga(i) = 0.5d0 * gwe(ngauss,jmax+1-i)
        end do
        do i = 2, jmax
        aga(jmax+i-1) = (gabs(ngauss,i) + 1.d0) * 0.5d0
        wga(jmax+i-1) = 0.5d0 * gwe(ngauss,i)
        end do

      end if

      write(603,*)'Gauss abscissae and weights on (0,1):'
      write(603,*)'Order: ', ngauss
      write(603,*)(aga(i),i=1,ngauss)
      write(603,*)(wga(i),i=1,ngauss)

      return
      end
