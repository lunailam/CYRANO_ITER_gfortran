      subroutine bafaga

      implicit none

c     Chooses gaussian abs. and weights for given order of integration.
c     Computes basis functions and deriv. at gaussian points
c             (reduced abscissa on [0,1], calls bafadn and old bafad)
c
c     argument:  ngauss, order of integration
c     result:    bafgn
c     both in comin2

      include 'pardim.copy'
      include 'comin2.copy'

      integer i, it
c     sizes must be changed in comin2 if higher degree than 5 introduced!
c     abscissa 0 is prepended to normal gaussian abscissae in aga (for tricks
c     at magnetic axis)

      write(603,*)'bafaga: ngauss,ngauss= ',ngauss,ngauss
      call gausaw(ngauss, aga(1), wga)
      aga(0) = 0.d0
      aga(ngauss+1) = 1.d0
      
      do i = 0, ngauss + 1
        do it = 1, ntyp
        call bafadn(it, aga(i), bafgn(1,0,it,i), bafsin(1,0,it,i))
        end do
      call bafad(aga(i), bafg(1,i))
      end do

      return
      end
