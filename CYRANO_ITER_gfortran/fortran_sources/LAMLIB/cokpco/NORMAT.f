      subroutine normat(lv, lx)
      
      implicit none
      double precision lv, lx
      
c     Builds auxiliary matrices depending on finite element lengths
      
      include 'pardim.copy'
      include 'comgrd.copy'
      
      integer i, j, ii, jj
      
      norma1(1,1) = 1.d0
      norma1(1,2) = lx
      norma1(2,1) = lv
      norma1(2,2) = lx * lv

      norma2(1,1) = 1.d0 / lv
      norma2(1,2) = lx * norma2(1,1)
      norma2(2,1) = 1.d0
      norma2(2,2) = lx

      norma3(1,1) = 1.d0 / lx
      norma3(1,2) = 1.d0
      norma3(2,1) = lv * norma3(1,1)
      norma3(2,2) = lv

      norma4(1,1) = norma2(1,1) * norma3(1,1)
      norma4(1,2) = norma2(1,1)
      norma4(2,1) = norma3(1,1)
      norma4(2,2) = 1.d0
      
      do ii = 0, 2, 2
      do jj = 2-ii, 2, 2
      do i = 1, 2
      do j = 1, 2
      norma1(ii+i,jj+j) = norma1(i,j)
      norma2(ii+i,jj+j) = norma2(i,j)
      norma3(ii+i,jj+j) = norma3(i,j)
      norma4(ii+i,jj+j) = norma4(i,j)
      end do
      end do
      end do
      end do
      
      return
      end
