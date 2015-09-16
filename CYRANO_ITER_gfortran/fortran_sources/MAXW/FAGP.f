      subroutine fagp(iev, jex, ivn, jxn, isp)
      
      implicit none
      integer iev, jex, ivn, jxn, isp

c     Computes f0, df0/dv, df0/dx at Gauss points
c     using values in F0i obtained by interpolation (earlier call to INTF0)
c     iev, iex: element indices
c     ivn, jxn: indices of 'top right' element node in mesh
c     isp: species index (among gdr species)
      
      include 'pardim.copy'
      include 'comgdr.copy'
       
      integer i, j, igv, igx, ivm1, jxm1
      
      double precision bn(4,4), t1, t2, t3
      
      call dset(5*5, 0.d0, f0gp, 1)
      call dset(5*5, 0.d0, f0vgp, 1)
      call dset(5*5, 0.d0, f0xgp, 1)
      
      ivm1 = ivn - 1
      jxm1 = jxn - 1
      
      bn(1,1) = F0i(jxm1, ivm1, 1, isp)
      bn(2,1) = F0i(jxm1, ivm1, 2, isp)
      bn(1,2) = F0i(jxm1, ivm1, 3, isp)
      bn(2,2) = F0i(jxm1, ivm1, 4, isp)
      bn(1,3) = F0i(jxn, ivm1, 1, isp)
      bn(2,3) = F0i(jxn, ivm1, 2, isp)
      bn(1,4) = F0i(jxn, ivm1, 3, isp)
      bn(2,4) = F0i(jxn, ivm1, 4, isp)
      bn(3,1) = F0i(jxm1, ivn, 1, isp)
      bn(4,1) = F0i(jxm1, ivn, 2, isp)
      bn(3,2) = F0i(jxm1, ivn, 3, isp)
      bn(4,2) = F0i(jxm1, ivn, 4, isp)
      bn(3,3) = F0i(jxn, ivn, 1, isp)
      bn(4,3) = F0i(jxn, ivn, 2, isp)
      bn(3,4) = F0i(jxn, ivn, 3, isp)
      bn(4,4) = F0i(jxn, ivn, 4, isp)
      
      do i = 1, 4
      do j = 1, 4
      t1 = bn(i,j) * norma1(i,j)
      t2 = bn(i,j) * norma2(i,j)
      t3 = bn(i,j) * norma3(i,j)
        do igv = 1, ngauv
         do igx = 1, ngaux
         f0gp(igv,igx)  = f0gp(igv,igx)  + t1 * bibj(i,j,igv,igx)
         f0vgp(igv,igx) = f0vgp(igv,igx) + t2 * bipbj(i,j,igv,igx)
         f0xgp(igv,igx) = f0xgp(igv,igx) + t3 * bibjp(i,j,igv,igx)
         end do
        end do
      end do
      end do
      
      return
      end
