      subroutine falp(jv, ix, irlow, jside, iev, jex, ivn, jxn, isp)
      
      implicit none
      integer jv, ix, irlow, jside, iev, jex, ivn, jxn, isp
c      double precision v, x

c     Computes f0, df0/dv, df0/dx, d2f0/dxdv at mesh point (jv,ix) of
c     interpolated grid, at radius of index irlow + jside - 1.
c     iev, iex: element indices
c     ivn, jxn: indices of 'top right' element node in mesh
c     isp: species index (among gdr species)

      include 'pardim.copy'
      include 'comgdr.copy'

      integer i, j, igv, igx, ivm1, jxm1, kr
      
      double precision bn(4,4), t1, t2, t3, t4
      
      kr = irlow + jside - 1

      f0l(kr,1) = 0.d0
      f0l(kr,2) = 0.d0
      f0l(kr,3) = 0.d0
      f0l(kr,4) = 0.d0
      
      ivm1 = ivn - 1
      jxm1 = jxn - 1
      
      bn(1,1) = F0(jxm1, ivm1, kr, 1, isp)
      bn(2,1) = F0(jxm1, ivm1, kr, 2, isp)
      bn(1,2) = F0(jxm1, ivm1, kr, 3, isp)
      bn(2,2) = F0(jxm1, ivm1, kr, 4, isp)
      bn(1,3) = F0(jxn, ivm1, kr, 1, isp)
      bn(2,3) = F0(jxn, ivm1, kr, 2, isp)
      bn(1,4) = F0(jxn, ivm1, kr, 3, isp)
      bn(2,4) = F0(jxn, ivm1, kr, 4, isp)
      bn(3,1) = F0(jxm1, ivn, kr, 1, isp)
      bn(4,1) = F0(jxm1, ivn, kr, 2, isp)
      bn(3,2) = F0(jxm1, ivn, kr, 3, isp)
      bn(4,2) = F0(jxm1, ivn, kr, 4, isp)
      bn(3,3) = F0(jxn, ivn, kr, 1, isp)
      bn(4,3) = F0(jxn, ivn, kr, 2, isp)
      bn(3,4) = F0(jxn, ivn, kr, 3, isp)
      bn(4,4) = F0(jxn, ivn, kr, 4, isp)
      
        do i = 1, 4
          do j = 1, 4
          t1 = bn(i,j) * norma1(i,j)
          t2 = bn(i,j) * norma2(i,j)
          t3 = bn(i,j) * norma3(i,j)
          t4 = bn(i,j) * norma4(i,j)

          f0l(kr,1) = f0l(kr,1) + t1 * bflnv(i,0,jv,jside) * bflnx(j,0,ix,jside)
          f0l(kr,2) = f0l(kr,2) + t2 * bflnv(i,1,jv,jside) * bflnx(j,0,ix,jside)
          f0l(kr,3) = f0l(kr,3) + t3 * bflnv(i,0,jv,jside) * bflnx(j,1,ix,jside)
          f0l(kr,4) = f0l(kr,4) + t4 * bflnv(i,1,jv,jside) * bflnx(j,1,ix,jside)
          end do
        end do
      
      return
      end
