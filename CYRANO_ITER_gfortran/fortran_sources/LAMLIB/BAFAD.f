      subroutine bafad(x, baf)

      implicit none
      double precision x, baf(21)
c
c     Computes basis functions and derivatives in local reduced coordinate (0,1).
c
      double precision x2, omx
      
      x2 = x * x
      omx = 1.d0 - x
      
c     C1 Hermite cubic element:

c 4 basis functions :
        baf(1) = ((2.d0 * x - 3.d0) * x2) + 1.d0
        baf(2) = ((x - 2.d0) * x + 1.d0) * x
        baf(3) = x2 * (3.d0 - 2.d0 * x)
        baf(4) = - x2 * omx
c 4 local derivatives:
        baf(5) = - 6.d0 * x * omx
        baf(6) = omx * (1.d0 - 3.d0 * x)
        baf(7) = - baf(5)
        baf(8) = x * (3.d0 * x - 2.d0)

c 4 local 2nd derivatives:
        baf(15) = 6.d0 * (2.d0 * x - 1.d0)
        baf(16) = 6.d0 * x - 4.d0
        baf(17) = - baf(15)
        baf(18) = 6.d0 * x - 2.d0

c     C0 Lagrange quadratic element:

c 3 basis functions:
        baf(9) = (2.d0 * x - 3.d0) * x + 1.d0
        baf(10) = 4.d0 * x * omx
        baf(11) = x * (2.d0 * x - 1.d0)

c 3 local derivatives:
        baf(12) = 4.d0 * x - 3.d0
        baf(13) = 4.d0 - 8.d0 * x
        baf(14) = 4.d0 * x - 1.d0

c 3 local 2nd derivatives:
        baf(19) = 4.d0
        baf(20) = - 8.d0
        baf(21) = 4.d0

        return
        end
