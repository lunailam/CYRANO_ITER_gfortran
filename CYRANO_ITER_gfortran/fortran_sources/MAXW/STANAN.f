      double precision function stanan(alpha)

      implicit none
      double precision alpha, pi
 
C     Reduces angle alpha to ]-pi,pi]
 
      data pi/3.14159265359d0/

      stanan = dmod(alpha, 2.d0*pi)

      if(dabs(stanan).gt.pi)
     ;stanan = stanan - dsign(2.d0*pi, stanan)
 
      return
      end
 
