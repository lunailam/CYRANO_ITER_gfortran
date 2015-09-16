      subroutine aqwrite(AQP, VEC, IBLK, NBLK, IWRI, IQUEUE, ISTATU)

      implicit none
      integer aqp, iblk, nblk, iwri, iqueue, istatu
      double precision vec(*)

      INCLUDE 'pardim.copy'

      integer irec, i, j
	  
C     Normal WRITE:
	  IREC = IBLK
      DO I = 1, NBLK
	  WRITE(UNIT=aqp, REC=IREC, IOSTAT=ISTATU)
     ;     (VEC(J), J = 1 + IOBLL * (I-1), I*IOBLL)
	  IREC = IREC + 1
      END DO    
        
      return
      end
