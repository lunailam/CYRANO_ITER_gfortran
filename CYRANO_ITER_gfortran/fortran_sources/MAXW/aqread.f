      subroutine aqread(aqp, vec, iblk, nblk, irea, iqueue, istatu)
      
      integer aqp, iblk, nblk, irea, iqueue, istatu
      double precision vec(*)

      include 'pardim.copy'

	integer i, j, irec
            
c     Normal READ; n.b. aqp has different use than on Cray.
      irec = iblk
      do i = 1, nblk
      read(unit=aqp, rec=irec, iostat=istatu)
     ;(vec(j), j = 1 + iobll * (i-1), i*iobll)
	  irec = irec + 1
      end do
      
      return
      end