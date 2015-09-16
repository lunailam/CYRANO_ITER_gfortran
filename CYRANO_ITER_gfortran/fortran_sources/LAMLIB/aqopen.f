      subroutine aqopen(aqp, naqp, lunit, istatu)
      
      integer aqp, naqp, lunit, istatu, lrecl, lbl

      include 'pardim.copy'
      
	lrecl = iobll * wordlzz
	  
c 19/12/03: added buffered='yes'
cJAC      open(unit=lunit,access='direct',status='scratch',recl=lrecl,buffered='yes',iostat=istatu)
      open(unit=lunit, access='direct', status='scratch', recl=lrecl , iostat=istatu)


c 19/12/03:
cJAC	inquire(unit=lunit, blocksize=lbl)
	inquire(unit=lunit)
	write(603,1000)lunit, lbl, lrecl
 1000 format('Opening solver unit', i4, ', block size =', i6
     ;     , ', record length =', i6)      
      return
      end
