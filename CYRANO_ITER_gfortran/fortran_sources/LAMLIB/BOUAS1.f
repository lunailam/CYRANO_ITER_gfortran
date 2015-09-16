      subroutine bouas1(rc, ni, npopi, axl, ipop1, ipop2, in)

      implicit none
      character*1 rc
      integer ni, npopi, ipop1(3), ipop2(3), in(6), m2
      complex*16 axl(3)      
c
c     Assembly list for first element, either at magnetic axis
c                                      or at inner metal shell.

c     Input:
C     RC: 'R' or 'C', for row or column.

C     Output:
C     NI: number of dof kept
C     IN(NI): list of indices of dof kept
C     NPOPI: number of equations (or unknowns) to eliminate.
C     AXL: row (or column) multiplication factor.
C     IPOP1: list of destination indices for elimination of equations
C            (or unknowns).
C     IPOP2: list of origin indices (equations or unknowns to be eliminated)

      include 'pardim.copy'
      include 'comgeo.copy'
      include 'compla.copy'
      include 'comsub.copy'
      include 'commod.copy'
      include 'comfin.copy'
      include 'comphy.copy'

      character*3 eletyp

      logical hot
      
      integer mode, icolo, icpblo, idllo, ibulo

      if(rc.eq.'R' .or. rc.eq.'r')then
      mode = m
      else if(rc.eq.'C' .or. rc.eq.'c')then
      mode = mpk
      end if
      
      eletyp = styp(1)
      icolo = iconn(1)
      ibulo = ibub(1)
      idllo = iddl(1)
      icpblo = icolo + ibulo

      hot = .not.vacuum(1) .and. .not.coldpl(1) .and. flrops(1)

      ni = 0
      npopi = 0
                        if(crown)then
c                       -------------
      if(eletyp.eq.'HEC')then
c     Eliminate E- in favour of E+: (equation E- = E+)
      npopi = 1
      axl(1) = cun
      ipop1(1) = 1
      ipop2(1) = 3
      ni = 4
      if(ndof.lt.6)ni = 3
c     Keep E+, E+', E-', E//':
      in(1) = 1
      in(2) = 2
      in(3) = 4
      in(4) = 6

      else if(eletyp.eq.'M23')then
      ni = 4
      if(ndof.lt.6)ni = 3
c     Keep Erho, Etheta', Ephi', Erho at element centre (bubble dof):
      in(1) = 1
      in(2) = 3
      in(3) = 5
      in(ni) = icolo+1

      else if(eletyp.eq.'CAR')then
	print *, 'New CAR element: must implement convolution for metallic condition'
	stop
      end if
                        else  ! Magnetic axis:
c                       ----
      if(eletyp.eq.'HEC')then

c     Odd modes, ||>1: all components and first derivative vanish:
      if(iabs(mode).gt.1.and. (mode/2)*2.ne.mode) return

      	if(hot)then
	
      		if(mode.eq.0)then  
      		ni = 3
      		if(ndof.lt.5)ni = 2
c           Keep  E+', E-', E//:
      		in(1) = 2
      		in(2) = 4
      		in(3) = 5

      		else if(iabs(mode).eq.1)then  
      		ni = 1
      		if(ndof.lt.6)ni = 0
c           Keep  E//':
      		in(1) = 6

      		else if(mode.lt.-1)then
      		ni = 1
c           Keep  E+':
      		in(1) = 2
c           Eliminate E-' in favour of E+': (equation (m+2) E+' + (m-2) E-' = 0)
      		npopi = 1
      		axl(1) = -(mode+2.d0)/(mode-2.d0)
      		ipop1(1) = 2
      		ipop2(1) = 4

      		else if(mode.gt.1)then
      		ni = 1
C           Keep  E-':
      		in(1) = 4
C           Eliminate E+' in favour of E-': (equation (m+2) E+' + (m-2) E-' = 0)
      		npopi = 1
      		axl(1) = -(mode-2.d0)/(mode+2.d0)
      		ipop1(1) = 4
      		ipop2(1) = 2
      		end if

      	else
	
      		if(mode.eq.0)then
      		ni = 3
      		if(ndof.lt.5)ni = 2
c           Keep  E+', E-', E//:
      		in(1) = 2
      		in(2) = 4
      		in(3) = 5

      		else if(mode.eq.1)then
      		ni = 2
      		if(ndof.lt.6)ni = 1
c           Keep  E-, E//':
      		in(1) = 3
      		in(2) = 6

      		else if(mode.eq.-1)then
      		ni = 2
      		if(ndof.lt.6)ni = 1
c           Keep  E+, E//':
      		in(1) = 1
      		in(2) = 6

      		else if(mode.lt.-1)then
c     		check!
      		ni = 1
c           Keep  E+':
      		in(1) = 2
c           Eliminate E-' in favour of E+': (equation (m+2) E+' + (m-2) E-' = 0)
      		npopi = 1
      		axl(1) = -(mode+2.d0)/(mode-2.d0)
      		ipop1(1) = 2
      		ipop2(1) = 4

      		else if(mode.gt.1)then
c     		check!
      		ni = 1
c           Keep  E-':
      		in(1) = 4
c           Eliminate E+' in favour of E-': (equation (m+2) E+' + (m-2) E-' = 0)
      		npopi = 1
      		axl(1) = -(mode-2.d0)/(mode+2.d0)
      		ipop1(1) = 4
      		ipop2(1) = 2
      		end if
	
      	end if

      else if(eletyp.eq.'M23')then

      	if(mode.eq.0)then  
      	ni = 3
      	if(ndof.lt.5)ni = 2
c       keep etheta', ephi, erho bubble:
      	in(1) = 3
      	in(2) = 4
      	in(ni) = icolo+1

      	else if(iabs(mode).eq.1)then 
      	ni = 2
      	if(ndof.lt.6)ni = 1
c       Keep Erho, Ephi':
      	in(1) = 1
      	in(2) = 5
      	npopi = 3
c       Eliminate Etheta in favour of Erho: (equation Etheta = i*m*Erho)
      		if(rc.eq.'R' .or. rc.eq.'r')then
      		axl(1) = - ci * mode
      		else if(rc.eq.'C' .or. rc.eq.'c')then
      		axl(1) =  ci * mode
      		end if
      	ipop1(1) = 1
      	ipop2(1) = 2
C       Eliminate Erho bubble in favour of Erho (left): (eq. Erho'(0) = 0)
      	axl(2) = 0.75d0
      	ipop1(2) = 1
      	ipop2(2) = icolo+1
c       ... and in favour of Erho (right): (same eq.)
      	axl(3) = 0.25d0
      	ipop1(3) = icpblo+1
      	ipop2(3) = icolo+1

      	else if(iabs(mode).gt.1)then
      	m2 = mode / 2
      	     if(2*m2 .eq. mode)then
c            Even modes: keep Etheta', Erho bubble:
      	     ni = 2
      	     in(1) = 3
      	     in(2) = icolo+1
      	     else
c     	     Odd modes:
C check!
      	     npopi = 2
c            Eliminate Erho bubble in favour of Erho (left): (eq. Erho'(0) = 0)
      	     axl(1) = 0.75d0
      	     ipop1(1) = 1
      	     ipop2(1) = icolo+1
c            ... and in favour of Erho (right): (same eq.)
      	     axl(2) = 0.25d0
      	     ipop1(2) = icpblo+1
      	     ipop2(2) = icolo+1
      	     end if
      	end if


      else if(eletyp.eq.'CAR')then

      	if(mode.eq.0)then  
      	ni = 3
      	if(ndof.lt.5)ni = 2
c       keep eR, eY, ephi:
      	in(1) = 1
      	in(2) = 3
      	in(3) = 5

      	else if((mode/2)*2.ne.mode)then 
c       keep eR', eY', ephi':
      	ni = 3
      	if(ndof.lt.6)ni = 2
      	in(1) = 2
      	in(2) = 4
      	in(3) = 6
      	end if

      end if  ! Element type selection
                        end if  ! To be or not to be a magnetic axis 
c                       ------

      return
      end
