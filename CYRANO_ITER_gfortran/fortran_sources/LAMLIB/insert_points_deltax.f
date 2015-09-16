
      real*8 function insert_points_deltax(Nin, vecin, value, SIDE, range, Npts, delta, vecout)
      
	implicit none

      integer, intent(in) :: Nin, Npts
      real*8, intent(in)  :: vecin(Nin), value, delta, range
      real*8, intent(out) :: vecout(Nin+Npts)
      character(1), intent(in) :: SIDE

	integer :: j, indx, find_nearest, ilow, ihig, iaux
	real*8 :: auxfine(Npts), delaux, wid, ddd
	
	integer isrchfge
	external isrchfge
	external find_nearest

c     Insert equally spaced Npts points in an original vector near point given by 'value'

c     MODIFIED from insert_points2.f: 
c     range is the actual total width in which to insert the points

c     Nin: size of original vector
c     vecin: origial vector
c     value: value near which the points will be added
c     SIDE: insert points to the left ('L') or to the right ('R') of value
c     range: range in which the new points will be inserted 
c     Npts: number of points to be inserted
c     delta: distance from the (sing.) point 'value', where to begin to insert the new points
c     vecout: output vector of size Nin+Npts with inserted points


c	BEGIN:
      auxfine(1:Npts)=0.d0
	vecout(1:Nin) = vecin(1:Nin)
	ddd = delta

c	Find first point >= value	
	ihig = isrchfge(Nin, vecin(1:Nin), 1, value)
	ilow = ihig-1 
	
	if(SIDE .eq. 'L')then   ! Insert points to the LEFT side of 'value'

		indx = ilow
	
	    if(indx .lt. 1)goto 4444
		
		! Create fine grid
c		wid = dabs(value-vecin(indx-range+1))
            wid = range !!!!!!!!!!!!!!!!!!!!!!!!
		delaux = wid / dfloat(Npts+1)
c	    if(delta.eq.0.d0 .or. delta>delaux)ddd=delaux

		do j = 1,Npts
			auxfine(Npts-j+1) = (value-delta) - dfloat(j-1) * delaux
	            if(auxfine(Npts-j+1)<0)print*,'oops:',value
		end do
		
		! Shift original points that are inside fine region
	      iaux = isrchfge(Nin, vecin(1:Nin), 1, value-delta-wid)
		vecout(iaux:indx) = value - delta - wid

	    ! Shift values Npts to the right
		vecout(indx+Npts:Nin+Npts) = vecout(indx:Nin)

	    ! Fill Npts empty points with fine grid
		vecout(indx+1:indx+Npts) = auxfine(1:Npts)

4444	continue

	else	! Insert points to the RIGHT side of 'value'

		indx = ihig

	    if(indx .ge. Nin)goto 5555

		! Create fine grid
c		wid = dabs(vecin(indx+range-1)-value)
	      wid = range !!!!!!!!!!!!!!!!!!!!!!!
		delaux = wid / dfloat(Npts+1)
c	    if(delta.eq.0.d0 .or. delta>delaux)ddd=delaux

		do j = 1,Npts
			auxfine(j) = (value+delta) + dfloat(j-1) * delaux
		end do

		! Shift original points that are inside fine region
		iaux = isrchfge(Nin, vecin(1:Nin), 1, value + delta + wid)
		vecout(indx:iaux-1) = value + delta + wid

		! Shift values Npts to the right
		vecout(indx+Npts:Nin+Npts) = vecout(indx:Nin)

	    ! Fill Npts empty points with fine grid
		vecout(indx:indx+Npts-1) = auxfine(1:Npts)

5555	continue

	end if ! left or right
      
	insert_points_deltax=0

      end