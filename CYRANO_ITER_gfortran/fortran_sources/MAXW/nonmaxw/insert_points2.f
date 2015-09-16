
      real*8 function insert_points2(Nin, vecin, value, SIDE, range, Npts, delta, vecout)
      
	implicit none

      integer, intent(in) :: Nin, Npts, range
      real*8, intent(in)  :: vecin(Nin), value, delta
      real*8, intent(out) :: vecout(Nin+Npts)
      character(1), intent(in) :: SIDE

	integer :: j, indx, find_nearest, ilow, ihig
	real*8 :: auxfine(Npts), delaux, wid, ddd
	
	integer isrchfge
	external isrchfge
	external find_nearest

c     Insert equally spaced Npts points in an original vector near point given by 'value'

c     Nin: size of original vector
c     vecin: origial vector
c     value: value near which the points will be added
c     SIDE: insert points to the left ('L') or to the right ('R') of value
c     range: range of points to be substituted by the Npts new points
c     Npts: number of points to be inserted
c     delta: distance from the (sing.) point 'value', where to begin to insert the new points
c     vecout: output vector of size Nin+Npts with inserted points


c	BEGIN:

	vecout(1:Nin) = vecin(1:Nin)
	ddd = delta

c	Find first point >= value	
	ihig = isrchfge(Nin, vecin(1:Nin), 1, value)
	ilow = ihig-1 
	
	if(SIDE .eq. 'L')then   ! Insert points to the LEFT side of 'value'

		indx = ilow
	
	    if(indx .lt. 3)goto 4444
		
		! Create fine grid
		wid = dabs(value-vecin(indx-range+1))
		delaux = wid / dfloat(Npts+1)
	    if(delta.eq.0.d0 .or. delta>delaux)ddd=delaux

		do j = 1,Npts
			auxfine(Npts-j+1) = (value-ddd) - dfloat(j-1) * delaux
		end do
		
		! Shift original points that are inside fine region
		vecout(indx-range+1:indx) = value - wid

	    ! Shift values Npts to the right
		vecout(indx+Npts:Nin+Npts) = vecout(indx:Nin)

	    ! Fill Npts empty points with fine grid
		vecout(indx+1:indx+Npts) = auxfine(1:Npts)

4444	continue

	else	! Insert points to the RIGHT side of 'value'

		indx = ihig

	    if(indx .ge. Nin)goto 5555

		! Create fine grid
		wid = dabs(vecin(indx+range-1)-value)
		delaux = wid / dfloat(Npts+1)
	    if(delta.eq.0.d0 .or. delta>delaux)ddd=delaux

		do j = 1,Npts
			auxfine(j) = (value+ddd) + dfloat(j-1) * delaux
		end do

		! Shift original points that are inside fine region
		vecout(indx:indx+range-1) = value + wid

		! Shift values Npts to the right
		vecout(indx+Npts:Nin+Npts) = vecout(indx:Nin)

	    ! Fill Npts empty points with fine grid
		vecout(indx:indx+Npts-1) = auxfine(1:Npts)

5555	continue

	end if ! left or right
      
	insert_points2=0

      end