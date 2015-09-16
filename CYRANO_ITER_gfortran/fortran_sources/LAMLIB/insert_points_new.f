
      real*8 function insert_points_new(Nin, vecin, ind1, ind2, Npts, vecout)
      
	implicit none

      integer, intent(in) :: Nin, Npts, ind1, ind2
      real*8, intent(in)  :: vecin(Nin)
      real*8, intent(out) :: vecout(Nin+Npts)
      
	integer :: j, indx, range
	real*8 :: auxfine(Npts), delaux, wid, ddd
	


c     Simplified version of insert_points2.f:
c     Insert equally spaced Npts points between indexes 
c     ind1 and ind2 of original vector

cccccccccc Version that allows range==1  !!!!!!!

c     Nin: size of original vector
c     vecin: origial vector
c     ind1: first index of interval to be refined
c     ind2: second index of interval to be refined
c     Npts: number of points to be inserted
c     vecout: output vector of size Nin+Npts with inserted points


c	BEGIN:

	vecout(1:Nin) = vecin(1:Nin)
      range = ind2 - ind1  

		! Create fine grid
		wid = dabs(vecin(ind2)-vecin(ind1))
		delaux = wid / dfloat(Npts+1)

		do j = 1,Npts
			auxfine(j) = vecin(ind1) + dfloat(j) * delaux
		end do
		
		! Shift original points that are inside fine region
		if(range>1)vecout(ind1+1:ind1+range-1) = vecin(ind1)

		! Shift values Npts to the right
		vecout(ind2+Npts:Nin+Npts) = vecout(ind2:Nin)

	    ! Fill Npts empty points with fine grid
		vecout(ind2:ind2+Npts-1) = auxfine(1:Npts)

	insert_points_new=0

      end 