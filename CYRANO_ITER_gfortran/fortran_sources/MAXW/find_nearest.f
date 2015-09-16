      integer function find_nearest(n, x, targ)
      
	implicit none

      integer n, i
      double precision x(n), targ
      
c     search for nearest element of ordered array x ~ targ
c     return its index

c	Cases where targ is outside table 
      if(targ.le.x(1))then
         find_nearest = 1
         return
      end if

      if(targ.ge.x(n))then
         find_nearest = n
         return
      end if

     
	do i = 2, n

         if(targ .eq. x(i))then
	      find_nearest = i
	      return
         end if

	   if(x(i).gt.targ)then
	   ! 'i' is the index of the first x element
	   ! larger then targ: Need to check if it is the nearest
	     
		 if( dabs((x(i)-targ)) .le. dabs(targ-x(i-1)) )then
	          find_nearest = i	 ! x(i) is the nearest point
		    return
		 else
		    find_nearest = i-1 ! x(i-i) is closer to xtarg
		    return
		 end if
	   	
	   end if

	end do

c      i = 1
c      j = 1
c   1  continue
c      if(x(i).ge.targ)then
c      isrchfge = j
c      return
c      else if(j.lt.n)then
c      j = j + 1
c      i = i + 1
c      go to 1
c      end if


      
      end