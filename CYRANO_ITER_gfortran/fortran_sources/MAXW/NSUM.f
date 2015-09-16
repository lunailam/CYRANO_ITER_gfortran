      integer function nsum(n,ia)                                       

      implicit none

      integer n, ia(n), i

      nsum=0         
	                                                   
      if(n.gt.0)then                                                  
	  do i = 1, n                                                        
	  nsum = nsum + ia(i)
	  end do
	end if

      return                                                            
      end                                                               
