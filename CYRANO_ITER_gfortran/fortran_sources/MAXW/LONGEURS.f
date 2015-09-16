      subroutine longeurs(n, x, xl, i)                                         

      implicit none
      integer n, i
      double precision x(0:n), xl(n)                                     

c -construit longeurs intervalles a partir des abscisses-               
c  (i=1),ou inversement(i=2)                                            

      integer j
      
      go to (11,22)i 
                                                         
  11  do 1 j = 1, n                                                        
      xl(j) = x(j) - x(j-1)                                                 
  1   continue                                                          
      return                                                            

   22 do 2 j = 1, n                                                        
   2  x(j) = x(j-1) + xl(j)                                                 

      return                                                            
      end                                                               
c                                                                       
c**********************************************************             
c                                                                       
