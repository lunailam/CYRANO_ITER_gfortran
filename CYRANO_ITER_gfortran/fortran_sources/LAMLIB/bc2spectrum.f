      subroutine bc2spectrum(co, si, co1, si1, k)                                                

      implicit none

      integer k
      complex*16 co, si, co1, si1                                         

c     Computes poloidal harmonic k of matrix 'bc2' of coord.transfo.
c     with 5 rows,6 cols, used to convert (+,+',-,-',//,//') into (rho, theta, theta', phi, phi')                                              
c     Results are sent through comro2       
c     For use in D-shaped and general geometry.
c     Input: a poloidal harmonic index (k), and the associated Fourier 
c     coefficient of cos and sin THETA (resp. co, si)
c     and their radial (y) derivative (resp. co1, si1)
c     NB: these variables are complex!                        
                                                                        
      include 'pardim.copy'
      include 'comgeo.copy'
c      include 'comma2.copy'
      include 'commag.copy'
      include 'comro2.copy'
      include 'comfin.copy'
      include 'comin2.copy'
      include 'comswe.copy'
      include 'comphy.copy'

      integer i, j
      double precision t                                                            
      complex*16 ist, ict, is1t, ic1t                                         
                                                                        
      t = sqrt2i                                                          
                                                                        
      ist = ci * si * t                                              
      ict = ci * co * t                                              
      is1t = ci * si1 * t                                            
      ic1t = ci * co1 * t                                            
                                                                        
      call zset(30, czero, bc2, 1)                                    
        if(k.eq.0)then                                                               
        bc2(1,1) = dcmplx(t, 0.d0)                                                      
        bc2(1,3) = dcmplx(t, 0.d0)     
        end if                                                 
                                                                        
      bc2(2,1) = - ict                                                  
      bc2(2,3) =   ict                                                  
      bc2(2,5) =   si                                         
                                                                        
      bc2(3,1) = -ic1t                                                  
      bc2(3,2) = -ict                                                   
      bc2(3,3) =  ic1t                                                  
      bc2(3,4) =  ict                                                   
      bc2(3,5) = si1                                          
      bc2(3,6) = si                                         
                                                                        
      bc2(4,1) = ist                                                    
      bc2(4,3) = - ist                                                  
      bc2(4,5) = co                                           
                                                                        
      bc2(5,1) = is1t                                                   
      bc2(5,2) = ist                                                    
      bc2(5,3) = - is1t                                                 
      bc2(5,4) = - ist                                                  
      bc2(5,5) = co1                                         
      bc2(5,6) = co                                           
                                                                        
        do i = 1, 5                                                     
          do j = 1, 6                                                     
          bc2h(j,i) = dconjg(bc2(i,j))                                     
          end do
        end do
                                                                        
      return                                                            
      end                                                               
