      subroutine loctr2                                                 

      implicit none

c     Computes matrix 'bc2' of coord.transfo. at a radial table point.
c     with 5 rows,6 cols.                                               
c     Trigono. data of magnetic field angle provided                    
c     by tables. Index in tables intab input through comswe.                  
c     Results are sent through comro2       
c     Not for use in dshape!                            
                                                                        
      include 'pardim.copy'
      include 'comgeo.copy'
      include 'comma2.copy'
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
                                                                        
      si = sitab(intab)                                                 
      co = cotab(intab)                                                 
      si1 = si1tab(intab)                                               
      co1 = co1tab(intab)                                               
                                                                        
      ist = dcmplx(0.d0,si*t)                                              
      ict = dcmplx(0.d0,co*t)                                              
      is1t = dcmplx(0.d0,si1*t)                                            
      ic1t = dcmplx(0.d0,co1*t)                                            
                                                                        
      call zset( 30, czero, bc2, 1 )                                    
                                                                        
      bc2(1,1) = t                                                      
      bc2(1,3) = t                                                      
                                                                        
      bc2(2,1) = - ict                                                  
      bc2(2,3) =   ict                                                  
      bc2(2,5) =   dcmplx(si,0.d0)                                         
                                                                        
      bc2(3,1) = -ic1t                                                  
      bc2(3,2) = -ict                                                   
      bc2(3,3) =  ic1t                                                  
      bc2(3,4) =  ict                                                   
      bc2(3,5) = dcmplx(si1,0.d0)                                          
      bc2(3,6) = dcmplx(si,0.d0)                                           
                                                                        
      bc2(4,1) = ist                                                    
      bc2(4,3) = - ist                                                  
      bc2(4,5) = dcmplx(co,0.d0)                                           
                                                                        
      bc2(5,1) = is1t                                                   
      bc2(5,2) = ist                                                    
      bc2(5,3) = - is1t                                                 
      bc2(5,4) = - ist                                                  
      bc2(5,5) = dcmplx(co1,0.d0)                                          
      bc2(5,6) = dcmplx(co,0.d0)                                           
                                                                        
      do 1 i = 1, 5                                                     
      do 1 j = 1, 6                                                     
      bc2h(j,i) = dconjg( bc2(i,j) )                                     
    1 continue                                                          
                                                                        
      return                                                            
      end                                                               
