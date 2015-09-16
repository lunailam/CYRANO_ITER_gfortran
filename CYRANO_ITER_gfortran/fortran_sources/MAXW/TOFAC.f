      complex*16 function tofac(ia)                                        

      implicit none
      integer ia                                                        

      include 'pardim.copy'
      include 'comgeo.copy'
      include 'commod.copy'
      include 'comant.copy'
      include 'comphy.copy'
      
      real*8 xtra, dsinc
      
      complex*16 cdeid
      
      external cdeid, dsinc
                                                                        
      if(cyl)then                                                       
      tofac = cdeid(- kphi * zaa(ia))                                
      xtra = kphi * dza(ia) * 0.5d0                                           

      else                                                              
      tofac = cdeid(- n * phiaa(ia))                                 
      xtra = n * dphia(ia) * 0.5d0                                            

      end if                                                            

      tofac = tofac * dsinc(xtra)        

cPL26/8/04
      if(cokpco)tofac = tofac * cdeid(- n * phibarac(ia))
      	                
      return                                                            
      end                                                               
