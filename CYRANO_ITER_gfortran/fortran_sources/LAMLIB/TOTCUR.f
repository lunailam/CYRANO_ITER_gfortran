      double precision function totcur(x)                               
c
c     Total toroidal current inside magnetic surface rho = x
      implicit none
      double precision x

      include 'pardim.copy'
      include 'comreg.copy'
      include 'comgeo.copy'
      include 'commag.copy'
                                                                        
      if(x.ge.ap)then
      totcur = ipl                                                        
      else 
      totcur = ipl * (1.d0-(1.d0-(x*api)**2)**(alpha+1.d0))         
c NB: approx: remains the same for analytical D-shape; 
c     current density being then divided by kappa.
      end if   
        
      return                                                            
      end                                                               
