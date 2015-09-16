      double precision function totcun(x)                               

      implicit none
      double precision x
c
c     Normalized (/rho**2) total toroidal current inside rho=x.
c     NB: approx. in noncircular; 
c     current density being then divided by kappa.
c
      include 'pardim.copy'
      include 'comreg.copy'
      include 'comgeo.copy'
      include 'commag.copy'
                                                                        
      ap = rhocur  
      api = 1.d0 / ap
                                                           
      if(x.ge.ap)then
      totcun = ipl / (x*x)                                                     
      else if(x .gt. 1.d-4*ap)then
      totcun = ipl * (1.d0-(1.d0-(x*api)**2)**(alpha+1.d0)) / (x*x)
      else       
      totcun = ipl * (alpha+1.d0) * api*api
      end if   
        
      return                                                            
      end                                                               
