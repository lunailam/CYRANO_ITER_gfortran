      double precision function dgamma(x)
      
      double precision x, s14aaf
      integer ifail
      external s14aaf
      
c     gamma function, using Nag:    
      ifail = 0  
      dgamma = s14aaf(x, ifail)
      return
      end