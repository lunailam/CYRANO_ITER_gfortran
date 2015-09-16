      subroutine intnod
      
      implicit none
c
c     Used for gigdr=.F.
c     Must follow a call to NORMAT
c     Called inside x and v region loops when the meshes are piecewise uniform.
c      

      include 'pardim.copy'
      include 'comgdr.copy'
      include 'comin2.copy' 

      integer ib, jb, ipov, ipox
      
      do ipov = 0, 1
      do ipox = 0, 1
c     Loop over v basis functions
      do ib = 1, 4
c       Loop over x basis functions
        do jb = 1, 4
        prods(ib,0,ipov,jb,0,ipox) = norma1(ib,jb) * 
     ;                               bafint(ib,0,1,ipov) * bafint(jb,0,1,ipox)
        prods(ib,1,ipov,jb,0,ipox) = norma2(ib,jb) * 
     ;                               bafint(ib,1,1,ipov) * bafint(jb,0,1,ipox)
        prods(ib,0,ipov,jb,1,ipox) = norma3(ib,jb) * 
     ;                               bafint(ib,0,1,ipov) * bafint(jb,1,1,ipox)
        prods(ib,1,ipov,jb,1,ipox) = norma4(ib,jb) * 
     ;                               bafint(ib,1,1,ipov) * bafint(jb,1,1,ipox)
        end do
      end do
      end do
      end do
      
            
      return
      end
