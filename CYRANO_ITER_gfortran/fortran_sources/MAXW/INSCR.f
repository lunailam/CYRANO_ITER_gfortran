      logical function inscr(ireg)

      implicit none
      integer ireg

C     Determines whether current point in region IREG is inside a screen
C     and whether this point deserves special treatment in
C     density and temperature calculation.
c     PHI through COMGEO!

      include 'pardim.copy'
      include 'comgeo.copy'
      include 'comant.copy'
      
      integer js
 
      if(nscree.gt.0)then
        do js = 1, nscree
        inscr = inscr .or.
     ;	     (ireg.gt.irbscr(js) .and.
     ;	      vscrsh(js) .and. (phi-thes1(js))*(phi-thes2(js)).le.0.d0)
        end do
      else
      inscr = .false.
      end if
 
      return
      end
