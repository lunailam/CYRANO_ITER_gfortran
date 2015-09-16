      subroutine filmao(fac, ldfac)

      implicit none
      integer ldfac
      complex*16 fac(ldfac,*)

c     Inscribes '1.' on matrix diagonal for all dof which vanish 
c     or were eliminated at the inner boundary.
c     FAC is matrix of first (innermost) finite element.

c     Implementations of this routine and BOUAS1 are not independent; should be
c     grouped in one; circular concentric case.

      include 'pardim.copy'
      include 'comequ.copy'
      include 'comreg.copy'
      include 'comsub.copy'
      include 'commod.copy'
      include 'comant.copy'
      include 'comfic.copy'
      include 'comfin.copy'
      include 'compla.copy'
      include 'comgeo.copy'
      include 'comphy.copy'

      logical lotr, hot
      
      character*3 eletyp
      
      integer 
     ;  i, j, iel, idllo, icolo, ibulo, idcatl, m1, m2, nie, ien(6)
     ;, ideca, idecam, idec

      write(6,*)'enter filmao'

c     nb: IEN = 'liste des absents', NIE = 'nombre d'absents'

      iel = 1
      eletyp = styp(1)
      idllo = iddl(1)
      icolo = iconn(1)
      ibulo = ibub(1)
      hot = .not.vacuum(1) .and. .not.coldpl(1) .and. flrops(1)
c      idcatl = ishibc
      idcatl = 0

      m1 = minf(0)
      m2 = msup(0)

      if(crown)then
c     *************
c     Inner metal wall: same conditions for all modes.
      nie = 2
      if(ndof.lt.5)nie = 1
        if(eletyp.eq.'HEC')then
c       E- has been eliminated (E- = E+); E//=0:
        ien(1) = 3
        ien(2) = 5
        else if(eletyp.eq.'M23')then
c       Etheta=0, Ephi=0:
        ien(1) = 2
        ien(2) = 4
c       else if(eletyp.eq.'CAR')then
c       CAR type: metallic condition expressed as convolution
        end if
      end if
c     ******

      do 42 m = m1, m2
c     ----------------
      mr = m + 1 - minf(0)
      ideca = idcatl + (mr-1) * icolo
      idecam = idcatl + nmode(0) * icolo + (mr-1) * ibulo

      if(.not.crown)then
c     ------------------
c     Magnetic axis is a boundary: conditions are mode-dependent.

          if(eletyp.eq.'HEC')then
c         +++++++++++++++++++++++
	
      		if(hot)then

      			if(m.eq.0)then
      			nie=3
      			if(ndof.lt.6)nie=2
c               E+=0, E-=0, E//'=0:
      			ien(1)=1
      			ien(2)=3
      			ien(3)=6
      			else if(m.eq.-1)then
      			nie=5
      			if(ndof.lt.5)nie=4
c               All gone except  E//':
      			ien(1)=1
      			ien(2)=2
      			ien(3)=3
      			ien(4)=4
      			ien(5)=5
      			else if(m.lt.-1)then
      				if( (m/2)*2.ne.m )then
c                   Odd modes, ||>1: all components and first derivative vanish:
      				nie=ndof
      				ien(1)=1
      				ien(2)=2
      				ien(3)=3
      				ien(4)=4
      				ien(5)=5
      				ien(6)=6
      				else
      				nie=5
      				if(ndof.lt.5)nie=3
c                   All gone except  E+':
      				ien(1)=1
      				ien(2)=3
      				ien(3)=4
      				ien(4)=5
      				ien(5)=6
      				end if
      			else if(m.eq.1)then
      			nie=5
      			if(ndof.lt.5)nie=4
c               All gone except  E//':
      			ien(1)=1
      			ien(2)=2
      			ien(3)=3
      			ien(4)=4
      			ien(5)=5
      			else if(m.gt.1)then
      				if( (m/2)*2.ne.m )then
c                   Odd modes, ||>1: all components and first derivative vanish:
      				nie=ndof
      				ien(1)=1
      				ien(2)=2
      				ien(3)=3
      				ien(4)=4
      				ien(5)=5
      				ien(6)=6
      				else
      				nie=5
      				if(ndof.lt.5)nie=3
c                   All gone except  E-':
      				ien(1)=1
      				ien(2)=2
      				ien(3)=3
      				ien(4)=5
      				ien(5)=6
      				end if
      			end if
		
      		else
c		    HEC, not HOT
	
      			if(m.eq.0)then
      			nie = 3
      			if(ndof.lt.6)nie = 2
      			ien(1) = 1
      			ien(2) = 3
      			ien(3) = 6
      			else if(m.eq.-1)then
      			nie = 4
      			if(ndof.lt.5)nie = 3
      			ien(1) = 2
      			ien(2) = 3
      			ien(3) = 4
      			ien(4) = 5
      			else if(m.lt.-1)then
c      			check!
            	if( (m/2)*2.ne.m )then
      				nie = ndof
      				ien(1) = 1
      				ien(2) = 2
      				ien(3) = 3
      				ien(4) = 4
      				ien(5) = 5
      				ien(6) = 6
      				else
      				nie = 5
      				if(ndof.lt.5)nie = 3
      				ien(1) = 1
      				ien(2) = 3
      				ien(3) = 4
      				ien(4) = 5
      				ien(5) = 6
      				end if
      			else if(m.eq.1)then
      			nie = 4
      			if(ndof.lt.5)nie = 3
      			ien(1) = 1
      			ien(2) = 2
      			ien(3) = 4
      			ien(4) = 5
      			else if(m.gt.1)then
c     			check!
      				if( (m/2)*2.ne.m )then
      				nie = ndof
      				ien(1) = 1
      				ien(2) = 2
      				ien(3) = 3
      				ien(4) = 4
      				ien(5) = 5
      				ien(6) = 6
      				else
      				nie = 5
      				if(ndof.lt.5)nie = 3
      				ien(1) = 1
      				ien(2) = 2
      				ien(3) = 3
      				ien(4) = 5
      				ien(5) = 6
      				end if
      			end if
		
      		end if
	
          else if(eletyp.eq.'M23')then
c         ++++++++++++++++++++++++++++
      		if(m.eq.0)then
      		nie = 3
      		if(ndof.lt.6)nie = 2
      		ien(1) = 1
      		ien(2) = 2
      		ien(3) = 5
      		else if(iabs(m).eq.1)then
c		    Bubble d.o.f. at axis eliminated
      		nie = 4
      		if(ndof.lt.5)nie = 3
      		ien(1) = 2
      		ien(2) = 3
      		ien(3) = 4
      		ien(nie) = icolo+1
      		else if(iabs(m).gt.1)then
      			if( (m/2)*2.ne.m )then
c			     bubble d.o.f. at axis eliminated
C check!
      			nie = icolo+1
      			ien(1) = 1
      			ien(2) = 2
      			ien(3) = 3
      			ien(4) = 4
      			ien(5) = 5
      			ien(nie) = icolo + 1
      			else
      			nie = 4
      			if(ndof.lt.5)nie = 2
      			ien(1) = 1
      			ien(2) = 2
      			ien(3) = 4
      			ien(4) = 5
      			end if
      		end if
	
          else if(eletyp.eq.'CAR')then
c         ++++++++++++++++++++++++++++
          nie = 0
	    idec = 0

      		if(m.ne.0)then
      		  if(ndof.lt.6)then
      		  ien(1) = 1
      		  ien(2) = 3
			  nie = 2
	              idec = 2
	              else
      		  nie = 3
	              idec = 3
      		  ien(1) = 1
      		  ien(2) = 3
      	  	  ien(3) = 5
	              end if
	            end if
      		if((m/2)*2 .eq. m)then
c		      Derivatives of even modes vanish
      		  if(ndof.lt.5)then
			  nie = nie + 2
      		  ien(idec+1) = 2
      		  ien(idec+2) = 4
	              else
      		  nie = nie + 3
      		  ien(idec+1) = 2
      		  ien(idec+2) = 4
      		  ien(idec+3) = 6
	              end if
      		end if	
          end if
c         ++++++
      end if
c     ------

      lotr = .not.crown .and. eletyp.eq.'M23' .and. (m/2)*2.ne.m

      do 43 j = 1, nie
      i = ideca + ien(j)
      if(lotr .and. j.eq.nie)i = idecam + 1
      fac(i,i) = cun
  43  continue

  42  continue
c     --------

      write(6,*)'Exit filmao'
      return
      end
