      SUBROUTINE correct_densities(inow)

      IMPLICIT NONE

c	NEW (Mar 2007):
c	Allow for species to have non parabolic densities
c	Species indicated by flag i_special (.ne. 0) in Jet_data_tight
c	Correct majority species density to maintain quasi-neutrality
c	(elec. density fixed!)																		


      include 'pardim.copy'
      include 'compla.copy'
      include 'coequi.copy'
      include 'comgeo.copy'	 
      include 'comgdr.copy'  ! logical 'DIRK'
c      include 'comgeo.copy'  ! plasma radius 	 rx0m(1)


	integer, intent(in) :: inow
	integer :: i_maj, j, NP, OpenStat
	real*8  :: int_elec, int_sp, int_maj, dx, xext(1000),yext(1000), auxdat(1000), dummyr
	character*100 :: FNAME
	character charaux

c	Majority species index ---------------------------------------------
	i_maj = maxloc(spefra(2:nspec,1),DIM=1) + 1


c	Incorporate original species(inow) density in majority --------
	dentab(1:i_rp,i_maj) = dentab(1:i_rp,i_maj) + 
     ;                        zch(inow)/zch(i_maj)*dentab(1:i_rp,inow) 

ccc	spefra(i_maj,1)      = spefra(i_maj,1)      + spefra(inow,1)



c	Read data from file	------------------------------------------------

	if((inow.eq.ispgdr(1) .or. inow.eq.my_ispgdr(1)) .and. DIRK)then
	   ! Read file from QLFP_data/.../ folder
	    charaux = char(48+inow)
	    FNAME = trim(QLFPDIR) // "fort.333ALL"
            print*,'Species',inow 
            print*,'Reading density from',trim(FNAME)
	else
	   ! Read file from INPUT/ folder
	    charaux = char(48+inow)
		FNAME = "./INPUT/dens_spec" // charaux // ".dat"
            print*,'Species',inow
            print*,'Reading density from',trim(FNAME)
	endif


	open(UNIT = 99, FILE = trim(FNAME), STATUS = "OLD",
     ;     IOSTAT = OpenStat, ACTION = "READ")
		 if(OpenStat>0) then
			 print*, trim(FNAME), 'not found!'
			 stop ! goto 666
		 else
c			read ( 99, *, END = 666), NP
		     do j = 1,1000
		        read ( 99, *, END = 666), xext(j), yext(j)
		     end do
		 end if	


666	continue
	close (99)

	NP = j - 1 

c		Interpolate to CYRANO abscissa (xext=rho(m))
c	    Assure data extend to r_p
		if(xext(NP).lt. ap)then
		   xext(NP+1) = ap
		   yext(NP+1) = 0.d0
		   NP = NP+1
	    end if
		call interp2(xext, 1, CONCFAC*yext, 1, NP, 
     ;                 abscis(1:i_rp), auxdat(1:i_rp), i_rp )




c	New density for 'special' species (10^20 in CYRANO def.)
	dentab(1:i_rp,inow) = 1.d-20 * auxdat(1:i_rp)
c	dentab(1:i_rp,inow) = 1.0*0.1*(1.d0-abscno(1:i_rp)**2)**7

c	Correct majority density (ne fixed)
	dentab(1:i_rp,i_maj) = dentab(1:i_rp,i_maj) - 
     ;                       zch(inow)/zch(i_maj)*dentab(1:i_rp,inow) 


c	New species fraction	 
	int_elec=0;
	int_sp=0;
	int_maj=0;

	do j=2,i_rp

		dx = abscis(j) - abscis(j-1)
		int_elec = int_elec + 
     ;               0.5 * (dentab(j,1)*abscis(j) + dentab(j-1,1)*abscis(j-1)) * dx
		int_sp = int_sp + 
     ;               0.5 * (dentab(j,inow)*abscis(j) + dentab(j-1,inow)*abscis(j-1)) * dx
		int_maj= int_maj + 
     ;               0.5 * (dentab(j,i_maj)*abscis(j) + dentab(j-1,i_maj)*abscis(j-1)) * dx
	end do

c	New 'special' species fraction 
	spefra(inow,1) = int_sp/int_elec !/ zch(inow)
c	Correct majority fraction
c	spefra(i_maj,1)      = spefra(i_maj,1)      - spefra(inow,1)
	spefra(i_maj,1) = int_maj/int_elec !/ zch(i_maj)

      return
 
      END SUBROUTINE correct_densities
