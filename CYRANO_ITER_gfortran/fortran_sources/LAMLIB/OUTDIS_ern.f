      SUBROUTINE OUTDIS_ern

      IMPLICIT NONE

C     BUILDS DISPERSION OUTPUT TABLES, to be plotted elsewhere.
cERN	ALL equilibrium profiles OUTPUT was transferred to OUTTAB.f

       
      include 'pardim.copy'
      include 'dynou2.copy'
      include 'comgeo.copy'
      include 'comant.copy'
      include 'commod.copy'
      include 'comreg.copy'
      include 'comswe.copy'
      include 'complo.copy'
      include 'complp.copy'
      include 'compla.copy'
      include 'comfin.copy'
      include 'comfic.copy' 

	character(100) :: GGs, GGs2 
	character(2) :: charaux, charaux2
	 
      integer :: i, j, l, ipol, npldis, 
     ;           nrom, nmod, ncols, nppp, nradp

      parameter(nrom=3)

	real*8 :: X1(4,npfft+1)

      DOUBLE PRECISION 
     ;  SECOND
      EXTERNAL SECOND
	
c	----------------------- MAIN ------------------------------

c	NB : The (R,Z) grid for the 2D OUTPUT was already built and
c	     written to Rcyr.dat and Zcyr.dat in subroutine OUTGRID.f.
c		 All relevant variables are passed from OUTGRID through 
c		 COMMONS (nploth, polplo, istp(ireg), ist11, etc... )
 
      WRITE(NOFILE,*)'ENTER OUTDIS ; time=',SECOND()-TIMIN

c     1) Compute the cut-offs and resonances in central region (reg.1)
c        Subroutine CHARSU returns: (1)L-n//2, (2)R-n//2 and (3)S-n//2 
c	   using real parts of the tensors (S,R and L) and considering 
c	   n// = nphi (refractive index of m=0 mode).

      if(.not.vacuum(1)) then

	   x1 = 0		! initialize table
	   tab = 0.0d0  ! initialize table
	   ireg = 1		! consider only central region in CHARSU
	   nradp = istp(ireg)  ! Number of radial points

c        1.1) Toroidal mode
         imoto = 1
         if(cyl)then
            kphi = ktoan(imoto)
         else
            n = motoan(imoto)
            kphi = n * r0i
         end if

c	   1.2) Number of poloidal points (CHECK consistent with OUTNPFT)
         npp  = npfft  + 1
         nppp = nploth + 1
         if(updsym)then
            npp  = npfft/2  + 1
            nppp = nploth/2 + 1
         end if

c        1.3) Compute resonance and cut-off's and store in TAB

	   do i = 1, nradp	! radial loop - - - - - - - - - 
		  intab = intaplot(i)  ! should be =intaplot(i) to take into account OUTGAU = F
c           print *, intab, eqt(intab,1,1), vttab(intab,1)
		  do ipol = 1, npp  ! poloidal loop
		     intabp = ipol
               call charsu(x1(1,ipol)) ! Compute resonances and cut-off's
            end do 

		  do l = 1, 3  ! Loop over res/cutoff + + + +       
	         
			 if(OUTNPFT)then  ! original poloidal GRID
		        tab(1:npp,i,l) = x1(l,1:npp)
                  if(updsym)then
                     do j = 1, npp-1
                        tab(npp+j,i,l) = tab(npp-j,i,l)
                     end do
                  end if
		     else	 ! poloidal interpolation
		        call interp1(polang, polang(npp), x1(l,1), 4,
     ;                         npp, polplo, tab(1,i,l), nppp, 'R')
                  if(updsym)then
                     do j = 1, nppp-1
                        tab(nppp+j,i,l) = tab(nppp-j,i,l)
                     end do
                  end if
		     end if ! (OUTNPFT = true)         
		  
		  end do ! end of res/cutoff loop + + + + + +

         end do ! end of radial loop - - - - - - - - - - - - - 

	end if ! ( not vacuum(1) )
 
	 
c	############### WRITING OUTPUT FILES ################# 

c     1) 1D dispersions for all (nmod) poloidal wavenumbers 
c        NB: abscissae (absdis) is R at Z=Zmag (built in groots.f)
c		   Does not depend on OUTGAU

c	Number of radial points
      npldis = nele + nreg
      if(.not.polsym)npldis = 2 * (nele + nreg)
c	Number of poloidal modes considered for dispersion
	nmod = (mstud2-mstud1) / mstust + 1
c	Number of columns (except first one -> R)
	ncols = 2*nrom*nmod

	call INT_TO_STRING(ncols,  charaux)
	     GGs  = "(A16," // charaux // "i16)"
	call INT_TO_STRING(ncols+1,  charaux2)
	     GGs2  = "(" // charaux2 // "g16.7)"

	open (UNIT = 999, FILE = paplda // 'dispersion_1D.dat', 
     ;      STATUS = "REPLACE", ACTION = "WRITE")
         write(999,*), '1D dispersions vs. R at Z=Zmag (first column)'
	   write(999,"(A7, A7, A7)"), 'rad', 'roots', 'modes'  
	   write(999,"(i7, i7, i7)"), npldis, nrom, nmod
	   write(999,GGs), 'pol.mode', mroots(1:ncols)
	   write(999,"(8A16)"), 'majR','Re','Im','Re','Im','Re','Im','...'          
	   do j = 1, npldis
            write(999,GGs2), absdis(j), roots(j,1:ncols)
         end do
	close(999)


c     2) 2D Fast wave cut-offs and perpendicular resonance 
c	   (Only in central region and for m=0 mode).
c	   NB: Write 3 different files for L-n//2, R-n//2 and S-n//2.
c		   Format: radial(lines) vs. poloidal(columns) to use
c				   with Rcyr.dat and Zcyr.dat (OUTGRID.f)

      if(.not.vacuum(1)) then

c	   2.1) L - n//^2
	   open (UNIT = 999, FILE = paplda // 'L-npar2.dat', 
     ;         STATUS = "REPLACE", ACTION = "WRITE")
	      write(999,*), 'L - n//^2 (m = 0)'
            write(999,*), nradp
            write(999,*), nploth+1
            do j = 1, nradp
               write(999,2000)(tab(ipol,j,1), ipol = 1, nploth+1)
            end do
	   close(999)

c	   2.2) R - n//^2
	   open (UNIT = 999, FILE = paplda // 'R-npar2.dat', 
     ;         STATUS = "REPLACE", ACTION = "WRITE")
	      write(999,*), 'R - n//^2 (m = 0)'
            write(999,*), nradp
            write(999,*), nploth+1
            do j = 1, nradp
               write(999,2000)(tab(ipol,j,2), ipol = 1, nploth+1)
            end do
	   close(999)

c	   2.3) S - n//^2
	   open (UNIT = 999, FILE = paplda // 'S-npar2.dat', 
     ;         STATUS = "REPLACE", ACTION = "WRITE")
	      write(999,*), 'S - n//^2 (m = 0)'
            write(999,*), nradp
            write(999,*), nploth+1
            do j = 1, nradp
               write(999,2000)(tab(ipol,j,3), ipol = 1, nploth+1)
            end do
	   close(999)

      end if  ! ( not vacuum(1) )

      return

 2000 format(200(g14.6))
      end
