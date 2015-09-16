      subroutine groots

      implicit none

C     Generates roots of the local dispersion at element nodes,
C     in equatorial plane, except at magnetic axis.
C     Generates equilibrium profiles in equatorial plane (REMOVED by ERN).
C     Generates AKPERP, kperp of fast wave on magn. axis.

c     19/7/98: changed storage in table ROOTS

      include 'pardim.copy'
      include 'comfic.copy'
      include 'comin2.copy'
      include 'comgeo.copy'
      include 'commag.copy'
      include 'comma2.copy'
      include 'comreg.copy'
      include 'comsub.copy'
      include 'comfin.copy'
      include 'comswe.copy'
      include 'comnud.copy'
      include 'compla.copy'
      include 'comber.copy'
      include 'commod.copy'
      include 'comant.copy'
      include 'complo.copy'
      include 'complp.copy'
      include 'comphy.copy'

      logical inscr
      
      integer j, npldis, mst1m, mst2m, mstsm, nrootp, iii, nkup, iplo, iphi, niphi
     ;, ipo, ipro, incplo, ie1, ie2, ipath, idro, idro2, idro3, iroot
     ;, nrom
     ;, izamax

	character(100) :: GGs, GGs2 
	character(2) :: charaux, charaux2
     
      real second
      
      double precision omc, rlar
      
      complex*16 krho2(3,maxpom)
      
      external inscr, izamax, second
      
c     Max. number of roots at each point, for each mode:
      data nrom/3/ 
      write(nofile,*)'enter dispersion study ; time=',second()-timin
      call dset(9*maxpom*2*(maxnel+maxreg), 0.d0, roots, 1)

cERN	Number of poloidal modes considered for dispersion
	nmoddisp = (mstud2-mstud1) / mstust + 1

c	Number of radial points considered 
      npldis = nreg + nele
 
C     Fast wave kperp on inner boundary, outer equatorial plane:

cERN  1) First element (mag. axis)

      ireg = 1
      iel = 1
      ig = 1
      intab = 1
      y = rx0(0)
      rho = rx0m(0)
      phi = 0.
      intabp = 1
c      r0or = r0orta(intab,1)
c      si = eqt(intab,1,15)
c      co = eqt(intab,1,14)
c      bmodul = bmotab(intab,1)
 
      mst1m = mstud1
      mst2m = mstud2
      mstsm = mstust
      mstud1 = 0
      mstud2 = 0
      mstust = 1
 
      if(vacuum(ireg))then
         nrootp = 1
      else if(coldpl(ireg) .or. .not.flrops(ireg))then
         nrootp = 2
      else
         nrootp = 3
      end if

      write(603,*)'ompe=',omptab(1,1),' vte=',vttab(1,1)
      call disper(1, krho2, nrootp, .false.)

      if(crown)then
      write(nofile,*)'Kperp on crown inner boundary (lfs): (m**-1)'
      else
      write(nofile,*)'Kperp on magnetic axis: (m**-1)'
      end if
      write(nofile,*) (cdsqrt(krho2(iii,1)),iii=1,nrootp)
      mstud1 = mst1m
      mstud2 = mst2m
      mstust = mstsm
       
      nkup = 0
        if(polsym)then
        iplo = 0
        niphi = 1
        else
        iplo = npldis
        niphi = 2
        end if
      write(81,*)npldis
      write(81,*)12

cERN  2) Other elements
      
	do 1 iphi = 1, niphi ! poloidal loop  ++++++++++++++++++++++++++++++
         ipo = 1 + npfft/2 * (iphi - 1)
         phi = polang(ipo)     ! poloidal angle [0 or pi]
         intabp = ipo
         incplo = 3 - 2 * iphi ! radial step: +1 for phi=0  (outer R), 
c											-1 for phi=pi (inner R)

         do 10 ireg = 1, nreg ! radial loop ---------------------------

c	      Number of roots: 	        
	      if(vacuum(ireg))then
               nrootp = 1	 ! only 1 root
            else if(coldpl(ireg) .or. .not.flrops(ireg))then
               nrootp = 2	 ! 2 roots (quadratic)
            else
               nrootp = 3	 ! 3 roots (cubic)
            end if
        
		  ie1 = ifiel(ireg) - 1
		  ie2 = ilael(ireg)

		  do 10 iel = ie1, ie2 
			 
			 intab = iel * (ngauss + 1) + ireg	! radial index
			 rho = abscis(intab)				! minor radius
			 y = abscno(intab)					! norm. minor radius y=rho/a
			 yinv = abscni(intab)				! 1/y = a/rho
			 rhoinv = yinv * rnori				! 1/rho
			 r0or = r0orta(intab,ipo)			! Ro/R
			 si = eqt(intab,ipo,15)				! sin(THETA) - pitch angle
			 co = eqt(intab,ipo,14)				! cos(THETA)
			 bmodul = bmotab(intab,ipo)			! Total magnetic field |Btot|

c			 Toroidal wavenumber: 
			 if(cyl)then
				kphilo = kphi			! kphi = n/Ro (const.)
			 else
				kphilo = kphi * r0or	! kphi = n/R
			 end if

c			 Larmor radius:
			 if(nrootp.eq.3)then		
				omc = qom(2) * bmodul	  ! 2nd species cyclotron freq.
				if(inscr(ireg))then
				   vt(2) = vttab2(intab,2)
				else
				   vt(2) = vttab(intab,2) ! 2nd species therm. velocity
				end if
                  rlar = vt(2) / dabs(omc)  ! 2nd species Larmor radius
			 end if
      
			 iplo = iplo + incplo	! radial increment (incplo = +1 or -1)	
      
c			 To plot vs. R:
			 absdis(iplo) = eqt(intab,ipo,1)  ! major radius R at Z=Zmag
c			 absdis(iplo) = incplo * rho
cERNc			 B0 profile:
cERN			 profil(iplo,2*nspec+1) = bmotab(intab,ipo)

			 ipath = 1
c			 Routine to compute the dispersion roots (krho^2)
			 call disper(ipath, krho2, nrootp, .true.)
			 idro3 = 2 * nrom * ((mstud2 - mstud1) / mstust + 1)
 
			 do m = mstud1, mstud2, mstust
				mr = (m - mstud1) / mstust + 1
				idro = 2 * nrom * (mr - 1)
				do iroot = 1, nrootp
				   roots(iplo,idro+2*iroot-1) = dreal(krho2(iroot,mr))
				   roots(iplo,idro+2*iroot) = dimag(krho2(iroot,mr))
				   mroots(idro+2*iroot-1) = m	
				   mroots(idro+2*iroot) = m	
				end do
			 end do ! (m)

c			 This is not always Bernstein root!!:
c			 Krho max. * Larmor:
			 if(nrootp.eq.3)then
				idro2 = 2 * mr - 1 + idro3
				roots(iplo,idro2) = rlar * dsqrt( cdabs( 
     ;				  krho2(izamax(nrootp, krho2(1,mr), 1), mr)  )) ! real
c				Keta * Larmor:
				roots(iplo,idro2+1) = dabs(kper * rlar)  ! imaginary
			 end if

  10	   continue  ! End of radial loop ------------------------
     
	   iplo = npldis + 1  ! correct radial index 'iplo' for inner radius

   1  continue  ! End of poloidal loop ++++++++++++++++++++++++++


cERN	######### Saving disp. to files (previously in OUTDIS.f) ###########



cERN  1) 2D dispersions



cERN	Equilibrium profiles already written in TABLES.F

c     Equilibrium profiles:
cERN      npldis = nreg + nele
cERN      nkup = 0
cERN        if(polsym)then
cERN        iplo = 0
cERN        niphi = 1
cERN        else
cERN        iplo = npldis
cERN        niphi = 2
cERN        end if

cERN      do 11 iphi = 1, niphi
cERN      ipo = 1 + npfft/2 * (iphi - 1)
cERN      intabp = ipo
cERN      incplo = 3 - 2 * iphi

cERN      do 12 ireg = 1, nreg
cERN      ie1 = ifiel(ireg) - 1
cERN      ie2 = ilael(ireg)
cERN        do iel = ie1, ie2
cERN        intab = iel * (ngauss + 1) + ireg
cERN        rho = abscis(intab)
cERN        y = abscno(intab)
cERN        iplo = iplo + incplo
cERNc       Plot vs. ±rho:
cERN        abspro(iplo) = incplo * rho
           
cERN          if(inscr(ireg))then
cERN            do ipro = 1, nspec
cERN            profil(iplo,ipro) = dfac * dentab2(intab,ipro)
cERN            profil(iplo,ipro+nspec) = temtab2(intab,ipro)
cERN            end do
cERN          else
cERN            do ipro = 1, nspec
cERN            profil(iplo,ipro) = dfac * dentab(intab,ipro)
cERN            profil(iplo,ipro+nspec) = temtab(intab,ipro)
cERN            end do
cERN          end if
cERNc       not appropriate: B should be vs R, not rho:
cERNc       profil(iplo,2*nspec+1) = bmotab(intab,ipo)
cERN        profil(iplo,2*nspec+2) = qtab(intab)
cERN        end do
cERN  12  continue
cERN      iplo = npldis + 1
cERN  11  continue
 
      write(nofile,*)'Dispersion study ended ; time=',second()-timin
      return
      end
