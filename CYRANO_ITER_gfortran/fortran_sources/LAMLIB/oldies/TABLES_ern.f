      subroutine tables

	use vec_functions

      implicit none
C
C     Generates tables of equilibrium profiles
C
C-------------------------------------------------------------------------------
C   N.B.: choosing the type of geometry:
C
C   CYL=.T.: 
C     cylindrical (arbitrary section possible) i.e. no toroidal curvature;
C              
C   POLSYM=.T.: 
C     geometry is poloidally symmetrical; active only if CYL=.T., DSHAPE=.F. and
C     GENEQ=.F.
C
C   DSHAPE=.T.:
C     analytical description with elongation, Shafranov shift, triangularity
C     if set to .F., concentric circular magnetic surfaces;
C
C   GENEQ=.T.:
C     equilibrium defined by radial tables of poloidal Fourier coefficients of 
C     cartesian coordinates. Overrules DSHAPE. To be implemented.
C                
C-------------------------------------------------------------------------------
C
C   Tabulated coefficients:
C   -----------------------
C
C    2D tables ti: in array EQT; parity in theta is given for up-down symmetric
C    ------------  cases.
C   
C    #   Parity     Contains
C    1     1        R            
C    2    -1        Y            
C    3     1        dR/dr            
C    4    -1        1/r dR/dtheta            
C    5    -1        dY/dr            
C    6     1        1/r dY/dtheta            
C    7     1        Nn = Sqrt(t4**2+t6**2)            
C    8    -1        Gn = t3*t4 + t5*t6            
C    9     1        Jn = t3*t6 - t5*t4 
C   10     1        mu = (Ytt*Rt-Rtt*Yt) / Nt**2         
C   11    -1        Cn = C/r**2 = (Rt*Yrt-Yt*Rrt) / r**2            
C   12     1        sin(THETA)/r            
C   13     1        tan THETA           
C   14     1        cos THETA
C   15     1        sin THETA
C   16     1        d(sin THETA)/dr
C   17    -1        1/r d(sin THETA)/dtheta           
C
C    1D (radial) tables si: array EQTA1D
C    ----------------------
C   
C    #    Contains
C    1    fn            
C    2    g            
C    3    In = I(r) / r**2            
C    4    q            
C    5    LAMBDAn            
C    6    1/r dPSI/dr
C    7    PSI (values at element nodes only)
C    8    V'/((2pi)**2 rho RA) (1-(rho/ap)**2)**alpha
C    9    dV/dr = (2pi)^2 <RJ> (torus) or 2pi <J> (cylinder)
C   10	Surface perimeter : L(r) = integ [r.Nn(r,theta) dtheta] !ERN
C
C    Add. tables for constant k// coords: 
C    2D in array ckt; parity is given for up-down symmetric cases.
C    ------------------------------------------------------------
C   
C    #   Parity     Contains
C    1     1        dThetabar/dTheta
C    2     1        dPhibar/dTheta
C    3    -1        Thbar (0->2pi)
C    4    -1        Phibar minus phi (+/-)
C    5    -1        Khi (magnetic field strength based poloidal angle) (0->2pi)
C    6    -1        dThbar/drho  
C    7    -1        dPhibar/drho           
C
cPL26/8/04 To do for cokpco noncircular: make 2D tables of R, Y on uniform thetabar mesh
c          will allow easy 2D field plots
C    1D
C    --
C    hachi
C-------------------------------------------------------------------------------

      include 'pardim.copy'
      include 'comfic.copy'
      include 'comgeo.copy'
      include 'commag.copy'
      include 'comma2.copy'
      include 'comreg.copy'
      include 'comsub.copy'
      include 'comfin.copy'
      include 'comswe.copy'
      include 'comin2.copy'
      include 'comrot.copy'
      include 'compla.copy'
      include 'comfou.copy'
      include 'commod.copy'
      include 'comant.copy'
      include 'complo.copy'
      include 'complp.copy'
      include 'comphy.copy'
	include 'coequi.copy'
c
      logical inscr, wr
      
      integer i, im, j, ipo, isp, imi, ima, is1, is2, iel1, iel2
     ;, index
     ;, idamin, idamax, i1st, nsum
      
      double precision polst, t, rle
     ;, gauss, totcun, g
     ;, jnor(npfft+1), jnorav, gfac

	complex*16 :: cdeid

      external inscr, gauss, totcun, g, idamin, idamax, i1st, nsum, gfac, cdeid

cERNcccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	character(120) :: GGs, dummychar
	character(3) :: auxchar
	integer :: Openstat 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


C     Symmetries in 2D tables (used with updsym=.T.):     
      uds(1) =  1
      uds(2) = -1
      uds(3) =  1
      uds(4) = -1
      uds(5) = -1
      uds(6) =  1
      uds(7) =  1
      uds(8) = -1
      uds(9) =  1
      uds(10) =  1
      uds(11) = -1
      uds(12) =  1
      uds(13) =  1
      uds(14) =  1
      uds(15) =  1
      uds(16) =  1
      uds(17) = -1
      udsck(1) =  1
      udsck(2) =  1
      udsck(3) = -1
      udsck(4) = -1
      udsck(5) = -1
      udsck(6) = -1
      udsck(7) = -1
c
      phi = 0.
      i = 0
cERN  ccccccccccccccccccccccccccccccc 
	if (READ_FLUX) updsym = .false.
c     ccccccccccccccccccccccccccccccc	

c     2) Radial coordinate ------------------------------------------- 

cPL2004 revoir points radiaux pour gicurl et/ou giplas false
      
      do ireg = 1, nreg
c     -----------------
	   is1 = nsum(ireg-1, ns) + 1
         is2 = is1 + ns(ireg) - 1
         i = i + 1
         iregoa(i) = ireg
         isuboa(i) = is1
C        List of abscissae (element boundaries + element Gauss points)
C        + region boundaries twice:
         abscno(i) = rx0(ireg-1)
         do isubr = is1, is2
		  iel1 = i1st(isubr) + 1
		  iel2 = i1st(isubr) + iele(isubr)
		  do iel = iel1, iel2
			 ifiabs(iel) = i
			 do j = 1, ngauss
				i = i + 1
				iregoa(i) = ireg
				isuboa(i) = isubr
				abscno(i) = fx0(iel-1) + aga(j) * fl(iel)
			 end do
			 i = i + 1
			 iregoa(i) = ireg
			 isuboa(i) = isubr
			 abscno(i) = fx0(iel)
cERN  		 ccccccccccccccccccccccccccc
     			 if (ireg .eq. 1)  i_rp  = i  ! index for plasma radius
			 if (ireg .eq. 2)  i_ran = i  ! index for antenna radius 
			 if (ireg .eq. 3)  i_rw  = i  ! index for wall radius (=nabsci) 
c			 ccccccccccccccccccccccccccc
		  end do
	   end do
      end do
c     ------

      nabsci = i
      do i = 1, nabsci
	   abscis(i) = abscno(i) * rnorm	! rho(m)
         abscni(i) = 1.d0 / dmax1(abscno(i), 1.d-10)
      end do
	
C     3) Poloidal angles for 2D tables --------------------------------

      if(updsym)then
C        Take up-down symmetry into account (circular and D-shape cases):
         npp = npfft / 2 + 1
      else
         npp = npfft + 1
      end if
      
      polst = twopi / dfloat(npfft)		! poloidal step
      do ipo = 1, npfft
         polang(ipo) = (ipo - 1) * polst  ! poloidal angle
      end do
      polang(npfft+1) = twopi				! overlap with first point


c     4) Profile tables at Gauss points and region boundaries: -------
c
      phi = 0.
      eqta1d(1,6) = 0.
cERN  cccccccccccccccccccccccccccc
      if(READ_FLUX)then
c		Read flux surfaces(R,Z), Btor, pitch angle(THETA) and safety factor 
c		given by FLUSH. The output is written in common variables (COEQUI).
		call read_equi(Nrad_flu)  
	end if
c	cccccccccccccccccccccccccccc

	do 30 i = 1, nabsci		! Radial loop
c     ====================================
         ireg = iregoa(i)
         y = abscno(i)		! rho/ap
         rho = abscis(i)		! rho

cERN  NEW option for kappa profile (force kappa = 1 near magnetic axis) 
	if(dshape .and. kapprof) kappa = kappa_orig - (kappa_orig-1) * exp(-50*y*y)
c	---------------------------------------------------------------------------

C        Out of screens:
	   if(READ_PROFILES)then
            call ypro_read(.false.)	! Read n and T profiles from file
	   else
            call ypro(.false.)
         end if
         do isp = 1, nspec
		  vttab(i,isp)  = vt(isp)		! Thermal velocity
		  omptab(i,isp) = omp(isp)		! Plasma frequency
		  temtab(i,isp) = temper(isp)	! Temperature
		  dentab(i,isp) = densit(isp)	! Density
         end do

         if(nscree.gt.0)then
C        Inside screens:
            call ypro(.true.)
            do isp = 1, nspec
			 vttab2(i,isp)  = vt(isp)
			 omptab2(i,isp) = omp(isp)
			 temtab2(i,isp) = temper(isp)
			 dentab2(i,isp) = densit(isp)
		  end do
         end if

cERN	   Seems to be at wrong place! Why not call together with DSHAPE2?		 
	   if(circ)then
C        Magnetic field angle sine, cosine and radial (y) derivative:
C        (1D table for circular cross-section)
           call comag2
		 sitab(i) = si
		 cotab(i) = co
		 si1tab(i) = si1
		 co1tab(i) = co1
	   end if
c
c        5) Geometrical coefficients: 2D tables
c
         do 40 ipo = 1, npp	 ! Poloidal loop
c	   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            phi = polang(ipo)  ! poloidal angle

cERN        cccccccccccccccccccccc
	      if(READ_FLUX)then
			 call genshap(i) ! Compute geom. coef. for general equlibrium
c		     NB: 'read_equi' must have been called first (outside pol. loop)	
			 goto 888
		  end if
c		  cccccccccccccccccccccc
            if(geneq)then
cERN		  All things already done in previous call (READ_FLUX = true)
c	      !!!Need to decide the link between geneq and READ_FLUX!!!
            else if(dshape)then
                 call dshap1(2) ! Compute geom. coef. for D-shape equlibrium
            else
			   call comag1(2) ! Compute geom. coef. for circular equlibrium
			   cn = 0.d0
			   newmu = 1.d0
			   jacn = 1.d0
			   g12n = 0.d0
			   ntn = 1.d0
		  end if

888	      continue
            r0orta(i,ipo) = r0or	 ! = Raxis/R

c	      Geom. quantities READY to be written to 2D tables 'eqt'.
c	      Pitch angle related quantities will be computed later.
c           N.B.: suffix 'n' below means normalized by radial coordinate,
c                 to have nonsingular coefficients at magnetic axis.

            eqt(i,ipo,1) = r			! R(rho,theta) coordinate
            eqt(i,ipo,2) = z			! Z(rho,theta) coordinate
            eqt(i,ipo,3) = drrho		! dR/drho
            eqt(i,ipo,4) = drthn		! 1/rho dR/dtheta
            eqt(i,ipo,5) = dzrho		! dZ/drho
		  eqt(i,ipo,6) = dzthn		! 1/rho dZ/dtheta
		  eqt(i,ipo,7) = ntn		! sqrt(drthn*drthn + dzthn*dzthn)
		  eqt(i,ipo,8) = g12n		! drrho*drthn + dzrho*dzthn
		  eqt(i,ipo,9) = jacn		! drrho*dzthn - dzrho*drthn
		  eqt(i,ipo,10) = newmu		! (d2zth2*drdth - d2rth2*dzth) / ntn**2 
		  eqt(i,ipo,11) = cn		! d2Z/drhodth*dRdth - d2R/drhodth*dZth

  40     continue	! end of poloidal loop (ipo = 1:npp)
c	   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
cERN	   Surface perimeter = integ [r.ntn(r,theta) dtheta] (0,2pi)	
	   eqta1d(i,10) = rho * sum(eqt(i,1:npp,7)) * 2*pi / dfloat(npp)

c	   Reset 1D tables
         eqta1d(i,4) = 0.
         eqta1d(i,5) = 0.
         eqta1d(i,8) = 0.
         eqta1d(i,9) = 0.
      
         if(cyl)then
c           Omit factor 1/R in coeff. Lambda and in q:
            do ipo = 1, npp-1
		     eqta1d(i,4) = eqta1d(i,4) + eqt(i,ipo,9)
		     eqta1d(i,5) = eqta1d(i,5) + eqt(i,ipo,7) ** 2 / eqt(i,ipo,9)
		     eqta1d(i,9) = eqta1d(i,9) + eqt(i,ipo,9)
            end do
	      eqta1d(i,5) = eqta1d(i,5) / dfloat(npp-1)
            eqta1d(i,8) = eqta1d(i,9) / dfloat(npp-1)
            eqta1d(i,9) = eqta1d(i,8) * twopi * rho
         else
            do ipo = 1, npp-1
c              For q:
               eqta1d(i,4) = eqta1d(i,4) + eqt(i,ipo,9) / eqt(i,ipo,1)
c              For LAMBDAn:
               eqta1d(i,5) = eqta1d(i,5) + eqt(i,ipo,7) ** 2
     ;                           / (eqt(i,ipo,1) * eqt(i,ipo,9))
c              For integrand of I(rho) and V'(rho):
               eqta1d(i,9) = eqta1d(i,9) + eqt(i,ipo,1) * eqt(i,ipo,9)
            end do
            eqta1d(i,5) = eqta1d(i,5) / dfloat(npp-1)
	      eqta1d(i,9) = eqta1d(i,9) / dfloat(npp-1)
c           Pol. average of (R/Ra Jn):
		  eqta1d(i,8) = eqta1d(i,9) / ra
c		  V'(rho):
		  eqta1d(i,9) = eqta1d(i,9) * twopi**2 * rho
         end if
         if(ireg.eq.1)then
            eqta1d(i,8) = eqta1d(i,8) * dmax1(1.d0-(rho*api)**2, 0.d0)**alpha
     ;                    * rho
         else
            eqta1d(i,8) = 0.d0
         end if
         eqta1d(i,2) = g(rho)
         if(circ)then
		  qtab(i) = q
         end if
 
  30  continue	! End of radial loop
c     ===============================
c	 5) Internal quantities (1D) 
c		NB: These 1D profiles are necessary to compute some 2D tables in D-shape 
c	    geometry, such as the pitch angle, toroidal magnetic field, Psi, etc...) 

c     5.1) Trapeze integration for I(rho):
c     NB if crown=.true. this fails to get the total current inside rho!
      eqta1d(1,3) = 0.d0
      do i = 2, nabsci
         eqta1d(i,3) = eqta1d(i-1,3) + (abscis(i) - abscis(i-1)) * 0.5d0 *
     ;   (eqta1d(i-1,8) + eqta1d(i,8))
      end do

c     5.2) I(rho) / rho ** 2:
      t = ipl / eqta1d(nabsci,3)
      do i = 2, nabsci
         eqta1d(i,3) = eqta1d(i,3) * t / abscis(i) ** 2
      end do
      eqta1d(1,3) = eqta1d(2,3)

c     5.3) dPsi/(r dr):  ('divided by Ra' in straight cyl.)
      do i = 1, nabsci
         eqta1d(i,6) = mu0 * eqta1d(i,3) / eqta1d(i,5)
      end do

c     5.4) Trapeze integration for psi(rho):
      eqta1d(1,7) = 0.
      do i = 2, nabsci
         eqta1d(i,7) = eqta1d(i-1,7) + (abscis(i) - abscis(i-1)) * 0.5d0 *
     ; (eqta1d(i-1,6)*abscis(i-1) + eqta1d(i,6)*abscis(i))
      end do

c     5.5) Volume of each finite element shell (surface in cylinder case):
      do iel = 1, nele
	   i = ifiabs(iel)
	   eltvol(iel) = gauss(fl(iel), abscis(i), eqta1d(i,9), 1, 0)
	end do


c      6) Magnetic field modulus and pitch angle (2D tables)
c		NB: THETA = atan(Bpol/Btor) is the pitch angle.
        
      do i = 1, nabsci  ! Radial loop	
c	================         
	   fn = 2.d-7 * eqta1d(i,3) / (b0 * eqta1d(i,5))
         if(.not.cyl)fn = fn / ra
         eqta1d(i,1) = fn		 ! f/rho
         eqta1d(i,4) = eqta1d(i,4)*eqta1d(i,2) / (eqta1d(i,1)*dfloat(npp-1))
	   qtab(i) = eqta1d(i,4) ! q-factor
         intab = i  
         do ipo = 1, npp  ! Poloidal loop
c	   ~~~~~~~~~~~~~~~	
		  intabp = ipo
cERNcccccccccccccccccccccc
		  if(geneq)then
c		     call genshap2  ! Not ready yet !!!!!
		  else
		     call dshap2	! Compute |B| and THETA in D-shape geometry
		  end if
cccccccccccccccccccccccccc		  
		  eqt(i,ipo,13) = tgth		  ! tan(THETA)
		  eqt(i,ipo,14) = co		  ! cos(THETA)
		  eqt(i,ipo,15) = si		  ! sin(THETA)
		  eqt(i,ipo,12) = co*tgthor   ! sin(THETA)/rho
		  bmotab(i,ipo) = bmodul	  ! |B| = sqrt(Btor^2 + Bpol^2)
         end do
c	   ~~~~~~
      end do
c	======      

c      7) Use Gaussian integration to build table of poloidal flux function:
c         for extra accuracy at element nodes ONLY.
        
      eqta1d(1,7) = 0.
      i = 0
      do ireg = 1, nreg
	   i = i + 1
         if(i.gt.1)eqta1d(i,7) = eqta1d(i-1,7)
C           Gauss quadrature of EQTA1D(,6)*RHO; 
C		  add up successive results at element nodes.
		  do iel = ifiel(ireg), ilael(ireg)
			 rle = fl(iel) * rnorm
			 im = i
			 i = i + ngauss + 1
			 eqta1d(i,7) = eqta1d(im,7) + 
     ;		 gauss(rle, abscis(im+1), eqta1d(im+1,6), 1, 1)
		  end do
      end do

c     8) Numerical derivatives:
c     8.1) Radial: caring for discontinuities that may occur at region boundaries.
      i = 0
      do ireg = 1, nreg
         i = i + 1
         t = 1.d0 / (abscis(i+1) - abscis(i))
         do ipo = 1, npp
            eqt(i,ipo,16) = (eqt(i+1,ipo,15) - eqt(i,ipo,15)) * t
         end do       
         do ig = 1, (ilael(ireg)-ifiel(ireg)+1) * (ngauss+1) - 1
            i = i + 1
            t = 1.d0 / (abscis(i+1) - abscis(i-1))
            do ipo = 1, npp
			 eqt(i,ipo,16) = (eqt(i+1,ipo,15) - eqt(i-1,ipo,15)) * t
		  end do       
         end do
         i = i + 1
         t = 1.d0 / (abscis(i) - abscis(i-1))
         do ipo = 1, npp
		  eqt(i,ipo,16) = (eqt(i,ipo,15) - eqt(i-1,ipo,15)) * t
         end do       
      end do
c     8.2) Poloidal:
      t = 1.d0 / (2.d0 * polst)
      do i = 1, nabsci
        do ipo = 2, npp-1
		 eqt(i,ipo,17) = (eqt(i,ipo+1,12) - eqt(i,ipo-1,12)) * t
        end do       
      end do

      if(updsym)then
        do i = 1, nabsci
		 eqt(i,1,17) = 0.d0
		 eqt(i,npp,17) = 0.d0
        end do
      else
        do i = 1, nabsci
		 eqt(i,1,17) = (eqt(i,2,12) - eqt(i,npp-1,12)) * t
		 eqt(i,npp,17) = eqt(i,1,17)
        end do
      end if ! (updsym)

c     9) Using symmetries:
      if(updsym)then
        do ipo = 1, npp-1
          do i = 1, nabsci
		   bmotab(i,npp+ipo) = bmotab(i,npp-ipo)
             r0orta(i,npp+ipo) = r0orta(i,npp-ipo)
          end do
        end do
        do index = 1, 17
           if(uds(index).eq.1)then
C     9.1) Even coefficients:
              do ipo = 1, npp-1
                 do i = 1, nabsci
			      eqt(i,npp+ipo,index) = eqt(i,npp-ipo,index)
			   end do
              end do
		 else
C     9.2) Odd coefficients:
              do ipo = 1, npp-1
			   do i = 1, nabsci
				  eqt(i,npp+ipo,index) = - eqt(i,npp-ipo,index)
			   end do
			end do
		 end if
        end do
      end if ! (updsym)

c     10) Find bmin, bmax on each magn. surface:

      if(.not.updsym) then
        do i = 1, nabsci
		 imi = idamin(npp, bmotab(i,1), nabplo)
		 ima = idamax(npp, bmotab(i,1), nabplo)
cERN		 cccccccccccccccccccccccccccccccc
c		 imi = minloc(bmotab(i,1:npp),1)
c		 ima = maxloc(bmotab(i,1:npp),1)
c		 cccccccccccccccccccccccccccccccc
		 bmin(i) = bmotab(i,imi)
		 bmax(i) = bmotab(i,ima)
		 write(nofile,*)imi,ima,bmin(i),bmax(i)
		 bbar(i) = 0.5d0 * (bmax(i)+bmin(i))
		 delb(i) = 0.5d0 * (bmax(i)-bmin(i))
           if(delb(i).lt.1.d-6)then
		 do ipo = 1, npp
		 ckt(i,ipo,5) = polang(ipo)
		 end do
           else
c			2D table of mu, on full section:
			do ipo = 1, npp
			   t = (bbar(i) - bmotab(i,ipo)) / delb(i)
			   ckt(i,ipo,5) = dacos(dsign(dmin1(1.d0,dabs(t)),t))
     ;	     * dfloat(isign(1, (ima-imi)*(ima-ipo)*(ipo-imi)))
			end do
			do ipo = ima+1, npp
c			This gives a range of mu slightly larger than [0,2pi] but which should
c			be suitable for contour plots:
			   ckt(i,ipo,5) = ckt(i,ipo,5) + twopi
c			   if(ckt(i,ipo,5).lt.0.)ckt(i,ipo,5) = ckt(i,ipo,5) + twopi
			end do
c			More accurate? 
CC			Caution with 2 pi jumps!
			if(imi.lt.ima)then
			   ckt(i,imi,5) = 0.d0
			else
			   ckt(i,imi,5) = twopi
			end if
			ckt(i,ima,5) = pi
          end if
        end do      

      else	! Up-down sym. case
        gloext = .false.
        wr = .true.
        do i = 1, nabsci
           mulext(i) = .false.
c          Magn. axis:
           if(i.eq.1)then
              imi = 1
              ima = npp
           else
c             Other absc:
              imi = idamin(npp, bmotab(i,1), nabplo)
              ima = idamax(npp, bmotab(i,1), nabplo)
cERN			cccccccccccccccccccccccccccccccccc
c			imi = minloc(bmotab(i,1:npp),1)
c			ima = maxloc(bmotab(i,1:npp),1)
c			cccccccccccccccccccccccccccccccccc
          end if
          if(imi.ne.1 .or. ima.ne.npp)then
             mulext(i) = .true.
             gloext = .true.
             if(wr)then
                write(6,*)'i=',I,': bmin and bmax at theta=',polang(imi)
     ;          ,' ',polang(ima),'; imi, ima:',imi,ima,' (first occurence)'
                wr = .false.
             end if
          end if
		bmin(i) = bmotab(i,imi)
		bmax(i) = bmotab(i,ima)
		bbar(i) = 0.5d0 * (bmax(i) + bmin(i))
		delb(i) = 0.5d0 * (bmax(i) - bmin(i))
          if(delb(i).lt.1.d-6)then
            do ipo = 1, npp
			 ckt(i,ipo,5) = polang(ipo)
            end do
          else
            do ipo = 1, npp
c			 2D table of mu, on top half section:
			 t = (bbar(i) - bmotab(i,ipo)) / delb(i)
			 ckt(i,ipo,5) = dacos(dsign(dmin1(1.d0,dabs(t)),t))
     ;		 * dfloat(isign(1, (ima-ipo)*(ipo-imi)))
C     ;		 * dfloat(isign(isign(1,ima-imi), (ima-ipo)*(ipo-imi)))
c			 More accurate?
c			 t1 = bmotab(i,npp) - bmotab(i,ipo)
c			 t2 = bmotab(i,ipo) - bmotab(i,1)
c			 ckt(i,ipo,5) = dacos( (t1-t2)/(t1+t2) )
            end do
            ckt(i,1,5) = 0.d0
            ckt(i,npp,5) = pi
          end if
          do ipo = 1, npp - 1
             ckt(i,npp+ipo,5) = twopi - ckt(i,npp-ipo,5)
          end do
        end do
      
	end if !(.not.updsym)

c     11) k// = constant coordinates ----------------------------------

c      Use rf field Fourier expansion in generalized poloidal and toroidal 
C      angles making k// constant on flux surfaces. Additional tables:

	if(cokpco)then

	   do i = 1, nabsci
c     11.1) dThbar/dTheta and 'H' coefficient: 
         hachi(i) = 0.d0
           do ipo = 1, npp
           ckt(i,ipo,1) = eqt(i,ipo,7) / eqt(i,ipo,12)
	     hachi(i) = hachi(i) + ckt(i,ipo,1)
	     end do       
		  hachi(i) = hachi(i) / dfloat(npp)
c		  Makes hachi the inverse of cursive H in paper:
		  if(hachi(i).ne.0.d0)hachi(i) = 1.d0 / hachi(i)
		  jnorav = 0.d0
		  do ipo = 1, npp
			 ckt(i,ipo,1) = ckt(i,ipo,1) * hachi(i)
c			 dPhibar/dTheta: 
			 jnor(ipo) = eqt(i,ipo,9) / eqt(i,ipo,1)
			 jnorav = jnorav + jnor(ipo)
		  end do       
		  jnorav = jnorav / dfloat(npp)
c           if(i.eq.nabsci)then
c           do ipo = 1, npp
c           dphibar/dtheta: first method
c           ckt(i,ipo,2) = ckt(i,ipo,1) * eqta1d(i,4)
c     ;     - eqta1d(i,2)*eqt(i,ipo,9)/(eqta1d(i,1)*eqt(i,ipo,1))
c           write(nofile,*)ckt(i,ipo,2)   
c           end do   
c           end if

c     11.2) dPhibar/dTheta: 
		  do ipo = 1, npp
			 ckt(i,ipo,2) = eqta1d(i,4) * (ckt(i,ipo,1) - jnor(ipo) / jnorav)
		  end do       
		  ckt(i,1,3) = 0.d0
		  ckt(i,1,4) = 0.d0
c     11.3) Thbar and Phibar, using trapeze integration: 
		  do ipo = 2, npp
		  ckt(i,ipo,3) = ckt(i,ipo-1,3) + polst * 0.5 * (ckt(i,ipo-1,1) + ckt(i,ipo,1))
		  ckt(i,ipo,4) = ckt(i,ipo-1,4) + polst * 0.5 * (ckt(i,ipo-1,2) + ckt(i,ipo,2))
		  end do
		  if(updsym)then
			 ckt(i,npp,3) = pi
		  else
			 ckt(i,npp,3) = 2.d0 * pi
		  end if
		  ckt(i,npp,4) = 0.d0
         end do

c     11.4) dThbar/drho and dPhibar/drho: 
         i = 0
         do ireg = 1, nreg
		  i = i + 1
		  t = 1.d0 / (abscis(i+1) - abscis(i))
		  do ipo = 1, npp
			 ckt(i,ipo,6) = (ckt(i+1,ipo,3)-ckt(i,ipo,3)) * t
			 ckt(i,ipo,7) = (ckt(i+1,ipo,4)-ckt(i,ipo,4)) * t
		  end do       
		  do ig = 1, (ilael(ireg)-ifiel(ireg)+1)*(ngauss+1) - 1
			 i = i + 1
			 t = 1.d0 / (abscis(i+1) - abscis(i-1))
			 do ipo = 1, npp
				ckt(i,ipo,6) = (ckt(i+1,ipo,3)-ckt(i-1,ipo,3)) * t
				ckt(i,ipo,7) = (ckt(i+1,ipo,4)-ckt(i-1,ipo,4)) * t
			 end do       
		  end do
		  i = i + 1
		  t = 1.d0 / (abscis(i) - abscis(i-1))
		  do ipo = 1, npp
			 ckt(i,ipo,6) = (ckt(i,ipo,3)-ckt(i-1,ipo,3)) * t
			 ckt(i,ipo,7) = (ckt(i,ipo,4)-ckt(i-1,ipo,4)) * t
		  end do       
	   end do ! (ireg)

cERN 24Jan05
c     11.5) 2D table of exp[-i.n.(Phibar-Phi)] in thetabar mesh 
	   do i = 1, nabsci
		do ipo = 1, npp
	      einphb(i,ipo) = cdeid(- n * ckt(i,ipo,4))
	      end do
         call interp2(ckt(i,1:npp,3), 1, einphb(i,1:npp), 1, npp, polang, einphb(i,1:npp), npp)  		  
	   end do

	print *, dcmplx(dcos(ckt(50,1:10,4)), dsin(ckt(50,1:10,4)))
	print *
	print*, cdeid(ckt(50,1:10,4))
	print *
	print*, cdeid_vec(ckt(50,1:10,4))
	stop

c    	11.6) Up-down symmetry
         if(updsym)then
		  do ipo = 1, npp-1
			 do i = 1, nabsci
				ckt(i,npp+ipo,1) = ckt(i,npp-ipo,1)
				ckt(i,npp+ipo,2) = ckt(i,npp-ipo,2)
				ckt(i,npp+ipo,3) = twopi - ckt(i,npp-ipo,3)
				ckt(i,npp+ipo,4) = - ckt(i,npp-ipo,4)
				ckt(i,npp+ipo,6) = - ckt(i,npp-ipo,6)
				ckt(i,npp+ipo,7) = - ckt(i,npp-ipo,7)
	                  einphb(i,npp+ipo) = dconjg(einphb(i,npp-ipo))
			 end do
		  end do
	   end if
      
c     Generation of Hkl coefficients: see call in GENERA
      
	end if !(cokpco = true)


ccccccccccccccccccccccccccccccccc ERN ccccccccccccccccccccccccccccccccccccc
c     12) Updating COMMONS that will be used later in other CYRANO routines
c	    NB: All other 1D tables (eqta1d) are only used internally in 
c		    TABLES.f or in OUTTAB.f 

	qfactor(1:nabsci)    = eqta1d(1:nabsci,4)	! Safety factor q(r)
c	used in: OUTPOW, FOUCOG2, FOUCOJ, M12COKPCO, OUTDIS (only save)  
	dPsidr_n(1:nabsci)   = eqta1d(1:nabsci,6)	! 1/r dPsi/dr
c	used in: DIRESP, POWGDR, M12COKPCO 
	surf_perim(1:nabsci) = eqta1d(1:nabsci,10)	! Surface perimeter
c	used in: PRESO2 (dynmod = T)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cERN	All 1D and 2D TABLES are written to files in OUTTAB.f
cERN  The (R,Z) grid used for PLOTTING are written in OUTGRID.f

c     13) Write (R,Z) grid used in CALCULATIONS (Gauss pts + el. boundaries)
c	    NB: These files are NOT used in the plot routines (see OUTGRID.f)
c	    FORMAT: (nabsci) lines x (npfft+1) columns

c     13.1) R(rho,theta) coordinate
	open (UNIT = 777, FILE = paplda // 'Rall.dat', 
     ;      STATUS = "REPLACE", ACTION = "WRITE")
		write(777,*), 'R coordinate'
		write(777,*), ' rad ' , ' pol ', ' pla ' 
		write(777,"(i5, i5, i5)"), nabsci, npfft+1, i_rp 
		do j = 1, nabsci	
		   write(777,2222)(eqt(j,ipo,1), ipo = 1, npfft+1)
		end do
	close(777)

c    	13.2) Z(rho,theta) coordinate
	open (UNIT = 777, FILE = paplda // 'Zall.dat', 
     ;      STATUS = "REPLACE", ACTION = "WRITE")
		write(777,*), 'Z coordinate'
		write(777,*), ' rad ' , ' pol ', ' pla ' 
		write(777,"(i5, i5, i5)"), nabsci, npfft+1, i_rp 
		do j = 1, nabsci	
		   write(777,2222)(eqt(j,ipo,2), ipo = 1, npfft+1)
		end do
	close(777)


      return

 2222 format(200(g14.6))   

      end



