      program cyrano

cJAC	USE DFPORT		! ccccccccccc ERNESTO cccccccccc
cJAC	USE DFLIB		! ccccccccccc ERNESTO cccccccccc

ccc      This causes a lot of trouble (variable name conflicts):
ccc      USE numerical_libraries  ! PL (for IMSL compilation checks)

      implicit none 

c     This is the main program with calls to various routines performing
c     reading of data, solver initialization, system assembly, direct solution
c     using out-of-core memory and numerical output.

C     "Parameter" statements:
      include 'pardim.copy'

C     "Common" statements:
      include 'comusr.f'
      include 'comfic.copy'
c      include 'comerr.copy'
      include 'comgeo.copy'
      include 'commag.copy'
      include 'comma2.copy'
      include 'comreg.copy'
      include 'comsub.copy'
      include 'comfin.copy'
      include 'comswe.copy'
      include 'comin2.copy'
      include 'comrot.copy'
      include 'comnud.copy'
      include 'compri.copy'
      include 'compla.copy'
      include 'comgdr.copy'
      include 'comfou.copy'
      include 'comber.copy'
      include 'commod.copy'
      include 'comant.copy'
      include 'commmk.copy'
      include 'comequ.copy'
      include 'com3di.copy'
      include 'complo.copy'
      include 'complp.copy'
      include 'comwrr.copy'
      include 'comphy.copy'
      include 'comfft.copy'
      include 'cokpco.copy'  ! ERN
      include 'coequi.copy'  ! ERN
      include 'compow.copy'  ! ERN

      character*4 machin
      logical wridat
      integer dobtyp(2), antyn(maxant), anetyn(maxant), equidn(maxspe+1)
      double precision ktoans

C     "Namelist" statements:
      include 'namgeo.copy'
      include 'nammag.copy'
      include 'nampla.copy'
      include 'namant.copy'
      include 'namreg.copy'
      include 'namsub.copy'
      include 'nampom.copy'
      include 'namsys.copy'
      include 'namout.copy'
      include 'namplo.copy'
      include 'namnud.copy'
      include 'nammes.copy'

      character*(*) iofnam
	parameter(iofnam = 'Data_Files/JET_data_tight.txt')
c	parameter(iofnam = '..\..\Data_Files\JET_dataL.txt')
      character*500 title
	character*3 run_index_s

cccccccccccccccccc ERNESTO cccccccccccccccccccccccc
	real :: T1(2)							! Variables for evaluating
	real :: tempo, tempo0, tempoini			! the elapsed CPU time
	real*8 :: dummy
	integer :: resp
	character(120) :: dummychar
	integer :: OpenStat
ccccccccccccccccccccccccccccccccccccccccccccccccccc

      logical unmame, yes, use_outpow
 
      integer i, i1, i2, j, kr, kt, l, motovs, mst1m, mst2m, mstep, mstsm
     ;, nifil2, nofil2, ispe, ip1, it
     ;, lres
     ;, ith, lbl2, nrhss
     ;, nsum, i1st, isrchfge, ioerr
     ;, idum
     ;, run_index
c     ;, ceilq, nwr

      double precision second
      integer system
      
      double precision 
     ;  rdum
     ;, fcy
     ;, chneu(maxreg), rhoplama

      complex*16 vdum

      external second, nsum, i1st
cERN     ;, long, isrchfge
     ;, longeurs, isrchfge   
       
      data
     ;  eel/1.60217733d-19/,   me/9.1093897d-31/,     mh/1.672623d-27/
     ;, clight/2.99792458d08/, eps0/8.854187817d-12/, mu0/1.2566370614d-06/

c      Conversion factors: input densities are in 10**20 m**-3,
c                          temperatures in eV,
c                          other inputs are in MKSA units.
      data      dfac/1.d20/, tfac/1.60217733d-19/
 
      data      
     ;  czero/(0.d0,0.d0)/, cun/(1.d0,0.d0)/, ci/(0.d0,1.d0)/, ctwo/(2.d0,0.d0)/
     ;, pi/3.1415926535897932385d0/
     ;, sqrt2/1.4142135623730950488d0/, sqrt2i/0.70710678118654752440d0/

      sqrtpi = dsqrt(pi)
      oopi = 1.d0 / pi
      twopi = 2.d0 * pi
 
ccccccccccccccccccccc Begin Cyrano  ccccccccccccccccccccccccc

	print *, '       CYRANO'
	print *, '       ======'
	print *

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


cccccccccccccccccccccccc QLFP ccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        DIRK=.false.                        ! Import data from BATCH
	QLFPDIR = './QLFP_data/runtest/'  ! Where to find the data
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccc Testing M12 matrices for standard coordinates ccccccc
	WRITE_M12STD=.true.
	if(WRITE_M12STD)then
	   resp = SYSTEM('rm -rf ./M12std/*')
c	   resp = SYSTEM('mkdir ./M12std/')
	end if
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     A few default values (in case omitted in namelists):

cERN	All species parabolic (NEW March 2007)
	i_special(:)=0

c     NEW: Don't use artificial damping in HOTTEN
      SWcoef = 0.d0

c     Default: don't use k//=ct coordinates:
      cokpco = .false.
c     D-shape parameters: initialized to circular concentric surfaces
      geneq = .false.
      dshape = .false.
      kappa = 1.
      delta = 0.
      shsh0 = 0.
      z0 = 0.
      updsym = .true.
c     Number of poloidal Fourier coeffs. for smooth (nonresonant) plasma terms:
      nftnr = 4
c     Small k// value to avoid trouble when truly 0:
      kparze = 1.d-8
c     Default: no output for TRANSP code:
      transp = .false.
      notrsp = 25
c     Default: print solver progress information:
      soldoc = .true.
c     Default: don't store element matrices:
      stomat = .false.
c     And don't re-compute element matrices afterwards:
      recalc_total = .false.
cPL      recalc = .false.
      rawel1 = .false.
      onlyab = .false.
      onlyab_now = .false.
      recalc_by_species = .false.
	recalc_by_species_now = .false.
	use_outpow = .true.
       
ccccccccccccc ERN: Non-Maxwellians
      nspgdr=0
	ispgdr(1)=0
      my_nspgdr=0
	my_ispgdr(1)=0
cccccccccccccccccccccccccccccccccc

c     Default: use Maxwellian distributions 
C     (older option related to Raymond's nonmaxwellian dielectric tensor
C      routines NEWTEN and NEWTENI):
      usenma = .false.
C     Idem: Initialize equilibrium distributions to Maxwellian:
      call iset(maxspe+1, 0, equidn, 1)
      v0alph = 0.
      udrift = 0.
        do i = 1, maxexc
        udsant(i) = .false.
        end do
c     Default: antennae assumed different
      samant = .false.
c     Initialize integration schemes to numerical Gaussian quadratures:
c     (respectively for plasma, vacuum and general dielectric response contrib.)
c     Under development:
c     When .false., use analytical integrals of products of basis functions
c     times linear interpolations of equilibrium coefficients.
      giplas=.true.
      gicurl=.true.
      gigdr=.true.

cERN  Default OUTPUT grid (for plotting only)
      OUTGAU  = .true.	! radial GRID at Gauss points + element boundaries
      OUTNPFT = .true.	! poloidal GRID at npfft points

cERNccccccccccccccccccc NEW OPTIONS cccccccccccccccccccccc

c	 1) Options for importing equilibrium profiles -------------------------

	READ_GENERAL  = .false.	! Import general quantities (overwrite NAMELISTS)
	READ_FLUX     = .false.	! Import flux surfaces (R,Z) (FLUSH/ITER)
	READ_JTOR	  = .false.	! Import current density profile Jtor(rho)
	READ_PROFILES = .false. ! Import exp. temperature and density profiles 
							! (profiles for ions not often available)
	SAME_TEMP = .true.		! Use same temperature profiles for all species
	SAME_DENS = .true.		! Use same density profiles for all species

c     N.B.: READ_FLUX is only active if READ_GENERAL = .TRUE.
	if( .not. READ_GENERAL) READ_FLUX = .false.

c      2) COKPCO general options ---------------------------------------------

	WRITE_REPORT = .false.   ! Write progress to the FILE cokpco_report.txt
	WRITE_SCREEN = .false.	 ! Write progress to the SCREEN
	WRITE_OUTPUT = .true.    ! Write computed M12 values to files in ../../M12run/

	WRITE_TABLES = .false.   ! Write M12 tables for further reading (3D runs)
	READ_TABLES  = .false.   ! Read stored tables instead of calc. M12 matrix	

cERN	14/03/05: Option for standard coord. in vacuum (if cokpco = TRUE)
c			 Already included in datafile (NAMGEO)
	STDVAC = .false.		 ! Default: use cokpco coord. in vacuum region

c     N.B.: READ_TABLES is only active if WRITE_TABLES = .TRUE.
	if ( .not. WRITE_TABLES) READ_TABLES  = .false.

c	 3) COKPCO TABLES options ----------------------------------------------
c		NB: Only active for WRITE_TABLES = .true.

	ALLRAD_TABLES = .true.  ! FALSE : calc. M12 matrix only at NODE elements 
	                        ! TRUE  : calc. M12 matrix at ALL radial elements
c      N.B.: If ALLRAD_TABLES = FALSE, then further radial interpolation is needed
		 
	ALLCOS_TABLES = .false. ! FALSE : calc. M12 matrix at 2*NSPEC cosXo values
						    ! TRUE  : calc. M12 matrix at ALL cosXo values (0-100)
c      N.B.: ALLCOS_TABLES = TRUE is only used for generating complete cosXo tables
c	       to compare with theory (Mathematica). The cosXo interpolation is NOT 
c            implemented and the program STOPS after generating the tables.

	FINEAA_TABLES = .true. ! FALSE : calc. M12 matrix at ~200 A values (FAST)
							! TRUE  : calc. M12 matrix at ~1500 A values (SLOW)
c      N.B.: In ANY case interpolation in A values takes place. For cosXo~1, the
c	       finer grid (TRUE) gives better results.

c      4) 'Dynamic' poloidal mode couplings -------------------------------------
c		 Re-compute number of poloidal modes in each element (nmode(iel))

	 dynmod = .false.	! default

c	 Options:
c	 dynopt = 1: Compute nmode(iel) based on surface perimeters and NEW variable 
c				 poloidal length 'pollength' [meters] 
c	 dynopt = 2: Compute nmode(iel) based on exponential function. Maximum value at 
c				 antenna (modva1/modva2) and minimum at axis, given by NEW variables
c				 modax1 and modax2 
	 dynopt = 2
	
	 pollength = 0.02	! Desired poloidal length to be resolved [m]
	 modax1 = -5		! Minimum poloidal mode at axis
	 modax2 = +10		! Maximum poloidal mode at axis
c	 NB: modax1/2 = modva1/2 gives the same as dynmod=FALSE (constant nmode)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
     
C     Messages and short output in standard output file:
cJAC  nofile = 6
      nofile = 603
      open(unit=nofile, file='STDOutput', status='unknown')

cJAC  NEW: File for timing studies (remove some time prints from screen)
      open(unit=605, file = 'Timings.dat', status = 'unknown')
        write(605,*), '       CYRANO (Timings)'
	write(605,*), '       ================'
	write(605,*)

c     Run data are in standard input file:
cJAC  NIFILE = 5
      NIFILE = 512
c     For Mac and PC versions; this file is empty on Cray:
      include 'openin.copy'


c Not in use:     One-dimensional plots are stored in following file:
      no1dpl=10
      write(nofile,*) 'enter main '

C     Replay can only be true after a first run:
      replay = .false.

      run_index = 0
 

 1000 continue
 
c     Time reference at current run start:
c     [NB: time measurement does not work yet on Mac version; 
c     see UNIX extensions in ABSOFT FORTRAN manual]
      timin = second()

c     Run index, used to create separate folders in Plot_data when the input file contains more than one run:
      run_index = run_index + 1
      call INT_TO_STRING2(run_index, run_index_s)


c     Run title and data input:
      read(nifile,*,end=1001,err=1002)title
      write(603,*)'reading namout ; nifile= ',nifile
      read(nifile,namout,end=1001,err=1002)
      if(wridat)write(nofile,namout)
c

c	String coding the type of computer on which you are running:
c	valid values are 'Mac', 'PC' (supersedes namelist input)
	machin = 'PC'
c       Path for plot files: "paplda"
cJAC    paplda = '                  Plot_data/run' // run_index_s // '/'
        paplda = 'Plot_data/run' // run_index_s // '/'
c       ! paplda is 41 characters, right-justified, machine-dependent!
cJAC     paplda modified do 17 characters in complp.copy


c
      read(nifile,namplo,end=1001,err=1002)
      if(wridat)write(nofile,namplo)
      nifil2=nifile
      if( ploold .and. .not.replay )then
c     Read a solution previously stored on unit NOLIFI:
c     Do not store it again:
      keepso=.false.
      nifil2=nolifi
c     read(nifil2,*)oldjon
      end if


ccc   write(nofile,*)'jobsiw=',jobsiw
ccc   write(nofile,*)'dycole=',dycole
c

cERN	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c	Link number of poloidal plot points to number of fft points:
c	(this supersedes namelist input from NAMPLO, and avoids 
c	 poloidal interpolation in OUTPUT routines!)
cERN	 nploth = npfft

cERN	New: introduction of logical 'outnpft' to control poloidal 
c		 OUTPUT grid (similar to 'outgau' for radial points)
c		 if OUTNPFT = TRUE, then nploth = npfft (see OUTGRID.f)
c	For the moment, impose OUTNPFT = TRUE (never interpolate)
	OUTNPFT = .true.
c	ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c 
      write(603,*)'Reading namelist data on file ',nifil2
      read(nifil2,namgeo,end=1001,err=1002)
      write(nofile,*)'namgeo ok'
      read(nifil2,nammag,end=1001,err=1002)
      write(nofile,*)'nammag ok'
      read(nifil2,nampla,end=1001,err=1002)
      write(nofile,*)'nampla ok'
      read(nifil2,namant,end=1001,err=1002)
      write(nofile,*)'namant ok'
      read(nifil2,namreg,end=1001,err=1002)
      write(nofile,*)'namreg ok'
      read(nifil2,namsub,end=1001,err=1002)
      write(nofile,*)'namsub ok'
      read(nifil2,nampom,end=1001,err=1002)
      write(nofile,*)'nampom ok'
      read(nifil2,namsys,end=1001,err=1002)
      write(nofile,*)'namsys ok'


 
cERNccccccccccccccccccccccccccc ERN ccccccccccccccccccccccccccccccccccc

c	Read general equilibrium quantities. 
c	This overwrites many NAMELIST variables given in JetdataL.txt, 
c	in particular r_p = rx0m(1) and r_wall = rx0m(nreg), which are used to 
c     build the abcissa vectors 'abscis' and 'abscno' in Cyrano (TABLES.f).

      if(READ_GENERAL)then
	   
	   open(UNIT = 11, FILE = 'Profiles/general.dat', 
     ;        STATUS = "OLD", IOSTAT = OpenStat, ACTION = "READ")
			if (Openstat > 0) then
				print *, '--> General equilibrium data not found!' 
				print *, '--> Using Namelist information.'
				print *
				geneq = .FALSE.
				READ_FLUX = .false.
				goto 4444 
			end if
			write(nofile,*) 
			write(nofile,*), 
     ;              '----> Importing experimental data (READ_GENERAL = T)'
			write(nofile,*) 
			print *, '--> Importing general data.'
			print *
	        read ( 11,'(A100)'), dummychar
			read ( 11,'(A100)'), dummychar
			read ( 11,'(A100)'), dummychar
			read ( 11,'(A13,I13)')  , dummychar, SHOT	! ERN
			read ( 11,'(A13,I13)')  , dummychar, Nrad_flu	! ERN
			read ( 11,'(A13,F13.8)'), dummychar, instant	! ERN
			read ( 11,'(A13,F13.8)'), dummychar, rx0m(1)	! Namelist
			read ( 11,'(A13,F13.8)'), dummychar, r0	        ! Namelist (axis)  
			read ( 11,'(A13,F13.8)'), dummychar, z0	        ! Namelist (axis)
			read ( 11,'(A13,F13.8)'), dummychar, Rgeo	! ERN		
			read ( 11,'(A13,F13.8)'), dummychar, Zgeo	! ERN
			read ( 11,'(A13,F13.8)'), dummychar, rx0m(nreg) ! Namelist
			read ( 11,'(A13,F13.8)'), dummychar, shsh0	! Namelist
			read ( 11,'(A13,F13.8)'), dummychar, ipl	! Namelist
			read ( 11,'(A13,F13.8)'), dummychar, b0		! Namelist
			read ( 11,'(A13,F13.8)'), dummychar, rfpow	! Namelist
			read ( 11,'(A13,F13.8)'), dummychar, fregag	! Namelist
			read ( 11,'(A13,F13.8)'), dummychar, n0(1,1)	! Namelist
			read ( 11,'(A13,F13.8)'), dummychar, nb(1,1)	! Namelist
			read ( 11,'(A13,F13.8)'), dummychar, t0(1,1)	! Namelist 
			read ( 11,'(A13,F13.8)'), dummychar, tb(1,1)	! Namelist 
	   close (11)
	  	
	   Raxis = r0 ! cERN: in Cyrano r0 is the axis radius	  	
           Zaxis = z0 ! cERN: in Cyrano z0 is the axis elevation

c	   Correction for density values (Cyrano uses 1e+20/m3 instead 1e+19/m3)	
	   n0 = n0 * 1.0d-1
	   nb = nb * 1.0d-1
c	   Update temperatures for all species
	   if (SAME_TEMP) then	
		   t0(2:nspec,1) = t0(1,1)
		   tb(2:nspec,1) = tb(1,1)
	   end if

	end if ! READ_GENERAL = TRUE

4444	continue

	   print *, 'SHOT number   =', SHOT
	   print *, 'N.flux surf.  =', Nrad_flu
	   print *, 'instant(s)    =', instant
	   print *
	   print *, 'Rgeo(m)       =', Rgeo 
	   print *, 'Zgeo(m)       =', Zgeo 
	   print *, 'Raxis(m)      =', r0
	   print *, 'Zaxis(m)      =', z0
	   print *, 'r_p(m)        =', rx0m(1)
	   print *, 'r_wall(m)     =', rx0m(nreg)
	   print *, 'ShafShift(m)  =', shsh0
	   print *
	   print *, 'Ip(MA)        =', ipl*1.d-6
	   print *, 'Bo(T)         =', b0
	   print *, 'P_icrh(MW)    =', rfpow*1.d-6
	   print *, 'f_icrh(MHz)   =', fregag*1.d-6 
	   print *, 'ne_0(1e19/m3) =', n0(1,1)*10 
	   print *, 'ne_b(1e19/m3) =', nb(1,1)*10
	   print *, 'Te_0(keV)     =', t0(1,1)*1.d-3
	   print *, 'Te_b(keV)     =', tb(1,1)*1.d-3
           print *

c	Update subregion limits (and ensure consistent with region boundaries):
	sx0m(0) = rx0m(0)	
	  do i = 1, nreg
	  sx0m(sum(ns(1:i))) = rx0m(i)
	  end do
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cERN	NEW option for using a KAPPA(elongation) profile
c	If kapprof = TRUE, CYRANO uses an elongation profile of the form
c	kappa(rho) = kappa_orig - (kappa_orig-1) * exp(-50*rho^2/ap2),
c	where kappa_orig = kappa is the NAMELIST value and kappa(0) = 1 
c	near the magnetic axis
c	Changes made in: TABLES.f(L.237), GENERA3.f(L.258)

c	Default
	kapprof = .true.
	kappa_orig = kappa

cERN	15/03/05 Store original cokpco value 
c	(cokpco may now be changed inside some subroutines 
c	 like GENERA3, FILANT, etc...)
	cokpco_orig = cokpco

      if(replay)then
c     ~~~~~~~~~~~~~~
c     Repeat the output for another antenna excitation or different
c     plot parameters:
      write(nofile,*)'Repeat the output with different'
     ;            // ' antenna currents:',(tancur(i),i=1,ntoant)
c2004      call repou2

c2004      else
C2004     ~~~~
c2004: new replay logic to compute power deposition etc:
c2004 NB: only works for single toroidal mode (monoto=.true.)!
cPL10/4/2004: Create folder \runxxx\ for current run index (3-digit xxx padded with zeros if required) in Plot_data folder:
	resp = SYSTEM("mkdir " // paplda)

	goto 3003
      end if
c     ~~~~~~
c     Normal problem solving.
c     Define problem structure using raw data.

      circ = .not. (dshape .or. geneq)
c      circ = .not.(dshape .or. geneq .or. cokpco)
c     'circ' stands for circular concentric magn. surfaces:
      if(circ)shsh0 = 0.d0
c     Reference for plots in cylinder (where actual R0 is infinite):
      if(cyl)R0 = 0.d0
c     Major radius at magnetic axis:
      Ra = R0 + shsh0
c     To deal with Shafranov shift and discontinuities, impose rhocur = 
c     first region boundary. All plasma current is taken inside rhocur!
      rhocur = rx0m(1)

c      updsym = .not.geneq .or. (dshape .and. udsk=0.d0)
      updsym = .not.geneq .and. (.not.dshape .or. (Z0.eq.0.d0))

        if(dshape)then
C         Check elongation and Shafranov shift are admissible:
          if(kappa.le.0.)then
          write(nofile,*)'Elongation KAPPA must be positive; run aborted'
          stop
          end if
          if(shsh0 .ge. 0.5*rhocur)then
          write(nofile,*)'Shafranov shift larger than half 
     ;current minor radius; singular geometry; run aborted'
          print *, 'STOP at Cyrano.f line 516!'
          stop
          end if
        end if
      
        if(nspgdr.gt.0)then
c       Impose k//=ct coordinates if general diel. response is used:      
        cokpco = .true.
c       find which species do not require general diel. response:
        j = 0
          do i = 1, nspec
          yes = .false.
            do L = 1, nspgdr
            if(ispgdr(L).eq.i)yes = .true.
            end do
          spegdr(i) = yes
            if(.not.yes)then
            j = j + 1
            ispnod(j) = i
            end if
          end do
        else
          do i = 1, nspec
          ispnod(i) = i
          spegdr(i) = .false.
          end do
        end if
      
c     Impose consistent polsym switch:
      if(.not.(circ .and. cyl))polsym = .false.
c     N.B.: polsym is input from namgeo. 
c           polsym=.false. is allowed in circular cylinder for testing purposes.
c     n.b. B0 is the toroidal induction on magnetic axis (R=RA).
      crown = rx0m(0) .ne. 0.d0
      rcrown = rx0m(0)
      rhowal = rx0m(nreg)

      ap = rhocur
      ap2 = ap * ap
      api = 1.d0 / ap

c     (now J0 is only used in circular: normaliz. of current density profile)
        if(geneq)then
cERN        stop 'geneq not allowed yet'
        j0 = ipl * (alpha + 1.) / (kappa * pi * rhocur ** 2)
c	WARNING: Wrong expression for geneq=T

        else if(dshape)then
        j0 = ipl * (alpha + 1.) / (kappa * pi * rhocur ** 2)
        else
        j0 = ipl * (alpha + 1.) / (pi * rhocur ** 2)
        end if

c       This sets characteristic length rnorm for vacuum cases:
        if(crown)then
        rnorm = rhowal - rcrown
        dobtyp(1) = - 2
        else
        rnorm = rhowal
        dobtyp(1) = - 1
        end if
 
      glovac = .true.
      glocol = .true.
      glo0or = .true.
	rhoplama = 0.d0
c     Compute densities from species fractions:
c     Assume first species is electrons, 
c     and no other negatively charged particle type is present.
c     CHNEU checks neutrality in each plasma region:
      ilaplreg = 0
        do i = 1, nreg
c       get index of last plasma region:
        if(.not.vacuum(i))ilaplreg = i
c       glovac: remains true when vacuum everywhere.
        glovac = glovac .and. vacuum(i)
c       glocol: remains true when all plasma regions are cold.
        glocol = glocol .and. ( coldpl(i) .or. vacuum(i) )
c       glo0or: remains true when all plasma regions are 0 order flr.
        glo0or = glo0or .and. ( vacuum(i) .or. coldpl(i) .or.
     ;                          .not.flrops(i) )
        if(.not.vacuum(i))rhoplama = rx0m(i)
        chneu(i) = 0.d0
          if(.not.vacuum(i))then
          spefra(1,i) = 1.d0
          chneu(i) = zch(1) * spefra(1,i)
            do ispe = 2, nspec
            n0(ispe,i) = spefra(ispe,i) * n0(1,i)
            nb(ispe,i) = spefra(ispe,i) * nb(1,i)
            chneu(i) = chneu(i) + zch(ispe) * spefra(ispe,i)
            end do
          end if
        end do
	  if(.not.glovac)rnorm = rhoplama
      write(nofile,*)'Normalisation length rnorm =', rnorm

      write(nofile,*)'Neutrality check, for each region (0=OK):' 
        do i = 1, nreg
        write(nofile,*)i, chneu(i)
        end do

      rnori = 1.d0 / rnorm
        if(cyl)then
        rnor0 = 0.d0
        r0i = 0.d0
        else
        r0i = 1.d0 / ra
        rnor0 = rnorm * r0i
        r0orn = 1.d0 / rnor0
        end if
      rnor02 = rnor0 ** 2

      glomax= .true.
        do ispe = 1, nspec
c       Charge-to-mass ratios:
        qom(ispe) = (eel * zch(ispe)) / (mh * amass(ispe))
          if(.not.spegdr(ispe))then
c         Initialize type of equilibrium distribution: 
c         NB: the following variables GLOMAX, EQUIDN, EQUIDF exclusively control
c         the use of R.Koch's routines NEWTEN and NEWTENI, *NOT* the species for
c         which a general dielectric response is required.
            if(equidn(ispe).eq.0)then
c           Usual Maxwellian equilibrium: don't use RK's routines
            equidf(ispe) = 'MAXWE'
            else
            glomax= .false.
              if(equidn(ispe).eq.1)then
c             Maxwellian alpha in D (NEWTEN)
              equidf(ispe) = 'CASE1'
              else if(equidn(ispe).eq.2)then
c             Slowing down alpha in D:
              equidf(ispe) = 'CASE2'
              else if(equidn(ispe).eq.3)then
c             Maxwellian D beam in D (NEWTEN)
              equidf(ispe) = 'CASE3'
              else if(equidn(ispe).eq.4)then
c             Slowing down D beam in D (NEWTEN)
              equidf(ispe) = 'CASE4'
              end if
            end if
          end if
        end do

c     When this switch is true, inner boundary conditions are treated as part of first equation block:   
      bcinb1 = bcinb1 .or. .not.(circ.or.crown)
c     Will only be possibly false in circular geom. or with crown.
      
        do 1 i = 0, nreg
c       ----------------
        ip1 = i + 1
        coljum = .true.
        if(i.gt.0)coljum = coldpl(i) .or. vacuum(i) .or. .not.flrops(i)
        if(ip1.le.nreg)coljum = coljum .and.
     ;  ( coldpl(ip1) .or. vacuum(ip1) .or. .not.flrops(ip1) )
C       When equilibrium is not usual Maxwellian,
C       I assume that a FLR expansion never takes place in other species
C       -> cold-like boundary jumps:
        if(.not.glomax)coljum=.true.
 
        rx0(i) = rx0m(i) * rnori
          if(i.eq.0)then
C         Conditions at axis; actual number of constraints depends 
C         on pol. mode index!
            if(bcinb1)then
c           For m nonzero; in circular, for |m|>1:
            ncstr(i) = 3
            end if
          else
          ncstr(i) = 5
          end if
          if((i.eq.0.and.crown) .or. i.eq.nreg)then
            if(coljum)then
            ncstr(i) = 2
            else
            ncstr(i) = 3
            end if
          end if
        passiv(i)=.true.
          if(i.eq.0 .or. i.eq.nreg)then
          i2 = min0(i+1,2)
            if(dobtyp(i2) .eq. -1)then
            rbtyp(i) = 'EAX'
            else if(dobtyp(i2) .eq. -2)then
            rbtyp(i) = 'MET'
            end if
          else
          rbtyp(i)='IPA'
          end if
        betaed(i) = 0.d0
        pedge(i) = 0.d0
          if( i.ne.0 .and. i.ne.nreg )then
            if(plstep(i) .and. .not.(coldpl(i).and.coldpl(ip1)) )then
c         Edge pressure and beta: (or inward pressure jumps at boundaries)
              do ispe = 1, nspec
              if(.not.vacuum(i) .and. .not.coldpl(i)) pedge(i) = pedge(i)
     ;      + dfac * tfac * nb(ispe,i) * tb(ispe,i)
              if(.not.vacuum(ip1) .and. .not.coldpl(ip1)) pedge(i) = pedge(i)
     ;      - dfac * tfac * n0(ispe,ip1) * t0(ispe,ip1)
              end do
            betaed(i) =  2.d0 * mu0 * pedge(i) / b0**2
            end if
          end if
   1    continue
c       --------
 
cERN      call long(nreg,rx0,rl,1)
          call longeurs(nreg,rx0,rl,1)
 
      omegag = 2.d0 * pi * fregag
      k0 = omegag / clight
      k02 = k0 ** 2
      k0rn2 = (k0*rnorm) ** 2
c     Global factor to normalize weak form to (i*Poynting flux):
      glofac = dcmplx(1.d0 / (2.d0 * omegag * mu0), 0.d0)
      
c     Wall skin depth:
      skiwal = dsqrt(2.d0 * walres / (omegag * mu0))
 
c       Antenna characteristics:
        any_falen = .false.
        do 9 i = 1, ntoant
C       Identify active internal boundaries: (antennae)
        rbtyp(irbant(i)) = 'IAC'
        passiv(irbant(i)) = .false.
c       Antenna end type: (short or open circuit)
          if( anetyn(i) .eq. 1 )then
          anetyp(i) = 'SHC'
          else if( anetyn(i) .eq. 2 )then
          anetyp(i) = 'OPC'
		else
		write(nofile,*)'Antenna #', i, ': unknown type of end connection'
		stop
          end if
c       Antenna type :
          if( antyn(i) .eq. 1 )then
c		General poloidal antenna
          antyp(i) = 'POL'
          udsant(i) = udsant(i) .or. abs(thea1(i)+thea2(i)).lt.1.d-8
          else if( antyn(i) .eq. 2 )then
c		General toroidal antenna - never used - must debug!
          antyp(i) = 'TOR'
          udsant(i) = udsant(i) .or. abs(thea1(i)+thea2(i)).lt.1.d-8
          else if( antyn(i) .eq. 3 )then
c         Textor antennae case A2T [Koch, Schliersee 1986]
          antyp(i) = 'A2T'
          anetyp(i) = 'SHC'
          thea1(i) = 0.D0
          thea2(i) = 2.8D0
c		@Check the toroidal location!!!:
          zaa(i) = 0.D0
          dza(i) = 0.17D0
          udsant(i) = .false.
 
          raprio(i) = 0.1d0
          laprio(i) = 147.d-9
          caprio(i) = 279.d-12
   
          else if(antyn(i) .eq. 4)then
c         One of Textor's fast wave antennae [Durodie]
c         These must be antennae #1 and 2.
          antyp(i) = 'FRE'
          anetyp(i) = 'SHC'
          thea1(i) = 0.d0
          thea2(i) = 1.25d0
          dza(i) = 0.13d0
            if(i .eq. 1)then
            zaa(i) = -0.205
            else
            zaa(i) = 0.205
            end if
          udsant(i) = .false.
   
          raprio(i) = 0.d0
          laprio(i) = 230.d-9
          caprio(i) = 260.d-12
   
          else if(antyn(i) .eq. 5)then
c         One of JET's A2 straps, highly idealized.
c         Angular aperture estimated for total strap heigth 1.135m,
c         elongation 1.7 and minor radius 1.3m at strap.
c         These must be antennae #1 to 4.
          antyp(i) = 'JET'
          anetyp(i) = 'SHC'
          thea1(i) = - 0.26d0
          thea2(i) = 0.26d0
          udsant(i) = .true.
          dza(i) = 0.16d0
            if(i .eq. 1)then
            zaa(i) = - 0.682
            else if(i .eq. 2)then
            zaa(i) = - 0.265
            else if(i .eq. 3)then
            zaa(i) = 0.265
            else if(i .eq. 4)then
            zaa(i) = 0.682
            end if
          end if
        udsant(i) = udsant(i) .and. .not.falen(i) .and. updsym
	  if(falen(i))any_falen = .true.
   9    continue

c       SAMANT (input): 
c	  .true. when all the straps in the problem only differ by
c       their toroidal location. In this case, calculation for a single wire
c       antenna.
c       NEFFAN: 'effective' number of different antennae for system solving.
c       Currently, 1 if all straps the same, NTOANT otherwise. 
c	  NEFFAN is used in PROINI, FILRHS, output routines. 
c	  ADDSCR (Faraday shields): to be updated
          if(samant)then
          if(ntoant.gt.1)write(nofile,*)
     ;    'The antenna straps only differ by their toroidal location.'
          write(nofile,*)
     ;   'This will reduce the number of right-hand sides in the linear system.'
          neffan = 1
          else
          neffan = ntoant
          end if

c       Screens are considered active (for mesh generation):
        do i = 1, nscree
        passiv(irbscr(i)) = .false.
        end do
 
      if(.not. glovac)then
c     Particle identification:
      DO ISPE = 1, NSPEC
      if( zch(ispe).eq.-1.d0 .and. amass(ispe).le.5.5d-4 )then
      paname(ispe) = 'Electrons'
      else if( zch(ispe).eq.1.d0 .and. amass(ispe).eq.1.d0 )then
      paname(ispe) = 'Hydrogen'
      else if( zch(ispe).eq.1.d0 .and. amass(ispe).eq.1.000001d0 )then
      paname(ispe) = 'Hbeam'
      else if( zch(ispe).eq.1.d0 .and. amass(ispe).eq.2.d0 )then
      paname(ispe) = 'Deuterium'
      else if( zch(ispe).eq.1.d0 .and. amass(ispe).eq.2.000001d0 )then
      paname(ispe) = 'Dbeam'
      else if( zch(ispe).eq.1.d0 .and. amass(ispe).eq.3.d0 )then
      paname(ispe) = 'Tritium'
      else if( zch(ispe).eq.2.d0 .and. amass(ispe).eq.3.d0 )then
      paname(ispe) = 'Helium^3'
      else if( zch(ispe).eq.2.d0 .and. amass(ispe).eq.4.d0 )then
      paname(ispe) = 'Alpha'
      else if( zch(ispe).eq.10.d0 .and. amass(ispe).eq.20.d0 )then
      paname(ispe) = 'Neon'
      else if( zch(ispe).eq.4.d0 .and. amass(ispe).eq.9.d0 )then
      paname(ispe) = 'Beryllium'
	else if( zch(ispe).eq.15.d0 .and. amass(ispe).eq.40.d0 )then
      paname(ispe) = 'Argon15'
	else if( zch(ispe).eq.16.d0 .and. amass(ispe).eq.40.d0 )then
      paname(ispe) = 'Argon16'
	else
	paname(ispe) = 'Species #' // char(48+ispe)
c	(needs upgrade if nspec > 9...)
      end if
      if( zch(ispe).eq.1.d0 .and. amass(ispe).eq.2.d0 .and.
     ; (equidf(ispe).eq.'case3' .or. equidf(ispe).eq.'case4') )
     ;paname(ispe) = 'D beam'
      end do
      paname(nspec+1) = 'total'
      end if
 
c     Max.order of coupled neighbour poloidal modes:
        if(polsym)then
        ncrot = 0
        klim = 0
        end if
        if(ncrot .gt. npfft/2)then
        write(nofile,*)'Warning: ncrot=',ncrot,' > npfft/2=',npfft/2
	  print*, 'Warning: ncrot=',ncrot,' > npfft/2=',npfft/2
	  print*, 'Stop at CYRANO line 874!'
        end if
	klim = min0(klim, modva2-modva1)
cPL2004: not used globally      ncou = max0(klim, ncrot)
C     Number of retained poloidal couplings: 
c	KLIM from plasma and NCROT from vacuum terms
 
c     Label current run in log file :
      nofil2 = nofile
      write(nofil2,*)title
        if(wridat)then
        write(nofil2,namout)
        write(nofil2,namplo)
        write(nofil2,namgeo)
        write(nofil2,nammag)
        write(nofil2,nampla)
        write(nofil2,namant)
        write(nofil2,namreg)
        write(nofil2,namsub)
        write(nofil2,nampom)
        write(nofil2,namsys)
        end if
      write(nofil2,9000)(betaed(i),i=1,nreg-1)
      write(nofil2,9001)(pedge(i),i=1,nreg-1)
 9000 format('Region edge beta =         ', 10(1x, g13.5))
 9001 format('Region edge pressure (Pa) =', 10(1x, g13.5), 1h )

        if(keepso)then
c       Solution is to be kept on unit NOLOFI:
        nofil2 = nolofi
c       write(nofil2,*) title
c       write(nofil2,namout)
c       write(nofil2,namplo)
        write(nofil2,namgeo)
        write(nofil2,nammag)
        write(nofil2,nampla)
        write(nofil2,namant)
        write(nofil2,namreg)
        write(nofil2,namsub)
        write(nofil2,nampom)
        write(nofil2,namsys)
        end if
 
        if(monoto)then
c       Run with a single toroidal mode:
        ntotor=1
          if(.not.cyl)then
          motov1 = motoan(1)
          motov2 = motoan(1)
	    motovs = 1
          end if
c         NB: in the cylindrical case, axial wavenumber is in KTOAN(1)
        else
c       Define toroidal mode set for 3D run: [not used yet!]
          if(.not.cyl)then
          i = 0
c		Assumes mode range and mode step give exact integer division:
		ntotor = (motov2 - motov1) / motovs
	      if(ntotor.gt.maxtom)then
	      write(nofile,*)'Max. allowed number of toroidal modes exceeded'
	      write(nofile,*)'NTOTOR=', ntotor, ',but MAXTOM=', maxtom, ';I stop.'
	      stop
	      end if
            do j = motov1, motov2, motovs
            i = i + 1
            motoan(i) = j
            end do
          ntotor = i
          else if(ktoans .ne. 0.)then
c         Cylindrical case: ktoan(1), ktoans and ntotor are input
            do j = 2, ntotor
            ktoan(j) = ktoan(1) + (j - 1) * ktoans
            end do
          end if

c       Initialize 3D outputs:
        call dset(nabplo*(3*nficom+4), 0.d0, yout, 1)
cPL not used        call dset(nabplo*ndof, 0.d0, yana, 1)
ccc     call dset(nabplo*(maxthp+1)*8, 0.d0, tab2, 1)
        end if
c     Total (3D) impedance:
      call zset( maxant**2, czero, zimpto, 1 )
C     COEIMP: weight factor for numerical integration yielding
C     3D impedance matrix.
      coeimp = 1.d0
 
c     Open solver files:
        if(soldoc)then
        i = 1
        else
        i = 0
        end if
      call fbtini(nofile, unitgl, 1, rdum, 0
     ;                  , unitrh, 1, rdum, 0
     ;                  , unitem, 1, rdum, 0, i, stomat)

c     Element type identification:
      do it = 1, neltyp
        if(it.eq.1)then
c       Hermite cubic for each field component (+,-,//):
        ibafty(1,it) = 1
        ibafty(2,it) = 1
        ibafty(3,it) = 1
        else if(it.eq.2)then
c       Lagrange quadratic for Erho, 
c       Hermite cubic other field components:
        ibafty(1,it) = 2
        ibafty(2,it) = 1
        ibafty(3,it) = 1
        else if(it.eq.3)then
c       Hermite cubic for each field component (R,Y,phi):
        ibafty(1,it) = 1
        ibafty(2,it) = 1
        ibafty(3,it) = 1
        end if
      end do

c     Basis function identification: 
c     IDERIV is derivative index, =0 for field dof, =1 for first derivative dof
c     inotyp is node type index of each basis function,
c     =1 for left, 2 for inside, 3 for right node.
c     NG, NM, ND: for each basis function set, number of basis functions resp.
c     associated with left, inside and right nodes.
c     Hermite cubic:
      nbaf(1) = 4
      deg(1) = 3
      ideriv(1,1) = 0
      ideriv(2,1) = 1
      ideriv(3,1) = 0
      ideriv(4,1) = 1
      inotyp(1,1) = 1
      inotyp(2,1) = 1
      inotyp(3,1) = 3
      inotyp(4,1) = 3
      ng(1) = 2
      nm(1) = 0
      nd(1) = 2
c     Lagrange quadratic:
      nbaf(2) = 3
      deg(2) = 2
      ideriv(1,2) = 0
      ideriv(2,2) = 0
      ideriv(3,2) = 0
      inotyp(1,2) = 1
      inotyp(2,2) = 2
      inotyp(3,2) = 3
      ng(2) = 1
      nm(2) = 1
      nd(2) = 1

        do i = 1, nficom
          do j = 1, neltyp
          nbafel(i,j) = nbaf(ibafty(i,j))
          end do
        end do

c     Analytical integrals of products of basis functions:
      call intbfp

ccccccccccccccccccccccccccccccccccccc
c     Loop over toroidal modes:     c
ccccccccccccccccccccccccccccccccccccc

      do 3000 imoto = 1, ntotor
c     -------------------------
      i1 = max0( imoto-1, 1 )
      i2 = min0( imoto+1, ntotor )
 
        if(cyl)then
        kphi = ktoan(imoto)
        if(.not.monoto) coeimp = ( ktoan(i2) - ktoan(i1) ) * 0.5
        else
        n = motoan(imoto)
        kphi = n * r0i
        ktoan(imoto) = kphi
        if(.not.monoto) coeimp = ( motoan(i2) - motoan(i1) ) * 0.5
        end if

      kprn = kphi * rnorm
      kprn2 = kprn ** 2
 
      if(mstust.eq.0 .or. iabs(mstust).ge.1.d10)mstust=1
 
        if(autome)then
        write(nofile,*)'Automatic mesh generation currently not working and '
        write(nofile,*)'has been disabled; input mesh is used!'
        autome = .false.
        end if
      
        if(autome)then
C       =============c       Automatic mesh generation - disabled.
          if(ploold)then
c         Mesh was previously generated
          read(nifil2,nammes)
          write(nofile,nammes)
          else 
          mst1m = mstud1
          mst2m = mstud2
          mstsm = mstust
C         This routine to revise for each new application!
          call autcut
          mstud1 = mst1m
          mstud2 = mst2m
          mstust = mstsm
          end if
        if(keepso)write(nofil2,nammes)
	
        else if(imoto.eq.1)then
c       ======================c       Get mesh from input data. The same mesh is used for all toroidal modes!
        sx0(0) = sx0m(0) * rnori
        nsreg = nsum(nreg, ns)
c       Loop over subregions:
          do isubr = 1, nsreg
          sx0(isubr) = sx0m(isubr) * rnori
          istart = i1st(isubr)
          kt = iele(isubr)
c         Subregion is cut in elts.:
          call cuteq( sx0(isubr-1), sx0(isubr), kt, fx0(istart) )
cERN          call long( kt, fx0(istart), fl(istart+1), 1 )
	        call longeurs( kt, fx0(istart), fl(istart+1), 1 )
          end do
cERN        call long( nsreg, sx0, sl, 1 )
            call longeurs( nsreg, sx0, sl, 1 )
        end if
c       =====c     Number of subregions:
      nsreg = nsum(nreg,ns)
c	Number of radial finite elements:
      nele = nsum(nsreg, iele)
      if(nele.gt.maxnel)stop 'nele too large'
      if(nsreg.gt.maxsub)stop 'nsreg too large'
      if(nreg.gt.maxreg)stop 'nreg too large'


cPL29/5/04 element and basis function identification moved up from here     
 
        do isubr = 1, nsreg
        kt = iele(isubr)
          if(istyp(isubr) .eq. 1)then
c	    Use Hermite cubic basis functions on E+, E-, E// electric field components
          styp(isubr) = 'HEC'
          iddl(isubr) = 2 * ndof
          iconn(isubr) = ndof
          else if(istyp(isubr) .eq. 2)then
c	    Use Lagrange quadratic on Erho and Hermite cubic basis functions on Etheta, Ephi
          styp(isubr) = 'M23'
          iddl(isubr) = 2 * ndof - 1
          iconn(isubr) = ndof - 1
          else if(istyp(isubr) .eq. 3)then
	    print *, 'New element type CAR under construction'
	    stop
c	    Use Hermite cubic basis functions on ER, EY, Ephi electric field components
          styp(isubr) = 'CAR'
          iddl(isubr) = 2 * ndof
          iconn(isubr) = ndof
          end if
        ibub(isubr) = iddl(isubr) - 2 * iconn(isubr)
        write(nofile,2001) isubr, kt, styp(isubr)
        end do
 
c     Generate some index lists for the run:
cERN
c	NB: The input option .false. makes the presol routine run as
c	    in previous version, without re-computing the poloidal mode couplings 
      call presol(.false.,0)

	if(.not.glovac)i_last_plasma_block = ilael(ilaplreg) + 2*(ilaplreg-1)

      if(imoto.eq.1)then  ! Set of actions only required at first toroidal mode of current run
c     =================c

cPL9/4/2004: Erase folder Plot_data and all contents if first run of this data file:
	  if(run_index .eq. 1)then
          resp = SYSTEM("rm -rf Plot_data")
cPL9/4/2004: Create new Plot_data folder:
	  resp = SYSTEM("mkdir Plot_data")
	  end if
cPL9/4/2004: Create folder /runxxx/ for current run index (3-digit xxx padded with zeros if required) in Plot_data folder:
	  resp = SYSTEM("mkdir " // paplda)

cJAC    Open LOGFILE 
	open (UNIT = 7771, FILE = paplda // 'CYRANO_resume.dat', 
     ;      STATUS = "REPLACE", ACTION = "WRITE")


ccccccccccccccccccc ERNESTO cccccccccccccccccccccccccc
c      Initialize sine and cosine tables for various FFT using npfft points:
c      write(6,*)'Before init fft: npfft=', npfft
cJAC   Don't use IMSL in JAC (Extracted routines re-included)

c       call cfft2(1, 1, npfft, vdum, cworkp, vdum) 
c       call cfft2(1, -1, npfft, vdum, cworkm, vdum) 
	call dfftci(npfft, cwork2) ! -- IMSL (complex FFT) --	
c       call rcfft2(1, 1, npfft, rdum, workp, vdum) 
c       call rcfft2(1, -1, npfft, rdum, workm, vdum)
	call dfftri(npfft, work2) ! -- IMSL (real FFT) --	
       write(603,*)'After init fft: npfft=', npfft
		 
c      For non resonant terms, will use (input) nftnr points:
       nftnr = min0(nftnr, npfft)
c       call cfft2(1, 1, nftnr, vdum, worknr, vdum) 
	call dfftci(nftnr, worknr2) ! -- IMSL (complex FFT) --	

ccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     Generate equilibrium tables:
c     NB: 
c     - Don't call any geometrical or magnetic equil.routine before this call!
c     - Automatic mesh generation does not work now; found too hard to be 
c       reliable! 

          tempo0 = second()-timin
c	  print *,'--> (CYRANO.f) Begin TABLES:', tempo0, 'sec' 
          print*, 'Enter TABLES'
          write(605,*), '--> (CYRANO.f) Begin TABLES:', tempo0, 'sec'

          call tables	! Generate equilibrium tables

          print*, 'Enter OUTGRID'
	  call outgrid  ! Generate 2D OUTPUT grid (accord. to OUTGAU and OUTNPFT)
			! NB: This grid will be used in ALL 2D output routines
          print*, 'Enter OUTTAB'
	  call outtab	! Save tables to files

          tempo = second()-timin

c         print *, '--> (CYRANO.f)   End TABLES:', tempo, 'sec' 
c	  print *, '     Time used:', tempo-tempo0, 'sec'
          print*, 'Exit  TABLES'
	  print *
          write(605,*), '--> (CYRANO.f)   End TABLES:', tempo, 'sec' 
	  write(605,*), '     Time used:', tempo-tempo0, 'sec'
	  write(605,*)

cc	  print *, 'STOP at Cyrano.f line 1100'
cc	  stop
	  	
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cERN
c     Characteristic poloidal length (mesh) in meters:
c	pollength = (surface perimeter) / (number pol. modes)
        if(.not.giplas .and. dynmod)then
           write(nofile,*)'dynmod=T not yet compatible with giplas=F. I stop.'
           print *, 'STOP at Cyrano line 1210!'
           stop
	end if

c       Re-generate the index lists computing the necessary number of 
c	poloidal modes for each element nmode(iel):
c	The second argument (dynopt) sets the type of nmode(iel) dependence:
c	    1 - Compute nmode(iel) based on surface perimeters and NEW variable 
c		poloidal length 'pollength' [meters] (Cyrano line 255)
c	    2 - Compute nmode(iel) based on exponential function. Maximum value at 
c		antenna (modva1/modva2) and minimum at axis, given by NEW variables
c		modax1 and modax2 (Cyrano line 256)

        if(dynmod) then
	   call presol(dynmod,dynopt)
	end if

! ---------- Reset output folders for k// = const. coordinates ----------------

cJACc     1) Tables (M12tables)	
cJAC	if (cokpco .and. WRITE_TABLES) then
cJAC	      resp = SYSTEM("rm -rf M12tables")
cJAC	      resp = SYSTEM("mkdir M12tables")
cJAC	end if
cJACc     2) Runs (M12run) 
cJAC	if (cokpco .and. WRITE_OUTPUT) then
cJAC	      resp = SYSTEM("rm -rf M12run")
cJAC	      resp = SYSTEM("mkdir M12run")
cJAC	end if

! ---------- Generate M12 tables for k// = const. coordinates ---------------- !

	if (cokpco .and. WRITE_TABLES) then	
           print*,'(Cyrano) COKPCO not reqdy!'
cJAC	! Generate and write the M12 tables	to folder ../../M12tables/
cJACc	  tempo0 = etime (T1)
cJAC	  tempo0 = second()-timin
cJAC          print *,'--> (CYRANO.f) Begin M12COKPCO_TABLES:', tempo0, 'sec' 
cJAC
cJAC	  call M12COKPCO_TABLES
cJAC
cJACc	  tempo  = etime (T1)
cJAC         tempo  = second()-timin        
cJAC          print *, '--> (CYRANO.f)   End M12COKPCO_TABLES:', tempo, 'sec' 
cJAC	  print *, '     Time used:', tempo-tempo0, 'sec'
cJAC	  print *
cJAC	end if ! WRITE_TABLES = TRUE
cJAC	WRITE_TABLES = .false.
cJAC	if (ALLCOS_TABLES) then
cJAC           print *, 'Stop at Cyrano line 1258.'
cJAC           stop
        end if
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cERNc       Write data to plot coordinate lines:
cERN        call coordp
cERN		  This call has been substituded by OUTGRID and OUTTAB 
cERN		  (called just after TABLES.f)

cccccccccccccccccccccc QLFP input cccccccccccccccccccc

cNEW	  logical DIRK variable for QLFP input
c	  DIRK = .true.

        if(.not.glovac)then
c       Read equilibrium distribution(s), output from QLFP code:
c       For Mac and PC versions; this file is empty on the Cray:
c        include 'openfp.copy'

         if(nspgdr.gt.0 .or. my_nspgdr.gt.0)then
	
     	        if(DIRK)then
cERN		   call read_QLFPheader  ! use Dirk's data in QLFPDIR folder
		else
                   open(unit=unqlfp, file='QLFP_data.txt', status='old')
cERN		   call read_QLFP        ! Generate numerical data in Cyrano (QLFPdata.txt)
		end if

          end if

cccccccccccccccccccccccccccccccccccccccccccccccccccccccc








c       Detect presence of a fundamental cyclotron resonance or a '2nd' harmonic
c       in the vessel (i.e. not just inside the plasma).
c       (Limited validity: toroidal geometry, and Doppler shift ignored here)
c       Does not detect Landau-type resonance, much less peaked than cyclotron.
c       'spres' and 'res' are later used to select the type of poloidal fft for
c        corresponding contributions to dielectric tensor 
c        (smaller number of points is used for nonresonant terms)
          do lres = -2, 2
          res(lres) = .false.
            do ispe = 1, nspec
	      fcy = (omegag - LRES * qom(ispe) * bmin(nabsci)) 
     ;          * (omegag - LRES * qom(ispe) * bmax(nabsci))
            spres(lres,ispe) = fcy .le. 0.
            res(lres) = res(lres) .or. spres(lres, ispe)
c	        if(lres.ne.0 .and. spres(lres,ispe))then
c	        lev = omegag / (dfloat(lres)*qom(ispe))
c             Approximate location of resonance layer
c                do i = 1, nabsci
c	          
c	          end do
c	        end if
            end do
          end do
        end if
c       Locate QLFP radial abscissae in Cyrano table ABSCIS:
c       NB: one can have nspgdr=0 but nraddr > 0 (for tests or initialization of
c       wave-QLFP iteration). When nspgdr=0, radial abscissae for QLFP are read
c       in input file (NAMPLA).
          if(nraddr.gt.0)then
            do kr = 1, nraddr
            irafp(kr) = isrchfge(nabsci, abscis, 1, raqlfp(kr))
              if(irafp(kr).gt.nabsci)then
              write(nofile,*)
     ;        'Warning: illicit radius for QLFP output: index=', kr
     ;,       '; radius=', raqlfp(kr)
     ;,       '. Using authorized maximum instead:', abscis(nabsci)  
              irafp(kr) = nabsci
              end if
            end do
          end if

c     Study the local dispersion for Maxwellians and first toroidal mode only:
      open(unit=81
     ;, file = paplda // 'Tensor elements', status='unknown')
c     [See what to do with general diel response!]
c     First switch off non-Maxwellian behaviour (for control of RK's routines):
      unmame = usenma
      usenma = .false.
      call groots
      usenma = unmame
      end if
C     =====
      if(distes)goto 3001

c     Compute antennae spectra: 
c     ------------------------

cERN	14/03/05: If STDVAC = true, use standard coordinates
c                 for antenna spectrum (restore value at the end)

          tempo0 = second()-timin
	  print *,'Enter FILANT' 
          write(605,*), '--> (CYRANO.f) Begin FILANT:', tempo0, 'sec'

	  if(STDVAC)cokpco=.false.
        
          call filant
	
          cokpco = cokpco_orig
  
          tempo = second()-timin
          print *, 'Exit  FILANT' 
c	  print *, '     Time used:', tempo-tempo0, 'sec'
	  print *
          write(605,*), '--> (CYRANO.f)   End FILANT:', tempo, 'sec' 
	  write(605,*), '     Time used:', tempo-tempo0, 'sec'
	  write(605,*)

       
        if(ploold)then
C@ TO SEE AGAIN!
c  Probably obsolete, now that all plots are dealt with in Matlab...
c       Case of plot a previously stored run:
        write(603,*)'old solution to be read in file ',nifil2
c       read(nifil2,1003)
c    ;   ((x(j,k), j=1,ngelim ), k=1,ntoant )
c       write(6,*)'old solution read in file ',nifil2
        goto 3001
        end if
 
c     Jumps across antennae:
      print *,'Enter ANJUMP'
      call anjump
      print *,'Exit  ANJUMP'
      print *


ccccccccccccccccccccccccc GENERA3 cccccccccccccccccccccccccccccc

        tempo0 = second()-timin
	print *,'--> Enter GENERA3 (gensys = T)'
	print *
	write(605,*),'--> (CYRANO.f) Begin GENERA3 (gensys = T):', tempo0, 'sec' 
	write(605,*)
c       Ensure both active and reactive terms are computed:
	onlyab_now = .false.
	recalc_by_species_now = .false.
c       Generate and factorize global matrix and right-hand side; solve system:
        gensys = .true.      

        call genera3

	tempo  = second()-timin
        print *
        print *, '--> Exit GENERA3 (gensys = T)'
c	print *, '     Time used:', tempo-tempo0, 'sec'
	print *
        write(605,*), '--> (CYRANO.f)   End GENERA3 (gensys = T):', second()-timin, 'sec' 
	write(605,*), '     Time used:', tempo-tempo0, 'sec'
	write(605,*)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 3001 continue  ! jump point for case of dispersion calculation only
 3000 continue  ! imoto (toroidal mode loop)
c     --------
 
      if(distes)goto 3002

c     Close solver matrix file:
      call fbtclo('MAT')


      if(wrisol)then
c     Solution vector:
      write(nofile,*)' '
      write(nofile,*)'Solution vector:'
      write(nofile,*)'----------------'
      imoto = 1
      i2 = 0
        do ith = 1, nthoma
        call rdsol(imoto, ith, bel, ldbel, lbl2, nrhss)
c        nwr = ceilq(lbl2,4)
c        do i = 1, nwr
          do i = 1, lbl2
          i2 = i2 + 1
          write(nofile,100)i2, bel(i,1)
 100      format(1h ,4(i4,2x,g12.4,2x,g12.4))

c          write(nofile,100)(j+i-1,bel(j+i-1,1),j=1,lbl2+1-i,nwr)
c 100      format(1h ,4(i4,2x,g11.3,2x,g11.3,))
c 100      format(1h ,4(i4,', (',g11.3,';',g11.3,') '))

          end do
        end do
      end if
 
 3002 continue  ! jump point for case of dispersion calculation only



c     Generate output and plot files:
c     Open instructions (computer dependent):
cERN      include 'openpl.copy'

c     Dispersion plots (eq. profiles in OUTTAB):
cERN	call outdis   ! OUTPUT FILES also commented in OPENPL.COPY
	print *,'Enter OUTDIS'
	call outdis_ern
	print *,'Exit OUTDIS'
        print *
c        print *, '    --> after outdis OK'
c        write(605,*), '    --> after outdis:', second()-timin 

 3003 continue  ! jump point for case of 'replay' the previous run 
c				with different antenna currents

        if(.not.distes)then



c       Element power balance data:
c       [NB: solution vector is modified in OUTRFF when there are Faraday shields
c       This case to be dealt with later!]
          if(stomat)then
	    onlyab_now = onlyab
	    recalc_by_species_now = .false.
	    call qufetc2(.false., idum, .true., .false., .true., idum)
          end if



          if(recalc_total)then
c         This re-computes all finite element matrices to get total power balance
c         on each element:
ccccccccccccccccccccccccc ERNESTO cccccccccccccccccccccccccccccc
c	    tempo0 = etime (T1)
c	    print *,'--> (CYRANO.f) Begin GENERA3 (gensys = F, ' //
c     ;            'recalc total element powers):', tempo0, 'sec'
c	    print *

          gensys = .false.
c         Select whether only active terms are recomputed (onlyab is input):
	    onlyab_now = onlyab
	    recalc_by_species_now = .false.
cERN		Don't write results in this ASPLAS call
          WRITE_OUTPUT = .false.
          call genera3

c	    tempo  = etime (T1)
c          print *, '--> (CYRANO.f)   End GENERA3 (gensys = F):', tempo, 'sec' 
c	    print *, '     Time used:', tempo-tempo0, 'sec'
c	    print *
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          end if

          if(.not.glovac .and. recalc_by_species)then
c         This re-computes all finite element matrices species by species
c         to get power balance species by species on each element:
ccccccccccccccccccccccccc ERNESTO cccccccccccccccccccccccccccccc
c	    tempo0 = etime (T1)
c	    print *,'--> (CYRANO.f) Begin GENERA3 (gensys = F, ' //
c     ;            'recalc element powers for each species):', tempo0, 'sec'
c	    print *,'--> (CYRANO.f) Begin GENERA3 (gensys = F):', tempo0, 'sec' 
c	    print *

          gensys = .false.
c         Select whether only active terms are recomputed (onlyab is input):
	    onlyab_now = onlyab
	    recalc_by_species_now = .true.
cERN		Don't write results in this ASPLAS call
          WRITE_OUTPUT = .false.
          call genera3

c	    tempo  = etime (T1)
c          print *, '--> (CYRANO.f)   End GENERA3 (gensys = F):', tempo, 'sec' 
c	    print *, '     Time used:', tempo-tempo0, 'sec'
c	    print *
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          end if

ccccccccccccccccccccccccc ERNESTO cccccccccccccccccccccccccccccc
c	  tempo0 = etime (T1)
          tempo0 = second()-timin
	  print *,'--> Begin OUTPUT routines:'
          print *
	  write(605,*),'--> (CYRANO.f) Begin OUTPUT routines:', tempo0, 'sec'

c       Write OUTPUT data:
cERN	NB: the 2D (R,Z) grid was already built in OUTGRID.f (line ~1165) 

c       RF fields:
	if(STDVAC)cokpco=.false.
	print *,'Enter OUTRFF'
        call outrff
	cokpco = cokpco_orig 
        print *,'Exit OUTRFF'
        print *
        write(605,*), '    --> after outrff:', second()-timin 

c       Power:
        if(use_outpow)then
           print *,'Enter OUTPOW (this may take some time ....)'
           call outpow
           print *, 'Exit OUTPOW'
           print *
           write(605,*), '    --> after outpow:', second()-timin 
        end if

ccccccccccccccccccc
cERN	  if(WRITE_M12STD)call outpow_m12std
c	  if(WRITE_M12STD)call outpow_m12std_2
cccccccccccccccccccc

c       2D fields: Must call after outpow if TRANSP output required (cf normalization factor norfac):
        if(plot2d)then
	  print *,'Enter OUTRF2'
          call outrf2
	  print *, 'Exit OUTRF2'
          print*
          write(605,*), '    --> after outrf2:', second()-timin 
        end if


c       Misc.:
        print *,'Enter OUTBIT'
        call outbit
	print *, 'Exit OUTBIT'
        print *
        write(605,*), '    --> after outbit:', second()-timin       

cccccccccccccccccccccccc
        if((nspgdr.ne.0 .or.my_nspgdr.ne.0 ) .and. DIRK)then
            print*, 'Enter OUT_RFCOEF'
cx            call OUT_RFCOEF(ispgdr(1))
c            call OUT_RFCOEF_2(ispgdr(1))
            print*, 'Exit OUT_RFCOEF'	    
            write(605,*), '     --> after OUT_RFCOEF' , second()-timin  
        end if
cccccccccccccccccccccc

        tempo = second()-timin
        print *, '--> End of OUTPUT routines!' 
c	print *, '     Time used:', tempo-tempo0, 'sec'
	print *
        write(605,*), '--> (CYRANO.f)   End of OUTPUT:', tempo, 'sec' 
	write(605,*), '     Time used:', tempo-tempo0, 'sec'
	write(605,*)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        end if


c2004 Keep solver files open in case a 'replay'case follows in input file
c2004c     Close solver system:
c2004      call fbtclo('ALL')
 
c     if( .not. monoto )
c    ;write(nofile,*)'Global (3D) impedance matrix: ',
c    ; ((zimpto(i,j),j=1,ntoant),i=1,ntoant)
 
c2004      end if
c     ~~~~~~
      write (nofile,*) 'Exit time= ',second()-timin
 


cJAC    Write a LOGFILE of the simulation -----------------------------
         
c     1.2) Species density	
c	open (UNIT = 7771, FILE = paplda // 'CYRANO_resume.dat', 
c     ;      STATUS = "REPLACE", ACTION = "WRITE")
		write(7771,*) '   '
		write(7771,*) 'CYRANO: Main parameters'
		write(7771,*) '-----------------------' 		
		write(7771,*) '   '
                if(READ_FLUX)then
                   write(7771,*) 'SHOT number   =', SHOT
                   write(7771,*) 'instant(s)    =', instant
                   write(7771,*) 'N.flux surf.  =', Nrad_flu
                write(7771,*)
                end if
                write(7771,"(A16, f12.3)") ' Rgeo(m)       =', Rgeo 
                write(7771,"(A16, f12.3)") ' Zgeo(m)       =', Zgeo 
                write(7771,"(A16, f12.3)") ' Raxis(m)      =', r0
                write(7771,"(A16, f12.3)") ' Zaxis(m)      =', z0
                write(7771,"(A16, f12.3)") ' r_p(m)        =', rx0m(1)
                write(7771,"(A16, f12.3)") ' r_wall(m)     =', rx0m(nreg)
                write(7771,"(A16, f12.3)") ' ShafShift(m)  =', shsh0
                write(7771,"(A16, f12.3)") ' Ip(MA)        =', ipl*1.d-6
                write(7771,"(A16, f12.3)") ' Bo(T)         =', b0
                write(7771,"(A16, f12.3)") ' P_icrh(MW)    =', rfpow*1.d-6
                write(7771,"(A16, f12.3)") ' f_icrh(MHz)   =', fregag*1.d-6 
                write(7771,"(A16, f12.0)") ' N             =', motoan(1)
                write(7771,"(A16, f12.3)") ' ne_0(1e19/m3) =', n0(1,1)*10 
                write(7771,"(A16, f12.3)") ' ne_b(1e19/m3) =', nb(1,1)*10
                write(7771,"(A16, f12.3)") ' Te_0(keV)     =', t0(1,1)*1.d-3
                write(7771,"(A16, f12.3)") ' Te_b(keV)     =', tb(1,1)*1.d-3
                write(7771,*) ' '  
                write(7771,*) ' ' 
                write(7771,*) '0D results:'
                write(7771,*) '-----------'
                write(7771,*) ' ' 
                write(7771,"(A16, 6A15)")   ' Species:', (paname(j), j = 1,nspec)
                write(7771,"(A16, 6f14.3)") ' n(0) [1e19/m3]:',(10.0*dentab(1,j), j=1,nspec) !10*dentab(1,2), 10*dentab(1,3) 
                write(7771,"(A16, 6f14.1)") ' FLR0 power [%]:',(100.0*parpow(j,1)/parpow(nspec+1,1), j=1,nspec)
                write(7771,"(A16, 6f14.1)") ' FLR2 power [%]:',(100.0*parpow(j,2)/parpow(nspec+1,2), j=1,nspec)
                write(7771,*) ' ' 
                write(7771,"(A19, f8.2, A1)") 'Total time used  =', second()-timin, 's'

      close(7771)


c       ---------------------------------------------------------------------




cccccccccccccc Copy some files to /Plot_data/runXXX ccccccccccccc
c	resp = SYSTEM("copy ..\..\STDOutput " // paplda // 
c     ;              "STDOutput_" // run_index_s // ".txt")
c	resp = SYSTEM("copy " // iofnam // paplda // "datafile_" // 
c     ;               run_index_s // ".txt")
c	if (READ_GENERAL) then
c	resp = SYSTEM("copy ..\..\Profiles\general.dat " // paplda // 
c     ;              "general_" // run_index_s // ".txt")
c	end if
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


ccccccccccccccccccccccccc END message cccccccccccccccccccccccccccccc          
          
          tempo  = second()-timin
	  print *
c	  call beepqq(1000,200)
	  print *, '---------------> Cyrano FINISHED!'
	  print *, '                 Total time used:', tempo, 'sec'
	  print *
	  write(605,*)	  
          write(605,*), '---------------> (CYRANO.f) Cyrano FINISHED!'
	  write(605,*), '                  Total time used:', tempo, 'sec'
	  write(605,*)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

          close(605)

      goto 1000

 1001 continue
      write(nofile,*)'End of namelist data'
c2004     Close solver system:
      call fbtclo('ALL')
      stop

 1002 continue
      write(nofile,*)'Error in namelist data'
c2004     Close solver system:
      call fbtclo('ALL')
      stop
 
 200  format(1h ,(12(g10.2,1x)))
 300  format(1h ,3g16.8)
 222  format(1h ,5(1x,i4,1x,g16.8))
 2001 format(1h ,'Subregion ',i1,': ',i3,' elements of type ',a3)
 1003 format(10g23.15)

      end

