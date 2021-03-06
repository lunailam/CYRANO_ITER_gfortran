      subroutine powabs(irad, onespe, isp, gdrtoo, temp_total)

      implicit none
      logical onespe, gdrtoo
      integer irad, isp
	     real temp_total

C     Symmetrized version, with resonant absorption only.
C     Should work in general geom. as well.
c     Uses IELM in comswe; yinv, si, co, rhoinv, phi, bmodul, copisi
c     Computes plasma absorption. this routine is consistent
c     with routine hotten for "cold" options.
c     Collisional effects added as perturbation of coll.less calculation
c     Absorption on slowing-down distribution by calls to newteni
c     Does not deal with general dielectric response!
      
      include 'pardim.copy'
      include 'dynou2.copy'
      include 'compow.copy'
      include 'comgeo.copy'
      include 'commag.copy'
      include 'comma2.copy'
      include 'compla.copy'
      include 'commod.copy'
      include 'comswe.copy'
      include 'compri.copy'
      include 'comfin.copy'
      include 'comfic.copy'
      include 'comphy.copy'
 
      logical inscr, insc

      integer 
     ;  ispe, im, nmim, nmi2, ipm, kstop1, kstop2, m1, m2, imu
     ;, spli(maxspe), is, ncl, i

      double precision 
     ;  fa(14,2*maxpom-1), kpar(2*maxpom-1), ketan(maxpom), tr, second
     ;, xmu(-2:2), ex(-2:2), omc, coeff, coef1, coef2, small, dam2, ntni
     ;, cabs2, drpc, aux, myfa(2*maxpom-1), myfa2(2*maxpom-1)

      complex*16 
     ;  pmder(15,maxpom), zdotc, zdotu, tra(15)
     ;, ks(3,3), efield(3), ytra(3), zzz, z1, z2, mypmder(maxpom)
c	 ;, omcol
ccccccccccccccccc teste do tempo ccccccccccccccccccccc
    real :: start, finish
ccccccccccccccccc teste do tempo ccccccccccccccccccccc
ccccccccccccccccccc ERNESTO (cokpco) cccccccccccccccccccccccccccccc
	real :: T1(2)	          !   Variables for evaluating
	real :: tempo1, tempo0       !   the elapsed CPU time
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      external second, inscr, zdotc, zdotu

      cabs2(zzz) = dreal(zzz*dconjg(zzz))
	drpc(z1, z2) = dreal(z1 * dconjg(z2))

c      write(6,*)'enter powab23, irad= ',irad,'; intabp= ',intabp
      
C     Init. all contrib. to 0:       
c     inicio contagem de tempo
      call cpu_time(start)
	  
      call dset(9, 0.d0, powvec, 1)

      if(gdrtoo)then
        if(onespe)then
        ncl = 1
        spli(1) = isp
        else
        ncl = nspec
          do i = 1, ncl
          spli(i) = i
          end do
        end if
      else
        if(onespe)then
        if(spegdr(isp))return
        ncl = 1
        spli(1) = isp
        else
        ncl = nspec - nspgdr
        call icopy(ncl, ispnod, 1, spli, 1)
        end if
      end if
	  
ccccccccccccccccccccc LUENNE ccccccccccccccccccccc      
      call cpu_time(finish)
	  temp_total = temp_total + (finish - start)
cccccccccccccccccccccccccccccccccccccccccccccccccc
	  
      insc = inscr(ireg)
c      write(6,*)'ntnt= ',ntnt(intabp),'; yinv= ',yinv,'; kprn= ',kprn
c      write(6,*)'r0or= ',r0or,'; si= ',si
      ntni = 1.d0 / ntnt(intabp)
        if(.not.cokpco)then      
          do im = 1, nmode(ielm)
          m = im - 1 + minf(ielm)
c         OK in CYL too:
cERN	12/05/05 : Problems at axis: yinv = 1e+10 -> ketan = m*1e+10
c			  Trying 2nd abscissa point
		if(irad .eq. 1) yinv = 1.0d0/(abscis(2))

          ketan(im) = m * yinv * ntni * co - kprn * r0or * si
          end do
	  else
          do im = 1, nmode(ielm)
          m = im - 1 + minf(ielm)
c          ketan(im) = ((m + n * ) * hachi() * co - n * ri) / si
          end do
	  end if

      do im = 1, nmode(ielm)
c     Modal quantities:
 
C     Nabla+ E+
      pmder(1,im) = sqrt2i*(xpmp(2,irad,im) - ketan(im) * xpmp(1,irad,im))
C     Nabla- E+
      pmder(2,im) = sqrt2i*(xpmp(2,irad,im) + ketan(im) * xpmp(1,irad,im))
C     Nabla+ E-
      pmder(3,im) = sqrt2i*(xpmp(4,irad,im) - ketan(im) * xpmp(3,irad,im))
C     Nabla- E-
      pmder(4,im) = sqrt2i*(xpmp(4,irad,im) + ketan(im) * xpmp(3,irad,im))
C     Nabla+ E//
      pmder(5,im) = sqrt2i*(xpmp(6,irad,im) - ketan(im) * xpmp(5,irad,im))
C     Nabla- E//
      pmder(6,im) = sqrt2i*(xpmp(6,irad,im) + ketan(im) * xpmp(5,irad,im))
 
C     E+
      pmder(7,im) = xpmp(1,irad,im)
C     E+'
      pmder(8,im) = xpmp(2,irad,im)
C     E-'
      pmder(9,im) = xpmp(4,irad,im)
C     E//
      pmder(10,im) = xpmp(5,irad,im)
C     E//'
      pmder(11,im) = xpmp(6,irad,im)
C     E+''
      pmder(12,im) = xpmp2(1,irad,im)
C     Nabla+ E-'
      pmder(13,im) = xpmp2(2,irad,im) - ketan(im) * xpmp(4,irad,im)
C     E//''
      pmder(14,im) = xpmp2(3,irad,im)
C     E-
      pmder(15,im) = xpmp(3,irad,im)

c     B//
      mypmder(im) = hrtp(3,irad,im)   !! Bphi for NOW!!!

c      write(6,*)im
c      write(6,*)(j, pmder(j,im),j=1,15) 
      end do
 
      if(ipl.eq.0.d0)then
C     Case of no poloidal field
C     ~~~~~~~~~~~~~~~~~~~~~~~~~
      kpar(1) = kphi * co * r0or
      if(dabs(kpar(1)).lt.kparze)kpar(1) = dsign(kparze,kpar(1))
      nmim = 1
      nmi2 = 1
c     Local fields and derivatives, stored in PMDER:
        do ipm = 1, 15
        tra(ipm) = zdotu(nmode(ielm), pmder(ipm,1), 15, copisi, 1)
        end do
      call zcopy(15, tra, 1, pmder, 1)
      efield(1) = pmder(7,1)
      efield(2) = pmder(15,1)
      efield(3) = pmder(10,1)
 
c     Combine local fields for absorption calc:
c     do im = 1, nmim
      im = 1
      fa(1,im) = cabs2( pmder(7,im) )
      fa(2,im) = cabs2( pmder(10,im) )
      fa(3,im) = dimag( pmder(7,im) * dconjg(pmder(5,im)) )
      fa(4,im) = dimag( conjg(pmder(10,im))
     ;                  * (pmder(3,im)-pmder(2,im)) )
      fa(5,im) = cabs2( pmder(1,im) )
      fa(6,im) = cabs2( pmder(2,im) )
      fa(7,im) = cabs2( pmder(3,im) )
      fa(8,im) = drpc(pmder(2,im), pmder(3,im))
      fa(9,im) = cabs2( pmder(5,im) )
      fa(10,im) = cabs2( pmder(6,im) )
      fa(11,im) = cabs2(pmder(8,im)) + drpc(pmder(12,im), pmder(7,im))
      fa(12,im) = sqrt2 * drpc(pmder(3,im), pmder(8,im))
     ;                  + drpc(pmder(13,im), pmder(7,im))
      fa(13,im) = cabs2(pmder(11,im)) + drpc(pmder(14,im), pmder(10,im))
      fa(14,im) = cabs2(pmder(3,im)-pmder(2,im))
c     end do
 
      else
C     With poloidal field:
C     ~~~~~~~~~~~~~~~~~~~
      nmim = 2 * nmode(ielm) - 1
      nmi2 = nmode(ielm)

        do ipm = 1, 15
          do im = 1, nmode(ielm)
          pmder(ipm,im) = pmder(ipm,im) * copisi(im)
          end do
        end do
ccccccccccccccccccccc
          do im = 1, nmode(ielm)
          mypmder(im) = mypmder(im) * copisi(im)
          end do
ccccccccccccccccccccc


        if(.not.cokpco)then
          do im = 1, nmim
          kpar(im) = kphi * co * r0or
     ;               + si * rhoinv * ntni * (minf(ielm) - 1. + 0.5 * (im + 1))
c         Avoids k// exactly 0:
          if(dabs(kpar(im)).lt.kparze)kpar(im) = dsign(kparze, kpar(im))
          end do
	  else
c Ernesto:
c k// best defined for current surface in calling outpow
	  end if

      call dset(14*(2*maxpom-1), 0.d0, fa, 1)

        do im = 1, nmim
        mav2 = im - 1 + 2 * minf(ielm)
c       Here I group k and -k: is this ok with no up-down sym?
c        kstop1 = mod(mav2, 2)
c       Mod. 22/09/2003: fortran mod function result may be < 0!
        kstop1 = iabs(mod(mav2, 2))
          if(mav2 .le. minf(ielm)+msup(ielm))then
          kstop2 = min(klim, mav2-2*minf(ielm))
          else
          kstop2 = min(klim, 2*msup(ielm)-mav2)
          end if
        if(mod(mav2-kstop2,2).ne.0)kstop2 = kstop2 - 1
 
          do k = kstop1, kstop2, 2
          m1 = (mav2 + k) / 2 + 1 - minf(ielm)
          m2 = m1 - k
            if(k.eq.0)then
            tr = 1.d0
            else
            tr = 2.d0
            end if

c		All quantities summed over k (* = complex conjugate)

c		-> Re { E+(m1).E+(m2)* }
          fa(1,im) =  fa(1,im) + tr * drpc(pmder(7,m1), pmder(7,m2))

c		-> Re { E//(m1).E//(m2)* }
          fa(2,im) =  fa(2,im) + tr * drpc(pmder(10,m1), pmder(10,m2))

c		-> Im { E+(m1).(Nabla+.E//(m2))* } + Im { E+(m2).[Nabla+.E//(m1)]* }
          fa(3,im) =  fa(3,im) + tr *
     ;     ( dimag(pmder(7,m1) * dconjg(pmder(5,m2)))
     ;     + dimag(pmder(7,m2) * dconjg(pmder(5,m1))) )

c		-> Im { [Nabla+.E-(m1) - Nabla-.E+(m1)] . E//(m2)*} +
c		   Im { [Nabla+.E-(m2) - Nabla-.E+(m2)] . E//(m1)*}
          fa(4,im) =  fa(4,im) + tr *
     ;     ( dimag((pmder(3,m1)-pmder(2,m1)) * dconjg(pmder(10,m2)))
     ;     + dimag((pmder(3,m2)-pmder(2,m2)) * dconjg(pmder(10,m1))) )

c		-> Re { Nabla+.E+(m1).[Nabla+.E+(m2)]* }
          fa(5,im) =  fa(5,im) + tr * drpc(pmder(1,m1), pmder(1,m2))

c		-> Re { Nabla-.E+(m1).[Nabla-.E+(m2)]* }
          fa(6,im) =  fa(6,im) + tr * drpc(pmder(2,m1), pmder(2,m2))

c		-> Re { Nabla+.E-(m1).[Nabla+.E-(m2)]* }
          fa(7,im) =  fa(7,im) + tr * drpc(pmder(3,m1), pmder(3,m2))

c		-> Re { Nabla+.E//(m1).[Nabla+.E//(m2)]* }
          fa(9,im) =  fa(9,im) + tr * drpc(pmder(5,m1), pmder(5,m2))

c		-> Re { Nabla-.E//(m1).[Nabla-.E//(m2)]* }
          fa(10,im) = fa(10,im) + tr * drpc(pmder(6,m1), pmder(6,m2))

c		-> Re { [Nabla+.E-(m1)-Nabla-.E+(m1)].[Nabla+.E-(m2)-Nabla-.E+(m2)]* }
          fa(14,im) = fa(14,im) + tr * drpc(pmder(3,m1)-pmder(2,m1), 
     ;    pmder(3,m2) - pmder(2,m2))

c         Re (B//(m1).B//(m2)*)
          myfa(im)= drpc(mypmder(m1),mypmder(m2))
c         Re (B//(m1)E//(m1) . B//(m2)E//(m2)*)
          myfa2(im)= drpc(mypmder(m1)*pmder(10,m1),mypmder(m2)*pmder(10,m2))


          end do
        end do
      end if
c     ~~~~~~
 
      do 2 is = 1, ncl
c     ~~~~~~~~~~~~~~~~
      ispe = spli(is)
        if(insc)then
c	    (inside faraday shields)
      	vt(ispe) = vttab2(intab,ispe)
      	omp(ispe) = omptab2(intab,ispe)
      	else
c	    (outside)
      	vt(ispe) = vttab(intab,ispe)
      	omp(ispe) = omptab(intab,ispe)
      	end if

      omc = qom(ispe) * bmodul
 
      if(.not.cokpco)then
c     +++++++++++++++++++
        if(equidf(ispe).eq.'MAXWE')then
 
c       Maxwellian equilibrium:
c       ~~~~~~~~~~~~~~~~~~~~~~
 
          if(coldpl(ireg))then
          small = 0.d0
 
c         Cold plasma with collisional damping:
c         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
          coeff = damp(ispe,ireg) * 0.5 * eps0 * omp(ispe) ** 2
          dam2  = damp(ispe,ireg) ** 2
            if(dam2 .ne. 0)then
            pic10 = pic10 + coeff * (
     ;      dreal(zdotc(nmi2, pmder(7,1), 15, pmder(7,1), 15))
     ;      / ((omegag - omc) ** 2 + dam2 )
     ;    + dreal(zdotc(nmi2, pmder(15,1), 15, pmder(15,1), 15))
     ;      / ((omegag + omc) ** 2 + dam2 )
     ;      )
            plan = plan + coeff * (
     ;      dreal(zdotc(nmi2, pmder(10,1), 15, pmder(10,1), 15))
     ;      / (omegag ** 2 + dam2 )
     ;      )
            end if
 
c     if( small.ne.0.d0 .and. flrops(ireg) ) then
c     coef2 = eps0*( small * omp(ispe) )**2
c     sigod(1,1) = sigod(1,1) + coef2 *0.5d0 * (
c    ;   1.d0 / omegag
c    ; - 2.d0 /(omegag - omc + ci*damp(ispe,ireg))
c    ; + 1.d0 /(omegag - 2.d0*omc + ci*damp(ispe,ireg))   )
c
c     sigod(2,2) = sigod(2,2) + coef2 *0.5d0 * (
c    ;   1.d0 /(omegag + 2.d0*omc)
c    ; - 2.d0 /(omegag + omc)
c    ; + 1.d0 / omegag                         )
c     sigod(1,2) = sigod(1,2) + coef2 *0.25d0 * (
c    ;   1.d0 /(omegag + omc)
c    ; - 2.d0 / omegag
c    ; + 1.d0 /(omegag - omc + ci*damp(ispe,ireg))        )
c     sigod(3,3) = sigod(3,3) + coef2 *0.25d0 * (
c    ;   1.d0 /(omegag + omc)
c    ; - 2.d0 / omegag
c    ; + 1.d0 /(omegag - omc + ci*damp(ispe,ireg))        )
c     end if
 
          else
c         Warm plasma:
c         ~~~~~~~~~~~~
	pttmp2 = 0.0d0
	aux = 0
cERN		 cccccccccccccccccccccccccccc     dabs() here ??
          small = vt(ispe) / (dabs(omc) * rnorm)
          coeff = 0.5 * sqrtpi * eps0 * omp(ispe) ** 2
 
          do 6 im = 1, nmim
c           The following loop excludes the possibility of 'anomalous' ion cyclotron 
c           damping, and of electron cyclotron damping (both for imu<0):
c           To enable such effects, loop imu = -2 to 2 and export contributions:
            do imu = 0, 2
cERN		 cccccccccccccccccccccccccccc     dabs() here ??
            xmu(imu) = (omegag - dfloat(imu) * omc) / (kpar(im) * vt(ispe))
              if(abs(xmu(imu)) .le. 13.04)then
c             (13.04 = sqrt(170), threshold value in rif2 function)
              ex(imu) = dexp(- xmu(imu) ** 2) / dabs(kpar(im) * vt(ispe))
              else
              ex(imu) = 0.d0
              end if
            end do

c	    Order 0:
c	    -------
 
c			0.1) Fundamental ion cyclotron damping:
			pic10 = pic10 + coeff * ex(1) * fa(1,im)

c			0.2) Landau damping:
			plan = plan + coeff * 2.d0 * xmu(0) ** 2 * ex(0) * fa(2,im)
 
CC        IF(FLROPS(IREG))THEN
c         Finite Larmor radius terms
c         ~~~~~~~~~~~~~~~~~~~~~~~~~~
C         Evaluated in all 'not cold' situations;
C         NB: M=1 mode provides a sing. contribution (KETA) when fields
C             are evaluated to zeroth order in Larmor radius.
 
c         Order 1:
c	    -------
          coef1 = coeff * small

c			1.1) Fundamental ion cyclotron damping:
			pic11 = pic11 - 2.d0 * coef1 * xmu(1) * ex(1) * fa(3,im)

c			1.2) Mixed Landau-TTMP:
c			pttmp = pttmp + dsqrt(2.d0) * coef1 * xmu(0) * ex(0) * fa(4,im)
			pttmp = pttmp + dsqrt(2.d0) * coef1 * xmu(0) * ex(0) * myfa2(im) 


c         Order 2:
c		-------	
          coef2 = coef1 * small

c			2.1) First harmonic ion cyclotron damping (omega=2 omc):
c			pic2  = pic2 + coef2 * ex(2) * fa(5,im)
			pic2  = pic2 + 2.d0*coef2 * ex(2) * fa(5,im)


cERN	------------- MISSING ---------------------
c			2.2) Fundamental ion cyclotron damping:
CC			PIC12 = PIC12+COEF2*EX(1)*(
CC   ;		XMU(1)**2*FA(9,IM)
CC   ;		)
cc			A TERMINER!!!
CC   ;		-FA(5,IM)-FA(6,IM)
CC   ;		+2.*FA(8,IM)+2.*FA(11,IM)-FA(12,IM) )

c	-------------------------------------------

cERN 05/03/05: Checking FLR2 terms

c			2.3) TTMP:
cERN			pttmp = pttmp + coef2 * ex(0) * fa(14,im)
cERN			New variable pttmp2 (COMPOW.copy) for 2nd order FLR TTMP
c			pttmp2 = pttmp2 + coef2 * ex(0) * fa(14,im)
			pttmp2 = pttmp2 + coef2 * ex(0) * myfa(im)


c			2.4) Landau damping:
			pel2 = pel2 - coef2 * ex(0) * xmu(0) ** 2 * (fa(9,im) + fa(10,im))
 




CC    ELSE
CC    END IF
 
cERN
cERN	

   6      continue
          end if
	
	


	

      else
c     Nonmaxwellian equilibria (Raymond's routines):
c     ~~~~~~~~~~~~~~~~~~~~~~~~~
      call newteni( equidf(ispe),
     ; omc, omp(ispe), omegag,
     ; akperp, -kphi,
     ; vt(1), vt(ispe), vt(ispe), v0alph, udrift, ks )
 
      call mucrvz(3, 3, ks, 3, 3, efield, 1, 3, ytra)
c     'Slowing down' power:
      psd = psd + 0.5 * omegag * eps0 * dimag(zdotc(3, efield, 1, ytra, 1))
 
      end if
c     ~~~~~~
      else
c     ++++
c     Constant k// coordinates
c Ernesto:
      end if
c     ++++++
 
   2  continue
c     ~~~~~~~~

      return
 
      end
