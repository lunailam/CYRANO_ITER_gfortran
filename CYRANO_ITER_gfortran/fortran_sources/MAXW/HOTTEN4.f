      subroutine hotten(kpar, sigsd, sigod, onespe, isp, gdrtoo, onlyab)

      implicit none
      integer isp
      logical  onespe, gdrtoo, onlyab
      double precision kpar
      complex*16  sigsd(3,3), sigod(3,3)

c     computes +-// components of maxwellian susceptibility tensor using
c     plasma dispersion function.
c     calls newteni routine for species described by R.Koch's nonmaxwellian
c     dielectric tensor.
c     arguments:
c     input:  kpar: component of wavevector along b0
c             onespe: if true, compute tensor for species isp only
c             gdrtoo: if false, only compute tensor for maxwellian species (i.e.
c                     not listed as having a general dielectric response;
c                     if true, compute maxwellian tensor irrespective of the
c                     type of the distribution function.
c             onlyab: if true, only compute active terms (damping calculation)
c                     if false (normal use), compute active and reactive.
c     input through commons: intab, intabp
c     output:  3*3 matrices
c              sigsd : non differential part of tensor
c              sigod : coefficients of differential terms
c             (REFERS TO POLOIDAL MODE EXPANSION OPERATOR)

      include 'pardim.copy'
      include 'comgeo.copy'
      include 'commag.copy'
      include 'compla.copy'
      include 'commod.copy'
      include 'comswe.copy'
      include 'comfic.copy'
      include 'comphy.copy'
      include 'comfin.copy'  ! ERN: Added for SWcoef

      logical  inscr, insc

      integer iwr

      integer ispe, imu, imu1, imu2
     ;, spli(maxspe), is, ncl, i

      double precision coeff, coef2, xmu(-2:2), omc, kpars, rlar, small

      complex*16  zfc(-2:2), zt, zif4, omcol, ks(3,3)
      
      complex*16 :: sstix, krb  ! ERN: Added for SWcoef
      real*8 :: npar2  ! ERN: Added for SWcoef


      external inscr, zif4

      iwr = nofile
      call zset(9, czero, sigsd, 1)
      call zset(9, czero, sigod, 1)
        if(flrops(ireg))then
        imu1 = -2
        imu2 = 2
        else
        imu1 = -1
        imu2 = 1
        end if

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

      bmodul = bmotab(intab, intabp)
c     phi must be defined to find whether we are inside a screen box:
      phi = polang(intabp)
      insc = inscr(ireg)

      do 2 is = 1, ncl
c     ~~~~~~~~~~~~~~~~
      ispe = spli(is)
        if(insc)then
      	vt(ispe) = vttab2(intab,ispe)
      	omp(ispe) = omptab2(intab,ispe)
      	else
      	vt(ispe) = vttab(intab,ispe)
      	omp(ispe) = omptab(intab,ispe)
      	end if
c      write(6,*)'vt(',ispe,')= ',vt(ispe)
      omc = qom(ispe) * bmodul
c      write(6,*)'omc= ',omc

      if(.not.usenma .or. equidf(ispe).eq.'MAXWE')then
c     ------------------------------------------------
      	if(coldpl(ireg))then
      	small = 0.d0
      	rlar = 0.d0
      	else
      	small = vt(ispe) / (dabs(omc) * rnorm)
      	if(ispe.eq.ispres)rlar = small * rnorm
      	end if

      if(coldpl(ireg))then
c     limit smallness of resonance width with damping coeff:
cc    if ( coldpl(ireg) .or.
cc   ; abs(kpar*vt(ispe)) .le. 1.08*damp(ispe,ireg) ) then
cc    if ( kpar*small .eq. 0.d0 ) then
c     ================================
      coeff = omp(ispe) ** 2 / omegag
      omcol = omegag + ci * damp(ispe,ireg)

        if(onlyab)then
c     Absorption only:
      sigsd(1,1) = sigsd(1,1) - dimag(coeff / (omcol - omc))
      sigsd(2,2) = sigsd(2,2) - dimag(coeff / (omcol + omc))
      sigsd(3,3) = sigsd(3,3) - dimag(coeff /  omcol)

      	if( small.ne.0.d0 .and. flrops(ireg) ) then
      	coef2 = (small * omp(ispe)) ** 2 /  omegag
      	sigod(1,1) = sigod(1,1) + coef2 * 0.5d0 * dimag(
     ;	   1.d0 / omcol
     ;	 - 2.d0 /(omcol - omc)
     ;	 + 1.d0 /(omcol - 2.d0*omc)                    )

      	sigod(2,2) = sigod(2,2) + coef2 * 0.5d0 * dimag(
     ;	   1.d0 /(omcol + 2.d0*omc)
     ;	 - 2.d0 /(omcol + omc)
     ;	 + 1.d0 / omcol                                )
      	sigod(1,2) = sigod(1,2) + coef2 * 0.25d0 * dimag(
     ;	   1.d0 /(omcol + omc)
     ;	 - 2.d0 / omcol
     ;	 + 1.d0 /(omcol - omc)                          )
      	sigod(3,3) = sigod(3,3) + coef2 * 0.25d0 * dimag(
     ;	   1.d0 /(omcol + omc)
     ;	 - 2.d0 / omcol
     ;	 + 1.d0 /(omcol - omc)                          )
      	end if
        else
c     Full tensor:
      sigsd(1,1) = sigsd(1,1) - coeff / (omcol - omc)
      sigsd(2,2) = sigsd(2,2) - coeff / (omcol + omc)
      sigsd(3,3) = sigsd(3,3) - coeff /  omcol

      	if( small.ne.0.d0 .and. flrops(ireg) ) then
      	coef2 = (small * omp(ispe)) ** 2 /  omegag
      	sigod(1,1) = sigod(1,1) + coef2 * 0.5d0 * (
     ;	   1.d0 / omcol
     ;	 - 2.d0 /(omcol - omc)
     ;	 + 1.d0 /(omcol - 2.d0*omc)               )

      	sigod(2,2) = sigod(2,2) + coef2 * 0.5d0 * (
     ;	   1.d0 /(omcol + 2.d0*omc)
     ;	 - 2.d0 /(omcol + omc)
     ;	 + 1.d0 / omcol                           )
      	sigod(1,2) = sigod(1,2) + coef2 * 0.25d0 * (
     ;	   1.d0 /(omcol + omc)
     ;	 - 2.d0 / omcol
     ;	 + 1.d0 /(omcol - omc)                     )
      	sigod(3,3) = sigod(3,3) + coef2 * 0.25d0 * (
     ;	   1.d0 /(omcol + omc)
     ;	 - 2.d0 / omcol
     ;	 + 1.d0 /(omcol - omc)                     )
      	end if
          end if

      else
c     ====
      if(dabs(kpar).lt.kparze)then
c      	if(kpar.eq.0.d0)then
      	kpars = dsign(kparze, kpar)
      	else
      	kpars = kpar
      	end if
      coeff = omp(ispe) ** 2 / (omegag * kpars * vt(ispe))
c     (singular at n=0 for m=0 or for all m if ipl=0)
c     if(wripoi)
c    ;write(iwr,400)ispe,vt(ispe),omp(ispe),omc,small,coeff

      	do imu = imu1, imu2
      	xmu(imu) = (omegag - float(imu) * omc) / (kpars * vt(ispe))
c       Calls to plasma dispersion function:
      	zfc(imu) = zif4(xmu(imu), kpars, onlyab)
c	    write(iwr,100)ispe,imu,xmu(imu),zfc(imu)
        end do
      
      sigsd(1,1) = sigsd(1,1) + zfc(1) * coeff
      sigsd(2,2) = sigsd(2,2) + zfc(-1) * coeff
      sigsd(3,3) = sigsd(3,3) + 2.d0 * xmu(0) * (cun + xmu(0)*zfc(0)) * coeff

      	if(flrops(ireg))then
      	coef2 = coeff * small**2
      	sigod(1,1) = sigod(1,1) - coef2 *0.5d0 * (
     ;	   zfc(0)
     ;	 - 2.d0*zfc(1)
     ;	 + zfc(2)                                )
      	sigod(2,2) = sigod(2,2) - coef2 *0.5d0 * (
     ;	   zfc(-2)
     ;	 - 2.d0*zfc(-1)
     ;	 + zfc(0)                                )
      	sigod(3,3) = sigod(3,3) - coef2 *0.5d0 * (
     ;	          xmu(-1)**2 * zfc(-1)
     ;	 - 2.d0 * xmu( 0)**2 * zfc( 0)
     ;	 +        xmu( 1)**2 * zfc( 1)           )

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cERN  NEW: If species=1 (elec) then add IBW damping to gradE+.gradE+
c     CHECK sign of correction (+/- i)
      if(ispe.eq.1)then

c      sstix = 0.5 * (sigsd(2,2) + cun + sigsd(1,1) + cun)
c      npar2 = kpar / k0 * kpar / k0

c      krb =  k02 * 2. * (npar2-sstix)/( k0rn2 * coef2*0.5d0 * 
c     ;                              ( zfc(0)-2.d0*zfc(1)+zfc(2) ) )

c	 if(iel.lt.1)iel=1
c        sigod(1,1) = sigod(1,1) + i * SWcoef*(fl(max0(iel,1)))**2
c	 sigsd(1,1) = sigsd(1,1) + i * SWcoef * dreal(krb)**2  
         sigod(1,1) = sigod(1,1) + i * SWcoef

      end if
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



      	sigod(1,2) = sigod(1,2) - coef2 *0.25d0 * (
     ;	 + zfc(-1)
     ;	 - 2.d0*zfc(0)
     ;	 + zfc(1)                                )
      	zt = ci * coeff * small * sqrt2i
      	sigod(1,3) = sigod(1,3) - zt * (
     ;	 - xmu(0)*zfc(0)
     ;	 + xmu(1)*zfc(1)               )
      	sigod(2,3) = sigod(2,3) + zt * (
     ;	   xmu(-1)*zfc(-1)
     ;	 - xmu( 0)*zfc( 0)             )
      	end if
	
      end if
c     ======

      else
c     ----
c     Nonmaxwellian species described by Raymond Koch's tensor 
c     (damping only - full tensor takes A LOT of time):
      call newteni( equidf(ispe), omc, omp(ispe), omegag,
     ;             akperp, -kphi,
     ;             vt(1), vt(ispe), vt(ispe),
     ;             v0alph, udrift, ks )
c     N.B.: this neglects the reactive contribution of the species!
      call zaxpy( 9, cun, ks, 1, sigsd, 1 )
      end if
c     ------

   2  continue
c     ~~~~~~~~
      sigod(2,1) = sigod(1,2)
      sigod(3,1) = sigod(1,3)
      sigod(3,2) = sigod(2,3)

      return

c  100   format(1h ,'resonant species n0',i4/1h ,1x,i4,2x,g23.15,2(2x,
c     ;   g23.15))
ccc  300   format(1h ,3(6(1x,g20.12)/1h )//)
c  400   format(1h ,'spec.',i4,1x,'vth=',e15.7,1x,'omp=',e15.7,1x,
c     ;    'omc=',e15.7,1x,'larm/rnorm=',e15.7,1x,'coef=',
c     ;      e15.7)

      end
