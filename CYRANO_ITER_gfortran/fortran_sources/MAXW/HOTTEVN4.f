      subroutine hottev(npoi, kparv, sigsd, sigod, onespe, isp, gdrtoo, onlyab)

      implicit none
      
      logical onespe, gdrtoo, onlyab
      
      integer npoi, isp
            
      double precision kparv(npoi)
      
c     19/3/92: ITAINC (in COMSWE) gives increment of index in poloidal tables.
c     
c     Computes +-// components of susceptibility tensor using
c     plasma dispersion function.
c     N.B.: this routine only deals with Maxwellian plasmas.
c     Vectorized version with non-axisymmetric N and T profiles.
c     Radial position specified with INTAB; NPOI poloidal points starting
c     from #1 with increments of ITAINC (given in COMSWE).
 
c     RESULTS:  3*3 MATRICES
c                SIGSD : NON DIFFERENTIAL PART OF TENSOR
c                SIGOD : COEFFICIENTS OF DIFFERENTIAL TERMS
c               (REFERS TO POLOIDAL MODE EXPANSION OPERATOR)
c 5/9/98: gdrtoo: if true, also compute the Maxwellian dielectric tensor for
c species involved in the general dielectric response. If false, omit these
c contributions.
c 28/1/2000: onlyab = .true. to compute damping only.
         
      include 'pardim.copy'
      integer nfc, ncalm
      parameter (nfc=npfft+1, ncalm=5*nfc*maxspe)
c      parameter (nfc=npfft/2+1, ncalm=5*nfc*maxspe)
      include 'comgeo.copy'
      include 'commag.copy'
      include 'comma2.copy'
      include 'commod.copy'
      include 'compla.copy'
      include 'comant.copy'
      include 'comswe.copy'
      include 'compri.copy'
      include 'comfic.copy'
      include 'comphy.copy'
      include 'comfin.copy'  ! ERN: Added for SWcoef
       
      logical allk0, inscr, insc

      integer i, itr(ncalm), ncl, isp1, isp2, itab, ipo, ipoi, ispe, imu
     ;, imu1, imu2, itast
     ;, spli(maxspe), is

      double precision tr(ncalm), tr2(ncalm), vtloc(maxspe,nfc)
     ;, omploc(maxspe,nfc)

      double precision 
     ;  coeff, coef2, small(maxspe,nfc)
     ;, kpar(5*nfc*maxspe), xmuv(5*nfc*maxspe), kparlo
     ;, xmu(-2:2), omc(maxspe,nfc), kpars

      complex*16 
     ;  sigsd(3,3,nfc), sigod(3,3,nfc)
     ;, zfc(-2:2), zt, zif4, zfcv(ncalm), omcol

      complex*16 :: sstix, krb  ! ERN: Added for SWcoef
      real*8 :: npar2  ! ERN: Added for SWcoef

      external inscr, zif4
 
      call zset(9*npoi,czero,sigsd,1)
      call zset(9*npoi,czero,sigod,1)
        if(flrops(ireg))then
        imu1 = -2
        imu2 = 2
        else
        imu1 = -1
        imu2 = 1
        end if
      itast = imu2 - imu1 + 1

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
 
      call dset(5*npoi*ncl,czero,xmuv,1)
 
      allk0 = .true.
      itab = 0
      ipo = 1

      do 3 ipoi = 1, npoi
c     ~~~~~~~~~~~~~~~~~~~
      phi = polang(ipo)
c     Requires phi defined:
      insc = inscr(ireg)
        if(insc)then
          do is = 1, ncl
          ispe = spli(is)
          vtloc(ispe,ipoi) = vttab2(intab,ispe)
          omploc(ispe,ipoi) = omptab2(intab,ispe)
          end do
        else
          do is = 1, ncl
          ispe = spli(is)
          vtloc(ispe,ipoi) = vttab(intab,ispe)
          omploc(ispe,ipoi) = omptab(intab,ispe)
          end do
        end if
 
c      kpars = mpk*si*rhoinv + kphi*co*r0orta(intab,ipo)
      kpars = kparv(ipoi)
      allk0 = allk0 .and. kpars.eq.0.d0

      do 1 is = 1, ncl
c     ~~~~~~~~~~~~~~~~
      ispe = spli(is)
      omc(ispe,ipoi) = qom(ispe) * bmotab(intab,ipo)
        if (coldpl(ireg)) then
        small(ispe,ipoi) = 0.d0
        else
        small(ispe,ipoi) = vtloc(ispe,ipoi) / (dabs(omc(ispe,ipoi)) * rnorm)
        end if
 
        do imu = imu1, imu2
        itab = itab + 1
        kpar(itab) = kpars
        if(kpars*small(ispe,ipoi) .ne. 0.d0)xmuv(itab) = 
     ;  (omegag - imu * omc(ispe,ipoi)) / (kpars * vtloc(ispe,ipoi))
        end do
  1   continue
      ipo = ipo + itainc
  3   continue

      if(.not.coldpl(ireg) .and. .not.allk0)
     ;call riftab(itab, xmuv, kpar, zfcv, tr, tr2, itr, onlyab)
 
      itab = imu1
      ipo = 1
      do 20 ipoi = 1, npoi

      do 2 is = 1, ncl
c     ~~~~~~~~~~~~~~~~
      ispe = spli(is)
      itab = itab + itast

      if(coldpl(ireg))then
c     ====================
c     Cold model with damping coeff:
      coeff = omploc(ispe,ipoi) ** 2 / omegag
      omcol = omegag + ci * damp(ispe,ireg)

        if(onlyab)then
c     Compute damping only:
      sigsd(1,1,ipoi) = sigsd(1,1,ipoi) - dimag(coeff/(omcol - omc(ispe,ipoi)))
      sigsd(2,2,ipoi) = sigsd(2,2,ipoi) - dimag(coeff/(omcol + omc(ispe,ipoi)))
      sigsd(3,3,ipoi) = sigsd(3,3,ipoi) - dimag(coeff/ omcol)
 
      if( small(ispe,ipoi).ne.0.d0 .and. flrops(ireg) ) then
      coef2 = ( small(ispe,ipoi) * omploc(ispe,ipoi) )**2 /  omegag

      sigod(1,1,ipoi) = sigod(1,1,ipoi)
     ; + coef2 *0.5d0 * dimag(
     ;   1.d0 / omcol
     ; - 2.d0 /(omcol - omc(ispe,ipoi))
     ; + 1.d0 /(omcol - 2.d0*omc(ispe,ipoi)) )
 
      sigod(2,2,ipoi) = sigod(2,2,ipoi) + coef2 *0.5d0 * dimag(
     ;   1.d0 /(omcol + 2.d0*omc(ispe,ipoi))
     ; - 2.d0 /(omcol + omc(ispe,ipoi))
     ; + 1.d0 / omcol                         )

      sigod(1,2,ipoi) = sigod(1,2,ipoi) + coef2 *0.25d0 * dimag(
     ;   1.d0 /(omcol + omc(ispe,ipoi))
     ; - 2.d0 / omcol
     ; + 1.d0 /(omcol - omc(ispe,ipoi)) )

      sigod(3,3,ipoi) = sigod(3,3,ipoi) + coef2 *0.25d0 * dimag(
     ;   1.d0 /(omcol + omc(ispe,ipoi))
     ; - 2.d0 / omcol
     ; + 1.d0 /(omcol - omc(ispe,ipoi)) )
      end if

        else 
c     Active and reactive terms
      sigsd(1,1,ipoi) = sigsd(1,1,ipoi) - coeff / (omcol - omc(ispe,ipoi))
      sigsd(2,2,ipoi) = sigsd(2,2,ipoi) - coeff / (omcol + omc(ispe,ipoi))
      sigsd(3,3,ipoi) = sigsd(3,3,ipoi) - coeff /  omcol
 
      if( small(ispe,ipoi).ne.0.d0 .and. flrops(ireg) ) then
      coef2 = ( small(ispe,ipoi) * omploc(ispe,ipoi) )**2 /  omegag

      sigod(1,1,ipoi) = sigod(1,1,ipoi)
     ; + coef2 *0.5d0 * (
     ;   1.d0 / omcol
     ; - 2.d0 /(omcol - omc(ispe,ipoi))
     ; + 1.d0 /(omcol - 2.d0*omc(ispe,ipoi)) )
 
      sigod(2,2,ipoi) = sigod(2,2,ipoi) + coef2 *0.5d0 * (
     ;   1.d0 /(omcol + 2.d0*omc(ispe,ipoi))
     ; - 2.d0 /(omcol + omc(ispe,ipoi))
     ; + 1.d0 / omcol                         )

      sigod(1,2,ipoi) = sigod(1,2,ipoi) + coef2 *0.25d0 * (
     ;   1.d0 /(omcol + omc(ispe,ipoi))
     ; - 2.d0 / omcol
     ; + 1.d0 /(omcol - omc(ispe,ipoi)) )

      sigod(3,3,ipoi) = sigod(3,3,ipoi) + coef2 *0.25d0 * (
     ;   1.d0 /(omcol + omc(ispe,ipoi))
     ; - 2.d0 / omcol
     ; + 1.d0 /(omcol - omc(ispe,ipoi)) )
      end if
        end if
 
      else
c     ====
        do 4 imu = imu1, imu2
        i = imu + itab
        zfc(imu) = zfcv(i)
        xmu(imu) = xmuv(i)
  4     continue
      kparlo = kpar(itab)
      if(dabs(kparlo).lt.kparze)kparlo = dsign(kparze,kparlo)
      coeff = omploc(ispe,ipoi) ** 2 / (omegag * kparlo * vtloc(ispe,ipoi))
c    (singular at n=0 for m=0; or for all m if ipl=0)
 
      sigsd(1,1,ipoi) = sigsd(1,1,ipoi) + zfc(1) * coeff
      sigsd(2,2,ipoi) = sigsd(2,2,ipoi) + zfc(-1) * coeff
      sigsd(3,3,ipoi) = sigsd(3,3,ipoi) + 2.d0 * xmu(0) * 
     ;                                    (cun + xmu(0)*zfc(0)) * coeff
 
      if(flrops(ireg))then
      coef2 = coeff * small(ispe,ipoi)**2

      sigod(1,1,ipoi) = sigod(1,1,ipoi) - coef2 *0.5d0 * (
     ;   zfc(0)
     ; - 2.d0*zfc(1)
     ; + zfc(2)                                )

      sigod(2,2,ipoi) = sigod(2,2,ipoi) - coef2 *0.5d0 * (
     ;   zfc(-2)
     ; - 2.d0*zfc(-1)
     ; + zfc(0)                                )

      sigod(3,3,ipoi) = sigod(3,3,ipoi) - coef2 *0.5d0 * (
     ;          xmu(-1)**2 * zfc(-1)
     ; - 2.d0 * xmu( 0)**2 * zfc( 0)
     ; +        xmu( 1)**2 * zfc( 1)           )

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cERN  NEW: If species=1 (elec) then add IBW damping to gradE+.gradE+
c     CHECK sign of correction (+/- i)
      if(ispe.eq.1)then

c      sstix = 0.5 * (sigsd(2,2,ipoi) + cun + sigsd(1,1,ipoi) + cun)
c      npar2 = kparv(ipoi) / k0 * kparv(ipoi) / k0

c      krb =  k02 * 2. * (npar2-sstix)/( k0rn2 * coef2*0.5d0 * 
c     ;                              ( zfc(0)-2.d0*zfc(1)+zfc(2) ) )

c         if(iel.lt.1)iel=1
c         sigod(1,1,ipoi) = sigod(1,1,ipoi) + i * SWcoef*(fl(max0(iel,1)))**2
c         sigsd(1,1,ipoi) = sigsd(1,1,ipoi) + i * SWcoef * dreal(krb)**2 
          sigod(1,1,ipoi) = sigod(1,1,ipoi) + i * SWcoef 

      end if
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      sigod(1,2,ipoi) = sigod(1,2,ipoi) - coef2 *0.25d0 * (
     ; + zfc(-1)
     ; - 2.d0*zfc(0)
     ; + zfc(1)                                )

      zt = ci*coeff*small(ispe,ipoi)*sqrt2i

      sigod(1,3,ipoi) = sigod(1,3,ipoi) - zt * (
     ; - xmu(0)*zfc(0)
     ; + xmu(1)*zfc(1)               )

      sigod(2,3,ipoi) = sigod(2,3,ipoi) + zt * (
     ;   xmu(-1)*zfc(-1)
     ; - xmu( 0)*zfc( 0)             )
 
      end if
      end if
c     ======
 
   2  continue
c     ~~~~~~~~
      sigod(2,1,ipoi) = sigod(1,2,ipoi)
      sigod(3,1,ipoi) = sigod(1,3,ipoi)
      sigod(3,2,ipoi) = sigod(2,3,ipoi)
      ipo = ipo + itainc
  20  continue
c     ~~~~~~~~
 
      return
 
      end
