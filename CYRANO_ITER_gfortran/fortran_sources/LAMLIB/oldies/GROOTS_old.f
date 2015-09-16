      subroutine groots

      implicit none

C     Generates roots of the local dispersion at element nodes,
C     in equatorial plane, except at magnetic axis.
C     Generates equilibrium profiles in equatorial plane.
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
      
      integer npldis, mst1m, mst2m, mstsm, nrootp, iii, nkup, iplo, iphi, niphi
     ;, ipo, ipro, incplo, ie1, ie2, ipath, idro, idro2, idro3, iroot
     ;, nrom
     ;, izamax
     
      real second
      
      double precision omc, rlar
      
      complex*16 krho2(3,maxpom)
      
      external inscr, izamax, second
      
c     Max. number of roots at each point, for each mode:
      data nrom/3/
      
      write(nofile,*)'enter dispersion study ; time=',second()-timin
      call dset(9*maxpom*2*(maxnel+maxreg), 0.d0, roots, 1)
 
      npldis = nreg + nele
 
C     Fast wave kperp on inner boundary, outer equatorial plane:
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
      nrootp=2
      else
      nrootp=3
      end if

      write(6,*)'ompe=',omptab(1,1),' vte=',vttab(1,1)
c	Compute dispersion roots (krho^2)      
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
      do 1 iphi = 1, niphi
      ipo = 1 + npfft/2 * (iphi - 1)
      phi = polang(ipo)
      intabp = ipo
      incplo = 3 - 2 * iphi

      do 10 ireg = 1, nreg
        if(vacuum(ireg))then
        nrootp = 1
        else if(coldpl(ireg) .or. .not.flrops(ireg))then
        nrootp = 2
        else
        nrootp = 3
        end if
      ie1 = ifiel(ireg) - 1
c      ie1 = ifiel(ireg)
c      if(ireg.gt.1 .or. crown)ie1 = ie1 - 1
      ie2 = ilael(ireg)

      do 10 iel = ie1, ie2
      intab = iel * (ngauss + 1) + ireg
      rho = abscis(intab)
      y = abscno(intab)
      yinv = abscni(intab)
      rhoinv = yinv * rnori
      r0or = r0orta(intab,ipo)
      si = eqt(intab,ipo,15)
      co = eqt(intab,ipo,14)
      bmodul = bmotab(intab,ipo)
 
        if(cyl)then
        kphilo = kphi
        else
        kphilo = kphi * r0or
        end if

      if(nrootp.eq.3)then
c     2nd species Larmor radius:
      omc = qom(2) * bmodul
        if(inscr(ireg))then
        vt(2) = vttab2(intab,2)
        else
        vt(2) = vttab(intab,2)
        end if
      rlar = vt(2) / dabs(omc)
      end if
      
      iplo = iplo + incplo
      
c     To plot vs. R:
      absdis(iplo) = eqt(intab,ipo,1)
c      absdis(iplo) = incplo * rho
c     B0 profile:
      profil(iplo,2*nspec+1) = bmotab(intab,ipo)
      ipath = 1
      call disper(ipath, krho2, nrootp, .true.)
      idro3 = 2 * nrom * ((mstud2 - mstud1) / mstust + 1)
 
      do m = mstud1, mstud2, mstust
      mr = (m - mstud1) / mstust + 1
      idro = 2 * nrom * (mr - 1)
        do iroot = 1, nrootp
        roots(iplo,idro+2*iroot-1) = dreal(krho2(iroot,mr))
        roots(iplo,idro+2*iroot) = dimag(krho2(iroot,mr))
        end do
      end do

c     This is not always Bernstein root!!:
c     Krho max. * Larmor:
        if(nrootp.eq.3)then
        idro2 = 2 * mr - 1 + idro3
        roots(iplo,idro2) = rlar * dsqrt(cdabs(krho2(
     ;  izamax(nrootp, krho2(1,mr), 1), mr)))
c       Keta * Larmor:
        roots(iplo,idro2+1) = dabs(kper * rlar)
        end if

  10  continue
      iplo = npldis + 1
   1  continue
 
c     Equilibrium profiles:
      npldis = nreg + nele
      nkup = 0
        if(polsym)then
        iplo = 0
        niphi = 1
        else
        iplo = npldis
        niphi = 2
        end if

      do 11 iphi = 1, niphi
      ipo = 1 + npfft/2 * (iphi - 1)
      intabp = ipo
      incplo = 3 - 2 * iphi

      do 12 ireg = 1, nreg
      ie1 = ifiel(ireg) - 1
      ie2 = ilael(ireg)
        do iel = ie1, ie2
        intab = iel * (ngauss + 1) + ireg
        rho = abscis(intab)
        y = abscno(intab)
        iplo = iplo + incplo
c       Plot vs. ±rho:
        abspro(iplo) = incplo * rho
           
          if(inscr(ireg))then
            do ipro = 1, nspec
            profil(iplo,ipro) = dfac * dentab2(intab,ipro)
            profil(iplo,ipro+nspec) = temtab2(intab,ipro)
            end do
          else
            do ipro = 1, nspec
            profil(iplo,ipro) = dfac * dentab(intab,ipro)
            profil(iplo,ipro+nspec) = temtab(intab,ipro)
            end do
          end if
c       not appropriate: B should be vs R, not rho:
c       profil(iplo,2*nspec+1) = bmotab(intab,ipo)
        profil(iplo,2*nspec+2) = qtab(intab)
        end do
  12  continue
      iplo = npldis + 1
  11  continue
 
      write(nofile,*)'Dispersion study ended ; time=',second()-timin
      return
      end
