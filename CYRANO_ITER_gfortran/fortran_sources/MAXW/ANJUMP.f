      subroutine anjump

      implicit none

c     Computes jumps at antennae for all mode dofs.
c     Result in JUANT(NDOF,NMODE(NA),NEFFAN), in COMANT

      include 'pardim.copy'
      include 'comgeo.copy'
      include 'compla.copy'
      include 'comfou.copy'
      include 'commod.copy'
      include 'comfic.copy'
      include 'comfin.copy'
      include 'comreg.copy'
      include 'comsub.copy'
      include 'comant.copy'
      include 'commag.copy'
      include 'comma2.copy'
      include 'comin2.copy'
      include 'comswe.copy'
      include 'comphy.copy'

      logical itis
c     ;, coljum

      integer i, j, ja, jm, imod, isc, iregne, nmodc
     ;, nsum
c     ;, ireg, iel

      complex*16 cucup, cucut, z4

      write(603,*)'enter anjump'

      do ja = 1, neffan + nscree * nmoscr
c     ===================================

        if(ja .le. neffan)then
        ireg = irbant(ja)
        else
        imod = mod(ja-neffan, nmoscr)
        isc = (ja-neffan-imod) / nmoscr
        if(imod.ne.0)isc = isc + 1
        ireg = irbscr(isc)
        end if

      isubr = nsum(ireg,ns)
c      intab = ireg + (ngauss + 1) * nsum(isubr,iele)
      intab = ifiabs(ifiel(ireg+1)) - 1
      iel = ilael(ireg)
      nmodc = nmode(iel)
        do i = 1, ncstr(ireg)
          do j = 1,nmodc
          juant(i,j,ja) = czero
          end do
        end do

      coljum = coldpl(ireg) .or. vacuum(ireg) .or. .not.flrops(ireg)
      iregne = ireg + 1
      if(iregne.ge.1 .and. iregne.le.nreg)coljum = coljum .and. (coldpl(iregne) .or. vacuum(iregne) .or. .not.flrops(iregne))
c     When equilibrium is not usual Maxwellian,
c     I assume here that a FLR expansion never takes place in other species -> cold-like boundary jumps:
      if(.not.glomax)coljum = .true.

      do i = 1, nmode(iel)
c     ====================
      m = i - 1 + minf(iel)

c     Find whether M belongs to MOANT:
      itis = .false.
        do j = 1, nmoant
          if(m.eq.moant(j))then
          itis = .true.
          jm = j
          end if
        end do
      mr = m + 1 - minf(iel)
      z4 = czero
      cucup = czero
      cucut = czero
      
        if(itis)then
        cucup = cura(jm,ja)
        cucut = rcurat(jm,ja)
c       Erho jump: charge mode / epsilon0
        if(coljum)z4 = chaa(jm,ja) / eps0
        end if

c     Jumps in
c     Erho: 
      juant(1,mr,ja) = z4
c     Etheta: always 0
      juant(2,mr,ja) = czero
c     Etheta': due to poloidal antenna current
      juant(3,mr,ja) = cucup
c     Ephi: always 0
      juant(4,mr,ja) = czero
c     Ephi': due to toroidal antenna current (this is currently 0!)
      juant(5,mr,ja) = cucut
c     write(6,*)'juant:',(juant(j2,mr,ja),j2=1,ncstr(ireg))

      end do
      end do
c     ======
      write(603,*)'Exit anjump'
      return
      end
