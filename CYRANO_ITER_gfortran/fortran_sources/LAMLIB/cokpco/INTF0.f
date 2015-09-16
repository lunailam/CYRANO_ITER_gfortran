      subroutine intf0
c      subroutine intf0(icase)
      
      implicit none

c      integer icase
      
      include 'PARDIM.COPY'
      include 'COMFIC.COPY'
      include 'COMFIN.COPY'
      include 'COMMAG.COPY'
      include 'COMPLA.COPY'
      include 'COMGEO.COPY'
      include 'COMSWE.COPY'
      include 'COMGDR.COPY'
      include 'COMIN2.COPY'
      include 'COMPHY.COPY' 
      
c     Interpolates the meshes and the equilibrium distribution radially.
c     xn mesh: impose separatrix as a region boundary; always uniform mesh
c     in each region. 
c     v mesh: interpolate region boundaries; impose uniform mesh on each region      
c     F0 and derivatives are computed at interpolated abscissae on 2 radial
c     meshes, using the finite element basis. Then a linear radial interpolation
c     is performed.

C      logical needrd
      
      integer irql, kr, jv, ix, isp, i, j, i1, iv1, ix1, j1, jf, ja, jo, ipo
     ;, irv, irx, jev, iex, L, i2, igx, igv
     ;, isrchfge

      double precision s(2), dri, ksi, ag
     
      irql = isrchfge(nraddr, raqlfp, 1, abscis(intab))

c     Local xn mesh:
c     Normalized x at separatrix (co-passing):
      xnsep = - 2.d0 * delb(intab) / bmax(intab)
      xnbl(0) = -1.d0
      xnbl(1) = xnsep
      xnbl(2) = 0.d0
      xnbl(3) = - xnsep
      xnbl(4) = 1.d0
        do i = 1, nxgreg
c       Region-wise list of element lengths (for uniform meshes):
        xlelr(i) = xnbl(i) - xnbl(i-1)
        call cuteq(xnbl(i-1), xnbl(i), ielrx(i), xngrid(ifiax(i)))
c       Powers of element length:
          do ipo = -1, 2
          xlpo(ipo,i) = xlelr(i) ** ipo
          end do
        end do

c  for an eventual x mesh, to declare...
c      xmax = B0 / bmin(intab)
c        do i = 1, nabsx
c        xgta(i,kr) = xmax * (1.d0 - dabs(xngta(i,kr)))
c        end do

c       Cases of 1 radial point or Cyrano abscissa outside of qlfp abscissae:
c       Use the available grid data at a different radius, on the same v mesh
c       but evaluated on a different xn mesh.
        if(nraddr.eq.1 .or. irql.eq.1 .or. irql.gt.nraddr)then
          if(nraddr.eq.1)then
          i1 = 1
          kr = 1
          ja = 1
          jo = 1
          s(1) = 1.d0
          s(2) = 0.d0
          else if(irql.eq.1)then
          i1 = 0
          kr = 1
          ja = 2
          jo = 2
          s(1) = 0.d0
          s(2) = 1.d0
          else if(irql.gt.nraddr)then
          i1 = nraddr
          kr = nraddr
          ja = 1
          jo = 1
          s(1) = 1.d0
          s(2) = 0.d0
          end if

c       Local v mesh:
        call dcopy(nvgreg, vrb(0,kr), 1, vrbl(0), 1)
        call dcopy(nabsv, vgta(1,kr), 1, vgrid, 1)
c       In this case the v nodes are identical; basis functions:
        call bafadn(1, 0.d0, bflnv(1,0,1,ja), bflvs(1,0,1,ja))
        call bafadn(1, 1.d0, bflnv(1,0,2,ja), bflvs(1,0,2,ja))
          do irv = 1, nvgreg
          call dcopy(maxbaf*3, bflnv(1,0,1,ja), 1, bflnv(1,0,ifiav(irv),ja), 1)
          ivt(ifiav(irv),ja) = ifiav(irv) + 1
            do jv = ifiav(irv)+1, ilaav(irv)
            ivt(jv,ja) = jv
            call dcopy(maxbaf*3, bflnv(1,0,2,ja), 1, bflnv(1,0,jv,ja), 1)
            end do
          end do

        else
        ja = 1
        jo = 2
c       Linear interpolation factors:
        i = irql
        i1 = i - 1
        dri = 1.d0 / (raqlfp(i) - raqlfp(i1))
        s(1) = (raqlfp(i) - abscis(intab)) * dri
        s(2) = (abscis(intab) - raqlfp(i1)) * dri
c       Local v mesh: piecewise uniform, with interpolated boundaries.
          vrbl(0) = 0.d0
          do j = 1, nvgreg
          vrbl(j) =  vrb(j,i1) * s(1) + vrb(j,i) * s(2)
          call cuteq(vrbl(j-1), vrbl(j), ielrv(1), vgrid(ifiav(j)))
          end do
c       Find corresp. nodes at radii rho1, rho2; compute basis functions in 
c       both grids:
c       v nodes loop:
        call dset(6*maxbaf*maxvg, 0.d0, bflnv, 1)
          do j = 1, 2
          j1 = i1 + j - 1
            do jv = 1, nabsv
            iv1 = isrchfge(nabsv, vgta(1,j1), 1, vgrid(jv))
            iv1 = max0(iv1, 2)
            ivt(jv,j) = iv1
              if(iv1.le.nabsv)then
              irv = irvoa(iv1)
              jev = iv1 - irv
              ksi = (vgrid(jv) - vgta(iv1-1,j1)) / vlelg(irv,j1)
              call bafadn(1, ksi, bflnv(1,0,jv,j), bflvs(1,0,jv,j))
              end if
            end do
          end do
        end if
c       Region-wise list of element lengths (for uniform meshes):
        do j = 1, nvgreg
        vlelr(j) = vrbl(j) - vrbl(j-1)
c       Powers of element length:
          do ipo = -1, 2
          vlpo(ipo,j) = vlelr(j) ** ipo
          end do
        end do

c     Find corresp. nodes at radii rho1, rho2; compute basis functions in 
c     both grids:
c     xn nodes loop:
      call dset(6*maxbaf*maxxg, 0.d0, bflnx, 1)
        do j = ja, jo
        j1 = i1 + j - 1
          do ix = 1, nabsx
          ix1 = isrchfge(nabsx, xngta(1,j1), 1, xngrid(ix))
          ix1 = max0(ix1, 2)
          ixt(ix,j) = ix1
            if(ix1.le.nabsx)then
            irx = irxoa(ix1)
            iex = ix1 - irx
            ksi = (xngrid(ix) - xngta(ix1-1,j1)) / xlelg(irx,j1)
            call bafadn(1, ksi, bflnx(1,0,ix,j), bflxs(1,0,ix,j))
            end if
          end do
        end do

      call dset(maxxg*maxvg*4*nspgdr, 0.d0, F0i, 1)
        do j = ja, jo
        j1 = i1 + j - 1
          do jv = 1, nabsv
          iv1 = ivt(jv,j)
          jev = iv1 - irvoa(jv)
            do ix = 1, nabsx
            ix1 = ixt(ix,j)
            iex = ix1 - irxoa(ix)
            call normat(vlelg(irvoa(jv),j1), xlelg(irxoa(ix),j1))
              do isp = 1, nspgdr
c             f0 and derivatives on interpolated grid, but at rho1 and rho2:
              call falp(jv, ix, i1, j, jev, iex, iv1, ix1, isp)
c               Linear interpolation:
                do jf = 1, 4
                F0i(ix,jv,jf,isp) = F0i(ix,jv,jf,isp) + s(j) * f0l(j1,jf)
                end do
              end do
            end do
          end do
        end do

c     Compute element lengths and grids in x and v:
      kr = 0
      L = 0
      i2 = 0
      if(gigdr)then
c     Use Gaussian abscissae:
        do i = 1, nxgreg
        L = L + 1
          do j = 1, ielrx(i)
          ag = xngrid(L)
          kr = kr + 1
          xlel(kr) = xngrid(L+1) - ag
            do igx = 1, ngaux
            i2 = i2 + 1
            xngaug(i2) = ag + xlel(kr) * agax(igx)
            end do
          L = L + 1
          end do
        end do
      lxg = ngaux * nelrx
      kr = 0
      L = 0
      i2 = 0
        do i = 1, nvgreg
        L = L + 1
          do j = 1, ielrv(i)
          ag = vgrid(L)
          kr = kr + 1
          vlel(kr) = vgrid(L+1) - vgrid(L)
            do igv = 1, ngauv
            i2 = i2 + 1
            vgaug(i2) = ag + vlel(kr) * agav(igv)
            end do
          L = L + 1
          end do
        end do
      lvg = ngauv * nelrv

      else
c     Use element nodes:
        do i = 1, nxgreg
        L = L + 1
          do j = 1, ielrx(i)
          ag = xngrid(L)
          kr = kr + 1
          xlel(kr) = xngrid(L+1) - ag
          L = L + 1
          end do
        end do
      lxg = nabsx
      call dcopy(nabsx, xngrid, 1, xngaug, 1)

      kr = 0
      L = 0
        do i = 1, nvgreg
        L = L + 1
          do j = 1, ielrv(i)
          ag = vgrid(L)
          kr = kr + 1
          vlel(kr) = vgrid(L+1) - vgrid(L)
          L = L + 1
          end do
        end do
      lvg = nabsv
      call dcopy(nabsv, vgrid, 1, vgaug, 1)

      end if

      return
      end