      subroutine regjum(iregb, side)

      implicit none
      character*1 side
      integer iregb

c     General jump constraints between 2 media are assembled for one
c     side of boundary ('L' or 'R').
c     Index of concerned region,subregion and element given in COMSWE.
c     Block-Thomas version with compact matrix storage as a vector.
 
      include 'pardim.copy'
      include 'comusr.f'
      include 'comreg.copy'
      include 'comsub.copy'
      include 'comfin.copy'
      include 'comswe.copy'
      include 'comequ.copy'
      include 'commod.copy'
      include 'comphy.copy'

      integer icolo, ibulo, icpblo, iell, ielb, ielr, nmodc, ncstrt
     ;, idum
      double precision cabs1, rownor(maxbll)
      complex*16 zdum
      cabs1(zdum) = dabs(dreal(zdum)) + dabs(dimag(zdum))
 
      if(iregb .gt. nreg) return
      icolo = iconn(isubr)
      ibulo = ibub(isubr)
      icpblo = icolo + ibulo

      write(603,*)'                    in regjum: ithoma=',ithoma,side

        if(side .eq. 'L' .or. side .eq. 'l')then
        iell = iel
        ielb = iel
        ielr = min0(iel+1,nele)
        nmodc = nmode(ielb)
        ncstrt = nmodc * ncstr(iregb)
c	  RHS: surface term in weak form and jump constraint
        call filrhs(.true.)
	 
        else if(side .eq. 'R' .or. side .eq. 'r')then
        if(iregb .eq. nreg) return
        iell = max0(1,iel-1)
        ielb = iel - 1
        ielr = iel
        nmodc = nmode(ielb)
        ncstrt = nmodc * ncstr(iregb)
        end if
c       write(603,*),'ael:', ael(lblock(ithoma+1)+6,1:10)
c       write(603,*),'btom:', btom(6,1:10)
c             write(603,*),'ael:', ael(42+6,1:10)
      call intbou(side, lblock(ithoma))

c      write(603,*),'ael(regjum):', ael(42+6,1:10)

c       print*, btom(6,1:30)
c        print*,'before fbstep'
c       print*, ael(306+6,1:10)
        if(gensys)then
c       Case routine is called during generation of linear system:
        call fbstep(0, rownor, rownor)
c        else if(.not.recalc_by_species_now)then
        else
c       Case routine is called after solution of linear system:
c       Compute block quadratic form, residual, ...:
        call qufetc2(.true., ithoma, .false., .true., .false., 1)
c        call qufetc2(.true., ithoma, .false., .true., .false., idum)
        end if
c        print*,'after fbstep'
c       print*, btom(6,1:30)

      ithoma = ithoma + 1
      write(603,*)'exit regjum'

      return
      end
