      subroutine interp_QLFP_noout(input_rad)
c      subroutine intf0(icase)
      
      implicit none

c      integer icase
      
      include 'pardim.copy'
      include 'comfic.copy'
      include 'comfin.copy'
      include 'commag.copy'
      include 'compla.copy'
      include 'comgeo.copy'
      include 'comswe.copy'
      include 'comgdr.copy'
      include 'comin2.copy'
      include 'comphy.copy'
      include 'cokpco.copy' ! ERN
      include 'compow.copy' ! ERN
      
c     Interpolates the the QLFP distribution given at the
c	raqlfp surfaces to current CYRANO radius abscis(intab)
c     xn mesh: use QLFP grid 
c     v mesh:  use QLFP mesh (interpolated in read_QLFPheader.f)       
     
      integer, intent(in) :: input_rad

      integer irql, kr, jv, ix, isp, i, j, i1, iv1, ix1, j1, jf, ja, jo, ipo
     ;, irv, irx, jev, iex, L, i2, igx, igv, ispa, isrchfge

      double precision s(2), dri, ksi, ag, vel, s1, s2, t1, t2, FFAC, vo, Abeam
     ;                , alp, bet
     
	character(100) :: FILE_NAME
	integer :: OpenStat,k


	isp = 1
	ispa = ispgdr(isp)

        intab = input_rad


c	Find closest QLFP radius (irql > abscis > irql-1)
      irql = isrchfge(nraddr, raqlfp, 1, abscis(intab))
	if(irql.eq.1)irql=2
	if(irql.gt.nraddr)irql=nraddr


c       Linear interpolation of fo, dfo/dv, dfo/dx to current CYRANO radius 
c	  - original QLFP data: F0  (:,:,radqlfp,:,:)
c	  - interpolated  data: F0i (:,:,:,:) at rho=abscis(intab)

        i = irql
        i1 = irql - 1

        dri = 1.d0 / (raqlfp(i) - raqlfp(i1))
        alp = (raqlfp(i) - abscis(intab)) * dri
        bet = (abscis(intab) - raqlfp(i1)) * dri
	
		F0i(1:nabsx,1:nabsv,1:3,isp) = alp * F0(1:nabsx,1:nabsv,i1,1:3,isp)  
     ;                                 + bet * F0(1:nabsx,1:nabsv,i ,1:3,isp)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	  Normalize fo to CYRANO definition (QLFP normalization: fo=1 at v=0)
cERN	  F0i(1:nabsx,1:nabsv,1:3,isp) =  
cERN     ;  F0i(1:nabsx,1:nabsv,1:3,isp) * 1.d20*dentab(intab,ispa) / (sqrtpi * vttab(intab,ispa)) ** 3
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c	x-mesh (no interp, use QLFP values at irql)
	lxg = nabsx
	xngaug(1:lxg)= xgta(1:nabsx,i)

cERN  Avoid xn=+-1 cccccccccccccccccccccccccccccccccccccccccc
	xngaug(1)  = -0.999d0
	xngaug(lxg)=  0.999d0
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c	v-mesh (no interp, use values at irql) --> see read_QLFPheader
	lvg = nabsv
	vgaug(1:lvg)= vgta(1:nabsv,i)




      return

2222  format(200(g18.6))   

      end
