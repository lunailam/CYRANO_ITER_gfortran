      subroutine interp_QLFP
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
     
      integer irql, kr, jv, ix, isp, i, j, i1, iv1, ix1, j1, jf, ja, jo, ipo
     ;, irv, irx, jev, iex, L, i2, igx, igv, ispa, isrchfge

      double precision s(2), dri, ksi, ag, vel, s1, s2, t1, t2, FFAC, vo, Abeam
     ;                , alp, bet
     
	character(100) :: FILE_NAME
	integer :: OpenStat,k


	isp = 1
	ispa = ispgdr(isp)


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




cERN	Output added
	if (WRITE_OUTPUT) then

          if(cokpco)then
             FILE_NAME = trim(COKFOLDER)  // "/VXgrid.dat"
          else
             FILE_NAME = trim(STDFOLDER)  // "/VXgrid.dat"
	  end if

			open (UNIT = 66, FILE = FILE_NAME, STATUS = "REPLACE",
c     ;			  BUFFERED='YES', BLOCKSIZE = 32*512, BUFFERCOUNT = 50, 
     ;              IOSTAT = OpenStat, ACTION = "WRITE")
	              if (OpenStat > 0) then
		              print *, 'Error writing file: ', FILE_NAME
	                  stop
				  end if
	  
				  write(66,*),"(v,xn) grid for Non-Maxwellians (interp_QLFP.f)"
c					  write(66,*),lvg,lxg
				  write(66,*) ' radius ', ' density(1/m3) ', ' vT (m/s) '
			  write(66,*) abscis(intab), 1.d20*dentab(intab,ispgdr(1)), vttab(intab,ispgdr(1))

				write(66,*)'v'				
				write(66,*),lvg
				do k = 1,lvg
				  write(66,*) vgaug(k) 
				end do

				write(66,*)'xn'				
				write(66,*),lxg
				do k = 1,lxg
				  write(66,*) xngaug(k) 
				end do

				write(66,*)'F0(v:lines,xn:cols)'				
c				write(66,*),lvg, lxg
				do k = 1,lvg
				  write(66,2222) vgaug(k), F0i(1:lxg,k,1,1)
				end do

				write(66,*)'dF0/dv(v:lines,xn:cols)'				
c				write(66,*),lvg, lxg
				do k = 1,lvg
				  write(66,2222) vgaug(k), F0i(1:lxg,k,2,1)
				end do

	close (66)

	end if





      return

2222  format(200(g18.6))   

      end
