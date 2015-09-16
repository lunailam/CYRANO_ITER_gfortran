      subroutine plasmb(pla, plas, ipoint, k, v, iv1, ldv, switch, add)

c 3/4/2004: removed one argument before last (nfft), obsolete

      implicit none
      logical add
      integer ipoint, k, iv1, ldv, switch
c     ;, nfft
      complex*16 pla(6,6), plas(6,6), v(iv1:iv1+ldv-1,28)

c     Builds volume plasma term matrix PMA
c     from matrix of poloidal harmonics V.
c     Harmonics of current order K are stored in Kth row of V.
c     Circular plasma cross section (magnetic field angle independent
c     of poloidal angle).

cERN	04/03/05 : New terms added for FLR2 corrections in DSHAPE

      include 'pardim.copy'
      include 'comswe.copy'
      include 'comgeo.copy'
      include 'comsub.copy'
      include 'comma2.copy'
      include 'compla.copy'
      include 'comphy.copy'

      character*3 eletyp
      integer kp, i, j
      double precision tem
      complex*16 tem2, pma(6,6), pmas(6,6)

      call zset(36, czero, pma, 1)

      eletyp = styp(isubr)

c      kp = iabs(ipoint)
      kp = ipoint
c      if(ipoint.lt.0)kp = ipoint + nfft
c     2003 fix in polfft:
      if(ipoint.lt.0)kp = ipoint + npfft

      if(switch.eq.1)then
c     ~~~~~~~~~~~~~~~~~~~
c     For volume terms

      if(glomax)then
c     --------------
      pma(1,1) = - y * v(kp,1)
      pma(3,3) = - y * v(kp,2)
      pma(5,5) = - y * v(kp,3)

        if(.not.circ .and. eletyp.eq.'M23')then
c		 M23 in plasma: must still check correct ldif matrix in genera
		 pma(1,3) = - y * v(kp,4)
		 pma(1,5) = - y * v(kp,5)
		 pma(3,5) = - y * v(kp,6)
		 pma(3,1) = - pma(1,3)
		 pma(5,1) = - pma(1,5)
		 pma(5,3) = pma(3,5)
        end if

        if( .not.coldpl(ireg) .and. flrops(ireg) )then
c to see: this is for circular!
cERN	04/03/05  looking now ...
        call zset(36, czero, pmas, 1)
        pmas(1,1) = - y * v(kp,17)
        pmas(2,2) = - y * v(kp,18)
        pmas(3,3) = - y * v(kp,19)
        pmas(1,2) = - y * v(kp,14)
        pmas(1,3) = - y * v(kp,7)
        pmas(1,4) = - y * v(kp,10)
        pmas(1,5) = - y * v(kp,8)
        pma(1,6)  = - y * v(kp,11)
        pma(2,2)  = - y * v(kp,4)
        pma(2,4)  = - y * v(kp,9)
        pmas(3,4) = - y * v(kp,15)
        pmas(3,5) = - y * v(kp,12)
        pma(3,6)  = - y * v(kp,13)
        pma(4,4)  = - y * v(kp,5)
        pma(6,6)  = - y * v(kp,6)

        pmas(2,3) = - pmas(1,4)
        pmas(2,1) =   pmas(1,2)
        pmas(3,1) =   pmas(1,3)
        pmas(4,1) =   pmas(1,4)
        pmas(5,1) = - pmas(1,5)
        pma(6,1)  = - pma(1,6)
        pma(3,2)  =   pma(2,3)
        pma(4,2)  =   pma(2,4)
        pmas(4,3) =   pmas(3,4)
        pmas(5,3) = - pmas(3,5)
        pma(6,3)  = - pma(3,6)

cERN  04/03/05: New option for DSHAPE -----------------------------

          if(k.ne.0)then

			if (circ) then	! Standard (original) procedure

				tem = - k * co
				pmas(1,1) = pmas(1,1) + tem*(-v(kp,14))
				pma(1,2)  = pma(1,2) +  tem*(-v(kp,4) )
				pmas(1,3) = pmas(1,3) + tem*(-v(kp,10))
				pma(1,4)  = pma(1,4) +  tem*  v(kp,9)
				pmas(3,1) = pmas(3,1) + tem*(-v(kp,10))
				pma(3,2)  = pma(3,2) +  tem*(-v(kp,9) )
				pmas(3,3) = pmas(3,3) + tem*  v(kp,15)
				pma(3,4)  = pma(3,4) +  tem*  v(kp,5)
				pma(5,1)  = pma(5,1) +  tem*( v(kp,11))
				pma(5,3)  = pma(5,3) +  tem*(-v(kp,13))
				pmas(5,5) = pmas(5,5) + tem*  v(kp,16)

cERN	------------------------------------------------------
			else		! DSHAPE  (NEW coefficients: 20 ... 28)
				
				tem = - k 	! coefs. already multiplyed by cos(THETA)		
				pmas(1,1) = pmas(1,1) + tem*(-v(kp,26))		!  v14 * cos
				pma(1,2)  = pma(1,2) +  tem*(-v(kp,20) )	!  v4  * cos
				pmas(1,3) = pmas(1,3) + tem*(-v(kp,23))		!  v10 * cos
				pma(1,4)  = pma(1,4) +  tem*  v(kp,22)		!  v9  * cos
				pmas(3,1) = pmas(3,1) + tem*(-v(kp,23))		!  v10 * cos
				pma(3,2)  = pma(3,2) +  tem*(-v(kp,22) )	!  v9  * cos
				pmas(3,3) = pmas(3,3) + tem*  v(kp,27)		!  v15 * cos
				pma(3,4)  = pma(3,4) +  tem*  v(kp,21)		!  v5  * cos
				pma(5,1)  = pma(5,1) +  tem*( v(kp,24))		!  v11 * cos
				pma(5,3)  = pma(5,3) +  tem*(-v(kp,25))		!  v13 * cos
				pmas(5,5) = pmas(5,5) + tem*  v(kp,28)		!  v16 * cos

			end if
          end if
c	-------------------------------------------------------



c       (Real scalar times complex vector:)
        call zdscal(36, k0rn2, pmas, 1)
        end if

      else
c     ----
c to see in non circ! (Non-Maxwellian)
      pma(1,1) = - y * v(kp,1)
      pma(3,3) = - y * v(kp,2)
      pma(5,5) = - y * v(kp,3)
      pma(1,3) = - y * v(kp,4)
      pma(1,5) = - y * v(kp,5)
      pma(3,5) = - y * v(kp,6)
      pma(3,1) = pma(1,3)
      pma(5,1) = pma(1,5)
      pma(5,3) = pma(3,5)
      end if
c     ------
      call zdscal(36, k0rn2, pma, 1)

      else if(switch.eq.2)then
c     ~~~~~~~~~~~~~~~~~~~~~~~~
c     For boundary conditions
        if(coljum)then

          if(glomax)then
          pma(1,1) = v(kp,1)
          pma(3,3) = v(kp,2)
          pma(5,5) = v(kp,3)
            if(.not.circ .and. eletyp.eq.'M23')then
            pma(1,3) = v(kp,4)
            pma(1,5) = v(kp,5)
            pma(3,5) = v(kp,6)
            pma(3,1) = - pma(1,3)
            pma(5,1) = - pma(1,5)
            pma(5,3) = pma(3,5)
            end if
          else
          pma(1,1) = v(kp,1)
          pma(3,3) = v(kp,2)
          pma(5,5) = v(kp,3)
          pma(1,3) = v(kp,4)
          pma(1,5) = v(kp,5)
          pma(3,5) = v(kp,6)
          pma(3,1) = pma(1,3)
          pma(5,1) = pma(1,5)
          pma(5,3) = pma(3,5)
          end if

        else
c       Jump condition with flr contributions
cERN	  Is this OK for DSHAPE ?????????????????????????????????
        tem2 = ci * omegag * eps0 * rnorm * 0.5d0
        pma(1,1) = tem2 * v(kp,1)
        pma(1,2) = tem2 * v(kp,2)
        pma(1,3) = tem2 * v(kp,3)
        pma(1,4) = tem2 * v(kp,4)
        pma(3,1) = tem2 * (-v(kp,3))
        pma(3,2) = tem2 * v(kp,4)
        pma(3,3) = tem2 * v(kp,5)
        pma(3,4) = tem2 * v(kp,6)
        pma(5,1) = tem2 * v(kp,7)
        pma(5,3) = tem2 * v(kp,8)
        pma(5,6) = tem2 * v(kp,9)
        end if

      else if(switch.eq.3)then
c     ~~~~~~~~~~~~~~~~~~~~~~~~
c     Case SWITCH=1 *(-1/Y)

        if(glomax)then
        pma(1,1) = v(kp,1)
        pma(3,3) = v(kp,2)
        pma(5,5) = v(kp,3)
          if(.not.circ .and. eletyp.eq.'M23')then
          pma(1,3) = v(kp,4)
          pma(1,5) = v(kp,5)
          pma(3,5) = v(kp,6)
          pma(3,1) = - pma(1,3)
          pma(5,1) = - pma(1,5)
          pma(5,3) = pma(3,5)
          end if

          if( .not.coldpl(ireg) .and. flrops(ireg) )then
          call zset( 36, czero, pmas, 1 )
          pmas(1,1) = v(kp,17)
          pmas(2,2) = v(kp,18)
          pmas(3,3) = v(kp,19)
          pmas(1,2) = v(kp,14)
          pmas(1,3) = v(kp,7)
          pmas(1,4) = v(kp,10)
          pmas(1,5) = v(kp,8)
          pmas(1,5) = v(kp,8)
          pma(1,6)  = v(kp,11)
          pma(2,2)  = v(kp,4)
          pma(2,4)  = v(kp,9)
          pmas(3,4) = v(kp,15)
          pmas(3,5) = v(kp,12)
          pma(3,6)  = v(kp,13)
          pma(4,4)  = v(kp,5)
          pma(6,6)  = v(kp,6)
    
          pmas(2,3) = - pmas(1,4)
          pmas(2,1) =   pmas(1,2)
          pmas(3,1) =   pmas(1,3)
          pmas(4,1) =   pmas(1,4)
          pmas(5,1) = - pmas(1,5)
          pma(6,1)  = - pma(1,6)
          pmas(3,2) =   pmas(2,3)
          pma(4,2)  =   pma(2,4)
          pmas(4,3) =   pmas(3,4)
          pmas(5,3) = - pmas(3,5)
          pma(6,3)  = - pma(3,6)
    
            if(k.ne.0)then
            yinv = abscni(intab)
            tem =   k * co * yinv
            pmas(1,1) = pmas(1,1) + tem*(-v(kp,14))
            pmas(1,2) = pmas(1,2) + tem*(-v(kp,4) )
            pmas(1,3) = pmas(1,3) + tem*(-v(kp,10))
            pmas(1,4) = pmas(1,4) + tem*  v(kp,9)
            pmas(3,1) = pmas(3,1) + tem*(-v(kp,10))
            pmas(3,2) = pmas(3,2) + tem*(-v(kp,9) )
            pmas(3,3) = pmas(3,3) + tem*  v(kp,15)
            pmas(3,4) = pmas(3,4) + tem*  v(kp,5)
            pmas(5,1) = pmas(5,1) + tem*( v(kp,11))
            pmas(5,3) = pmas(5,3) + tem*(-v(kp,13))
            pmas(5,5) = pmas(5,5) + tem*  v(kp,16)
            end if
          call zdscal(36, k0rn2, pmas, 1)
          end if  ! FLR

        else  ! .not. glomax:
        pma(1,1) = v(kp,1)
        pma(3,3) = v(kp,2)
        pma(5,5) = v(kp,3)
        pma(1,3) = v(kp,4)
        pma(1,5) = v(kp,5)
        pma(3,5) = v(kp,6)
        pma(3,1) = pma(1,3)
        pma(5,1) = pma(1,5)
        pma(5,3) = pma(3,5)
        end if  ! glomax

      call zdscal(36, k0rn2, pma, 1)

      end if
c     ~~~~~~

      if(.not.add)then
      call zcopy(36, pma, 1, pla, 1)
      else 
        do i = 1, 6
          do j = 1, 6
          pla(i,j) = pla(i,j) + pma(i,j)
          end do
        end do
      end if

      if(.not.coldpl(ireg) .and. flrops(ireg))then
        if(.not.add)then
        call zcopy(36, pmas, 1, plas, 1)
        else 
          do i = 1, 6
            do j = 1, 6
            plas(i,j) = plas(i,j) + pmas(i,j)
            end do
          end do
        end if
      end if

      return
      end

