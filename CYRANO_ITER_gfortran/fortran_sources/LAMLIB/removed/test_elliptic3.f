      ! *************************************************************** !
      !                                                                 !
      ! SUBROUTINE test_elliptic3(cx1, cx2)				      !
      !                                                                 !
      ! Routine for testing the various functions that compute          !
      ! the Elliptic Integral of 3rd kind with complex parameter        !
      ! alpha = cx1 + i*cx2 and modulus 0<k<1 (in Carlson.f)            !
      !                                                                 !                       !
      ! (1) cdelpc(-alpha,k)     : Byrd & Friedman                      !
	! (2) cdellpic (-alpha,k)  : Carlson                              !  
      ! (3) elit3 (+alpha,k)     : S. Zhang & J. Jin                    !
      !                                                                 ! 
      ! Output to file: ...../M12run/TESTPI.dat                         !     
      ! Run MATLAB routine testpi.m to see results                      !
      !                                                                 !                                                      
      ! *************************************************************** !
                                               
      SUBROUTINE test_elliptic3(cx1, cx2)

      implicit none

      ! Input
      real*8, intent(in) :: cx1, cx2

      integer    :: j, OpenStat
	character(100) :: FILE_NAME
      real*8     :: kaux(100)
	complex*16 :: TESTPI(100), TESTPI2(100), TESTPI3(100) 
      complex*16 :: cdelpc, cdellpic

      external cdelpc, cdellpic 

cc    complex alpha parameter (cx1 + i*cx2)
c      cx1 = -1.8d0
c	 cx2 = -1.6d0

      do j =1,100
         kaux(j) = dfloat(j-1)/100.0d0
         TESTPI(j) = cdelpc(-dcmplx(cx1,cx2), kaux(j))
         TESTPI2(j) = cdellpic(-dcmplx(cx1,cx2), kaux(j))
         call elit3(+dcmplx(cx1,cx2), kaux(j), TESTPI3(j))
	end do

         FILE_NAME = "../../M12run/TESTPI.dat" 
	   open (UNIT = 111, FILE = FILE_NAME, STATUS = "REPLACE",
     ;         IOSTAT = OpenStat, ACTION = "WRITE")
	         if (OpenStat > 0) then
	             print *, 'Error writing file: ', FILE_NAME
	             stop
	         end if
               write(111,"(7G16.8)"), cx1, cx2, '0.0d0', '0.0d0', '0.0d0', '0.0d0', '0.0d0'
               do j =1,100
                  write(111,"(7G16.8)"), kaux(j) 
     ;                                 , dreal(TESTPI(j)), dimag(TESTPI(j))
     ;					   , dreal(TESTPI2(j)), dimag(TESTPI2(j))
     ;                                 , dreal(TESTPI3(j)), dimag(TESTPI3(j))
                  print *, kaux(j), dreal(TESTPI(j)), dimag(TESTPI(j))
               end do
          close (111)

	                                                              
      END SUBROUTINE test_elliptic3                                     
                                                                  
! *************************************************************** !

