	MODULE vec_functions

	IMPLICIT NONE

	CONTAINS

c	-------------------------------------------------
c     1) Complex exponential of real*8 vector: exp(i.x)
      FUNCTION cdeid_vec(X)

      real*8, dimension(:), intent(in) :: X
	complex*16, dimension(size(X)) :: cdeid_vec

	cdeid_vec = dcmplx(dcos(X), dsin(X))

      END FUNCTION cdeid_vec
c	-------------------------------------------------






	END MODULE vec_functions
