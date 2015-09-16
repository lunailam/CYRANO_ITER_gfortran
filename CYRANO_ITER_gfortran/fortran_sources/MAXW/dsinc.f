      double precision function dsinc(x)

      double precision x

	dsinc = 1.d0

	if(dabs(x) .gt. 1.d-16)dsinc = dsin(x) / x

	return
	end
