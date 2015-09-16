      complex*16 function cdsinc(z)

      complex*16 z

	cdsinc = (1.d0, 0.d0)

	if(z .ne. (0.d0,0.d0))cdsinc = cdsin(z) / z

	return
	end
