      complex*16 function cdeid(x)
      
      double precision x

c     Exponential of i times a double precision number
      
      cdeid = dcmplx(dcos(x), dsin(x))
      
      return
      end