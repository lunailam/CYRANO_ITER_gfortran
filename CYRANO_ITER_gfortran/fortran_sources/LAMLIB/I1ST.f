      INTEGER FUNCTION I1ST(ISUB)                                       

      IMPLICIT NONE
      INTEGER ISUB

c     GIVES INDEX OF LEFT ABSCISSA OF 1ST ELT.IN SUBREGION ISUB    
c     in the finite element mesh.        
                                                                        

      include 'pardim.copy'
      include 'comfin.copy'
      include 'comsub.copy'

      INTEGER N, I

      N = 0                                                               
        DO I = 1, ISUB - 1                                                   
        N = N + IELE(I) 
        END DO 
                                                     
      I1ST = N                                                            

      RETURN                                                            
      END                                                               
