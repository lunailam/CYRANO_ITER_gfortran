      SUBROUTINE LONG(N, X, XL, I)                                         

      IMPLICIT NONE
      INTEGER N, I
      DOUBLE PRECISION X(0:N), XL(N)                                     

C -CONSTRUIT LONGEURS INTERVALLES A PARTIR DES ABSCISSES-               
C  (I=1),OU INVERSEMENT(I=2)                                            

      INTEGER J
      
      GO TO (11,22)I 
                                                         
  11  DO 1 J=1,N                                                        
      XL(J)=X(J)-X(J-1)                                                 
  1   CONTINUE                                                          
      RETURN                                                            

   22 DO 2 J=1,N                                                        
   2  X(J)=X(J-1)+XL(J)                                                 

      RETURN                                                            
      END                                                               
C                                                                       
C**********************************************************             
C                                                                       
