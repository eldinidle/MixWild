! ************************************************
!                                                                   
!                   **********************                         
!                   *  SUBROUTINE  GEND  *                          
!                   **********************                          
!                                                                  
!    GENERATE A SCALAR MATRIX ( A SPECIAL CASE OF WHICH IS THE    
!    IDENTITY MATRIX), I.E., A SQUARE MATRIX ALL DIAGONAL        
!    ELEMENTS OF WHICH ARE EQUAL TO A SPECIFIED CONSTANT, WITH THE  
!    OFF-DIAGONAL ELEMENTS EQUAL TO ZERO                           
!                                                                 
!    CALL GEND (A,X,N,MS)                                        
!                                                               
!    A .......... OUTPUT MATRIX, N BY N                        
!    X .......... INPUT CONSTANT OR VARIABLE, SPECIFIES VALUE FOR   
!                 DIAGONAL ELEMENTS IN A                           
!    N .......... NUMBER OF ROWS IN A, EQUAL TO NUMBER OF COLUMNS 
!    MS ......... DESIRED STORAGE MODE OF A                      
!                                                               
! ************************************************
SUBROUTINE GEND(A,C,N,MS)                                     
   DOUBLE PRECISION A,C                                         
   DIMENSION A(1)                                              
   
   SELECT CASE(MS)
   CASE (0)        ! rectangular
      L = N*N                                                  
      DO J=1,L                                                   
         A(J) = 0                                                    
      END DO
      K = N + 1                                                  
      DO J=1,L,K                                              
         A(J) = C                                                 
      END DO
   
   CASE (1,3,4)    ! packed symm, upper or lower triangle
      L = (N*(N+1))/2                                        
      DO J=1,L                                            
         A(J) = 0                                             
      END DO
      K = 1                                               
      DO J=1,N                                         
         A(K) = C                                          
         K = K + J + 1                                    
      END DO
      
   CASE (2)        ! diagonal
      DO J=1,N                                     
         A(J) = C                                      
      END DO
   END SELECT

   RETURN                               
END SUBROUTINE GEND


! ************************************************
!                   **********************
!                   *  SUBROUTINE  GRAM  *                          
!                   **********************                         
!                                                                 
!    OBTAIN THE GRAMIAN MATRIX OF PRODUCTS OF ROW VECTORS OF A   
!    SPECIFIED MATRIX BY POSTMULTIPLYING THE MATRIX BY ITS TRANS-   
!    POSE                                                          
!                                                                 
!    CALL GRAM (A,C,M,N)                                         
!                                                               
!    A .......... INPUT MATRIX, M BY N, GENERAL RECTANGULAR (MSA=0) 
!    C .......... OUTPUT MATRIX, M BY M, GRAMIAN, SYMMETRIC (MSC=1)
!    M .......... NUMBER OF ROWS IN A                             
!    N .......... NUMBER OF COLUMNS IN A                         
!                                                               
! ************************************************
SUBROUTINE GRAM(A,C,M,N)                                      
   DOUBLE PRECISION A,C,X                                       
   DIMENSION A(M,N),C(1)                                       

   IC=0                                                      
   DO I=1,M                                             
      DO J=1,I                                            
         X=0.0D0                                                
         DO K=1,N                                          
            X=X+A(I,K)*A(J,K)                                    
         END DO
         IC=IC+1                                             
         C(IC)=X                                            
      END DO
   END DO
   
   RETURN                                        
END SUBROUTINE GRAM

! ************************************************
!                   **********************                         
!                   *  SUBROUTINE  LOCAT *                        
!                   **********************                       
!                                                               
!    LOCATE A SINGLE ELEMENT IN A MATRIX BY CALCULATING ITS ONE-   
!    DIMENSIONAL SUBSCRIPT FROM GIVEN ROW AND COLUMN SUBSCRIPTS   
!                                                                
!    CALL LOCAT (I,J,LSUB,M,MS)                                 
!                                                              
!    I .......... ROW SUBSCRIPT OF ELEMENT TO BE LOCATED      
!    J .......... COLUMN SUBSCRIPT OF ELEMENT TO BE LOCATED  
!    LSUB ....... ONE-DIMENSIONAL SUBSCRIPT OF ELEMENT (OUTPUT)    
!    M .......... NUMBER OF ROWS IN THE MATRIX THAT CONTAINS THE  
!                 ELEMENT                                        
!    MS ......... STORAGE MODE OF THIS MATRIX                   
!
!     Storage modes: 0 SQUARE, 1 PACKED SYMMETRIC, 2 DIAGONAL, 
!               3 PACKED LOWER TRIANGLE
!               4 PACKED UPPER TRIANGLE
!                                                              
! ************************************************
SUBROUTINE LOCAT(I,J,LSUB,M,MS)                              
   USE MIXLIB, ONLY: POST_ERROR
   
   IF(M .LT. I) THEN
      CALL POST_ERROR('ERROR MESSAGE FROM SUBROUTINE LOCAT: INDEX OUT OF RANGE')
      STOP
   ELSE
   
      SELECT CASE (MS)
      CASE (0)   ! Square
         LSUB=(J-1)*M+I                                           
         
      CASE (1)   ! Packed symmetric
         IF(I.LT.J) THEN
            LSUB=(I*(I-1))/2+J                                    
         ELSE
            LSUB = (J*(J-1))/2 + I                                         
         ENDIF
         
      CASE (2)   ! Diagonal
         IF(I .EQ. J) THEN
            LSUB=I                                             
         ELSE
            CALL POST_ERROR('ERROR MESSAGE FROM SUBROUTINE LOCAT: INDEX OUT OF RANGE')
            STOP
         ENDIF
         
      CASE (3)   ! Packed Lower Triangle
         IF(I.LT.J) THEN
            CALL POST_ERROR('ERROR MESSAGE FROM SUBROUTINE LOCAT: INDEX OUT OF RANGE')
            STOP
         ELSE
            LSUB=(I*(I-1))/2+J                                    
         ENDIF
         
      CASE (4)   ! Packed Upper Triangle
         IF(I.LE.J) THEN
            LSUB = (J*(J-1))/2 + I                                         
         ELSE
            CALL POST_ERROR('ERROR MESSAGE FROM SUBROUTINE LOCAT: INDEX OUT OF RANGE')
            STOP
         ENDIF
      END SELECT
   ENDIF
   
   RETURN                                                        
END SUBROUTINE LOCAT

! ************************************************
!                   **********************                  
!                   *  SUBROUTINE   TRP  *                 
!                   **********************                
!                                                                    
!    TRANSPOSE A GENERAL RECTANGULAR MATRIX (MS=0)                  
!                                                                  
!    CALL TRP (A,B,MA,NA)                                         
!                                                                
!    A .......... INPUT MATRIX, M BY N, GENERAL RECTANGULAR (MS=0)   
!    B .......... OUTPUT MATRIX, N BY M, TRANSPOSE OF A, GENERAL    
!                 RECTANGULAR (MS=0)                               
!    MA ......... NUMBER OF ROWS IN A                             
!    NA ......... NUMBER OF COLUMNS IN A                         
!                                                               
!                                                              
! ************************************************
SUBROUTINE TRP(A,B,MA,NA)                                    
   DOUBLE PRECISION A,B,C                                      
   DIMENSION B(NA,MA),A(MA,NA)                                
   
   IF(MA .EQ. NA) THEN
      B(1,1)=A(1,1)                                            
      IF(NA.NE.1) THEN
         DO J=2,NA                                            
            L = J - 1                                             
            B(J,J)=A(J,J)                                        
            DO I=1,L                                          
               C = A(I,J)                                         
               B(I,J) = A(J,I)                                   
               B(J,I) = C                                       
            END DO
         END DO
      ENDIF
   ELSE
      DO J=1,MA                                    
         DO I=1,NA                                   
            B(I,J) = A(J,I)                              
         END DO
      END DO
   ENDIF
   
   RETURN                                 
END SUBROUTINE TRP

! ************************************************
!                                                                   
!                   **********************                         
!                   *  SUBROUTINE  EXTC  *                        
!                   **********************                       
!                                                               
!    EXTRACT ONE COLUMN OF A MATRIX AND STORE IT IN A SEPARATE 
!    VECTOR                                                   
!                                                            
!    CALL EXTC (A,J,V,M,MS)                                 
!                                                          
!    A .......... INPUT MATRIX, M BY (J OR LARGER)        
!    J .......... INDEX OF COLUMN TO BE EXTRACTED        
!    V .......... OUTPUT VECTOR OF LENGTH M, EXTRACTED COLUMN       
!    M .......... NUMBER OF ROWS IN A                              
!    MS ......... STORAGE MODE OF A                               
!                                                                
! ************************************************
SUBROUTINE EXTC (A,J,B,M,MS)                                   
   USE MIXLIB, ONLY: POST_ERROR
   DOUBLE PRECISION A,B                                          
   DIMENSION A(1),B(1)                                          
   
   IF(MS.GT.0 .AND. M .LT. J) THEN
      CALL POST_ERROR('ERROR MESSAGE FROM SUBROUTINE EXTC: INDEX OUT OF RANGE')  
   ELSE
      SELECT CASE (MS)
      
      CASE (0)   ! Square
         K = M*(J-1)                                              
         DO I=1,M                                             
            K = K + 1                                              
            B(I) = A(K)                                           
         END DO
         
      CASE (3)   ! Packed Lower Triangle
         K = (J*(J+1))/2                                                   
         DO L=J,M                                                       
            B(L) = A(K)                                                     
            K = K + L                                                      
         END DO
         L=J-1                                                        
         IF(L.GT.0) THEN
            DO I=1,L                                                
               B(I)=0.0D0                                                
            END DO
         ENDIF
   
      CASE (2)   ! Diagonal
         DO I=1,M                                                   
            B(I) = 0.0D0                                                
         END DO
         B(J) = A(J)                                                
   
      CASE (1,4)   ! Symmetric or Packed Upper Triangle
         IF(MS .EQ. 1) THEN
            K = (J*(J+1))/2                                                   
            DO L=J,M                                                       
               B(L) = A(K)                                                     
               K = K + L                                                      
            END DO
         ENDIF
         
         K=(J*(J-1))/2                                           
         DO I=1,J                                             
            KZ = K + I                                            
            B(I) = A(KZ)                                         
         END DO
         
         IF(MS .EQ. 4) THEN
            L=J+1                                              
            IF(M.GE.L) THEN
               DO I=L,M                                      
                  B(I)=0.0D0                                      
               END DO
            ENDIF
         ENDIF
      END SELECT
   ENDIF
      
   RETURN                                    
END SUBROUTINE EXTC
