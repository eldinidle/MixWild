!*********************************************************************
!                   **********************                        
!                   *  SUBROUTINE PARTR  *                       
!                   **********************                      
!
!*********************************************************************

SUBROUTINE PARTR(A,B,C,MB,N,MC)
   DOUBLE PRECISION A,B,C
   DIMENSION A(1),B(1),C(1)
   
   L  = 0
   L1 = 0 
   L2 = 0
   
   DO J = 1,N    
   
      DO I = 1,MB   
         L1 = L1+1
         L  = I + (J-1)*(MB+MC)
         B(L1) = A(L)
      END DO
 
      DO I = 1,MC   
         L2 = L2+1
         L  = I + (J-1)*(MB+MC)
         C(L2) = A(L)
      END DO
      
   END DO

   RETURN  
END SUBROUTINE PARTR

!*********************************************************************
!                   **********************                        
!                   *  SUBROUTINE PRTRC  *                       
!                   **********************                      
!
!*********************************************************************
SUBROUTINE PRTRC(A,B,C,D,MC,NC)
   DOUBLE PRECISION A,B,C,D
   DIMENSION A(1),B(1),C(1),D(1)
   
   L = 0
   L1= 0
   L2= 0
   
   DO I = 1,MC+NC
      DO J = 1,I
         L = L+1
         NC2 = (NC*(NC+1))/2
         IF(L .LE. NC2) THEN
            B(L) = A(L)
         ELSE IF (L .GT. NC2 .AND. J .LE. NC) THEN
            L1   = L1+1
            LL   = (MC *(J-1)) + ((L1-1)/NC + 1)
            C(LL)= A(L)
         ELSE
            L2   = L2+1
            D(L2)= A(L)
         END IF
      END DO 
   END DO 

   RETURN  
END SUBROUTINE PRTRC
   
!*********************************************************************
!                                                                   
!                   **********************                         
!                   *  SUBROUTINE  TRAM  *                        
!                   **********************                       
!                                                               
!    OBTAIN THE TRACE (SUM OF DIAGONAL ELEMENTS) OF A SQUARE MATRIX 
!                                                                  
!    CALL TRAM (A,TR,N,MS)                                        
!                                                                
!    A .......... INPUT MATRIX, N BY N                          
!    TR ......... OUTPUT SCALAR, TRACE                         
!    N .......... NUMBER OF ROWS IN A, EQUAL TO NUMBER OF COLUMNS   
!    MS ......... STORAGE MODE OF A                                
!                                                                 
!*********************************************************************
SUBROUTINE TRAM(A,C,N,MS)                                       
   DOUBLE PRECISION A,C,SUM                                       
   DIMENSION A(1)                                                
   
   SUM=0.0D0                                                    
  
   SELECT CASE (MS)
   CASE (0)         ! Square
      J=N*N                                                     
      K=N+1                                                    
      DO I=1,J,K                                           
         SUM=SUM+A(I)                                           
      END DO
   CASE (1,3,4)     ! Packed symm, upper, or lower triangle
      J=0                                                  
      DO K=1,N                                         
         J=J+K                                              
         SUM=SUM+A(J)                                      
      END DO
   CASE (2)         ! diagonal 
      DO I=1,N                                     
         SUM=SUM+A(I)                                   
      END DO
   END SELECT
   C=SUM                                         
   
   RETURN                                   
END SUBROUTINE TRAM

!*********************************************************************
!                                       
!                   **********************                         
!                   *  SUBROUTINE ADJRC  *                        
!                   **********************                       
!                                                               
!    ADJOIN ONE RECTANGULAR MATRIX (MS=0) AND TWO TRIANGULAR   
!    (MS=3 OR 4) OR SYMMETRIC (MS=1) MATRICES TO OBTAIN A LARGE    
!    TRIANGULAR OR SYMMETRIC MATRIX                               
!                                                                
!    CALL ADJRC (A,B,C,D,MB,NB)                                 
!                                                              
!    A .......... INPUT MATRIX, NB BY NB, SYMMETRIC OR TRUE   
!                 TRIANGULAR (MSA=1,3,4)                     
!    B .......... INPUT MATRIX, MB BY NB, GENERAL RECTANGULAR      
!                 (MSB=0)                                         
!    C .......... INPUT MATRIX, MB BY MB, SYMMETRIC OR TRUE      
!                 TRIANGULAR (MSC=1,3,4)                        
!    D .......... OUTPUT MATRIX, (NB+MB) BY (NB+MB), SYMMETRIC OR  
!                 TRUE TRIANGULAR (MSD=1,3,4)                     
!    MB ......... NUMBER OF ROWS IN B, EQUAL TO NUMBER OF ROWS AND 
!                 NUMBER OF COLUMNS IN C                          
!    NB ......... NUMBER OF COLUMNS IN B, EQUAL TO NUMBER OF ROWS
!                 AND NUMBER OF COLUMNS IN A                    
!                                                              
!*********************************************************************
SUBROUTINE ADJRC (A,B,C,D,MB,NB)                             
   DOUBLE PRECISION A,B,C,D                                          
   DIMENSION A(1),B(MB,NB),C(1),D(1)                                
   
   K = 0                                                           
   L = (NB*(NB+1))/2                                              
   
   DO J=1,L                                            
      D(J) = A(J)                                         
   END DO
   
   DO I=1,MB                                        
      DO J=1,NB                                       
         L = L + 1                                        
         D(L) = B(I,J)                                   
      END DO
      
      DO J=1,I                                     
         L = L + 1                                     
         K = K + 1                                    
         D(L) = C(K)                                 
      END DO
   END DO
   
   RETURN                                  
END SUBROUTINE ADJRC

!*********************************************************************
!                                      
!                   **********************                          
!                   *  SUBROUTINE  INVT  *                         
!                   **********************                        
!                                                                
!    INVERT A TRUE TRIANGULAR MATRIX (MS=3,4) IN PLACE AND      
!    CALCULATE ITS DETERMINANT                                 
!                                                             
!    CALL INVT (A,N,DET)                                     
!                                                           
!    A .......... INPUT-OUTPUT MATRIX, N BY N, TRUE TRIANGULAR      
!                 (MS=3,4)                                         
!    N .......... NUMBER OF ROWS IN A, EQUAL TO NUMBER OF COLUMNS 
!    DET ........ OUTPUT SCALAR, DETERMINANT OF A                
!                                                               
!*********************************************************************
SUBROUTINE INVT(A,N,C)                                        
USE MIXLIB, ONLY: POST_ERROR
   DOUBLE PRECISION D,X,Y,Z,U,A,C                               
   DIMENSION A(1)                                              
   INTEGER DIAG,ROWNO,COLNO,ROWCOL                            
   
   D=1.0D0                                                   
   DIAG=0                                                   
   
   DO I=1, N
   
      DIAG=DIAG+I                                           
      X=A(DIAG)                                                          
      D=D*X                                                             
      IF(X .EQ. 0)  EXIT    ! escape out of  loop
      A(DIAG)=1.0D0/X                                                 
      COLNO=DIAG                                                     
      
      DO K=1,I-1
         COLNO=COLNO-1                                               
         Y=A(COLNO)                                                 
         Y=-Y/X                                                    
         A(COLNO)=Y                                               
         ROWNO=DIAG                                              
         ROWCOL=COLNO                                           
         DO J=I,N-1
            ROWNO=ROWNO+J                                       
            ROWCOL=ROWCOL+J                                    
            Z=A(ROWCOL)                                       
            U=A(ROWNO)                                       
            A(ROWCOL)=Z+Y*U                                 
         END DO    ! loop j = I,N-1
      END DO   ! loop k=1,i-1
      ROWNO=DIAG                                 
      
      DO J=I,N-1
         ROWNO=ROWNO+J                           
         Y=A(ROWNO)                             
         A(ROWNO)=Y/X                          
      END DO    ! loop j = I,N-1
      
   END DO  ! loop I=1,N
   
   C=D                              
   
   IF(C .EQ. 0) THEN
      CALL POST_ERROR('ERROR MESSAGE FROM SUBROUTINE INVT: MATRIX IS SINGULAR')
   ENDIF
   
   RETURN                                                      
END SUBROUTINE INVT

!*********************************************************************
!                                                                  
!                   **********************                        
!                   *  SUBROUTINE  INVD  *                       
!                   **********************                      
!                                                              
!    INVERT A DIAGONAL MATRIX (MS=2) IN PLACE AND CALCULATE ITS     
!    DETERMINANT                                                   
!                                                                 
!    CALL INVD (A,N,DET)                                         
!                                                               
!    A .......... INPUT-OUTPUT MATRIX, N BY N, DIAGONAL (MS=2) 
!    N .......... NUMBER OF ROWS IN A, EQUAL TO NUMBER OF COLUMNS   
!    DET ........ OUTPUT SCALAR, DETERMINANT OF A                  
!                                                                 
!*********************************************************************
SUBROUTINE INVD(A,N,C)                                          
USE MIXLIB, ONLY: POST_ERROR
   DOUBLE PRECISION D,X,A,C                                       
   DIMENSION A(N)                                                
   
   D=1.0D0                                                      
   
   DO I=1,N                                                 
      X=A(I)                                                     
      D=D*X                                                     
      IF (X .EQ. 0)  EXIT  ! leave loop
      A(I)  =1.0D0/X                                          
   END DO
   
   C=D                                                    
   IF(C .EQ. 0) THEN
      CALL POST_ERROR('ERROR MESSAGE FROM SUBROUTINE INVD: MATRIX IS SINGULAR')
   ENDIF
   
   RETURN                                                      
END SUBROUTINE INVD
   
!*********************************************************************
!                                                          
!                   **********************                
!                   *  SUBROUTINE  HMPY  *               
!                   **********************              
!                                                      
!    CALCULATE THE HADAMARD (TERM BY TERM) PRODUCT OF TWO MATRICES  
!                                                                  
!    CALL HMPY (A,B,C,M,N,MS)                                     
!                                                                   
!    A .......... INPUT MATRIX, M BY N                             
!    B .......... INPUT MATRIX, M BY N                            
!    C .......... OUTPUT MATRIX, M BY N, HADAMARD PRODUCT        
!    M .......... NUMBER OF ROWS IN A, EQUAL TO NUMBER OF ROWS IN B 
!    N .......... NUMBER OF COLUMNS IN A, EQUAL TO NUMBER OF       
!                 COLUMNS IN B                                    
!    MS ......... STORAGE MODE OF A, EQUAL TO STORAGE MODE OF B  
!                                                               
!*********************************************************************
SUBROUTINE HMPY(A,B,C,M,N,MS)                                 
   DOUBLE PRECISION A,B,C                                       
   DIMENSION A(1),B(1),C(1)                                    
   
   SELECT CASE (MS)
   CASE (0)
      MN=M*N                                                   
   CASE (1,3,4)
      MN=(M*(M+1))/2                                         
   CASE (2)
      MN=M                                                 
   END SELECT
   
   DO I=1,MN                                       
      C(I)=A(I)*B(I)                                     
   END DO
   RETURN                                        
END SUBROUTINE HMPY
   
!*********************************************************************
!                                            
!                   **********************  
!                   *  SUBROUTINE  CHOL  * 
!                   **********************
!                                        
!    COMPUTE THE CHOLESKY (SQUARE ROOT) FACTOR AND THE DETERMINANT  
!    OF A SYMMETRIC (MS=1) MATRIX                                  
!                                                                 
!    CALL CHOL (A,B,N,DET)                                       
!                                                               
!    A .......... INPUT MATRIX, N BY N, SYMMETRIC (MSA=1)      
!    B .......... OUTPUT MATRIX, N BY N, CHOLESKY FACTOR, TRUE
!                 LOWER TRIANGULAR (MSB=3)                   
!    N .......... NUMBER OF ROWS IN A, EQUAL TO NUMBER OF COLUMNS   
!    DET ........ OUTPUT VARIABLE CONTAINING DETERMINANT OF A      
!                                                                 
!*********************************************************************
SUBROUTINE CHOL(A,B,N,C)                                        
   DOUBLE PRECISION A,B,C,X,Y                                     
   DIMENSION A(1),B(1)                                           
   
   X=DSQRT(A(1))                                                
   C=X                                                         
   B(1)=X                                                     
   
   IF(N .GT. 1) THEN
      N1=N-1                                                   
      KC=1                                                    
      IFIR=1                                                 
      
      DO J=1,N1                                         
         KC=KC+J                                              
         B(KC)=A(KC)/X                                       
      END DO

      DO I=1,N1                                      
         IFIR=IFIR+I                                       
         KC=IFIR                                          
         X=0.0D0                                         
         
         DO J=1,I                                   
            X=X+B(KC)*B(KC)                               
            KC=KC+1                                      
         END DO
         
         X=DSQRT(A(KC)-X)                            
         C=C*X                                      
         B(KC)=X                                   
         II=I+1                                   
         IF (II.EQ.N) EXIT   ! Escape out of loop
         JC=IFIR                                
         
         DO J=II,N1                        
            JC=JC+J                              
            IC=JC                               
            KC=IFIR                            
            Y=0.0D0                           
            
            DO K=1,I                     
               Y=Y+B(IC)*B(KC)                 
               KC=KC+1                        
               IC=IC+1                       
            END DO
            
            B(IC)=(A(IC)-Y)/X            
         END DO
      END DO
   ENDIF

   RETURN                                                           
END SUBROUTINE CHOL
   
!*********************************************************************
!                                                               
!                   **********************                     
!                   *  SUBROUTINE  EXTD  *                    
!                   **********************                   
!                                                           
!    EXTRACT THE DIAGONAL ELEMENTS OF A SQUARE MATRIX AND STORE     
!    THEM IN A SEPARATE VECTOR                                     
!                                                                 
!    CALL EXTD (A,V,N,MS)                                        
!                                                               
!    A .......... INPUT MATRIX, N BY N                         
!    V .......... OUTPUT VECTOR OF LENGTH N, CONTAINS DIAGONAL OF A 
!    N .......... NUMBER OF ROWS IN A, EQUAL TO NUMBER OF COLUMNS  
!    MS ......... STORAGE MODE OF A                               
!                                                                
!*********************************************************************
SUBROUTINE EXTD(A,B,N,MS)                                      
   DOUBLE PRECISION A,B                                          
   DIMENSION A(1),B(1)                                          
   
   SELECT CASE (MS)
   CASE (0)       ! square
      L = 1                                                     
      K = N + 1                                                
      DO J=1,N                                              
         B(J) = A(L)                                            
         L = L + K                                             
      END DO
   
   CASE (1,3,4)   ! Packed, symm, upper, or lower triangle
      L = 1                                               
      DO J=1,N                                         
         B(J) = A(L)                                       
         L = L + J + 1                                    
      END DO
   
   CASE (2)       ! diagonal
      DO J=1,N                                     
         B(J) = A(J)                                   
      END DO
   END SELECT
   
   RETURN                                  
END SUBROUTINE EXTD
   
!*********************************************************************
!                                      
!                   **********************                          
!                   *  SUBROUTINE  ADJR  *                         
!                   **********************                        
!                                                                
!    ADJOIN THE ROWS OF TWO MATRICES TO PRODUCE A RESULTANT MATRIX  
!                                                                  
!    CALL ADJR (A,B,C,MA,N,MB)                                      
!                                                                  
!    A .......... INPUT MATRIX, MA BY N, GENERAL RECTANGULAR (MS=0) 
!    B .......... INPUT MATRIX, MB BY N, GENERAL RECTANGULAR (MS=0)
!    C .......... OUTPUT MATRIX, (MA+MB) BY N, RESULT OF ADJOINING, 
!                 GENERAL RECTANGULAR (MS=0)                       
!    MA ......... NUMBER OF ROWS IN A                             
!    N .......... NUMBER OF COLUMNS IN A, EQUAL TO NUMBER OF     
!                 COLUMNS IN B                                  
!    MB ......... NUMBER OF ROWS IN B                               
!                                                                  
!*********************************************************************
SUBROUTINE ADJR(A,B,C,MA,N,MB)                                   
   DOUBLE PRECISION A,B,C                                          
   DIMENSION A(MA,N),B(MB,N),C(1)                                 
   
   MAB = MA + MB                                                 
   LL = N + 1                                                   
   LK = MA + 1                                                 
   
   DO JJ=1,N                                                
      J = LL-JJ                                                 
      L = MAB*J-MB                                             
      
      DO II=1,MA                                            
         I = LK-II                                              
         C(L) = A(I,J)                                         
         L = L - 1                                            
      END DO
   END DO
   
   DO J=1,N                                          
      L = MAB*J-MB                                       
      DO I=1,MB                                       
         L = L + 1                                        
         C(L) = B(I,J)                                   
      END DO
   END DO
   
   RETURN                                    
END SUBROUTINE ADJR                                     
   
!*********************************************************************
!                                        
!                   **********************                          
!                   *  SUBROUTINE  ADJC  *                         
!                   **********************                        
!                                                                
!    ADJOIN THE COLUMNS OF TWO MATRICES TO PRODUCE A RESULTANT  
!    MATRIX                                                    
!                                                             
!    CALL ADJC (A,B,C,M,NA,NB)                               
!                                                                   
!    A .......... INPUT MATRIX, M BY NA, GENERAL RECTANGULAR (MS=0)
!    B .......... INPUT MATRIX, M BY NB, GENERAL RECTANGULAR (MS=0) 
!    C .......... OUTPUT MATRIX, M BY (NA+NB), RESULT OF ADJOINING,
!                 GENERAL RECTANGULAR (MS=0)                      
!    M .......... NUMBER OF ROWS IN A, EQUAL TO NUMBER OF ROWS IN B 
!    NA ......... NUMBER OF COLUMNS IN A               
!    NB ......... NUMBER OF COLUMNS IN B              
!                                                    
!*********************************************************************
SUBROUTINE ADJC(A,B,C,M,NA,NB)                     
   DOUBLE PRECISION A,B,C                            
   DIMENSION A(1),B(1),C(1)                         
   
   K = M*NA                                        
   KK = M*NB                                      
   
   DO L=1,K                                    
      C(L) = A(L)                                  
   END DO
   
   DO L=1,KK                                 
      JZ = L + K                                 
      C(JZ) = B(L)                              
   END DO

   RETURN                              
END SUBROUTINE ADJC

!*********************************************************************
!                   **********************                         
!                   *  SUBROUTINE MPYRT  *                        
!                   **********************                       
!                                                               
!    MULTIPLY TWO MATRICES, THE SECOND ONE ENTERING IN TRANSPOSED  
!    FORM                                                         
!                                                                
!    CALL MPYRT (A,B,C,MA,NA,MSA,MB)                            
!                                                              
!    A .......... INPUT MATRIX, MA BY NA, FIRST FACTOR IN MULTI-   
!                 PLICATION                                       
!    B .......... INPUT MATRIX, MB BY NA, TRANSPOSED SECOND FACTOR 
!                 IN MULTIPLICATION, GENERAL RECTANGULAR (MSB=0)  
!    C .......... OUTPUT MATRIX, MA BY MB, RESULT OF MULTIPLI-   
!                 CATION                                        
!    MA ......... NUMBER OF ROWS IN A                          
!    NA ......... NUMBER OF COLUMNS IN A, EQUAL TO NUMBER OF  
!                 COLUMNS IN B                               
!    MSA ........ STORAGE MODE OF A                         
!    MB ......... NUMBER OF ROWS IN B                              
!                                                                 
!*********************************************************************
SUBROUTINE MPYRT(A,B,C,MA,NA,MSA,MB)                            
   DOUBLE PRECISION A,B,C,X                                       
   DIMENSION A(1),B(1),C(1)                                      
   
   K=0                                                         
   ICA=0                                                      
   ICCA=0                                                    
   DO J=1,MB                                            
      INB=J                                                   
      JNB=INB                                                
      
      DO I=1,MA                                         
      
         SELECT CASE(MSA)
         
         CASE (0)           ! square
            INA=I                                               
            INCA=MA                                            
            LLA=1                                             
            LHA=NA                                           
            
         CASE (1,3)         ! Symmetric, packed lower
            INA=I*(I-1)/2+1                           
            INCA=1                                   
            LLA=1                                   
            IF(MSA.EQ.1) THEN
               LHA=NA                              
               ICCA=0                             
               ICA=0                             
            ELSE
               LHA=I                                 
            ENDIF
         
         CASE (2)         ! Diagonal
            INA=I                                          
            INCA=0                                        
            LLA=I                                        
            LHA=I                                       
         
         CASE (4)         ! packed upper
            INA=I*(I+1)/2                   
            INCA=I                         
            ICCA=1                        
            ICA=0                        
            LLA=I                       
            LHA=NA                     
         END SELECT
         
         K=K+1                     
         X=0.0D0                  
         
         IF(LLA .LT. 1) CYCLE   ! restart loop
         
         IF (LLA .GT. 1) INB=INB+(LLA-1)*MB    
         
         IF(NA .GE. LHA) THEN
            DO M=LLA,LHA     
               X=X+A(INA)*B(INB)  
               IF(MSA .EQ. 1 .AND. M .EQ. I) THEN
                  INCA=I 
                  ICCA=1
               ENDIF
               INA=INA+INCA+ICA                                              
               ICA=ICA+ICCA                                                 
               INB=INB+MB                                                  
            END DO
            INB=JNB                                                    
         ENDIF
         C(K)=X                                                    
      END DO
   END DO

   RETURN                                               
END SUBROUTINE MPYRT

!*********************************************************************
!                   **********************                         
!                   *  SUBROUTINE GRAMT  *                        
!                   **********************                       
!
!*********************************************************************
SUBROUTINE GRAMT(A,C,M,N)
   REAL*8 A,C,X
   DIMENSION A(M,N),C(1)
   
   IC=0
   
   DO I=1,N
      DO J=1,I
         X=0.0D0
         DO K=1,M
            X=X+A(K,I)*A(K,J)
         END DO
         IC=IC+1
         C(IC)=X
      END DO
   END DO
   
   RETURN
END SUBROUTINE GRAMT
