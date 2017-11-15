MODULE RRMSET
USE MIXLIB
CONTAINS

! *****************************************************************
!   SUBROUTINE VECS (CVEC,RVEC,NROW,NCOL)                 
!                                                         
!   EXHANGES A TWO-DIMENSIONAL ARRAY STORED BY ROWS       
!   FOR A TWO-DIMENSIONAL ARRAY STORED BY COLUMNS         
!                                                         
!   Parameters Sent                                       
!   RVEC = NROW x NCOL ARRAY STORED BY ROWS               
!   NROW = NUMBER OF ROWS                                 
!   NCOL = NUMBER OF COLUMNS                              
!                                                         
!   Parameters Returned
!   CVEC = NROW x NCOL ARRAY STORED BY COLUMNS            
! *****************************************************************
SUBROUTINE VECS(CVEC,RVEC,NROW,NCOL)
   IMPLICIT REAL*8 (A-H,O-Z)

   REAL*8  CVEC(*),RVEC(*)

   IC = 0
   DO I = 1,NCOL
      DO J = 1,NROW
         IC = IC+1
         IC2= NCOL * (J-1) + I
         CVEC(IC) = RVEC(IC2)
      END DO
   END DO
   RETURN
END SUBROUTINE VECS
   
! *****************************************************************
!   SUBROUTINE OMEGAI (NS,MAXNI,TIME,AUTO,NP,NI,NII2,NI2, 
!                      XTIME,RC,OM,OMAD,W1,W2)            
!                                                         
!   BUILDS OM (OMEGA) OMAD (AUTOCORR DERIVATIVES)
!   ASSUMING XTIME CONTAINS THE TIMEPTS OBSERVED AND      
!   TIME CONTAINS ALL POSSIBLE TIMEPOINTS                 
!                                                         
!   Parameters Sent                                       
!      NS = 0 - Stationary AR1  1 - Non Stationary AR1    
!         = 2 - Stationary MA1  3 - Stationary ARMA(1,1)  
!         = 4 or 5 - Toeplitz(NP) matrix                  
!   MAXNI = TOTAL NUMBER OF TIMEPOINTS POSSIBLE           
!   TIME  = VECTOR (MAXNI) OF ALL POSSIBLE TIMEPOINTS     
!   AUTO  = VECTOR (NP) OF AUTOCORRELATION PARAMETERS     
!    NP   = NUMBER OF AUTOCORRELATION PARAMETERS          
!    NI   = NUMBER OF TIMEPOINTS OBSERVED                 
!   NII2  = NI * NI                                       
!   NI2   = NI * NI+1 / 2                                 
!   XTIME = VECTOR (NI) OF OBSERVED TIMEPOINTS            
!                                                         
!   Parameters Returned                                   
!   OMAD  = VECTOR (NI2) OF OMEGA DERIVATIVES FOR RHO     
!   OM    = VECTOR (NI2) OF OMEGA ELEMENTS                
!   RC    = RETURN CODE                                   
!           -1 = TOO MANY MISSING TIMEPOINTS              
!            0 = OK                                       
!            1 = TOO FEW MISSING TIMEPOINTS               
!                                                         
!   Local Vectors                                         
!   OMFLAG = VECTOR (MAXNI) INDICATING WHETHER TIMEPOINT  
!            IS MISSING OR NOT FOR ALL MAXNI ELEMENTS     
!                                                         
!   Work Vectors                                          
!       W1 = VECTOR (NII2)                                
!       W2 = VECTOR (NII2)                                
!                                                         
!   Maximum Values                                        
!   MAXNI =  200   NI =  200                              
! *****************************************************************

SUBROUTINE OMEGAI(NS,MAXNI,TIME,AUTO,NP,NI,NII2,NI2,XTIME,RC, &
                  OM,OMAD,W1,W2)
   USE MIXLIB
   IMPLICIT REAL*8 (A-H,O-Z)

   REAL*8 TIME(MAXNI),XTIME(NI),AUTO(NP),OM(NI2),OMAD(NI2), &
          W1(NII2),W2(NII2)
   INTEGER RC,COUNT,COUNT2
   LOGICAL FLAG,OMFLAG(200)


   ! DETERMINE WHETHER OMEGA SHOULD EQUAL THE IDENTITY OR NOT

   IF (AUTO(1) .EQ. 0.0D0) THEN

      ! BUILD IDENTITY MATRIX AND EXIT SUBROUTINE IF AUTO(1) = 0

      COUNT = 0
      RC    = 0
      DO K1 = 1,NI
         DO K2 = 1,K1
            COUNT = COUNT+1
            IF (K1 .EQ. K2) THEN
               OM  (COUNT) = 1.0D0
               OMAD(COUNT) = 1.0D0
            ELSE
               OM  (COUNT) = 0.0D0
               OMAD(COUNT) = 0.0D0
            ENDIF
         END DO
      END DO
   ELSE
   
      ! INITIALIZE FOR BUILDING OMEGA

      MISS  = 0
      COUNT = 0
      RC    = 0
      K2    = 1

      ! CREATE OMFLAG WHICH FLAGS MISSING TIMEPOINTS

      DO K1 = 1,MAXNI
         IF (K2 .GT. NI) THEN
            FLAG = .FALSE.
            MISS = MISS + 1
         ELSE
           IF (FP_EQUAL(TIME(K1), XTIME(K2))) THEN
              FLAG = .TRUE.
              K2   = K2 + 1
           ELSE
              FLAG = .FALSE.
              MISS = MISS + 1
           ENDIF
         ENDIF
         OMFLAG(K1) = FLAG
      END DO

      ! CHECK TO SEE IF NUMBER OF MISSING TIMEPOINTS IS CORRECT

      IF (MISS - (MAXNI - NI) .LT. 0) THEN
 
         ! TOO MANY MISSING TIMEPOINTS
         RC = -1
 
      ELSEIF (MISS - (MAXNI - NI) .GT. 0) THEN
 
         ! TOO FEW MISSING TIMEPOINTS
         RC =  1

      ELSE
      
         ! CORRECT NUMBER OF MISSING TIMEPOINTS
         ! PROCEED TO BUILD VECTORS OM & OMAD
   
         COUNT  = 0
         COUNT2 = 0
         
         DO K1 = 1,MAXNI
            DO K2 = 1,K1
               COUNT = COUNT + 1
               IF (OMFLAG(K1) .AND. OMFLAG(K2)) THEN
                  COUNT2 = COUNT2 + 1
                  KK = K1 - K2
                  KK1= KK + 1 
                  KK2= KK - 1 
                  IF (NS .EQ. 0) THEN
                     AUTOSCAL = 1.0D0 / (1.0D0 - AUTO(1)**2)
                     OM  (COUNT2) = (AUTO(1)**KK) * AUTOSCAL
                     OMAD(COUNT2) = ((DBLE(KK)) * AUTO(1)**(KK2)- &
                                  (DBLE(KK-2)) * AUTO(1)**(KK1)) * &
                                  (AUTOSCAL**2)
                  ELSEIF (NS .EQ. 1) THEN
                     OM  (COUNT2) = (AUTO(1)**KK)
                     OMAD(COUNT2) = (DBLE(KK)) * AUTO(1)**(KK2)
                  ELSEIF (NS .EQ. 2) THEN
                     IF (K1 .EQ. K2) THEN
                        OM  (COUNT2) = 1.0D0 + AUTO(1)*AUTO(1)
                        OMAD(COUNT2) = 2.0D0 * AUTO(1)
                     ELSEIF (K1 .EQ. K2+1) THEN
                        OM  (COUNT2) = - AUTO(1)
                        OMAD(COUNT2) = - 1.0D0
                     ELSE
                        OM  (COUNT2) = 0.0D0
                        OMAD(COUNT2) = 0.0D0
                     ENDIF
                  ELSEIF (NS .EQ. 3) THEN
                     AUTOSCAL = 1.0D0 / (1.0D0 - AUTO(1)**2)
                     AUTOPROD = AUTO(1)*AUTO(2)
                     IF (K1 .EQ. K2) THEN
                        OM(COUNT2) = (1.0D0 + AUTO(2)**2  &
                                     - 2.0D0*AUTOPROD) * AUTOSCAL
                     ELSE
                        OM(COUNT2) = (AUTO(1)**(KK2))  &
                                  * (1.0D0 - AUTOPROD) &
                                  * (AUTO(1) - AUTO(2)) * AUTOSCAL
                     ENDIF
                     OMAD(COUNT2) = KK
                  ELSE   
                     IF (K1 .EQ. K2) THEN
                        OM  (COUNT2) = 1.0D0
                        OMAD(COUNT2) = 0.0D0
                     ELSE
                        IF(KK.LE.NP) OM(COUNT2) = AUTO(KK)
                        IF(KK.GT.NP) OM(COUNT2) = 0.0D0
                        OMAD(COUNT2) = KK
                     ENDIF
                  ENDIF
               ELSE
               ENDIF
            END DO
         END DO
         
         IF (NS .NE. 0 .AND. NS .LE. 1) THEN
            CALL MPYM(OMAD,OM,W1,NI,NI,3,4,NI)
            CALL MPYM(OM,OMAD,W2,NI,NI,3,4,NI)
            CALL ADDM(W1,W2,W1,NI,NI,0)
            CALL CHAMS(W1,OMAD,NI,0,1)
            CALL CHAMS(OM,W1,NI,3,0)
            CALL GRAM(W1,OM,NI,NI)
         ENDIF
      ENDIF
   ENDIF

   RETURN
END SUBROUTINE OMEGAI
   
! *****************************************************************
!   SUBROUTINE TRANDI(TRAN,N,N3)                          
!                                                         
!   BUILDS THE TRANSFORMATION MATRIX PSI_(n) n x n*n      
!   FROM MAGNUS (1988) page 109 WHICH TRANSFORMS          
!   w(A) (the vector containing just the diagonal elements
!   of A) into  vec A     PSI_(n) ' w(A)  = vec A         
!                                                         
!   Parameters Sent                                       
!      N = number of rows ( = # of columns) of A          
!     N3 = N * N*N                                        
!                                                         
!   Parameters Returned
!   TRAN = N x N*N VECTOR  PSI_(n)                        
! *****************************************************************
SUBROUTINE TRANDI(TRAN,N,N3)
   IMPLICIT REAL*8 (A-H,O-Z)

   REAL*8  TRAN(N3)

   N2 = N * (N+1)
   I  = 0
   
   DO
      I  = I+1
      
      TRAN(I) = 1.0D0
      IF (I .EQ. N3) EXIT
      
      DO J = 1,N2
         I = I+1
         TRAN(I) = 0.0D0
      END DO
   END DO
   
   RETURN
END SUBROUTINE TRANDI
   
! *****************************************************************
!   SUBROUTINE TRANLT(TRAN,N,NN)                          
!                                                         
!   BUILDS THE NUN x N*N TRANSFORMATION MATRIX            
!   L tilde _(n)                                          
!   FROM MAGNUS (1988) page 96 WHICH TRANSFORMS           
!   v tilde (A) (the vector containing the unique elements
!   of a lower triangular matrix A) into vec A            
!   L tilde _(n) ' v tilde (A) = vec A                    
!                                                         
!   Parameters Sent                                       
!      N = number of rows ( = # of columns) of A          
!     NN = 1/2 N (N-1) *  N*N                             
!                                                         
!   Parameters Returned
!   TRAN = NN VECTOR  L tilde _(n)                        
! *****************************************************************
SUBROUTINE TRANLT(TRAN,N,NN)
   IMPLICIT REAL*8 (A-H,O-Z)

   REAL*8  TRAN(NN)

   NUN = (N * (N-1)) / 2
   I  = 0
   DO J=1,N
      IF (J .EQ. N) THEN
          I0 = (J-1) * NUN
      ELSE
          I0 = J * NUN
      ENDIF
 
      DO J2=1,I0
         I = I+1
         TRAN(I) = 0.0D0
      END DO
 
      IF (I .EQ. NN) CYCLE  ! go back to the top of the loop
 
      I1 = N-J
      DO J3=1,I1
         I = I+1
         TRAN(I) = 1.0D0
         DO J4=1,NUN
            I = I+1
            TRAN(I) = 0.0D0
         END DO
      END DO
   END DO
   RETURN
END SUBROUTINE TRANLT
   
! *****************************************************************
!   SUBROUTINE TRANNN(TRAN,N,N2,N2N2)                     
!                                                         
!   BUILDS THE N*N x N*N TRANSFORMATION MATRIX            
!   N_(n) = 1/2 ( I_(n*n) + K_(nn) )                      
!   FROM MAGNUS (1988) page 48 WITH THE PROPERTY          
!   N_(n) vec A  =  vec 1/2( A + A' )                     
!                                                         
!   Parameters Sent                                       
!      N = number of rows ( = # of columns) of A          
!     N2 = N*N                                            
!   N2N2 = N2 * N2                                        
!                                                         
!   Parameters Returned
!   TRAN = N2N2 VECTOR  N_(n)                             
! *****************************************************************
SUBROUTINE TRANNN(TRAN,N,N2,N2N2,EI,RIN,EIN)
   USE MIXLIB
   IMPLICIT REAL*8 (A-H,O-Z)

   REAL*8  TRAN(N2N2),RIN(N2),EI(N),EIN(N2N2)

   CALL GEND(RIN,1.0D0,N,0)
   CALL GEND(TRAN,0.0D0,N2,0)
 
   DO I=1,N
      DO J=1,N
         IF (I .EQ. J) THEN
            EI(J) = 1.0D0
         ELSE
            EI(J) = 0.0D0
         ENDIF
      END DO
 
      CALL KMPY(EI,RIN,EIN,N,1,0,N,N)
      CALL KMPY(EIN,EI,EIN,N2,N,0,1,N)
      CALL ADDM(TRAN,EIN,TRAN,N2,N2,0)
   END DO

   L = 0
   DO I=1,N2
      DO J=1,N2
         L = L+1
         IF (I .EQ. J) TRAN(L) = TRAN(L) + 1.0D0
         TRAN(L) = 0.5D0 * TRAN(L)
      END DO
   END DO
   
   RETURN
END SUBROUTINE TRANNN
   
! *****************************************************************
!   SUBROUTINE TRANTO(GAM,R,N,N2,NN,N2R,EIN,EI3N,IDN,     
!                     U,E,ETEMP,GTEMP)                    
!                                                         
!   BUILDS THE N*N x R  TRANSFORMATION MATRIX  GAM_(n,r)  
!   FROM MAGNUS (1988) page 117 WITH THE PROPERTY         
!   GAM_(n,r) * gamma_(r) A  =  vec A                     
!                                                         
!   where gamma_(r) A  contains the r unique elements     
!   of the symmetric r-Toeplitz n x n matrix              
!                                                         
!   Parameters Sent                                       
!      N = number of rows ( = # of columns) of A          
!      R = number of unique elements of A                 
!     N2 = N*N                                            
!     NN = (N * (N-1)) / 2                                
!    N2R = N2 * R                                         
!                                                         
!   Work Vectors                                          
!    EIN(N)  EI3N(N)                                      
!     U(N2)    E(N2)   ETEMP(N2)                          
!    IDN(NN)           GTEMP(N2R)                         
!                                                         
!   Parameters Returned                                   
!   GAM  = N2R  TRANSFORMATION VECTOR                     
! *****************************************************************
SUBROUTINE TRANTO(GAM,R,N,N2,NN,N2R,EIN,EI3N,IDN,U,E,ETEMP,GTEMP)
   USE MIXLIB
   IMPLICIT REAL*8 (A-H,O-Z)

   REAL*8  GAM(N2R),GTEMP(N2R),IDN(NN),EIN(N),EI3N(N),U(N2),E(N2), &
           ETEMP(N2)
   INTEGER R

   CALL GEN(GAM,0.0D0,N2,R,0)
   CALL GEND(IDN,1.0D0,N,1)

   DO I=1,R
      I1  = I-1
      
      IF (I1 .EQ. 0) THEN
         CALL GEND(U,1.0D0,N,0)
      ELSE
         N1  = N  - I1 
         CALL GEN(E,0.0D0,N,N,0)
         
         DO I2=1,N1
            I3  = I2 + I1
            CALL EXTC(IDN,I2,EIN,N,1)       
            CALL EXTC(IDN,I3,EI3N,N,1)       
            CALL MPYRT(EIN,EI3N,ETEMP,N,1,0,N)
            CALL ADDM(E,ETEMP,E,N,N,0)
         END DO
         
         CALL TRP(E,ETEMP,N,N)
         CALL ADDM(E,ETEMP,U,N,N,0)
      ENDIF  
 
      CALL EXTC(IDN,I,EIN,N,1)
      CALL MPYRT(U,EIN,GTEMP,N2,1,0,R)
      CALL ADDM(GAM,GTEMP,GAM,N2,R,0)
   
   END DO
   
   RETURN
END SUBROUTINE TRANTO
   
! *****************************************************************
!   SUBROUTINE TRANGG(TRAN,N,NNN2,N2N2,TRANN)             
!                                                         
!   BUILDS THE TRANSFORMATION MATRIX  G_(n) n*n x nn      
!   FROM McCULLOCH (1982) JASA pp 679 WHICH TRANSFORMS    
!   vech(A) into vec(A)      G_(n) vech(A) = vec(A)       
!   FOR THE SYMMETRIC MATRIX A                            
!                                                         
!   Parameters Sent                                       
!      n = number of rows ( = # of columns) of A          
!   nnn2 = n*n  x  (n * (n+1)) / 2                        
!   n2n2 = n * n  x  n * n                                
!  TRANN = transformation matrix N_(nn) of size (n2n2)    
!                                                         
!   Parameters Returned
!   TRAN = n*n x (n * (n+1))/2  vector  G_(n)             
! *****************************************************************
SUBROUTINE TRANGG(TRAN,N,NNN2,N2N2,TRANN)
   IMPLICIT REAL*8 (A-H,O-Z)

   REAL*8  TRAN(NNN2),TRANN(N2N2)

   N2 = N * N
   IC = 0
   DO I1 = 1,N
      IND = (I1 - 1) * (N**3)
      
      DO I2 = 1,I1
         DO I3 = 1,N2
            IC = IC +1
            IND= IND+1
            TRAN(IC) = TRANN(IND)
       
            ! trapping values of TRAN which equal 0.5 and setting them
            ! equal to the correct value of 1.0
       
            IF (TRAN(IC) .GT. 0.1D0 .AND. TRAN(IC) .LT. 0.9D0) &
                TRAN(IC) = 1.0D0
         END DO
      END DO
   END DO
   
   RETURN
END SUBROUTINE TRANGG
   
! *****************************************************************
!   SUBROUTINE TRANGA(TRAN,R,N2,N2R,TRANI)                
!                                                         
!   BUILDS A TRANSFORMATION MATRIX TRAN_(n,r) n*n x r     
!                                          OR nn  x r     
!                                  where nn = n(n+1)/2    
!   FOR THE FOLLOWING CASE                                
!   SUPPOSE A VECTOR v(A) CONTAINS ALL POSSIBLE ELEMENTS  
!   FOR THE SYMMETRIC MATRIX A (where in general the order
!   of v(A) = r  will be less than (n * (n+1))/2 )        
!   THEN   TRAN_(n,r) v(A) = vec A  IF N2 = N*N           
!                          = vechA  IF N2 = (N*(N+1))/2   
!                                                         
!   Parameters Sent                                       
!      n = number of rows ( = # of columns) of A          
!      r = number of unique elements of A                 
!     n2 = n * n  OR  n(n+1)/2                            
!    n2r = n2* r                                          
!  TRANI = matrix (n2) which contains the real index      
!          i = 1 .. r  indicating which of the r          
!          parameters occupies each element of vec A      
!                                                         
!   Parameters Returned                                   
!   TRAN = n2  x  r vector                                
! *****************************************************************
SUBROUTINE TRANGA(TRAN,R,N2,N2R,TRANI)
   IMPLICIT REAL*8 (A-H,O-Z)

   REAL*8  TRAN(N2R),TRANI(N2)
   INTEGER R

   IC = 0
   DO I1 = 1,R
      IC2 = 0
      DO I2 = 1,N2
         IC = IC + 1
         IC2= IC2+1
         IF (INT(TRANI(IC2)) .EQ. I1) THEN
             TRAN(IC) = 1.0D0
         ELSE
             TRAN(IC) = 0.0D0
         ENDIF
      END DO
   END DO
   RETURN
END SUBROUTINE TRANGA
   
! *****************************************************************
!   SUBROUTINE ARMADER(D,NS,N2,N2S,LAG,A)                 
!                                                         
!   BUILDS THE DERIVATIVE MATRIX     D_(n,ns) n*n x ns    
!                                    or n*(n+1)/2 x ns    
!   FOR THE ARMA(1,1) CASE                                
!                                                         
!   Parameters Sent                                       
!      n = number of timepoints                           
!     ns = number of ARMA parameters ( = 2)               
!          the first one is for AR and the second for MA  
!     n2 = n * n  OR  n(n+1)/2                            
!    n2s = n2* ns                                         
!    LAG = matrix (n2) which contains the lag index       
!                                                         
!   Parameters Returned                                   
!      D = n*n x ns  derivative vector                    
!          or n(n+1)/2 x ns (depending on what n2 is)     
! *****************************************************************
SUBROUTINE ARMADER(D,NS,N2,N2S,LAG,A)
   IMPLICIT REAL*8 (A-H,O-Z)

   REAL*8 D(N2S),LAG(N2),A(NS),AP,AS,AR,MA,AR2,MA2,AS2,MAR

   IF (NS .EQ. 2) THEN

      AR  = A(1)
      MA  = A(2)
      AP  = AR*MA
      AR2 = AR*AR
      MA2 = MA*MA
      MAR = MA-AR
      AS  = 1.0D0 /  (1.0D0 - AR2)
      AS2 = 1.0D0 / ((1.0D0 - AR2)**2)

      IC = 0
      DO I1 = 1,N2
         IC  =  IC + 1
         IC2 =  IC + N2
         IF (INT(LAG(IC)) .EQ. 0) THEN
             D(IC)  = 2.0D0*AS2*(MAR * (AP-1.0D0) - (AP * AR))
             D(IC2) = 2.0D0 * AS  *  MAR
         ELSE
             K = INT(LAG(IC))
             K1= K-1
             K3= K-3
             ! K2= K-2
             ! D1= AR**(K1) * AS2
             ! D2= (1.0D0 + MA2) * (DBLE(K) - DBLE(K2)*AR2)
             ! D3= AP*(4.0D0 - DBLE(K1)*((1.0D0/AR2) - AR2))
             ! D(IC)  = D1 * (D2 - D3)
             D1= AR**(K1) * AS2
             D2= (1.0D0 - AR2) * (1.0D0 - 2.0D0*AP + MA2)
             D3= (DBLE(K1)/AR - DBLE(K3)*AR)
             D4= (AR - AP*AR - MA + AP*MA)
             D(IC)  = D1 * (D2 + D3*D4)
             D(IC2) = (AR**(K1) * AS) * (2.0D0*AP - AR2 - 1.0D0)
         ENDIF
      END DO
   ELSE 
      WRITE(6,*)'Wrong number of parameters for ARMA(1,1)'
   ENDIF
   
   RETURN
END SUBROUTINE ARMADER
   
!*************************************************************
!  SUBROUTINE PARTIV (A,NA,NA2,B,NAC,C,NC,NC2,       
!                     D,ND2,W1,W2,W3,IER)            
!                                                    
!  assumes A and C symmetric    B rectangular        
!                                                    
!  sends back D which is the symmetrix inverse of    
!              -------------                         
!             |  A   |  B ' |                        
!             |------|------|                        
!             |  B   |  C   |                        
!              -------------                         
!                                                    
!  IF PROBLEMS OCCURED IN THE INVERSION ==> IER = 1  
! ************************************************************

SUBROUTINE PARTIV(A,NA,NA2,B,NAC,C,NC,NC2,D,ND2,W1,W2,W3,IER)

   USE MIXLIB
   IMPLICIT REAL*8 (A-H,O-Z)

   REAL*8 A(NA2),C(NC2),D(ND2),B(NAC),W1(NA2),W2(NAC),W3(NC2)

   IER = 0
   CALL INVS(A,NA,DET,W1,IER)
   CALL GRAMM(B,A,W3,NC,NA,1,W1)
   CALL SUBM(C,W3,C,NC,NC,1)
   CALL INVS(C,NC,DET,W3,IER)

   CALL MPYM(B,A,W2,NC,NA,0,1,NA)
   CALL MPYM(C,W2,B,NC,NC,1,0,NA)
   UNON = -1.0D0
   CALL SCM(B,UNON,B,NC,NA,0)

   CALL GRMMT(W2,C,W1,NC,NA,1,W3)
   CALL ADDM(A,W1,A,NA,NA,1)

   CALL ADJRC(A,B,C,D,NC,NA)

   RETURN
END SUBROUTINE PARTIV
   
!*************************************************************
!  SUBROUTINE GFACT  (A,L,D,N,N2)                    
!                                                    
!  assumes A symmetric in packed form  N2 = N(N+1)/2 
!                                                    
!  sends back L lower triangle  D diagonal           
!  performing the Gaussian Factorization of the      
!  symmetric matrix A                                
!   A = L D L'                                       
! ************************************************************
SUBROUTINE GFACT(A,L,D,N,N2)

   IMPLICIT REAL*8 (A-H,O-Z)
   REAL*8 A(N2),L(N2),D(N)

   CALL CHOL(A,L,N,DET)

   ! get diagonal D from the Cholesky factor L

   IC = 0
   DO I = 1,N
      DO J = 1,I
         IC = IC + 1
         IF (I .EQ. J) D(I) = L(IC)*L(IC)
      END DO
   END DO

   ! rescale the Cholesky Factor L to be strictly lower triangular

   IC = 0
   DO I = 1,N
      DO J = 1,I
         IC = IC + 1
         L(IC) = L(IC) / DSQRT(D(J))
      END DO
   END DO

   RETURN
END SUBROUTINE GFACT
   
! *****************************************************************
!     Subroutine AUTOSPEC (modified from Jenkins & Watts   
!                          pg. 310-311 to use matrices)    
!                                                        
!     calculates smoothed autospectrum frequency points 
!     and resulting smoothed autocorrelations          
!                                                     
!       COR = Vector (M) inputing autocorrelations   
!             and returning smoothed autocorrelations     
!      SPEC = Vector (NF) returns smoothed autospectral  
!             estimate                                  
!         M = Number of smoothed autocorrelations to be
!             computed                                
!        NF = Number of frequency points for spectrum
!             NF = (2 or 3 * (M+1)) + 1             
!         W = Vector (N * M) of weights            
!                                                 
! *****************************************************************
SUBROUTINE AUTOSPEC (COR,SPEC,W,M,NF,MNF)
   USE MIXLIB

   IMPLICIT REAL*8 (A-H, O-Z)
   DOUBLE PRECISION COR(M),W(MNF),SPEC(NF),PI
   PARAMETER (PI = 3.141592654D0)

   !   calculate the weight matrix and the smoothed autocorrelations

   IC = 0
   DO K = 1,M
      WEIGHT  = 0.5D0 * (1.0D0 + DCOS((PI * DBLE(K)) / DBLE(M+1)))
      COR(K)  = COR(K)* WEIGHT
      DO I = 1,NF
         I1      = I - 1
         IC      = IC+ 1
         W(IC)   = DCOS((PI*DBLE(K)*DBLE(I1))/DBLE(NF-1))
      END DO
   END DO

   !   calculate smoothed autospectral estimate at each frequency point

   CALL MPYM(W,COR,SPEC,NF,M,0,0,1)

   DO I = 1,NF
      SPEC(I) = 2.0D0 + 4.0D0 * SPEC(I)
   END DO

   RETURN
END SUBROUTINE AUTOSPEC
   
!  *****************************************************************
!     Subroutine TOECHECK (AUTO,NP,MAXP,TOEMAT,TOECHO,NONPOS) 
!                                                          
!     checks to see whether TOEMAT is positive definite    
!                                                          
!      AUTO = Vector (NP) of non-zero autocorrelations     
!      MAXP = number of total timepoints                   
!     TOEMAT= Symmetric Toeplitz matrix (maxp * maxp+1 /2) 
!             which is formed from AUTO and checked for    
!             positive definiteness                        
!     TOECHO= Cholesky of Toeplitz matrix TOEMAT           
!     NONPOS= 0 if not positive definite                   
!             1 if positive definite                       
!                                                          
!  *****************************************************************
SUBROUTINE TOECHECK (AUTO,NP,MAXP,TOEMAT,TOECHO,NONPOS, &
                     iun,nd,xlab,head)
   USE MIXLIB
   IMPLICIT REAL*8 (A-H, O-Z)
   DOUBLE PRECISION AUTO(NP),TOEMAT(*),TOECHO(*)
   character*4 xlab(*),head(*)

   IC = 0
   DO K1 = 1,MAXP
      DO K2 = 1,K1
         IC = IC + 1
         KK    = K1 - K2
         IF (K1 .EQ. K2) THEN
            TOEMAT(IC) = 1.0D0
         ELSE
            IF(KK.LE.NP) TOEMAT(IC) = AUTO(KK)
            IF(KK.GT.NP) TOEMAT(IC) = 0.0D0
         ENDIF
      END DO
   END DO

   CALL PRNT(IUN,toemat,maxp,maxp,1,XLAB,XLAB,ND,HEAD,1,80,7,1,1, &
      "toeplitz autocorrelation matrix")
   CALL CHOL(TOEMAT,TOEMAT,MAXP,DET)
   WRITE(IUN,90)DET
   90 FORMAT(1X,' The determinant is ',F12.6)
   CALL CHSKY(TOEMAT,TOECHO,MAXP,NONPOS)

   RETURN
END SUBROUTINE TOECHECK
   
!  *****************************************************************
!     Subroutine SYMTGG (A,N,NCOL)                         
!                                                          
!     performs A' G'G                                      
!     ==> leaves A in original form  (nrow x ncol)         
!                                                          
!      AUTO = NROW x NCOL INPUT-OUTPUT MATRIX              
!      NROW = (N*(N+1))/2 number of rows and columns      
!             in G'G                                       
!      NCOL = number of columns in A                       
!                                                          
!  *****************************************************************
SUBROUTINE SYMTGG (A,N,NCOL)                          

   IMPLICIT REAL*8 (A-H, O-Z)
   DOUBLE PRECISION A(*)

   NROW = (N*(N+1))/2
   IROW = 0 

   DO I = 1,N 
      DO I2= 1,I
    
         IROW = IROW + 1
   
         IF (I2.LT.I) THEN
            DO J=1,NCOL
               IND = (J-1)*NROW + IROW
               A(IND) = 2.0D0 * A(IND)
            END DO
         ENDIF
      END DO
   END DO

   RETURN
END SUBROUTINE SYMTGG
   
!  *****************************************************************
!     Subroutine AVECHA (A,B,C,NRA,NCA,NCB,WNCA,NOC)       
!                                                          
!     for vech[ A B_1 A'  A B_2 A'   ....  A B_ncb A ' ]   
!                                                          
!     ==> each column of B is assumed to be a nca x nca    
!     ==> symmetric matrix                                 
!                                                          
!     ==> A is assumed to be a nra x nca rectangular matrix
!                                                          
!                                                          
!      A    = nra x nca rectangular input matrix           
!      B    = (nca*(nca+1)/2) x ncb input matrix           
!             with a vech in each column                   
!      C    = (nra*(nra+1)/2) x ncb output matrix          
!             with a vech in each column                   
!             UNLESS NOC=1, AND THEN                       
!      C    = nra x ncb output matrix                      
!             the off-diagonal elements are removed        
!      WNCA = work vector of size nca                      
!                                                          
!  *****************************************************************
SUBROUTINE AVECHA (A,B,C,NRA,NCA,NCB,WNCA,NOC,NOCOV)                          

   USE MIXLIB
   IMPLICIT REAL*8 (A-H, O-Z)
   DOUBLE PRECISION A(*),B(*),C(*),WNCA(NCA)

   NRB  = (NCA*(NCA+1))/2
   IF (NOC .EQ. 0) THEN
      NRC  = (NRA*(NRA+1))/2
   ELSEIF (NOC .EQ. 1) THEN
      NRC  = NRA
   ENDIF

   DO I = 1,NCB 

      IND = (I-1)*NRB + 1
      IN2 = (I-1)*NRC + 1
      CALL GRAMM(A,B(IND),C(IN2),NRA,NCA,1,WNCA)
      IF (NOCOV .EQ. 1) CALL CHAMS(C(IN2),C(IN2),NRA,1,2)
   
   END DO

   RETURN
END SUBROUTINE AVECHA
   
!  *****************************************************************
!     Subroutine GRAMD  (A,B,NRA,NCA)                      
!                                                          
!     B = diag(AA')                                        
!                                                          
!                                                          
!      A    = nra x nca rectangular input  matrix (MS=0)   
!      B    = (nra*nra) diagonal    output matrix (MS=2)   
!                                                          
!  *****************************************************************
SUBROUTINE GRAMD (A,B,NRA,NCA)                          

   IMPLICIT REAL*8 (A-H, O-Z)
   DOUBLE PRECISION A(*),B(NRA)

   IC = 0 
   DO I = 1,NRA          
      B(I) = 0.0D0 
      DO J = 1,NCA           
 
         IC = IC + 1
         B(I) = B(I) + A(IC)*A(IC)
 
      END DO
   END DO

   RETURN
END SUBROUTINE GRAMD
   
!  *****************************************************************
!     Subroutine GRAMTD (A,B,NRA,NCA)                      
!                                                          
!     B = diag(A'A)                                        
!                                                          
!                                                          
!      A    = nra x nca rectangular input  matrix (MS=0)   
!      B    = (nca*nca) diagonal    output matrix (MS=2)   
!                                                          
!  *****************************************************************
SUBROUTINE GRAMTD (A,B,NRA,NCA)                          

   IMPLICIT REAL*8 (A-H, O-Z)
   DOUBLE PRECISION A(*),B(NCA)

   IC = 0 
   DO I = 1,NCA          
      B(I) = 0.0D0 
      DO J = 1,NRA           
 
         IC = IC + 1
         B(I) = B(I) + A(IC)*A(IC)
         
      END DO
   END DO

   RETURN
END SUBROUTINE GRAMTD
   
!  *****************************************************************
!     Subroutine GRMMTD (A,B,C,NRA,NCA)                    
!                                                          
!     C = diag(A'BA)                                       
!                                                          
!                                                          
!      A    = nra x nca rectangular input  matrix (MS=0)   
!      B    = nra x nra symmetric   input  matrix (MS=1)   
!      C    = (nca*nca) diagonal    output matrix (MS=2)   
!                                                          
!  *****************************************************************
SUBROUTINE GRMMTD (A,B,C,NRA,NCA)                          

   IMPLICIT REAL*8 (A-H, O-Z)
   DOUBLE PRECISION A(*),B(*),C(NCA)

   DO IN = 1,NCA          
      C(IN) = 0.0D0 
      IC  = 0 
      IC2 = (IN-1)*NRA
      DO I = 1,NRA           
         DO J = 1,I           
    
           IC = IC + 1
           C(IN) = C(IN) + B(IC)*A(IC2+I)*A(IC2+J)
         
         END DO
      END DO
   END DO

   RETURN
END SUBROUTINE GRMMTD
   
END MODULE RRMSET
