!     Last change:  N     3 Nov 2003    2:06 pm
MODULE MixLib

   CHARACTER (len=80) ERROR_STACK(20)
   INTEGER NEXT_FREE,NEXT_RETRIEVED
   
CONTAINS

!
! subroutine to produce contrast matrices for
! Helmert contrasts (IH=2)
! or repeated contrasts (IH=1)
!  
! ************************************************
! SUBROUTINE PHIFN
! Calculate the probability distribution function (Intercept) for
! various distributions:
! NORMAL, LOGISTIC, Complementary Log-Log, OR Log-Log 
!   = 0       = 1                    = 2          =3 
! ************************************************
REAL(kind=8) FUNCTION PHIFN(Z,nfn)
    implicit none
    real(kind=8),intent(in)::z
    integer,intent(in)::nfn
    real(kind=8)::z2,ord,e,g
   
   SELECT CASE (nfn)
   
   CASE (0)    ! Normal distribution
     IF(Z.LT.-8.0D0) THEN
        PHIFN=0.000000000000001D0
     ELSEIF(Z.GT.8.0D0) THEN
        PHIFN=0.999999999999999D0
     ELSE
        Z2 = 0.0D0 - ((Z*Z)/2.0D0)
        ORD=DEXP(Z2)/2.506628275D0
        E=1.0D0/(1.0D0+0.2316418D0*DABS(Z))
        G=((((1.330274429D0*E-1.821255978D0)* &
           E+1.781477937D0)*E-0.356563782D0)*E+0.319381530D0)*E
        G=G*ORD
        IF(Z.LE.0.0D0)PHIFN=G
        IF(Z.GT.0.0D0)PHIFN=1.0D0-G
     ENDIF
     
   CASE (1)    ! Logistic distribution
     IF(Z.LT.-34.0D0) THEN
        PHIFN=0.000000000000001D0
        RETURN
     ELSEIF(Z.GT. 34.0D0) THEN
        PHIFN=0.999999999999999D0
        RETURN
     ELSE
        PHIFN = 1.0D0 / ( 1.0D0 + DEXP(0.0D0 - Z))
     ENDIF
     
   CASE (2)    ! Complementary Log-Log distribution
     PHIFN = 1.0D0 - DEXP(0.0D0 - DEXP(Z))
     
   CASE (3)    ! Log-Log distribution
     PHIFN = DEXP(0.0D0 - DEXP(Z))
   END SELECT
   
END FUNCTION PHIFN

! ************************************************
! SUBROUTINE PHIY
! Calculate the probability distribution function (ordinant) for
! various distributions:
! NORMAL, LOGISTIC, Complementary Log-Log, OR Log-Log 
!   = 0       = 1                    = 2          =3 
! ************************************************
real(kind=8) FUNCTION PHIY(Z,nfn)
    implicit none
    real(kind=8),intent(in)::z
    integer,intent(in)::nfn
    real(kind=8)::fn,az
   
   SELECT CASE (nfn)
   
   CASE (0)    ! Normal distribution
     AZ=Z*Z
     IF(AZ.GT.360.0D0) THEN
        PHIY=0.0D0
   !       ELSEIF(Z.LT.-10.0D0) THEN
   !          PHIY=0.0D0
     ELSE
        PHIY=(DEXP(-Z*Z/2.0D0))/2.506628275D0
     ENDIF
     
   CASE (1)    ! Logistic distribution
     FN    = 1.0D0 / ( 1.0D0 + DEXP(0.0D0 - Z))
     PHIY  = FN * (1.0D0 - FN)
     
   CASE (2)    ! Complementary Log-Log distribution
     FN    = 1.0D0 - DEXP(0.0D0 - DEXP(Z))
     PHIY  = (FN - 1.0D0 ) * (0.0D0 - DEXP(Z))
     
   CASE (3)    ! Log-Log distribution
     FN    = DEXP(0.0D0 - DEXP(Z))
     PHIY  = FN * (0.0D0 - DEXP(Z))
     
   END SELECT
   
END FUNCTION PHIY

! *****************************************************************
!     FUNCTION FP_EQUAL(A,B)
! *****************************************************************
!     This short routine does a 'safe compare' between two floating
!     point numbers.  This will get around the compiler message:
! 
! "Floating-point comparisons for equality may produce inconsistent results."
! 
!     which IS a LEGITIMATE message - it warns of a HIGH DEGREE of 
!     susceptibility to floating point roundoff errors, which should
!     be fixed!  For a quick introduction to the background, read 
!     http://www.lahey.com/float.htm   Tony Gray 12/18/00
! *****************************************************************
LOGICAL FUNCTION FP_EQUAL(A,B)
   DOUBLE PRECISION A,B, Epsilon
   PARAMETER (EPSILON = .0000005) ! required closeness of comparands
   IF (ABS(B - A) .LE. (ABS(B+A)*EPSILON)) THEN
      FP_EQUAL = .True.
   ELSE
      FP_EQUAL = .False.
   ENDIF
   
END FUNCTION FP_EQUAL

   ! Standardize variables - in each column in which the STD is not 0
SUBROUTINE STANDZ(XDAT,NR,NC,MEANX,STDX)
    implicit none
        REAL(KIND=8) :: XDAT(NR,NC),MEANX(NC),STDX(NC)
        INTEGER :: NR,NC,I,J

           do i=1,NR
               do j=1,NC
                  if (stdx(j) > 0.0d0) then 
                      XDAT(i,j) = (XDAT(i,j) - meanx(j)) / stdx(j) 
                  end if
               end do
           end do

END SUBROUTINE STANDZ

subroutine unique(IDs,n,ID)
    implicit none
    !
    ! Purpose:
    ! To sort and remove the duplicate ID numbers
    !
        INTEGER,INTENT(IN):: n
        INTEGER,INTENT(IN),DIMENSION(n):: IDs
        INTEGER,INTENT(OUT),DIMENSION(n)::ID
        INTEGER:: i,j,k,n1,swap

        ID(:)=0

        !remove the duplicate ID
        j=1
        ID(j)=IDs(j)
        do i=1,n
              if (ID(j)==IDs(i)) then
                 continue
              else
                 j=j+1
                 ID(j)=IDs(i)
              end if
        end do

        !sort the ID
        n1=COUNT(ID>0)
        do i=1,n1-1

              !find the location of minimum value
               k=i
               do j=i+1,n1
                  IF(ID(j)<ID(k)) THEN
                    k=j
                  END if
               end do

              !swap the minimum value
              IF(i/=k) THEN
                swap=ID(i)
                ID(i)=ID(k)
                ID(k)=swap
              END if
        end do

end subroutine unique

subroutine meanc(a,nobs,cols,m,s)
    !
    !purpose:
    !to calculate the mean and std in each column
    !
        INTEGER,INTENT(IN)::nobs,cols               ! cols is the # of columns
        REAL(KIND=8),INTENT(IN),DIMENSION(nobs,cols)::a   ! input matrix
        REAL(KIND=8),INTENT(OUT),dimension(cols):: m,s ! m is mean vector,s is std vector
        REAL(KIND=8),DIMENSION(nobs)::t             ! temporary array
        INTEGER:: i

        ! calculate mean
        t(:)=0.
        do i=1,cols
            t(:)=a(:,i)
            m(i)=SUM(t)/nobs
        end do

        ! calculate std
        t(:)=0
        do i=1,cols
            t(:)=(a(:,i)-m(i))**2
            s(i)=SQRT(SUM(t)/(nobs-1))
        end do

end subroutine meanc

subroutine descripc(a,nobs,cols,m,s,minv,maxv)
        !
        !purpose:
        !to calculate the mean, std, min, and max in each column
        !
        INTEGER,INTENT(IN)::nobs,cols                         ! cols is the # of columns
        REAL(KIND=8),INTENT(IN),DIMENSION(nobs,cols)::a             ! input matrix
        REAL(KIND=8),INTENT(OUT),dimension(cols):: m,s,minv,maxv ! m is mean vector,s is std vector
        REAL(KIND=8),DIMENSION(nobs)::t                       ! temporary array
        INTEGER:: i

        ! calculate mean, min, max
        t(:)=0.
        do i=1,cols
        t(:)=a(1:nobs,i)
        m(i)=SUM(t)/nobs
        minv(i)=minval(t)
        maxv(i)=maxval(t)
        end do

        ! calculate std
        t(:)=0
        do i=1,cols
        t(:)=(a(1:nobs,i)-m(i))**2
        s(i)=SQRT(SUM(t)/(nobs-1))
        end do

end subroutine descripc

SUBROUTINE CONTRAST(IH,C,R,P,NOMU,IPRIOR,NPAR,BIGH)
!
   IMPLICIT REAL*8(A-H,O-Z),INTEGER*4(I-N)
   INTEGER C, R, NOMU, RR, P, NPAR, INMU, IND, J, J2, I, IPRIOR, IH, &
      JJ, JI, IJ, II, JIND
   REAL (KIND=8), ALLOCATABLE :: H(:), CONTR(:), ZPR(:),             &
      ZRRR(:), CONTP(:), ZRRP(:), CONTRR(:), BIGHH(:)
   DOUBLE PRECISION BIGH
   DIMENSION BIGH(1)

!  for three categories 
!  gets Helmert contrast matrix (IH=2) in form like
!  .5 .5
!  -1  1
!  or repeated contrast matrix (IH=1) in form like
!  1  0
! -1  1
!  note: C = number of categories - 1
   !
   RR = R * (R+1) / 2 
   ALLOCATE(H(C*C))
   IF (IH .EQ. 1) THEN
      IND = 0
      DO J = 1,C
         J2 = J+1
         DO I = 1,C
            IND = IND+1
            IF (I.EQ.J) THEN
               H(IND) = 1.0D0
            ELSE IF (I.EQ.J2) THEN
               H(IND) = -1.0D0
            ELSE 
               H(IND) = 0.0D0 
            END IF
         END DO
      END DO
   ELSE IF (IH .EQ. 2) THEN
      IND = 0
      DO J = 1,C
         J2 = J+1
         DO I = 1,C
            IND = IND+1
            IF (I.GT.J2) THEN
               H(IND) = 0.0D0
            ELSE IF (I.EQ.J2) THEN
               H(IND) = -1.0D0
            ELSE IF (I.LT.J2) THEN
               H(IND) = 1.0D0 / DBLE(C+1 - I)
            END IF
         END DO
      END DO
   END IF
!
! get the part for MU
!
   INMU = 0
   IF (R .GT. 0 .AND. NOMU .EQ. 0 .AND. IPRIOR .EQ. 0) THEN
      ALLOCATE(CONTR(R*C*R*C))
      INMU = 1
      IND = 0
      DO JJ = 1,C
         DO IJ = 1,R
            DO JI = 1,C
               DO II = 1,R
                  EYE = 0.0D0
                  IF (II .EQ. IJ) EYE = 1.0D0
                  IND = IND+1
                  JIND = C*(JJ-1)+JI
                  CONTR(IND) = EYE * H(JIND)
               END DO
            END DO
         END DO
      END DO
      IF (P .GT. 0) THEN
         ALLOCATE(ZPR(P*C*R*C))
         CALL GEN(ZPR,0.0D0,P*C,R*C,0)
      END IF
      ALLOCATE(ZRRR(RR*C*R*C))
      CALL GEN(ZRRR,0.0D0,RR*C,R*C,0)
   END IF
!
! get the part for ALPHA
!
   IF (P .GT. 0) THEN
      ALLOCATE(CONTP(P*C*P*C))
      IND = 0
      DO JJ = 1,C
         DO IJ = 1,P
            DO JI = 1,C
               DO II = 1,P
                  EYE = 0.0D0
                  IF (II .EQ. IJ) EYE = 1.0D0
                  IND = IND+1
                  JIND = C*(JJ-1)+JI
                  CONTP(IND) = EYE * H(JIND)
               END DO
            END DO
         END DO
      END DO
      IF (R .GT. 0) THEN
          ALLOCATE(ZRRP(RR*C*P*C))
          CALL GEN(ZRRP,0.0D0,RR*C,P*C,0)
      END IF
   END IF
!
! get the part for SIGMA (assumes that IVCAT=1)
!
   IF (R .GT. 0) THEN
      ALLOCATE(CONTRR(RR*C*RR*C))
      IND = 0
      DO JJ = 1,C
         DO IJ = 1,RR
            DO JI = 1,C
               DO II = 1,RR
                  EYE = 0.0D0
                  IF (II .EQ. IJ) EYE = 1.0D0
                  IND = IND+1
                  JIND = C*(JJ-1)+JI
                  CONTRR(IND) = EYE * H(JIND)
               END DO
            END DO
         END DO
      END DO
   END IF
   DEALLOCATE(H)
! 
! produce the big transformation matrix
!

   ALLOCATE(BIGHH(NPAR*NPAR))

!
   IF       (INMU .EQ. 0 .AND. P .GT. 0 .AND. R .EQ. 0) THEN
!       matrix = contp
        CALL RELOC(CONTP,BIGH,P*C,P*C,0)
   ELSE IF  (INMU .EQ. 0 .AND. P .GT. 0 .AND. R .GT. 0) THEN
!       matrix = contp + zrrp + contrr
        CALL ADJR(CONTP,ZRRP,BIGH,P*C,P*C,RR*C)
        CALL ADJR(ZRRP,CONTRR,BIGHH,P*C,RR*C,RR*C)
        CALL ADJC(BIGH,BIGHH,BIGH,P*C+RR*C,P*C,RR*C)
        DEALLOCATE(ZRRP)
   ELSE IF  (INMU .EQ. 0 .AND. P .EQ. 0 .AND. R .GT. 0) THEN
!       matrix = contrr
        CALL RELOC(CONTRR,BIGH,RR*C,RR*C,0)
   ELSE IF  (INMU .GT. 0 .AND. P .GT. 0) THEN
!       matrix = contr + zpr + contp + zrrr + zrrp + contrr
        CALL ADJR(CONTR,ZPR,BIGH,R*C,R*C,P*C)
        CALL ADJR(ZPR,CONTP,BIGHH,R*C,P*C,P*C)
        CALL ADJC(BIGH,BIGHH,BIGH,R*C+P*C,R*C,P*C)
        CALL ADJR(ZRRR,ZRRP,BIGHH,R*C,RR*C,P*C)
        CALL ADJC(BIGH,BIGHH,BIGH,R*C+P*C,R*C+P*C,RR*C)
        CALL ADJC(BIGHH,CONTRR,BIGHH,RR*C,R*C+P*C,RR*C)
        CALL ADJR(BIGH,BIGHH,BIGH,R*C+P*C,R*C+P*C+RR*C,RR*C)
        DEALLOCATE(ZPR,ZRRR,ZRRP)
   ELSE IF  (INMU .GT. 0 .AND. P .EQ. 0) THEN
!       matrix = contr + zrrr + contrr
        CALL ADJR(CONTR,ZRRR,BIGH,R*C,R*C,RR*C)
        CALL ADJR(ZRRR,CONTRR,BIGHH,R*C,RR*C,RR*C)
        CALL ADJC(BIGH,BIGHH,BIGH,R*C+RR*C,R*C,RR*C)
        DEALLOCATE(ZRRR)
   ENDIF

   DEALLOCATE(BIGHH)
   IF (ALLOCATED(CONTR))  DEALLOCATE(CONTR)
   IF (ALLOCATED(CONTP))  DEALLOCATE(CONTP)
   IF (ALLOCATED(CONTRR)) DEALLOCATE(CONTRR)

END SUBROUTINE CONTRAST

! ************************************************
! SUBROUTINE GAMMAS
!     Implement the Gamma function, used by HERMIT
!     when calculating Gauss-Hermite quadrature. 
! ************************************************
REAL*8 FUNCTION GAMMAS(X)
   
   IMPLICIT REAL*8(A-H,O-Z),INTEGER*4(I-N)
   
   GAM(Y)=(((((((.035868343*Y-.193527818)*Y+.482199394)*Y- &
   .756704078)*Y+.918206857)*Y-.897056937)*Y+.988205891)*Y &
   -.577191652)*Y+1.0
   Z=X
   IF(Z .LE. 0.0D0 .OR. Z .GE. 70.0D0) THEN
      GAMMAS=0.
      WRITE(6,*) 'Function GAMMAS(X), X outside 1-70 range', Z
      STOP
   ENDIF
   
   IF(Z .LT. 1.0D0) THEN
      GAMMAS=GAM(Z)/Z
   ELSE IF (FP_EQUAL(Z , 1.0D0)) THEN
      GAMMAS=1.
   ELSE
      ZA=1.
      DO WHILE(Z > 1)
         Z=Z-1.
         ZA=ZA*Z
      END DO
      IF(FRACTION(Z) .EQ. 0.00D0) THEN
         GAMMAS=ZA
      ELSE
         GAMMAS=ZA*GAM(Z)
      ENDIF
   END IF
   RETURN
   
END FUNCTION GAMMAS

! ************************************************
!                        **********************                         
!                        *  SUBROUTINE CHAMS  *                         
!                        **********************                         
!                                                                       
!         CHANGE MODE OF STORAGE OF A SQUARE MATRIX FROM ONE 
!         PACKED FORM TO ANOTHER
!                                                                       
!         CALL CHAMS (A,B,N,MSA,MSB)                                    
!                                                                       
!         A .......... INPUT MATRIX, N BY N                             
!         B .......... OUTPUT MATRIX, N BY N                            
!         N .......... NUMBER OF ROWS IN A, EQUAL TO NUMBER OF COLUMNS  
!         MSA ........ STORAGE MODE OF A, 0..4                               
!         MSB ........ DESIRED STORAGE MODE OF B, 0..4, UNEQUAL TO MSA        
!                                                                       
!     Storage modes: 0 SQUARE, 1 PACKED SYMMETRIC, 2 DIAGONAL, 
!                    3 PACKED LOWER TRIANGLE
!                    4 PACKED UPPER TRIANGLE
!
! About Matrix storage modes:
! 
! Matrices are stored in one of three formats in memory:
! Square, diagonal, or packed.  For matrices with some symmetry,
! (upper or lower triangular, symmetric, or diagonal), you don't
! have to store every element because some of the elements are known
! by virtue of the symmetry.  This routine handles conversion between
! the different storage modes, taking advantage of the symmetry to 
! fill in the missing elements.  
!
! 0 - Square matrix, every element included, stored in row major order.
! 1 - Packed symmetric - Mat(I,J) = Mat(J,I), so only the upper triangle
!     is actually stored, and the elements are copied across the diagonal
!     on expansion.
! 2 - Diagonal: linear array of N elements, representing Mat(I,I), where
!     all nondiagonal elements are 0
! 3 - Packed lower triangle - Mat(I,J) = 0 for I < J, others are stored 
!     as an upper triangle (like packed symmetric) and are transposed 
!     on expansion.  Upper triangle is zeroed on expansion.
! 4 - Packed upper triangle - Mat(I,J) = 0 for I > J, and is stored as
!     a packed array.  On expansion no transpose is required, but lower
!     triangle is zeroed. 
!
! The 'packed' formats (1,3,4) are all stored in the same way, as an 
! upper triangle contiguously recorded in memory, so for 1 and 4 you get
! elements (1,1) (1,2) (2,2) (1,3) (2,3) (3,3) (1,4) (2,4) (3,4) (4,4) etc.
! and for 3 you get 
!          (1,1) (2,1) (2,2) (3,1) (3,2) (3,3) (4,1) (4,2) (4,3) (4,4) etc.
! though to pack and unpack type 3, it is first transposed and then handled
! like an upper triangle.  1/2/01 Tony Gray
!
! Modified 1/30/01 to always use array indices within bounds to -chk 
! option can be enabled in compiler.  TG
! ************************************************
SUBROUTINE CHAMS(A,B,N,MSA,MSB)                                   
   DOUBLE PRECISION A,B
   DIMENSION A(N,N),B(N,N)
   NAMELIST /indices/ N, IR, JR, K   ! For debug dumps
   
   ! To select among the conversions, take 5* the source format
   ! and add the destination format.  Some of the conversions take
   ! more than one step, which is why there's a loop.  For the 
   ! ones that take more than one step, the select variable is 
   ! modified to give the next step in the process
   MMS=5*MSA+MSB
   
   DO  ! keep looping till all conversions are done
   
      SELECT CASE (MMS)
      
      CASE (1,4) ! SQUARE TO PACKED UPPER TRIANGLE : 0->1, 0->4
         K = 0                                                             
         DO J=1,N                                                         
            DO I=1,J              
               K = K + 1         
               ! B(K,1) = A(I,J)  
               B(MOD(K-1,N) + 1,(K-1)/N +1) = A(I,J)  
            END DO
         END DO
         EXIT
         
      CASE (2)   ! SQUARE TO DIAGONAL  : 0->2
         DO J=1,N                                                      
            B(J,1) = A(J,J) 
         END DO
         EXIT
      
      CASE (3)   ! first step of SQUARE TO PACKED LOWER TRIANGLE : 0->3
                 ! or the last step of 3 -> 0
         DO J=2,N                                                        
            DO I=1,J-1
               IF(MSA.EQ.3) THEN
                  ! complete 3-> 0 by copying upper tri to lower and 
                  ! zeroing the upper.
                  B(J,I)=B(I,J)       
                  B(I,J)=0.0D0       
               ELSE
                  ! start 0-> 3 by copying lower tri to upper and 
                  ! zeroing the lower.
                  A(I,J)=A(J,I)     
                  A(J,I)=0.0D0     
               ENDIF
            END DO
         END DO
         
         IF (MSA .LT. 3) THEN
            MMS = 1  ! Finish 0-> 3 by packing upper triangle into B
         ELSEIF (MSA .EQ. 3) THEN
            EXIT
         ELSE
            MMS = -1 ! Illegal case, error
         ENDIF
      
      CASE (5)   ! SYMMETRIC TO SQUARE : 1->0
         L = N + 1                                                         
         ! Take a triangular array stored in A() w/o its empty elements
         ! and unpack it into B() and then copy elements across the 
         ! diagonal to make B() symmetric NxN.  TG 
         K = (N*(N+1))/2 + 1 ! count of elements in triang. less one 
         DO JR=N,1,-1                                                      
            DO IR=JR,1,-1          
               K = K - 1            
               ! B(IR,JR) = A(K,1)     
               B(IR,JR) = A(MOD(K-1,N) + 1,(K-1)/N +1)
               ! IF(ir > N .or. jr > n .OR. K > N) then
               !    write(6,indices) 
               ! END IF
            END DO
         END DO
         
         ! Now copy across the diagonal to make a complete array
         DO J=2,N                                                      
            DO I=1,J                                                      
               B(J,I) = B(I,J) 
            END DO
         END DO
         EXIT
         
      CASE (10)  ! DIAGONAL TO SQUARE : 2->0
         DO J=1,N                                                      
            DO I=1,N                                                      
               B(I,J) = 0.0D0         
            END DO
         END DO
         DO J=1,N                                                      
            B(J,J) = A(J,1)          
         END DO
         
         EXIT
      
      CASE (7,17,22) ! SYMMETRIC TO DIAGONAL: 1->2, 3->2, 4->2 
         DO J=1,N                                                      
            K = J*(J+1)/2
            B(J,1) = A(MOD(K-1,N) + 1,(K-1)/N +1)
         END DO
         EXIT
      
      CASE (11,13,14) ! DIAGONAL TO SYMMETRIC : 2->1, 2->3, 2->4
         L = N + 1                                                         
         K = (N*L)/2                                                       
         DO J=1,N                                                      
            M = N + 1 - J           
            L = L - 1           
            !B(K,1) = A(L,1)    
            B(MOD(K-1,N) + 1,(K-1)/N +1) = A(L,1)    
            K = K - M         
         END DO
         L = 2                                                             
         DO J=2,N                                                      
            LL = L + J - 2   
            DO I=L,LL                                                     
               !B(I,1) = 0.0D0 
               B(MOD(I-1,N) + 1,(I-1)/N +1) = 0.0D0 
            END DO
            L = L + J        
         END DO
     
         EXIT
      
      CASE (15,20) ! LOWER/UPPER TRUE TRIANGLE TO SQUARE : 3->0, 4->0
         L = N + 1                                                         
         K = (L*N)/2 + 1                                                   
         DO JR=N,1,-1
            DO IR=JR,1,-1
               K = K - 1      
               !B(IR,JR) = A(K,1) 
               B(IR,JR) = A(MOD(K-1,N) + 1,(K-1)/N +1)
            END DO
         END DO
         
         DO J=2,N                                                      
            L = J - 1           
            DO I=1,L                                                      
               B(J,I) = 0.0D0  
            END DO
         END DO
         IF(MSA.EQ.3) THEN
            MMS = 3 ! If upper true triangle, transpose to lower
         ELSE
            EXIT
         ENDIF
      
      CASE DEFAULT
         CALL POST_ERROR('ERROR MESSAGE FROM SUBROUTINE CHAMS: ' // &
                'ILLICIT COMBINATION OF STORAGE MODES')
      END SELECT
   END DO
   
END SUBROUTINE CHAMS

! ************************************************
!                  **********************                         
!                  *  SUBROUTINE  CHSKY *                         
!                  **********************                         
!   CALL CHSKY(A,B,N,NONPOS)                                      
!                                                                 
!   A .......... INPUT MATRIX, N BY N, SYMMETRIC (MSA=1)          
!   B .......... OUTPUT MATRIX, N BY N, CHOLESKY FACTOR, TRUE     
!                LOWER TRIANGULAR (MSB=3)                         
!   N .......... NUMBER OF ROWS IN A, EQUAL TO NUMBER OF COLUMNS  
!   NONPOS ..... OUTPUT VARIABLE, EQUALS  1  IF  A  IS POSITIVE-  
!                DEFINITE, EQUALS  0  OTHERWISE                   
!                                                                       
! ************************************************
SUBROUTINE CHSKY(A,B,N,NONPOS)                                    
   IMPLICIT REAL*8 (A-H,O-Z)                                         
   REAL*8 A(1),B(1)                                                  
   
   IF(A(1).LT.1.D-30) THEN
      NONPOS=0                                                          
   ELSE
      NONPOS=1                                                          
      X=DSQRT(A(1))                                                     
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
            
            X = A(KC) - X              
            
            IF(X.LT.1.D-30) THEN
               NONPOS=0               
               EXIT
            END IF
            
            X=DSQRT(X)             
            B(KC)=X               
            II=I+1               
            
            IF (II.EQ.N) EXIT
            
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
      END IF
   END IF
   
END SUBROUTINE CHSKY
   
! ************************************************
!                                                                       
!                  **********************                         
!                  *  SUBROUTINE   GEN  *                         
!                  **********************                         
!                                                                       
!   GENERATE A MATRIX ALL STORED ELEMENTS OF WHICH ARE EQUAL TO A 
!   SPECIFIED CONSTANT                                            
!                                                                 
!   CALL GEN (A,X,M,N,MS)                                         
!                                                                 
!   A .......... OUTPUT MATRIX, M BY N                            
!   X .......... INPUT CONSTANT OR VARIABLE, SPECIFIES VALUE FOR  
!                ELEMENTS OF A                                    
!   M .......... NUMBER OF ROWS IN A                              
!   N .......... NUMBER OF COLUMNS IN A                           
!   MS ......... DESIRED STORAGE MODE OF A                        
!                                                                       
! ************************************************
SUBROUTINE GEN (A,C,M,N,MS)                                       
   DOUBLE PRECISION A,C                                              
   DIMENSION A(1)                                                    
   
   SELECT CASE (MS)
   CASE (0)
      L=M*N                                                             
   CASE (1,3,4)
      L=(N*(N+1))/2                                                     
   CASE (2)
      L=N                                                               
   END SELECT
   
   DO J=1,L                                                        
      A(J)=C                                                            
   END DO
END SUBROUTINE GEN

! ************************************************
!                                                                       
!                  **********************                         
!                  *  SUBROUTINE GRMCV  *                         
!                  **********************                         
!                                                                 
!   CALCULATE THE MATRIX PRODUCT OF A VECTOR, MULTIPLY IT BY A    
!   SCALAR, AND ADD IT TO A SYMMETRIC (MS=1) MATRIX ALREADY IN    
!   MEMORY                                                        
!                                                                 
!   CALL GRMCV (A,B,X,C,N)                                        
!                                                                 
!   A .......... INPUT MATRIX, N BY N, SYMMETRIC (MSA=1)          
!   B .......... OUTPUT MATRIX, N BY N, SYMMETRIC (MSB=1), RESULT 
!                OF ADDITION                                      
!   X .......... INPUT VECTOR OF LENGTH N                         
!   C .......... INPUT VARIABLE OR CONSTANT (SCALAR)              
!   N .......... NUMBER OF ROWS IN A, EQUAL TO NUMBER OF COLUMNS, 
!                EQUAL TO LENGTH OF X                             
!                                                                       
! ************************************************
SUBROUTINE GRMCV(A,B,X,C,N)                                       
   DOUBLE PRECISION A,B,C,X                                          
   DIMENSION A(1),B(1),X(1)                                          
   IC=0                                                              
   DO I=1,N                                                      
      DO J=1,I                                                      
         IC=IC+1                                                           
         B(IC)=A(IC)+C*X(I)*X(J)                                           
      END DO
   END DO
END SUBROUTINE GRMCV

! ************************************************
!                                                                 
!                  **********************                         
!                  *  SUBROUTINE GRMMT  *                         
!                  **********************                         
!                                                                 
!   OBTAIN THE GRAMIAN MATRIX OF GENERALIZED PRODUCTS OF COLUMN   
!   VECTORS OF A SPECIFIED MATRIX BY POST- AND PREMULTIPLYING A   
!   METRIC MATRIX (SYMMETRIC OR DIAGONAL, MS=1,2) BY THE MATRIX   
!   AND ITS TRANSPOSE
!                                                                 
!   CALL GRMMT (A,B,C,M,N,MSB,W)                                  
!                                                                 
!   A .......... INPUT MATRIX, M BY N, GENERAL RECTANGULAR (MSA=0)
!   B .......... INPUT MATRIX, M BY M, METRIC MATRIX, SYMMETRIC OR
!                DIAGONAL                                         
!   C .......... OUTPUT MATRIX, N BY N, GENERALIZED GRAMIAN,      
!                SYMMETRIC (MSC=1)                                
!   M .......... NUMBER OF ROWS IN A                              
!   N .......... NUMBER OF COLUMNS IN A                           
!   MSB ........ STORAGE MODE OF B                                
!   W .......... WORKING VECTOR OF LENGTH M                       
!                                                                 
! ************************************************
SUBROUTINE GRMMT(A,B,C,M,N,MSB,W)                                 
   DOUBLE PRECISION X,A,B,C,W                                        
   DIMENSION A(M,N),B(1),C(1),W(1)
   
   IF (MSB .GT. 1) THEN
   
      IC=0                                                              
      DO I=1,N                                                      
         DO J=1,I                                                      
            X=0.0D0          
            DO K=1,M        
               X=X+A(K,I)*A(K,J)*B(K) 
            END DO
            IC=IC+1     
            C(IC)=X    
         END DO
      END DO
   
   ELSEIF (MSB .EQ. 1) THEN
   
      KC=0                                                              
      DO I=1,N                                                      
         KK=0                                                              
         
         DO K=1,M                                                      
            X=0.0D0               
            IC=KK                
            DO L=1,K                                                      
               IC=IC+1          
               X=X+A(L,I)*B(IC)
            END DO
            IF (K .NE. M) THEN
               J=K+1          
               IC=IC+K       
               DO L=J,M     
                  X=X+A(L,I)*B(IC)   
                  IC=IC+L           
               END DO
            END IF
            W(K)=X                 
            KK=KK+K               
         END DO
      
         DO K=1,I                 
            KC=KC+1             
            X=0.0D0            
            DO L=1,M          
               X=X+A(L,K)*W(L)  
            END DO
            C(KC)=X            
         END DO
         
      END DO
   
   ELSE  ! MSB = 0
      CALL POST_ERROR('ERROR MESSAGE FROM SUBROUTINE GRMMT: ' // &
             'METRIC IS GENERAL RECTANGULAR')
   ENDIF
   
END SUBROUTINE GRMMT

!*************************************************
! SUBROUTINE HRECUR
!
!      PURPOSE: CALCULATE GAUSS-HERMITE QUADRATURE
! ************************************************
SUBROUTINE HRECUR(PN,DPN,PN1,X,NN)
   IMPLICIT REAL*8(A-H,O-Z),INTEGER*4(I-N)
   
   P1=1.
   P=X
   DP1=0.
   DP=1.
   
   DO J=2,NN
      FJ=J
      FJ2=(FJ-1.)/2.
      Q=X*P-FJ2*P1
      DQ=X*DP+P-FJ2*DP1
      P1=P
      P=Q
      DP1=DP
      DP=DQ
   END DO
   
   PN=P
   DPN=DP
   PN1=P1
   
END SUBROUTINE HRECUR

!*************************************************
! SUBROUTINE HERMIT
!
!      PURPOSE: CALCULATE GAUSS-HERMITE QUADRATURE
! ************************************************
SUBROUTINE HERMIT(X,A,NN,EPSQ)
   IMPLICIT REAL*8(A-H,O-Z),INTEGER*4(I-N)
   DIMENSION X(1),A(1)
   DATA PI/3.141592654D0/
   
   FN=DBLE(NN)
   N1=NN-1
   N2=(NN+1)/2
   CC=1.7724538509*GAMMAS(FN)/(2.**N1)
   S=(2.*FN+1.)**.16667
   
   DO I=1,N2
      SELECT CASE (I)
      CASE (1)
         XT=S**3-1.85575/S
      CASE (2)
         XT=XT-1.14*FN**.426/XT
      CASE (3)
         XT=1.86*XT-.86*X(1)
      CASE (4)
         XT=1.91*XT-.91*X(2)
      CASE DEFAULT
         XT=2.*XT-X(I-2)
      END SELECT
      
      CALL HROOT(XT,NN,DPN,PN1,EPSQ)
      X(I)=XT
      A(I)=CC/DPN/PN1
      NI=NN-I+1
      X(NI)=-XT
      A(NI)=A(I)
   END DO
   
   DO I=1,NN
      X(I)=-X(I)*DSQRT(2.D0)
      A(I)=A(I)*(1./DSQRT(PI))
   END DO
   
END SUBROUTINE HERMIT

!*************************************************
! SUBROUTINE HROOT
!
!      PURPOSE: CALCULATE GAUSS-HERMITE QUADRATURE
! ************************************************
SUBROUTINE HROOT(X,NN,DPN,PN1,EPSQ)
   IMPLICIT REAL*8(A-H,O-Z),INTEGER*4(I-N)
   
   DO ITER=1,10
      CALL HRECUR(P,DP,PN1,X,NN)
      D=P/DP
      X=X-D
      IF(DABS(D) .LE. EPSQ) EXIT
   END DO
   
   DPN=DP
   RETURN
END SUBROUTINE HROOT

! ************************************************
!                                                                 
!                  **********************                         
!                  *  SUBROUTINE  KMPY  *                         
!                  **********************                         
!                                                                 
!   FORM THE KRONECKER (DIRECT) PRODUCT OF TWO MATRICES           
!                                                                 
!   CALL KMPY (A,B,C,MA,NA,MS,MB,NB)                              
!                                                                 
!   A .......... INPUT MATRIX, MA BY NA, FIRST FACTOR IN          
!                KRONECKER MULTIPLICATION                         
!   B .......... INPUT MATRIX, MB BY NB, SECOND FACTOR IN         
!                KRONECKER MULTIPLICATION                         
!   C .......... OUTPUT MATRIX, (MA*MB) BY (NA*NB), KRONECKER     
!                PRODUCT                                          
!   MA ......... NUMBER OF ROWS IN A                              
!   NA ......... NUMBER OF COLUMNS IN A                           
!   MS ......... STORAGE MODE OF A, EQUAL TO STORAGE MODE OF B    
!   MB ......... NUMBER OF ROWS IN B                              
!   NB ......... NUMBER OF COLUMNS IN B                           
!                                                                 
! ************************************************
SUBROUTINE KMPY(A,B,C,P,Q,MS,R,S)                                 
   DOUBLE PRECISION A,B,C                                            
   INTEGER P,Q,R,S,EA,EB                                             
   DIMENSION A(P,Q),B(R,S),C(1)                                      
   EQUIVALENCE (EB,KK),(NC,LL)                                       
   
   IF (MS .LT. 1) THEN
      IC=P*Q*R*S                                                        
      II=Q+1                                                            
      
      DO I=1,Q                                                      
         II=II-1                                                           
         JJ=S+1                                                            
         
         DO J=1,S             
            JJ=JJ-1          
            KK=P+1          
            
            DO K=1,P       
               KK=KK-1    
               LL=R+1    
               
               DO L=1,R 
                  LL=LL-1  
                  C(IC)=A(KK,II)*B(LL,JJ) 
                  IC=IC-1                
               END DO
            END DO
         END DO
      END DO
      
   ELSEIF (MS .EQ. 1) THEN
   
      EA=(P*(P+1)/2)+1                                                  
      NEB=(R*(R+1)/2)+1                                                 
      IC=P*R                                                            
      IC=(IC*(IC+1))/2                                                  
      KR=R-1                                                            
      II=P+1                                                            
      
      DO I=1,P                                                      
         II=II-1                                                           
         KEA=EA                                                            
         EB=NEB                                                            
         KJJ=R+1                                                           
         
         DO J=1,R                                                      
            KJJ=KJJ-1
            EA=KEA 
            
            DO K=1,II
               EA=EA-1 
               JJ=KR  
               NC=R  
               IF (K .EQ. 1) NC=KJJ  
               
               DO L=1,NC            
                  IF (K .EQ. 1) THEN
                     EB=EB-1        
                  ELSEIF (L .EQ. 1) THEN
                     EB=NEB-J      
                  ELSEIF (L .LE. J) THEN
                     EB=EB-JJ     
                     JJ=JJ-1     
                  ELSE
                     EB=EB-1
                  ENDIF
                  ! C(IC)=A(EA,1)*B(EB,1)  
                  C(IC)=A(MOD(EA-1,P) + 1,(EA-1)/P +1) &
                        * B(MOD(EB-1,R) + 1,(EB-1)/R +1)
                  IC=IC-1               
               END DO
            END DO
         END DO
      END DO
   
   ELSE
      IC=P*R                                                            
      II=P+1                                                            
      
      DO I=1,P                                                      
         II=II-1                                                           
         JJ=R+1                                                            
         
         DO J=1,R                                                      
            JJ=JJ-1      
            C(IC)=A(II,1)*B(JJ,1)  
            IC=IC-1     
         END DO
      END DO
   ENDIF
   
END SUBROUTINE KMPY

! ************************************************
!                                                                 
!                  **********************                         
!                  *  SUBROUTINE   ADDM *                         
!                  **********************                         
!                                                                 
!   ADD TWO MATRICES                                              
!                                                                 
!   CALL ADDM(A,B,C,M,N,MS)                                       
!                                                                 
!   A .......... INPUT MATRIX, M BY N                             
!   B .......... INPUT MATRIX, M BY N                             
!   C .......... OUTPUT MATRIX, M BY N, RESULT OF ADDITION        
!   M .......... NUMBER OF ROWS IN A, EQUAL TO NUMBER OF ROWS IN B
!   N .......... NUMBER OF COLUMNS IN A, EQUAL TO NUMBER OF       
!                COLUMNS IN B                                     
!   MS ......... STORAGE MODE OF A, EQUAL TO STORAGE MODE OF B    
!                                                                       
! ************************************************
SUBROUTINE ADDM(A,B,C,M,N,MS)                                     
   DOUBLE PRECISION A,B,C                                            
   DIMENSION A(1),B(1),C(1)                                          
   
   SELECT CASE (MS)
   CASE (0)         ! rectangular
      K=M*N                                                            
   CASE (1,3,4)     ! packed: symmetric, upper or lower triangular
      K=(M*(M+1))/2                                                    
   CASE (2)         ! diagonal
      K=M                                                              
   END SELECT
   
   DO I=1,K                                                        
      C(I)=A(I)+B(I)                                                    
   END DO

   RETURN                                                            
END SUBROUTINE ADDM


! ************************************************
!                  **********************                         
!                  *  SUBROUTINE   SUBM *                         
!                  **********************                         
!                                                                 
!   SUBTRACT TWO MATRICES                                         
!                                                                 
!   CALL SUBM(A,B,C,M,N,MS)                                       
!                                                                 
!   A .......... INPUT MATRIX, M BY N, MINUEND                    
!   B .......... INPUT MATRIX, M BY N, SUBTRAHEND                 
!   C .......... OUTPUT MATRIX, M BY N, RESULT OF SUBTRACTION     
!   M .......... NUMBER OF ROWS IN A, EQUAL TO NUMBER OF ROWS IN B
!   N .......... NUMBER OF COLUMNS IN A, EQUAL TO NUMBER OF       
!                COLUMNS IN B                                     
!   MS ......... STORAGE MODE OF A, EQUAL TO STORAGE MODE OF B    
!                                                                       
! ************************************************
SUBROUTINE SUBM(A,B,C,M,N,MS)                                     
   DOUBLE PRECISION A,B,C                                            
   DIMENSION A(1),B(1),C(1)                                          
   
   SELECT CASE (MS)
   CASE (0)         ! rectangular
      K=M*N                                                            
   CASE (1,3,4)     ! packed: symmetric, upper or lower triangular
      K=(M*(M+1))/2                                                    
   CASE (2)         ! diagonal
      K=M                                                              
   END SELECT
   
   DO I=1,K                                                        
      C(I)=A(I)-B(I)                                                    
   END DO

   RETURN                                                            
END SUBROUTINE SUBM


! ************************************************
!                                                                 
!                  **********************                         
!                  *  SUBROUTINE   SCM  *                         
!                  **********************                         
!                                                                 
!   MULTIPLY A MATRIX BY A SCALAR                                 
!                                                                 
!   CALL SCM (A,X,B,M,N,MS)                                       
!                                                                 
!   A .......... INPUT MATRIX, M BY N                             
!   X .......... SCALAR CONSTANT OR VARIABLE, FACTOR IN MULTI-    
!                PLICATION                                        
!   B .......... OUTPUT MATRIX, M BY N, RESULT OF MULTIPLICATION  
!   M .......... NUMBER OF ROWS IN A                              
!   N .......... NUMBER OF COLUMNS IN A                           
!   MS ......... STORAGE MODE OF A                                
!                                                                       
! ************************************************
SUBROUTINE SCM(A,C,B,M,N,MS)                                      
   DOUBLE PRECISION A,B,C                                            
   DIMENSION A(1),B(1)                                               
   
   SELECT CASE (MS)
   CASE (0)         ! rectangular
      K=M*N                                                            
   CASE (1,3,4)     ! packed: symmetric, upper or lower triangular
      K=(M*(M+1))/2                                                    
   CASE (2)         ! diagonal
      K=M                                                              
   END SELECT
   
   DO I=1,K                                                     
      B(I)=C*A(I)                                                       
   END DO
   RETURN                                                            
END SUBROUTINE SCM
   
! ************************************************
!                                                                 
!                  **********************                         
!                  *  SUBROUTINE   MPYM *                         
!                  **********************                         
!                                                                 
!   MULTIPLY TWO MATRICES                                         
!                                                                 
!   CALL MPYM(A,B,C,MA,NA,MSA,MSB,NB)                             
!                                                                 
!   A .......... INPUT MATRIX, MA BY NA, FIRST FACTOR IN MULTI-   
!                PLICATION                                        
!   B .......... INPUT MATRIX, NA BY NB, SECOND FACTOR IN MULTI-  
!                PLICATION                                        
!   C .......... OUTPUT MATRIX, MA BY NB, RESULT OF MULTIPLI-     
!                CATION, GENERAL RECTANGULAR (MS=0)               
!   MA ......... NUMBER OF ROWS IN A                              
!   NA ......... NUMBER OF COLUMNS IN A, EQUAL TO NUMBER OF ROWS  
!                IN B                                             
!   MSA ........ STORAGE MODE OF A                                
!   MSB ........ STORAGE MODE OF B                                
!   NB ......... NUMBER OF COLUMNS IN B                           
!                                                                       
! ************************************************
SUBROUTINE MPYM (A,B,C,MA,NA,MSA,MSB,NB)                          
   DOUBLE PRECISION A,B,C,X                                          
   DIMENSION A(1),B(1),C(1)                                          
   
   ICA=0                                                             
   ICB=0                                                             
   ICCA=0                                                            
   ICCB=0                                                            
   LOOP=1                                                            
   II=0                                                              
   JJ=0
                                            
   ! Abort if either array is zero-sized
   IF(MA .EQ. 0 .OR. NA .EQ. 0 .OR. NB .EQ. 0) RETURN
   
   IF((MSA .GE. 3)  .AND. (MSB .EQ. 0)) LOOP=2 
   
   IF (LOOP .EQ. 1) THEN
      CALL SYNCB()   !   DEFINITION OF PARAMETERS FOR MATRIX B
      CALL SYNCA()   !   DEFINITION OF PARAMETERS FOR MATRIX A
   ELSE
      CALL SYNCA()   !   DEFINITION OF PARAMETERS FOR MATRIX A
      CALL SYNCB()   !   DEFINITION OF PARAMETERS FOR MATRIX B
   ENDIF
   
   DO
      !
      !     SYNCHRONIZATION OF PARAMETERS  
      !
      K=I+(J-1)*MA  ! linear Index in destination array
      X=0.0D0                                                           
  
      IF(.NOT.(( LLA .LT. LLB) .AND. (LHA .LT. LLB)) &
       .AND. .NOT.((LLA .GT. LLB) .AND. (LHB .LT. LLA))) THEN
     
         IF (LLA .EQ. LLB) THEN
            K1=LLA  
         ELSEIF (LLA .LT. LLB) THEN
            K1=LLB 
            K3=K1-1 
            DO M=LLA,K3                                                    
               IF((MSA .EQ. 1) .AND. (M .EQ. I)) THEN
                  INCA=I 
                  ICCA=1
               ENDIF 
               INA=INA+INCA+ICA  
               ICA=ICA+ICCA     
            END DO
            
         ELSE ! lla .GT. llb
            K1=LLA 
            K3=K1-1  
            DO M=LLB,K3                                                    
               IF((MSB .EQ. 1) .AND. (M .EQ. J)) THEN
                  INCB=J  
                  ICCB=1 
               ENDIF
               INB=INB+INCB+ICB
               ICB=ICB+ICCB   
            END DO
         ENDIF
      
         IF(LHA .LT. LHB) THEN
            K2=LHA  
         ELSE
            K2=LHB 
         ENDIF
         ! 
         ! VECTOR MULTIPLICATION AND RESETTING OF PARAMETERS 
         !
         DO M=K1,K2                                                     
            X=X+A(INA)*B(INB) 
            IF((MSA .EQ. 1) .AND. (M .EQ. I)) THEN
               INCA=I
               ICCA=1 
            ENDIF
            IF((MSB .EQ. 1) .AND. (M .EQ. J)) THEN
               INCB=J
               ICCB=1
            ENDIF
            INA=INA+INCA+ICA    
            ICA=ICA+ICCA       
            INB=INB+INCB+ICB  
            ICB=ICB+ICCB     
         END DO
         
         IF(MSB .EQ. 1) THEN
            INCB=1          
            ICCB=0         
         ENDIF
         INB=JNB                                                           
         ICB=0                                                             
         IF(MSA .EQ. 1) THEN
            INCA=1        
            ICCA=0       
         ENDIF
         INA=JNA                                                           
         ICA=0                                                             
      ENDIF
  
      C(K)=X  ! store the accumulated row/column product 
              ! in the result array
      
      IF (LOOP .EQ. 1) THEN
         IF(II .LT. MA) THEN
            CALL SYNCA()
         ELSEIF (II .EQ. MA) THEN
            IF(JJ .LT. NB) THEN
               II=0   
               CALL SYNCB()
               CALL SYNCA()
            ELSE
               EXIT
            ENDIF
         ELSE
            EXIT
         ENDIF
      ELSE             ! LOOP .EQ. 2
         IF(JJ .LT. NB) THEN
            CALL SYNCB()
         ELSEIF (JJ .EQ. NB) THEN
            IF(II .LT. MA) THEN
               JJ=0
               CALL SYNCA()
               CALL SYNCB()
            ELSE
               EXIT
            ENDIF
         ELSE
            EXIT
         ENDIF
      ENDIF
   ENDDO
   
   RETURN                                                            

   CONTAINS
   
   SUBROUTINE SYNCA
      II=II+1                                                           
      IF ((MSA .EQ. 3) .AND. (MSB .EQ. 0)) THEN
         ! Storage mode of B is 0 and A is 3.  Note that this
         ! DOES NOT MATCH the equivalent test for B in SYNCB, where
         ! A is 0 and B is 4  TG 1/2/01
         I=MA-II+1                                                         
      ELSE
         I=II                                                              
      ENDIF
      
      SELECT CASE (MSA)
      CASE (0)
         INA=I                                                             
         INCA=MA                                                           
         LLA=1                                                             
         LHA=NA                                                            
      
      CASE (1,3)
         INA=I*(I-1)/2+1                                                   
         INCA=1                                                            
         LLA=1                                                             
         IF(MSA .NE. 1) THEN
            LHA=I                                                             
         ELSE
            LHA=NA                                                            
            ICCA=0                                                            
            ICA=0                                                             
         ENDIF
      
      CASE (2)
         INA=I                                                             
         INCA=0                                                            
         LLA=I                                                             
         LHA=I                                                             
      
      CASE(4)
         INA=I*(I+1)/2                                                     
         INCA=I                                                            
         ICCA=1                                                            
         ICA=0                                                             
         LLA=I                                                             
         LHA=NA                                                            
      
      END SELECT
      
      JNA=INA                                                           
   END SUBROUTINE SYNCA
   
   SUBROUTINE SYNCB
      JJ=JJ+1                                                           
      IF ((MSA .EQ. 0) .AND. (MSB .EQ. 4)) THEN
         ! Storage mode of A is 0 and B is 4
         J=NB-JJ+1                                                         
      ELSE
         J=JJ                                                              
      ENDIF
      
      SELECT CASE (MSB)
      CASE (0)
         INB=(J-1)*NA+1                                                    
         INCB=1                                                            
         LLB=1                                                             
         LHB=NA                                                            
      
      CASE (2)
         INB=J                                                             
         INCB=0                                                            
         LLB=J                                                             
         LHB=J                                                             
      
      CASE (3)
         INB=J*(J+1)/2                                                     
         INCB=J                                                            
         ICCB=1                                                            
         LLB=J                                                             
         LHB=NA                                                            
      
      CASE(1,4)
         INB=J*(J-1)/2+1                                                   
         INCB=1                                                            
         LLB=1                                                             
         IF(MSB .NE. 1) THEN
            LHB=J                                                             
         ELSE
            LHB=NA                                                            
         ENDIF
         
      END SELECT
      JNB=INB                                                           
   END SUBROUTINE SYNCB
   
END SUBROUTINE MPYM
                                                                 
! ************************************************
!                  **********************                         
!                  *  SUBROUTINE MPDSD  *                         
!                  **********************                         
!                                                                 
!   RESCALE A SYMMETRIC (GRAMIAN) MATRIX (MS=1) BY PRE- AND POST- 
!   MULTIPLYING IT BY A DIAGONAL MATRIX (MS=2)                    
!                                                                 
!   CALL MPDSD (A,B,C,N)                                          
!                                                                 
!   A .......... INPUT MATRIX, N BY N, DIAGONAL (MSA=2), CONTAINS 
!                SCALING FACTORS                                  
!   B .......... INPUT MATRIX, N BY N, SYMMETRIC (MSB=1)          
!   C .......... OUTPUT MATRIX, N BY N, SYMMETRIC (MSC=1), RESULT 
!                OF RESCALING                                     
!   N .......... NUMBER OF ROWS AND NUMBER OF COLUMNS IN A AND B  
!                                                                 
! ************************************************
SUBROUTINE MPDSD(A,B,C,N)                                         
   DOUBLE PRECISION A,B,C                                            
   DIMENSION A(1),B(1),C(1)                                          
   IC=0                                                              
   DO I=1,N                                                      
      DO J=1,I                                                      
         IC=IC+1                                                           
         C(IC)=B(IC)*A(I)*A(J)                                             
      END DO
   END DO
END SUBROUTINE  MPDSD        
                                                                  
! ************************************************
!                  **********************                         
!                  *  SUBROUTINE MPYTR  *                         
!                  **********************                         
!                                                                 
!   MULTIPLY TWO MATRICES, THE FIRST ONE ENTERING IN ITS TRANS-   
!                POSED FORM                                       
!                                                                 
!   CALL MPYTR (A,B,C,MA,NA,MSB,NB)                               
!                                                                 
!   A .......... INPUT MATRIX, MA BY NA, TRANSPOSED FIRST FACTOR  
!                IN MULTIPLICATION, GENERAL RECTANGULAR (MSA=0)   
!   B .......... INPUT MATRIX, MA BY NB, SECOND FACTOR IN MULTI-  
!                PLICATION                                        
!   C .......... OUTPUT MATRIX, NA BY NB, RESULT OF MULTIPLICATION
!                GENERAL RECTANGULAR (MSC=0)                      
!   MA ......... NUMBER OF ROWS IN A, EQUAL TO NUMBER OF ROWS IN B
!   NA ......... NUMBER OF COLUMNS IN A                           
!   MSB ........ STORAGE MODE OF B                                
!   NB ......... NUMBER OF COLUMNS IN B                           
!                                                                       
! ************************************************
SUBROUTINE MPYTR(A,B,C,MA,NA,MSB,NB)                              
   DOUBLE PRECISION A,B,C,X                                          
   DIMENSION A(1),B(1),C(1)                                          
   K=0                                                               
   ICB=0                                                             
   ICCB=0                                                            
   DO J=1,NB                                                     
   
      SELECT CASE (MSB)
      
      CASE (0)
         INB=(J-1)*MA+1                                                    
         INCB=1                                                            
         LLB=1                                                             
         LHB=MA                                                            
      
      CASE (2)
         INB=J                                                             
         INCB=0                                                            
         LLB=J                                                             
         LHB=J                                                             
      
      CASE (3)
         INB=J*(J+1)/2                                                     
         INCB=J                                                            
         ICCB=1                                                            
         LLB=J                                                             
         LHB=MA                                                            
      
      CASE (1,4)
         INB=J*(J-1)/2+1                                                   
         INCB=1                                                            
         LLB=1                                                             
         IF(MSB .NE. 1) THEN
            LHB=J
         ELSE
            LHB=MA 
         ENDIF
      
      END SELECT
      
      JNB=INB                                                           
      
      DO I=1,NA                                                     
         INA=(I-1)*MA+1                                                    
         K=K+1                                                             
         X=0.0D0                                                           
         
         IF (LLB .GE. 1) THEN
            IF (LLB .GT. 1) INA=INA+LLB-1
            
            IF(MA .GE. LHB) THEN
               DO M=LLB,LHB             
                  X=X+A(INA)*B(INB)    
                  IF((MSB .EQ. 1) .AND. (M .EQ. J)) THEN
                     INCB=J           
                     ICCB=1          
                  ENDIF
                  INA=INA+1         
                  INB=INB+INCB+ICB 
                  ICB=ICB+ICCB    
               END DO
               IF(MSB .EQ. 1) THEN
                  INCB=1         
                  ICCB=0        
               ENDIF
               INB=JNB 
               ICB=0  
            ENDIF
         
         ENDIF
         C(K)=X                                                            
      END DO
   END DO
   RETURN                                                            
END SUBROUTINE  MPYTR


!*******************************************************
!
!             **********************
!             *  SUBROUTINE PRNT  *
!             **********************
!
!   PRINT MATRICES IN ANY STORAGE MODE WITH VARIABLE FORMAT
!
!   CALL PRNT(A,M,N,MS,RLAB,CLAB,ND,HEAD,IH,MCW,NSIG,NR,NC,
!  *60H         MATRIX TITLE                                   )
!     * IN COLUMN 6                   ) IN COLUMN 70,71, OR 72
!
!   A ......INPUT MATRIX, M BY N
!   M ......NUMBER OF ROWS IN A
!   N ......NUMBER OF COLUMNS IN A
!   MS .....STORAGE MODE OF A
!   RLAB ...INPUT STRING OF M ROW LABELS, EACH CONTAINING
!           4*NR CHARACTERS
!   CLAB ...INPUT STRING OF N COLUMN LABELS, EACH CONTAINING
!           4*NC CHARACTERS
!   ND .....INPUT-OUTPUT SCALAR, DISPLAY NUMBER
!   HEAD ...INPUT STRING CONTAINING 120 CHARACTERS OF PAGE
!           HEADING
!  IH .....=0  DO NOT PAGE EJECT OR PRINT PAGE HEADING
!          =1  PAGE EJECT AND PRINT HEADING
!          =2  ??
!  MCW.....MAXIMUM COLUMN WIDTH OF PRINTOUT(LE.132.AND.GE.40)
!  NSIG ...NUMBER OF SIGNIFICANT FIGURES TO BE PRINTED (LE.16)
!          IF NSIG.EQ.0, ONLY THE INTEGRAL PART OF THE NUMBER
!          WILL BE PRINTED AND NSIG WILL BE SET TO MAGNITUDE OF
!          GREATEST ABSOLUTE ELEMENT.
!  NR .....NUMBER OF 4-CHARACTER SEGMENTS IN ROW LABELS (LE.6)
!  NC .....NUMBER OF 4-CHARACTER SEGMENTS IN COLUMN LABELS (LE.6)
!     .....ANY TITLE UP TO 60 CHARACTERS
!
!  NOTES:  IF MCW.LE.130, THE HEADING WILL PRINT IN TWO
!          60-CHARACTER LINES
!
!          IF MCW.LT.76, THE TITLE WILL PRINT IN LINES OF 32
!          AND 28 CHARACTERS
!
!            **********************
!            *  SUBROUTINE PRNTR  *
!            **********************
!
!  PRINT THE TRANSPOSE OF A MATRIX IN MEMORY.  THE CALLING
!  SEQUENCE IS THE SAME AS IN SUBROUTINE PRNT, EXCEPT THAT THE
!  MODE OF STORAGE (MS) PARAMETER IS OMITTED.   ROW AND COLUMN
!  PARAMETERS REFER TO THE UNTRANSPOSED MATRIX.
!
SUBROUTINE PRNT(LSTUN,B,NUMR,NUMC,MS,ROW,COL,ND, &
                   HEAD,IH,NCT,NSF,NR,NC,TITLE)
   REAL*8 B,BIG
   CHARACTER (len=*) TITLE
   CHARACTER*4 HEAD,TEMP,FMT,FORM,FLB,ROW,COL,BLANK,CH
   CHARACTER*40 NEWFMT,NEWFLB
   CHARACTER*36 NEWFORM
   CHARACTER*10 OneSpace
   LOGICAL*1 X,D1,D2,TR,ITT
   DIMENSION B(1),HEAD(1),CH(36),FMT(10),FORM(9), &
             FLB(10),TEMP(20),ROW(1),COL(1)
   DATA TEMP/'I3, ','X,  ','.   ','(   ','X,99','(   ', &
      ' A4,','(I6 ',',3X ','    ','(I6,','3X, ',' 2X,',' 99D', &
      '))  ',')   ','X)) ',' 99F','1A4,',' A4 '/
   DATA CH/'1','2','3','4','5','6','7','8','9','10','11','12', &
      '13','14','15','16','17','18','19','20','21','22', &
      '23','24','25','26', &
      '27','28','29','30','31','32','33','34','35','36'/
   DATA BLANK/'    '/
   
   BIG = DABS (B(1))
   TR=.FALSE.
   MMS=MS+1
   
   !     The following entry point is for printing the transpose of a 
   !     matrix.  To enable this entry point, uncomment the following 5 
   !     lines.  TG commented this code out 1/3/01 while getting rid of 
   !     the last of the GOTO's.  
   !      GO TO 6
   !      ENTRY PRNTR (B,NUMC,NUMR,COL,ROW,ND,HEAD,IH,NCT,NSF,NC,NR,TITLE)
   !      TR=.TRUE.
   !      MMS=1
   !    6 CONTINUE
   
   IF(NUMR*NUMC.GE.22501) THEN
      CALL POST_ERROR('ERROR MESSAGE FROM SUBROUTINE PRNT: OUTPUT  MATRIX TOO LARGE')
      RETURN
   ENDIF

   OneSpace = "(1H )"  ! used for small formatting output
   
   ND=ND+1
   DO I=1,10
      FMT(I)=BLANK
      FLB(I)=BLANK
   END DO
   DO I=1,9
      FORM(I)=BLANK
   END DO
   
   IHH=IH
   ITT=.TRUE.
   NUMROW=NUMR
   NRL=NR
   NCL=NC
   NSIG=NSF
   NCOL=NCT
   
   IF(NRL.LT.0)NRL=0
   IF(NRL.GT.6)NRL=6
   IF(NCL.LT.0)NCL=0
   IF(NCL.GT.6)NCL=6
   LSIZE=4*NCL
   D1=.TRUE.
   IF(NRL.EQ.0)D1=.FALSE.
   D2=.TRUE.
   IF(NCL.EQ.0)D2=.FALSE.
   IF(NCOL.EQ.0)NCOL=132
   IF(NCOL.LT.40)NCOL=40
   IF(NCOL.GT.132)NCOL=132
   IF(NSIG.LT.0)NSIG=1
   IF(NSIG.GT.16)NSIG=16
   
   NCOL=NCOL-4*NRL-12
   IF(D1)NCOL=NCOL+3
   IF(NCOL.LT.4)NCOL=4
   
   SELECT CASE (MMS)
   CASE (1)       ! square array, MS .EQ. 0
      J=NUMROW*NUMC
      NUMCOL=NUMC
      X=.TRUE.
   CASE (2,4,5)   ! Packed array: symmetric, lower or upper triangle
      J=(NUMROW*(NUMROW+1))/2
      NUMCOL=NUMROW
      TR=.FALSE.
      X=.FALSE.
   CASE (3)       ! diagonal
      J=NUMROW
      NUMCOL=1
      TR=.FALSE.
      X=.TRUE.
   END SELECT
   
   IF(J.GT.1) THEN
      DO I=2,J
         BIG = DMAX1 (BIG,DABS(B(I)))
      END DO
   ENDIF
   
   IF(BIG.EQ.0.0D0) THEN
      MODE=1
   ELSE
      MODE = DLOG10 (DABS (BIG))
      MODE =MODE+1
   ENDIF
   
   FMT(1)=TEMP(11)
   FMT(2)=TEMP(12)
   IF(NRL.GT.0) FMT(3)=CH(NRL)
   FMT(4)=TEMP(7)
   FMT(5)=TEMP(13)
   FORM(1)=TEMP(4)
   FORM(2)=CH(11+4*NRL)
   FORM(3)=TEMP(5)
   FORM(4)=TEMP(6)
   FLB(1)=TEMP(4)
   FLB(2)=CH(11+4*NRL)
   FLB(3)=TEMP(5)
   FLB(4)=TEMP(6)
   FMT( 8)=TEMP(3)
   FMT(10)=TEMP(16)
   FORM(6)=TEMP(2)
   FORM(7)=TEMP(1)
   FORM(9)=TEMP(17)
   FLB(6)=TEMP(2)
   IF(NCL.GT.0) FLB(7)=CH(NCL)
   FLB(8)=TEMP(7)
   FLB(10)=TEMP(17)
   
   IF((MODE.GT.16) .OR. .NOT. &
      ((MODE.LE.NSIG.AND.MODE.GE.1) .OR. &
      (NSIG.EQ.0.AND.MODE.GT.0))) THEN

   !****+D-FORMAT

      NCW=NSIG+10
      MAXCOL=NCOL/NCW
      IF(MAXCOL.LT.1)MAXCOL=1
      LCEN=(NCW-LSIZE)/2+1
      NCEN=(NCW-3)/2+1
      FMT(6)=TEMP(14)
      FMT(7)=CH(NCW)
      FMT(9)=CH(NSIG)
      FORM(8)=CH(NCEN)
      FORM(5)=CH(NCW-NCEN-3)
      FLB(5)=CH(LCEN)
      FLB(9)=CH(NCW-LCEN-LSIZE)
   ELSE

   !****+F-FORMAT

      NCW=NSIG+3
      IF(NCW.LT.LSIZE+1)NCW=LSIZE+1
      MAXCOL=NCOL/NCW
      IF(MAXCOL.LT.1)MAXCOL=1
      LCEN=(NCW-LSIZE)/2+1
      NCEN=(NCW-3)/2+1
      IF(NSIG.EQ.0)NSIG=MODE
      NDEC=NSIG-MODE
      FMT(6)=TEMP(18)
      FMT(7)=CH(NCW)
      IF(NDEC.GT.0)FMT(9)=CH(NDEC)
      IF(NDEC.LE.0)FMT(9)=TEMP(10)
      
      IF(NCEN.GT.0) THEN
         FORM(8)=CH(NCEN)
      ELSE
         FORM(8)=TEMP(10)
         FORM(9)=TEMP(15)
      ENDIF
      
      IF((NCW-NCEN-3).GT.0) THEN
         FORM(5)=CH(NCW-NCEN-3)
      ELSE
         FORM(5)=TEMP(10)
         FORM(6)=TEMP(10)
      ENDIF
      
      IF(LCEN.GT.0) THEN
         FLB(5)=CH(LCEN)
      ELSE
         FLB(5)=TEMP(10)
         FLB(6)=TEMP(10)
      ENDIF
      
      IF((NCW-LCEN-LSIZE).GT.0) THEN
         FLB(9)=CH(NCW-LCEN-LSIZE)
      ELSE
         FLB(8)=TEMP(20)
         FLB(9)=TEMP(10)
         FLB(10)=TEMP(15)
      ENDIF
   ENDIF
   
   IF(.NOT. D1) THEN
      FMT(3)=TEMP(10)
      FMT(4)=TEMP(10)
      FMT(5)=TEMP(10)
      FORM(2)=CH(9)
      FLB(2)=CH(9)
   ENDIF
   
   MCM1=MAXCOL-1


   !     create new character expressions NEWFMT NEWFLB & NEWFORM
   !     removing blanks from the corresponding FMT FLB & FORM

   IC1 = 0
   IC2 = 0
   DO I1 = 1,10
      DO I2 = 1,4
         IF (FMT(I1)(I2:I2) .NE. ' ') THEN
             IF (IC1 .EQ. 0) THEN
                 NEWFMT(1:1) = FMT(I1)(I2:I2)
             ELSE
                 NEWFMT = NEWFMT(1:IC1)//FMT(I1)(I2:I2)
             END IF
             IC1 = IC1 + 1
         ELSE
         END IF
         IF (FLB(I1)(I2:I2) .NE. ' ') THEN
             IF (IC2 .EQ. 0) THEN
                 NEWFLB(1:1) = FLB(I1)(I2:I2)
             ELSE
                 NEWFLB = NEWFLB(1:IC2)//FLB(I1)(I2:I2)
             END IF
             IC2 = IC2 + 1
         ELSE
         END IF
      END DO
   END DO

   IC2 = 0
   DO I1 = 1,9
      DO I2 = 1,4
         IF (FORM(I1)(I2:I2) .NE. ' ') THEN
             IF (IC2 .EQ. 0) THEN
                 NEWFORM(1:1) = FORM(I1)(I2:I2)
             ELSE
                 NEWFORM = NEWFORM(1:IC2)//FORM(I1)(I2:I2)
             END IF
             IC2 = IC2 + 1
         END IF
      END DO
   END DO


   !    done with character manipulations

   JCRANK=NUMROW
   IF(TR)JCRANK=1
   
   DO K=1,NUMCOL,MAXCOL
   
      KSTOP= MIN0 (NUMCOL,K+MAXCOL-1)
      KK=(K-1)*NUMROW
      IF(TR)KK=K-1
      LINE=(-1)
      KSM1=KSTOP-K
      ISTART=1
      IF(.NOT. X) THEN
         JCRANK=1
         JADD=K-1
         JSTART=((JADD*K)/2)+1
         ISTART=K
      ENDIF
      
      DO II=ISTART,NUMROW
         
         LINE=LINE+1
         IF(LINE .EQ. 0) THEN
            IF(IHH.GT.1) THEN
               IHH=1
               IF(IH.EQ.2)IHH=0
               MST=MMS-1
               WRITE(LSTUN,"(1X,4H@PRN/1X,9I4)") &
                     NUMR,NUMC,MST,ND,IHH,NCT,NSF,NR,NC
            ENDIF
     
            IF(IHH.EQ.1) THEN
               NW=30
               IF(NCT.LT.120)NW=15
               IF(NCT.LT.40)NW=10
               WRITE(LSTUN,"(1H-,30A4)")(HEAD(I),I=1,NW)
               IF(NCT.LT.120.AND.NCT.GE.40) &
                  WRITE(LSTUN,"(1X,15A4)")(HEAD(I),I=16,30)
               IF(NCT.LT.40) &
                  WRITE(LSTUN,"(1X,10A4)")(HEAD(I),I=11,30)
               WRITE(LSTUN,OneSpace)
               IHH=0
            ENDIF
     
            NW=15
            IF(NCT.LT.76)NW=8
            IF(ITT) THEN
               WRITE(LSTUN,"(/1X,A60)")(TITLE(1:4*NW))
               ITT=.FALSE.
            ENDIF
            IF(NCT.LT.76) &
               WRITE(LSTUN,fmt="(16X,A60)")(TITLE(36:60))
            WRITE(LSTUN,OneSpace)
            WRITE(LSTUN,NEWFORM)(I,I=K,KSTOP)
            
            IF(D2) THEN
               L1=1+(K-1)*NC
               L2=1+(KSTOP-1)*NC+NC-1
               WRITE(LSTUN,NEWFLB)(COL(L),L=L1,L2)
            ENDIF
            
            WRITE(LSTUN,OneSpace)
         ENDIF
         JJ=MCM1
   
         IF(X) THEN
            JSTART=II+KK
            IF(TR)JSTART=(II-1)*NUMCOL+KK+1
            JSTOP=JSTART+NUMROW* MIN0 (JJ,KSM1)
            IF(TR)JSTOP=JSTART+MIN0(JJ,KSM1)
         ELSE
            JSTART=JSTART+JADD
            JSTOP=JSTART+ MIN0 (LINE,MCM1)
            JADD=JADD+1
         ENDIF
         IF(D1) THEN
            L1=1+(II-1)*NR
            L2=L1+NR-1
            WRITE(LSTUN,NEWFMT)II,(ROW(L),L=L1,L2), &
                            (B(J),J=JSTART,JSTOP,JCRANK)
         ELSE
            WRITE(LSTUN,NEWFMT)II,(B(J),J=JSTART,JSTOP,JCRANK)
         ENDIF
      END DO
      WRITE(LSTUN,OneSpace)
   END DO
   RETURN
END SUBROUTINE  PRNT
                                                                 
! ************************************************
!                  **********************                         
!                  *  SUBROUTINE RELOC  *                         
!                  **********************                         
!                                                                 
!   DUPLICATE A MATRIX                                            
!                                                                 
!   CALL RELOC (A,B,M,N,MS)                                       
!                                                                 
!   A .......... INPUT MATRIX, M BY N                             
!   B .......... OUTPUT MATRIX, M BY N, EQUAL TO A                
!   M .......... NUMBER OF ROWS IN A                              
!   N .......... NUMBER OF COLUMNS IN A                           
!   MS ......... STORAGE MODE OF A                                
!                                                                       
! ************************************************
SUBROUTINE RELOC (A,B,M,N,MS)                                     
   DOUBLE PRECISION A,B                                              
   DIMENSION A(1),B(1)                                               
   SELECT CASE (MS)
   CASE (0)       ! square matrix
      L = M*N                                                           
   CASE (1,3,4)   ! packed matrix
      L = (M*(M+1))/2                                                   
   CASE (2)       ! diagonal matrix
      L=N                                                               
   END SELECT
   DO J=1,L                                                        
      B(J) = A(J)                                                       
   END DO
   RETURN                                                            
END SUBROUTINE  RELOC          

! ************************************************
! SUBROUTINE YSAME (YVEC,N,NSAME)                            
!                                                            
! INDICATES WHETHER ALL ELEMENTS OF YVEC ARE THE SAME        
!                                                            
! Parameters Sent                                            
! YVEC  = N x 1 vector of values                             
! N     = NUMBER OF ROWS IN YVEC                             
!                                                            
! Parameters Returned                                        
! NSAME = 0 if elements are not all the same                 
!       = 1 if elements are all the same                     
! ************************************************
SUBROUTINE YSAME(YVEC,N,NSAME)
   IMPLICIT REAL*8 (A-H,O-Z)

   REAL*8  YVEC(*)

   IF (N .EQ. 1) THEN 
       NSAME = 1        
   ELSE                 
       NSAME = 0
       NSUM  = 0
       YOLD  = YVEC(1)   
       DO I = 2,N    
          IF (FP_EQUAL(YVEC(I) , YOLD)) NSUM = NSUM + 1
          YOLD = YVEC(I)
       END DO
       IF (NSUM .EQ. (N-1)) NSAME = 1
   ENDIF
   RETURN
END SUBROUTINE YSAME

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

! **********************************************************************
!                                                                 
!                  *********************                         
!                  *  SUBROUTINE  INVS  *                         
!                  **********************                         
!                                                                 
!   INVERT A SYMMETRIC MATRIX (MS=1) IN PLACE AND CALCULATE THE   
!   DETERMINANT                                                   
!                                                                 
!   CALL INVS (A,N,DET,W)                                         
!                                                                 
!   A .......... INPUT-OUTPUT MATRIX, N BY N,SYMMETRIC (MS=1)     
!   N .......... NUMBER OF ROWS IN A, EQUAL TO NUMBER OF COLUMNS  
!   DET ........ OUTPUT SCALAR, DETERMINANT OF A                  
!   W .......... WORKING VECTOR OF LENGTH N                       
!   IER......... OPTIONAL ERROR FLAG
!                                                                       
! **********************************************************************

SUBROUTINE INVS(A,N,C,W,IER,VERBOSE)                                      
   DOUBLE PRECISION U,X,Y,Z,D,A,C,W                                  
   DIMENSION A(1), W(1)                                             
   INTEGER DIAGMK                                                    
   INTEGER DIAG,DIAG2,ROWNO,ROWCOL                                   
   INTEGER COLNO
   INTEGER, INTENT(in out),OPTIONAL:: IER
   LOGICAL, INTENT(in),OPTIONAL:: VERBOSE

   IF (PRESENT(IER)) THEN
      IER = 0
   END IF 

   D=A(1)                                                            
   
   IF(D .NE. 0)  THEN
   
      A(1)=1.0D0/D                                                      
      
      IF(N .GT. 1) THEN
         
         DIAG=1                                                            
         
         DO K=2,N                                                     
            KM1=K-1             
            DIAGMK=DIAG        
            DIAG=DIAG+K       
            U=A(DIAG)        
            COLNO=DIAGMK    
            DIAG2=0        
            
            DO I=1,KM1                                                     
               X=0.0D0    
               COLNO=COLNO+1         
               ROWNO=DIAGMK         
               J=1                 
               ROWCOL=DIAG2       
               DO WHILE (J .LT. I)
                  ROWCOL=ROWCOL+1    
                  ROWNO=ROWNO+1     
                  Y=A(ROWCOL)      
                  Z=A(ROWNO)      
                  X=X+Y*Z        
                  J=J+1         
               END DO
               
               ROWCOL=ROWCOL+1   
               
               DO WHILE (J .LT. K)
                  ROWNO=ROWNO+1        
                  Y=A(ROWCOL)         
                  Z=A(ROWNO)         
                  X=X+Y*Z           
                  ROWCOL=ROWCOL+J  
                  J=J+1           
               END DO
               
               W(I)=-X          
               Y=A(COLNO)      
               U = U-X*Y      
               DIAG2=DIAG2+I 
            END DO
            D=D*U
            
            IF(U .NE. 0) THEN
               A(DIAG)=1.0D0/U   
               ROWNO=DIAGMK     
               DIAG2=0         
               
               DO I=1,KM1     
                  ROWNO=ROWNO+1  
                  DIAG2=DIAG2+I 
                  X=W(I)       
                  X=X/U       
                  A(ROWNO)=X 
                  ROWCOL=DIAG2      
                  DO J=I,KM1     
                     Y=W(J)             
                     Z=A(ROWCOL)       
                     A(ROWCOL)=Z+X*Y  
                     ROWCOL=ROWCOL+J 
                  END DO
               END DO
            ENDIF
         END DO
      ENDIF
   ENDIF
   
   C = D
   IF(D .EQ. 0) THEN
      IF(PRESENT(VERBOSE)) THEN
         IF(VERBOSE) CALL POST_ERROR('ERROR MESSAGE FROM SUBROUTINE INVS: MATRIX IS SINGULAR')
      ELSE
         CALL POST_ERROR('ERROR MESSAGE FROM SUBROUTINE INVS: MATRIX IS SINGULAR')
      ENDIF
      IF (PRESENT(IER)) THEN
         IER = 1
      END IF 
   ENDIF
   RETURN                                                            
END SUBROUTINE INVS

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
 
! *****************************************************************
! SUBROUTINE READAT (FILEDAT,N,NTOT,MAXK,MAXCOL,R,P,ALLDAT,  
!                     IDNI,IDIND,YIND,XIND,WIND,
!                     MISS,YMISS,XMISS,WMISS,
!                     IWT,WTIND,IOFS,OFIND)    
!                                                            
! READ DATA FROM FILE  filedat  IN FREE FORMAT               
! ID     Y     X(R)     W(P)     WT   OFFSET                 
! in columns specified by                                    
! idind  yind  xind(r)  wind(p)  wtind  ofind                
!     
! The last four parameters are optional.  If you want to
! specify a weighting column, include IWT (flag 0|1) and 
! and WTIND.  If you want to use offset/censoring 
! which are mathematically identical,  include IOFS and 
! OFIND in the call.  Tony Gray 12/14/00
!                                                            
! Parameters Sent                                            
! R  = DIMENSION OF RANDOM BETA SPACE  H = 1 ... R           
! P  = DIMENSION OF FIXED ALPHA SPACE  L = 1 ... P           
! MAXCOL = max number of columns of variables to read        
! IDIND    = indicates which of the columns of TEMP          
!            to use for the ID variable                      
! YIND     = indicates which of the columns of TEMP          
!            to use for the dependent variable               
! XIND (R) = indicates which of the columns of TEMP          
!            to use in the analysis                          
! WIND (P) = indicates which of the columns of TEMP          
!            to use in the analysis                          
! WTIND    = indicates which of the columns of TEMP          
!            to use for the weighting variable               
!            if supplied.  may be left out by caller         
! OFIND    = indicates which of the columns of TEMP          
!            to use for the offset/censor    variable        
! MISS     = 1  then check for missing values                
! YMISS    = missing value code for Y variable               
! XMISS(R) = missing value codes for R X variables           
! WMISS(P) = missing value codes for P W variables           
!                                                            
! Parameters Returned                                        
! N    = TOTAL NUMBER OF SUBJECTS                            
! NTOT = TOTAL NUMBER OF OBSERVATIONS                        
! MAXK = MAX NUMBER OF OBS PER SUBJECT                       
!                                                            
! IDNI   = VECTOR (2 * MAXN) contains IDs and number of      
!          observations per ID                               
! FOR ORDER FOR EACH SUBJECT'S DATA ==>                      
! ID = SCALAR - SUBJECT'S ID                                 
! NI = SCALAR - NOBS FOR SUBJECT                             
!                                                            
! ALLDAT = VECTOR ((MAXNI*(R+P+1)) * MAXN) ALL DATA          
!          or equals (MAXNI*(R+P+2)) * MAXN) if iofs=1       
! FOR ORDER FOR EACH SUBJECT'S DATA ==>                      
! Y  = VECTOR - NI DEPENDENT VARIABLES                       
! X  = VECTOR - RNI (R*NI) DESIGN FOR RANDOM EFFECTS         
! W  = VECTOR - PNI (P*NI) FIXED COVARIATES                  
! C  = VECTOR - NI (optional) OFFSET/CENSOR CODES            
!                                                            
! Work Vectors                                               
! TEMP  = VECTOR (R+P+2) holds data temporarily              
!         or equals (R+P+3) if iofs=1                        
! XTEMP = VECTOR (R) HOLDS X OF EACH SUBJ TEMPORARILY        
! WTEMP = VECTOR (P) HOLDS W OF EACH SUBJ TEMPORARILY        
!                                                            
! *****************************************************************

SUBROUTINE READAT(FILEDAT,N,NTOT,MAXK,MAXCOL,R,P,ALLDAT, &
      IDNI,IDIND,YIND,XIND,WIND,MISS,YMISS,XMISS,WMISS, &
      IWT,WTIND,IOFS,OFIND)  ! Last 4 are optional
  
   IMPLICIT REAL*8 (A-H,O-Z)

   INTEGER H,R,P,XIND(*),WIND(*),IDIND,YIND,NTOT, &
           IDOLD,IDTEMP
   INTEGER ,POINTER ::IDNI(:)
   INTEGER PASS
   INTEGER, INTENT(in), OPTIONAL:: WTIND,IWT,IOFS,OFIND
   REAL (KIND = 8),POINTER ::ALLDAT(:)
   REAL*8  TEMP(MAXCOL),XTEMP(R),WTEMP(P),XMISS(R),WMISS(P)
   LOGICAL FIRST
   CHARACTER (len=*), INTENT(IN):: FILEDAT

   ! INITIALIZE 
   
   DO PASS = 1,2
   
   IF (PASS .EQ. 2) THEN
      ! write(6,"(' TG Allocating alldat(',I6,'), IDNI(',I6,')')")
      !+         ICOUNT,2*N
      ALLOCATE (ALLDAT(2*ICOUNT + 1))
      ALLDAT = 0.0D0
      ALLOCATE (IDNI(2*N))
      IDNI = 0
   ENDIF
   
   I     = 1
   K     = 1
   MAXK  = 0
   ICOUNT= 0
   NTOT  = 0
   FIRST  = .TRUE.

   ! READ IN DATA UNTIL END

   OPEN(1,ACTION='READ',FILE=FILEDAT)
   
   DO   ! loop forever
   
      READ(1,*,END=1999)(TEMP(H),H=1,MAXCOL)
    
      IDTEMP = INT(TEMP(IDIND))
      YTEMP  = TEMP(YIND)
      
      IF (PRESENT(IWT)) THEN
         IF (IWT  .EQ. 1) WTTEMP  = TEMP(WTIND)
      END IF
      IF (PRESENT(IOFS)) THEN
         IF (IOFS .EQ. 1) OFTEMP  = TEMP(OFIND)
      END IF
      
      ! Load up XTEMP with input col's indicated by XIND
      DO H = 1,R
         IN = XIND(H)
         XTEMP(H) = TEMP(IN)
      END DO

      ! Load up WTEMP with input col's indicated by WIND
      DO L = 1,P
         IN = WIND(L)
         WTEMP(L) = TEMP(IN)
      END DO

      ! Check for missing values if MISS=1
      MISSY = 0
      IF (MISS .EQ. 1) THEN
         IF (FP_EQUAL(YTEMP , YMISS)) THEN
            MISSY = 1
         END IF
            
         DO H = 1,R
            IF (FP_EQUAL(XTEMP(H) , XMISS(H))) THEN
               MISSY = 1
            END IF
         END DO
         DO L = 1,P
            IF (FP_EQUAL(WTEMP(L) , WMISS(L))) THEN
               MISSY = 1
            END IF
         END DO
      ELSE
      ENDIF
      IF (MISSY .NE. 0) THEN
         CYCLE  ! give up on current value and go read next one
      END IF

      ! QUERY FOR NEW ID AND SET PARAMETERS ACCORDINGLY

      IF (.NOT. FIRST) THEN 
         IF (R .GE. 1 .AND. IDTEMP .EQ. IDOLD) THEN
            K     = K+1
         ELSE
            IF (PRESENT(IWT)) THEN
               IF (IWT .EQ. 1) THEN
                  ICOUNT = ICOUNT+1
                  IF (PASS .EQ. 2) ALLDAT(ICOUNT) = WTOLD
                  WTOLD  = WTTEMP
               ENDIF
            END IF
            
            IC2 = 2*I
            IF (PASS .EQ. 2) THEN
               IDNI(IC2-1) = IDOLD
               IDNI(IC2) = K
            ENDIF
            NTOT = NTOT+K
            IF (K .GT. MAXK) MAXK = K
            K     = 1
            I     = I+1
         ENDIF
      ELSE
         WTOLD  = WTTEMP
      ENDIF
      

      ! PUT TEMPORARY VALUES INTO DATA VECTORS AND MATRICES

      IDOLD = IDTEMP
      
      ICOUNT = ICOUNT+1
      IF (PASS .EQ. 2) ALLDAT(ICOUNT) = YTEMP
      FIRST  = .FALSE.

      DO H=1,R
         ICOUNT = ICOUNT+1
         IF (PASS .EQ. 2) ALLDAT(ICOUNT) = XTEMP(H)
      END DO

      DO L=1,P
         ICOUNT = ICOUNT+1
         IF (PASS .EQ. 2) ALLDAT(ICOUNT) = WTEMP(L)
      END DO

      IF (PRESENT(IOFS)) THEN
         ! Caller is using optional offset/censoring in column OFIND
         IF (IOFS.EQ.1) THEN
             ICOUNT = ICOUNT+1
             IF (PASS .EQ. 2) ALLDAT(ICOUNT) = OFTEMP
         ELSE
         ENDIF
      END IF

   END DO   ! loop back to read next line

   ! cleanup final entry
   1999 IF (PRESENT(IWT)) THEN
      IF (IWT .EQ. 1) THEN
         ICOUNT = ICOUNT+1
         IF (PASS .EQ. 2) ALLDAT(ICOUNT) = WTOLD
      ENDIF
   END IF
   IC2 = 2*I
   IF (PASS .EQ. 2) THEN
      IDNI(IC2-1) = IDOLD
      IDNI(IC2) = K
   ENDIF
   NTOT = NTOT+K
   IF (K .GT. MAXK) MAXK = K
   
   N = I
   CLOSE(1)
   END DO   ! two passes, one to get size, second to read data
   
   !WRITE(6,"('first 10 of alldat: ',10F15.6)")(ALLDAT(J),J=1,10)
   !WRITE(6,"('first 10 of idni: ',10I6)")(IDNI(J),J=1,10)
   
   ! RETURN FROM SUBROUTINE
   RETURN
END SUBROUTINE READAT
   
   
! *****************************************************************
!   SUBROUTINTE QUADP (B,B1,A,NQ1,NQ,NDIM,IUNIF,WA1,WA2,sca) **  
!                                                         
!   DEFINES QUADRATURE POINTS AND WEIGHTS                 
!                                                         
!   Parameters Sent                                       
!   NQ1    = Number of Quadrature Points per Dimension    
!   NDIM   = Number of Dimensions                         
!   IUNIF  = 0 for Gaussian Distribution                  
!          = 1 for Uniform  Distribution                  
!          = 2 for log gamma Distribution                 
!   SCA    = value of scale parameter (only for log gamma)
!                                                         
!   Parameters Returned                                   
!   NQ     = Total Number of Quadrature Points            
!   B1     = NQ1 quadrature points (1 dimension)          
!   B      = NDIM*(NQ1**NDIM) quadrature points           
!   A      = NQ1**NDIM quadrature weights                 
!                                                         
!   Work Vectors                                          
!   WA1    = NQ1 work vector                              
!   WA2    = NDIM*(NQ1**NDIM) work vector                 
!                                                         
! *****************************************************************

SUBROUTINE QUADP(B,B1,A,NQ1,NQ,NDIM,IUNIF,WA1,WA2,SCA)
   IMPLICIT REAL*8(A-H,O-Z),INTEGER*4(I-N)
   DOUBLE PRECISION B,B1,A,WA1,WA2
   REAL (KIND = 8), INTENT(in out),OPTIONAL:: SCA
   DIMENSION B(1),B1(1),A(1),WA1(1),WA2(1)

   ! number of quadrature points, quadrature nodes & weights  

   IF (NDIM .EQ. 0) then
      NQ   = 1                           
      A(1) = 1.0d0
   ELSE
      NQ    = NQ1
      DO i=1,NQ
         B(I)=0.0D0
         A(I)=0.0D0
      END DO

      EPSQ=1.0D-8
      CALL HERMIT(B,A,NQ,EPSQ)
      IF (iunif .eq. 1) then
         DO i = 1,NQ
            A(I) = 1.0d0 / DBLE(NQ)
         END DO
      ELSEIF (iunif .eq. 2) then 
         IF (PRESENT(SCA)) THEN
            !CALL QLGAM(B,A,nq,SCA)
         ENDIF
      ENDIF
   ENDIF

   IF (NDIM .GE. 2 .AND. NDIM .LE. 8) then

      NQ  = NQ1**NDIM
      DO J=1,NQ1
         B1(J)  = B(J)
         WA1(J) = A(J)
      END DO

      J3=0
      DO J1 = 1,NDIM
         NDIV = NQ1**(J1-1)
         M2   = NQ/NDIV
         DO J2 = 1,NDIV
            DO J  = 1,M2
               K    = ((J-1)/(NQ1**(NDIM-J1)))+1
               J3   = J3 + 1
               B(J3)   = B1(K)
               WA2(J3) = WA1(K)
            END DO
         END DO
      END DO

      SUMW=0.0D0
      DO J = 1,NQ
         aqt = 1.0d0
         DO j1 = 2,NDIM
            J2  = NQ*(j1-1)+J
            aqt = aqt * WA2(J2)
         END DO
         A(J)  = WA2(J) * AQT
         SUMW  = SUMW   + A(J)
      END DO
      DO J = 1,NQ
         A(J) = A(J)/SUMW
      END DO
   ELSEIF (NDIM .gt. 8) then
      write(6,*)' Program can only have 8 random effects'
   ELSE
   ENDIF

   RETURN
END SUBROUTINE QUADP

!*********************************************************************
!                                                                 
!                  **********************                        
!                  *  SUBROUTINE GRAMM  *                       
!                  **********************                      
!                                                             
!   OBTAIN THE GRAMIAN MATRIX OF GENERALIZED PRODUCTS OF ROW 
!   VECTORS OF A SPECIFIED MATRIX BY PRE- AND POSTMULTIPLYING A   
!   METRIC MATRIX (SYMMETRIC OR DIAGONAL, MS=1,2) BY THE MATRIX  
!   AND ITS TRANSPOSE                                           
!                                                              
!   CALL GRAMM (A,B,C,M,N,MSB,W)                              
!                                                            
!   A .......... INPUT MATRIX, M BY N, GENERAL RECTANGULAR (MSA=0)
!   B .......... INPUT MATRIX, N BY N, METRIC MATRIX, SYMMETRIC OR
!                DIAGONAL                                        
!   C .......... OUTPUT MATRIX, M BY M, GENERALIZED GRAMIAN,    
!                SYMMETRIC (MSC=1)                             
!   M .......... NUMBER OF ROWS IN A                          
!   N .......... NUMBER OF COLUMNS IN A, EQUAL TO NUMBER OF ROWS  
!                AND NUMBER OF COLUMNS IN B                      
!   MSB ........ STORAGE MODE OF B                              
!   W .......... WORKING VECTOR OF LENGTH N                    
!                                                             
!     Storage modes: 0 SQUARE, 1 PACKED SYMMETRIC, 2 DIAGONAL, 
!                    3 PACKED LOWER TRIANGLE
!                    4 PACKED UPPER TRIANGLE
!
!*********************************************************************
SUBROUTINE GRAMM(A,B,C,M,N,MSB,W)                            
   DOUBLE PRECISION X,A,B,C,W                                  
   DIMENSION A(M,N),B(1),C(1),W(1)                            
   
   IF (MSB .EQ. 0) THEN
      CALL POST_ERROR('ERROR MESSAGE FROM SUBROUTINE GRAMM: METRIC IS GENERAL RECTANGULAR')
   
   ELSEIF (MSB .EQ. 1) THEN
      KC=0                                                          
      DO I=1,M                                                 
         KK=0                                                        
         DO K=1,N                                               
            X=0.0D0                                                   
            IC=KK                                                    
            DO L=1,K                                            
               IC=IC+1                                                
               X=X+A(I,L)*B(IC)                                      
            END DO
            IF (K .NE. N) THEN
               J=K+1                                               
               IC=IC+K                                            
               DO L=J,N                                      
                  X=X+A(I,L)*B(IC)                                 
                  IC=IC+L  
               END DO
            ENDIF
            W(K)=X 
            KK=KK+K 
         END DO
         
         DO K=1,I                                                   
            KC=KC+1                                                       
            X=0.0D0                                                     
            DO L=1,N                                               
               X=X+A(K,L)*W(L)                                           
            END DO
            C(KC)=X                                                  
         END DO
      END DO
   
   ELSE   ! MSB .GE. 2
      IC=0                                                 
      DO I=1,M                                        
         DO J=1,I                                       
            X=0.0D0                                           
            DO K=1,N                                     
               X=X+A(I,K)*A(J,K)*B(K)                          
            END DO
            IC=IC+1                                        
            C(IC)=X                                       
         END DO
      END DO
   ENDIF
   RETURN                                   
END SUBROUTINE GRAMM                                   

! SUBRUTINE TO READ LEVEL 3 DATA

 SUBROUTINE READLV3(FILEDAT,NOBS,NVAR,R3,ID3IND,IDIND,X3IND,N3,NI3,X3)
  INTEGER NTOT,MAXCOL,R3,ID3IND,X3IND(*)
  INTEGER,INTENT(OUT):: N3
  INTEGER,POINTER:: NI3(:)
  REAL(KIND=8),POINTER::X3(:,:)
  CHARACTER (len=*), INTENT(IN):: FILEDAT

  INTEGER I,J,IC1,IC2
  INTEGER,ALLOCATABLE:: IDS3(:),IDS2(:), ID3(:),NI3A(:),TEMP(:),TEMP1(:)
  REAL,ALLOCATABLE:: DATX(:,:)


  OPEN(UNIT=1,FILE=FILEDAT)
  ALLOCATE(datx(nobs,nvar))           ! read the file until the end
  READ(1,*) ((datx(i,j),j=1,nvar),i=1,nobs)
  CLOSE(1) 
  
 ALLOCATE(IDs3(nobs))
 IDs3(:)=INT(datx(:,id3ind))        ! vector of Level-3 ids
 ALLOCATE(temp(nobs))
 call unique(IDs3,nobs,temp)
 i=COUNT(temp>0)
 ALLOCATE(ID3(i))
 ID3(1:i)=temp(1:i)                 ! vector of unique level-3 ids
 DEALLOCATE(TEMP)

 ALLOCATE(x3(nobs,r3))
 if (r3>0) then
    do i=1,r3
       x3(:,i)=datx(:,x3ind(i))     ! Level-3 random effects
    end do
 end if

 if (r3==0) THEN
   n3=nobs
   ALLOCATE(ni3a(n3))
   ni3a(1:n3)=1
 else if (r3>0) then
   n3=COUNT(id3>0)                    ! number of level-3 clusters  
   ALLOCATE(ni3a(n3))
    DO i=1,n3
       ni3a(i)=COUNT(IDs3==id3(i))    ! n of each level-3 cluster
    END DO
 END if

 allocate(ids2(nobs))
 IDs2(:)=INT(datx(:,idind))        ! vector of Level-2 ids
 
! n level-2 units within each level-3 cluster
ALLOCATE(ni3(n3))
ni3(:)=0
ic1=0
ic2=0
do i=1,n3
   ic1=ic2+1
   ic2=ic1+ni3a(i)-1
   ALLOCATE(temp(ni3a(i)))
   ALLOCATE(temp1(ni3a(i)))
   temp(:)=0
   temp1(1:ni3a(i))=ids2(ic1:ic2)
   call unique(temp1,ni3a(i),temp)
   ni3(i)=COUNT(temp>0)
   DEALLOCATE(temp)
   DEALLOCATE(temp1)
end do

END SUBROUTINE READLV3



!*********************************************************************
!                                                                       
!                        **********************                        
!                        * FUNCTION ERROR_MESSAGE
!                        **********************                      
!                                                                   
! Either records an error message in the stack for retrieval by
! the DLL caller, or retrieves the message depending on the value
! of the PUSHPOP flag.  To register a message from FORTRAN code, 
! call with a string that contains the message.
! To pull a message out of the queue, supply a CHARACTER variable
! for the second string, and put a nonzero value in the PUSHPOP
! flag.  The function returns 0 if there is are no messages in the
! queue, 1 if there is one or more (including the one being pulled
! out).  
!                                                                   
!*********************************************************************
SUBROUTINE POST_ERROR(ErrorString)
   CHARACTER (len=*), INTENT(IN):: ErrorString
   IF(NEXT_FREE .LT. LBOUND(ERROR_STACK,1)) &
      NEXT_FREE = LBOUND(ERROR_STACK,1)
   IF(NEXT_FREE .GT. UBOUND(ERROR_STACK,1)) &
      NEXT_FREE = LBOUND(ERROR_STACK,1)
   IF(NEXT_RETRIEVED .LT. LBOUND(ERROR_STACK,1)) &
      NEXT_RETRIEVED = LBOUND(ERROR_STACK,1)
   IF(NEXT_RETRIEVED .GT. UBOUND(ERROR_STACK,1)) &
      NEXT_RETRIEVED = LBOUND(ERROR_STACK,1)

   ! Write the string to unit 6, std out or file for dll
   WRITE(6,"(/,1X,A)") ErrorString
   
   ! Store a string in the stack, and print it to unit 6
   ERROR_STACK(NEXT_FREE)(1:MIN(LEN(ErrorString),80)) = &
      ErrorString(1:MIN(LEN(ErrorString),80)) 
   NEXT_FREE = NEXT_FREE + 1
   IF(NEXT_FREE .GT. UBOUND(ERROR_STACK,1)) &
      NEXT_FREE = LBOUND(ERROR_STACK,1)
   IF(NEXT_FREE .EQ. NEXT_RETRIEVED) THEN
      ! We have overrun the buffer, possibly because
      ! the caller is not retrieving messages.  The oldest message
      ! has been overwritten, so bump up the queue tail.
      NEXT_RETRIEVED = NEXT_RETRIEVED +1
   ENDIF
END SUBROUTINE POST_ERROR
   
INTEGER FUNCTION RETRIEVE_ERROR(ErrorString)

   CHARACTER (len=*), INTENT(OUT):: ErrorString
   
   IF(NEXT_FREE .LT. LBOUND(ERROR_STACK,1)) &
      NEXT_FREE = LBOUND(ERROR_STACK,1)
   IF(NEXT_FREE .GT. UBOUND(ERROR_STACK,1)) &
      NEXT_FREE = LBOUND(ERROR_STACK,1)
   IF(NEXT_RETRIEVED .LT. LBOUND(ERROR_STACK,1)) &
      NEXT_RETRIEVED = LBOUND(ERROR_STACK,1)
   IF(NEXT_RETRIEVED .GT. UBOUND(ERROR_STACK,1)) &
      NEXT_RETRIEVED = LBOUND(ERROR_STACK,1)
   
   ! Retrieve the next string, if any
   IF(NEXT_RETRIEVED .NE. NEXT_FREE) THEN
      RETRIEVE_ERROR = 1
      ErrorString(1:MIN(LEN(ErrorString),LEN(ERROR_STACK(NEXT_RETRIEVED)))) &
         = ERROR_STACK(NEXT_RETRIEVED)(1:MIN(LEN(ErrorString),LEN(ERROR_STACK(NEXT_RETRIEVED))))
      NEXT_RETRIEVED = NEXT_RETRIEVED + 1
   ELSE
      RETRIEVE_ERROR = 0
   ENDIF
   
   RETURN
   
END FUNCTION RETRIEVE_ERROR
   
END MODULE MixLib
