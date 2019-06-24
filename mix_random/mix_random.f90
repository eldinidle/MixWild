PROGRAM MIX_random
    implicit none
    INTEGER :: I,nsubj,ndatasets,r,k,myseed,nonpos,mysize,j,ns,numre,numre2,mytotal,myid
    INTEGER, ALLOCATABLE :: tempseed(:)
    REAL(KIND=8) :: random_normal
    REAL(KIND=8),ALLOCATABLE:: chols(:,:),means(:,:),temp(:),myrand(:,:),temp2(:),cur_val(:)
    CHARACTER(LEN=80) :: fileebmeans,fileout

    OPEN(1, FILE='mix_random.def')
    READ(1,*)FILEebmeans
    READ(1,*)FILEout
    READ(1,*) Nsubj, r, ns, ndatasets, myseed
    close(1)

    numre = r+ns
    numre2 = numre*(numre+1)/2
    mytotal = numre + numre2
    allocate(chols(nsubj,numre2))
    allocate(means(nsubj,numre))
    allocate(temp(mytotal))
    allocate(temp2(numre2))
    allocate(cur_val(numre))

    open(1, file=fileebmeans)

    do k=1,nsubj
        read(1, *) myid,(temp(i),i=1,mytotal)
        means(k,:) = temp(1:numre)
        call chsky(temp(numre+1:mytotal),temp2,numre,nonpos)
        if(nonpos .eq. 0) write(*,*) "Problem with Cholesky for subject ",k
        chols(k,:) = temp2
        write(*,'(100f10.3)') (means(k,j), j=1,numre)
        write(*,'(100f10.3)') (chols(k,j), j=1,numre2)
    end do
    close(1)
    write(*,*) "Choleskys and means successfully saved!"
    call random_seed(size=mysize)
    allocate(tempseed(mysize))
    tempseed(:) = myseed
    call random_seed(put=tempseed)
    allocate(myrand(nsubj,numre))
    write(*,*) "Running seed number ", tempseed

    open(1, file=fileout)
    write(1,*) Nsubj, r, ns, ndatasets
    do i=1,ndatasets
        write(1,'(i4)',advance='no') i
        do j=1,nsubj
            do k=1,numre
                myrand(j,k) = random_normal()
            end do
        end do
        do k=1,nsubj
            cur_val = 0
            call mpytr(myrand(k,:),chols(k,:),cur_val,numre,1,3,numre)
            cur_val = cur_val + means(k,:)
            write(1,'(100f10.3)',advance='no') (cur_val(j),j=1,numre)
        end do
        write(1,*)
    end do

!   CALL MPYTR (A,B,C,MA,NA,MSB,NB)
!
!   A .......... INPUT MATRIX, MA BY NA, TRANSPOSED FIRST FACTOR
!                IN MULTIPLICATION, GENERAL RECTANGULAR (MSA=0)
!   B .......... INPUT MATRIX, MA BY NB, SECOND FACTOR IN MULTI-
!                PLICATION
!   C .......... OUTPUT MATRIX, NA BY NB, RESULT OF MULTIPLICATION
!                GENERAL RECTANGULAR (MSC=0)


end program mix_random

FUNCTION random_normal() RESULT(fn_val)

! Adapted from the following Fortran 77 code
!      ALGORITHM 712, COLLECTED ALGORITHMS FROM ACM.
!      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
!      VOL. 18, NO. 4, DECEMBER, 1992, PP. 434-435.

!  The function random_normal() returns a normally distributed pseudo-random
!  number with zero mean and unit variance.

!  The algorithm uses the ratio of uniforms method of A.J. Kinderman
!  and J.F. Monahan augmented with quadratic bounding curves.

REAL(kind=8) :: fn_val

!     Local variables
REAL(kind=8)     :: s = 0.449871, t = -0.386595, a = 0.19600, b = 0.25472,           &
            r1 = 0.27597, r2 = 0.27846, u, v, x, y, q

!     Generate P = (u,v) uniform in rectangle enclosing acceptance region

DO
  CALL RANDOM_NUMBER(u)
  CALL RANDOM_NUMBER(v)
  v = 1.7156 * (v - half)

!     Evaluate the quadratic form
  x = u - s
  y = ABS(v) - t
  q = x**2 + y*(a*y - b*x)

!     Accept P if inside inner ellipse
  IF (q < r1) EXIT
!     Reject P if outside outer ellipse
  IF (q > r2) CYCLE
!     Reject P if outside acceptance region
  IF (v**2 < -4.0*LOG(u)*u**2) EXIT
END DO

!     Return ratio of P's coordinates as the normal deviate
fn_val = v/u
RETURN

END FUNCTION random_normal


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
