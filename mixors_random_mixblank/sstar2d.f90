! *****************************************************************
! ****  SUBROUTINE GETDNPLUS(DNPLUS,N,N4)                      ****
! ****                                                         ****
! ****  BUILDS D_n+ which is THE N*N x N*(N+1)/2               ****
! ****  INVERTED DUPLICATION MATRIX FROM MAGNUS (1988) page 56 ****
! ****  WITH THE PROPERTY D_n+ vec(A) = v(A)                   ****
! ****                                                         ****
! ****  where A is n x n and symmetric                         ****
! ****                                                         ****
! ****  Parameters Sent                                        ****
! ****     N = number of rows ( = # of columns) of A           ****
! ****     N4 = n*n*n*(n+1)/2                                  ****
! ****  DNPLUS = N4  INVERTED DUPLICATION MATRIX               ****
! *****************************************************************
SUBROUTINE getdnplus(dnplus,n,nc2,n2)
    implicit none
    integer,intent(in)::n,n2,nc2
    integer::i,j,row,col1,col2
    real(kind=8)::dnplus(nc2,n2)
    
    dnplus = 0.0
    do j=1,n
        row = n*(j-1)+j-j*(j-1)/2
        col1 = j+n*(j-1)
        dnplus(row,col1) = 1    
        do i=j+1,n
            row = n*(j-1)+i-j*(j-1)/2
            col1 = i+n*(j-1)
            col2 = j+n*(i-1)
            dnplus(row,col1) = .5
            dnplus(row,col2) = .5
        end do
    end do
END SUBROUTINE getdnplus

! *****************************************************************
! ****  SUBROUTINE GETUNP(UNP,N,N4)                            ****
! ****                                                         ****
! ****  BUILDS U_n' which is THE revised N*(N+1)/2 x N*N       ****
! ****  ELIMINATION MATRIX FROM MAGNUS (1988) page 76          ****
! ****  WITH THE PROPERTY U_n' v(A) = vec A                    ****
! ****                                                         ****
! ****  where A is n x n and upper triangular                  ****
! ****                                                         ****
! ****  Parameters Sent                                        ****
! ****     N = number of rows ( = # of columns) of A           ****
! ****     N4 = n*n*n*(n+1)/2                                  ****
! ****  UNP = N4  ELIMINATION MATRIX                           ****
! *****************************************************************
SUBROUTINE getUnp(unp,n,nc2,n2)
    implicit none
    integer,intent(in)::n,n2,nc2
    integer::i,j,row,col
    real(kind=8),intent(out)::unp(n2,nc2)
    
    unp = 0
    do i=1,n
        do j=i,n
            col = i+j*(j-1)/2
            row = i+n*(j-1)
            unp(row,col) = 1
        end do
    end do
END SUBROUTINE getunp

subroutine getInKS(InKS,Svech,n,nc2,n2)
    implicit none
    integer,intent(in)::n,n2,nc2
    real(kind=8)::temp(n,n)
    integer::num,sstart,i,j,k
    real(kind=8),intent(out)::InKS(n2,n2)
    real(kind=8),intent(in)::Svech(nc2)
    
    temp = 0
    i=1
    do k=1,n
        do j=k,n
            temp(j,k) = svech(i)
            i=i+1
        end do
    end do
!    write(*,*) ((temp(k,j),k=1,n),j=1,n)
        
    InKS = 0
    do num = 1, n
        sstart=(num-1)*n
        InKS(sstart+1:sstart+n,sstart+1:sstart+n) = temp
    end do
end subroutine getInKS

subroutine getSStar(Svech, n, nc2, sStar)
    implicit none
    integer,intent(in)::n,nc2
    integer::n2,i,j
    real(kind=8),intent(in)::Svech(nc2)
    real(kind=8),intent(out)::sStar(nc2,nc2)
    real(kind=8),allocatable::dnplus(:,:),unp(:,:),InKS(:,:),work(:,:)

    n2 = n*n
    allocate(dnplus(nc2,n2))
    allocate(unp(n2,nc2))
    allocate(InKS(n2,n2))
    allocate(work(nc2,n2))

    call getdnplus(dnplus,n,nc2,n2)
!do j=1,nc2
!    write(69,'(6f5.2)') (dnplus(j,i), i=1,n2)
!end do
    call getInKS(InKS,Svech,n,nc2,n2)
!    write(*,*)
!do j=1,n2
!    write(69,'(6f5.2)') (inKS(j,i), i=1,n2)
!end do
    work = matmul(dnplus,InKS)
!do j=1,nc2
!    write(69,'(6f5.2)') (work(j,i), i=1,n2)
!end do
    call getunp(unp,n,nc2,n2)
!do j=1,n2
!    write(69,'(6f5.2)') (unp(j,i), i=1,nc2)
!end do
!write(*,*)
    sstar = matmul(work,unp)
!do j=1,nc2
!    write(69,'(6f5.2)') (sstar(j,i), i=1,nc2)
!end do
    deallocate(dnplus)
    deallocate(unp)
    deallocate(InKS)
    deallocate(work)
end subroutine getSStar


subroutine inverse(a,c,n)
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed 
! during the calculation
!===========================================================
    implicit none 
    integer n
    double precision a(n,n), c(n,n)
    double precision L(n,n), U(n,n), b(n), d(n), x(n)
    double precision coeff
    integer i, j, k
    
    ! step 0: initialization for matrices L and U and b
    ! Fortran 90/95 aloows such operations on matrices
    L=0.0
    U=0.0
    b=0.0
    
    ! step 1: forward elimination
    do k=1, n-1
        do i=k+1,n
            coeff=a(i,k)/a(k,k)
            L(i,k) = coeff
            do j=k+1,n
                a(i,j) = a(i,j)-coeff*a(k,j)
            end do
       end do
    end do
    
    ! Step 2: prepare L and U matrices 
    ! L matrix is a matrix of the elimination coefficient
    ! + the diagonal elements are 1.0
    do i=1,n
      L(i,i) = 1.0
    end do
    ! U matrix is the upper triangular part of A
    do j=1,n
      do i=1,j
        U(i,j) = a(i,j)
      end do
    end do
    
    ! Step 3: compute columns of the inverse matrix C
    do k=1,n
      b(k)=1.0
      d(1) = b(1)
    ! Step 3a: Solve Ld=b using the forward substitution
      do i=2,n
        d(i)=b(i)
        do j=1,i-1
          d(i) = d(i) - L(i,j)*d(j)
        end do
      end do
    ! Step 3b: Solve Ux=d using the back substitution
      x(n)=d(n)/U(n,n)
      do i = n-1,1,-1
        x(i) = d(i)
        do j=n,i+1,-1
          x(i)=x(i)-U(i,j)*x(j)
        end do
        x(i) = x(i)/u(i,i)
      end do
    ! Step 3c: fill the solutions x(n) into column k of C
      do i=1,n
        c(i,k) = x(i)
      end do
      b(k)=0.0
    end do
end subroutine inverse

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

! ----------------------------------------------------------------------------
!
! SUBROUTINE PHIder
! Calculate the derivative of the pdf for
! various distributions:
! NORMAL, LOGISTIC, Complementary Log-Log, OR Log-Log
!    = 0     = 1              = 2         =3
! ----------------------------------------------------------------------------
!@    phider - calculate the derivative of the pdf for Normal, Logistic, Compl. Log-Log, or Log-Log

real(kind=8) function PHIder(Z, nfn)
    implicit none
    integer, intent(IN) :: nfn
    real(kind=8), intent(in) :: Z

    real(kind=8) :: phifn,pi
    parameter (PI=3.14159265358979323846264338327950D+00)
  ! Z = no. of standard deviations from the mean.
    select case(nfn)

        case(0)  ! Normal distribution (dnorm)
            if(ABS(Z) .gt. 37.0D0) then
                PHIder = 0d0
            else
                PHIder =-Z*(DEXP(-Z*Z*.5d0)) / SQRT(2*pi)
            end if

        case(1)  ! Logistic distribution (dlogis)
            PHIFN = 1.0D0 /(1.0D0+DEXP(0.0D0-Z))
            phider = phifn*(1-phifn)*(1-2*phifn)

        case(2)  ! Complementary Log-Log distribution
            PHIder = DEXP(Z-dexp(z))*(1-dexp(z))

        case(3)  ! Log-Log distribution
        ! revise
        PHIder = DEXP(-Z-dexp(-z))*(-1+dexp(-z))
        
  end select
  return

end function PHIder

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

!Calculates ratio of Standard Normal pdfs PHI(newB)/PHI(origB)
subroutine GET_PHI_RATIO(newB,origB, phiRatio)
    implicit none
    REAL(KIND=8)::newB,origB,phiratio
    PHIRATIO = exp((origB*origB-newB*newB)/2)
END subroutine GET_PHI_RATIO

subroutine descripc(a,nobs,cols,m,s,minv,maxv)
    implicit none
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
