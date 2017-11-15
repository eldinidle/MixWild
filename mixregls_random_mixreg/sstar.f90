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
SUBROUTINE getdnplus(dnplus,n,n4)
    implicit none
    integer,intent(in)::n,n4
    integer::i,j,row,col1,col2,n2
    real(kind=8)::dnplus(n4)
    
    n2 = n*(n+1)/2
    dnplus = 0.0
    do j=1,n
        row = n*(j-1)+j-j*(j-1)/2
        col1 = j+n*(j-1)
        dnplus(row+n2*(col1-1)) = 1    
        do i=j+1,n
            row = n*(j-1)+i-j*(j-1)/2
            col1 = i+n*(j-1)
            col2 = j+n*(i-1)
            dnplus(row+n2*(col1-1)) = .5
            dnplus(row+n2*(col2-1)) = .5
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
SUBROUTINE getUnp(unp,n,n4)
    implicit none
    integer,intent(in)::n,n4
    integer::i,j,row,col
    real(kind=8)::unp(n4)
    
    unp = 0
    do i=1,n
        do j=i,n
            col = i+j*(j-1)/2
            row = i+n*(j-1)
            unp(row+n*n*(col-1)) = 1
        end do
    end do
END SUBROUTINE getunp

!subroutine getInKS(InKS,vechS,n,nStar,n4)
!    implicit none
!    integer,intent(in)::n,n4,nStar
!    integer::i,j,n2,countInKS,countS
!    real(kind=8)::InKS(n4)
!    real(kind=8),intent(in)::vechS(nStar)
!    real(kind=8),allocatable::Iden(:)
!    
!    allocate(Iden(n*n))
!    call gend(Iden,1.0D0,n,3)
!    call kmpy(Iden,vechS,InKS,n,n,3,3,n)
!    deallocate(Iden)
!end subroutine getInKS
!    

subroutine getInKSprime(InKS,S,n,nStar,n4)
    implicit none
    integer,intent(in)::n,n4,nStar
    integer::row,n2,num,indent,sstart,ksstart
    real(kind=8)::InKS(n4)
    real(kind=8),intent(in)::S(nStar)
    
    InKS = 0
    n2 = n*n
    do row = 1, n2
        num = mod(row-1, n)+1
        indent = (row - 1) / n
        sstart = num*(num-1)/2
        ksstart = n2*(row-1)+indent*n
        InKS(ksstart+1:ksstart+num) = S(sstart+1:sstart+num)
    end do
end subroutine getInKSprime

subroutine getSStar(vechS, n, nStar, sStar)
    implicit none
    integer,intent(in)::n,nStar
    integer::n2,n2nStar,n4!,i
    real(kind=8),intent(in)::vechS(nStar)
    real(kind=8)::sStar(nStar*nStar)
    real(kind=8),allocatable::dnplus(:),unp(:),InKS(:),InKSprime(:),work(:)

    n2 = n*n
    n2nStar = n2 * nStar
    n4 = n2 * n2
    allocate(dnplus(n2nStar))
    allocate(unp(n2nStar))
    allocate(InKS(n4))
    allocate(InKSprime(n4))
    allocate(work(n2nStar))

    call getdnplus(dnplus,n,n2nStar)
!    write(*,'(6f4.1)') (dnplus(i), i=1,n2nStar)
    call getInKSprime(InKSprime,vechS,n,nStar,n4)
!    write(*,*)
!    write(*,'(9f4.1)') (inKS(i), i=1,n4)
    call trp(InKSprime,InKS,n2,n2)
    call mpym(dnplus,InKS,work,nStar,n2,0,0,n2)
!    write(*,*)
!    write(*,'(16f4.1)') (work(i), i=1,n2nStar)
    call getunp(unp,n,n2nStar)
!    write(*,*)
!    write(*,'(9f4.1)') (unp(i), i=1,n2nStar)
    call mpym(work,unp,sStar,nStar,n2,0,0,nStar)
    deallocate(dnplus)
    deallocate(unp)
    deallocate(InKS)
    deallocate(InKSprime)
    deallocate(work)
end subroutine getSStar

SUBROUTINE TRP(A,B,MA,NA)    
    implicit none
    integer::ma,na,i,j
    real(kind=8)::a(ma,na),b(na,ma)
    do i=1,na
        do j=1,ma
            b(i,j) = a(j, i)
        end do
    end do
end subroutine trp

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