module lsboth
    implicit none
    save
    INTEGER :: NOBS,NVAR,NQ,AQUAD,ID2IND,YIND,P,R,rr,S,ns,MISS,MAXK,NC2,MAXIT,NCENT,&
                PNINT,RNINT,SNINT,POLD,ROLD,SOLD,ncov,num0,discard0,nv,nvar2,mls,&
                rv,sv,pv,numloc,chol,npar,npar2,ndim,ndim2,nqwr1,nqwr0,stage2,nreps,nors,&
                pfixed,ptheta,pomega,pto,readcats,nvar3,maxj,multi2nd,myseed,sepfile,nvarsep,&
                nc2sep,id2indsep,nobssep,maxksep
    INTEGER,ALLOCATABLE :: XIND(:),UIND(:),WIND(:),IDNI(:,:),varIND(:),ids(:),var2ind(:),idnisep(:,:),allzeros(:)
    REAL(KIND=8) :: RIDGEIN,CONV,YMISS,SDLV,RCORR,errv,cutoff
    REAL(KIND=8),ALLOCATABLE:: Y(:),X(:,:),U(:,:),W(:,:),BETA(:),TAU(:),SPAR(:), tempsums(:,:),mychol(:),&
                               alpha(:),thetas(:,:),thetavs(:,:),var(:,:),varavg(:,:),simvals(:,:,:),icode(:),&
                               data2(:,:)
                               
    CHARACTER(LEN=16) :: YLABEL
    CHARACTER(LEN=4) :: HEAD(36)
    CHARACTER(LEN=16),ALLOCATABLE :: BLAB(:),ALAB(:),TLAB(:),var2label(:)
    character(len=12),allocatable::varlabel(:)
    CHARACTER(LEN=80) :: FILEDAT, FILEOUT, FILEDEF,fileprefix,filedat2,ebfile
    CHARACTER(LEN=24),ALLOCATABLE :: intlabel(:)

end module lsboth

module procedures
!implicit none

    contains
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

    !SUBROUTINE INVS(A,N,C,W,IER,VERBOSE)
    SUBROUTINE INVS(A,N,C,IER)
        implicit none
        real(kind=8),intent(inout)::a(n*(n+1)/2),c
        integer,intent(in)::n
        integer,intent(out)::ier
        integer::k,j,i,km1
            double precision::x,u,y,z,d,w(n)
        INTEGER DIAGMK                                                    
        INTEGER DIAG,DIAG2,ROWNO,ROWCOL                                   
        INTEGER COLNO
        !   INTEGER, INTENT(in out),OPTIONAL:: IER
        !   LOGICAL, INTENT(in),OPTIONAL:: VERBOSE

        !   IF (PRESENT(IER)) THEN
              IER = 0
        !   END IF

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
                 IER = 1
        ENDIF
    END SUBROUTINE INVS

    !                        **********************
    !                        *  SUBROUTINE   MPYM *
    !                        **********************
    !
    !         MULTIPLY TWO MATRICES
    !
    !         CALL MPYM(A,B,C,MA,NA,MSA,MSB,NB)
    !
    !         A .......... INPUT MATRIX, MA BY NA, FIRST FACTOR IN MULTI-
    !                      PLICATION
    !         B .......... INPUT MATRIX, NA BY NB, SECOND FACTOR IN MULTI-
    !                      PLICATION
    !         C .......... OUTPUT MATRIX, MA BY NB, RESULT OF MULTIPLI-
    !                      CATION, GENERAL RECTANGULAR (MS=0)
    !         MA ......... NUMBER OF ROWS IN A
    !         NA ......... NUMBER OF COLUMNS IN A, EQUAL TO NUMBER OF ROWS
    !                      IN B
    !         MSA ........ STORAGE MODE OF A
    !         MSB ........ STORAGE MODE OF B
    !         NB ......... NUMBER OF COLUMNS IN B
    !
          SUBROUTINE MPYM (A,B,C,MA,NA,MSA,MSB,NB)
        implicit none
    real(kind=8),intent(in)::a(:),b(:)
    real(kind=8),intent(out)::c(ma*nb)
    integer,intent(in)::ma,na,msa,msb,nb
    integer::k,j,i,ica,icb,icca,iccb,iback,jback,jj,mma,mmb,loop,ii,jnb,inb,incb,llb,lhb,ina,inca,lla,lha,jna,k1,k2,k3,m
        real(kind=8)::x
          MMA=MSA+1
          MMB=MSB+1
          ICA=0
          ICB=0
          ICCA=0
          ICCB=0
          LOOP=1
          IBACK=1
          JBACK=1
          II=0
          JJ=0
          IF(MSA-1)99,101,97
    99 IF(MSB-1)101,101,96
    97 IF(MSB-1)95,101,101
    96 IF(MSB-3)101,101,93
    93 JBACK=2
          GOTO 101
    95 IF(MSA-3)101,91,90
    90 LOOP=2
          GOTO 101
    91 LOOP=2
          IBACK=2
    101 GOTO(102,105),LOOP
    !
    !     DEFINITION OF PARAMETERS FOR MATRIX B
    !
    102 JJ=JJ+1
          GOTO(110,111),JBACK
    110 J=JJ
          GOTO 112
    111 J=NB-JJ+1
    112 GOTO(11,15,13,14,15),MMB
    11 INB=(J-1)*NA+1
          INCB=1
          LLB=1
          LHB=NA
          GOTO 20
    13 INB=J
          INCB=0
          LLB=J
          LHB=J
          GOTO 20
    14 INB=J*(J+1)/2
          INCB=J
          ICCB=1
          LLB=J
          LHB=NA
          GOTO 20
    15 INB=J*(J-1)/2+1
          INCB=1
          LLB=1
          IF(MSB-1)35,75,35
    35 LHB=J
          GOTO 20
    75 LHB=NA
    20 JNB=INB
          GOTO(105,106),LOOP
    !
    !     DEFINITION OF PARAMETERS FOR MATRIX A
    !
    105 II=II+1
          GOTO(120,121),IBACK
    120 I=II
          GOTO 122
    121 I=MA-II+1
    122 GOTO(12,18,17,18,16),MMA
    12 INA=I
          INCA=MA
          LLA=1
          LHA=NA
          GOTO 21
    17 INA=I
          INCA=0
          LLA=I
          LHA=I
          GOTO 21
    18 INA=I*(I-1)/2+1
          INCA=1
          LLA=1
          IF(MSA-1)73,74,73
    73 LHA=I
          GOTO 21
    74 LHA=NA
          ICCA=0
          ICA=0
          GOTO 21
    16 INA=I*(I+1)/2
          INCA=I
          ICCA=1
          ICA=0
          LLA=I
          LHA=NA
    21 JNA=INA
          GOTO(106,102),LOOP
    106 K=I+(J-1)*MA
          X=0.0D0
    !
    !     SYNCHRONIZATION OF PARAMETERS
    !
          IF(LLA-LLB)23,25,40
    23 IF(LHA-LLB)100,24,24
    40 IF(LHB-LLA)100,26,26
    25 K1=LLA
          GOTO 41
    24 K1=LLB
          K3=K1-1
          DO 44 M=LLA,K3
          IF(MSA-1)244,45,244
    45 IF(M-I)244,46,244
    46 INCA=I
          ICCA=1
    244 INA=INA+INCA+ICA
    44 ICA=ICA+ICCA
          GOTO 41
    26 K1=LLA
          K3=K1-1
          DO 29 M=LLB,K3
          IF(MSB-1)129,28,129
    28 IF(M-J)129,55,129
    55 INCB=J
          ICCB=1
    129 INB=INB+INCB+ICB
    29 ICB=ICB+ICCB
    41 IF(LHA-LHB)19,22,22
    19 K2=LHA
          GOTO 27
    22 K2=LHB
    !
    !     VECTOR MULTIPLICATION AND RESETTING OF PARAMETERS
    !
    27 DO 30 M=K1,K2
          X=X+A(INA)*B(INB)
          IF(MSA-1)60,61,60
    61 IF(M-I)60,62,60
    62 INCA=I
          ICCA=1
    60 IF(MSB-1)65,63,65
    63 IF(M-J)65,64,65
    64 INCB=J
          ICCB=1
    65 INA=INA+INCA+ICA
          ICA=ICA+ICCA
          INB=INB+INCB+ICB
    30 ICB=ICB+ICCB
          IF(MSB-1)83,82,83
    82 INCB=1
          ICCB=0
    83 INB=JNB
          ICB=0
          IF(MSA-1)130,131,130
    131 INCA=1
          ICCA=0
    130 INA=JNA
          ICA=0
    100 C(K)=X
          GOTO(140,141),LOOP
    140 IF(II-MA)105,143,6666
    143 IF(JJ-NB)146,6666,6666
    141 IF(JJ-NB)102,145,6666
    145 IF(II-MA)147,6666,6666
    146 II=0
          GOTO 102
    147 JJ=0
          GOTO 105
    6666 RETURN
          END SUBROUTINE MPYM
                                                                      
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
        implicit none
    real(kind=8),intent(in)::a(ma*na),b(:)
    real(kind=8),intent(out)::c(na*nb)
    integer,intent(in)::ma,na,msb,nb
    integer::k,j,i,icb,iccb,jnb,inb,incb,llb,lhb,ina,m
        real(kind=8)::x
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
        implicit none
    real(kind=8),intent(in)::a(n*(n+1)/2),x(n),c
    real(kind=8),intent(out)::b(n*(n+1)/2)
    integer,intent(in)::n
    integer::ic,j,i
    IC=0                                                              
    DO I=1,N                                                      
          DO J=1,I                                                      
             IC=IC+1                                                           
             B(IC)=A(IC)+C*X(I)*X(J)                                           
          END DO
    END DO
    END SUBROUTINE GRMCV

    SUBROUTINE GRMCV_inout(A,X,C,N)                                       
        implicit none
    real(kind=8),intent(in)::x(n),c
    real(kind=8),intent(inout)::a(n*(n+1)/2)
    integer,intent(in)::n
    integer::ic,j,i
    IC=0                                                              
    DO I=1,N                                                      
          DO J=1,I                                                      
             IC=IC+1                                                           
             A(IC)=A(IC)+C*X(I)*X(J)                                           
          END DO
    END DO
    END SUBROUTINE GRMCV_inout
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
        implicit none
    real(kind=8),intent(in)::a(m*n)
    real(kind=8),intent(out)::c(m*(m+1)/2)
    integer,intent(in)::m,n
    integer::k,ic,j,i
        real(kind=8)::x

    IC=0                                                      
    DO I=1,M                                             
          DO J=1,I                                            
             X=0.0D0                                                
             DO K=1,N                                          
                X=X+A((I-1)*n+K)*A((J-1)*n+K)                                    
             END DO
             IC=IC+1                                             
             C(IC)=X                                            
          END DO
    END DO
    
    RETURN                                        
    END SUBROUTINE GRAM

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
        implicit none
    real(kind=8),intent(in)::a(:),b(mb*nb),c(:)
    real(kind=8),intent(out)::d(:)
    integer,intent(in)::mb,nb
    integer::k,l,j,i
    
    K = 0                                                           
    L = (NB*(NB+1))/2                                              
    
    DO J=1,L                                            
          D(J) = A(J)                                         
    END DO
    
    DO I=1,MB                                        
          DO J=1,NB                                       
             L = L + 1                                        
             D(L) = B((j-1)*mb+i)                                   
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
        implicit none
    real(kind=8),intent(in)::a(m*na),b(m*nb)
    real(kind=8),intent(out)::c(m*(na+nb))
    integer,intent(in)::m,na,nb
    integer::k,kk,l,jz
    
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
        implicit none
        integer,intent(in)::n
        real(kind=8),intent(in)::a(n*(n+1)/2)
        real(kind=8),intent(out)::b(n*(n+1)/2)
        integer,intent(out)::nonpos
        integer::n1,kc,ifir,i,j,ic,jc,ii,k
        real(kind=8)::x,y
    
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
        integer::n2!,i,k
        real(kind=8),intent(in)::Svech(nc2)
        real(kind=8),intent(out)::sStar(nc2,nc2)
        real(kind=8),allocatable::dnplus(:,:),unp(:,:),InKS(:,:),work(:,:)

        n2 = n*n
        allocate(dnplus(nc2,n2))
        allocate(unp(n2,nc2))
        allocate(InKS(n2,n2))
        allocate(work(nc2,n2))
        

        call getdnplus(dnplus,n,nc2,n2)
        call getInKS(InKS,Svech,n,nc2,n2)
        work = matmul(dnplus,InKS)
        call getunp(unp,n,nc2,n2)
        sstar = matmul(work,unp)
        deallocate(dnplus)
        deallocate(unp)
        deallocate(InKS)
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
        DOUBLE PRECISION A,B, Eps
        PARAMETER (EPS = .0000005) ! required closeness of comparands
        IF (ABS(B - A) .LE. (ABS(B+A)*EPS)) THEN
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

    !Pulls necessary Standard Deviation from Packed Variance matrix
    subroutine getSDev(n, length, Var, sdev)
        implicit none
        integer::n, length, i, j, place
        real(kind=8)::Var(length*(length+1)/2), sdev
        place = 1
        do i=1,length
            do j=1,i
                if(i==j .and. i==n) sdev = dsqrt(Var(place))
                place = place + 1
            end do
        end do
    end subroutine getSDev

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

    SUBROUTINE SUBMANDV()
            use lsboth, only: y,rcorr,nc2,idni,sdlv
            implicit none
            REAL(KIND=8) :: RTEMP, meanlv
            REAL(KIND=8),ALLOCATABLE:: IDMV(:,:),TEMPR(:),meanMV(:),stdMV(:)
            INTEGER :: IC,I,J

            ALLOCATE (IDMV(NC2,2))
            IDMV = 0.0D0
            IC=0
            DO I=1,NC2
                DO J=1,IDNI(I,2)
                     IC=IC+1
                     IDMV(I,1) = IDMV(I,1) + Y(IC)
                     IDMV(I,2) = IDMV(I,2) + (Y(IC))**2
                END DO
                IDMV(I,1) = IDMV(I,1)/DBLE(IDNI(I,2))
                IF (IDNI(I,2) == 1) THEN
                     IDMV(I,2) = (IDMV(I,2) - DBLE(IDNI(I,2))*(IDMV(I,1))**2) / (DBLE(IDNI(I,2)))
                ELSE IF (IDNI(I,2) > 1) THEN 
                     IDMV(I,2) = (IDMV(I,2) - DBLE(IDNI(I,2))*(IDMV(I,1))**2) / (DBLE(IDNI(I,2)) - 1.0D0)
                END IF
                IF (IDMV(I,2) .LT. 0.0D0) IDMV(I,2) = 0.0D0
            END DO
            
            MEANLV=SUM(DLOG(IDMV(:,2)+1.0D0))/DBLE(NC2)
                 ALLOCATE(TEMPR(NC2))
                 tempR(:)=0.0D0
                 tempR(:)=((DLOG(IDMV(:,2)+1.0D0))-MEANLV)**2
                 RTEMP = SUM(tempR)/DBLE(NC2-1)
                 SDLV=DSQRT(RTEMP)
                 DEALLOCATE(TEMPR)

            CALL CORR(IDMV(:,1),IDMV(:,2),NC2,RCORR)
            !WRITE(*,*)SDLV,RCORR

            ! MAKE Z VALUES
            ALLOCATE(meanMV(2))
            ALLOCATE(stdMV(2))
            call meanc(IDMV,NC2,2,meanMV,stdMV)
            IDMV(:,1) = (IDMV(:,1)-MEANMV(1)) / STDMV(1)
            IDMV(:,2) = (IDMV(:,2)-MEANMV(2)) / STDMV(2)
            DEALLOCATE(meanMV)
            DEALLOCATE(stdMV)

    END SUBROUTINE SUBMANDV 

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

    subroutine CORR(X,Y,N,R)
        !purpose:
        !to calculate the correlation between two vectors
        !
            INTEGER,INTENT(IN)::N                       ! number of rows of two vectors
            REAL(KIND=8),INTENT(IN),DIMENSION(n)::X,Y   ! input vectors
            REAL(KIND=8),INTENT(OUT):: R                ! correlation of X and Y
            REAL(KIND=8),DIMENSION(N)::XT,YT            ! temporary arrays
            REAL(KIND=8):: AX,AY,SXX,SYY,SXY

            ! calculate correlation
            AX = SUM(X)/DBLE(N)
            AY = SUM(Y)/DBLE(N)

            XT(:) = X(:)-AX
            YT(:) = Y(:)-AY

            SXX = DOT_PRODUCT(XT,XT)
            SYY = DOT_PRODUCT(YT,YT)
            SXY = DOT_PRODUCT(XT,YT)

            R = SXY/DSQRT(SXX*SYY)

    end subroutine CORR

    FUNCTION random_normal() RESULT(fn_val)
        implicit none

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
                r1 = 0.27597, r2 = 0.27846, u, v, x, y, q, half=.5

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
    END FUNCTION random_normal

    subroutine cholesky_sub(A,n,zeros)
        implicit none

        ! formal vars
        integer,intent(in) :: n      ! number of rows/cols in matrix
        real(kind=8),intent(inout)   :: A(n,n) ! matrix to be decomposed
        integer,intent(out) :: zeros

        ! local vars
        integer :: j      ! iteration counter

        ! begin loop
        zeros = 0
        A(1,1) = sqrt(A(1,1))
        if(.not. fp_equal(A(1,1), 0.0d0)) then
            A(2:n,1) = A(2:n,1) / A(1,1)
            chol: do j = 2,n

                ! perform diagonal component
                A(j,j) = sqrt(A(j,j) - dot_product(A(j,1:j-1),A(j,1:j-1)))
                ! perform off-diagonal component
                if (j < n) A(j+1:n,j) = (A(j+1:n,j) - matmul(A(j+1:n,1:j-1),A(j,1:j-1))) / &
            &           A(j,j)

            end do chol

                do j=1,n
                    A(1:j-1,j) = 0
                end do
        else
            zeros = 1
            A = 0
        end if
    end subroutine cholesky_sub
end module procedures