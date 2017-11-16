PROGRAM repeat_mixor
    implicit none
    INTEGER :: I,NOBS,NVAR,ID2IND,numRE,r,MISS,MAXK,NC2,j,&
                nvar2,pfixed,ptheta,pomega,pto,k,nvar3,diffk,&
                nsubj,ndatasets,nreps,n,maxj,maxits,numits,numgood,discard,ns
    INTEGER,ALLOCATABLE :: var2IND(:),ids(:),idni(:,:),jcodes(:),isgood(:)
    REAL(KIND=8) :: YMISS,zval,pval,phifn,logl,loglsd,cutoff
    REAL(KIND=8),ALLOCATABLE:: tempsums(:,:),temp3(:),&
                                betas(:,:),meanbetas(:),totalvarhats(:),varbetas(:),&
                               thetatemp(:),thetas(:,:),tempdata(:),liks(:),dev(:)
    CHARACTER(LEN=16),ALLOCATABLE :: var2label(:)
    CHARACTER(LEN=24),ALLOCATABLE :: intlabel(:)
    CHARACTER(LEN=100) :: FILEDAT,fileeb,filetempdat,FILEprefix,temp80
    character(len=160) :: fileout
    character(len=24) :: mystr
    logical :: fp_equal

    OPEN(1, FILE='repeat_mixor.def')
    read(1,*) filedat
    read(1,*) fileeb
    read(1,*) fileprefix
     read(1,*) nvar,nreps,ymiss,pfixed,ptheta,pomega,pto,maxj,cutoff
     discard = 0
     if(.not. fp_equal(cutoff, 0.0d0)) discard = 1
    miss = 1
    maxits = 100
    if(fp_equal(ymiss, 0.0d0)) miss = 0
    nvar2 = 1+max(pfixed,0)+max(pomega,0)+max(ptheta,0)+max(pto,0)
    allocate(var2ind(nvar2))
    allocate(var2label(nvar2))
    allocate(jcodes(maxj))
    read(1,*) id2ind,var2ind(1)
    k = 1
    IF (Pfixed .GE. 1) THEN
        READ(1,*) (var2ind(k+I), I=1,Pfixed)
        k = k + pfixed
     END IF
    IF (Ptheta .GE. 1) THEN
        READ(1,*) (var2IND(k+I), I=1,Ptheta)
        k = k + ptheta
    END IF
    IF (Pomega .GE. 1) THEN
        READ(1,*) (var2IND(k+I), I=1,Pomega)
        k = k + pomega
     END IF
    IF (Pto .GE. 1) THEN
        READ(1,*) (var2IND(k+I), I=1,Pto)
    END IF
    read(1,*) (jcodes(k),k=1,maxj)
! read in the labels
     READ(1,*) var2label(1)
    k = 1
    IF (Pfixed .GE. 1) THEN
        READ(1,*) (var2label(k+I), I=1,Pfixed)
        k = k + pfixed
     END IF
    IF (Ptheta .GE. 1) THEN
        READ(1,*) (var2label(k+I), I=1,Ptheta)
        k = k + ptheta
    END IF
    IF (Pomega .GE. 1) THEN
        READ(1,*) (var2label(k+I), I=1,Pomega)
        k = k + pomega
     END IF
    IF (Pto .GE. 1) THEN
        READ(1,*) (var2label(k+I), I=1,Pto)
     END IF
    CLOSE(1)
! write out the DEF file 
    OPEN(1,FILE=trim(fileprefix)//".def")
    WRITE(1,5)FILEDAT
    write(1,*) fileeb
    WRITE(1,5)FILEprefix
     write(1,*) pfixed,ptheta,pomega,pto,maxj,cutoff
    write(1,*) var2ind(1)
    k=1
    IF (Pfixed .GE. 1) THEN
        write(1,*) (var2ind(k+I), I=1,Pfixed)
        k = k + pfixed
     END IF
    IF (Ptheta .GE. 1) THEN
        write(1,*) (var2IND(k+I), I=1,Ptheta)
        k = k + ptheta
    END IF
    IF (Pomega .GE. 1) THEN
        write(1,*) (var2IND(k+I), I=1,Pomega)
        k = k+pomega
     END IF
    IF (Pto .GE. 1) THEN
        write(1,*) (var2IND(k+I), I=1,Pto)
     END IF
    write(1,*) (jcodes(k),k=1,maxj)

     write(1,*) var2label(1)
     k=1
     IF (Pfixed .GE. 1) THEN
        write(1,*) (var2Label(k+I), I=1,Pfixed)
        k = k + pfixed
     END IF 
     IF (Ptheta .GE. 1) THEN
        write(1,*) (var2Label(k+I), I=1,Ptheta)
        k = k + ptheta
     END IF
     IF (Pomega .GE. 1) THEN
        write(1,*) (var2Label(k+I), I=1,Pomega)
        k = k + pomega
     END IF
     IF (Pto .GE. 1) THEN
        write(1,*) (var2Label(k+I), I=1,Pto)
     END IF
    CLOSE(1)

    open(5*nreps+6,file=fileeb)
    read(5*nreps+6,*) nsubj,r,ns,ndatasets
    numRE = r+ns
    nvar3 = 4+R+pfixed+pomega+R*ptheta+pto
    write(*,*) nsubj,nvar2,nvar3
    allocate(thetatemp(nsubj*numRE))
    allocate(thetas(nsubj,numRE))
    allocate(tempdata(nvar3))
    allocate(betas(nreps,nvar3))
    allocate(meanbetas(nvar3))
    allocate(varbetas(nvar3))
    allocate(totalvarhats(nvar3))
    allocate(temp3(nvar3))
    allocate(liks(nreps))
    allocate(isgood(nreps))
    allocate(dev(nreps))
    totalvarhats = 0
    isgood = 1

    CALL READAT(FILEDAT,NC2,NOBS,MAXK,NVAR,nvar2,tempsums,IDNI,&
                ID2IND,var2ind,MISS,YMISS,ids)

        allocate(intLabel(nvar3))
        intlabel(1) = var2label(1)
        intlabel(2) = "Intercept             "
        intlabel(3:pfixed+3) = var2label(2:1+pfixed)
        do j=1,R
            write(mystr, '(A6, I1, A17)') "Locat_",j,"                 "
            intLabel(1+Pfixed+1+(Ptheta+1)*(j-1)+1) = mystr
            do i=1,pTheta
                write(mystr, '(A6, I1, A17)') "Locat_",j,"*"//var2label(1+pfixed+i)
                intlabel(1+pfixed+1+(Ptheta+1)*(j-1)+1+i) = mystr
            end do
        end do
        if(pOmega >= 0) intlabel(1+pfixed+1+(1+Ptheta)*R+1) = "Scale                 "
        if(pOmega > 0) intlabel(1+pfixed+1+(1+Ptheta)*R+2:1+pfixed+1+(1+Ptheta)*R+1+pomega) = &
            "Scale*"//var2label(1+pfixed+ptheta+1:1+pfixed+ptheta+pomega)
        if(pTO >= 0) intlabel(1+pfixed+1+(1+Ptheta)*R+1+pomega+1) = "Locat_1*Scale           "
        if(pTO >= 1) intlabel(1+pfixed+1+(1+Ptheta)*R+1+pomega+2:1+pfixed+1+(1+Ptheta)*R+1+pomega+pto) = &
            "L*S*"//var2label(1+pfixed+ptheta+pomega+1:1+pfixed+ptheta+pomega+pto)

        filetempdat = trim(fileprefix) // "_.dat"
    OPEN(1,FILE="mixor.def")
        WRITE(1,9) "                                                                               "
        WRITE(1,9) "                                                                               "
    9 FORMAT(18A4)
        WRITE(1,5) adjustl(filetempdat)
        if(len(trim(fileprefix)) > 34) then
            temp80 = fileprefix(1:34) // "_x.out"
        else
            temp80 = trim(fileprefix) // "_x.out"
        end if
    5 FORMAT(A80)
    WRITE(1,5) adjustl(temp80)
    WRITE(1,'("temp_.def")')

    if(discard .eq. 1) then
        WRITE(1,'(4I3,E10.1E3,17I3,i4)') 0, nvar3+1, 0, nvar3-1, .0005, maxj, 2, 0, 0, 0, 11, 1, 1, 0,0, 0, 0, 0, 0, 0, 0, 0, maxits
    else
        WRITE(1,'(4I3,E10.1E3,17I3,i4)') 0, nvar3+1, 0, nvar3-1, .0005, maxj, 0, 0, 0, 0, 11, 1, 1, 0,0, 0, 0, 0, 0, 0, 0, 0, maxits
    end if
    WRITE(1,'(20I3)') 1, 2
    write(1,'(20I3)') (i+2, i=1, nvar3-1)
    write(1,'(20I3)')(jCODEs(J), J = 1,MAXJ)
    if(discard .eq. 1) then
        write(1,*) -999
        write(1,*) (-999, j=1, pfixed+2), cutoff, (-999, j=1,nvar3-4-pfixed)
    end if
    WRITE(1,*) intlabel(1)
    write(1,*) (intlabel(i), i=2, nvar3)
    CLOSE(1)

    do n=1, nreps
        read(5*nreps+6,*) j,(thetatemp(k),k=1,nsubj*numRE)
        do i=1,nsubj
            thetas(i,:) = thetatemp((i-1)*numRe+1:i*numRe)
        end do
        write(*,*) n!,(thetatemp(k),k=1,nsubj*r)
    ! READ THE DATA 

    ! THE DATA MUST BE SORTED BY THE ID VARIABLE

    ! NOTE THAT Y, X, U, W, IDNI ARE ALLOCATED IN READAT
        open(n*5,file=filetempdat)
        DO I=1,NC2  ! go over level-2 clusters
            tempdata(1) = tempsums(i,1)
            k=1
            diffk = 0
            if(pfixed .ge. 0) then
                tempdata(2) = 1d0
                diffk = diffk + 1
                if(pfixed .ge. 1) then
                    tempdata(k+diffk+1:k+diffk+pfixed) = tempsums(i,k+1:k+pfixed)
                    k = k + pfixed
                end if
            end if
            if(ptheta .ge. 0) then
                do j=1,R
                    tempdata(k+diffk+1) = thetas(i,j)
                    diffk = diffk + 1
                    if(ptheta .ge. 1) then
                        tempdata(k+diffk+1:k+diffk+ptheta) = thetas(i,j)*tempsums(i,k+1:k+ptheta)
                        if(j .ne. R) diffk = diffk + ptheta
                    end if
                end do
                k = k + ptheta
            end if
            if(pomega .ge. 0) then
                tempdata(k+diffk+1) = thetas(i,R+1)
                diffk = diffk + 1
                if(pomega .ge. 1) then
                    tempdata(k+diffk+1:k+diffk+pomega) = thetas(i,R+1)*tempsums(i,k+1:k+pomega)
                    k = k + pomega
                end if
            end if
            if(pto .ge. 0) then
                tempdata(k+diffk+1) = thetas(i,1)*thetas(i,R+1)
                diffk = diffk + 1
                if(pto .ge. 1) then
                    tempdata(k+diffk+1:k+diffk+pto) = thetas(i,1)*thetas(i,R+1)*tempsums(i,k+1:k+pto)
                    k = k + pto
                end if
            end if
            write(n*5,'(i16,25e15.6)') idni(i,1),(tempdata(k),k=1,nvar3)
        end do
        close(n*5)

        CALL SYSTEM("./mixor > temp_")
        open(n*5+4, file="mixor.lik")
        read(n*5+4,*) liks(n), numits
        close(n*5+4)
        if(numits >= maxits) then
            isgood(n) = 0
        else
            open(n*5+2, file="mixor.est")
            do j=1,nvar3-1
                read(n*5+2,*) mystr, betas(n,j)
            end do
            close(n*5+2)
            open(n*5+3, file="mixor.var")
            do j=1,nvar3-1
                read(n*5+3,*) (temp3(i),i=1,nvar3-1)
                totalvarhats(j) = totalvarhats(j) + temp3(j)
            end do
            close(n*5+3)
        end if
    end do
    numgood = sum(isgood)
    logl = dot_product(liks, isgood)/nreps
    loglsd = sqrt(dot_product((liks - logl)**2, isgood)/(numgood-1))
        write(mystr, '(I5)') nreps
        fileout = trim(fileprefix) // "_"//trim(adjustl(mystr)) //".out"
    open(2,file=fileout)
    write(2,*) "Number of successful replications =", numgood
    write(2,*)
    WRITE(2,'("-------------")')
    WRITE(2,'("Final Results")')
    WRITE(2,'("-------------")')
    
   WRITE(2,562)LOGL,loglsd,LOGL-nvar3+1, &
   LOGL-0.5D0*DBLE(nvar3-1)*DLOG(DBLE(NC2)),0-2*LOGL,0-2*(LOGL-nvar3+1), &
   0-2*(LOGL-0.5D0*DBLE(nvar3-1)*DLOG(DBLE(NC2)))
   562 FORMAT(1X,'Average Log Likelihood         = ',F12.3,' (sd=',F6.3,')'/, &
         1X,"Akaike's Information Criterion = ",F12.3,/, &
         1X,"Schwarz's Bayesian Criterion   = ",F12.3,//,&
         1X,"==> multiplied by -2             ",      /  &
         1X,'Log Likelihood                 = ',F12.3,/, &
         1X,"Akaike's Information Criterion = ",F12.3,/, &
         1X,"Schwarz's Bayesian Criterion   = ",F12.3,/)
            WRITE(2,57)
            57 FORMAT(/,'Variable',12x,'    Estimate',4X,'AsymStdError',4x, &
                      '     z-value',4X,'     p-value',/,'----------------',4x,  &
                      '------------',4X,'------------',4X,'------------',4X,'------------')
    varbetas = 0
    totalvarhats(:) = totalvarhats(:) / numgood
    do j=1,nvar3-1
        meanbetas(j) = dot_product(isgood,betas(:,j))/numgood
        do i=1,nreps
            varbetas(j) = varbetas(j) + isgood(i)*(betas(i,j)-meanbetas(j))**2
        end do
        varbetas(j) = varbetas(j) / (numgood - 1)
        ZVAL = meanbetas(j)/sqrt(totalvarhats(j)+varbetas(j))
        PVAL = 2.0D0 *(1.0D0 - PHIFN(DABS(ZVAL)))
        write(*,'(A16,4(4x,F10.4))') intlabel(j+1), meanbetas(j), sqrt(totalvarhats(j)+varbetas(j)),zval,pval
        write(2,804) intlabel(j+1), meanbetas(j), sqrt(totalvarhats(j)+varbetas(j)),zval,pval
    end do
    
     804 FORMAT(A16,4(4x,F12.5))
     call system("mkdir work")
     call system("mv mixor.est work/")
     call system("mv mixor.var work/")
     call system("mv mixor.lik work/")
     call system("mv mixor.its work/")
     call system("mv mixor.def work/")
     call system("mv temp_* work/")
     call system("mv *_.dat work/")
     call system("mv *_x.out work/")

contains

SUBROUTINE READAT(FILEDAT,NC2,NOBS,MAXK,NVAR,nvar2,tempsums,IDNI,ID2IND,var2ind,miss,YMISS,ids)
    implicit none
    CHARACTER(LEN=80),intent(in) :: FILEDAT
    integer,intent(in) :: nvar,miss,id2ind,nvar2
    integer,intent(out) :: nc2,nobs,maxk
    REAL(KIND=8),intent(in) :: YMISS
    REAL(KIND=8),ALLOCATABLE,intent(out):: tempsums(:,:)
    INTEGER,ALLOCATABLE,intent(in) :: var2ind(:)
    INTEGER,ALLOCATABLE,intent(out) :: idni(:,:),ids(:)

    INTEGER :: myPASS,I,K,ICOUNT,myindex,IDTEMP,IDOLD,hasmiss
    REAL(KIND=8),ALLOCATABLE:: TEMPR(:)
    LOGICAL FIRST, FP_EQUAL

    ALLOCATE (TEMPR(NVAR))

   ! INITIALIZE
    DO myPASS = 1,2
        IF (myPASS .EQ. 2) THEN
                allocate(ids(icount))
                allocate(tempsums(nc2,nvar2))

              tempsums = 0
        ! IDNI has IDs and Nobs per ID
              ALLOCATE (IDNI(NC2,2))
              IDNI = 0
        ENDIF
   
        I     = 1
        K     = 1
        MAXK  = 0
        ICOUNT= 0
        NOBS  = 0
        FIRST  = .TRUE.

        ! READ IN DATA UNTIL END
        OPEN(1,ACTION='READ',FILE=FILEDAT)

        DO   ! loop forever
        
              READ(1,*,END=1999)(TEMPR(myindex),myindex=1,NVAR)
                hasmiss = 0
                IF (MISS .EQ. 1) THEN
                    do myindex = 1,nvar2
                        IF (FP_EQUAL(tempr(var2IND(myindex)), YMISS)) THEN
                            hasMISS = 1
                            exit
                        END IF
                    end do
                end if
                IF (hasMISS .NE. 0) THEN
                    CYCLE  ! give up on current value and go read next one
                end if

      ! QUERY FOR NEW ID AND SET PARAMETERS ACCORDINGLY

              IDTEMP = INT(TEMPR(ID2IND))
              ! QUERY FOR NEW ID AND SET PARAMETERS ACCORDINGLY

              IF (.NOT. FIRST) THEN 
                  ! if r=0 and rnint=1 then NO random effects 
                 IF (R .GE. 1 .AND. IDTEMP .EQ. IDOLD) THEN
                    K     = K+1
                 ELSE
                    IF (myPASS .EQ. 2) THEN
                       IDNI(I,1) = IDOLD
                       IDNI(I,2) = K
                    ENDIF
                    NOBS = NOBS+K
                    IF (K .GT. MAXK) MAXK = K
                    K     = 1
                    I     = I+1
                 ENDIF
              ENDIF

              ! PUT TEMPORARY VALUES INTO DATA VECTORS AND MATRICES

              IDOLD = IDTEMP
              ICOUNT = ICOUNT+1

              FIRST  = .FALSE.
              IF (myPASS == 2) THEN
                    ids(icount) = idtemp
                    do myindex=1,nvar2
                        if(var2ind(myindex) .ne. 0) tempsums(i,myindex) = tempsums(i,myindex) + tempr(var2ind(myindex))
                    end do
              END IF

    END DO   ! loop back to read next line

    ! cleanup final entry
    1999   IF (myPASS .EQ. 2) THEN
              IDNI(I,1) = IDOLD
              IDNI(I,2) = K
           ENDIF
    NOBS = NOBS+K
    IF (K .GT. MAXK) MAXK = K

    NC2 = I
    CLOSE(1)
   END DO   ! two passes, one to get size, second to read data
        do i=1,nc2
            tempsums(i,1:nvar2) = tempsums(i,1:nvar2)/idni(i,2)
        end do
   DEALLOCATE(TEMPR)

   !  END OF READAT
   RETURN 
   END SUBROUTINE READAT

end program repeat_mixor

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

! ************************************************
! SUBROUTINE PHIFN
! Calculate the probability distribution function (Intercept) for
! Normal distributions:
! ************************************************
REAL(kind=8) FUNCTION PHIFN(Z)
    implicit none
    real(kind=8),intent(in)::z
    real(kind=8)::z2,ord,e,g
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
END FUNCTION PHIFN
