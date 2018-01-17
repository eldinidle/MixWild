! 2-level random-intercepts LOCATION SCALE REGRESSION MODEL
!  with variables influencing BS and WS variance
!  Newton-Raphson solution
!  puts error variance into W
!  casts random effect as VARIANCE (not sd)
!  adds in random scale, correlated with random location
!
! includes effects of the scale STD DEV in terms of a log-linear model

! ******************************************
! THE DATA MUST BY SORTED BY THE ID VARIABLE
! ******************************************

PROGRAM MIXREGLS_subject
    implicit none
    INTEGER :: I,NOBS,NVAR,NQ,AQUAD,ID2IND,YIND,P,R,S,PP,RR,SS,MISS,MAXK,NC2,&
                MAXIT,NCENT,PNINT,RNINT,SNINT,POLD,ROLD,SOLD,j,ll,h,ko,kv,nvar2,num0,&
                pfixed,ptheta,pomega,pto,k,nvar3,ncov,ns,pv,rv,sv,nv,nreps,no2nd,nors,discard0
    INTEGER,ALLOCATABLE :: XIND(:),UIND(:),WIND(:),IDNI(:,:),var2IND(:),varind(:)
    REAL(KIND=8) :: RIDGEIN,MEANY,STDY,CONV,YMISS,MINY,MAXY,SDLV,RCORR,RTEMP,&
                    temp,xb,errv,cutoff
    REAL(KIND=8),ALLOCATABLE:: TEMPR(:),Y(:),X(:,:),U(:,:),W(:,:), &
                               MEANX(:),STDX(:),MEANU(:),STDU(:),MEANW(:),STDW(:), &
                               MINX(:),MAXX(:),MINU(:),MAXU(:),MINW(:),MAXW(:), &
                               BETA(:),TAU(:),SPAR(:), tempsums(:,:),tempr2(:),&
                               IDMV(:,:),alpha(:),thetas(:,:),thetavs(:,:),&
                               tempdata(:),tempVector(:),var(:,:),varavg(:,:)
    character(len=12),allocatable::varlabel(:)
    CHARACTER(LEN=16) :: YLABEL
    CHARACTER(LEN=16),ALLOCATABLE :: BLAB(:),ALAB(:),tlab(:),var2label(:)
    CHARACTER(LEN=24),ALLOCATABLE :: intlabel(:)
    CHARACTER(LEN=4) :: HEAD(36)
    CHARACTER(LEN=80) :: FILEDAT, FILEprefix
    character(len=86) :: fileout
    CHARACTER(LEN=160) :: FILEOUT2
    character(len=24) :: mystr
    logical :: fp_equal,longdef

    OPEN(1, FILE='mixregls_random_mixreg.def')
    READ(1,9) HEAD
    9 FORMAT(18A4)
    READ(1,*)FILEDAT
    READ(1,*)FILEprefix
    5 FORMAT(A80)

    pv = 0
    rv = 0
    sv = 0
    longdef = .false.
    READ(1,*) NVAR, P, R, S, PNINT, RNINT, SNINT, CONV

   IF (conv > .9 .or. fp_equal(conv,0d0)) THEN
        BACKSPACE 1
        READ(1,*) NVAR,P,R,S,PNINT,RNINT,SNINT,pv,rv,sv,CONV,NQ,AQUAD,MAXIT,yMISS,NCENT,NCOV,RIDGEIN,&
                nreps,cutoff,nors,no2nd,discard0
        longdef = .true.
    else
        BACKSPACE 1
        READ(1,*) NVAR, P, R, S, PNINT, RNINT, SNINT, CONV,NQ,AQUAD,MAXIT,yMISS,NCENT,NCOV,RIDGEIN,nreps,cutoff,nors,no2nd,discard0
   ENDIF

    ! set SNINT=0 (for the error variance) if SNINT=1 AND S=0
    ! DO NOT ALLOW A MODEL TO BE FIT WITHOUT AN ERROR VARIANCE
    IF (S==0 .AND. SNINT==1) SNINT=0
    miss = 1
    if(fp_equal(ymiss, 0.0d0)) miss = 0
        
    !   SCALEP=0.50D0
    ! nvar     =  number of variables
    ! ridgein  =  initial value of the ridge
    ! scalep   =  scale parm starting value (as proportion of the error variance starting value)
    ! nq       =  number of quad pts (2 to 20)
    ! aquad    =  adaptive quadrature 0=no 1=yes
    ! maxit    =  maximum number of iterations
    ! idIND    =  field for Level-2 ID
    ! yIND     =  field for continuous outcome
    ! p        =  number of covariates
    ! r        =  number of random effect variance terms
    ! s        =  number of error variance terms
    ! pnint    =  1 if no intercept for mean model (0 otherwise)
    ! rnint    = 1 if no intercept for BS variance model 
    ! snint    = 1 if no intercept for WS variance model
    ! ncent    =  1 for standardizing of all RHS variables (0 for no standardizing)
   
    READ(1,*) ID2IND, YIND

! intercepts are needed for all
    IF (P .GE. 1) THEN
        ALLOCATE(XIND(P))
        READ(1,*) (XIND(I), I=1,P)
    END IF
        IF (R .GE. 1) THEN
        ALLOCATE(UIND(R))
        READ(1,*) (UIND(I), I=1,R)
    END IF
    IF (S .GE. 1) THEN
        ALLOCATE(WIND(S))
        READ(1,*) (WIND(I), I=1,S)
    END IF

    nv = pv + rv + sv
     if(nv > 0) then
        allocate(varind(nv))
        j = 0
        allocate(varlabel(nv))
    end if
     IF (Pv .GE. 1) THEN
        READ(1,*) (varind(I), I=1,Pv)
        j = pv
     END IF
    IF (Rv .GE. 1) THEN
        READ(1,*) (varIND(I+j), I=1,Rv)
        j = j + rv
     END IF
     IF (Sv .GE. 1) THEN
        READ(1,*) (varIND(I+j), I=1,Sv)
     END IF

    POLD=P
    ROLD=R
    SOLD=S

    p = p + 1 - pnint + 2*pv
    r = r + 1 - rnint + 2*rv
    s = s + 1 - snint + 2*sv

! read in the labels
     READ(1,*) YLABEL
     IF (P .GE. 1) THEN
        ALLOCATE(BLAB(P))
        if(pold > 0) READ(1,*) (BLAB(I+1-pnint), I=1,Pold)
     END IF
     IF (R .GE. 1) THEN
        ALLOCATE(ALAB(R))
        if(rold > 0) READ(1,*) (ALAB(I+1-rnint), I=1,Rold)
     END IF
     IF (S .GE. 1) THEN
        ALLOCATE(TLAB(S))
        if(sold > 0) READ(1,*) (TLAB(I+1-snint), I=1,Sold)
     END IF
     j = 0
     IF (Pv .GE. 1) THEN
        READ(1,*) (varlabel(I), I=1,Pv)
        j = pv
     END IF
     IF (Rv .GE. 1) THEN
        READ(1,*) (varlabel(I+j), I=1,Rv)
        j = j + rv
     END IF
     IF (Sv .GE. 1) THEN
        READ(1,*) (varlabel(I+j), I=1,Sv)
     END IF
     if(no2nd .ne. 1) then
        read(1,*) pfixed,ptheta,pomega,pto
         if(nors .ne. 0) then
            pomega = -1
            pto = -1
        end if
        nvar2 = 1+max(pfixed,0)+max(pomega,0)+max(ptheta,0)+max(pto,0)
        nvar3 = 5+pfixed+pomega+ptheta+pto
        allocate(var2ind(nvar2))
        allocate(var2label(nvar2))
        read(1,*) var2ind(1)
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
    else
        pfixed = -1
        ptheta = -1
        pomega = -1
        pto = -1
        nvar2 = 0
        nvar3 = 0
    end if
    CLOSE(1)
! write out the DEF file 
    OPEN(1,FILE=trim(fileprefix)//".def")
    WRITE(1,9) HEAD
    WRITE(1,5)FILEDAT
    WRITE(1,5)FILEprefix
!    WRITE(1,5)FILEDEF
     if(longdef) then
        WRITE(1,'(10I3,E10.1E3, i4, i2, i5, f12.3, 2i2, f6.3, i4, f10.3, 3i2)') NVAR, Pold, Rold, Sold, PNINT, RNINT, SNINT, &
                pv, rv, sv, CONV, NQ, AQUAD, MAXIT, yMISS, NCENT, NCOV, RIDGEIN, nreps, cutoff, nors, no2nd, discard0
    else
        WRITE(1,'(7I3,E10.1E3, i4, i2, i5, f12.3, 2i2, f6.3, i4, f10.3, 3i2)') NVAR, Pold, Rold, Sold, PNINT, RNINT, SNINT, &
                CONV, NQ, AQUAD, MAXIT, yMISS, NCENT, NCOV, RIDGEIN, nreps, cutoff, nors, no2nd, discard0
    end if    
    WRITE(1,'(20I3)') ID2IND, YIND
    IF (P .GE. 1) THEN
        WRITE(1,'(20I3)') (XIND(I), I=1,Pold)
    END IF
    IF (R .GE. 1) THEN
        WRITE(1,'(20I3)') (UIND(I), I=1,Rold)
    END IF
    IF (S .GE. 1) THEN
        WRITE(1,'(20I3)') (WIND(I), I=1,Sold)
    END IF
    j=0
     IF (Pv .GE. 1) THEN
        write(1,'(20I3)') (varind(I), I=1,Pv)
        j = pv
     END IF
    IF (Rv .GE. 1) THEN
        write(1,'(20I3)') (varIND(I+j), I=1,Rv)
        j = j + rv
     END IF
     IF (Sv .GE. 1) THEN
        write(1,'(20I3)') (varIND(I+j), I=1,Sv)
     END IF
     
    WRITE(1,*) YLABEL
    IF (P .GE. 1) THEN
        WRITE(1,*) (BLAB(I+1-pnint), I=1,Pold)
    END IF
    IF (R .GE. 1) THEN
        WRITE(1,*) (ALAB(I+1-rnint), I=1,Rold)
    END IF
    IF (S .GE. 1) THEN
        WRITE(1,*) (TLAB(I+1-snint), I=1,Sold)
    END IF
     j = 0
     IF (Pv .GE. 1) THEN
        write(1,*) (varlabel(I), I=1,Pv)
        j = pv
     END IF
    IF (Rv .GE. 1) THEN
        write(1,*) (varlabel(I+j), I=1,Rv)
        j = j + rv
     END IF
     IF (Sv .GE. 1) THEN
        write(1,*) (varlabel(I+j), I=1,Sv)
     END IF
     if(no2nd .ne. 1) then
         write(1,*) pfixed,ptheta,pomega,pto
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

    ! write the labels
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
    end if
    CLOSE(1)

    fileout = trim(fileprefix) // "_1.out"
! READ THE DATA 

! THE DATA MUST BE SORTED BY THE ID VARIABLE

! NOTE THAT Y, X, U, W, IDNI ARE ALLOCATED IN READAT
    CALL READAT(FILEDAT,NC2,NOBS,MAXK,NVAR,R,P,S,nv,nvar2,&
                Y,X,U,W,var,varavg,tempsums,IDNI,&
                ID2IND,YIND,XIND,UIND,WIND,varind,var2ind,&
                MISS,YMISS,pold,rold,sold,pnint,rnint,snint,&
                discard0,num0,fileprefix)

    if(pnint .ne. 1) then
        x(:,1) = 1
        blab(1) = "Intercept   "
    end if
    if(rnint .ne. 1) then
        u(:,1) = 1
        alab(1) = "Intercept   "
    end if
    if(snint .ne. 1) then
        w(:,1) = 1
        tlab(1) = "Intercept   "
    end if
    do j=1, pv
        ll = pold + 1 - pnint
        blab(ll+j*2-1) = trim(varlabel(j)) // "_BS"
        blab(ll+j*2) = trim(varlabel(j)) // "_WS"
    end do
    do j=1, rv
        ll = rold + 1 - rnint
        alab(ll+j*2-1) = trim(varlabel(j+pv)) // "_BS"
        alab(ll+j*2) = trim(varlabel(j+pv)) // "_WS"
    end do
    do j=1, sv
        ll = sold + 1 - snint
        tlab(ll+j*2-1) = trim(varlabel(j+pv+rv)) // "_BS"
        tlab(ll+j*2) = trim(varlabel(j+pv+rv)) // "_WS"
    end do
    ko = 0
    do i=1,nc2
        do h = 1,idni(i,2)
            ko = ko + 1
            do j=1, pv
                ll = pold + 1 - pnint
                x(ko,ll+j*2-1) = varavg(i,j)
                x(ko,ll+j*2) = var(ko,j) - varavg(i,j)
            end do
            kv = pv
            do j=1, rv
                ll = rold + 1 - rnint
                u(ko,ll+j*2-1) = varavg(i,j+kv)
                u(ko,ll+j*2) = var(ko,j+kv) - varavg(i,j+kv)
            end do
            kv = kv + rv
            do j=1, sv
                ll = sold + 1 - snint
                w(ko,ll+j*2-1) = varavg(i,j+kv)
                w(ko,ll+j*2) = var(ko,j+kv) - varavg(i,j+kv)
            end do
        end do
    end do

    PP = ((P+1)*P)/2
    SS = ((S+1)*S)/2
    RR = ((R+1)*R)/2

    ! get the mean of the subject means and Variances
    CALL SUBMANDV(Y,IDNI,NC2,IDMV,RCORR,SDLV)


! get descriptive stats on all of the variables
    meany=SUM(y(1:nobs))/DBLE(nobs)
    miny=minval(y(1:nobs))
    maxy=maxval(y(1:nobs))
    ALLOCATE(TEMPR(nobs))
    tempR(:)=0.0D0
    tempR(:)=(y(1:nobs)-meany)**2
    RTEMP=SUM(tempR)/DBLE(nobs-1)
    stdy=DSQRT(RTEMP)
    DEALLOCATE(TEMPR)
    if (p>0) then
    ALLOCATE(meanx(p))
        ALLOCATE(stdx(p))
        ALLOCATE(minx(p))
        ALLOCATE(maxx(p))
        call descripc(x,nobs,p,meanx,stdx,minx,maxx)
        ! standardize x if ncent=1
        if (ncent == 1) CALL STANDZ(X,NOBS,P,MEANX,STDX)
    end if
    if (r>0) then
        ALLOCATE(meanu(r))
        ALLOCATE(stdu(r))
        ALLOCATE(minu(r))
        ALLOCATE(maxu(r))
        call descripc(u,nobs,r,meanu,stdu,minu,maxu)
        ! standardize u if ncent=1
        if (ncent == 1) CALL STANDZ(U,NOBS,R,MEANU,STDU)
    end if
    if (s>0) then
        ALLOCATE(meanw(s))
        ALLOCATE(stdw(s))
        ALLOCATE(minw(s))
        ALLOCATE(maxw(s))
        call descripc(w,nobs,s,meanw,stdw,minw,maxw)
        ! standardize w if ncent=1
        if (ncent == 1) CALL STANDZ(W,NOBS,S,MEANW,STDW)
    end if

 ! get starting values
     if (p==0) then
        xb=0
     else if (p>0) then
         ALLOCATE(BETA(P))
         CALL STARTBETA(Y,X,NOBS,P,BETA)
     end if
     if (s>0) then
         ALLOCATE(TAU(S))
         ALLOCATE(TEMPR2(NOBS))
         NS=NCOV+1
         ALLOCATE(SPAR(ncov+1))
         CALL STARTTAU(Y,X,W,NOBS,P,S,SDLV,tAU,TEMPR2,ERRV,SPAR)
     END IF 
     IF (R>0) then
         ! note that startalpha uses tempr2 and errv which are obtained from starttau above 
        ALLOCATE(ALPHA(R))
        call STARTALPHA(U,TEMPR2,ERRV,NOBS,R,ALPHA)
     END IF
     DEALLOCATE(TEMPR2)
    allocate(thetas(nc2,2))
    allocate(thetavs(nc2,3))
       CALL RANDINTEM(Y,X,IDNI,NC2,NOBS,P,BETA,ALPHA(1),TAU(1),thetas,thetavs)

       ! PRINT OUT THE descriptives and starting values
       ! writes results out to mixREGLS1.OUT
       CALL PRINTDESC(HEAD,FILEprefix,FILEOUT,CONV,NQ,AQUAD,MAXIT,NOBS,NC2,IDNI,YLABEL,meany,miny,maxy,stdy, &
       NCENT,P,R,S,BLAB,meanx,minx,maxx,stdx,ALAB,meanu,minu,maxu,stdu,TLAB,meanw,minw,maxw,stdw,num0)

    CALL mixREGLSEST(IDNI,Y,X,U,W,BLAB,ALAB,TLAB,NC2,P,R,S,CONV,NQ,AQUAD,MAXIT,NCENT,ncov,RIDGEIN, &
                       BETA,TAU,SPAR,alpha,thetas,thetavs, maxk,nors)

#if defined(_WIN32)
     FILEOUT2 = "COPY mixREGLS51.OUT+mixREGLS52.OUT " // FILEOUT
     CALL SYSTEM(FILEOUT2)
     CALL SYSTEM("DEL mixREGLS51.OUT mixregls52.out")
     call system("mkdir work")
     call system("move mixregls5* work")
#else
     FILEOUT2 = "cat mixREGLS51.OUT mixREGLS52.OUT >> " // FILEOUT
     CALL SYSTEM(FILEOUT2)
     CALL SYSTEM("rm mixREGLS51.OUT mixREGLS52.out")
     call system("mkdir work")
     call system("mv mixREGLS5* work")
#endif

    if(no2nd .ne. 1) then
        allocate(tempdata(nvar3))
        open(32,file=trim(fileprefix)//'_level2.dat')
        DO I=1,NC2  ! go over level-2 clusters
            write(32,'(i16,25f15.6)') idni(i,1),(tempsums(i,k),k=1,nvar2)
        end do
        close(32)
        allocate(intLabel(nvar3))
        intlabel(1) = var2label(1)
        intlabel(2) = "Intercept             "
        intlabel(3:pfixed+3) = var2label(2:1+pfixed)
        j=1
            write(mystr, '(A6, I1, A17)') "Locat_",j,"                 "
            intLabel(1+Pfixed+1+(Ptheta+1)*(j-1)+1) = mystr
            do i=1,pTheta
                write(mystr, '(A6, I1, A17)') "Locat_",j,"*"//var2label(1+pfixed+i)
                intlabel(1+pfixed+1+(Ptheta+1)*(j-1)+1+i) = mystr
            end do
        if(pomega .ne. -1) intlabel(1+pfixed+1+(1+Ptheta)+1) = "Scale                 "
        if(pOmega > 0) intlabel(1+pfixed+1+(1+Ptheta)+2:1+pfixed+1+(1+Ptheta)+1+pomega) = &
            "Scale*"//var2label(1+pfixed+ptheta+1:1+pfixed+ptheta+pomega)
        if(pTO >= 0) intlabel(1+pfixed+1+(1+Ptheta)+1+pomega+1) = "Locat_1*Scale           "
        if(pTO >= 1) intlabel(1+pfixed+1+(1+Ptheta)+1+pomega+2:1+pfixed+1+(1+Ptheta)+1+pomega+pto) = &
            "L*S*"//var2label(1+pfixed+ptheta+pomega+1:1+pfixed+ptheta+pomega+pto)

        OPEN(1, FILE='mix_random.def')
        write(1,*) trim(fileprefix)//'_ebvar.dat'
        write(1,*) trim(fileprefix)//'_ebrandom.dat'
        if(nors .ne. 0) then
            write(1,*) nc2, 1, 0, nreps, 123
        else
            write(1,*) nc2, 1, 1, nreps, 123
        end if
        close(1)

#if defined(_WIN32)
        call system("mix_random")
#else
        call system("./mix_random")
#endif

        open(1, file="repeat_mixreg.def")
        write(1,*) trim(fileprefix)//'_level2.dat'
        write(1,*) trim(fileprefix)//'_ebrandom.dat'
        write(1,*) trim(fileprefix)//"_random"
        write(1,*) nvar2+1,nreps,pfixed,ptheta,pomega,pto,cutoff
        write(1,*) 1, 2
            k=2
            if(pfixed .ge. 1) then
                write(1,*) (k+i, i=1,pfixed)
                k = k + pfixed
            end if
            if(ptheta .ge. 1) then
                write(1,*) (k+i, i=1,ptheta)
                k = k + ptheta
            end if
            if(pomega .ge. 1) then
                write(1,*) (k+i, i=1,pomega)
                k = k + pomega
            end if
            if(pto .ge. 1) then
                write(1,*) (k+i, i=1,pto)
                k = k + pto
            end if
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
        call system("cp repeat_mixreg.def "//trim(fileprefix)//"_repeat_mixreg.def")
        call system("./repeat_mixreg")
        open(3, file=trim(fileprefix)//'_desc2.out')
        
    write(3,9) head
    write(3,*)
    write(3,*) "Level 2 obervations =",nc2
        !CALL SYSTEM("DEL mixREGLS52.OUT")
         ALLOCATE(tempVector(nc2))
    200  FORMAT(1x,A16,4F12.4)
        WRITE(3,'("------------")')
         WRITE(3,'("Descriptives")')
         WRITE(3,'("------------")')
         WRITE(3,*)

        meany=SUM(tempsums(1:nc2,1))/DBLE(nc2)
         miny=minval(tempsums(1:nc2,1))
         maxy=maxval(tempsums(1:nc2,1))
         tempVector(:)=(tempsums(1:nc2,1)-meany)**2
         temp=SUM(tempVector)/DBLE(nc2-1)
         stdy=DSQRT(TEMP)
     WRITE(3,'(" Dependent variable")')
     WRITE(3,'("                         mean         min         max     std dev")') 
     WRITE(3,'(" ----------------------------------------------------------------")')
     WRITE(3,200) var2Label(1),meany,miny,maxy,stdy
     WRITE(3,*)

    write(3,*)
     WRITE(3,'(" Independent variables")')
     WRITE(3,'("                         mean         min         max     std dev")') 
     WRITE(3,'(" ----------------------------------------------------------------")')

    do i=1,pfixed
            meany=sum(tempsums(1:nc2,1+i))/dble(nc2)
            miny=minval(tempsums(1:nc2,1+i))
            maxy=maxval(tempsums(1:nc2,1+i))
            tempVector(:)=(tempsums(1:nc2,1+i)-meany)**2
            TEMP=SUM(tempVector)/DBLE(nc2-1)
            stdy=DSQRT(TEMP)
            WRITE(3,200) var2label(i+1),meany,miny,maxy,stdy
    end do

         WRITE(3,*)
         WRITE(3,*)
         WRITE(3,'(" Random Location and Scale EB mean estimates")')
         WRITE(3,'("                         mean         min         max     std dev")') 
         WRITE(3,'(" ----------------------------------------------------------------")')

        do j=1,2
            meany=sum(thetas(:,j))/dble(nc2)
            miny=minval(thetas(:,j))
            maxy=maxval(thetas(:,j))
            tempVector(:)=(thetas(:,j)-meany)**2
            TEMP=SUM(tempVector)/DBLE(nc2-1)
            stdy=DSQRT(TEMP)
            write(mystr, '(A6, I1, A17)') "Locat_",j,"                 "
            if(j .eq. 2) mystr = "Scale"
            if(j .eq. 1 .or. (j .eq. 2 .and. nors .ne. 1)) WRITE(3,200) mystr,meany,miny,maxy,stdy
        end do

        if(pto >= 0) then
            meany=sum(thetas(:,1)*thetas(:,2))/dble(nc2)
            miny=minval(thetas(:,1)*thetas(:,2))
            maxy=maxval(thetas(:,1)*thetas(:,2))
            tempVector(:)=(thetas(:,1)*thetas(:,2)-meany)**2
            TEMP=SUM(tempVector)/DBLE(nc2-1)
            stdx=DSQRT(TEMP)
            WRITE(3,200) "Locat_1*Scale     ",meany,miny,maxy,stdy
        end if

         close(3)
            write(mystr, '(I5)') nreps

        CALL SYSTEM("cat "//trim(fileprefix)//"_desc2.out "//trim(fileprefix) &
                    //"_random_"//trim(adjustl(mystr))//".out >> "//trim(fileprefix)//"_2.out")
        call system("rm "//trim(fileprefix)//"_desc2.out")
        call system("mv mix_random.def work")
        call system("mv "//trim(fileprefix)//"_ebvar.dat work")
        call system("mv "//trim(fileprefix)//"_random* work")

    end if
    !deallocate(tempsums,tempdata,tempvector)
CONTAINS

! READ THE DATA 

! THE DATA MUST BE SORTED BY THE ID VARIABLE


! NOTE THAT Y, X, U, W, IDNI ARE ALLOCATED IN READAT
SUBROUTINE READAT(FILEDAT,NC2,NOBS,MAXK,NVAR,R,P,S,nv,nvar2,Y,X,U,W,var,varavg,tempsums, &
    IDNI,ID2IND,YIND,XIND,UIND,WIND,varind,nsind,miss,YMISS,pold,rold,sold,pnint,rnint,snint,&
    discard0,num0,fileprefix)
    implicit none
        
    CHARACTER(LEN=80),intent(in) :: FILEDAT,fileprefix
    integer,intent(in) :: nvar,r,p,s,yind,miss,pold,rold,sold,nv,pnint,rnint,snint,id2ind,nvar2,discard0
    integer,intent(out) :: nc2,nobs,maxk,num0
    REAL(KIND=8),intent(in) :: YMISS
    REAL(KIND=8),ALLOCATABLE,intent(out):: Y(:),X(:,:),U(:,:),W(:,:),var(:,:),varavg(:,:),tempsums(:,:)
    INTEGER,ALLOCATABLE,intent(in) :: XIND(:),UIND(:),WIND(:),varind(:),nsind(:)
    INTEGER,ALLOCATABLE,intent(out) :: idni(:,:)

    INTEGER :: myPASS,I,K,ICOUNT,myindex,IDTEMP,IDOLD,hasmiss,nvartotal,discardi,m
    REAL(KIND=8) :: firsty
    REAL(KIND=8),ALLOCATABLE:: TEMPR(:)
    INTEGER,ALLOCATABLE :: allvarsind(:)
    LOGICAL FIRST, FP_EQUAL
    

        ALLOCATE (TEMPR(NVAR))
        nvarTotal = 1+pold+rold+sold+nv+nvar2
        allocate(allvarsIND(nvarTotal))
        allvarsIND(1) = yind
        allvarsIND(2:pold+1) = xind(1:pold)
        allvarsIND(pold+2:pold+rold+1) = uind(1:rold)
        allvarsIND(pold+rold+2:pold+rold+sold+1) = wind(1:sold)
        allvarsIND(pold+rold+sold+2:pold+rold+sold+1+nv) = varind(1:nv)
        allvarsIND(pold+rold+sold+2+nv:nvarTotal) = nsind(1:nvar2)
        num0 = 0
        if(discard0 .ne. 0) open(16, file=trim(fileprefix)//"_removed.dat")

        ! INITIALIZE
        DO myPASS = 1,2
            IF (myPASS .EQ. 2) THEN
               ALLOCATE (Y(ICOUNT))
                ALLOCATE (X(ICOUNT,P))
                ALLOCATE (U(ICOUNT,R))
                ALLOCATE (W(ICOUNT,S))
                allocate (var(icount,nv))
                allocate (varavg(nc2,nv))
                allocate(tempsums(nc2,nvar2))
                
                Y = 0.0D0
                X = 0.0D0
                U = 0.0D0
                W = 0.0D0
              var = 0.0D0
              varavg = 0
        ! IDNI has IDs and Nobs per ID
                ALLOCATE (IDNI(NC2,2))
                IDNI = 0
                tempsums = 0
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
                hasMISS = 0
                IF (MISS .EQ. 1) THEN
                    do myindex = 1,nvartotal
                        IF (FP_EQUAL(tempr(allvarsIND(myindex)), YMISS)) THEN
                            hasMISS = 1
                            exit
                        END IF
                    end do
                end if
                IF (hasMISS .NE. 0) THEN
                    CYCLE  ! give up on current value and go read next one
                end if

                IDTEMP = INT(TEMPR(ID2IND))

                  ! QUERY FOR NEW ID AND SET PARAMETERS ACCORDINGLY

                IF (.NOT. FIRST) THEN 
                      ! if r=0 and rnint=1 then NO random effects 
                    IF ((R .GE. 1 .OR. RNINT .NE. 1) .AND. IDTEMP .EQ. IDOLD) THEN
                        K     = K+1
                    ELSE
                        discardi = 0
                        if(r > 0 .and. discard0 == 1 .and. mypass .eq. 2) then
                            discardi = 1
                            firsty = y(icount-k+1)
                            do j=2,k
                                if(.not. fp_equal(firsty, y(icount-k+j))) discardi = 0
                            end do
                        end if
                        if(discardi .eq. 0) then
                            IF (myPASS .EQ. 2) THEN
                               IDNI(I,1) = IDOLD
                               IDNI(I,2) = K
                            ENDIF
                            NOBS = NOBS+K
                            IF (K .GT. MAXK) MAXK = K
                            I     = I+1
                        else
                            write(*,*) "REMOVED"
                            num0 = num0 + 1
                            icount = icount - k
                            do j=1,k
                                write(16,'(i9,16f10.3)') idold, y(icount+j), (x(icount+j,m),m=1,pold), &
                                        (u(icount+j,m),m=1,rold), (w(icount+j,m), m=1,sold), &
                                        (var(icount+j,m), m=1,nv), (tempsums(i,m)/k,m=1,nvar2)
                            end do
                        tempsums(i,:) = 0
                        varavg(i,:) = 0
                        end if
                        K     = 1
                    ENDIF
                ENDIF

                  ! PUT TEMPORARY VALUES INTO DATA VECTORS AND MATRICES

                IDOLD = IDTEMP
                ICOUNT = ICOUNT+1

                FIRST  = .FALSE.
                IF (myPASS == 2) THEN
                    Y(icount)  = TEMPR(YIND)
                    do myindex=1,pold
                        x(icount,myindex+1-pnint) = tempr(xind(myindex))
                    end do
                    do myindex=1,rold
                        u(icount,myindex+1-rnint) = tempr(uind(myindex))
                    end do
                    do myindex=1,sold
                        w(icount,myindex+1-snint) = tempr(wind(myindex))
                    end do
                    do myindex=1,nv
                        var(icount,myindex) = tempr(varind(myindex))
                        varavg(i,myindex) = varavg(i,myindex) + tempr(varind(myindex))
                    end do
                    do myindex=1,nvar2
                        if(nsind(myindex) .ne. 0) tempsums(i,myindex) = tempsums(i,myindex) + tempr(nsind(myindex))
                    end do
                END IF

            END DO   ! loop back to read next line

            ! cleanup final entry
        1999  discardi = 0
                        if(r > 0 .and. discard0 == 1 .and. mypass .eq. 2) then
                            discardi = 1
                            firsty = y(icount-k+1)
                            do j=2,k
                                if(.not. fp_equal(firsty, y(icount-k+j))) discardi = 0
                            end do
                        end if
                        if(discardi .eq. 0) then
                            IF (myPASS .EQ. 2) THEN
                               IDNI(I,1) = IDOLD
                               IDNI(I,2) = K
                            ENDIF
                            NOBS = NOBS+K
                            IF (K .GT. MAXK) MAXK = K
                        else
                            write(*,*) "REMOVED"
                            num0 = num0 + 1
                            icount = icount - k
                            do j=1,k
                                write(16,'(i9,16f10.3)') idold, y(icount+j), (x(icount+j,m),m=1,pold), &
                                        (u(icount+j,m),m=1,rold), (w(icount+j,m), m=1,sold), &
                                        (var(icount+j,m), m=1,nv), (tempsums(i,m)/k,m=1,nvar2)
                            end do
                            i = i -1
                        end if

            NC2 = I
            CLOSE(1)
        END DO   ! two passes, one to get size, second to read data
        DEALLOCATE(TEMPR)

        !  END OF READAT
        do i=1,nc2
            tempsums(i,1:nvar2) = tempsums(i,1:nvar2)/idni(i,2)
            varavg(i,1:nv) = varavg(i,1:nv)/idni(i,2)
        end do
        RETURN 
    END SUBROUTINE READAT
! nobs     =  number of total observations
   ! NC2      =  number of level-2 clusters


    SUBROUTINE SUBMANDV(Y,IDNI,NC2,IDMV,RCORR,SDLV)

        REAL(KIND=8) :: Y(:),RCORR,SDLV,MEANLV,RTEMP
        REAL(KIND=8),ALLOCATABLE:: IDMV(:,:),TEMPR(:),meanMV(:),stdMV(:)
        INTEGER :: NC2,IDNI(:,:),IC,I,J

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
        REAL(KIND=8),INTENT(IN),DIMENSION(:,:)::a   ! input matrix
        REAL(KIND=8),INTENT(OUT),dimension(:):: m,s ! m is mean vector,s is std vector
        REAL(KIND=8),DIMENSION(nobs)::t             ! temporary array
        INTEGER:: i

        ! calculate mean
        t(:)=0.
        do i=1,cols
            t(:)=a(1:nobs,i)
            m(i)=SUM(t)/nobs
        end do

        ! calculate std
        t(:)=0
        do i=1,cols
            t(:)=(a(1:nobs,i)-m(i))**2
            s(i)=SQRT(SUM(t)/(nobs-1))
        end do

    end subroutine meanc

    subroutine descripc(a,nobs,cols,m,s,minv,maxv)
        !
        !purpose:
        !to calculate the mean, std, min, and max in each column
        !
        INTEGER,INTENT(IN)::nobs,cols                         ! cols is the # of columns
        REAL(KIND=8),INTENT(IN),DIMENSION(:,:)::a             ! input matrix
        REAL(KIND=8),INTENT(OUT),dimension(:):: m,s,minv,maxv ! m is mean vector,s is std vector
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

    subroutine STARTBETA(Y,X,NOBS,P,BETA)
!
!purpose:
!to calculate starting values for the mean regression coefficients 
!
        INTEGER,INTENT(IN)::nobs,p                             ! cols is the # of columns in the X matrix
        REAL(KIND=8),INTENT(IN),DIMENSION(:)::Y                ! input matrix Y
        REAL(KIND=8),INTENT(IN),DIMENSION(:,:)::X              ! input matrix X
        REAL(KIND=8),INTENT(OUT),dimension(:):: BETA           ! starting values for regression coefficients
        REAL(KIND=8),DIMENSION(((P+1)*P)/2)::TEMPR             ! temporary array
        REAL(KIND=8),DIMENSION(P)::TEMPR2,TEMPR3               ! temporary arrays
        REAL(KIND=8)::DET
        INTEGER:: ier
        CALL GRAMT(X,TEMPR,NOBS,P)
        CALL INVS(TEMPR,P,DET,TEMPR3,IER)
        CALL MPYTR(X,Y,TEMPR2,NOBS,P,0,1)
        CALL MPYM(TEMPR,TEMPR2,BETA,P,P,1,0,1)
    END SUBROUTINE STARTBETA

    SUBROUTINE STARTTAU(Y,X,W,NOBS,P,S,SDLV,TAU,TEMPR2,ERRV,SPAR)
    !
    !purpose:
    !to calculate starting values for the ERROR VARIANCE regression coefficients 
    ! and the random scale parameters
        implicit none
        INTEGER,INTENT(IN)::nobs,p,s                          
        REAL(KIND=8),INTENT(IN),DIMENSION(:)::Y                ! input matrix Y
        REAL(KIND=8),INTENT(IN),DIMENSION(:,:)::X              ! input matrix X
        REAL(KIND=8),INTENT(IN),DIMENSION(:,:)::W              ! input matrix W
        REAL(KIND=8),INTENT(OUT),dimension(:):: TAU            ! starting values for regression coefficients
        REAL(KIND=8),INTENT(OUT),dimension(:):: SPAR           ! starting values for random scale parameters
        REAL(KIND=8),INTENT(OUT),dimension(:):: TEMPR2         ! LN of squared residuals
        REAL(KIND=8)::TP1(NOBS),TP2(((S+1)*S)/2),TP3(S)        ! temporary array
        REAL(KIND=8)::DET,ERRV,SDLV,RTEMP,SMALL
        INTEGER:: ier, i
        RTEMP = 1.0D0
        SMALL  = TINY(RTEMP)
        CALL MPYM(X,BETA,TP1,NOBS,P,0,0,1)
        tempR2=(Y(1:nobs) - TP1)**2
        ! add a little number for zero deviations
        do i = 1,nobs
            if (tempR2(i) .LE. SMALL) TEMPR2(i) = SMALL 
        end do
        ERRV = (SUM(TEMPR2)/DBLE(NOBS))  ! error variance ~ SSE/n 
        TEMPR2 = DLOG(TEMPR2)
        CALL GRAMT(W,TP2,NOBS,S)
        CALL INVS(TP2,S,DET,TP3,IER)
        CALL MPYTR(W,TEMPR2,TP3,NOBS,S,0,1)
        CALL MPYM(TP2,TP3,TAU,S,S,1,0,1)
        TAU(1) = DLOG(.6*ERRV)  ! ln of error variance SSE

! number of parameters associated with the random scale part

        SPAR(1) = SDLV      ! random scale STD DEV
    END SUBROUTINE STARTTAU

    SUBROUTINE STARTALPHA(U,TEMPR2,ERRV,NOBS,R,ALPHA)
        !
    !purpose:
    !to calculate starting values for the BS VARIANCE regression coefficients 
    ! note that TEMPR2 and ERRV are input and result from previous call startalpha
        
    INTEGER,INTENT(IN)::nobs,r                         
    REAL(KIND=8),INTENT(IN),DIMENSION(:,:)::U              ! input matrix W
    REAL(KIND=8),INTENT(IN),dimension(:):: TEMPR2          ! LN of squared residuals
    REAL(KIND=8),INTENT(OUT),dimension(:):: ALPHA          ! starting values for regression coefficients
    REAL(KIND=8)::TEMPR(((R+1)*R)/2),TEMPR3(R)             ! temporary arrays
    REAL(KIND=8)::DET,ERRV
    INTEGER:: ier,i
            CALL GRAMT(U,TEMPR,NOBS,R)
            CALL INVS(TEMPR,R,DET,TEMPR3,IER)
            CALL MPYTR(U,TEMPR2,TEMPR3,NOBS,R,0,1)
            CALL MPYM(TEMPR,TEMPR3,ALPHA,R,R,1,0,1)
            ALPHA(1) = DLOG(.4*ERRV)
            ! make sure that starting values are not too big or small
            do i=1,r
                if (alpha(i) .gt.  1.0d0) alpha(i) =  1.0d0
                if (alpha(i) .lt. -1.0d0) alpha(i) = -1.0d0
            end do
    END SUBROUTINE STARTALPHA


    SUBROUTINE RANDINTEM(Y,X,IDNI,NC2,NOBS,P,BETA,LBSVAR,LWSVAR,THETAEM,PVAREM)
!
!purpose:
! EM Algorithm for a random intercept model - 20 ITERATIONS
! this provides excellent starting values for beta, BS and WS variances, and EB estimates (mean and post var)

    
        INTEGER,INTENT(IN)::nobs,p,NC2                          
        INTEGER,INTENT(IN),DIMENSION(:,:)::IDNI                ! ids and nobs per cluster
        REAL(KIND=8),INTENT(IN),DIMENSION(:)::Y                ! input matrix Y
        REAL(KIND=8),INTENT(IN),DIMENSION(:,:)::X              ! input matrix X
        REAL(KIND=8),INTENT(INOUT),dimension(:):: BETA         ! starting values for regression coefficients
        REAL(KIND=8),INTENT(INOUT)::LBSVAR,LWSVAR              ! starting values for Ln(BSvar) and Ln(WSvar)
        REAL(KIND=8),INTENT(OUT),dimension(NC2):: THETAEM,PVAREM ! posterior mean and variance

        REAL(KIND=8)::TEMPR(NC2),TEMPR2(NOBS),THETAEM2(NOBS),TEMPPP((P*(P+1)/2)),TEMPP(P),TEMPP2(P),TEMPNOBS(NOBS),TEMPNOBS2(NOBS) ! temporary arrayS
        REAL(KIND=8)::RTEMP,RHO,RTEMP2,ERRV,DET
        INTEGER:: ier,ITER,I,IC,J

       THETAEM(:) = 0.0D0
       PVAREM(:)  = 0.0D0
        ITERLOOP:DO ITER=1,20
           ! E STEP
           TEMPR(:) = 0.0D0
           RHO = DEXP(LBSVAR) / (DEXP(LBSVAR) + DEXP(LWSVAR))
           DO I=1,NC2
              TEMPR(I) = (DBLE(IDNI(I,2))*RHO)/ (1.0D0 + (DBLE(IDNI(I,2)) - 1.0D0)*RHO)
           END DO
           PVAREM(:) = DEXP(LBSVAR)*(1.0D0 - TEMPR(:))
           TEMPR = TEMPR / DBLE(IDNI(:,2))
           CALL MPYM(X,BETA,TEMPR2,NOBS,P,0,0,1)
           IC=0
           DO I=1,NC2
               RTEMP = 0.0D0
               DO J=1,IDNI(I,2)
                  IC=IC+1
                  RTEMP = RTEMP + (Y(IC) - TEMPR2(IC))
               END DO
               THETAEM(I) = TEMPR(I)*RTEMP
           END DO

       ! create version of THETAEM that is NOBS long
           IC=0
           DO I=1,NC2
               DO J=1,IDNI(I,2)
                  IC=IC+1
                  THETAEM2(IC) = THETAEM(I)
               END DO
           END DO
           ! M STEP
           ! beta
           CALL GRAMT(X,TEMPPP,NOBS,P)
           CALL INVS(TEMPPP,P,DET,TEMPP2,IER)
           TEMPNOBS = Y(1:nobs) - THETAEM2
           CALL MPYTR(X,TEMPNOBS,TEMPP,NOBS,P,0,1)
           CALL MPYM(TEMPPP,TEMPP,BETA,P,P,1,0,1)
        
           ! cluster var
           RTEMP = DOT_PRODUCT(THETAEM,THETAEM)       
           RTEMP2 = SUM(PVAREM)
           RTEMP = (RTEMP + RTEMP2)/DBLE(NC2)
           LBSVAR = DLOG(RTEMP)
           ! error var
           CALL MPYM(X,BETA,TEMPNOBS2,NOBS,P,0,0,1)
           ERRV=0.0D0
           IC=0
           DO I=1,NC2
              DO J=1,IDNI(I,2)
                 IC=IC+1
                 ERRV = ERRV + (TEMPNOBS(IC) - TEMPNOBS2(IC))**2 + PVAREM(I)
              END DO
           END DO
           LWSVAR = DLOG(ERRV/DBLE(NOBS))
      END DO ITERLOOP
!     OPEN(UNIT=3,FILE='MIXREGB5.RES')
!     DO I=1,NC2
!        WRITE(3,'(2F15.6)')THETAEM(I),PVAREM(I)
!     END DO
!     CLOSE(3)
    END SUBROUTINE RANDINTEM


    subroutine CORR(X,Y,N,R)
    !purpose:
    !to calculate the correlation between two vectors
    !
        INTEGER,INTENT(IN)::N                       ! number of rows of two vectors
        REAL(KIND=8),INTENT(IN),DIMENSION(:)::X,Y   ! input vectors
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


! PRINT OUT THE descriptives and starting values
    SUBROUTINE PRINTDESC(HEAD,FILEDAT,FILEprefix,CONV,NQ,AQUAD,MAXIT,NOBS,NC2,IDNI,YLABEL,meany,miny,maxy,stdy, &
               NCENT,P,R,S,BLAB,meanx,minx,maxx,stdx,ALAB,meanu,minu,maxu,stdu,TLAB,meanw,minw,maxw,stdw,num0)

        CHARACTER(LEN=16),INTENT(IN):: YLABEL
        CHARACTER(LEN=16),INTENT(IN),dimension(:):: BLAB,ALAB,TLAB
        CHARACTER(LEN=4),INTENT(IN),DIMENSION(:):: HEAD
        CHARACTER(LEN=80),INTENT(IN):: FILEDAT, FILEprefix
        INTEGER,INTENT(IN)::NQ,AQUAD,MAXIT,NOBS,NC2,NCENT,P,R,S,num0
        INTEGER,INTENT(IN),DIMENSION(:,:)::IDNI                ! ids and nobs per cluster
        REAL(KIND=8),INTENT(IN),dimension(:):: MEANX,MINX,MAXX,STDX,MEANU,MINU,MAXU,STDU,MEANW,MINW,MAXW,STDW
        REAL(KIND=8),INTENT(IN)::CONV,MEANY,MINY,MAXY,STDY
        INTEGER:: IUN,I

             IUN    = 16
             OPEN(UNIT=IUN,FILE="mixREGLS51.OUT")
             
             WRITE(IUN,'("MIXREGLS: Mixed-effects Location Scale Model with BS and WS variance models")')
             write (IUN,*)
             WRITE(IUN,'("-----------------------------")')
             WRITE(IUN,'("mixREGLS.DEF specifications")')
             WRITE(IUN,'("-----------------------------")')
             WRITE(IUN,"(1x,18a4)")HEAD
             WRITE(IUN,*)
             WRITE(IUN,'(" data and output files:")')
             WRITE(IUN,"(1x,a80)")FILEDAT
             WRITE(IUN,"(1x,a80)")trim(FILEprefix)//"_1.out"
             WRITE(IUN,*)
             WRITE(IUN,"(' CONVERGENCE CRITERION = ',F11.8)")CONV
             WRITE(IUN,"(' RIDGEIN    = ',F8.4)")RIDGEIN
            !WRITE(IUN,"('SCALEP     = ',F8.4)")SCALEP
             WRITE(IUN,"(' NQ         = ',I4)")NQ
             WRITE(IUN,"(' QUADRATURE = ',I4,' (0=non-adaptive, 1=adaptive)')")AQUAD
             WRITE(IUN,"(' MAXIT      = ',I4)")MAXIT

             WRITE(IUN,*)
             WRITE(IUN,*)
             WRITE(IUN,'("------------")')
             WRITE(IUN,'("Descriptives")')
             WRITE(IUN,'("------------")')
             WRITE(IUN,*)
             WRITE(IUN,'(" Number of level-1 observations = ",I8)') nobs
             WRITE(IUN,*)
             write(IUN,'(" Number of level-2 clusters     = ",I8)') nc2
             write(IUN,*)
             write(IUN,'(" Number of level-1 observations for each level-2 cluster")')
             write(IUN,'(1x,13I6)') (IDNI(i,2), i=1,nc2)

        200  FORMAT(1x,A8,4F12.4)

             WRITE(IUN,*)
             WRITE(IUN,'(" Dependent variable")')
             WRITE(IUN,'("                 mean         min         max     std dev")') 
             WRITE(IUN,'(" --------------------------------------------------------")')
             WRITE(IUN,200) YLABEL,meany,miny,maxy,stdy
             WRITE(IUN,*)

             if (ncent==1) then
                WRITE(IUN,'(" ==> Standardization of covariates has been selected")') 
                WRITE(IUN,'(" ==> All covariates have mean=0 and std dev=1")') 
                WRITE(IUN,'(" ==> Means and std devs listed below are pre-standardization")') 
                WRITE(IUN,*)
             end if
             if (pold>0) then
                WRITE(IUN,'(" Mean model covariates")')
                WRITE(IUN,'("                 mean         min         max     std dev")') 
                WRITE(IUN,'(" --------------------------------------------------------")')
                do i=1,p
                   WRITE(IUN,200) BLAB(i),meanx(i),minx(i),maxx(i),stdx(i)
                end do
                WRITE(IUN,*)
             end if

             if (rold>0) then
                WRITE(IUN,'(" BS variance model covariates")')
                WRITE(IUN,'("                 mean         min         max     std dev")') 
                WRITE(IUN,'(" --------------------------------------------------------")')
                do i=1,r
                   WRITE(IUN,200) ALAB(i),meanu(i),minu(i),maxu(i),stdu(i)
                end do
                WRITE(IUN,*)
             end if

             if (sold>0) then
                WRITE(IUN,'(" WS variance model covariates")')
                WRITE(IUN,'("                 mean         min         max     std dev")') 
                WRITE(IUN,'(" --------------------------------------------------------")')
                do i=1,s
                   WRITE(IUN,200) TLAB(i),meanw(i),minw(i),maxw(i),stdw(i)
                end do
                WRITE(IUN,*)
             end if
        if(discard0 .ne. 0) WRITE(IUN,508)num0
        508 FORMAT(//,1x,'==> The number of level 2 observations removed because of non-varying responses =', I6)
        if(discard0 .ne. 0) write(IUN,*) '(see '//trim(fileprefix)//'_removed.dat for information about those clusters)'

        CLOSE(IUN)
    END SUBROUTINE PRINTDESC

! estimation subroutine 
    SUBROUTINE mixREGLSEST(IDNI,Y,X,U,W,BLAB,ALAB,TLAB,NC2,P,R,S,CONV,NQ,AQUAD,MAXIT,NCENT,ncov,RIDGEIN, &
        BETA,TAU,SPAR,alpha,thetas,thetavs,maxk,nors)
        
        INTEGER,INTENT(IN)::NC2,P,R,NQ,AQUAD,MAXIT,NCENT,maxk,nors
        INTEGER,INTENT(IN),DIMENSION(:,:)::IDNI                                     ! ids and nobs per cluster
        REAL(KIND=8),INTENT(IN),DIMENSION(:)::Y                      ! input matrix Y, and random location mean and posterior variance
        REAL(KIND=8),INTENT(INOUT),DIMENSION(:)::BETA,TAU,alpha,spar         ! parameters (mean, bs var, ws var, random scale)
        REAL(KIND=8),INTENT(INOUT),DIMENSION(:,:)::thetas,thetavs               ! for adaptive quadrature
        REAL(KIND=8),INTENT(IN),DIMENSION(:,:)::X,U,W                             ! data matrices
        REAL(KIND=8),INTENT(IN)::CONV,RIDGEIN
        ! note that S and ncov vary between the cycles and so cannot be specified as intent(in)
        INTEGER :: S,I,J,L,LL,l2,k,PP,SS, & !RR
                   NPAR,IUN,CYCLES,NCYCLE,RP2,RPS2,NQ1,S0,NS,NS2,IFIN,ITER,NN,NPAR2, &
                   Q,NOB,IER,RIDGEIT,IUNS,NQwR0,NQwR1,ndim,ndim2,myqdim,totalqR0,totalqR1,mytotalq,&
                   ncov0,ncov
        REAL(KIND=8) :: RIDGE,LOGLP,PSUM,DET,LOGL,PI,XB,WT,WSVAR,ERRIJ,LPROB,PRA,LOGDIFF,MAXCORR,ONE,SCAL,&
                        ORIDGE,PVAL,ZVAL,RTEMP,RTEMP2,SMALL,BIG,LOGBIG,MAXDER,PHIFN,sdev,phiRatio,&
                        ua,bsvar,myz,tauhat,tauhatlow,tauhatup
        REAL(KIND=8),ALLOCATABLE:: DER2A(:),DER2B(:),DER(:),LIK(:,:),H(:), &
                                   DERP(:),DERP2(:),DERQ(:),DERQ2(:),DERQQ(:),DZ(:),DERPS(:),COREC(:), &
                                   BDER2(:),ADER2(:),TDER2(:),SDER2(:), &
                                   ABDER2(:),TBDER2(:),TADER2(:),SBDER2(:),SADER2(:),STDER2(:), &
                                   TBADER2(:),SBADER2(:),SBATDER2(:),DZ2(:),SE(:), &
                                   THETA(:),thetav(:),WORK2(:),WORK3(:), &
                                   myweights(:),mypoints(:,:), myweights0(:),&
                                   mypoints0(:,:),pointsR0(:,:),weightsR0(:),pointsR1(:,:),weightsR1(:),&
                                   vaug(:)
        CHARACTER(LEN=16),INTENT(IN),DIMENSION(:):: BLAB(:),ALAB(:),TLAB(:)

            ! parameters
        PI = 3.141592653589793238462643d0
        ONE = 1.0D0
        RTEMP = 1.0D0
        SMALL  = TINY(RTEMP)
        BIG    = HUGE(RTEMP)
        LOGBIG   = DLOG(BIG)

        PP = ((P+1)*P)/2
!        RR = ((R+1)*R)/2
        ndim = 2
        ndim2 = (ndim+1)*ndim/2
        IUN    = 16
        OPEN(UNIT=IUN,FILE="mixREGLS52.OUT")
        IUNS    = 17
        OPEN(UNIT=IUNS,FILE="mixREGLS5_details.ITS")
        open(unit=18, file="MixRegls5.its")
        NQwR1 = NQ**(2)
        nqwR0 = nq
        allocate(weightsR1(NQwr1))
        allocate(pointsR1(NQwr1,2))
        allocate(weightsR0(NQwR0))
        allocate(pointsR0(NQwR0,1))
        call getquad(1, nq, nqwr0, pointsR0, weightsR0)
        call getquad(2, nq, nqwr1, pointsR1, weightsR1)
    ! number of quadrature points, quadrature nodes & weights

        ! ALLOCATE VECTORS FOR NR (these don't vary by CYCLES)
        RP2 = ((R+P+1)*(R+P))/2
        ALLOCATE(BDER2(PP))
        ALLOCATE(ABDER2(P*R))
        ALLOCATE(ADER2(R*(R+1)/2))
        allocate(vaug(ncov+1))
        !Establishing arrays at maximum size needed
        S0 = S
        SS = (S*(S+1))/2
        RPS2 = (RR+P+S)*(RR+P+S+1)/2
        ! NS = number of additional var cov parameters due to random SCALE
        NS = ncov+1
        NS2 = NS*(NS+1)/2
        ncov0 = ncov
        NPAR = P+RR+S+NS ! number of parameters
        NPAR2 = NPAR*(NPAR+1)/2
        totalqR0 = nq ** 1
        totalqR1 = nq ** (2)
        mytotalq = totalqR1
        myqdim = ndim
        ALLOCATE (DER2A(NPAR2))
        ALLOCATE (DER2B(NPAR2))
        ALLOCATE (DER(NPAR))
        ALLOCATE (LIK(nc2,nqwr1))
        ALLOCATE (H(NC2))
        ALLOCATE (derp(npar))
        ALLOCATE (derp2(NPAR2))
        ALLOCATE (DERQ(NPAR))
        ALLOCATE (DERQ2(NPAR2))
        ALLOCATE (DERQQ(NPAR2))
        ALLOCATE (DZ(NPAR))
        dz(:) = 0D0
        ALLOCATE (DERPS(NPAR))
        ALLOCATE (COREC(NPAR))
        corec(:) = 0D0
        ALLOCATE (TBDER2(S*P))
        ALLOCATE (TADER2(S*RR))
        ALLOCATE (TBADER2(S*(P+RR)))
        ALLOCATE (TDER2(SS))
        ALLOCATE (DZ2(NPAR2))
        ALLOCATE (SE(NPAR))
        ALLOCATE (THETA(ndim))
        ALLOCATE (thetav(ndim2))
        ALLOCATE (SBDER2(NS*P))
        ALLOCATE (SADER2(NS*RR))
        ALLOCATE (STDER2(NS*S))
        ALLOCATE (SBADER2(NS*(P+RR)))
        ALLOCATE (SBATDER2(NS*(P+RR+S)))
        ALLOCATE (SDER2(NS2))
        ALLOCATE (WORK2(ndim))
        ALLOCATE (WORK3(ndim2))
        allocate(myweights(nqwR1))
        allocate(mypoints(nqwR1,2))

        ! start cycles
        ! cycles = 1: random intercept model with BS variance terms
        ! cycles = 2: add in scale (WS) variance terms
        ! cycles = 3: add in random scale 
        ! cycles = 4: use NS = R+1

        ncycle = 4
        if(nors .eq. 1) ncycle = 2
        CYCLELOOP:do cycles=1,ncycle
            if(cycles==2 .AND. s0==1) cycle
            if(cycles==4 .and. ncov0==0) cycle
             ! don't do cycle=2 if sold=0 (no WS covariates)
                     
            if (cycles==1) then
                WRITE(IUN,*)
                WRITE(IUN,*)
                WRITE(IUN,'("------------------------------")')
                WRITE(IUN,'("Model without Scale Parameters")')
                WRITE(IUN,'("------------------------------")')
                WRITE(*,*)
                WRITE(*,*) "------------------------------"
                WRITE(*,*) "Model without Scale Parameters"
                WRITE(*,*) "------------------------------"
                WRITE(IUNS,*)
                WRITE(IUNS,*) "------------------------------"
                WRITE(IUNS,*) "Model without Scale Parameters"
                WRITE(IUNS,*) "------------------------------"
                S  = 1
                SS = (S*(S+1))/2
                RPS2 = (R+P+S)*(R+P+S+1)/2
                NS = 0 ! number of additional var cov parameters due to random SCALE
                NPAR = P+R+S ! number of parameters
                NPAR2 = NPAR*(NPAR+1)/2
                myqdim = 1
                mytotalq = totalqR0
                mypoints0 = pointsR0
                myweights0 = weightsR0
                ncov=0
            else if (CYCLES==2) THEN
                WRITE(IUN,*)
                WRITE(IUN,*)
                WRITE(IUN,'("---------------------------")')
                WRITE(IUN,'("Model WITH Scale Parameters")')
                WRITE(IUN,'("---------------------------")')
                WRITE(IUNS,*)
                WRITE(IUNS,'("---------------------------")')
                WRITE(IUNS,'("Model WITH Scale Parameters")')
                WRITE(IUNS,'("---------------------------")')
                WRITE(*,*)
                WRITE(*,*) "---------------------------"
                WRITE(*,*) "Model WITH Scale Parameters"
                WRITE(*,*) "---------------------------"
                S = S0
                SS = (S*(S+1))/2
                RPS2 = (R+P+S)*(R+P+S+1)/2
                tau(2:S) = 0
                NS = 0 ! number of additional var cov parameters due to random SCALE
                NPAR = P+R+S ! number of parameters
                NPAR2 = NPAR*(NPAR+1)/2
                myqdim = 1
                mytotalq = totalqR0
                mypoints0 = pointsR0
                myweights0 = weightsR0
                ncov=0
            else if (CYCLES==3) THEN
                WRITE(IUN,*)
                WRITE(IUN,*)
                WRITE(IUN,'("-----------------------")')
                WRITE(IUN,'("Model WITH RANDOM Scale")')
                WRITE(IUN,'("-----------------------")')
                WRITE(IUNS,*)
                WRITE(IUNS,'("-----------------------")')
                WRITE(IUNS,'("Model WITH RANDOM Scale")')
                WRITE(IUNS,'("-----------------------")')
                WRITE(*,*)
                WRITE(*,*) "-----------------------"
                WRITE(*,*) "Model WITH RANDOM Scale"
                WRITE(*,*) "-----------------------"
                S = S0
                SS = (S*(S+1))/2
                RPS2 = (r+P+S)*(r+P+S+1)/2
               ! NS = number of additional var cov parameters due to random SCALE
                ncov = 0
                NS = 1
                NS2 = NS*(NS+1)/2
                NPAR = P+r+S+NS ! number of parameters
                NPAR2 = NPAR*(NPAR+1)/2
                myqdim = 2
                mytotalq = totalqR1
                mypoints0 = pointsR1
                myweights0 = weightsR1
            else if (CYCLES==4) THEN
                S = S0
                SS = (S*(S+1))/2
                RPS2 = (r+P+S)*(r+P+S+1)/2
                ! NS = number of additional var cov parameters due to random SCALE
                ncov = ncov0
                NS = 1+ncov
                NS2 = NS*(NS+1)/2
                NPAR = P+r+S+NS ! number of parameters
                NPAR2 = NPAR*(NPAR+1)/2
                myqdim = 2
                mytotalq = totalqR1
                spar(ns) = spar(1)
                spar(1:(ns-1)) = 0
                mypoints0 = pointsR1
                myweights0 = weightsR1
            end if

        !
        ! start iterations
        !

             ! ifin = 1 for iterations using NR (or BHHH depending on the inversion of DER2)
             !      = 2 for final iteration
             !

             ridge  = ridgein
             RIDGEIT= 0

             ! set the ridge up for the models with scale parameters
             if (CYCLES<=2) then
                ridge = ridgein+.05D0
             END IF
             IF (CYCLES>=3) THEN
                ridge = RIDGEIN+.2D0
             END IF
             if(ridge < 0) ridge = ridgein

             loglp  = -999999999999999.0

            ITER=1
        ! START WITH NEWTON RAPHSON (MIXREG for starting values)
            IFIN=1
            IFINLOOP:DO WHILE (ifin < 3)

                 ! set ifin=2 if on the last iteration
                if (ifin == 2) then
                     ifin = 3
                 end if

                 ! put the ridge back to its initial value after the first 5 & 10 iterations
                IF (CYCLES>=3 .AND. ITER==11) THEN
                    ridge = ridge - .2
                 else if (CYCLES==2 .and. iter==6) then
                    ridge = ridge - .05
                END IF
                 if(ridge < 0) ridge = 0 !To get rid of -0 values
                 
            !
            ! calculate the derivatives and information matrix
            !
                 LOGL= 0.0D0
                 NN=0
                 DER2A(:)=0.0D0
                 DER2b(:)=0.0D0
                 DER(:)=0.0D0
                 LIK(:,:)=0.0D0
                 H(:)=0.0D0

                ILOOP:DO I=1,NC2  ! go over level-2 clusters
                     if (I > 1) THEN
                ! augment nn by the last level-2 ID's number of observations
                        nn = nn + IDNI(I-1,2)
                     end if
                     DERP(:)=0.0D0
                     DERP2(:)=0.0D0

                    NQ1=NQ

                   ! modify the points and weights for adaptive quadrature
            
                    mypoints = mypoints0
                    myweights = myweights0
                    if(aquad .ne. 0 .and. (iter > 8 .or. cycles .eq. 2 .or. cycles .eq. 4)) then
            !        if(aquad .ne. 0 .and. (iter >= 8 .or. cycles .eq. 4)) then
                         do k=1,myqdim
                             call getSDev(k,myqdim,thetavs(i,:), sdev)
                             do q=1,mytotalq
                                 temp = mypoints(q,k)
                                 mypoints(q,k) = thetas(i,k) + sdev*mypoints(q,k)
                                 call get_phi_ratio(mypoints(q,k),temp, phiRatio)
                                 myweights(q) = sdev*phiRatio*myweights(q)
                             end do
                         end do
                    end if
                   
                    QLOOP:DO Q=1, mytotalq  ! go over quadrature points
                           
                        PSUM     = 0.0D0
                        DERQ(:)  = 0.0D0
                        DERQ2(:) = 0.0D0
                        DZ(:)    = 0.0D0

                        JLOOP:DO J=1,IDNI(I,2)  ! loop over level-1 observations
                            NOB=NN+J
                            XB = DOT_PRODUCT(BETA,X(NOB,:))        ! X BETA for the current LEVEL-1 obs
                            UA = DOT_PRODUCT(ALPHA,U(NOB,:))       ! U ALPHA for the current LEVEL-1 obs
                            WT = DOT_PRODUCT(TAU(1:S),W(NOB,1:S))  ! W TAU for the current LEVEL-1 obs
                                                                   ! note that S changes over CYCLES
                            vaug = 0                                       
                            !vaug is the derivative vector for the s parameters
                            if(cycles >= 3) then
                                if(ncov>=1) vaug(1) = mypoints(q,1)
                                if(ncov==2) vaug(2) = vaug(1)*vaug(1)
                                vaug(ncov+1) = mypoints(q,2)
                            end if
                            !mytheta is the quadrature points vector for spar
!                            if(cycles >= 3) then
!                                mytheta = mypoints(q,2)
!                            end if
                            rtemp = wt
                            if(ns > 0) RTEMP = rtemp + spar(ncov+1)*mypoints(q,2)
                            if(ns > 1) rtemp = rtemp + spar(1)*mypoints(q,1)
                            if(ns > 2) rtemp = rtemp + spar(2)*mypoints(q,1)**2
                            IF (RTEMP .GE. LOGBIG) THEN
                              WSVAR = BIG
                            ELSE
                              WSVAR = DEXP(RTEMP)
                            END IF
                                                       
                            IF (UA .GE. LOGBIG) THEN
                                  BSVAR = BIG
                            ELSE 
                                  BSVAR = DEXP(UA)
                            END IF
                            IF (BSVAR .LE. SMALL) BSVAR = SMALL
                            ERRIJ = Y(NOB) - (XB + DSQRT(BSVAR)*mypoints(q,1))
                            IF (WSVAR .LE. SMALL) WSVAR = SMALL   
                            LPROB = -.5D0*(DLOG(2.0D0*PI)+DLOG(WSVAR) + ERRIJ**2/WSVAR)
                            PSUM  = PSUM + LPROB
                            ! GET FIRST DERIVATIVES
                            DZ(1:P) = (-2.0D0/WSVAR)*ERRIJ*X(NOB,:)                    ! beta
                            DZ(P+1:P+R) = (-ERRIJ/WSVAR)*mypoints(q,1)*DSQRT(BSVAR)*U(NOB,:) ! alpha
                            DZ(P+R+1:P+r+S) = (1.0D0 - (ERRIJ**2/WSVAR))*W(NOB,1:S)    ! tau
                            if(ns > 0) dz(p+r+s+1:p+r+s+ns) = (1.0D0 - (ERRIJ**2/WSVAR))*vaug(1:ns) ! sparam

            !write(IUNS,'(17F8.3)') i*1.0,q*1.0,j*1.0,mypoints(q,1),lprob,psum,errij,y(nob),xb,uchth,wt,(dz(k), k=1,npar)                
                            DERQ = DERQ + (-.5D0)*DZ

                            ! 2ND DERIVATIVE MATRIX  (NEWTON RAPHSON)

                               ! beta         (p x p symmetric)
                            CALL GRAM(X(NOB,:),BDER2,P,1)
                            BDER2 = (-1.0D0/wsvar)*BDER2

                            ! alpha, beta  (r rows & p columns)
                            CALL MPYM(U(nob,:),X(NOB,:),ABDER2,r,1,0,0,P)
                            ABDER2 = (-.5D0/wsvar)*mypoints(q,1)*DSQRT(bsvar)*ABDER2

                            ! alpha (r x r symmetric)
                            CALL GRAM(u(nob,:),ADER2,r,1)
                            ADER2 = (.25D0/wsvar)*mypoints(q,1)*DSQRT(bsvar)*(errij-mypoints(q,1)*DSQRT(bsvar))*ADER2

                            CALL ADJRC(BDER2,ABDER2,ADER2,DZ2(1:RP2),r,P)

                            ! tau, beta    (s rows & p columns)
                            CALL MPYM(W(NOB,1:S),X(NOB,:),TBDER2,S,1,0,0,P)
                            TBDER2 = (-1.0D0/wsvar)*errij*TBDER2

                            ! tau, alpha   (s rows & r columns)
                            CALL MPYM(W(NOB,1:S),u(nob,:),TADER2,S,1,0,0,r)
                            TADER2 = (-.5D0/wsvar)*errij*mypoints(q,1)*DSQRT(bsvar)*TADER2
                            ! tau
                            CALL GRAM(W(NOB,1:S),TDER2,S,1)
                            TDER2 = (-.5D0/wsvar)*errij*errij*TDER2

                            CALL ADJC(TBDER2,TADER2,TBADER2,S,P,R)
                            CALL ADJRC(DZ2(1:RP2),TBADER2,TDER2,DZ2(1:RPS2),S,P+R)

                               ! random scale
                            IF (NS > 0) THEN
                                ! scale, beta    (NS rows & p columns)
                                 CALL MPYM(vaug,X(NOB,:),SBDER2,NS,1,0,0,P)
                                  SBDER2 = (-1.0D0/wsvar)*errij*SBDER2
                                  ! scale, alpha   (NS rows & r columns)
                                  CALL MPYM(vaug,U(NOB,:),SADER2,NS,1,0,0,R)
                                  SADER2 = (-.5D0/wsvar)*errij*mypoints(q,1)*DSQRT(bsvar)*SADER2
                                  ! scale, tau     (NS rows & s columns)
                                  CALL MPYM(vaug,W(NOB,1:S),STDER2,NS,1,0,0,S)
                                  STDER2 = (-.5D0/wsvar)*errij*errij*STDER2
                                  ! scale          (2 x 2 symmetric)
                                  CALL GRAM(vaug,SDER2,NS,1)
                                  SDER2 = (-.5D0/wsvar)*errij*errij*SDER2
                                  CALL ADJC(SBDER2,SADER2,SBADER2,NS,P,R)
                                  CALL ADJC(SBADER2,STDER2,SBATDER2,NS,P+R,S)
                                  CALL ADJRC(DZ2(1:RPS2),SBATDER2,SDER2,DZ2(1:NPAR2),NS,P+R+S)
                            ENDIF
                            DERQ2 = DERQ2 + DZ2

                        END DO JLOOP

                        IF (PSUM .GE. LOGBIG) THEN
                             LIK(I,Q) = BIG
                        ELSE 
                             LIK(I,Q) = DEXP(PSUM)
                        END IF

                        rtemp2 = myweights(q)
                        IF (RTEMP2 .LE. SMALL) RTEMP2 = SMALL            
                        PRA   = DEXP(PSUM + DLOG(RTEMP2))
                        H(I)     = H(I) + PRA
                        DERP     = DERP + DERQ*PRA

                        CALL GRAM(DERQ,DERQQ,NPAR,1)
                        DERP2 = DERP2 + ((DERQQ + DERQ2)*PRA)

            !if(ns > 0) write(IUNS,'(18F10.3)') i*1.0, q*1.0, psum, rtemp2,pra,h(i),(derq(k), k=1,npar),derp2(npar2),derq2(npar2)
            !if(ns == 0) write(IUNS,'(18F10.3)') i*1.0, q*1.0, psum, rtemp2,pra,h(i),(derq(k), k=1,npar)
                    END DO QLOOP

                       ! calculate the empirical Bayes estimate
                    if(h(i) .le. small) h(i) = small
                    theta = 0
                    thetav = 0
                    do k=1, myqdim
                        do q=1,mytotalq
                               theta(k) = theta(k) + lik(i,q)*myweights(q)*mypoints(q,k)
                        end do
                    end do
                    theta = theta/h(i)
                    do q=1,mytotalq
                        temp = lik(i,q)*myweights(q)
                        do k=1, myqdim
                           work2(k) = mypoints(q,k) - theta(k)
                        end do
                        call grmcv(thetav, thetav, work2, temp, myqdim)
                    end do
                    thetav = thetav/h(i)

                    thetas(i,:) = theta
                    thetavs(i,:) = thetav

                    LOGL  = LOGL + DLOG(H(I))
                    SCAL  = DEXP(0.0D0 - DLOG(H(I)))
                    DERPS = SCAL*DERP
                    DER   = DER + DERPS
                    CALL GRMCV(DER2A,DER2A,DERPS,ONE,NPAR)
                    DER2B = DER2B + SCAL*DERP2
                END DO ILOOP

                 LOGDIFF = LOGL-LOGLP
                 LOGLP   = LOGL

                 ! determine if an NR iteration is bad and increase the ridge
                 ! take the ridge off after 10 good iterations
     IF (LOGDIFF/LOGLP > .005 .AND. ITER < MAXIT) THEN
        RIDGEIT = 0
        RIDGE = RIDGE + .1D0
        WRITE(IUN,'("==> BAD NR ITERATION ",I5," with NEW ridge = ",F8.4,/)') ITER,RIDGE
        GO TO 99
     END IF
     IF (LOGDIFF/LOGLP <= .000001 .AND. RIDGEIT < 10) THEN
        RIDGEIT = RIDGEIT+1
     ELSE IF (LOGDIFF/LOGLP <= .000001 .AND. RIDGEIT >= 10 .and. ifin==1) then
        ridge = ridgein
!        if(ridge < ridgein) RIDGE = ridgein
     END IF
!                 IF (LOGDIFF/LOGLP > .05 .AND. ITER < MAXIT) THEN
!                    RIDGEIT = 0
!                    RIDGE = RIDGE + .05D0
!                    WRITE(IUN,'("==> BAD NR ITERATION ",I5," with NEW ridge = ",F8.4,/)') ITER,RIDGE
!                    GO TO 99
!                 END IF
!                 IF (LOGDIFF/LOGLP <= .000001 .AND. RIDGEIT < 10) THEN
!                    RIDGEIT = RIDGEIT+1
!                 ELSE IF (LOGDIFF/LOGLP <= .000001 .AND. RIDGEIT >= 10 .and. ifin==1) then
!                        ridge = ridge - .05D0
!                        if ( ridge < RIDGEIN) ridge = RIDGEIN
!                 END IF
                 if(maxder < 1 .and. iter > 10) ridge = 0

                if(ifin < 2) then
                     ! ridge adjustment - diagonal elements only
                     DO L=1,NPAR
                        LL=(L*(L+1))/2
                        DER2A(LL)=DER2A(LL)+RIDGE*DER2A(LL)
                     END DO
                     DER2B = DER2A - DER2B
                     DO L=1,NPAR
                        LL=(L*(L+1))/2
                        DER2B(LL)=DER2B(LL)+RIDGE*DER2B(LL)
                     END DO
                else
                    DER2B = DER2A - DER2B
                end if

                write(IUNS,*)"2nd Derivatives"
                ll = 0
                do k=1,npar
                    write(IUNS,'(25f11.2)') (der2b(ll+l), l=1,k)
                    ll = ll + k
                end do

                 ! INVERT 2nd DERIVATIVE MATRIX - COREC is a working vector for INVS
                 IF (IFIN == 1) THEN
                    CALL INVS(DER2B,NPAR,DET,COREC,IER)
                    IF (IER == 0) THEN
                        WRITE(*,'(" Newton-Raphson Iteration ",I5," with ridge ",F8.4)') ITER,RIDGE
                        WRITE(IUNS,'(" Newton-Raphson Iteration ",I5," with ridge ",F8.4)') ITER,RIDGE
                        WRITE(18,'(" Newton-Raphson Iteration ",I5," with ridge ",F8.4)') ITER,RIDGE
                       CALL MPYM(DER2B,DER,COREC,NPAR,NPAR,1,0,1)
                    ELSE IF (IER == 1) THEN
                       WRITE(*,'(" BHHH Iteration ",I5," with ridge ",F8.4)') ITER,RIDGE
                       WRITE(IUNS,'(" BHHH Iteration ",I5," with ridge ",F8.4)') ITER,RIDGE
                       WRITE(18,'(" BHHH Iteration ",I5," with ridge ",F8.4)') ITER,RIDGE
                       CALL INVS(DER2A,NPAR,DET,COREC,IER)
                       CALL MPYM(DER2A,DER,COREC,NPAR,NPAR,1,0,1)
                    END IF
                 ELSE
                    CALL INVS(DER2B,NPAR,DET,COREC,IER)
                    CALL MPYM(DER2B,DER,COREC,NPAR,NPAR,1,0,1)
                    DO L=1,NPAR
                       LL=(L*(L+1))/2
                       SE(L)=DSQRT(DABS(DER2B(LL)))
!                       if(der2b(ll) < 0) then
!                            if(der(l) < 0 .and. corec(l) < 0) corec(l) = -corec(l)
!                            if(der(l) > 0 .and. corec(l) > 0) corec(l) = -corec(l)
!                            corec(l) = -corec(l)
!                        end if
                    END DO
                 END IF
                if(iter<=5 .and. cycles > 1) corec = corec*.5

                write(IUNS,*)"Corrections"
                write(IUNS,'(20f11.3)') (corec(k), k=1,npar)
                write(IUNS,*)"Derivatives"
                write(IUNS,'(20f11.3)') (der(k), k=1,npar)
                write(IUNS,*)"Beta"
                write(IUNS,'(20f11.3)') (beta(k), k=1,p)
                write(IUNS,*)"Alpha"
                write(IUNS,'(20f11.3)') (alpha(k), k=1,r)
                write(IUNS,*)"Tau"
                write(IUNS,'(20f11.3)') (tau(k),k=1,s)
                if(ns > 0) then
                    write(IUNS,*)"Spar"
                    write(IUNS,'(20f11.3)') (spar(k),k=1,ns)
                end if

                 MAXDER=MAXVAL(ABS(DER))
                 MAXCORR=MAXVAL(ABS(COREC))
                 WRITE(*,*) '  maximum correction and derivative'
                 WRITE(*,*) '  ',maxcorr,MAXDER
                  WRITE(IUNS,*) '  maximum correction and derivative'
                 WRITE(IUNS,*) '  ',MAXCORR,MAXDER
                  WRITE(18,*) '  maximum correction and derivative'
                 WRITE(18,*) '  ',MAXCORR,MAXDER
                 ! For MODELS WITH SCALE, do step-halving for first FIVE iterations */
            !     IF ((ns > 0) .AND. ITER < 6) THEN
            !        COREC = 0.5D0*COREC
            !     END IF


                 ! done with NR and onto last iteration
                 IF (IFIN==1 .AND. (MAXCORR <= CONV .OR. ITER >= MAXIT)) THEN
                     IFIN=2
                     ORIDGE=RIDGE
                     RIDGE=0.0D0
                 END IF

                 ! UPDATE PARAMETERS
                 BETA  = BETA  + COREC(1:P)
                 if(iter >= 10 .or. cycles == 1) alpha = alpha + corec(p+1:p+r)
                do k=1,ns
                    spar(k) = spar(k) + corec(npar-ns+k)
                end do
                 IF (S==1) THEN 
                     TAU(1)= TAU(1)   + COREC(P+R+1)
                 ELSE IF (S>1) THEN
                     TAU   = TAU      + COREC(P+R+1:P+R+S)
                 END IF 

                

            99  WRITE(*,'("   -2 Log-Likelihood = ",F14.5)') -2*LOGL
                 WRITE(IUNS,'("   -2 Log-Likelihood = ",F14.5)') -2*LOGL
                 WRITE(18,'("   -2 Log-Likelihood = ",F14.5)') -2*LOGL
                 ITER = ITER+1
            END DO IFINLOOP
                 
            if((cycles .ne. 3) .or. (cycles == 3 .and. ncov0 == 0)) then
                 ! WRITE RESULTS
               WRITE(IUN,562)ITER-1,ORIDGE,LOGL,LOGL-NPAR, &
               LOGL-0.5D0*DBLE(NPAR)*DLOG(DBLE(NC2)),0-2*LOGL,0-2*(LOGL-NPAR), &
               0-2*(LOGL-0.5D0*DBLE(NPAR)*DLOG(DBLE(NC2)))
               562 FORMAT(1X,'Total  Iterations =',I4,/, &
                          1X,'Final Ridge value =',F4.1,//, &
                     1X,'Log Likelihood                 = ',F12.3,/, &
                     1X,"Akaike's Information Criterion = ",F12.3,/, &
                     1X,"Schwarz's Bayesian Criterion   = ",F12.3,//,&
                     1X,"==> multiplied by -2             ",      /  &
                     1X,'Log Likelihood                 = ',F12.3,/, &
                     1X,"Akaike's Information Criterion = ",F12.3,/, &
                     1X,"Schwarz's Bayesian Criterion   = ",F12.3,/)
                WRITE(IUN,57)
 57 FORMAT(/,'Variable',12x,'    Estimate',4X,'AsymStdError',4x, &
          '     z-value',4X,'     p-value',/,'----------------',4x,  &
          '------------',4X,'------------',4X,'------------',4X,'------------')
            
                 PVAL=0.0D0
                 ZVAL=0.0D0
            
                 IF (NCENT==1 .AND. P>1) THEN
                     WRITE(IUN,'("STANDARDIZED BETA (regression coefficients)")')
                 ELSE
                     WRITE(IUN,'("BETA (regression coefficients)")')
                 END IF
                 DO L=1,P
                    ZVAL = BETA(L)/SE(L)
                    PVAL = 2.0D0 *(1.0D0 - PHIFN(DABS(ZVAL),0))
                        WRITE(IUN,804)BLAB(L),BETA(L),SE(L),ZVAL,PVAL
                 END DO
            
                 IF (NCENT==1 .AND. R>1) THEN
                     WRITE(IUN,'("STANDARDIZED ALPHA (BS variance parameters: log-linear model)")')
                 ELSE
                     WRITE(IUN,'("ALPHA (BS variance parameters: log-linear model)")')
                 END IF
                 DO L=1,R
                    L2 = P+L
                    ZVAL = ALPHA(L)/SE(L2)
                    PVAL = 2.0D0 *(1.0D0 - PHIFN(DABS(ZVAL),0))
                        WRITE(IUN,804)ALAB(L),ALPHA(L),SE(L2),ZVAL,PVAL
                 END DO
                 IF (NCENT==1 .AND. S>1) THEN 
                     WRITE(IUN,'("STANDARDIZED TAU (WS variance parameters: log-linear model)")')
                 ELSE
                     WRITE(IUN,'("TAU (WS variance parameters: log-linear model)")')
                 END IF 
                 DO L=1,S
                    L2 = P+R+L
                    ZVAL = TAU(L)/SE(L2)
                    PVAL = 2.0D0 *(1.0D0 - PHIFN(DABS(ZVAL),0))
                        WRITE(IUN,804)TLAB(L),TAU(L),SE(L2),ZVAL,PVAL
                 END DO
            
                 IF (cycles==4 .or. (cycles==3 .and. ncov0==0)) THEN
                    WRITE(IUN,'("Random scale standard deviation")')
                    L2=P+R+S+ncov+1
                    if(spar(ncov+1) < 0) spar(ncov+1)=abs(spar(ncov+1))
                    ZVAL = SPAR(ncov+1)/SE(L2)
                    PVAL = 2.0D0 *(1.0D0 - PHIFN(DABS(ZVAL),0))
                    WRITE(IUN,804)'Std Dev         ',SPAR(ncov+1),SE(L2),ZVAL,PVAL
                 end if
                 IF (NS>0) THEN
                    if(ncov >= 1) then
                        L2=P+R+S+1
                        ZVAL = SPAR(1)/SE(L2)
                        PVAL = 2.0D0 *(1.0D0 - PHIFN(DABS(ZVAL),0))
                        WRITE(IUN,'("Random location (mean) effect on WS variance")')
                        WRITE(IUN,804)'Loc  Eff        ',SPAR(1),SE(L2),ZVAL,PVAL
                    end if
                    if(ncov==2) then
                        WRITE(IUN,'("Random quadratic location (mean) effects on WS variance")')
                        L2=P+R+S+2
                        ZVAL = SPAR(2)/SE(L2)
                        PVAL = 2.0D0 *(1.0D0 - PHIFN(DABS(ZVAL),0))
                        WRITE(IUN,804)'Quad Loc        ',SPAR(2),SE(L2),ZVAL,PVAL
                    end if
                end if
             804 FORMAT(A16,4(4x,F12.5))
            
                         IF (MAXDER > CONV .AND. ITER >= MAXIT) THEN
                            WRITE(IUN,'("NOTE: CONVERGENCE CRITERION WAS NOT ACHIEVED")')
                            WRITE(IUN,'("FINAL FIRST DERIVATIVE AND CORRECTION VALUES")')
                            DO L=1,NPAR
                              WRITE(IUN,*)DER(L),COREC(L)
                            END DO
                         END IF
    
                 ! write out the deviance, estimates, and standard errors
                 IF (CYCLES.EQ.1) OPEN(2, FILE='mixREGLS5.EST')
                    WRITE(2,'(F15.6,2I8)') -2*LOGL,ITER-1,MAXIT
                     WRITE(2,'(5F15.8)')(BETA(L),L=1,P)
                     WRITE(2,'(5F15.8)')(alpha(L),L=1,R)
                     WRITE(2,'(5F15.8)')(TAU(L),L=1,S)
                     write(2,'(5F15.8)')(SPAR(L),L=1,ns)
                     WRITE(2,'(5F15.8)')(SE(L),L=1,NPAR)
            end if
        END DO CYCLELOOP
     write(iun,*)
     write(iun,*)
         myz = 1.959964
     WRITE(IUN,'("BS variance ratios and 95% CIs")')
     write(iun,'("------------------------------")')
     write(iun,*)
    WRITE(IUN,808) 'Variable        ','Ratio','Lower','Upper'
    write(iun,808)'---------------------','------------------','------------','------------'
         WRITE(IUN,'("ALPHA (BS variance parameters: log-linear model)")')
     DO L=1,R
        L2 = P+L
        tauhat = exp(alpha(l))
        tauhatlow = exp(alpha(l)-myz*se(l2))
        tauhatup = exp(alpha(l)+myz*se(l2))
        WRITE(IUN,804)aLAB(L),tauhat, tauhatlow, tauhatup
     END DO
     write(iun,*)
     write(iun,*)
808 FORMAT(A16,3(4x,A12))
     WRITE(IUN,'("WS variance ratios and 95% CIs")')
     write(iun,'("------------------------------")')
     write(iun,*)
    WRITE(IUN,808) 'Variable        ','Ratio','Lower','Upper'
    write(iun,808)'---------------------','------------------','------------','------------'

         WRITE(IUN,'("TAU (WS variance parameters: log-linear model)")')
     DO L=1,S
        L2 = P+R+L
        tauhat = exp(tau(l))
        tauhatlow = exp(tau(l)-myz*se(l2))
        tauhatup = exp(tau(l)+myz*se(l2))
        WRITE(IUN,804)TLAB(L),tauhat, tauhatlow, tauhatup
     END DO
         if(ncov > 0) then
        WRITE(IUN,'("Random location (mean) effect on WS variance")')
             l2 = p+r+s+1
            tauhat = exp(spar(1))
            tauhatlow = exp(spar(1)-myz*se(l2))
            tauhatup = exp(spar(1)+myz*se(l2))
            WRITE(IUN,804)"Location Effect ",tauhat, tauhatlow, tauhatup
            if(ncov > 1) then
                 l2 = p+r+s+2
                tauhat = exp(spar(2))
                tauhatlow = exp(spar(2)-myz*se(l2))
                tauhatup = exp(spar(2)+myz*se(l2))
                WRITE(IUN,804)"Quad Location   ",tauhat, tauhatlow, tauhatup
            end if
        end if
    write(iun,*)
        WRITE(IUN,'("Random scale standard deviation")')
        L2=P+R+S+ns
        tauhat = exp(spar(ns))
        tauhatlow = exp(spar(ns)-myz*se(l2))
        tauhatup = exp(spar(ns)+myz*se(l2))
        WRITE(IUN,804)'Std Dev                 ',tauhat, tauhatlow, tauhatup
         
           
         CLOSE(IUN)
         CLOSE(IUNS)
         CLOSE(1)
         CLOSE(2)
         CLOSE(3)
         close(18)
    open(23,file=trim(fileprefix)//"_ebvar.dat")
    k=myqdim*(myqdim+1)/2
    DO I=1,NC2  ! go over level-2 clusters
        write(23,'(20F15.8)') (THETAs(I,j), j=1,myqdim), (thetavs(i,j), j=1,k)
    end do
    close(23)
!    open(29,file=trim(fileprefix)//"_mvinfo.dat")
!    write(29,*) nobs, maxk, nc2, 2, 1, ns
!    write(29,*) (spar(j),j=1,ns)
!    nob = 0
!    do i=1,nc2
!        write(29,*) idni(i,2), h(i)
!        do j=1,IDNI(I,2)  ! loop over level-1 observations
!           nob=nob+1 
!                            XB = DOT_PRODUCT(BETA,X(NOB,:))        ! X BETA for the current LEVEL-1 obs
!                            UA = DOT_PRODUCT(ALPHA,U(NOB,:))       ! U ALPHA for the current LEVEL-1 obs
!                            WT = DOT_PRODUCT(TAU(1:S),W(NOB,1:S))  ! W TAU for the current LEVEL-1 obs
!                                                                   ! note that S changes over CYCLES
!            write(29,*) y(nob)-xb, exp(.5*ua), wt
!        end do
!    end do
!    close(29)

        deallocate(weightsR1,pointsR1,weightsR0,pointsR0,myweights,mypoints)
        deallocate(bder2,abder2,ader2,der2a,der2b,der,derp,derp2,derq,derq2,derqq,dz,derps,tader2, &
                    tbder2,tbader2,tder2,dz2,sbder2,sader2,stder2,sbader2,sbatder2,sder2)
        deallocate(vaug,lik,h,corec,se,theta,thetav,work2,work3)
    END SUBROUTINE mixREGLSEST


END PROGRAM mixREGLS_subject


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
SUBROUTINE INVS(A,N,C,W,IER)
   DOUBLE PRECISION U,X,Y,Z,D,A,C,W
   DIMENSION A(1), W(1)                                             
   INTEGER DIAGMK                                                    
   INTEGER DIAG,DIAG2,ROWNO,ROWCOL                                   
   INTEGER COLNO, IER
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
!      IF(PRESENT(VERBOSE)) THEN
!         IF(VERBOSE) CALL POST_ERROR('ERROR MESSAGE FROM SUBROUTINE INVS: MATRIX IS SINGULAR')
!          WRITE(*,*) 'ERROR MESSAGE FROM SUBROUTINE INVS: MATRIX IS SINGULAR'
!      ELSE
!         CALL POST_ERROR('ERROR MESSAGE FROM SUBROUTINE INVS: MATRIX IS SINGULAR')
!          WRITE(*,*) 'ERROR MESSAGE FROM SUBROUTINE INVS: MATRIX IS SINGULAR'
!      ENDIF
!      IF (PRESENT(IER)) THEN
         IER = 1
!      END IF
   ENDIF
   RETURN                                                            
END SUBROUTINE INVS


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



!
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
      DOUBLE PRECISION A,B,C,X
      DIMENSION A(1),B(1),C(1)
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


! ************************************************
! SUBROUTINE PHIFN
! Calculate the probability distribution function (Intercept) for
! various distributions:
! NORMAL, LOGISTIC, Complementary Log-Log, OR Log-Log 
!   = 0       = 1                    = 2          =3 
! ************************************************
REAL*8 FUNCTION PHIFN(Z,nfn)
   IMPLICIT REAL*8(A-H,O-Z)
   
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
!         CALL POST_ERROR('ERROR MESSAGE FROM SUBROUTINE CHAMS: ' // &
!                'ILLICIT COMBINATION OF STORAGE MODES')
      END SELECT
   END DO
   
END SUBROUTINE CHAMS

SUBROUTINE CHAMS0to1(A,B,N)                                   
   DOUBLE PRECISION A,B
   DIMENSION A(N,N),B(N*(N+1)/2)
   
         K = 0                                                             
         DO J=1,N                                                         
            DO I=1,J              
               K = K + 1         
               B(K) = A(i,j)  
            END DO
         END DO

END SUBROUTINE CHAMS0to1

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

