program mixregls_both
    use lsboth,only:mls,fileout,stage2
    implicit none
    call readdef_both()
    call writedef_both()
    call readat()
    call adjustdata()
    call printdesc()
    call callmixreg()

    if(mls .eq. 1) then
        call mixregmlsest()
#if defined(_WIN32)
        CALL SYSTEM("COPY mixregls_both1.OUT+mixregls_both3.OUT+mixregls_both2.OUT " // FILEOUT)
#else
        CALL SYSTEM("cat mixregls_both1.OUT mixregls_both3.OUT mixregls_both2.OUT >> " // FILEOUT)
#endif
        call mixreglsest()
#if defined(_WIN32)
        CALL SYSTEM("COPY mixregls_both1.OUT+mixregls_both2.OUT " // FILEOUT)
#else
        CALL SYSTEM("cat mixregls_both1.OUT mixregls_both2.OUT >> " // FILEOUT)
#endif
    end if
    if(stage2 .ne. 0) then
        call run_stage2()
    end if
#if defined(_WIN32)
    CALL SYSTEM("DEL mixregls_both1.OUT mixregls_both2.OUT mixreg.est mixreg.var mixreg.def mixreg.lik")
    if(mls .eq. 1) call system("del mixregls_both3.OUT")
#else
    CALL SYSTEM("rm mixregls_both1.OUT")
    CALL SYSTEM("rm mixregls_both2.OUT")
    CALL SYSTEM("rm mixreg.est")
    CALL SYSTEM("rm mixreg.var")
    CALL SYSTEM("rm mixreg.def")
    CALL SYSTEM("rm mixreg.lik")
    if(mls .eq. 1) call system("rm mixregls_both3.OUT")
#endif

#if defined(_WIN32)
    call system("mkdir work")
    call system("move mixregls_both_* work")
#else
    call system("mkdir work")
    call system("mv mixregls_both_* work")
#endif
end program mixregls_both

subroutine readdef_both()
    use lsboth
    use procedures
    implicit none
    INTEGER :: I,j,k

    OPEN(1, FILE='lsboth_random_mixblank.def')
    READ(1,'(18A4)') HEAD
    READ(1,'(A80)')FILEDAT
    READ(1,'(A80)')FILEprefix
    filedef = trim(fileprefix)//".def"
    fileout = trim(fileprefix)//".out"
    filedat = adjustl(filedat)

    READ(1,*) NVAR,P,R,S,PNINT,RNINT,SNINT,pv,rv, sv, CONV, NQ, AQUAD, MAXIT, yMISS, NCENT, ncov, &
        RIDGEIN, discard0, mls, chol, nreps, cutoff, nors, myseed, stage2, multi2nd, sepfile

    if(mls .ne. 1) mls = 0
    if(ncov > 2) ncov = 2
    if(ncov < 0) ncov = 0
    if(mls .eq. 1) ncov = 1
    if(nors .ne. 1) nors = 0
    if(multi2nd .ne. 1) multi2nd = 0
    if(sepfile .ne. 1) sepfile = 0
! set SNINT=0 (for the error variance) if SNINT=1 AND S=0
! DO NOT ALLOW A MODEL TO BE FIT WITHOUT AN ERROR VARIANCE
    IF (S==0 .AND. SNINT==1) SNINT=0
    miss = 1
    if(fp_equal(ymiss, 0.0d0)) miss = 0

!   SCALEP=0.50D0
! nvar     =  number of variables
! ridgein  =  initial value of the ridge
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

    p = pold + 1 - pnint + 2*pv
    r = rold + 1 - rnint + 2*rv
    s = sold + 1 - snint + 2*sv

    !Special case, random slope term requires occasion-level variables
    if(mls .eq. 1) r = rold + 1 - rnint + rv

! read in the labels
     READ(1,*) YLABEL
     IF (P .GE. 1) THEN
        ALLOCATE(BLAB(P))
        if(pold > 0) READ(1,*) (BLAB(I+1-pnint), I=1,Pold)
        if(pnint .ne. 1) blab(1) = "intercept"
     END IF
     IF (R .GE. 1) THEN
        ALLOCATE(ALAB(R))
        if(rold > 0) READ(1,*) (ALAB(I+1-rnint), I=1,Rold)
        if(rnint .ne. 1) alab(1) = "intercept"
     END IF
     IF (S .GE. 1) THEN
        ALLOCATE(TLAB(S))
        if(sold > 0) READ(1,*) (TLAB(I+1-snint), I=1,Sold)
        if(snint .ne. 1) tlab(1) = "intercept"
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
     ns = 0
        rr = (r+1)*r/2
     if(nors .ne. 1) then
        ns = 1+ncov
        if(mls .eq. 1) ns = 1+ncov*r
    end if
    npar = p+r+s+ns
    npar2 = npar*(npar+1)/2
    numloc = max(r*mls,1)
    ndim = numloc + 1
    ndim2 = ndim*(ndim+1)/2
    NQwR1 = NQ**ndim
    nqwR0 = nq**numloc
    select case(stage2)
        case(1,3)
            readcats = 0
        case(2,4)
            readcats = 1
        case default
            stage2 = 0
            nvar2 = 0
            nvar3 = 0
    end select
    if(stage2 .ne. 0) then
        if(nors .eq. 1) then
            pomega = -1
            pto = -1
            read(1,*) pfixed, ptheta
        else
            read(1,*) pfixed,ptheta,pomega,pto
        end if
        nvar2 = 1+max(pfixed,0)+max(pomega,0)+max(ptheta,0)+max(pto,0)
        nvar3 = 3+pfixed+pomega+numloc*(ptheta+1+pto+1)

        allocate(var2ind(nvar2))
        allocate(var2label(nvar2))
        if(readcats .eq. 1) then
            read(1,*) maxj
            allocate(icode(maxj))
            READ(1,*)(ICODE(J), J = 1,MAXJ)
        end if
        if(sepfile .eq. 1) read(1,*) filedat2
        if(sepfile .eq. 1) read(1,*) nvarsep,id2indsep
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
    end if
    CLOSE(1)

end subroutine readdef_both

subroutine writedef_both
    use lsboth
    implicit none
    INTEGER :: I,j,k
    OPEN(1,FILE=FILEDEF)
    WRITE(1,'(18A4)') HEAD
    WRITE(1,'(A80)') FILEDAT
    WRITE(1,'(A80)') FILEprefix
    WRITE(1,'(10I3,E10.1E3, i4, i2, i5, f12.3, 2i2, f6.3,3i2,i5,f6.3,i2,i12,3i2)') NVAR, Pold, Rold, Sold, PNINT, RNINT, SNINT, &
        pv,rv,sv, CONV, NQ, AQUAD, MAXIT,yMISS,NCENT,ncov,ridgein,discard0,MLS,chol,nreps,cutoff,nors,myseed,stage2,multi2nd,sepfile
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

    if(stage2 .ne. 0) then
        if(nors .eq. 1) then
            write(1,'(2i3)') pfixed, ptheta
        else
            write(1,'(4i3)') pfixed,ptheta,pomega,pto
        end if
        if(readcats .eq. 1) then
            write(1,*) maxj
            write(1,*)(ICODE(J), J = 1,MAXJ)
        end if
        if(sepfile .eq. 1) write(1,*) filedat2
        if(sepfile .eq. 1) write(1,*) nvarsep, id2indsep
        write(1,'(20I3)') var2ind(1)
        k = 1
        IF (Pfixed .GE. 1) THEN
            write(1,'(20I3)') (var2ind(k+I), I=1,Pfixed)
            k = k + pfixed
         END IF
        IF (Ptheta .GE. 1) THEN
            write(1,'(20I3)') (var2IND(k+I), I=1,Ptheta)
            k = k + ptheta
        END IF
        IF (Pomega .GE. 1) THEN
            write(1,'(20I3)') (var2IND(k+I), I=1,Pomega)
            k = k + pomega
         END IF
        IF (Pto .GE. 1) THEN
            write(1,'(20I3)') (var2IND(k+I), I=1,Pto)
         END IF

         write(1,*) var2label(1)
        k = 1
        IF (Pfixed .GE. 1) THEN
            write(1,*) (var2label(k+I), I=1,Pfixed)
            k = k + pfixed
         END IF
        IF (Ptheta .GE. 1) THEN
            write(1,*) (var2label(k+I), I=1,Ptheta)
            k = k + ptheta
        END IF
        IF (Pomega .GE. 1) THEN
            write(1,*) (var2label(k+I), I=1,Pomega)
            k = k + pomega
         END IF
        IF (Pto .GE. 1) THEN
            write(1,*) (var2label(k+I), I=1,Pto)
        END IF
    end if
    CLOSE(1)
end subroutine writedef_both

SUBROUTINE READAT()
    use lsboth
    use procedures
    implicit none

    INTEGER :: myPASS,I,K,ICOUNT,myindex,IDTEMP,IDOLD,hasmiss,nvartotal,discardi,m,j
    real(kind=8):: firsty
    REAL(KIND=8),ALLOCATABLE:: TEMPR(:)
    INTEGER,ALLOCATABLE :: allvarsind(:)
    LOGICAL FIRST

    ALLOCATE (TEMPR(NVAR))
        nvarTotal = 1+pold+rold+sold+nv+nvar2
        if(sepfile .eq. 1) nvarTotal = 1+pold+rold+sold+nv
        allocate(allvarsIND(nvarTotal))
        allvarsIND(1) = yind
        allvarsIND(2:pold+1) = xind(1:pold)
        allvarsIND(pold+2:pold+rold+1) = uind(1:rold)
        allvarsIND(pold+rold+2:pold+rold+sold+1) = wind(1:sold)
        allvarsIND(pold+rold+sold+2:pold+rold+sold+1+nv) = varind(1:nv)
        if(sepfile .ne. 1) allvarsIND(pold+rold+sold+2+nv:nvarTotal) = var2ind(1:nvar2)
    if(discard0 .ne. 0) open(16, file=trim(fileprefix)//"_removed.dat")
    num0 = 0

   ! INITIALIZE
    DO myPASS = 1,2
        IF (myPASS .EQ. 2) THEN
               ALLOCATE (Y(ICOUNT))
                ALLOCATE (X(ICOUNT,P))
                ALLOCATE (U(ICOUNT,R))
                ALLOCATE (W(ICOUNT,S))
                allocate (var(icount,nv))
                allocate (varavg(nc2,nv))
                allocate(ids(icount))
                if(sepfile .ne. 1) then
                    allocate(tempsums(nc2,nvar2))
                    allocate(data2(icount,nvar2))
                    tempsums = 0
                    data2 = 0
                end if

              Y = 0.0D0
              X = 0.0D0
              U = 0.0D0
              W = 0.0D0
              var = 0.0D0
              varavg = 0
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

      ! QUERY FOR NEW ID AND SET PARAMETERS ACCORDINGLY

              IDTEMP = INT(TEMPR(ID2IND))
              ! QUERY FOR NEW ID AND SET PARAMETERS ACCORDINGLY

              IF (.NOT. FIRST) THEN
                  ! if r=0 and rnint=1 then NO random effects
                 IF (R .GE. 1 .AND. IDTEMP .EQ. IDOLD) THEN
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
                        icount = icount - k
                        write(*,*) "REMOVED"
                           num0 = num0 + 1
                        do j=1,k
                            if(sepfile .ne. 1) then
                                write(16,'(i9,16f10.3)') ids(icount+j), y(icount+j), (x(icount+j,m),m=1,pold), &
                                                    (u(icount+j,m),m=1,rold), (w(icount+j,m), m=1,sold), &
                                        (var(icount+j,m), m=1,nv), (data2(i,m),m=1,nvar2)
                            else
                                write(16,'(i9,16f10.3)') ids(icount+j), y(icount+j), (x(icount+j,m),m=1,pold), &
                                                    (u(icount+j,m),m=1,rold), &
                                        (w(icount+j,m), m=1,sold), (var(icount+j,m), m=1,nv)
                            end if
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
                    ids(icount) = idtemp
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
                    if(sepfile .ne. 1) then
                        do myindex=1,nvar2
                            if(var2ind(myindex) .ne. 0) then
                                tempsums(i,myindex) = tempsums(i,myindex) + tempr(var2ind(myindex))
                                data2(icount,myindex) = tempr(var2ind(myindex))
                            end if
                        end do
                    end if
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
                        icount = icount - k
                            write(*,*) "REMOVED"
                           num0 = num0 + 1
                        do j=1,k
                            if(sepfile .ne. 1) then
                                write(16,'(i9,16f10.3)') ids(icount+j), y(icount+j), (x(icount+j,m),m=1,pold), &
                                                    (u(icount+j,m),m=1,rold), (w(icount+j,m), m=1,sold), &
                                        (var(icount+j,m), m=1,nv), (data2(i,m),m=1,nvar2)
                            else
                                write(16,'(i9,16f10.3)') ids(icount+j), y(icount+j), (x(icount+j,m),m=1,pold), &
                                                    (u(icount+j,m),m=1,rold), &
                                        (w(icount+j,m), m=1,sold), (var(icount+j,m), m=1,nv)
                            end if
                        end do
                        i = i -1
                    end if

    NC2 = I
    CLOSE(1)
   END DO   ! two passes, one to get size, second to read data
        do i=1,nc2
            if(sepfile .ne. 1) tempsums(i,1:nvar2) = tempsums(i,1:nvar2)/idni(i,2)
            varavg(i,1:nv) = varavg(i,1:nv)/idni(i,2)
        end do
   DEALLOCATE(TEMPR)
END SUBROUTINE READAT

subroutine adjustdata()
    use lsboth
    use procedures
    implicit none
    integer:: j,ll,h,ko,i,kv
    if(pnint .ne. 1) x(:,1) = 1
    if(rnint .ne. 1) u(:,1) = 1
    if(snint .ne. 1) w(:,1) = 1

    do j=1, pv
        ll = pold+1-pnint
        blab(ll+j*2-1) = trim(varlabel(j)) // "_BS"
        blab(ll+j*2) = trim(varlabel(j)) // "_WS"
    end do
    do j=1, rv
        ll = rold + 1 - rnint
        if(mls .eq. 1) then
            alab(ll+j) = trim(varlabel(j+pv)) // "_WS"
        else
            alab(ll+j*2-1) = trim(varlabel(j+pv)) // "_BS"
            alab(ll+j*2) = trim(varlabel(j+pv)) // "_WS"
        end if
    end do
    do j=1, sv
        ll = sold+1-snint
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
                if(mls .eq. 1) then
                    u(ko,ll+j) = var(ko,j+kv) - varavg(i,j+kv)
                else
                    u(ko,ll+j*2-1) = varavg(i,j+kv)
                    u(ko,ll+j*2) = var(ko,j+kv) - varavg(i,j+kv)
                end if
            end do
            kv = kv + rv
            do j=1, sv
                ll = sold + 1 - snint
                w(ko,ll+j*2-1) = varavg(i,j+kv)
                w(ko,ll+j*2) = var(ko,j+kv) - varavg(i,j+kv)
            end do
        end do
    end do
    CALL SUBMANDV()
end subroutine adjustData

SUBROUTINE PRINTDESC()
        use lsboth
    use procedures
        implicit none
        INTEGER:: IUN,I

    real(kind=8)::tempr(nobs),meany,miny,maxy,stdy,meanx(p),stdx(p),minx(p),maxx(p),meanu(r),stdu(r),minu(r),maxu(r),&
                    meanw(s),stdw(s),minw(s),maxw(s)

     meany=SUM(y(1:nobs))/DBLE(nobs)
     miny=minval(y(1:nobs))
     maxy=maxval(y(1:nobs))
     tempR(:)=0.0D0
     tempR(:)=(y(1:nobs)-meany)**2
     stdy=dsqrt(SUM(tempR)/DBLE(nobs-1))
     if (p>0) then
        call descripc(x,nobs,p,meanx,stdx,minx,maxx)
        ! standardize x if ncent=1
        if (ncent == 1) CALL STANDZ(X,NOBS,P,MEANX,STDX)
     end if
     if (r>0) then
        call descripc(u,nobs,r,meanu,stdu,minu,maxu)
        ! standardize u if ncent=1
        if (ncent == 1) CALL STANDZ(U,NOBS,R,MEANU,STDU)
     end if
     if (s>0) then
        call descripc(w,nobs,s,meanw,stdw,minw,maxw)
        ! standardize w if ncent=1
        if (ncent == 1) CALL STANDZ(W,NOBS,S,MEANW,STDW)
     end if
             IUN    = 16
             OPEN(UNIT=IUN,FILE="mixREGLS_both1.OUT")

             WRITE(IUN,'("MIXREGLS_both: Mixed-effects Location Scale Model")')
             write (IUN,*)
             WRITE(IUN,'("-----------------------------")')
             WRITE(IUN,'("mixREGLS_both.DEF specifications")')
             WRITE(IUN,'("-----------------------------")')
             WRITE(IUN,"(1x,18a4)")HEAD
             WRITE(IUN,*)
             WRITE(IUN,'(" data and output files:")')
             WRITE(IUN,"(1x,a80)")FILEDAT
             WRITE(IUN,"(1x,a80)")FILEOUT
             WRITE(IUN,*)
             WRITE(IUN,"(' MULITPLE LOCATION EFFECTS  = ',L1)")mls .eq. 1
             WRITE(IUN,"(' CONVERGENCE CRITERION = ',F11.8)")CONV
             WRITE(IUN,"(' RIDGEIN    = ',F8.4)")RIDGEIN
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
        if(discard0 .ne. 0) WRITE(IUN,508)num0
        508 FORMAT(//,1x,'==> The number of level 2 observations removed because of non-varying responses =', I6)
        if(discard0 .ne. 0) write(IUN,*) '(see '//trim(fileprefix)//'_removed.dat for information about those clusters)'
             write(IUN,'(" Number of level-1 observations for each level-2 cluster")')
             write(IUN,'(1x,13I6)') (IDNI(i,2), i=1,nc2)

        200  FORMAT(1x,A16,4F12.4)

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
             if (p>0) then
                WRITE(IUN,'(" Mean model covariates")')
                WRITE(IUN,'("                 mean         min         max     std dev")')
                WRITE(IUN,'(" --------------------------------------------------------")')
                do i=1,p
                   WRITE(IUN,200) BLAB(i),meanx(i),minx(i),maxx(i),stdx(i)
                end do
                WRITE(IUN,*)
             end if

             if (r>0) then
                WRITE(IUN,'(" BS variance model covariates")')
                WRITE(IUN,'("                 mean         min         max     std dev")')
                WRITE(IUN,'(" --------------------------------------------------------")')
                do i=1,r
                   WRITE(IUN,200) ALAB(i),meanu(i),minu(i),maxu(i),stdu(i)
                end do
                WRITE(IUN,*)
             end if

             if (s>0) then
                WRITE(IUN,'(" WS variance model covariates")')
                WRITE(IUN,'("                 mean         min         max     std dev")')
                WRITE(IUN,'(" --------------------------------------------------------")')
                do i=1,s
                   WRITE(IUN,200) TLAB(i),meanw(i),minw(i),maxw(i),stdw(i)
                end do
                WRITE(IUN,*)
             end if

        CLOSE(IUN)
END SUBROUTINE PRINTDESC


SUBROUTINE mixregmlsEST()
        use lsboth
        use procedures
        implicit none

        INTEGER :: I,J,LL,l2,k,kk,PP,SS,ii,jj, &
                IUN,CYCLES,NCYCLE,RP2,RPS2,NQ1,NS2,IFIN,ITER,NN, &
                   Q,NOB,IER,RIDGEIT,IUNS,myqdim,totalqR0,totalqR1,mytotalq,&
                   counter,cholamt,jj1,myorder(ndim2),myqdim2,rr2,myorder2(ndim2),cholcol
        REAL(KIND=8) :: RIDGE,LOGLP,PSUM,DET,LOGL,PI,XB,WT,WSVAR,ERRIJ,LPROB,PRA,LOGDIFF,MAXCORR,ONE,SCAL,&
                        ORIDGE,PVAL,ZVAL,RTEMP,RTEMP2,SMALL,BIG,LOGBIG,MAXDER,sdev,phiRatio,uchth,&
                        tauhat, tauhatlow, tauhatup, myz, temp
        REAL(KIND=8),ALLOCATABLE:: DER2A(:),DER2B(:),DER(:),LIK(:,:),H(:), &
                                   DERP(:),DERP2(:),DERQ(:),DERQ2(:),DERQQ(:),DZ(:),DERPS(:),COREC(:), &
                                   BDER2(:),ADER2(:),TDER2(:),SDER2(:), &
                                   ABDER2(:),TBDER2(:),TADER2(:),SBDER2(:),SADER2(:),STDER2(:), &
                                   TBADER2(:),SBADER2(:),SBATDER2(:),DZ2(:),SE(:), &
                                   THETA(:),thetav(:),WORK2(:),WORK3(:), myder2sq(:,:),&
                                   mytheta(:), myweights(:),mypoints(:,:), myweights0(:),&
                                   mypoints0(:,:),pointsR0(:,:),weightsR0(:),pointsR1(:,:),weightsR1(:),&
                                   cholTheta(:),uth(:),uthfull(:),mycholspar(:),mycholspar2(:),&
                                   sstar(:,:), work(:,:), sigma(:),asstar2(:,:),adjVar(:,:)
        character(len=80)::templabel



            ! parameters
        PI = 3.141592653589793238462643d0
        ONE = 1.0D0
        RTEMP = 1.0D0
        SMALL  = TINY(RTEMP)
        BIG    = HUGE(RTEMP)
        LOGBIG   = DLOG(BIG)
        !Requiring correlation between random scale and random location effect(s)
        ncov = 1

        PP = ((P+1)*P)/2
        rr2 = (rr+1)*rr/2
        rp2 = (rr+p)*(rr+p+1)/2
        IUN    = 16
        IUNS    = 17
        OPEN(UNIT=IUNS,FILE="mixREGLS_both_details_.ITS")
        open(unit=18, file="MixRegls_both_.its")

    ! number of quadrature points, quadrature nodes & weights

        ! ALLOCATE VECTORS FOR NR (these don't vary by CYCLES)
        RP2 = ((RR+P+1)*(RR+P))/2
        allocate(cholTheta(R))

        !Establishing arrays at maximum size needed
        SS = (S*(S+1))/2
        RPS2 = (RR+P+S)*(RR+P+S+1)/2
        ! NS = number of additional var cov parameters due to random SCALE
        NS = ndim
        NS2 = NS*(NS+1)/2
        NPAR = P+RR+S+NS ! number of parameters
        NPAR2 = NPAR*(NPAR+1)/2
        totalqR0 = nq ** R
        totalqR1 = nq ** (R+1)
        mytotalq = totalqR1
        myqdim = ndim
        allocate(mytheta(ns))
        ALLOCATE (THETA(ndim))
        ALLOCATE (thetav(ndim2))
        ALLOCATE (SE(NPAR))
        ALLOCATE (COREC(NPAR))
        corec(:) = 0D0
        allocate(uthfull(r*r))
        allocate(uth(rr))
        ALLOCATE (WORK2(ndim))
        ALLOCATE (WORK3(ndim2))

        ALLOCATE(BDER2(PP))
        ALLOCATE(ABDER2(P*RR))
        ALLOCATE(ADER2(RR*(RR+1)/2))
        ALLOCATE (DER2A(NPAR2))
        ALLOCATE (DER2B(NPAR2))
        allocate(myder2sq(npar,npar))
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
        ALLOCATE (TBDER2(S*P))
        ALLOCATE (TADER2(S*RR))
        ALLOCATE (TBADER2(S*(P+RR)))
        ALLOCATE (TDER2(SS))
        ALLOCATE (DZ2(NPAR2))
        ALLOCATE (SBDER2(NS*P))
        ALLOCATE (SADER2(NS*RR))
        ALLOCATE (STDER2(NS*S))
        ALLOCATE (SBADER2(NS*(P+RR)))
        ALLOCATE (SBATDER2(NS*(P+RR+S)))
        ALLOCATE (SDER2(NS2))

        cholamt = rr+ns
        cholcol = r+1
        allocate(sstar(cholamt,cholamt))
        allocate(work(npar,npar))
        allocate(adjVar(npar,npar))
        allocate(asstar2(npar,npar))
        allocate(sigma(ndim2))
        allocate(mycholspar(cholamt))
        allocate(mycholspar2(cholamt))

        allocate(myweights(nqwR1))
        allocate(mypoints(nqwR1,r+1))
        allocate(myweights0(nqwR1))
        allocate(mypoints0(nqwR1,r+1))
        allocate(weightsR1(NQwr1))
        allocate(pointsR1(NQwr1,r+1))
        allocate(weightsR0(NQwR0))
        allocate(pointsR0(NQwR0,r))
        call getquad(r, nq, nqwr0, pointsR0, weightsR0)
        call getquad(r+1, nq, nqwr1, pointsR1, weightsR1)


        ! start cycles
        ! cycles = 1: random intercept model with BS variance terms
        ! cycles = 2: add in scale (WS) variance terms
        ! cycles = 3: add in random scale
        ! cycles = 4: use NS = R+1

        OPEN(UNIT=IUN,FILE="MIXREGLS_both2.OUT")
        ncycle = 4
        CYCLELOOP:do cycles=2,ncycle
            if(cycles > 2 .and. nors==1) cycle
        if (CYCLES==2) THEN
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
                SS = (S*(S+1))/2
                RPS2 = (RR+P+S)*(RR+P+S+1)/2
                NS = 0 ! number of additional var cov parameters due to random SCALE
                NPAR = P+RR+S ! number of parameters
                NPAR2 = NPAR*(NPAR+1)/2
                myqdim = R
                myqdim2 = myqdim*(myqdim+1)/2
                mytotalq = totalqR0
                mypoints0(1:mytotalq,1:myqdim) = pointsR0(1:mytotalq,1:myqdim)
                myweights0(1:mytotalq) = weightsR0(1:mytotalq)
                tau(2:S) = 0
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
                 SS = (S*(S+1))/2
                RPS2 = (RR+P+S)*(RR+P+S+1)/2
               ! NS = number of additional var cov parameters due to random SCALE
                NS = 1
                NS2 = NS*(NS+1)/2
                NPAR = P+RR+S+NS ! number of parameters
                NPAR2 = NPAR*(NPAR+1)/2
                myqdim = R+1
                myqdim2 = myqdim*(myqdim+1)/2
                mytotalq = totalqR1
                mypoints0 = pointsR1
                myweights0 = weightsR1
            else if (CYCLES==4) THEN
                SS = (S*(S+1))/2
                RPS2 = (RR+P+S)*(RR+P+S+1)/2
                ! NS = number of additional var cov parameters due to random SCALE
                NS = 1+R
                NS2 = NS*(NS+1)/2
                NPAR = P+RR+S+NS ! number of parameters
                NPAR2 = NPAR*(NPAR+1)/2
                myqdim = R+1
                myqdim2 = myqdim*(myqdim+1)/2
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

             ! set the ridge up for the models with random scale
             IF (CYCLES>=3) THEN
                ridge = RIDGEIN+.1D0
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
                IF (CYCLES>=3 .AND. ITER==13) THEN
                    ridge = ridge - .1
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
                    if(aquad .ne. 0 .and. (iter >= 8 .or. cycles .eq. 2 .or. cycles .eq. 4)) then
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
                        DERQ(1:npar)  = 0.0D0
                        DERQ2(1:npar2) = 0.0D0
                        DZ(1:npar)    = 0.0D0

                        JLOOP:DO J=1,IDNI(I,2)  ! loop over level-1 observations
                            NOB=NN+J
                            XB = DOT_PRODUCT(BETA,X(NOB,1:p))        ! X BETA for the current LEVEL-1 obs
                            ii=0
                            uchth = 0
                            do k=1, R
                                do jj=1, k
                                    ii = ii + 1
                                    uchth = uchth + U(nob,k) * mychol(ii) * mypoints(q, jj)
                                end do
                            end do

                            WT = DOT_PRODUCT(TAU(1:S),W(NOB,1:S))  ! W TAU for the current LEVEL-1 obs
                                                                   ! note that S changes over CYCLES

                            !mytheta is the quadrature points vector for spar
                            if(cycles >= 3) then
                                mytheta(1:ns) = mypoints(q,(R+1+1-ns):(R+1))
                            end if

                            if(ns > 0) then
                                RTEMP = WT + DOT_PRODUCT(mytheta(1:ns), spar(1:ns))
                                IF (RTEMP .GE. LOGBIG) THEN
                                  WSVAR = BIG
                                ELSE
                                  WSVAR = DEXP(RTEMP)
                                END IF
                            ELSE
                               IF (WT .GE. LOGBIG) THEN
                                  WSVAR = BIG
                               ELSE
                                  WSVAR = DEXP(WT)
                               END IF
                            END IF

                            !Uth is equivalent to JR* in the equations
!                            call mpytr(u(nob,:), mypoints(q,1:R), uthfull,1,R,0,R)
            !               call chams0to1(uthfull, uth, R)
!                            call chams(uthfull,uth,r,0,3)
                            ii=0
                            uth = 0
                            do k=1, R
                                do jj=1, k
                                    ii = ii + 1
                                    uth(ii) = + U(nob,k) * mypoints(q, jj)
                                end do
                            end do


                            ERRIJ = Y(NOB) - (XB + uchth) !uchth replaces ua
                            IF (WSVAR .LE. SMALL) WSVAR = SMALL
                            LPROB = -.5D0*(DLOG(2.0D0*PI)+DLOG(WSVAR) + ERRIJ**2/WSVAR)
                            PSUM  = PSUM + LPROB
                             DZ(1:P) = (-2.0D0/WSVAR)*ERRIJ*X(NOB,:)                    ! beta
                            dz(p+1:p+rr) = (-2.0D0/WSVAR)*ERRIJ*uth                         !cholesky
                            DZ(P+Rr+1:P+Rr+S) = (1.0D0 - (ERRIJ**2/WSVAR))*W(NOB,1:S)    ! tau
                            if(ns > 0) dz(p+rr+s+1:p+rr+s+ns) = (1.0D0 - (ERRIJ**2/WSVAR))*mytheta(1:ns) ! sparam

            !write(IUNS,'(17F8.3)') i*1.0,q*1.0,j*1.0,mypoints(q,1),lprob,psum,errij,y(nob),xb,uchth,wt,(dz(k), k=1,npar)
                            DERQ(1:npar) = DERQ(1:npar) + (-.5D0)*DZ(1:npar)
!write(iuns,*) errij,dz(1:npar),derq(1:npar),w(nob,1:s)

                               ! beta         (p x p symmetric)
                            CALL GRAM(X(NOB,:),BDER2,P,1)
                            BDER2 = (-1.0D0/wsvar)*BDER2

                            ! mychol, beta  (rr rows & p columns)
                            CALL MPYM(Uth,X(NOB,:),ABDER2,RR,1,0,0,P)
                            abder2 = -abder2/wsvar

                            ! Mychol is a stacked triangular R x R matrix
                            CALL GRAM(Uth,ADER2,RR,1)
                            ADER2 = -ADER2/wsvar

                            CALL ADJRC(BDER2,ABDER2,ADER2,DZ2(1:RP2),RR,P)

                            ! tau, beta    (s rows & p columns)
                            ii=0
                            do k=1, S
                                do jj=1, P
                                    ii = ii + 1
                                    tbder2(ii) = w(nob,k)*x(nob,jj)
                                end do
                            end do
!                            CALL MPYM(W(NOB,1:S),X(NOB,:),TBDER2,S,1,0,0,P)
                            TBDER2 = (-1.0D0/wsvar)*errij*TBDER2

                            ! tau, mychol   (s rows & Rr columns)
                            CALL MPYM(W(NOB,1:S),uth,TADER2,S,1,0,0,RR)

                           TADER2 = (-errij/wsvar)*TADER2
                            ! tau
                            CALL GRAM(W(NOB,1:S),TDER2,S,1)
                            TDER2 = (-.5D0/wsvar)*errij*errij*TDER2

                            CALL ADJC(TBDER2,TADER2,TBADER2,S,P,RR)
                            CALL ADJRC(DZ2(1:RP2),TBADER2,TDER2,DZ2(1:RPS2),S,P+RR)

                               ! random scale
                            IF (NS > 0) THEN
                                ! scale, beta    (NS rows & p columns)
                                CALL MPYM(mytheta,X(NOB,:),SBDER2,NS,1,0,0,P)
                                SBDER2 = (-1.0D0/wsvar)*errij*SBDER2
                                ! scale, alpha   (NS rows & Rr columns)
                                CALL MPYM(mytheta,Uth,SADER2,NS,1,0,0,RR)
                                SADER2 = (-errij/wsvar)*SADER2
                                ! scale, tau     (NS rows & s columns)
                                CALL MPYM(mytheta,W(NOB,1:S),STDER2,NS,1,0,0,S)
                                STDER2 = (-.5D0/wsvar)*errij*errij*STDER2
                                ! scale          (2 x 2 symmetric)

                                CALL GRAM(mytheta,SDER2,NS,1)
                                SDER2 = (-.5D0/wsvar)*errij*errij*SDER2

                                CALL ADJC(SBDER2,SADER2,SBADER2,NS,P,RR)
                                CALL ADJC(SBADER2,STDER2,SBATDER2,NS,P+RR,S)
                                CALL ADJRC(DZ2(1:RPS2),SBATDER2,SDER2,DZ2(1:NPAR2),NS,P+RR+S)
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
                    theta(1:myqdim) = 0
                    thetav(1:myqdim2) = 0
                    do k=1, myqdim
                        do q=1,mytotalq
                               theta(k) = theta(k) + lik(i,q)*myweights(q)*mypoints(q,k)
                        end do
                    end do
                    theta(1:myqdim) = theta(1:myqdim)/h(i)
                    do q=1,mytotalq
                        temp = lik(i,q)*myweights(q)
                        do k=1, myqdim
                           work2(k) = mypoints(q,k) - theta(k)
                        end do
                        call grmcv(thetav(1:myqdim2), thetav(1:myqdim), work2, temp, myqdim)
                    end do
                    thetav(1:myqdim2) = thetav(1:myqdim2)/h(i)

                    thetas(i,1:myqdim) = theta(1:myqdim)
                    thetavs(i,1:myqdim2) = thetav(1:myqdim2)

                    LOGL  = LOGL + DLOG(H(I))
                    SCAL  = DEXP(0.0D0 - DLOG(H(I)))
                    DERPS(1:npar) = SCAL*DERP(1:npar)
                    DER(1:npar)   = DER(1:npar) + DERPS(1:npar)
!write(iuns,*) derps(1:npar),der(1:npar)
!                    CALL GRMCV(DER2A(1:npar2),DER2A(1:npar2),DERPS(1:npar),ONE,NPAR)
                            ii=0
                            do k=1, npar
                                do jj=1, k
                                    ii = ii + 1
                                    der2a(ii) = der2a(ii) + derps(k)*derps(jj)
                                end do
                            end do
                    DER2B(1:npar2) = DER2B(1:npar2) + SCAL*DERP2(1:npar2)
                END DO ILOOP

                 LOGDIFF = LOGL-LOGLP
                 LOGLP   = LOGL

                 ! determine if an NR iteration is bad and increase the ridge
                 ! take the ridge off after 10 good iterations
                 IF (LOGDIFF/LOGLP > .05 .AND. ITER < MAXIT) THEN
                    RIDGEIT = 0
                    RIDGE = RIDGE + .05D0
                    WRITE(IUN,'("==> BAD NR ITERATION ",I5," with NEW ridge = ",F8.4,/)') ITER,RIDGE
                write(IUNS,*)"Corrections"
                write(IUNS,'(20f11.3)') (corec(k), k=1,npar)
                write(IUNS,*)"Derivatives"
                write(IUNS,'(20f11.3)') (der(k), k=1,npar)
                write(IUNS,*)"Beta"
                write(IUNS,'(20f11.3)') (beta(k), k=1,p)
                write(IUNS,*)"Chol"
                write(IUNS,'(20f11.3)') (mychol(k),k=1,rr)
                write(IUNS,*)"Tau"
                write(IUNS,'(20f11.3)') (tau(k),k=1,s)
                if(ns > 0) then
                    write(IUNS,*)"Spar"
                    write(IUNS,'(20f11.3)') (spar(k),k=1,ns)
                end if
                    GO TO 99
                 END IF
                 IF (LOGDIFF/LOGLP <= .000001 .AND. RIDGEIT < 10) THEN
                    RIDGEIT = RIDGEIT+1
                 ELSE IF (LOGDIFF/LOGLP <= .000001 .AND. RIDGEIT >= 10 .and. ifin==1) then
                        ridge = ridge - .05D0
                    if(maxder < 1 .or. (iter > 50 .and. maxder < 5) .or. (iter > 75 .and. maxder < 10)) then
                        ridge = 0
                    else
                        if ( ridge < RIDGEIN) ridge = RIDGEIN
                    end if
                 END IF

                if(ifin < 2) then
                     ! ridge adjustment - diagonal elements only
                     DO k=1,NPAR
                        LL=(k*(k+1))/2
                        DER2A(LL)=DER2A(LL)+RIDGE*DER2A(LL)
                     END DO
                     DER2B = DER2A - DER2B
                     DO k=1,NPAR
                        LL=(k*(k+1))/2
                        DER2B(LL)=DER2B(LL)+RIDGE*DER2B(LL)
                     END DO
                else
                    DER2B(1:npar2) = DER2A(1:npar2) - DER2B(1:npar2)
                end if

                write(IUNS,*)"2nd Derivatives"
                ll = 0
                do k=1,npar
                    write(IUNS,'(25f10.2)') (der2b(ll+kk), kk=1,k)
                    ll = ll + k
                end do

                 ! INVERT 2nd DERIVATIVE MATRIX - COREC is a working vector for INVS
                 IF (IFIN == 1) THEN
                    CALL INVS(DER2B,NPAR,DET,IER)
                    IF (IER == 0) THEN
                        WRITE(*,'(" Newton-Raphson Iteration ",I5," with ridge ",F8.4)') ITER,RIDGE
                        WRITE(IUNS,'(" Newton-Raphson Iteration ",I5," with ridge ",F8.4)') ITER,RIDGE
                        WRITE(18,'(" Newton-Raphson Iteration ",I5," with ridge ",F8.4)') ITER,RIDGE
                       CALL MPYM(DER2B,DER,COREC,NPAR,NPAR,1,0,1)
                    ELSE IF (IER == 1) THEN
                       WRITE(*,'(" BHHH Iteration ",I5," with ridge ",F8.4)') ITER,RIDGE
                       WRITE(IUNS,'(" BHHH Iteration ",I5," with ridge ",F8.4)') ITER,RIDGE
                       WRITE(18,'(" BHHH Iteration ",I5," with ridge ",F8.4)') ITER,RIDGE
                       CALL INVS(DER2A,NPAR,DET,IER)
                       CALL MPYM(DER2A,DER,COREC,NPAR,NPAR,1,0,1)
                    END IF
                 ELSE
                    CALL INVS(DER2B,NPAR,DET,IER)
                    CALL MPYM(DER2B,DER,COREC,NPAR,NPAR,1,0,1)
                    DO k=1,NPAR
                       LL=(k*(k+1))/2
                       SE(k)=DSQRT(DABS(DER2B(LL)))
                    END DO
                 END IF
             do k=1,npar
                 if(corec(k) > 1) corec(k) = 1
                 if(corec(k) < -1) corec(k) = -1
            end do
            if(iter<=5) corec(1:npar) = corec(1:npar)*.5
            if(ifin==1) then

                write(IUNS,*)"Corrections"
                write(IUNS,'(20f11.3)') (corec(k), k=1,npar)
                write(IUNS,*)"Derivatives"
                write(IUNS,'(20f11.3)') (der(k), k=1,npar)
                write(IUNS,*)"Beta"
                write(IUNS,'(20f11.3)') (beta(k), k=1,p)
                write(IUNS,*)"Chol"
                write(IUNS,'(20f11.3)') (mychol(k),k=1,rr)
                write(IUNS,*)"Tau"
                write(IUNS,'(20f11.3)') (tau(k),k=1,s)
                if(ns > 0) then
                    write(IUNS,*)"Spar"
                    write(IUNS,'(20f11.3)') (spar(k),k=1,ns)
                end if

                 MAXDER=MAXVAL(ABS(DER(1:npar)))
                 MAXCORR=MAXVAL(ABS(COREC(1:npar)))
                 WRITE(*,'(i4,a)') iter,' maximum correction and derivative and ridge'
                 WRITE(*,'(2g16.4,f12.2)') maxcorr,MAXDER,ridge
                  WRITE(IUNS,'(i4,a)') iter,' maximum correction and derivative and ridge'
                 WRITE(IUNS,'(2g16.4,f12.2)') MAXCORR,MAXDER,ridge
                  WRITE(18,'(i4,a)') iter,' maximum correction and derivative and ridge'
                 WRITE(18,'(2g16.4,f12.2)') MAXCORR,MAXDER,ridge

                 ! done with NR and onto last iteration
                 IF (IFIN==1 .AND. (MAXCORR <= CONV .OR. ITER >= MAXIT)) THEN
                     IFIN=2
                     ORIDGE=RIDGE
                     RIDGE=0.0D0
                 END IF

                 ! UPDATE PARAMETERS
                 BETA  = BETA  + COREC(1:P)
                 if(iter >= 10 .or. (iter>=5 .and. cycles .eq. 2)) then
                    mychol = mychol + COREC(P+1:P+RR)
                else
                    corec(p+1:p+rr) = 0
                end if
                    if(ns > 1) then
                        spar(1:ns-1) = spar(1:ns-1) + corec(npar-ns+1:npar-1)*.5
                    end if
                    if(ns > 0 .and. iter > 5) then
                        spar(ns) = spar(ns) + corec(npar)*.5
                    else if(ns > 0) then
                        corec(npar) = 0
                    end if
                     TAU   = TAU      + COREC(P+RR+1:P+RR+S)
                    kk = 1
                    do k=1, r
                        do LL=1, k
                            if(k==LL .and. mychol(kk) < 0) mychol(kk) = -mychol(kk)
                            kk = kk + 1
                        end do
                    end do
            end if

            99  WRITE(*,'("   -2 Log-Likelihood = ",F14.5)') -2*LOGL
                 WRITE(IUNS,'("   -2 Log-Likelihood = ",F14.5)') -2*LOGL
                 WRITE(18,'("   -2 Log-Likelihood = ",F14.5)') -2*LOGL
                 ITER = ITER+1
            END DO IFINLOOP

        if(cycles .ne. 3) then
write(iuns,*) "Original var/covar values"
                ll = 0
                do k=1,npar
                    write(IUNS,'(25f18.6)') (der2b(ll+kk), kk=1,k)
                    ll = ll + k
                end do
            cholamt = rr
            cholcol = r
            mycholspar(1:rr) = mychol
            sigma(1:rr) = mychol
            do j=1,ndim2
                myorder(j) = j
            end do
            myorder2 = myorder
            ii=1
            do i=1,npar
                do j=1,i
                    myder2sq(i,j) = der2b(ii)
                    myder2sq(j,i) = myder2sq(i,j)
                    ii=ii+1
                end do
            end do
            if(cycles >= 3) then
                mycholspar(rr+1:rr+ns) = spar(1:ns)
                sigma(rr+1:rr+ns) = spar(1:ns)
                if(chol .ne. 1 .and. chol .ne. 2) then
                    cholamt = rr+ns
                    cholcol = r+1
                    do j=0,r-1
                        jj1 = j*(j-1)/2
                        do i=1,r-j
                            myorder(j*r-jj1+i) = r*j-jj1+j+i
                        end do
                        myorder(rr+j+1)= r*(j+1)-jj1+1
                    end do
                end if
            end if
                    do j=1,ndim2
                        myorder2(j) = myorder(j)+p
                         if(myorder(j)>rr) myorder2(j) = myorder2(j)+s
                    end do
write(iuns,*) "Cholesky"
write(iuns,'(20f6.2)') (mycholspar(j),j=1,cholamt)
                    do i=1,cholamt
                        mycholspar2(i) = mycholspar(myorder(i))
                    end do
write(iuns,'(20i3)') (myorder(j),j=1,cholamt)
write(iuns,'(20f6.2)') (mycholspar2(j),j=1,cholamt)
            if(chol .ne. 2) then
!    open(7,file="sstarinfo.out",access="append")
                call getSStar(mycholspar2(1:cholamt),cholcol,cholamt,sstar(1:cholamt,1:cholamt))
write(iuns,*) "Matrix to transform cholesky back"
do i=1,cholamt
    write(iuns,'(20f6.2)') (sstar(i,j),j=1,cholamt)
                end do
                do i=1,cholamt
                    sigma(i) = dot_product(mycholspar(1:cholamt),sstar(i,1:cholamt))
                end do
!                sigma = matmul(sstar,mycholspar)
write(iuns,*) "Sigma"
write(iuns,'(20f6.2)') (sigma(j),j=1,cholamt)
                asstar2 = 0
                ii=1
                do i=1,npar
                    asstar2(i,i) = 1
                    jj=1
                    if((i > p .and. i <= p+rr) .or. (i>=npar+rr-cholamt+1)) then
                        do j=1,npar
                            if((j > p .and. j <= p+rr) .or. (j>=npar+rr-cholamt+1)) then
                                asstar2(i,j) = 2*sstar(ii,jj)
                                jj=jj+1
                            end if
                        end do
                        ii=ii+1
                    end if
                end do
write(iuns,*) "Matrix to transform var/covar"
do i=1,npar
    write(iuns,'(20f6.2)') (asstar2(i,j),j=1,npar)
end do
write(iuns,*) "Original var/covar values"
do i=1,npar
    write(iuns,'(20F15.8)') (myder2sq(i,j), j=1,npar)
end do
                work(1:npar,1:npar) = matmul(asstar2(1:npar,1:npar), myder2sq(1:npar,1:npar))
                adjVar(1:npar,1:npar) = matmul(work(1:npar,1:npar),&
                                                transpose(asstar2(1:npar,1:npar)))
            else
                sigma(1:rr) = mycholspar(1:rr)
                adjvar = myder2sq
            end if

            DO k=1,npar
                se(k) = dsqrt(dabs(adjvar(k,k)))
            END DO

                 ! WRITE RESULTS
            WRITE(IUN,562)ITER-1,ORIDGE,LOGL,LOGL-NPAR, &
            LOGL-0.5D0*DBLE(NPAR)*DLOG(DBLE(NC2)),0-2*LOGL,0-2*(LOGL-NPAR), &
            0-2*(LOGL-0.5D0*DBLE(NPAR)*DLOG(DBLE(NC2)))
            562 FORMAT(1X,'Total  Iterations =',I4,/, &
                          1X,'Final Ridge value =',F5.2,//, &
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
             DO k=1,P
                ZVAL = BETA(k)/SE(k)
                PVAL = 2.0D0 *(1.0D0 - PHIFN(DABS(ZVAL),0))
                WRITE(IUN,804)BLAB(k),BETA(k),SE(k),ZVAL,PVAL
             END DO
            write(iun,*)
            if(chol .ne. 2) then
                WRITE(IUN,'("Random (location) Effect Variances and Covariances")')
            else
                WRITE(IUN,'("Random (location) Effect Variances Cholesky terms")')
            end if
            counter = 1
             do i=1,r
                 do j=1,i
                    ZVAL = sigma(myorder(counter))/SE(myorder2(counter))
                    PVAL = 2.0D0 *(1.0D0 - PHIFN(DABS(ZVAL),0))
                    templabel = "Covariance"
                    if(i==j) templabel = alab(i)
                    WRITE(IUN,804) templabel,sigma(myorder(counter)),SE(myorder2(counter)),ZVAL,PVAL
                    counter = counter+1
                end do
             END DO

            804 FORMAT(A16,4(4x,F12.5))
                write(iun,*)
            IF (NCENT==1 .AND. S>1) THEN
                     WRITE(IUN,'("STANDARDIZED TAU (WS variance parameters: log-linear model)")')
                 ELSE
                     WRITE(IUN,'("TAU (WS variance parameters: log-linear model)")')
                 END IF
                 DO k=1,S
                    L2 = P+RR+k
                    ZVAL = TAU(k)/SE(L2)
                    PVAL = 2.0D0 *(1.0D0 - PHIFN(DABS(ZVAL),0))
                    WRITE(IUN,804)TLAB(k),TAU(k),SE(L2),ZVAL,PVAL
                 END DO
                 if(cycles == 4) then
                     write(iun,*)
write(iuns,*) (myorder(k),k=1,ndim2)
write(iuns,*) (myorder2(k),k=1,ndim2)
            if(chol .ne. 2 .and. chol .ne. 1) then
                WRITE(IUN,'("Random Scale Variance and Covariance")')
            else
                WRITE(IUN,'("Random Scale Variance Cholesky terms")')
            end if
                    do k=1,ns
                        ZVAL = sigma(myorder(rr+k))/se(myorder2(rr+k))
                        PVAL = 2.0D0 *(1.0D0 - PHIFN(DABS(ZVAL),0))
                        templabel = "Scale int var"
                        if(k<ns) templabel = "Cov "//ALAB(k)
                        WRITE(IUN,804) templabel,sigma(myorder(rr+k)),se(myorder2(rr+k)),ZVAL,PVAL
                    end do
            end if

         IF (MAXDER > CONV .AND. ITER >= MAXIT) THEN
            WRITE(IUN,'("NOTE: CONVERGENCE CRITERION WAS NOT ACHIEVED")')
            WRITE(IUN,'("FINAL FIRST DERIVATIVE AND CORRECTION VALUES")')
            DO k=1,NPAR
              WRITE(IUN,*)DER(k),COREC(k)
            END DO
         END IF

         ! write out the deviance, estimates, and standard errors
            OPEN(2, FILE='mixregls_both_.EST',access="append")
            WRITE(2,'(F15.6,2I8)') -2*LOGL,ITER-1,MAXIT
             WRITE(2,'(35F15.8)')(BETA(k),k=1,P)
             WRITE(2,'(35F15.8)')(sigma(myorder(k)),k=1,RR)
             WRITE(2,'(35F15.8)')(TAU(k),k=1,S)
             write(2,'(35F15.8)')(sigma(myorder(rr+k)),k=1,ns)
             WRITE(2,'(35F15.8)')(SE(k),k=1,p),(se(myorder2(k)),k=1,rr), &
             (se(k),k=p+rr+1,p+rr+s),(se(myorder2(k)),k=rr+1,rr+ns)
            close(2)
        end if
        close(IUN)
        OPEN(UNIT=IUN,FILE="MIXREGLS_both2.OUT",access="append")
     END DO CYCLELOOP

     write(iun,*)
     write(iun,*)
     write(iun,*)
808 FORMAT(A16,3(4x,A12))
     WRITE(IUN,'("WS variance ratios and 95% CIs")')
     write(iun,'("------------------------------")')
     write(iun,*)
    WRITE(IUN,808) 'Variable        ','Ratio','Lower','Upper'
    write(iun,808)'---------------------','------------------','------------','------------'
         WRITE(IUN,'("TAU (WS variance parameters: log-linear model)")')
         myz = 1.959964
     DO k=1,S
        L2 = P+RR+k
        tauhat = exp(tau(k))
        tauhatlow = exp(tau(k)-myz*se(l2))
        tauhatup = exp(tau(k)+myz*se(l2))
        WRITE(IUN,804)TLAB(k),tauhat, tauhatlow, tauhatup
     END DO
         CLOSE(IUN)
         CLOSE(IUNS)
         CLOSE(1)
         CLOSE(2)
         CLOSE(3)
         close(18)

    open(23,file=trim(fileprefix)//"_ebvar.dat")
    k=myqdim*(myqdim+1)/2
    DO I=1,NC2  ! go over level-2 clusters
        write(23,'(i8,20F15.8)') idni(i,1),(THETAs(I,j), j=1,myqdim), (thetavs(i,j), j=1,k)
    end do
    close(23)
    open(37,file="mixregls_both_.var")
    do i=1,npar
        write(37,'(20F15.8)') (adjvar(i,j), j=1,npar)
    end do
    close(37)
        deallocate(weightsR1,pointsR1,weightsR0,pointsR0,myweights,mypoints)
        deallocate(bder2,abder2,ader2,der2a,der2b,der,derp,derp2,derq,derq2,derqq,dz,derps,tader2, &
                    tbder2,tbader2,tder2,dz2,sbder2,sader2,stder2,sbader2,sbatder2,sder2)
        deallocate(cholTheta,mytheta,lik,h,corec,se,theta,thetav,work2,work3,uthfull,uth)
        deallocate(sstar,work,adjvar,asstar2,sigma,mycholspar)
END SUBROUTINE MIXREGmLSEST

SUBROUTINE mixREGLSEST()
    use lsboth
    use procedures
    implicit none
    INTEGER :: I,J,L,LL,l2,k,PP,SS, & !RR
               IUN,CYCLES,NCYCLE,RP2,RPS2,NQ1,S0,NS2,IFIN,ITER,NN, &
               Q,NOB,IER,RIDGEIT,IUNS,myqdim,totalqR0,totalqR1,mytotalq,ncov0
    REAL(KIND=8) :: RIDGE,LOGLP,PSUM,DET,LOGL,PI,XB,WT,WSVAR,ERRIJ,LPROB,PRA,LOGDIFF,MAXCORR,ONE,SCAL,&
                    ORIDGE,PVAL,ZVAL,RTEMP,RTEMP2,SMALL,BIG,LOGBIG,MAXDER,sdev,phiRatio,&
                    ua,bsvar,myz,tauhat,tauhatlow,tauhatup,temp
    REAL(KIND=8),ALLOCATABLE:: DER2A(:),DER2B(:),DER(:),LIK(:,:),H(:), &
                               DERP(:),DERP2(:),DERQ(:),DERQ2(:),DERQQ(:),DZ(:),DERPS(:),COREC(:), &
                               BDER2(:),ADER2(:),TDER2(:),SDER2(:), &
                               ABDER2(:),TBDER2(:),TADER2(:),SBDER2(:),SADER2(:),STDER2(:), &
                               TBADER2(:),SBADER2(:),SBATDER2(:),DZ2(:),SE(:), &
                               THETA(:),thetav(:),WORK2(:),WORK3(:), &
                               myweights(:),mypoints(:,:), myweights0(:),&
                               mypoints0(:,:),pointsR0(:,:),weightsR0(:),pointsR1(:,:),weightsR1(:),&
                               vaug(:),varSQ(:,:)

        ! parameters
    PI = 3.141592653589793238462643d0
    ONE = 1.0D0
    RTEMP = 1.0D0
    SMALL  = TINY(RTEMP)
    BIG    = HUGE(RTEMP)
    LOGBIG   = DLOG(BIG)

    PP = ((P+1)*P)/2
    RR = ((R+1)*R)/2
    ndim = 2
    ndim2 = (ndim+1)*ndim/2
    IUN    = 16
    OPEN(UNIT=IUN,FILE="mixREGLS_both2.OUT")
    IUNS    = 17
    OPEN(UNIT=IUNS,FILE="mixREGLS_both_details_.ITS")
    open(unit=18, file="MixRegls_both_.its")
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
        allocate(varSQ(npar,npar))

        ! start cycles
        ! cycles = 1: random intercept model with BS variance terms
        ! cycles = 2: add in scale (WS) variance terms
        ! cycles = 3: add in random scale
        ! cycles = 4: use NS = R+1

        ncycle = 4
        CYCLELOOP:do cycles=1,ncycle
            if(cycles==4 .and. ncov0==0) cycle

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
                        call grmcv_inout(thetav, work2, temp, myqdim)
                    end do
                    thetav = thetav/h(i)

                    thetas(i,:) = theta
                    thetavs(i,:) = thetav

                    LOGL  = LOGL + DLOG(H(I))
                    SCAL  = DEXP(0.0D0 - DLOG(H(I)))
                    DERPS = SCAL*DERP
                    DER   = DER + DERPS
                    CALL GRMCV_inout(DER2A,DERPS,ONE,NPAR)
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
     END IF
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
                    CALL INVS(DER2B,NPAR,DET,IER)
                    IF (IER == 0) THEN
                        WRITE(*,'(" Newton-Raphson Iteration ",I5," with ridge ",F8.4)') ITER,RIDGE
                        WRITE(IUNS,'(" Newton-Raphson Iteration ",I5," with ridge ",F8.4)') ITER,RIDGE
                        WRITE(18,'(" Newton-Raphson Iteration ",I5," with ridge ",F8.4)') ITER,RIDGE
                       CALL MPYM(DER2B,DER,COREC,NPAR,NPAR,1,0,1)
                    ELSE IF (IER == 1) THEN
                       WRITE(*,'(" BHHH Iteration ",I5," with ridge ",F8.4)') ITER,RIDGE
                       WRITE(IUNS,'(" BHHH Iteration ",I5," with ridge ",F8.4)') ITER,RIDGE
                       WRITE(18,'(" BHHH Iteration ",I5," with ridge ",F8.4)') ITER,RIDGE
                       CALL INVS(DER2A,NPAR,DET,IER)
                       CALL MPYM(DER2A,DER,COREC,NPAR,NPAR,1,0,1)
                    END IF
                 ELSE
                    CALL INVS(DER2B,NPAR,DET,IER)
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

                 ! done with NR and onto last iteration
                 IF (IFIN==1 .AND. (MAXCORR <= CONV .OR. ITER >= MAXIT)) THEN
                     IFIN=2
                     ORIDGE=RIDGE
                     RIDGE=0.0D0
                 END IF

                 ! UPDATE PARAMETERS
                 BETA  = BETA  + COREC(1:P)
                 if(iter >= 10 .or. cycles == 1) alpha = alpha + corec(p+1:p+r)
                 if(iter >= 5 .and. ns>0) then
                    do k=1,ns
                        spar(k) = spar(k) + corec(npar-ns+k)
                    end do
                end if
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
                 IF (CYCLES.EQ.1) OPEN(2, FILE='mixREGLS_both_.EST')
                    WRITE(2,'(F15.6,2I8)') -2*LOGL,ITER-1,MAXIT
                     WRITE(2,'(35F15.8)')(BETA(L),L=1,P)
                     WRITE(2,'(35F15.8)')(alpha(L),L=1,R)
                     WRITE(2,'(35F15.8)')(TAU(L),L=1,S)
                     write(2,'(35F15.8)')(SPAR(L),L=1,ns)
                     WRITE(2,'(35F15.8)')(SE(L),L=1,NPAR)
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
        write(23,'(i16,20F15.8)') idni(i,1), (THETAs(I,j), j=1,myqdim), (thetavs(i,j), j=1,k)
    end do
    close(23)
    ll=1
    do i=1,npar
        do j=1,i
            varSQ(i,j) = der2b(ll)
            varSQ(j,i) = varSQ(i,j)
            ll=ll+1
        end do
    end do
    open(37,file="mixregls_both_.var")
    do k=1,npar
        write(37,'(25f18.6)') (varSQ(k,j), j=1,npar)
    end do
    close(37)
    deallocate(weightsR1,pointsR1,weightsR0,pointsR0,myweights,mypoints)
    deallocate(bder2,abder2,ader2,der2a,der2b,der,derp,derp2,derq,derq2,derqq,dz,derps,tader2, &
                tbder2,tbader2,tder2,dz2,sbder2,sader2,stder2,sbader2,sbatder2,sder2)
    deallocate(vaug,lik,h,corec,se,theta,thetav,work2,work3,varsq)
END SUBROUTINE mixREGLSEST


subroutine callmixreg
    use lsboth
    use procedures
    implicit none
    integer::i,j,iun,counter,myio,k
    real(kind=8)::temp,logl,pval,zval,se(p+rr+1),temp3(2+rr+r+p),v(rr)
    CHARACTER(LEN=16)::templabel

    open(10,file="mixregls_both_temp_.dat")
    do i=1,nobs
        write(10,*) ids(i), y(i), (x(i,j),j=1,p), (u(i,j),j=1,numloc)
    end do
    close(10)

    OPEN(1,FILE="mixreg.def")
    WRITE(1,'(18A4)') HEAD
    WRITE(1,'("mixregls_both_temp_.dat")')
    write(1,'("mixregls_both_temp_.out")')
    write(1,'("mixregls_both_temp_.def")')
    write(1,'("mixregls_both_temp_.ebv")')
        WRITE(1,'(5I3,E10.1E3, 9I3)') 1, 20, numloc+p+2, numloc, p, CONV, 0, 0, 0, 0, 1, 0, 1, 0, 0
        WRITE(1,'(20I3)') 1, 2
        write(1,'(20I3)') (2+p+j,j=1,numloc)
        write(1,'(20I3)') (2+j,j=1,p)
        write(1,*) ylabel
        write(1,*) (alab(j), j=1,numloc)
        write(1,*) (blab(j), j=1,p)

    rr=numloc*(numloc+1)/2
        allocate(beta(p))
        allocate(tau(s))
        allocate(mychol(RR))
        allocate(spar(max(ndim,1+ncov)))
        tau = 0
        spar = 0
        spar(1) = abs(log(sdlv))
    ALLOCATE (THETAs(NC2,ndim))
    ALLOCATE (thetavs(NC2,ndim2))

#if defined(_WIN32)
       CALL SYSTEM("MIXREG.EXE > mixregls_both_temp_")
#else
       CALL SYSTEM("./mixreg > mixregls_both_temp_")
#endif
       open(unit=1,file="mixreg.lik")
       read(1,*) logl, npar
       close(1)
       if(npar < p+rr+1) stop
       open(unit=1,file="mixreg.est")
       do i=1,npar
            read(1,*) templabel, temp
            if(i .le. p) then
                beta(i) = temp
            else if (i .le. p+rr) then
                mychol(i-p) = temp
                if(mls .eq. 1) v(i-p) = temp
            else if (i .eq. p+rr+1) then
                tau(1) = log(temp)
            end if
        end do
       close(1)

        if(mls .ne. 1) then
            allocate(alpha(r))
            alpha = 0
            alpha(1) = log(mychol(1))
        else

           open(unit=1,file="mixreg.var")
           do i=1,p
                read(1,*) (temp3(j), j=1,p)
                se(i) = sqrt(temp3(i))
           end do
           do i=p+1,rr+p+1
                read(1,*) (temp3(j), j=1,rr+1)
                se(i) = sqrt(temp3(i-p))
           end do
           close(1)
!           deallocate(temp3)

!           allocate(temp3(2+r+rr))
           open(unit=1,file="mixregls_both_temp_.ebv")
           do i=1,nc2
               read(1,*) (temp3(j),j=1,2+r+rr)
               do j=1,r
                   thetas(i,j) = temp3(2+j)
               end do
               do j=1,rr
                   thetavs(i,j) = temp3(2+r+j)
               end do
           end do
           close(1)
!           deallocate(temp3)

        IUN    = 16
        OPEN(UNIT=IUN,FILE="mixregls_both3.OUT")
                WRITE(IUN,*)
                WRITE(IUN,*)
                WRITE(IUN,'("------------------------------")')
                WRITE(IUN,'("Model without Scale Parameters")')
                WRITE(IUN,'("------------------------------")')

         ! WRITE RESULTS
    WRITE(IUN,562) LOGL,LOGL-NPAR, &
    LOGL-0.5D0*DBLE(NPAR)*DLOG(DBLE(NC2)),0-2*LOGL,0-2*(LOGL-NPAR), &
    0-2*(LOGL-0.5D0*DBLE(NPAR)*DLOG(DBLE(NC2)))
    562 FORMAT(1X,'Log Likelihood                 = ',F12.3,/, &
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


        WRITE(IUN,'("BETA (regression coefficients)")')
         DO i=1,P
            ZVAL = BETA(i)/SE(i)
            PVAL = 2.0D0 *(1.0D0 - PHIFN(DABS(ZVAL),0))
                WRITE(IUN,804)BLAB(i),BETA(i),SE(i),ZVAL,PVAL
         END DO
         write(iun,*)
        WRITE(IUN,'("Random (location) Effect Variances and Covariances")')
        counter = 1
         do i=1,r
             do j=1,i
                ZVAL = mychol(counter)/SE(counter+p)
                PVAL = 2.0D0 *(1.0D0 - PHIFN(DABS(ZVAL),0))
                if(i==j) then
                    WRITE(IUN,804)alab(i),mychol(counter),SE(counter+p),ZVAL,PVAL
                else
                    WRITE(IUN,805)'Covariance',j,i,mychol(counter),SE(counter+p),ZVAL,PVAL
                end if
                counter = counter+1
            end do
         END DO
        ZVAL = temp/SE(1+p+rr)
        PVAL = 2.0D0 *(1.0D0 - PHIFN(DABS(ZVAL),0))
        write(iun,*)
        WRITE(IUN,804)'Error variance  ',temp,SE(1+p+rr),ZVAL,PVAL

        804 FORMAT(A16,4(4x,F12.5))
        805 FORMAT(A10,I0,I0,4X,4(4x,F12.5))
        close(iun)
           se(npar) = se(npar)/temp
        open(2,file="mixregls_both_.est")
        WRITE(2,'(F15.6,2I8)') -2*LOGL,20,MAXIT
         WRITE(2,'(35F15.8)')(BETA(k),k=1,P)
         WRITE(2,'(35F15.8)')(mychol(k),k=1,RR)
         WRITE(2,'(35F15.8)')(TAU(1))
         WRITE(2,'(35F15.8)')(SE(k),k=1,1+p+rr)
         close(2)

!        v = mychol
           call CHSKY(v,mychol,r,myio)
    end if
end subroutine callmixreg

subroutine dostage2
    use lsboth
    use procedures
    implicit none
    CHARACTER(LEN=80) :: filetempdat,temp80
    integer :: i,n,j,k,diffk,totalvar,catcount(maxj),ni(nc2),myobs,ii,nvalid,&
                jlabel,jval,nvarpercat,nobs2,ncolinvar,nc2a
    REAL(KIND=8),ALLOCATABLE:: temp3(:),liks(:),tempdata(:),&
                                betas(:,:),meanbetas(:),totalvarhats(:),varbetas(:)
    real(kind=8) :: logl,loglsd,zval,pval,meany,miny,maxy,temp,stdy,seval
    character(len=24) :: mystr,progname
    REAL(KIND=8),ALLOCATABLE:: tempVector(:)

    ni = 1
    nvalid = 0
    nc2a = nc2
    write(*,*) "There are ",sum(allzeros),"subjects with unestimable random effects"
    if(multi2nd .eq. 1) then
        write(*,*) "Multilevel second stage"
        ni = idni(:,2)
        nobs2 = nobs
        if(sepfile .eq. 1) then
            j=0
            do i=1,nc2
                if(idni(i,1) .eq. idnisep(i+j,1)) then
                    ni(i) = idnisep(i+j,2)
                else
                    ni(i) = 0
                    j = j - 1
                    nc2a = nc2a - 1
                end if
                write(*,*) i,i+j, idni(i,1), idnisep(i+j,1),ni(i)
            end do
        end if
        nobs2 = sum(ni)
    else
        write(*,*) "Second stage"
        nobs2 = nc2
        if(sepfile .eq. 1) then
            j=0
            do i=1,nc2
                if(idni(i,1) .ne. idnisep(i+j,1)) then
                    ni(i) = 0
                    j = j - 1
                end if
            end do
            nobs2 = sum(ni)
        end if
        data2(1:nobs2,1:nvar2) = tempsums(1:nobs2,1:nvar2)
        nc2a = nobs2
    end if
    nvarpercat = nvar3-1+multi2nd
    write(*,*) nobs,nc2,nobs2,nc2a
        allocate(tempdata(nvar3))
        allocate(intLabel(nvar3+2+max(maxj-1,0)))
        intlabel(1) = var2label(1)
        intlabel(2) = "Intercept             "
        intlabel(3:pfixed+3) = var2label(2:1+pfixed)
        do j=1,numloc
            write(mystr, '(A6, I1, A17)') "Locat_",j,"                 "
            intLabel(1+Pfixed+1+(Ptheta+1)*(j-1)+1) = mystr
            do i=1,pTheta
                write(mystr, '(A6, I1, A17)') "Locat_",j,"*"//var2label(1+pfixed+i)
                intlabel(1+pfixed+1+(Ptheta+1)*(j-1)+1+i) = mystr
            end do
        end do
        if(pOmega >= 0) intlabel(1+pfixed+1+(1+Ptheta)*numloc+1) = "Scale                 "
        if(pOmega > 0) intlabel(1+pfixed+1+(1+Ptheta)*numloc+2:1+pfixed+1+(1+Ptheta)*numloc+1+pomega) = &
            "Scale*"//var2label(1+pfixed+ptheta+1:1+pfixed+ptheta+pomega)
        if(pTO >= 0) then
            do j=1,numloc
                write(mystr, '(A6, I1, A17)') "Locat_",j,"*Scale                 "
                intlabel(1+pfixed+1+(1+Ptheta)*numloc+1+pomega+(j-1)*(pto+1)+1) = mystr
                do k=1, pto
                    write(mystr, '(A6, I1, A17)') "Locat_",j,"*S*"//trim(var2label(1+pfixed+ptheta+pomega+k))
                    intlabel(1+pfixed+1+(1+Ptheta)*numloc+1+pomega+(j-1)*numloc+1+k) = mystr
                end do
            end do
        end if
        intLabel(nvar3+1) = "Random.Int.Var"
        intLabel(nvar3+1+multi2nd) = "Residual.Variance"
        if(stage2 .eq. 2 .and. maxj > 2) then
            intLabel(nvar3+maxj) = "Random.Int.Var"
            do j=1,maxj-1
                write(mystr,'(a10,f8.3,a6)') "Threshold ",icode(j+1),"      "
                intLabel(nvar3+j) = mystr
            end do
        end if

        open(3, file=trim(fileprefix)//'_desc2.out')
         ALLOCATE(tempVector(nobs2))
    200  FORMAT(1x,A16,4F12.4)
!    write(3,9) head
!    write(3,*)
    write(3,*) "Level 2 units =",nc2a
    if(multi2nd .eq. 1) write(3,*) "Level 1 observations =",nobs2
        WRITE(3,'("------------")')
         WRITE(3,'("Descriptives")')
         WRITE(3,'("------------")')
         WRITE(3,*)

        if(stage2 .eq. 1 .or. stage2 .eq. 3) then
            meany=SUM(data2(1:nobs2,1))/DBLE(nobs2)
             miny=minval(data2(1:nobs2,1))
             maxy=maxval(data2(1:nobs2,1))
             tempVector(:)=(data2(1:nobs2,1)-meany)**2
             temp=SUM(tempVector)/DBLE(nobs2-1)
             stdy=DSQRT(TEMP)
             WRITE(3,'(" Dependent variable")')
             WRITE(3,'("                         mean         min         max     std dev")')
             WRITE(3,'(" ----------------------------------------------------------------")')
             WRITE(3,200) var2Label(1),meany,miny,maxy,stdy
        else
            catcount = 0
            do i=1,nobs2
                do j=1,maxj
                    if(fp_equal(data2(i,1),icode(j))) then
                        catcount(j) = catcount(j) + 1
                        cycle
                    end if
                end do
            end do
             WRITE(3,'(" Categories of the Dependent Variable")')
             WRITE(3,'(" ----------------------------------------------------------------")')
             WRITE(3,'(" Category             Frequency         Proportion")')
            DO J = 1,MAXJ
                WRITE(3,"(1X,F8.2,8X,F12.2,8x,f12.5)") iCODE(J),real(catcount(J)),real(catcount(J))/nobs2
            end do
        end if

        if(pfixed > 0) then
            write(3,*)
            write(3,*)
             WRITE(3,'(" Independent variables")')
             WRITE(3,'("                         mean         min         max     std dev")')
             WRITE(3,'(" ----------------------------------------------------------------")')

            do i=1,pfixed
                meany=sum(data2(1:nobs2,1+i))/dble(nobs2)
                miny=minval(data2(1:nobs2,1+i))
                maxy=maxval(data2(1:nobs2,1+i))
                tempVector(:)=(data2(1:nobs2,1+i)-meany)**2
                TEMP=SUM(tempVector)/DBLE(nobs2-1)
                stdy=DSQRT(TEMP)
                WRITE(3,200) var2label(i+1),meany,miny,maxy,stdy
            end do
        end if

         WRITE(3,*)
         WRITE(3,*)
         WRITE(3,'(" Random Location and Scale EB mean estimates")')
         WRITE(3,'("                         mean         min         max     std dev")')
         WRITE(3,'(" ----------------------------------------------------------------")')
        do j=1,numloc+1-nors
            meany=sum(thetas(:,j),ni>0)/dble(nc2a)
            miny=minval(thetas(:,j),ni>0)
            maxy=maxval(thetas(:,j),ni>0)
            tempVector(:)=(thetas(:,j)-meany)**2
            TEMP=SUM(tempVector,ni>0)/DBLE(nc2a-1)
            stdy=DSQRT(TEMP)
            write(mystr, '(A6, I1, A17)') "Locat_",j,"                 "
            if(j .eq. numloc+1) mystr = "Scale"
            if(j .le. numloc .or. (j .eq. numloc+1 .and. nors .ne. 1)) WRITE(3,200) mystr,meany,miny,maxy,stdy
        end do
        if(pto >= 0) then
            do j=1,numloc
                meany=sum(thetas(:,j)*thetas(:,numloc+1),ni>0)/dble(nc2a)
                miny=minval(thetas(:,j)*thetas(:,numloc+1),ni>0)
                maxy=maxval(thetas(:,j)*thetas(:,numloc+1),ni>0)
                tempVector(:)=(thetas(:,j)*thetas(:,numloc+1)-meany)**2
                TEMP=SUM(tempVector,ni>0)/DBLE(nc2a-1)
                stdy=DSQRT(TEMP)
                write(mystr, '(A6, I1, A17)') "Locat_",j,"*Scale                 "
                WRITE(3,200) mystr,meany,miny,maxy,stdy
            end do
        end if
         WRITE(3,*)
         write(3,*) "There are",sum(allzeros),"subjects with unestimable random effect values"
         WRITE(3,*)
         WRITE(3,*)

         close(3)
            write(mystr, '(I5)') nreps


    totalvar = nvar3-1+multi2nd
    select case(stage2)
        case(1)
            progname = "mixreg"
            totalvar = nvar3+multi2nd
        case(2)
            !progname = "mixor"
            !if(multi2nd .eq. 1)
            progname = "mixors"
            totalvar = nvar3-1+multi2nd+maxj-1
        case(3)
            progname = "mixpreg"
        case(4)
            progname = "mixno"
            totalvar = (maxj-1)*(nvar3-1+multi2nd)
    end select

        filetempdat = trim(fileprefix) // "_.dat"
        OPEN(1,FILE=trim(progname)//".def")
        WRITE(1,*) " x                                                                             "
        WRITE(1,*) " x                                                                              "
        WRITE(1,5) adjustl(filetempdat)
        if(len(trim(fileprefix)) > 34) then
            temp80 = fileprefix(1:34) // "_x.out"
        else
            temp80 = trim(fileprefix) // "_x.out"
        end if
        WRITE(1,5) adjustl(temp80)
    5 FORMAT(A80)

        if(stage2 .ne. 2) WRITE(1,'("temp_.def")')

    select case(stage2)
        case(1)
            write(1,'("temp_.res")')
            WRITE(1,'(5I3,E10.1E3, 9I3)') 1, 5, nvar3+2, multi2nd, nvar3-1, .0005, 0, 0, 0, 0, 0, 0, 1, 0, 0
        case(2)
        !if(multi2nd .ne. 1) then
        !    WRITE(1,'(4I3,E10.1E3,17I3,i4)') 0, nvar3+2, multi2nd, nvar3-1, .0005, maxj, 0, 0,0,0,11,1,1,0,0,0,0,0,1,0,0,0,100
        !else
            WRITE(1,'(9I3,E10.1E3,4I4,E10.1E3,2i2)') nvar3+2, maxj, nvar3-2, 0, 0, 0, 0, 0, 1-multi2nd, .0005, 11,1,100,0,.1,1,1
            write(1,*) 0, -1
        !end if
        case(3)
            WRITE(1,'(4I3,E10.1E3, 12I3)') 1, nvar3+2, multi2nd, nvar3-1, .0005, 0, 0, 0, 0, 0, 0, 11, 0, 0, 1, 0, 0
        case(4)
            WRITE(1,'(4I3,E10.1E3,17I3,i4)') 0, nvar3+2, multi2nd, nvar3-1, .0005, maxj, 0, 0,0,0,0,0,11,0,0,1,0,1,0

    end select
        WRITE(1,*) 1, 2
        if(stage2 .ne. 2) then! .or. multi2nd .ne. 1) then
            if(multi2nd .eq. 1) write(1,*) nvar3+2
            write(1,*) (i+2, i=1, nvar3-1)
        else
            write(1,*) (i+2, i=2, nvar3-1)
        end if
        if(stage2 .eq. 2 .or. stage2 .eq. 4)     write(1,'(20f8.3)')(iCODE(J), J = 1,MAXJ)
        WRITE(1,*) intlabel(1)
        if(stage2 .ne. 2) then! .or. multi2nd .ne. 1) then
            if(multi2nd .eq. 1) write(1,*) "Intercept"
            if(stage2 .eq. 4) then
                write(1,'(12a8)') (intlabel(i), i=2, nvar3)
            else
                write(1,*) (intlabel(i), i=2, nvar3)
            end if
        else
            write(1,*) (intlabel(i), i=3, nvar3)
        end if
        CLOSE(1)

        allocate(betas(nreps,totalvar))
        allocate(totalvarhats(totalvar))
        totalvarhats = 0
        allocate(liks(nreps))
        allocate(temp3(totalvar))
        write(*,*) trim(progname),nvar3,multi2nd,totalvar
        betas = 0
        liks = 0
        ncolinvar = totalvar
        if(stage2 .eq. 1) ncolinvar = nvar3-1
    do n=1, nreps
        myobs = 0
        if(mod(n,10) .eq. 0) write(*,*) n
        open(n*5,file=filetempdat)
        DO I=1,NC2  ! go over level-2 clusters
            do ii=1,ni(i)
                myobs = myobs + 1
                tempdata(1) = data2(myobs,1)
                k=1
                diffk = 0
                if(pfixed .ge. 0) then
                    tempdata(2) = 1d0
                    diffk = diffk + 1
                    if(pfixed .ge. 1) then
                        tempdata(k+diffk+1:k+diffk+pfixed) = data2(myobs,k+1:k+pfixed)
                        k = k + pfixed
                    end if
                end if
                if(ptheta .ge. 0) then
                    do j=1,numloc
                        tempdata(k+diffk+1) = simvals(n,i,j)
                        diffk = diffk + 1
                        if(ptheta .ge. 1) then
                            tempdata(k+diffk+1:k+diffk+ptheta) = simvals(n,i,j)*data2(myobs,k+1:k+ptheta)
                            if(j .ne. numloc) diffk = diffk + ptheta
                        end if
                    end do
                    k = k + ptheta
                end if
                if(pomega .ge. 0) then
                    tempdata(k+diffk+1) = simvals(n,i,numloc+1)
                    diffk = diffk + 1
                    if(pomega .ge. 1) then
                        tempdata(k+diffk+1:k+diffk+pomega) = simvals(n,i,numloc+1)*data2(myobs,k+1:k+pomega)
                        k = k + pomega
                    end if
                end if
                if(pto .ge. 0) then
                    do j=1,numloc
                        tempdata(k+diffk+1) = simvals(n,i,j)*simvals(n,i,numloc+1)
                        diffk = diffk + 1
                        if(pto .ge. 1) then
                            tempdata(k+diffk+1:k+diffk+pto) = tempdata(k+diffk+1)*data2(myobs,k+1:k+pto)
                            if(j .ne. numloc) diffk = diffk + ptheta
                        end if
                    end do
                    k = k + pto
                end if
                write(n*5,'(i16,25e18.6)') idni(i,1),(tempdata(k),k=1,nvar3),1.0
            end do
        end do
        close(n*5)
#if defined(_WIN32)
        CALL SYSTEM(progname//"> _temp")
#else
        CALL SYSTEM("./"//progname//"> _temp")
#endif
        nvalid = nvalid + 1
        open(n*5+2, file=trim(progname)//".est")
        do j=1,totalvar
            read(n*5+2,*) mystr, betas(nvalid,j)
        end do
        close(n*5+2)
        if(isnan(betas(nvalid,1))) then
            nvalid = nvalid - 1
            write(*,*) "INVALID REPLICATION"
        else
            open(n*5+3, file=trim(progname)//".var")
                do j=1,ncolinvar
                    read(n*5+3,*) (temp3(i),i=1,ncolinvar)
                    totalvarhats(j) = totalvarhats(j) + temp3(j)
                end do
                if(stage2 .eq. 1) then
                    do j=1,1+multi2nd
                        read(n*5+3,*) (temp3(i),i=1,1+multi2nd)
                        totalvarhats(nvar3-1+j) = totalvarhats(nvar3-1+j) + temp3(j)
                    end do
                end if
            close(n*5+3)
            open(n*5+4, file=trim(progname)//".lik")
            read(n*5+4,*) liks(n)
            close(n*5+4)
        end if
    end do

            write(mystr, '(I5)') nvalid
        fileout = trim(fileprefix) // "_"//trim(adjustl(mystr)) //".out"
    open(2,file=fileout)
    write(2,*)
!write(2,*) "Level 2 obervations =",nc2
    write(2,*) "Number of replications =", nvalid
    write(2,*)
    WRITE(2,'("-------------")')
    WRITE(2,'("Final Results")')
    WRITE(2,'("-------------")')
    logl = sum(liks)/nvalid
    loglsd = sqrt(sum((liks(1:nvalid)-logl)**2)/(nvalid-1))
   WRITE(2,562)LOGL,loglsd,LOGL-totalvar+1, &
   LOGL-0.5D0*DBLE(totalvar-1)*DLOG(DBLE(NC2)),0-2*LOGL,0-2*(LOGL-totalvar+1), &
   0-2*(LOGL-0.5D0*DBLE(totalvar-1)*DLOG(DBLE(NC2)))
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
    allocate(varbetas(totalvar+1))
    allocate(meanbetas(totalvar+1))
    varbetas = 0
    totalvarhats(:) = totalvarhats(:) / (nvalid-1)

    do j=1,totalvar
        meanbetas(j) = sum(betas(1:nvalid,j))/nvalid
    end do
     804 FORMAT(A16,4(4x,F12.5))
        681 FORMAT(1X,'Category',F8.3,' vs Category',F8.3,/)
    do j=1,totalvar !nvar-1 fixed effects plus the residual variance
        if(stage2 .eq. 4 .and. mod(j,nvar3-1+multi2nd) .eq. 1)   write(2,681) icode(j/(nvar3-1)+2),icode(1)

        do i=1,nvalid
            varbetas(j) = varbetas(j) + (betas(i,j)-meanbetas(j))**2
        end do
        varbetas(j) = varbetas(j) / (nvalid - 1)
        jlabel = j+1
        jval = j
        if(stage2 .eq. 4) then
            jlabel = mod(j,nvarpercat)+1
            jval = j - j/nvarpercat
            if(jlabel .eq. 1) then
                jlabel = nvarpercat+1
                jval = totalvar - maxj+1+j/nvarpercat
            end if
        end if
        if(intLabel(jlabel) .eq. "Random.Int.Var") write(2,*)
        seval = sqrt(totalvarhats(jval)+varbetas(jval)*(1+1/nvalid))
        ZVAL = meanbetas(jval)/seval
        PVAL = 2.0D0 *(1.0D0 - PHIFN(DABS(ZVAL),0))
        if(stage2 .ne. 2) then
            write(2,804) intlabel(jlabel), meanbetas(jval), seval,zval,pval
        else if(maxj .eq. 2) then
            if(j .ne. nvar3) write(2,804) intlabel(jlabel), meanbetas(jval), seval,zval,pval
        else
            if(j .ne. 1) write(2,804) intlabel(jlabel), meanbetas(jval), seval,zval,pval
        end if
        if(stage2 .eq. 4 .and. mod(j,nvar3-1+multi2nd) .eq. 0) write(2,*)
    end do

#if defined(_WIN32)
        CALL SYSTEM("copy "//trim(fileprefix)//"_desc2.out+"//trim(fileprefix) &
                    //"_"//trim(adjustl(mystr))//".out "//trim(fileprefix)//"_2.out")
#else
        CALL SYSTEM("cat "//trim(fileprefix)//"_desc2.out "//trim(fileprefix) &
            //"_"//trim(adjustl(mystr))//".out >> "//trim(fileprefix)//"_2.out")
#endif
!!    call system("del "//trim(fileprefix)//"_desc2.out")
!    CALL SYSTEM("DEL temp_.*")
end subroutine dostage2

subroutine make_random
    use lsboth
    use procedures
    implicit none

    INTEGER :: I,k,j,numre,c
    REAL(KIND=8),ALLOCATABLE:: chols(:,:,:),myrand(:,:)
    integer,allocatable :: tempseed(:)

    numre = numloc+(1-nors)
    allocate(chols(nc2,numre,numre))
    allocate(simvals(nreps,nc2,numre))
    allocate(myrand(nc2,numre))
    allocate(allzeros(nc2))

    call random_seed(size=k)
    allocate(tempseed(k))
    tempseed(:) = myseed
    call random_seed(put=tempseed)

    allzeros = 0
    do k=1,nc2
        c=1
        do i=1,numre
            do j=1,i
                chols(k,i,j) = thetavs(k,c)
                c = c+1
                chols(k,j,i) = chols(k,i,j)
            end do
        end do

        call cholesky_sub(chols(k,:,:),numre,allzeros(k))
!        write(*,*) ((chols(k,i,j),i=1,numre),j=1,numre)
    end do
!        write(*,'(20f8.3)') ((chols(1,i,j),i=1,numre),j=1,numre)
!        write(*,'(20f8.3)') (thetavs(1,c),c=1,numre*(numre+1)/2)
    do i=1,nreps
        do j=1,nc2
            do k=1,numre
                myrand(j,k) = random_normal()
            end do
        end do
        do k=1,nc2
            simvals(i,k,:) = matmul(myrand(k,:),chols(k,:,:)) + thetas(k,:)
!            write(*,*) k,(simvals(i,j,k),j=1,nc2)
        end do
    end do
end subroutine make_random

SUBROUTINE READAT2()
    use lsboth, only: tempsums,data2,nc2sep,nvar2,idnisep,id2indsep,miss,ymiss,nvarsep,maxksep,nobssep,var2ind,&
                      filedat2
    use procedures
    implicit none

    INTEGER :: myPASS,I,K,ICOUNT,myindex,IDTEMP,IDOLD,hasmiss
    REAL(KIND=8),ALLOCATABLE:: TEMPR(:)
    LOGICAL FIRST

    ALLOCATE (TEMPR(NVARsep))

   ! INITIALIZE
    DO myPASS = 1,2
        IF (myPASS .EQ. 2) THEN
                allocate(tempsums(nc2sep,nvar2))
                allocate(data2(icount,nvar2))

              tempsums = 0
              data2 = 0
        ! IDNI has IDs and Nobs per ID
              ALLOCATE (IDNIsep(NC2sep,2))
              IDNIsep = 0
        ENDIF

        I     = 1
        K     = 1
        MAXKsep  = 0
        ICOUNT= 0
        NOBSsep  = 0
        FIRST  = .TRUE.

        ! READ IN DATA UNTIL END
        OPEN(1,ACTION='READ',FILE=FILEDAT2)

        DO   ! loop forever

              READ(1,*,END=1999)(TEMPR(myindex),myindex=1,NVARsep)
                hasmiss = 0
                IF (MISS .EQ. 1) THEN
                    do myindex = 1,nvar2
                        IF (FP_EQUAL(tempr(var2ind(myindex)), YMISS)) THEN
                            hasMISS = 1
                            exit
                        END IF
                    end do
                end if
                IF (hasMISS .NE. 0) THEN
                    CYCLE  ! give up on current value and go read next one
                end if

      ! QUERY FOR NEW ID AND SET PARAMETERS ACCORDINGLY

              IDTEMP = INT(TEMPR(ID2INDsep))
              ! QUERY FOR NEW ID AND SET PARAMETERS ACCORDINGLY

              IF (.NOT. FIRST) THEN
                 IF (IDTEMP .EQ. IDOLD) THEN
                    K     = K+1
                 ELSE
                    if(mypass .eq. 2) then
                        IDNIsep(I,1) = IDOLD
                        IDNIsep(I,2) = K
                    ENDIF
                        NOBSsep = NOBSsep+K
                        IF (K .GT. MAXKsep) MAXKsep = K
                        I     = I+1
                    K     = 1
                 ENDIF
              ENDIF

              ! PUT TEMPORARY VALUES INTO DATA VECTORS AND MATRICES

              IDOLD = IDTEMP
              ICOUNT = ICOUNT+1

              FIRST  = .FALSE.
              IF (myPASS == 2) THEN
                    do myindex=1,nvar2
                        if(var2ind(myindex) .ne. 0) then
                            tempsums(i,myindex) = tempsums(i,myindex) + tempr(var2ind(myindex))
                            data2(icount,myindex) = tempr(var2ind(myindex))
                        end if
                    end do
              END IF

        END DO   ! loop back to read next line

    ! cleanup final entry
    1999                IF (myPASS .EQ. 2) THEN
                           IDNIsep(I,1) = IDOLD
                           IDNIsep(I,2) = K
                        ENDIF
                        NOBSsep = NOBSsep+K
                        IF (K .GT. MAXKsep) MAXKsep = K
    NC2sep = I
    CLOSE(1)
   END DO   ! two passes, one to get size, second to read data
        do i=1,nc2sep
            tempsums(i,1:nvar2) = tempsums(i,1:nvar2)/idnisep(i,2)
        end do
   DEALLOCATE(TEMPR)
END SUBROUTINE READAT2

subroutine run_stage2()
    use lsboth
    implicit none
    integer::i,j,k

        open(1,file="stage2only.def")
    WRITE(1,'(18A4)') HEAD
        if(sepfile .ne. 1) then
            filedat2 = filedat
            nvarsep = nvar
            id2indsep = id2ind
        end if
        write(1,*)FILEDAT2
        write(1,*)trim(fileprefix)//"_ebvar.dat"
        write(1,*)trim(FILEprefix)//"_stage2"
        write(1,*) NVARsep,nc2,numloc,nreps, nors, myseed, stage2, multi2nd,Ymiss
        write(1,*) pfixed,ptheta,pomega,pto
       if(readcats .eq. 1) then
            write(1,*) maxj
            write(1,*)(iCODE(J), J = 1,MAXJ)
        end if
        write(1,*) id2indsep,var2ind(1)
        k = 1
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
            k = k + pomega
         END IF
        IF (Pto .GE. 1) THEN
            write(1,*) (var2IND(k+I), I=1,Pto)
         END IF

         write(1,*) var2label(1)
        k = 1
        IF (Pfixed .GE. 1) THEN
            write(1,*) (var2label(k+I), I=1,Pfixed)
            k = k + pfixed
         END IF
        IF (Ptheta .GE. 1) THEN
            write(1,*) (var2label(k+I), I=1,Ptheta)
            k = k + ptheta
        END IF
        IF (Pomega .GE. 1) THEN
            write(1,*) (var2label(k+I), I=1,Pomega)
            k = k + pomega
         END IF
        IF (Pto .GE. 1) THEN
            write(1,*) (var2label(k+I), I=1,Pto)
        END IF
#if defined(_WIN32)
        call system("stage2only")
#else
        call system("./stage2only")
#endif
    CLOSE(1)
end subroutine run_stage2
