!  **************************************************************
!  MIXOR with Random Scale
!
! *************************************************************
!                                                                            

MODULE mixor_globals
    IMPLICIT NONE
    save
    INTEGER :: NVAR,NQ,AQUAD,ID2IND,YIND,P,R,S,RR,rnint,pv,rv,sv,nv, &
               MISS,MAXj,NS,MAXIT,myio,npar,myqdim,mytotalq,pold,sold,rold,nors,&
               rr1,chol,nfn,ngamma,qtotalR0,qtotalR1,nobs,nclust,iadd,s_cycle,ns_cycle,npar_cycle,rr_cycle,&
               myseed,stage2,multi2nd,pfixed,ptheta,pomega,pto,readcats,nvar3,maxj2,nvar2,nreps,&
               numloc,sepfile,nvarsep,id2indsep,nalpha,mls,ncov
    INTEGER,ALLOCATABLE :: XIND(:),UIND(:),WIND(:),Yasint(:),ids(:),ni(:), cluster_nob(:,:), catcount(:),&
                            varind(:),var2ind(:)
    REAL(KIND=8) :: RIDGEIN,CONV,YMISS
    REAL(KIND=8),ALLOCATABLE:: Y(:),X(:,:),U(:,:),W(:,:), mychol(:),mygamma(:),se(:),var(:,:),varavg(:,:),&
                               BETA(:),TAU(:),SPAR(:),thetas(:,:),thetavs(:,:,:),jcodes(:),pointsR0(:,:),&
                               pointsR1(:,:),weightsR0(:),weightsR1(:),mypoints(:,:),myweights(:),meanx(:),&
                               jcodes2(:), alpha(:)
    CHARACTER(LEN=24) :: YLABEL, templabel
    CHARACTER(LEN=24),ALLOCATABLE :: BLABel(:),tLABel(:),clabel(:),varlabel(:),var2label(:)
    CHARACTER(LEN=80) :: FILEDAT, head1, head2, fileprefix,filedat2
    REAL(kind=8), PARAMETER :: Pi = 3.141592653589793238462643
end module mixor_globals

PROGRAM MIXorS
    use mixor_globals
    implicit none

    call readDef()
    call writeDef()
    call readData()
    call adjustData()
    call printdesc()
    call startv2()    
    call mixorEst()
#if defined(_WIN32)
    CALL SYSTEM("COPY mixors1.OUT+mixors_both.OUT " // trim(fileprefix)//"_stage1.out")
    if(stage2 .ne. 0) call run_stage2()
    call system("mkdir work")
    call system("move mixors_both_* work")
#else
    CALL SYSTEM("cat mixors1.OUT mixors_both.OUT >> " // trim(fileprefix)//"_stage1.out")
    if(stage2 .ne. 0) call run_stage2()
    call system("mkdir work")
    call system("mv mixors_both_* work")
#endif
end program mixors

subroutine run_stage2()
    use mixor_globals
    implicit none
    integer::i,j,k
    
        open(1,file="stage2only.def")
    WRITE(1,"(A80)") HEAD1
    WRITE(1,"(A80)") HEAD2
        if(sepfile .ne. 1) then
            filedat2 = filedat
            nvarsep = nvar
            id2indsep = id2ind
        end if
        write(1,*)FILEDAT2
        write(1,*)trim(fileprefix)//"_ebvar.dat"
        write(1,*)trim(FILEprefix)//"_stage2"
        write(1,*) NVARsep,nclust,numloc,nreps, nors, myseed, stage2, multi2nd,Ymiss
        write(1,*) pfixed,ptheta,pomega,pto
       if(readcats .eq. 1) then
            write(1,*) maxj2
            write(1,*)(jCODEs2(J), J = 1,MAXJ2)
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
     

subroutine readDef
    use mixor_globals
    implicit none
    integer :: i,j,k
    logical :: fp_equal
    OPEN(1, FILE='mixors_random_mixblank.def')
    READ(1,"(A80)")head1
    READ(1,"(A80)")head2
    READ(1,*)FILEDAT
    READ(1,*)FILEprefix

    READ(1,*) NVAR, maxj, P, R, S, pv, rv, sv, rnint, CONV, NQ, AQUAD, MAXIT, yMISS, ridgein, nfn, nors,&
                myseed, stage2, multi2nd, nreps, sepfile, mls, ncov
    READ(1,*) chol, iadd
    IF (IADD .NE. 1 .AND. IADD .NE. -1) IADD = -1
    miss = 1
    if(fp_equal(ymiss, 0.0d0)) miss = 0
    if(nors .ne. 1) nors = 0
    if(mls .ne. 0) then
        mls = 1
        numloc = R + 1-rnint + rv
        nalpha = 0
        ns = numloc+1
        ncov = 1
    else
        numloc = 1
        nalpha = R
        ns = max(min(ncov+1,3),1)
    end if
    
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

    p = p + 2*pv
    !Special case, random slope term requires occasion-level variables
    if(mls .eq. 1) then
        r = r + 1 - rnint + rv
    else
        r = r + 1 - rnint
        nalpha = nalpha + 2*rv + 1 - rnint
    end if
    rr = numloc*(numloc+1)/2
    s = s + 2*sv
    allocate(jcodes(maxj))
    READ(1,*)(jCODEs(i), i = 1,MAXJ)
      
    READ(1,*) YLABEL
     IF (P .GE. 1) THEN
        ALLOCATE(BLABel(P))
        if(pold > 0) READ(1,*) (BLABel(I), I=1,Pold)
     END IF
     IF (R .GE. 1) THEN
        ALLOCATE(cLABel(R))
        if(rold > 0) READ(1,*) (cLABel(I+1-rnint), I=1,Rold)
     END IF
     IF (S .GE. 1) THEN
        ALLOCATE(TLABel(S))
        if(sold > 0) READ(1,*) (TLABel(I), I=1,Sold)
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
    ngamma = maxj-1
    NPAR = P+RR*mls+nalpha*(1-mls)+ngamma+S+NS ! number of parameters
    
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
            read(1,*) maxj2
            allocate(jcodes2(maxj2))
            READ(1,*)(jCODEs2(J), J = 1,MAXJ2)
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

end subroutine readDef

subroutine writeDef
    use mixor_globals
    implicit none
    integer :: i,j,k
    OPEN(1,FILE=trim(fileprefix)//".def")
    
    WRITE(1,"(A80)") HEAD1
    WRITE(1,"(A80)") HEAD2
    WRITE(1,"(A80)")FILEDAT
    WRITE(1,"(A80)")FILEprefix
    write(1,*) NVAR, maxj, Pold, Rold, Sold, rnint, pv, rv, sv, CONV, NQ, AQUAD, MAXIT, yMISS, ridgein, nfn, nors,&
                myseed, stage2, multi2nd, sepfile, mls, ncov
    write(1,*) chol, iadd
             
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

    WRITE(1,*)(jCODEs(i), i = 1,MAXJ)

    WRITE(1,*) YLABEL
    IF (P .GE. 1) THEN
        WRITE(1,*) (BLABel(I), I=1,Pold)
    END IF
    IF (R .GE. 1) THEN
        WRITE(1,*) (cLABel(I+1-rnint), I=1,Rold)
    END IF
    IF (S .GE. 1) THEN
        WRITE(1,*) (TLABel(I), I=1,Sold)
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
        write(1,*) pfixed,ptheta,pomega,pto

        if(readcats .eq. 1) then
            write(1,*) maxj2
            write(1,*)(jCODEs2(J), J = 1,MAXJ2)
        end if
        if(sepfile .eq. 1) write(1,*) filedat2
        if(sepfile .eq. 1) write(1,*) nvarsep,id2indsep
        write(1,*) var2ind(1)
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

    ! read in the labels
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
end subroutine writeDef

subroutine readData
    use mixor_globals
    implicit none
    
    INTEGER :: myPASS,I,K,ICOUNT,myindex,IDTEMP,IDOLD,hasmiss,nvartotal
    REAL(KIND=8),ALLOCATABLE:: TEMPR(:)
    INTEGER,ALLOCATABLE :: allvarsind(:)
    LOGICAL FIRST, FP_EQUAL

    ALLOCATE (TEMPR(NVAR))
        nvarTotal = 1+pold+rold+sold+nv
        allocate(allvarsIND(nvarTotal))
        allvarsIND(1) = yind
        allvarsIND(2:pold+1) = xind(1:pold)
        allvarsIND(pold+2:pold+rold+1) = uind(1:rold)
        allvarsIND(pold+rold+2:pold+rold+sold+1) = wind(1:sold)
        if(nv > 0) allvarsIND(pold+rold+sold+2:nvarTotal) = varind(1:nv)

   ! INITIALIZE
    DO myPASS = 1,2
        IF (myPASS .EQ. 2) THEN
               ALLOCATE (Y(ICOUNT))
                ALLOCATE (X(ICOUNT,P))
                ALLOCATE (U(ICOUNT,max(numloc,nalpha)))
                ALLOCATE (W(ICOUNT,S))

              Y = 0.0D0
              X = 0.0D0
              U = 0.0D0
              W = 0.0D0
                if(nv > 0) then
                allocate (var(icount,nv))
                allocate (varavg(nclust,nv))
              var = 0.0D0
              varavg = 0
            end if
        ! IDNI has IDs and Nobs per ID
              allocate(ni(nclust))
              ni = 0
              allocate(Ids(nclust))
        ENDIF
   
        I     = 1
        K     = 1
        ICOUNT= 0
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
                    IF (myPASS .EQ. 2) THEN
                        ids(i) = idold
                        ni(i) = k
                    ENDIF
!                    NOBS = NOBS+K
                    K     = 1
                    I     = I+1
                 ENDIF
              ENDIF

              ! PUT TEMPORARY VALUES INTO DATA VECTORS AND MATRICES

              IDOLD = IDTEMP
              ICOUNT = ICOUNT+1

              FIRST  = .FALSE.
              IF (myPASS == 2) THEN
                    Y(icount)  = TEMPR(YIND)
                    x(icount,1:p) = tempr(xind(1:pold))
                    do myindex=1,rold
                        u(icount,myindex+1-rnint) = tempr(uind(myindex))
                    end do
                    do myindex=1,sold
                        w(icount,myindex) = tempr(wind(myindex))
                    end do
                    do myindex=1,nv
                        var(icount,myindex) = tempr(varind(myindex))
                        varavg(i,myindex) = varavg(i,myindex) + tempr(varind(myindex))
                    end do
              END IF

    END DO   ! loop back to read next line

    ! cleanup final entry
    1999   IF (myPASS .EQ. 2) THEN
              IDs(i) = IDOLD
              NI(I) = K
           ENDIF
!    NOBS = NOBS+K
!    IF (K .GT. MAXK) MAXK = K

    NClust = I
    nobs = icount
    CLOSE(1)
    END DO   ! two passes, one to get size, second to read data
    if(nv > 0) then
        do i=1,nclust
            varavg(i,1:nv) = varavg(i,1:nv)/ni(i)
        end do
    end if
   DEALLOCATE(TEMPR)
END SUBROUTINE readdata

subroutine adjustData
    use mixor_globals
    implicit none
    integer:: j,ll,h,ko,i,kv
    logical :: fp_equal
    if(rnint .ne. 1) then
        u(:,1) = 1
        clabel(1) = "Intercept   "
    end if
    do j=1, pv
        ll = pold
        blabel(ll+j*2-1) = trim(varlabel(j)) // "_BS"
        blabel(ll+j*2) = trim(varlabel(j)) // "_WS"
    end do
    do j=1, rv
        ll = rold + 1 - rnint
        if(mls .eq. 1) then
            clabel(ll+j) = trim(varlabel(j+pv)) // "_WS"
        else
            clabel(ll+j*2-1) = trim(varlabel(j+pv)) // "_BS"
            clabel(ll+j*2) = trim(varlabel(j+pv)) // "_WS"
        end if        
    end do
    do j=1, sv
        ll = sold
        tlabel(ll+j*2-1) = trim(varlabel(j+pv+rv)) // "_BS"
        tlabel(ll+j*2) = trim(varlabel(j+pv+rv)) // "_WS"
    end do
    ko = 0
    do i=1,nclust
        do h = 1,ni(i)
            ko = ko + 1
            do j=1, pv
                ll = pold
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
                ll = sold
                w(ko,ll+j*2-1) = varavg(i,j+kv)
                w(ko,ll+j*2) = var(ko,j+kv) - varavg(i,j+kv)
            end do
        end do
    end do
    allocate(yasint(nobs))
    yasint = -1
    allocate(catcount(maxj))
    catcount = 0
    do i=1,nobs
        do j=1,maxj
            if(fp_equal(y(i),jcodes(j))) then
                yasint(i) = j
                catcount(j) = catcount(j) + 1
                cycle
            end if
        end do
        if(fp_equal(y(i),-1.0)) write(*,*) "Observation did not match any category!", y(i), (jcodes(j),j=1,maxj)
    end do
end subroutine adjustData

SUBROUTINE printDesc()
        use mixor_globals
        implicit none

        INTEGER:: I,j
        real(kind=8):: catpct(maxj),cumprob(maxj),stdx(p),minx(p),maxx(p),&
                        meanu(r),stdu(r),minu(r),maxu(r),meanw(s),stdw(s),minw(s),maxw(s)
    allocate(meanx(p))

! get descriptive stats on all of the variables
    if (p>0) then
        call descripc(x,nobs,p,meanx,stdx,minx,maxx)
    end if
    if (r>0) then
        call descripc(u,nobs,r,meanu,stdu,minu,maxu)
    end if
    if (s>0) then
        call descripc(w,nobs,s,meanw,stdw,minw,maxw)
    end if
     OPEN(UNIT=17,FILE="mixors1.OUT")
     
     WRITE(17,'("MixorS: Mixed-effects (multiple) Location Ordinal Scale Model")')
     write (17,*)
     WRITE(17,'("-----------------------------")')
     WRITE(17,'("MixorS.DEF specifications")')
     WRITE(17,'("-----------------------------")')
     WRITE(17,"(1x,a80)")HEAD1
     WRITE(17,"(1x,a80)")HEAD2
     WRITE(17,*)
     WRITE(17,'(" data and output files:")')
     WRITE(17,"(1x,a80)")FILEDAT
     WRITE(17,*)trim(fileprefix)//".out"
     WRITE(17,*)
     WRITE(17,"(' CONVERGENCE CRITERION = ',F11.8)")CONV
     WRITE(17,"(' RIDGEIN    = ',F8.4)")RIDGEIN
     WRITE(17,"(' NQ         = ',I4)")NQ
     WRITE(17,"(' QUADRATURE = ',I4,' (0=non-adaptive, 1=adaptive)')")AQUAD
     WRITE(17,"(' MAXIT      = ',I4)")MAXIT

     WRITE(17,*)
     WRITE(17,*)
     WRITE(17,'("------------")')
     WRITE(17,'("Descriptives")')
     WRITE(17,'("------------")')
     WRITE(17,*)
     WRITE(17,'(" Number of level-1 observations = ",I8)') nobs
     WRITE(17,*)
     write(17,'(" Number of level-2 clusters     = ",I8)') nclust
     write(17,*)
     write(17,'(" Number of level-1 observations for each level-2 cluster:")')
     write(17,'(1x,13I6)') (ni(i),i=1,nclust)

200  FORMAT(1x,A24,4F12.4)

    WRITE(17,258) YLabel
    258 FORMAT(//,1x,'Categories of the response variable ',A24, &
               /,1x,'-------------------------------------------------------------',/)
    WRITE(17,457)
    457 FORMAT(1X,'Category',5X,'   Frequency',5x,'  Proportion',/)
    cumprob = 0   
    DO J = 1,MAXJ
        CATpct(J) = real(catcount(J)) / nobs
        WRITE(17,"(1X,F8.2,5X,F12.2,5x,f12.5)") jCODEs(J),real(catcount(J)),catpct(J)
    end do
    write(17,*)

     if (p>0) then
        WRITE(17,'(" Mean model covariates")')
        WRITE(17,'("                                 mean         min         max     std dev")') 
        WRITE(17,'(" ------------------------------------------------------------------------")')
        do i=1,p
           WRITE(17,200) BLABel(i),meanx(i),minx(i),maxx(i),stdx(i)
        end do
        WRITE(17,*)
     end if

     if (r>0) then
        WRITE(17,'(" BS variance variables")')
        WRITE(17,'("                                 mean         min         max     std dev")') 
        WRITE(17,'(" ------------------------------------------------------------------------")')
        do i=1,max(r,nalpha)
           WRITE(17,200) cLABel(i),meanu(i),minu(i),maxu(i),stdu(i)
        end do
        WRITE(17,*)
     end if

     if (s>0) then
        WRITE(17,'(" WS variance model covariates")')
        WRITE(17,'("                                 mean         min         max     std dev")') 
        WRITE(17,'(" ------------------------------------------------------------------------")')
        do i=1,s
           WRITE(17,200) TLABel(i),meanw(i),minw(i),maxw(i),stdw(i)
        end do
        WRITE(17,*)
     end if
     close(17)
end subroutine printdesc

SUBROUTINE STARTV2()
    use mixor_globals
    implicit none
    integer:: i,j,k
    real(kind=8) :: cumprob,lncumodds,myint,tempparam
    real(kind=8), allocatable :: xpx(:,:), xpxi(:,:), catpct(:),ytemp(:,:),xpy(:,:),x1(:,:)
    
    allocate(mygamma(ngamma))
    allocate(beta(p))
    allocate(mychol(rr))
    allocate(alpha(nalpha))
    allocate(tau(s))
    allocate(spar(ns))

    allocate(x1(nobs,p+1))
    allocate(xpx(p+1,p+1))
    allocate(xpxi(p+1,p+1))
    allocate(ytemp(nobs,1))
    allocate(xpy(p+1,1))
    x1 = 1
    x1(:,1:p) = x
    xpx = matmul(transpose(x1(1:nobs,:)), x1(1:nobs,:))
    call inverse(xpx, xpxi, p+1)
    ytemp(:,1) = y(1:nobs)
    xpy = matmul(transpose(x1(1:nobs,:)), ytemp)
    xpy = matmul(xpxi, xpy)
    beta = xpy(1:p,1)
    myint = 0!dot_product(beta, meanx)
    !Need to calculate xb from linear regression, and use where currently says xb
    allocate(catpct(maxj))
    cumprob = 0   
    DO J = 1,MAXJ
        CATpct(J) = real(catcount(J)) / nobs
        if(j < maxj) then
            cumprob = cumprob + catpct(j)
            lncumodds = log(cumprob) - log(1-cumprob)
            IF (nfn .eq. 0) THEN
                tempparam = 0.625 * (lncumodds)
            ELSEIF (nfn .eq. 1) THEN
                tempparam = lncumodds
            ELSEIF (nfn .eq. 2) THEN
                tempparam = DLOG(0.0D0 - DLOG(1.0d0 - cumprob))
            ELSEIF (nfn .eq. 3) THEN
                tempparam = DLOG(0.0D0 - DLOG(cumprob))
            ENDIF
            mygamma(j) = tempparam - myint*iadd
        end if 
    END DO
    tau(:) = 0
    if(s>0) tau(1) = 1
    mychol(:) = 0
    if(mls .eq. 1) then
        k = 0
        DO I  = 1,R
            DO j = 1,I
                k = k + 1
                IF (j .EQ. I) THEN
                    IF (I .EQ. 1) THEN
                        mychol(k) = 1.0d0
                    ELSE
                        mychol(k) = 0.5d0
                    ENDIF
                ENDIF
            END DO
        END DO
        IF (nfn .EQ. 1) mychol(:) = mychol(:)*pi/DSQRT(3.0D0)
        IF (nfn .GE. 2) mychol(:) = mychol(:)*pi/DSQRT(6.0D0)
    else
        alpha = 0
        alpha(1) = 2*log(pi)-log(3.0)
    end if
    deallocate(xpx,ytemp,xpxi,catpct,xpy,meanx)
END SUBROUTINE STARTV2
  


! ----------------------------------------------------------------------------
!

SUBROUTINE MIXorEST()
    use mixor_globals
    implicit none

    INTEGER :: I,J,L,l2,k,t,m,k2,kk,LL,ii,jj,jj1,&
               CYCLES,NCYCLE,IFIN,ITER,iun, &
               Q,NOB,RIDGEIT,IUNS,ndim,mynob,h,cholamt,myorder(rr+r+1),myorder2(rr+r+1)
    REAL(KIND=8) :: RIDGE,LOGLP,LOGL, &
                    LOGDIFF,MAXCORR,myscale,tempC,tempx(npar,npar),&
                    ORIDGE,PVAL,ZVAL,phiY, cdf_upper, cdf_lower, pdf_upper, pdf_lower, &
                    MAXDER,PHIFN,sdev,phiRatio,phider,ddf_upper,ddf_lower,&
                    tauhat, tauhatlow, tauhatup, myz, lambda_upper, lambda_lower,iprob,qprob,qlogprob,jprob,&
                    COREC(npar), temp(npar,npar), myider2(npar,npar),minder2,&
                    myqderp2(nobs,npar,npar),myider(npar),myqder2(npar,npar),nopoints(1,1),noweights(1),sigmaalpha
    real(kind=8),allocatable::THETA1(:,:),thetav(:,:), myweights0(:),mypoints0(:,:),corec1(:,:),myqderp1(:,:,:),&
                                dlambda_upper1(:,:),dlambda_lower1(:,:),myqder1(:,:),myider1(:,:),myder2(:,:),&
                               sstar(:,:), work(:,:), sigma(:), asstar2(:,:), adjVar(:,:), &
                               qprobs(:),myder(:),myqlogprob(:),mycholspar(:)


    ! parameters
    ndim = numloc+1
    IUN    = 16
    IUNS    = 17
    OPEN(UNIT=IUNS,FILE="MIXorS_both_details.ITS")
    open(19,file="MixorS_both_all.var")
    open(13,file="MixorS_both_all.est")
    if(aquad .ne. 0) open(29,file=trim(fileprefix)//"_ebvar.dat")
    qtotalR0 = nq ** numloc
    qtotalR1 = nq ** (numloc+1)
    allocate(pointsR0(qtotalR0,numloc))
    allocate(pointsR1(qtotalR1,numloc+1))
    allocate(weightsR0(qtotalR0))
    allocate(weightsR1(qtotalR1))
    call getquad(numloc, nq, qtotalR0, pointsR0, weightsR0)
    call getquad(numloc+1, nq, qtotalR1, pointsR1, weightsR1)
    nopoints(1,1) = 0
    noweights(1) = 1

    ! ALLOCATE VECTORS FOR NR (these don't vary by CYCLES)

    !Establishing arrays at maximum size needed
    ! NS = number of additional var cov parameters due to random SCALE
    myqdim = ndim
    allocate(myder(npar))
    allocate(myder2(npar,npar))
    allocate(dlambda_upper1(npar,1))
    allocate(dlambda_lower1(npar,1))
    allocate(myider1(npar,1))
    allocate(corec1(npar,1))
    allocate(myqder1(npar,1))
    allocate(myqderp1(nobs,npar,1))
    ALLOCATE (SE(NPAR))
    ALLOCATE (THETA1(ndim,1))
    allocate(thetas(nclust,ndim))
    ALLOCATE (thetav(ndim,ndim))
    allocate(qprobs(qtotalR1))
    allocate(thetavs(nclust,ndim,ndim))
    allocate(myqlogprob(nobs))

    cholamt=rr+ns
    allocate(sstar(cholamt,cholamt))
    allocate(work(npar,npar))
    allocate(adjVar(npar,npar))
    allocate(asstar2(npar,npar))
    allocate(sigma(cholamt))
    allocate(mycholspar(cholamt))
    
! start cycles
! cycles = 1: random intercept model with BS variance terms
! cycles = 2: add in scale (WS) variance terms
! cycles = 3: add in random scale 
! cycles = 4: use NS = R+1
    open(unit=iun,file="mixors_both.out",status="replace")
    close(iun)

    ncycle = 5
    write(IUNS,*)'MULTIPLE LOCATION EFFECTS =',mls==1
    CYCLELOOP:do cycles=1,ncycle
        OPEN(UNIT=IUN,FILE="MIXorS_both.OUT",access="append")
        if(cycles==2 .AND. r==0) cycle
        if(cycles==3 .AND. s==0) cycle
        if(cycles > 3 .and. nors==1) cycle
        if(cycles .eq. 5 .and. ncov==0) cycle

            WRITE(IUN,*)
            WRITE(IUN,*)
            WRITE(*,*)
            WRITE(IUNS,*)
             
        if (cycles==1) then
            write(IUN,*)'MULTIPLE LOCATION EFFECTS =',mls==1
            WRITE(IUN,'("----------------------------")')
            WRITE(IUN,'("Model without Random Effects")')
            WRITE(IUN,'("----------------------------")')
            WRITE(*,'("----------------------------")')
            WRITE(*,'("Model without Random Effects")')
            WRITE(*,'("----------------------------")')
            WRITE(IUNS,'("----------------------------")')
            WRITE(IUNS,'("Model without Random Effects")')
            WRITE(IUNS,'("----------------------------")')
            s_cycle  = 0
            NS_cycle = 0 ! number of additional var cov parameters due to random SCALE
            RR_cycle = 0
            myqdim = 0
            mytotalq = 1
            mypoints0 = nopoints
            myweights0 = noweights
        else if (cycles==2) then
            WRITE(IUN,'("------------------------------")')
            WRITE(IUN,'("Model without Scale Parameters")')
            WRITE(IUN,'("------------------------------")')
            WRITE(*,*) "------------------------------"
            WRITE(*,*) "Model without Scale Parameters"
            WRITE(*,*) "------------------------------"
            WRITE(IUNS,*) "------------------------------"
            WRITE(IUNS,*) "Model without Scale Parameters"
            WRITE(IUNS,*) "------------------------------"
            s_cycle  = 0
            NS_cycle = 0 ! number of additional var cov parameters due to random SCALE
            myqdim = numloc
            RR_cycle = RR
            mytotalq = qtotalR0
            mypoints0 = pointsR0
            myweights0 = weightsR0
        else if (CYCLES==3) THEN
            WRITE(IUN,'("---------------------------")')
            WRITE(IUN,'("Model WITH Scale Parameters")')
            WRITE(IUN,'("---------------------------")')
            WRITE(IUNS,'("---------------------------")')
            WRITE(IUNS,'("Model WITH Scale Parameters")')
            WRITE(IUNS,'("---------------------------")')
            WRITE(*,*) "---------------------------"
            WRITE(*,*) "Model WITH Scale Parameters"
            WRITE(*,*) "---------------------------"
            s_cycle = S
            NS_cycle = 0 ! number of additional var cov parameters due to random SCALE
            myqdim = numloc
            RR_cycle = RR
            mytotalq = qtotalR0
            mypoints0 = pointsR0
            myweights0 = weightsR0
            tau(:) = 0
        else if (CYCLES==4) THEN
            WRITE(IUN,'("-----------------------")')
            WRITE(IUN,'("Model WITH RANDOM Scale")')
            WRITE(IUN,'("-----------------------")')
            WRITE(IUNS,'("-----------------------")')
            WRITE(IUNS,'("Model WITH RANDOM Scale")')
            WRITE(IUNS,'("-----------------------")')
            WRITE(*,*) "-----------------------"
            WRITE(*,*) "Model WITH RANDOM Scale"
            WRITE(*,*) "-----------------------"
            s_cycle = S
            ! NS = number of additional var cov parameters due to random SCALE
            NS_cycle = 1
            myqdim = numloc+1
            RR_cycle = RR
            mytotalq = qtotalR1
            mypoints0 = pointsR1
            myweights0 = weightsR1
            spar(1) = .1
        else if (CYCLES==5) THEN
!            WRITE(IUN,'("------------------------------------------------------")')
!            WRITE(IUN,'("Model WITH RANDOM Scale and Location-Scale Association")')
!            WRITE(IUN,'("------------------------------------------------------")')
!            WRITE(IUNS,'("------------------------------------------------------")')
!            WRITE(IUNS,'("Model WITH RANDOM Scale and Location-Scale Association")')
!            WRITE(IUNS,'("------------------------------------------------------")')
!            WRITE(*,*) "------------------------------------------------------"
!            WRITE(*,*) "Model WITH RANDOM Scale and Location-Scale Association"
!            WRITE(*,*) "------------------------------------------------------"
            s_cycle = S
           ! NS = number of additional var cov parameters due to random SCALE
            RR_cycle = RR
            NS_cycle = ns
            myqdim = numloc+1
            mytotalq = qtotalR1
            spar(ns) = spar(1)
            spar(1:(ns-1)) = 0
            mypoints0 = pointsR1
            myweights0 = weightsR1
        end if
        NPAR_cycle = P+RR_cycle*(mls+nalpha)+S_cycle+NS_cycle+ngamma ! number of parameters
!    write(IUNS,*) mls, numloc, RR, nalpha, RR_cycle, ns_cycle, NPAR_cycle
!    write(*,*) mls, numloc, RR, nalpha, RR_cycle, ns_cycle, NPAR_cycle
!
! start iterations
!
     ! ifin = 1 for iterations using NR (or BHHH depending on the inversion of DER2)
     !      = 2 for final iteration
     !
        ridge  = ridgein
        RIDGEIT= 0
    
         ! set the ridge up for the models with scale parameters
        if (CYCLES>1) then
            ridge = ridgein+.05D0
        END IF
        IF (CYCLES>=4) THEN
            ridge = RIDGEIN+.15D0
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
            minder2 = myder2(1,1)
            do k=2, npar_cycle
                if(minder2 > myder2(k,k)) minder2 = myder2(k,k)
            end do
        
             ! put the ridge back to its initial value after the first 5 & 10 iterations
            IF (ridge > ridgein .and. mod(iter, 10)==0) THEN
                ridge = ridge - .05
            END IF
            if(ridge < 0) ridge = 0 !To get rid of -0 values
            if(iter >= 10 .and. (maxval(abs(myder(1:npar_cycle))) < 2.5)) then! .or. maxval(abs(corec(1:npar_cycle)))<conv*10)) then
                ridge = 0
            else 
                if (cycles <= 3 .and. iter >= 5+cycles*5 .and. maxval(abs(myder(1:npar_cycle))) < 10 .and. minder2 > 100) then
                    ridge = 0
                else
                    if(ridge < ridgein) ridge = ridgein
                end if
            end if
        !
        ! calculate the derivatives and information matrix
        !
            LOGL= 0.0D0
            myDER(1:npar_cycle)=0.0D0
            myder2(1:npar_cycle,1:npar_cycle) = 0.0D0
            mynob = 0

            ILOOP:DO I=1,NClust  ! go over level-2 clusters
                theta1 = 0
                thetav = 0
                mypoints = mypoints0
                myweights = myweights0
                    if(myqdim>0 .and. aquad .ne. 0 .and. (iter >= 10 .or. cycles .eq. 3 .or. cycles .eq. 5)) then
                         do k=1,myqdim
                             sdev = sqrt(thetavs(i,k,k))
                             do q=1,mytotalq
                                 tempC = mypoints(q,k)
                                 mypoints(q,k) = thetas(i,k) + sdev*mypoints(q,k)
                                 call get_phi_ratio(mypoints(q,k),tempC, phiRatio)
                                 myweights(q) = sdev*phiRatio*myweights(q)
                             end do
                         end do
                    end if
                iprob = 0
                myider(1:npar_cycle) = 0
                myider2(1:npar_cycle,1:npar_cycle) = 0
                QLOOP:DO Q=1, mytotalq  ! go over quadrature points
                    qlogprob = 0
                    myqder1(1:npar_cycle,1) = 0
                    jloop: do j=1,ni(i)
                        nob = mynob + j
                        call lambdafun(nob, q, lambda_upper, lambda_lower,myscale)
                        cdf_upper = PhiFn(lambda_upper,nfn)
                        cdf_lower = PhiFn(lambda_lower,nfn)
                        pdf_upper = PhiY (lambda_upper,nfn)
                        pdf_lower = PhiY (lambda_lower,nfn)
                        ddf_upper = PhiDer(lambda_upper, nfn)
                        ddf_lower = PhiDer(lambda_lower, nfn)
                        jprob = cdf_upper - cdf_lower
                        if(jprob < 1e-32) jprob = 1e-32
                        myqlogprob(j) = log(jprob)
                        qlogprob = qlogprob + myqlogprob(j)

                        ! GET FIRST DERIVATIVES
                        dlambda_upper1(1:npar_cycle,1) = 0
                        dlambda_upper1(1:P,1) = -X(NOB,1:p)*myscale                    ! beta
                        t=0
                        do k=1, min(R,RR_cycle*mls) !to account for the case when RR_cycle=0
                            do m=1, k
                                t = t + 1
                                dlambda_upper1(p+ngamma+t,1) = -U(nob,k) * mypoints(q, m) * myscale !cholesky
                            end do
                        end do
                        if((1-mls)*RR_cycle > 0) then
                            sigmaalpha = exp(.5*dot_product(U(NOB,1:nalpha),alpha(1:nalpha)))
                            dlambda_upper1(p+ngamma+1:p+ngamma+nalpha,1) = -mypoints(q,1) * myscale * sigmaalpha*.5*U(nob,1:nalpha)
                        end if
                        dlambda_lower1(1:npar_cycle,1) = dlambda_upper1(1:npar_cycle,1)
                        t = ngamma+P+Rr_cycle*(mls+nalpha)
                        if(s_cycle > 0) then
                            dlambda_upper1(t+1:t+s_cycle,1) = -w(nob,1:s_cycle)*lambda_upper !tau
                            dlambda_lower1(t+1:t+s_cycle,1) = -w(nob,1:s_cycle)*lambda_lower !tau
                        end if

                        if(ns_cycle > 0) then
                            if(ns_cycle > 1) then
                                if(mls .eq. 1) then
                                    dlambda_upper1(npar_cycle-ns_cycle+1:npar_cycle,1) = -mypoints(q,1:ns_cycle)*lambda_upper
                                    dlambda_lower1(npar_cycle-ns_cycle+1:npar_cycle,1) = -mypoints(q,1:ns_cycle)*lambda_lower
                                else
                                    dlambda_upper1(npar_cycle,1) = -mypoints(q,numloc+1)*lambda_upper
                                    dlambda_lower1(npar_cycle,1) = -mypoints(q,numloc+1)*lambda_lower
                                    dlambda_upper1(npar_cycle-ns_cycle+1,1) = -mypoints(q,1)*lambda_upper
                                    dlambda_lower1(npar_cycle-ns_cycle+1,1) = -mypoints(q,1)*lambda_lower
                                    if(ncov > 1) then
                                        dlambda_upper1(npar_cycle-1,1) = -mypoints(q,1)**2*lambda_upper
                                        dlambda_lower1(npar_cycle-1,1) = -mypoints(q,1)**2*lambda_lower
                                    end if
                                end if                                
                            else                                
                                dlambda_upper1(npar_cycle,1) = -mypoints(q,numloc+1)*lambda_upper
                                dlambda_lower1(npar_cycle,1) = -mypoints(q,numloc+1)*lambda_lower
                            end if
                        end if

                        if(yasint(nob) <= ngamma) dlambda_upper1(p+yasint(nob),1) = myscale !threshold
                        if(yasint(nob) >= 2) dlambda_lower1(p+yasint(nob)-1,1) = myscale !threshold
                        myqderp1(j,1:npar_cycle,1) = pdf_upper/jprob*dlambda_upper1(1:npar_cycle,1) - &
                                        pdf_lower/jprob*dlambda_lower1(1:npar_cycle,1)
                        myqder1(1:npar_cycle,1) = myqder1(1:npar_cycle,1) + myqderp1(j,1:npar_cycle,1)
                        myqderp2(j,1:npar_cycle,1:npar_cycle) = ddf_upper*matmul(dlambda_upper1(1:npar_cycle,:),&
                                    transpose(dlambda_upper1(1:npar_cycle,:))) - &
                                    ddf_lower*matmul(dlambda_lower1(1:npar_cycle,:),&
                                    transpose(dlambda_lower1(1:npar_cycle,:)))
                        if((s_cycle > 0 .or. ns_cycle > 0)) then
                            do h=1,npar_cycle
                                do k=1,s_cycle+ns_cycle
                                    if(h <= t+k) then
                                        if(k<=s_cycle) then
                                            myqderp2(j,t+k,h) = myqderp2(j,t+k,h) - w(nob,k)* &
                                                (pdf_upper*dlambda_upper1(h,1)-pdf_lower*dlambda_lower1(h,1))
                                            myqderp2(j,h,k+t) = myqderp2(j,t+k,h)
                                        else
                                            kk = k - s_cycle
                                            if(kk==ns_cycle) then
                                                myqderp2(j,npar_cycle,h) = myqderp2(j,npar_cycle,h) - mypoints(q,numloc+1)* &
                                                        (pdf_upper*dlambda_upper1(h,1)-pdf_lower*dlambda_lower1(h,1))
                                            else if(mls .eq. 1) then 
                                                myqderp2(j,npar_cycle-ns_cycle+kk,h) = &
                                                    myqderp2(j,npar_cycle-ns_cycle+kk,h) - mypoints(q,kk)*&
                                                    (pdf_upper*dlambda_upper1(h,1)-pdf_lower*dlambda_lower1(h,1))
                                            else if(kk==1) then
                                                myqderp2(j,npar_cycle-ns_cycle+kk,h) = &
                                                    myqderp2(j,npar_cycle-ns_cycle+kk,h) - mypoints(q,1)*&
                                                    (pdf_upper*dlambda_upper1(h,1)-pdf_lower*dlambda_lower1(h,1))
                                            else if(kk==2) then
                                                myqderp2(j,npar_cycle-ns_cycle+kk,h) = &
                                                    myqderp2(j,npar_cycle-ns_cycle+kk,h) - mypoints(q,1)**2*&
                                                    (pdf_upper*dlambda_upper1(h,1)-pdf_lower*dlambda_lower1(h,1))
                                            end if
                                            myqderp2(j,h,npar_cycle-ns_cycle+kk) = myqderp2(j,npar_cycle-ns_cycle+kk,h)
                                        end if
                                    end if
                                end do
                            end do
                        end if
                        do h=1,rr_cycle*nalpha
                            do k=1,rr_cycle*nalpha
                                myqderp2(j,p+ngamma+k,p+ngamma+h) = myqderp2(j,p+ngamma+k,p+ngamma+h) - 0.5 * u(nob,k)* &
                                         (pdf_upper-pdf_lower)*dlambda_lower1(p+ngamma+h,1)
                                myqderp2(j,p+ngamma+h,ngamma+p+k) = myqderp2(j,p+ngamma+k,p+ngamma+h)
                            end do
                        end do
                    end do jloop
                    qprob = exp(qlogprob)*myweights(q)
                    qprobs(q) = qprob
                    myqder2(1:npar_cycle,1:npar_cycle) = 0
                    iprob = iprob + qprob
                    myider(1:npar_cycle) = myider(1:npar_cycle) + myqder1(1:npar_cycle,1)*qprob
                    do j=1,ni(i)
                        myqder2(1:npar_cycle,1:npar_cycle) = myqder2(1:npar_cycle,1:npar_cycle) &
                                                + myqderp2(j,1:npar_cycle,1:npar_cycle)*exp(-myqlogprob(j)) &
                                                - matmul(myqderp1(j,1:npar_cycle,:),transpose(myqderp1(j,1:npar_cycle,:)))
                    end do
                    myqder2(1:npar_cycle,1:npar_cycle) = myqder2(1:npar_cycle,1:npar_cycle) + &
                                            matmul(myqder1(1:npar_cycle,:), transpose(myqder1(1:npar_cycle,:)))
                    myider2(1:npar_cycle,1:npar_cycle) = myider2(1:npar_cycle,1:npar_cycle) + &
                                                    myqder2(1:npar_cycle,1:npar_cycle)*qprob
                    if(myqdim > 0) theta1(1:myqdim,1) = theta1(1:myqdim,1) + qprob*mypoints(q,1:myqdim)
                END DO qLOOP
                myider1(1:npar_cycle,1) = myider(1:npar_cycle)/iprob
                myder(1:npar_cycle) = myder(1:npar_cycle) + myider1(1:npar_cycle,1)
                logl = logl + log(iprob)
                !negative of actual hessian is used since corec is added to the current estimate
                do k=1,npar_cycle
                    do j=1,npar_cycle
                        tempx(k,j) = myider1(k,1)*myider1(j,1)
                    end do
                end do
!                tempx(1:npar_cycle,1:npar_cycle) = matmul(myider1(1:npar_cycle,:),transpose(myider1(1:npar_cycle,:)))
                myder2(1:npar_cycle,1:npar_cycle) = myder2(1:npar_cycle,1:npar_cycle) + tempx(1:npar_cycle,1:npar_cycle)
                myder2(1:npar_cycle,1:npar_cycle) = myder2(1:npar_cycle,1:npar_cycle) - myider2(1:npar_cycle,1:npar_cycle)/iprob
                if(myqdim > 0) then
                    theta1(1:myqdim,1) = theta1(1:myqdim,1) / iprob
                    do q=1,mytotalq
                        thetav(1:myqdim,1:myqdim) = thetav(1:myqdim,1:myqdim) + qprobs(q) * &
                            matmul(theta1(1:myqdim,1:1)-transpose(mypoints(q:q,1:myqdim)), &
                            transpose(theta1(1:myqdim,1:1))-mypoints(q:q,1:myqdim))
                    end do
                    thetav(1:myqdim,1:myqdim) = thetav(1:myqdim,1:myqdim) / iprob
                    thetas(i,1:myqdim) = theta1(1:myqdim,1)
                    thetavs(i,1:myqdim,1:myqdim) = thetav(1:myqdim,1:myqdim)
                end if
                mynob = mynob + ni(i)
            END DO ILOOP
            LOGDIFF = LOGL-LOGLP
            LOGLP   = LOGL
            ! determine if an NR iteration is bad and increase the ridge
            ! take the ridge off after 10 good iterations
            
            IF (LOGDIFF/LOGLP > .05 .AND. ITER < MAXIT) THEN
                RIDGEIT = 0
                RIDGE = RIDGE + .05D0
                WRITE(IUN,'("==> BAD NR ITERATION ",I5," with NEW ridge = ",F8.4,/)') ITER,RIDGE
                WRITE(IUNS,'("==> BAD NR ITERATION ",I5," with NEW ridge = ",F8.4,/)') ITER,RIDGE
                WRITE(*,'("==> BAD NR ITERATION ",I5," with NEW ridge = ",F8.4,/)') ITER,RIDGE
                corec(1:npar_cycle) = -.5 * corec(1:npar_cycle)
                GO TO 99
            END IF
        
            if(ifin < 2 .and. cycles > 1) then
                 ! ridge adjustment - diagonal elements only
                do k=1,npar_cycle
                    myder2(k,k) = abs(myder2(k,k))*(1 + ridge)
                end do
            end if
        
            temp(1:npar_cycle,1:npar_cycle) = myder2(1:npar_cycle,1:npar_cycle)
    write(IUNS,*)"2nd Derivatives"
                do k=1,npar_cycle
                    write(IUNS,'(25g15.4)') (myder2(k,i), i=1,k)
                end do
            call inverse(temp(1:npar_cycle,1:npar_cycle), temp(1:npar_cycle,1:npar_cycle), npar_cycle)
!    write(IUNS,*)"Information Matrix"
!                do k=1,npar_cycle
!                    write(IUNS,'(25g15.4)') (temp(k,i), i=1,npar_cycle)
!                end do
            corec1(1:npar_cycle,1) = myder(1:npar_cycle)
            corec1(1:npar_cycle,1) = matmul(temp(1:npar_cycle,1:npar_cycle), corec1(1:npar_cycle,1))
            corec(1:npar_cycle) = corec1(1:npar_cycle,1)
            if(iter<=5) corec = corec*.5
        
            write(IUNS,*)"Derivatives"
            write(IUNS,'(25g15.4)') (myder(k), k=1,npar_cycle)
!            write(*,*)"Derivatives"
!            write(*,'(25f9.3)') (myder(k), k=1,npar_cycle)
            write(IUNS,*)"Beta"
            write(IUNS,'(25g15.4)') (beta(k), k=1,p)
!            write(*,*)"Beta"
!            write(*,'(25f12.3)') (beta(k), k=1,p)
            if(RR_cycle*mls > 0) then
                write(IUNS,*)"Chol"
                write(IUNS,'(25g15.4)') (mychol(k),k=1,rr)
!                write(*,*)"Chol"
!                write(*,'(25f12.3)') (mychol(k),k=1,rr)
            end if
            if(RR_cycle*(1-mls) > 0) then
                write(IUNS,*)"Alpha"
                write(IUNS,'(25g15.4)') (alpha(k),k=1,nalpha)
            end if
            if(S_cycle > 0) then
!                write(*,*)"Tau"
!                 write(*,'(25f12.3)') (tau(k),k=1,s_cycle)
                write(IUNS,*)"Tau"
                write(IUNS,'(25g15.4)') (tau(k),k=1,s_cycle)
            end if
            write(IUNS,*)"Gamma"
            write(IUNS,'(25g15.4)') (mygamma(k), k=1,ngamma)
!            write(*,*)"Gamma"
!            write(*,'(25f8.3)') (mygamma(k), k=1,ngamma)
            if(ns_cycle > 0) then
                write(IUNS,*)"Spar"
                write(IUNS,'(25f15.4)') (spar(k),k=1,ns_cycle)
!                write(*,*)"Spar"
!                write(*,'(25f12.3)') (spar(k),k=1,ns_cycle)
            end if
        
            MAXDER=MAXVAL(ABS(myDER(1:npar_cycle)))
            MAXCORR=MAXVAL(ABS(COREC(1:npar_cycle)))
            WRITE(*,*) iter,'  maximum correction and derivative and ridge'
            WRITE(*,'(25g15.4)') maxcorr,MAXDER,ridge
            WRITE(IUNS,*) iter,'  maximum correction and derivative and ridge'
            WRITE(IUNS,'(25g15.4)') MAXCORR,MAXDER,ridge
        
             ! done with NR and onto last iteration
             IF (IFIN==1 .AND. (MAXCORR <= CONV .OR. ITER >= MAXIT)) THEN
                 IFIN=2
                 ORIDGE=RIDGE
                 RIDGE=0.0D0
             END IF
        do k=1,npar_cycle
            if(corec(k) > 1) corec(k) = 1
            if(corec(k) < -1) corec(k) = -1
        end do
        if(iter<=10) corec(1:npar_cycle) = corec(1:npar_cycle)/2
             ! UPDATE PARAMETERS
            IF (LOGDIFF/LOGLP <= .000001 .AND. RIDGEIT < 10) THEN
                RIDGEIT = RIDGEIT+1
            ELSE IF (LOGDIFF/LOGLP <= .000001 .AND. RIDGEIT >= 10 .and. ifin==1) then
                ridge = ridge - .05D0
                if ( ridge < RIDGEIN) ridge = RIDGEIN
            END IF
            do k=1,npar_cycle
                if(corec(k)*myder(k) < 0 .and. (k<= p .or. k>p+ngamma)) corec(k) = sign(corec(k),myder(k))
            end do

         99  WRITE(*,'("   -2 Log-Likelihood = ",F14.5)') -2*LOGL
             WRITE(IUNS,'("   -2 Log-Likelihood = ",F14.5)') -2*LOGL
!             WRITE(IUN,'("   -2 Log-Likelihood = ",F14.5)') -2*LOGL
            write(IUNS,*)"Corrections"
            write(IUNS,'(25f15.4)') (corec(k), k=1,npar_cycle)
!            write(*,*)"Corrections"
!            write(*,'(25f9.3)') (corec(k), k=1,npar_cycle)

            beta = beta + corec(1:p)
            mygamma(1) = mygamma(1) + corec(p+1)
            do k=2, ngamma
                mygamma(k) = mygamma(k) + corec(p+k)
                if(mygamma(k) < mygamma(k-1)) then 
                    write(*,*) k-1,k,mygamma(k-1), mygamma(k), "threshold value too small!"
                    write(iuns,*) k-1,k,mygamma(k-1), mygamma(k), "threshold value too small!"
                    mygamma(k) = mygamma(k-1)+.01
                end if
            end do
            if (s_cycle > 0) TAU   = TAU      + &
                COREC(P+RR_cycle*(mls+nalpha)+ngamma+1:P+RR_cycle*(mls+nalpha)+ngamma+S_cycle)
            if(RR_cycle*mls > 0 .and. iter > 14) then
                mychol = mychol + COREC(P+ngamma+1:P+ngamma+RR_cycle)/2
                kk = 1
                    do k=1, r
                        do LL=1, k
                            if(k==LL .and. mychol(kk) < 0) mychol(kk) = -mychol(kk)
                            kk = kk + 1
                        end do
                    end do
            else if(RR_cycle*(1-mls) > 0) then
                alpha = alpha + COREC(P+ngamma+1:P+ngamma+nalpha)
            end if
            if(ns_cycle > 0 .and. (iter > 6 .or. (iter > 1 .and. cycles >= 5))) then
                do k=1,ns_cycle
                        spar(k) = spar(k) + corec(npar_cycle-ns_cycle+k)/2
!                    end if
                end do
                if(spar(ns_cycle) < 0) spar(ns_cycle) = abs(spar(ns_cycle))
            end if
             ITER = ITER+1
        END DO IFINLOOP

        if(mls .eq. 1) then
            if(cycles .ne. 4) then
                if(rr_cycle*mls > 0 .and. chol .ne. 2) then         
                    cholamt = rr_cycle
                    mycholspar(1:rr_cycle) = mychol
                    sigma(1:rr_cycle) = mychol
                    do j=1,rr+ns
                        myorder(j) = j
                    end do
                    k=1
                    if(cycles >= 4) then 
                        mycholspar(rr+1:rr+ns) = spar(1:ns)
                        sigma(rr+1:rr+ns) = spar(1:ns)
                        if(chol .ne. 1) then
                            cholamt = rr+ns
                            do j=0,r-1
                                jj1 = j*(j-1)/2
                                do i=1,r-j
                                    myorder(j*r-jj1+i) = r*j-jj1+j+i
                                end do
                                myorder(rr+j+1)= r*(j+1)-jj1+1
                            end do
                        end if
                    end if
                    do j=1,cholamt
                        myorder2(j) = myorder(j)+p+ngamma
                        if(myorder(j)>rr) myorder2(j) = myorder2(j)+s
                    end do
!write(67,*) "Cholesky"
!write(67,'(20f6.2)') (mycholspar(j),j=1,cholamt)
                    if(chol .ne. 2) then
!write(67,*) r, numloc, ns_cycle, r+(1-chol)*min(1,ns_cycle),cholamt
                        call getSStar(mycholspar,r+(1-chol)*min(1,ns_cycle),cholamt,sstar(1:cholamt,1:cholamt))
!write(67,*) "Matrix to transform cholesky back"
!do i=1,cholamt
!    write(67,'(20f6.2)') (sstar(i,j),j=1,cholamt)
!end do
                        do i=1,cholamt
                            sigma(i) = dot_product(mycholspar(1:cholamt),sstar(1:cholamt,i))
                        end do
!write(67,*) "Sigma"
!write(67,'(20f6.2)') (sigma(j),j=1,cholamt)
                        asstar2 = 0
                        ii=1
                        do i=1,npar_cycle
                            asstar2(i,i) = 1
                            jj=1
                            if((i > p .and. i <= p+rr) .or. (i>=npar_cycle+rr-cholamt+1)) then
                                do j=1,npar_cycle
                                    if((j > p .and. j <= p+rr) .or. (j>=npar_cycle+rr-cholamt+1)) then
                                        asstar2(i,j) = 2*sstar(ii,jj)
                                        jj=jj+1
                                    end if
                                end do
                                ii=ii+1
                            end if
                        end do
!write(67,*) "Matrix to transform var/covar"
!do i=1,npar_cycle
!    write(67,'(20f6.2)') (asstar2(i,j),j=1,npar_cycle)
!end do
!write(67,*) "Original var/covar values"
!do i=1,npar_cycle
!    write(67,'(20F15.8)') (temp(i,j), j=1,npar_cycle)
!end do
                        work(1:npar_cycle,1:npar_cycle) = matmul(asstar2(1:npar_cycle,1:npar_cycle), &
                            temp(1:npar_cycle,1:npar_cycle))        
                        adjVar(1:npar_cycle,1:npar_cycle) = matmul(work(1:npar_cycle,1:npar_cycle),&
                                                        transpose(asstar2(1:npar_cycle,1:npar_cycle)))
                    else
                        sigma(1:rr) = mycholspar(1:rr)
                        adjvar = temp
                    end if
                else
                    adjvar = temp
                end if
            end if
        else
            adjvar = temp
        end if
        DO k=1,npar_cycle
            se(k) = dsqrt(dabs(adjvar(k,k)))
        END DO
                
        do k=1,npar_cycle
                write(19,'(25g15.4)') (adjvar(k,i), i=1,k)
        end do
            write(19,*)
             WRITE(13,'(35F15.8)')(BETA(k),k=1,P)
             if(mls .eq. 1) WRITE(13,'(35F15.8)')(sigma(myorder(k)),k=1,RR_cycle*mls)
             if(mls .eq. 0) WRITE(13,'(35F15.8)')(alpha(k),k=1,nalpha)
             WRITE(13,'(35F15.8)')(TAU(k),k=1,S_cycle)
             if(mls .eq. 1) write(13,'(35F15.8)')(sigma(myorder(rr+k)),k=1,ns_cycle)
             if(mls .eq. 0) write(13,'(35F15.8)')(spar(k),k=1,ns_cycle)
!             WRITE(13,'(35F15.8)')(SE(k),k=1,p),(se(myorder2(k)),k=1,rr_cycle), &
!             (se(k),k=p+rr+1,p+rr+s),(se(myorder2(k)),k=rr+1,rr+ns)

        !if(cycles .ne. 4 .or. ncov .eq. 0) then
        
             ! WRITE RESULTS
            WRITE(iun,562)ITER-1,ORIDGE,LOGL,LOGL-NPAR, &
            LOGL-0.5D0*DBLE(NPAR)*DLOG(DBLE(NClust)),0-2*LOGL,0-2*(LOGL-NPAR), &
            0-2*(LOGL-0.5D0*DBLE(NPAR)*DLOG(DBLE(NClust)))
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
         57 FORMAT(/,'Variable',20x,'    Estimate',4X,'AsymStdError',4x, &
                  '     z-value',4X,'     p-value',/,'----------------',12x,  &
                  '------------',4X,'------------',4X,'------------',4X,'------------')
        
            PVAL=0.0D0
            ZVAL=0.0D0
            WRITE(IUN,'("BETA (regression coefficients)")')
            DO L=1,P
                ZVAL = BETA(L)/SE(L)
                PVAL = 2.0D0 *(1.0D0 - PHIFN(DABS(ZVAL),0))
                WRITE(IUN,804)BLABel(L),BETA(L),SE(L),ZVAL,PVAL
             END DO
            write(iun,*)
            if(rr_cycle*mls > 0) then
                if(chol .ne. 1) then
                    WRITE(IUN,'("Random (location) Effect Variances and Covariances")')
                    k = 1
                    do i=1,r
                        do j=1,i
                            ZVAL = sigma(k)/SE(k+p+ngamma)
                            PVAL = 2.0D0 *(1.0D0 - PHIFN(DABS(ZVAL),0))
                            if(i==j) then
                                WRITE(IUN,804)clabel(i),sigma(k),SE(k+p+ngamma),ZVAL,PVAL
                            else
                                WRITE(IUN,805)'Covariance',j,i,sigma(k),SE(k+p+ngamma),ZVAL,PVAL
                            end if
                            k = k+1
                        end do
                    END DO
                else
                    WRITE(IUN,'("Random (location) Effect Cholesky Matrix")')
                    k = 1
                    do i=1,r
                        do j=1,i
                            ZVAL = mychol(k)/SE(k+p+ngamma)
                            PVAL = 2.0D0 *(1.0D0 - PHIFN(DABS(ZVAL),0))
                            if(i==j) then
                                WRITE(IUN,804)clabel(i),mychol(k),SE(k+p+ngamma),ZVAL,PVAL
                            else
                                WRITE(IUN,805)'Covariance',j,i,mychol(k),SE(k+p+ngamma),ZVAL,PVAL
                            end if
                            k = k+1
                        end do
                    END DO
                end if
            else if(rr_cycle*(1-mls) > 0) then
                WRITE(IUN,'("ALPHA (BS variance parameters: log-linear model)")')
                DO L=1,nalpha
                    L2 = P+ngamma+L
                    ZVAL = alpha(L)/SE(L2)
                    PVAL = 2.0D0 *(1.0D0 - PHIFN(DABS(ZVAL),0))
                    WRITE(IUN,804)cLABel(L),alpha(L),SE(L2),ZVAL,PVAL
                END DO
            end if    
            write(iun,*)
        804 FORMAT(A24,4(4x,F12.5))
        805 FORMAT(A10,I0,I0,12X,4(4x,F12.5))
            if(s_cycle > 0) then
                WRITE(IUN,'("TAU (WS variance parameters: log-linear model)")')
                DO L=1,S_cycle
                    L2 = P+RR_cycle*(mls+nalpha)+ngamma+L
                    ZVAL = TAU(L)/SE(L2)
                    PVAL = 2.0D0 *(1.0D0 - PHIFN(DABS(ZVAL),0))
                    WRITE(IUN,804)TLABel(L),TAU(L),SE(L2),ZVAL,PVAL
                END DO
                write(iun,*)
            end if
            IF (ngamma .GT. 0) THEN
                WRITE(IUN,584)
                584 FORMAT('Thresholds (for identification)')
                DO k = 1,ngamma
                    k2   = p+k 
                    ZVAL = mygamma(k)/SE(k2)
                    PVAL = 2.0D0 *(1.0D0 - PHIFN(DABS(ZVAL),0))
                    WRITE(IUN,'(I2,22X,4(4x,F12.5))') k,mygamma(k),SE(k2),ZVAL,PVAL
                END DO
            ENDIF

            if(cycles == 5) then
                write(iun,*)
                if(mls .eq. 1) then
                    WRITE(IUN,'("Random location effects on WS variance (log-linear model)")')
                    do k=1,ns_cycle-1
                       L2=P+RR_cycle+S_cycle+k+ngamma
                        ZVAL = SPAR(k)/SE(L2)
                        PVAL = 2.0D0 *(1.0D0 - PHIFN(DABS(ZVAL),0))
                        WRITE(IUN,804)cLABel(k),spar(k),SE(L2),ZVAL,PVAL
                    end do
                else
                    WRITE(IUN,'("Random location effects on WS variance (log-linear model)")')
                    do k=1,ns_cycle-1
                       L2=P+nalpha+S_cycle+k+ngamma
                        ZVAL = SPAR(k)/SE(L2)
                        PVAL = 2.0D0 *(1.0D0 - PHIFN(DABS(ZVAL),0))
                        if(k .eq. 1) WRITE(IUN,804)"Linear",spar(k),SE(L2),ZVAL,PVAL
                        if(k .eq. 2) WRITE(IUN,804)"Quadratic",spar(k),SE(L2),ZVAL,PVAL
                    end do
                end if
            end if
        
            IF (cycles>=4) THEN
                 write(iun,*)
                write(IUNS,'(25g15.4)') (se(k), k=1,npar_cycle)
                WRITE(IUN,'("Random scale standard deviation")')
                L2=P+RR_cycle*(mls+nalpha)+S_cycle+ns_cycle+ngamma
                if(spar(ns_cycle) < 0) spar(ns_cycle)=abs(spar(ns_cycle))
                ZVAL = SPAR(ns_cycle)/SE(L2)
                PVAL = 2.0D0 *(1.0D0 - PHIFN(DABS(ZVAL),0))
                WRITE(IUN,804)'Std Dev         ',SPAR(ns_cycle),SE(L2),ZVAL,PVAL
            end if
        
        
            IF (MAXDER > CONV .AND. ITER >= MAXIT) THEN
                WRITE(IUN,'("NOTE: CONVERGENCE CRITERION WAS NOT ACHIEVED")')
                WRITE(IUN,'("FINAL FIRST DERIVATIVE AND CORRECTION VALUES")')
                DO L=1,NPAR_cycle
                  WRITE(IUN,*)myDER(L),COREC(L)
                 END DO
         !   END IF
            end if
        close(IUN)
    END DO CYCLELOOP

    close(19)
    close(13)
    open(19,file="MixorS_both.var")
    open(13,file="MixorS_both.est")
     
                write(19,'(25g15.4)') adjvar(p+1,p+1),(adjvar(p+1,i), i=1,npar_cycle)
            do k=1,npar_cycle
                write(19,'(25g15.4)') adjvar(k,p+1),(adjvar(k,i), i=1,npar_cycle)
            end do

             WRITE(13,'(a24, F15.8)')"intercept",-1*mygamma(1)
             WRITE(13,'(a24, F15.8)')(blabel(k),BETA(k),k=1,P)
             WRITE(13,'(a24, F15.8)')("Threshold",mygamma(k),k=1,ngamma)             
             if(mls .eq. 1) WRITE(13,'(a24,F15.8)')("VarCov",sigma(myorder(k)),k=1,RR_cycle)
             if(mls .eq. 0) WRITE(13,'(a24,F15.8)')("Alpha",alpha(k),k=1,RR_cycle*nalpha)
             WRITE(13,'(a24,F15.8)')(tlabel(k),TAU(k),k=1,S_cycle)
             if(mls .eq. 1) write(13,'(a24,F15.8)')("VarCov",sigma(myorder(rr+k)),k=1,ns_cycle)
             if(mls .eq. 0) write(13,'(a24,F15.8)')("S",spar(k),k=1,ns_cycle)
    open(23,file="mixors.lik")
    write(23,*) logl
    close(23)

        OPEN(UNIT=IUN,FILE="MIXorS_both.OUT",access="append")
     write(iun,*)
808 FORMAT(A24,3(4x,A12))
     WRITE(IUN,'("WS variance ratios and 95% CIs")')
     write(iun,'("------------------------------")')
     write(iun,*)
    WRITE(IUN,808) 'Variable        ','Ratio','Lower','Upper'
    write(iun,808)'---------------------','------------------','------------','------------'
         WRITE(IUN,'("TAU (WS variance parameters: log-linear model)")')
         myz = 1.959964
     DO L=1,S_cycle
        L2 = P+ngamma+RR_cycle+L
        tauhat = exp(tau(l))
        tauhatlow = exp(tau(l)-myz*se(l2))
        tauhatup = exp(tau(l)+myz*se(l2))
        WRITE(IUN,804)TLABel(L),tauhat, tauhatlow, tauhatup
     END DO
         if(mls .eq. 0) WRITE(IUN,'("ALPHA (BS variance parameters: log-linear model)")')
     do k=1,(1-mls)*rr_cycle
        L2 = P+ngamma+k
        tauhat = exp(alpha(k))
        tauhatlow = exp(alpha(k)-myz*se(l2))
        tauhatup = exp(alpha(k)+myz*se(l2))
        WRITE(IUN,804)cLABel(L),tauhat, tauhatlow, tauhatup
     END DO
    write(iun,*)
    if(ns_cycle > 0) then
        WRITE(IUN,'("Random (location) Effect Variances")')
        do k=1,ns_cycle-1
            L2=P+RR_cycle+S_cycle+k+ngamma
            tauhat = exp(spar(k))
            tauhatlow = exp(spar(k)-myz*se(l2))
            tauhatup = exp(spar(k)+myz*se(l2))
            WRITE(IUN,804)cLABel(k),tauhat, tauhatlow, tauhatup
        end do

        write(iun,*)
        WRITE(IUN,'("Random scale standard deviation")')
        L2=P+RR_cycle*(mls+nalpha)+S_cycle+ns_cycle+ngamma
        tauhat = exp(spar(ns_cycle))
        tauhatlow = exp(spar(ns_cycle)-myz*se(l2))
        tauhatup = exp(spar(ns_cycle)+myz*se(l2))
        WRITE(IUN,804)'Std Dev         ',tauhat, tauhatlow, tauhatup
    end if   
    
    if(aquad .ne. 0) then
        do k=1,nclust
            write(29,'(i8,25g23.4)') ids(k),(thetas(k,i), i=1,ndim), ((thetavs(k,i,j),j=1,i),i=1,ndim)
        end do
    end if
     CLOSE(IUN)
     CLOSE(IUNS)
     CLOSE(1)
     CLOSE(2)
     CLOSE(3)
     close(19)
    close(29)
    close(13)
    DEALLOCATE(myDER2,dlambda_upper1,dlambda_lower1,myider1,corec1,myqder1,myqderp1,se,theta1,thetas,thetav,&
            qprobs,thetavs, sstar, work, adjvar, asstar2, sigma, mycholspar, pointsR0, &
            weightsR0, pointsR1, weightsR1)
END SUBROUTINE MIXorEST

! ----------------------------------------------------------------------------
!@    lambdafun - return lambda(m) and lambda(m-1), the linear predictors

subroutine lambdafun(nob,q,lambda_upper,lambda_lower,myscale)

    use mixor_globals
    implicit none

      ! subroutine arguments:
    integer, intent(IN) :: nob           ! index of current obs
    integer, intent(IN) :: q             ! index of quadrature node
    real(kind=8), intent(OUT)   :: lambda_upper  ! lambda ('z') value for upper
                                                !   category of ds_y%data(nob,1), say m
    real(kind=8), intent(OUT)   :: lambda_lower  ! lambda ('z') value for lower
                                                !   category below ds_y%data(nob,1), m-1
    real(kind=8), intent(OUT)   :: myscale       ! exp(wij'tau)
      
      ! Local variables:
    integer :: k,m,t
    real(kind=8) :: Xb, wt, Us


!     Location linear predictor
    Xb = 0  ! calculate X BETA for the current subject
    do k = 1, p
        Xb=Xb+beta(k)*X(nob,k)
    end do
!     Scale linear predictor
    wt = 0  ! calculate Z GAMMA  for the current subject
    do k = 1, S_cycle
        wt = wt + tau(k)*w(nob,k)
    end do
    if (ns_cycle > 0) then
        wt = wt + mypoints(q, numloc+1)*spar(ns_cycle)
        do k=1, ns_cycle-1
            if(mls .eq. 1 .or. k .eq. 1) then
                wt = wt + mypoints(q,k)*spar(k)
            else
                wt = wt + mypoints(q,k)**2*spar(k)
            end if
        end do
    end if
!     Random effect term (linear predictor?)
    Us = 0  ! calculate U SIGMA for the current subject
    if(rr_cycle*mls > 0) then
        t=0
        do k=1, R
            do m=1, k
                t = t + 1
                Us = Us + U(nob,k) * mychol(t) * mypoints(q, m)
            end do
        end do
    else if(rr_cycle*(1-mls) > 0) then
        Us = Us + mypoints(q,1) * exp(.5*dot_product(alpha(1:nalpha),U(nob,1:nalpha)))
    end if

    lambda_lower = -1.0d+10  ! -10 billion
    lambda_upper = +1.0d+10  ! +10 billion
    m = Yasint(nob)
    myscale = exp(-wt)

    if (m <= ngamma) then  
        lambda_upper = (mygamma(m) -Xb -Us)*myscale
    end if
    if (m > 1) then  
        lambda_lower = (mygamma(m-1) -Xb -Us)*myscale
    end if
end subroutine lambdafun

function isequal_dp(x,y)
    implicit none
    logical                    :: isequal_dp
    real(kind=8), intent(IN)      :: x,y
    real(kind=8),parameter        :: delta = .00001d0

    isequal_dp = .FALSE.
    if(x .gt. (y-delta) .and. x .lt. (y+delta)) isequal_dp = .TRUE.
end function isequal_dp
