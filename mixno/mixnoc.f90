!  **************************************************************
!  MIXNO        
!
!  modified on 
!  10/11/02 to correct AIC & BIC statistics
!           and to put in estimation of multiple random effects
!
!  modified on 8/5/04 to correct BIC statistics - use N not NTOT
!
!
! *************************************************************
!  nominal RRM MAIN PROGRAM                               
!                                                         
!  Model                                                  
!  Y = X BETA + W ALPHA + ERROR                           
!                                                         
!  Indices                                                
!     N  = TOTAL NUMBER OF SUBECTS         I = 1 .. N     
!     R  = DIMENSION OF RANDOM BETA SPACE  H = 1 .. R     
!     P  = DIMENSION OF FIXED ALPHA SPACE  L = 1 .. P     
!    JJ  = NUMBER OF CATEGORIES            J = 1 .. MAXJ  
!  NI(I) = NUMBER OF OBSERVED TIMEPTS      K = 1 .. NI(I) 
!          WITH MAXNI = MAXIMUM OF NI(I)                  
!  RNII  = SIZE OF EACH SUBJ' X VECTOR    HR = 1 .. RNII  
!  PNII  = SIZE OF EACH SUBJ' W VECTOR    LP = 1 .. PNII  
!                                                         
!  Parameters                                             
!  Y   = MATRIX (N, MAXNI) DEPENDENT VARIABLE             
!  X   = MATRIX (N, RNI) DESIGN FOR RANDOM EFFECTS        
!  W   = MATRIX (N, PNI) FIXED COVARIATES                 
! ALPHA= VECTOR (P) OF FIXED EFFECTS                      
!  MU  = VECTOR (R) OF MEAN BETA VALUES                   
! SIGMA= VECTOR (RR) OF BETA VARIANCE COVARIANCE          
!        WHERE HR = 1 .. RR    RR = (R * (R+1)) / 2       
!                                                         
! CONV = CONVERGENCE CRITERION                            
!  NPR = NUMBER OF SUBJECTS TO PRINT DATA OF              
!                                                         
! *************************************************************
! this program was modified on 1/19/95 
! a) add a ridge = .05 to the matrix of 2nd derivatives at all steps
!    until convergence is achieved, unless the log-likelihood doesn't
!    increase in which case change the ridge to 1.0
! b) comment out the unneccessary statements in the 20-30 lines below
!    (no need to have region.blg since MAT and CMAT are hard-coded)
! c) include code to calculate the % of non-varying response vectors
!    (see subroutine YSAME)
! d) use a modified INVS & CHSKY subroutines to detect when computational
!    difficulties are occuring and then to reduce the # of random effects
!    ==> if a variance term gets too small (ISIG=1)
!    ==> when INVS or CHSKY detects a singularity problem
!    ==> look at the flag variable IER    (=1 means trouble)
!    ==> look at the flag variable NONPOS (=0 means trouble)
! e) modify the PRNT subroutine to work with LAHEY
!    (declare TITLE as REAL*4)
! f) estimate the ln of the diagonal elements of ISIGMA (the cholesky
!    of the random-effects variance covariance matrix) if IFIN=0
!
! this program was modified on 7/25/95 
! a) to allow the user to leave out the specifications for IPRIOR and
!    IUNIF - to agree with the version of the manual sent to 
!    the biomedical computer programs journal 
!
! this program was modified on 9/28/95 
! a) to put in the complementary log-log response function           
!
! this program was modified on 9/27/96
! a) allow the variance covariance matrix of the random effects
!    to be a diagonal matrix (uncorrelated random effects)
!    (IDIAG=1)
! b) allow the mean of the random effects
!    to be zero 
!    (NOMU=1)
!
! this version was updated 9/27/96 to 
! a) write out estimates to MIXNO.EST (with labels)
! b) write out asymptotic variance-covariance matrix (in unpacked form)
!  to MIXNO.VAR (without labels - in the order specified in MIXNO.EST)
!
! this version was updated 10/9/96 to include
!
!    IVSEP  =  0 R random effects don't indicate separate groups of subjects
!          1 R random effects do    indicate separate groups of subjects
!            (or in IRT case, X variables for R random effects change
!             but there is only 1 random effect beta)
!
! this version was updated 5/15/97 to include
!
!    IVCAT  =  0 RR random effects don't vary across MAXJ-1 category contrasts
!          1 RR random effects do    vary across MAXJ-1 category contrasts
!
! this version was updated 8/10/01 to include
!
!    ICNT  = reference cell (=0), repeated (=1), or helmert (=2) contrasts
!            for the outcome variable                                      
!

program RRM_MODULE ! must have same name in all 4 programs
	IMPLICIT NONE
	
	! SHARED variables - these are used to pass data between the 
	! parts of the analysis. 
	
	CHARACTER*4 HEAD(30)
	CHARACTER (LEN=4), ALLOCATABLE :: IXLAB(:)
	CHARACTER (LEN=8), ALLOCATABLE :: IBLABEL(:),IALABEL(:) 
	
	REAL (KIND=8), ALLOCATABLE :: IWORKR(:),ISIGTAU(:),IAQEMP(:), &
	      ICODE(:), IBQ(:), IBQ1(:), IAQ(:), IALPHA0(:), IMU0(:), &
	      IPRJ(:), INUMJ(:), IDER(:), IDER2(:),ISIGMA0(:), ICON(:), &
	      ITHETA(:), ITHETAV(:)
	      
	REAL (KIND=8), ALLOCATABLE, TARGET:: ISIGMA(:),IMU(:),IALPHA(:)
	
	REAL (KIND = 8),POINTER ::ALLDAT(:)
	
	REAL*8    CONV, RIDGE, RLOGLP, RLOGL, AIC, SBC, &
	          DEV, AICD, SBCD, RNTOT, RNPAR
	
	INTEGER ICCY, IFIN, IGO, IRBAD, IRES, ISIG, IT, ITLAST, IUN, &
	      NCON, ND, NFN, NONPOS, NPR, P, R, RR, NQ, NQ1, IUNIF, MJ, &
	      RRJ, IER, IDIAG, IPRIOR, IVCAT, IVSEP, IWT, MAXJ, N, NOMU, &
	      NPAR, NPARR, NRP1,ICNT, IRT, IRTT
	
	INTEGER :: INITIALIZED = 0
	
	INTEGER ,POINTER ::IDNI(:)
	CHARACTER, PARAMETER:: PROGNAME*5 = 'MIXNO'
	
	   CALL INIT()  ! supply command line params to INIT to override
	                       ! filenames for DEF, DAT, etc.
	   DO WHILE(MainLoop() .NE. 0) 
	   ENDDO
	   CALL CLEANUP()
	
	CONTAINS
	
	
	
	! ************************************************
	!               **********************                         
	!               *  SUBROUTINE INIT(FILENAMES)
	!               **********************                         
	!                                                                       
	! Read in the specified DEF file, the indicated DAT file,
	! and set up to process the information.  This includes all 
	! of the one-time allocation of arrays, and initialization
	! of variables.  This is done in a separate routine so that 
	! the DLL caller can split up the execution into smaller
	! pieces so processing can be monitored and possibly cancelled
	! midway.  TG 3/7/01
	!
	! CALL INIT(FILENAMES)
	!                                                                       
	! FILENAMES .. ARRAY OF 4 x 40 character file names, to be used
	!              to open the input and output files.  This array
	!              is loaded either from the command line for the
	!              EXE version, or with a string passed through the
	!              DLL interface.  It is assumed to have been
	!              prepared by MIXLIB/SET_FILENAMES, and to contain
	!              in elements 1-4 the names: DEFFILE,DATFILE,OUTFILE,
	!              BASENAME.  If no values are supplied on the command
	!              line, then the default behaviour matches the historical
	!              behaviour, where the module name is used for the 
	!              DEF, EST,VAR, and RES files, and the DAT and OUT
	!              file names are taken from the DEF file.  TG 3/7/01
	!                                                                       
	! ************************************************
	SUBROUTINE INIT()
	USE MIXLIB
	IMPLICIT NONE
	
	   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	   ! LOCAL VARS ONLY - SHARED VARS ARE DECLARED
	   ! BELOW THE 'MODULE' STATEMENT AT THE TOP OF THE 
	   ! FILE
	   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	   
	   CHARACTER*8  YLabel,XLabel,TLABEL
	   CHARACTER*80 FILEDAT,FILEOUT
	   CHARACTER*160 LABTEMP,TEMPFORM
	   
	   INTEGER MAXCOL, NQR, NQRR, NCPAR, IDIND, &
	         L, J, ICATYX, MAXXJ, ICXW, ICX, NCATYX, &
	         MISS, ISTART, NALL, NTOT, MAXK, NXLAB, &
	         K, MAXJXJ, II, KI, K2, IR, IC, IC2, IC3, ICW, &
	         K3, ICAT, K4, I, IC4, ICA, J2, J3, L2, KIND, RJ,PJ, &
	         H,HR,YIND,P1,WTIND,XCIND,P1P1
	   
	   INTEGER , ALLOCATABLE:: IXIND(:), IWIND(:)
	  
	   REAL*8 YMISS, UNO, CATADD, WSUM, SUM2, &
	         DENUM, DET, TEMPSUM, WA, YDEV
	         
	   REAL (KIND=8), ALLOCATABLE :: IXMISS(:), IWMISS(:),  &
	      IAQ1(:),IAQ2(:),CODEX(:), IWP1P1(:), ICATFQ(:), &
	      ICATFQX(:),IWORKCAT(:),IWORKP1(:),IWRKP1(:)
	     
	   LOGICAL KEEP
	   
	   ! Start by checking if we are ready to run yet.  You 
	   ! can't call init twice in a row, and must call 
	   ! cleanup before you can start another session.
	   
	   IF(INITIALIZED .NE. 0) THEN
	      CALL POST_ERROR("Init() routine called out of sequence in " // PROGNAME)
	      RETURN
	   ENDIF
	   
	   ! You can override one or all of the filenames.  
	   ! If the DATFILE or OUTFILE name are included, 
	   ! they override the name given in the DEFFILE.  
	   ! If an extension is used with the override name, 
	   ! then it will be used, otherwise the default extensions 
	   ! of DEF,DAT, and OUT will be used for the three files.   
	   ! If the BASENAME is included, it overrides the base of 
	   ! the names used for the incidental output files (the 
	   ! extensions do not change).  The override names would 
	   ! be supplied either as arguments to the DLL routine
	   
	   ! READ IN TWO LINES (60 CHARS EACH) FOR THE TITLE 
	   ! AND THEN PARAMETERS AND STARTING VALUES FROM MIXOR.DEF
	
	    open(1, file="mixno.def")
	
	   READ(1,"(15A4)") HEAD
	   READ(1,"(A80)")FILEDAT
	   READ(1,"(A80)")FILEOUT
	   READ(1,"(80X)")   ! Skip over the DEF file name
	
	
	   OPEN(2, FILE= fileout)
	   OPEN(3, FILE= "mixno.est")
	   OPEN(4, FILE= "mixno.var")
	
	   READ(1,*) NPR,MAXCOL,R,P,CONV,MAXJ,MISS,ISTART,IWT, &
	             ICATYX,NQ1,NCON,IDIAG,NOMU,IVSEP,IVCAT,ICNT
	
	   IF (NQ1 .LE. 1) THEN
	      BACKSPACE 1
	      READ(1,*) NPR,MAXCOL,R,P,CONV,MAXJ,MISS,ISTART, &
	             IWT,ICATYX,IPRIOR,IUNIF,NQ1,NCON,IDIAG,NOMU,IVSEP,IVCAT,ICNT
	   ENDIF
	   
	   IPRIOR = 0
	
	   ! IRES = 0 DON'T PRINT OUT INDIVIDUAL PARAMETERS          
	   !      = 1 PRINT INDIVIDUAL PARAMETERS TO FILERES         
	
	   NFN =1
	   IRES=0
	   IF(R.GE.1) IRES=1
	   IF (IRES .EQ. 1) OPEN(5, file="mixno.res")
	
	   ! ALLOCATE  R 
	   ALLOCATE(IXMISS(R))
	   ALLOCATE(IXIND(R))
	   ALLOCATE(IWORKR(R))
	   ALLOCATE(ISIGTAU(R))
	
	   ! IRT = dimension of theta
	
	   IF (IVSEP .EQ. 0) THEN
	      IRT = R
	   ELSEIF (IVSEP .EQ. 1) THEN
	      IRT = 1
	   ENDIF          
	   ALLOCATE(ITHETA(IRT))
	
	   ! ALLOCATE  rr = (r * (r+1) / 2) 
	
	   IF (IDIAG .EQ. 0 .AND. IVSEP .EQ. 0) THEN
	      RR = (R * (R+1)) / 2
	   ELSE
	      RR = R
	   ENDIF
	   RRJ = RR
	   IF (IVCAT .EQ. 1) RRJ = RR*(MAXJ-1)
	   ALLOCATE(ISIGMA(RRJ))
	
	   IRTT = (IRT*(IRT+1))/2
	   ALLOCATE(ITHETAV(IRTT))
	
	   ! ALLOCATE  nquad points per dimension (nq1)
	   ALLOCATE(IBQ1(NQ1))
	
	   ! ALLOCATE  nq1**r
	   IF (IVSEP .EQ. 0) THEN
	      NQR = NQ1**R
	   ELSE
	      NQR = NQ1
	   ENDIF
	   ALLOCATE(IAQ(NQR))
	   ALLOCATE(IAQEMP(NQR))
	
	   ! ALLOCATE  R * (nq1**r)
	   IF (IVSEP .EQ. 0) THEN
	      NQRR = R*NQR
	   ELSE
	      NQRR = NQ1
	   ENDIF
	   ALLOCATE(IBQ(NQRR))
	
	   ! ALLOCATE  these vectors consectutively 
	   !           so that random-effects can 
	   !           be reassigned as fixed-effects
	   !           if numerical problems develop 
	
	   RJ = R*(MAXJ-1)
	   PJ = P*(MAXJ-1)
	   ALLOCATE(IMU(MAX(1,RJ)))   ! Lahey doesn't like 0 length allocs
	   IMU = 0.0D0
	   ALLOCATE(IALPHA(PJ))
	   ALLOCATE(IBLABEL(MAX(1,R)))
	   ALLOCATE(IALABEL(P))
	
	   ! ALLOCATE  p 
	   ALLOCATE(IWMISS(P))
	   ALLOCATE(IWIND(P))
	
	   ! ALLOCATE  maxj
	   ALLOCATE(ICODE(MAXJ))
	   ALLOCATE(IPRJ(MAXJ))
	   ALLOCATE(INUMJ(MAXJ))
	
	   ! ALLOCATE    (npar * (npar+1)) / 2
	   !    where    npar = r + rr + p + ngam
	
	   IF (IPRIOR .EQ. 0 .AND. NOMU .EQ. 0 .AND. IVCAT .EQ. 1) THEN
	      npar = (p+r)*(MAXJ-1)+rr*(MAXJ-1)
	   ELSEIF (IPRIOR .EQ. 0 .AND. NOMU .EQ. 0 .AND. IVCAT .EQ. 0) THEN
	      npar = (p+r)*(MAXJ-1)+rr 
	   ELSEIF ((IPRIOR .EQ. 1 .OR. NOMU .EQ. 1) .AND. IVCAT .EQ. 1) THEN
	      npar = p*(MAXJ-1)+rr*(MAXJ-1)
	   ELSEIF ((IPRIOR .EQ. 1 .OR. NOMU .EQ. 1) .AND. IVCAT .EQ. 0) THEN
	      npar = p*(MAXJ-1)+rr 
	   ENDIF
	
	   RNPAR = DBLE(NPAR)
	   
	   ALLOCATE(IDER(NPAR))
	
	   NCPAR = NPAR*NCON
	   NPARR = (NPAR * (NPAR+1)) / 2
	   ALLOCATE(ICON(NCPAR))
	   ALLOCATE(IDER2(NPARR))
	
	   WRITE(6,*) R,' random terms'
	   WRITE(6,*) P,' fixed  terms'
	
	   ! MAXCOL is the dimension of the data matrix to read in     
	   !
	   ! IDIND (integer) indicates which column contains the level 2 ID
	   !
	   ! YIND (integer) indicates which column contains the dependent var
	   !
	   ! R is the dimension of the X matrix to be used in the analysis
	   !
	   ! IXIND() (integer) indicates which columns are to be used for X 
	   !
	   ! P is the dimension of the W matrix to be used in the analysis
	   !
	   ! IWIND() (integer) indicates which columns are to be used for W 
	   !
	   ! WTIND (integer) indicates which column to be used for the 
	   !       weighting of level-2 units           
	   !
	   ! ISTART = 0 FOR AUTOMATIC STARTING VALUES
	   !          1 FOR STARTING VALUES READ IN
	   !
	   ! IWT    = 0 EACH PERSON (2ND LEVEL UNIT) IS WEIGHTED EQUALLY
	   !          1 FOR DIFFERENTIAL WEIGHTING FOR EACH PERSON (2ND LEVEL)
	   !
	   !
	   ! ICATYX = 0 NO CROSSTAB FOR ANY (X or W) BY Y                        
	   !          1 DO A CROSSTAB FOR A SPECIFIED (X or W) BY Y   
	   !          (the x or w variable is specified below 
	   !           as being in column XCIND)
	   !
	   ! IPRIOR = 0 specified PRIOR FOR RANDOM EFFECTS                 
	   !          1 EMPIRICAL PRIOR FOR RANDOM EFFECTS                     
	   !
	   ! IUNIF  = 0 NORMAL PRIOR FOR RANDOM EFFECTS                 
	   !          1 rectagular PRIOR FOR RANDOM EFFECTS                     
	   !
	   ! NQ1    = number of quadrature nodes per dimension          
	   !
	   ! NFN    = 0 probit                                          
	   !          1 logistic                                                
	   !          2 complementary log-log
	   !          3 log-log
	   !
	   ! NCON   =  number of transforms of estimated parameters 
	   !           (linear re-expressions)                          
	   !
	   ! IDIAG =  0 correlated random effects                       
	   !          1 independent random effects                              
	   !
	   ! NOMU  =  0 estimate the mean of the random effects         
	   !          1 don't estimate the mean of the random effects 
	   !            (assumed to be zero)
	   !
	   ! IVSEP =  0 R random effects don't indicate separate groups of subjects
	   !          1 R random effects do    indicate separate groups of subjects
	   !
	   ! ICNT  =  0 contrasts to 1st category                                   
	   !          1 repeated contrasts (2 vs 1, 3 vs 2, etc)                     
	   !          2 helmert contrasts (e.g., 1 vs 2+3, 2 vs 3)
	   !
	   !  IFIN = 0 ESTIMATE TAU (ln of diagonaol elements of cholesky)
	   !       = 1 ESTIMATE SIGMA (cholesky of variance-covariance
	   !                           matrix of random effects)
	   
	   READ(1,*) IDIND,YIND
	
	   IF (R .GE. 1) READ(1,*) (IXIND(H), H=1,R)
	   IF (P .GE. 1) READ(1,*) (IWIND(L), L=1,P)
	   IF (IWT .EQ. 1) READ(1,*) WTIND
	
	   READ(1,*)(ICODE(J), J = 1,MAXJ)
	   IF (ICATYX .EQ. 1) THEN
	      READ(1,*) XCIND,MAXXJ
	      ALLOCATE(CODEX(MAXXJ))
	      BACKSPACE(1)
	      READ(1,*) XCIND,MAXXJ,(CODEX(J),J=1,MAXXJ)
	      
	      DO H = 1,R
	         IF (XCIND .EQ. IXIND(H)) THEN
	            ICXW = 1
	            ICX  = H
	            NCATYX = H + 1
	         ENDIF
	      END DO
	      
	      DO L = 1,P
	         IF (XCIND .EQ. IWIND(L)) THEN
	            ICXW = 2
	            ICX  = L
	            NCATYX = L + R + 1
	         ELSE
	         ENDIF
	      END DO
	   ENDIF
	
	   IF (MISS .EQ. 1) THEN
	      READ(1,*) YMISS
	      IF (R .GE. 1) READ(1,*) (IXMISS(H), H=1,R)
	      IF (P .GE. 1) READ(1,*) (IWMISS(L), L=1,P)
	   ENDIF
	
	   READ(1,8) YLabel
	   IF (YLABEL .EQ. '        ') READ(1,8) YLABEL
	
	   IF (R .NE. 0) THEN
	      READ(1,8) TLABEL
	      IF (TLABEL .EQ. '        ') THEN
	         READ(1,8) (IBLABEL(H), H=1,R)
	      ELSE
	         BACKSPACE 1
	         READ(1,8) (IBLABEL(H), H=1,R)
	      ENDIF
	      IF (NOMU .EQ. 0 .AND. (ISTART .EQ. 1 .or. IPRIOR .EQ. 1)) &
	         READ(1,*) (IMU(H), H=1,R*(MAXJ-1))
	   ENDIF
	
	   IF (P .NE. 0) THEN
	      READ(1,8) TLABEL
	      IF (TLABEL .EQ. '        ') THEN
	         READ(1,8) (IALABEL(L), L=1,P)
	      ELSE
	         BACKSPACE 1
	         READ(1,8) (IALABEL(L), L=1,P)
	      ENDIF
	      
	      IF (ISTART.EQ.1)READ(1,*) (IALPHA(L), L=1,P*(MAXJ-1))
	      8 FORMAT(10A8)
	
	      IF (ICATYX .EQ. 1) THEN
	         IF (ICXW .EQ. 1) XLabel = IBLABEL(ICX)
	         IF (ICXW .EQ. 2) XLabel = IALABEL(ICX)
	      ENDIF
	   ENDIF
	
	   IF (R .NE. 0 .AND. ISTART .EQ. 1)READ(1,*)(ISIGMA(HR),HR=1,RRJ)
	
	   IF (NCON .NE. 0) READ(1,*) (ICON(L), L=1,NCPAR)
	
	   CLOSE(1)
	
	   NRP1 = 1 + R + P
	   UNO  = 1.0D0
	
	   ! get quadrature nodes & weights  
	
	   ALLOCATE(IAQ1(NQ1))
	   ALLOCATE(IAQ2(NQRR))
	
	   IF (IVSEP .EQ. 0) THEN
	      CALL QUADP(IBQ,IBQ1,IAQ,NQ1,NQ,R,IUNIF,IAQ1,IAQ2)
	   ELSE
	      CALL QUADP(IBQ,IBQ1,IAQ,NQ1,NQ,1,IUNIF,IAQ1,IAQ2)
	   ENDIF
	
	   ! GET THE DATA 
	   ! N    (total number of level 2 units)
	   ! NTOT (total number of observations)
	   ! 
	   
	   CALL READAT(filedat,N,NTOT,MAXK,MAXCOL,R,P,ALLDAT, &
	         IDNI,IDIND,YIND,IXIND,IWIND,MISS,YMISS,IXMISS,IWMISS, &
	         IWT,WTIND)
	
	   ! NALL is the number of elements for the data read in from filedat 
	   ! that ALLDAT contains 
	
	   RNTOT = DBLE(NTOT)
	
	   IF (IWT .EQ. 0) THEN
	      NALL = NRP1*NTOT  
	   ELSE
	      NALL = NRP1*NTOT + N
	   ENDIF
	
	   ! XLAB dimension is largest of maxni or npar (= p + ngam + r + rr)
	   ! or of ncon
	   ! NOTE - CAN'T CALL PRNT BEFORE THIS POINT!!
	
	   IF (MAXK .GE. NPAR) THEN
	      NXLAB = MAXK
	   ELSE
	      NXLAB = NPAR
	   ENDIF
	   
	   IF (NXLAB .LT. NCON) NXLAB = NCON
	   ALLOCATE(IXLAB(NXLAB))
	   DO K = 1,NXLAB
	      IXLAB(K) = '    '
	   END DO
	
	   ! calculate means, standard deviations, proportions and starting values
	
	   MAXJXJ = MAXJ*MAXXJ
	   ALLOCATE(ICATFQ(MAXJ))
	   ALLOCATE(ICATFQX(MAXJXJ))
	   ALLOCATE(IWORKCAT(MAXJ))
	
	   ICATFQ  = 0.0D0    ! Array assignment gets every element
	   ICATFQX = 0.0D0    ! Array assignment gets every element
	
	   II = 1
	   KI = 1
	   WSUM = 0.0D0
	   SUM2 = 0.0D0
	   
	   DO K  = 1,NALL
	      K2 = K - (II - 1)
	      IR     = MOD(K2,NRP1)
	      IF (IR .EQ. 0) IR = NRP1
	      
	      IF (IR .EQ. 1) THEN
	         IF (IWT .EQ. 1 .AND. KI .EQ. (IDNI(II*2)+1)) THEN
	            WSUM = WSUM + ALLDAT(K)*DBLE(IDNI(II*2))
	            SUM2 = SUM2 + ALLDAT(K)
	            KEEP = .FALSE.
	            KI = 1
	            II = II + 1
	         ELSE
	            KEEP = .TRUE.
	            KI = KI + 1
	         ENDIF
	      ENDIF
	      
	      IF (KEEP) THEN
	         IC     = NALL + IR
	         IC2    = NALL + NRP1 + IR
	         IC3    = NALL + NRP1 + NRP1 + IR
	
	         ! get the minimums
	
	         IF (K .LE. NRP1) THEN
	            ALLDAT(IC) = ALLDAT(K)
	         ELSE
	            IF (ALLDAT(K) .LT. ALLDAT(IC)) ALLDAT(IC) = ALLDAT(K)
	         ENDIF
	
	         ! get the maximums
	
	         IF (K2 .LE. NRP1) THEN
	            ALLDAT(IC2) = ALLDAT(K)
	         ELSE
	            IF (ALLDAT(K) .GT. ALLDAT(IC2)) ALLDAT(IC2) = ALLDAT(K)
	         ENDIF
	
	         ! get the sums
	         ! multiply by the weight if appropriate
	
	         IF (IWT .EQ. 0) THEN
	            ALLDAT(IC3) = ALLDAT(IC3) +ALLDAT(K)
	         ELSE
	            ICW = (K + (IDNI(II*2) * NRP1)-IR+1)-((KI-2)*NRP1)
	            ALLDAT(IC3) = ALLDAT(IC3)+ ALLDAT(ICW)*ALLDAT(K)
	         ENDIF
	
	         ! for Y, get the frequency in each category
	         ! and add the weight if appropriate
	
	         IF (IR .EQ. 1) THEN
	            IF (IWT .EQ. 0) THEN
	               CATADD = 1.0D0
	            ELSE
	               CATADD = ALLDAT(ICW)
	            ENDIF
	            
	            DO K3 = 1,MAXJ
	               IF (FP_EQUAL(ALLDAT(K) , ICODE(K3))) THEN
	                  ICATFQ(K3) = ICATFQ(K3) + CATADD
	                  ICAT      = K3
	               ELSE
	               ENDIF
	            END DO
	         ENDIF
	
	         ! for one of the Xs or Ws get the crosstab by Y
	         ! and add the weight if appropriate
	
	         IF (ICATYX .EQ. 1 .AND. IR .EQ. NCATYX) THEN
	            DO K4 = 1,MAXXJ
	               IF (ALLDAT(K) .GE. (CODEX(K4)-0.0001D0) &
	                  .AND. ALLDAT(K) .LE. (CODEX(K4)+0.0001D0)) THEN
	                   KIND = ((ICAT-1)*MAXXJ) + K4
	                   ICATFQX(KIND) = ICATFQX(KIND) + CATADD
	               ELSE
	               ENDIF
	            END DO
	         ENDIF
	      ELSE
	      ENDIF
	   END DO
	
	   ! calculate the means
	
	   ICCY = 0
	   
	   DO K  = 1,NRP1
	      IC3    = NALL + NRP1 + NRP1 + K
	      IF (IWT .EQ. 0) THEN
	         ALLDAT(IC3) = ALLDAT(IC3) / DBLE(NTOT)
	      ELSE
	         ALLDAT(IC3) = ALLDAT(IC3) / WSUM
	      ENDIF
	
	      !  see if the mean of the 1 random effect is 1
	      !  for intraclass correlation calculation later
	
	      IF (K .EQ. 2 .AND. R .EQ. 1 &
	         .AND. ALLDAT(IC3) .GT. 0.999D0 &
	         .AND. ALLDAT(IC3) .LT. 1.001D0) &
	         ICCY = ICCY + 1
	
	   END DO
	
	   ! get starting values for the p regression coefficient plus
	   ! the mean of the 2nd random effect (when r = 2)
	   !
	   ! use IWORKP1() and IWRKP1() as work vectors (P1 = p + r - 1)
	   !     IWP1P1() as a work vector (P1 * (P1+1) / 2)
	
	   IC = 0       
	   IF (NOMU .EQ. 0) THEN
	      P1 = P + R - 1
	   ELSE
	      P1 = P - 1
	   ENDIF
	   
	   P1P1 = (P1 * (P1+1)) / 2
	   ALLOCATE(IWORKP1(MAX(P1,NRP1)))
	   ALLOCATE(IWRKP1(P1))
	   ALLOCATE(IWP1P1(P1P1))
	   
	   IWRKP1 = 0.0D0    ! Array assignment gets every element
	   IWP1P1 = 0.0D0    ! Array assignment gets every element
	   
	   ! get the sums of squared deviations about the means
	   ! and work vectors for the regression coefficients
	
	   II = 1
	   KI = 1
	   
	   DO K  = 1,NALL
	      K2 = K - (II - 1)
	      IR     = MOD(K2,NRP1)
	      IF (IR .EQ. 0) IR = NRP1
	      
	      IF (IR .EQ. 1) THEN
	         IF (IWT .EQ. 1 .AND. KI .EQ. (IDNI(II*2)+1)) THEN
	            KEEP = .FALSE.
	            KI = 1
	            II = II + 1
	         ELSE
	            KEEP = .TRUE.
	            KI = KI + 1
	         ENDIF
	      ENDIF
	      
	      IF(KEEP) THEN
	         IC3    = NALL + NRP1 + NRP1 + IR
	         IC4    = NALL + NRP1 + NRP1 + NRP1 + IR
	         
	         IF (IWT .EQ. 0) THEN
	            ALLDAT(IC4) = ALLDAT(IC4) + ((ALLDAT(K)-ALLDAT(IC3))**2)
	            IF (IR .EQ. 1) YDEV = ALLDAT(K) - ALLDAT(IC3)
	
	            ! this assumes that the second column equals the intercept
	            ! so you don't need to get the regression coefficient for
	            ! ir eq 2
	
	            IF (IR .GT. 2) THEN
	               ! Don't allow IR2 to run off the end of the
	               ! work vector TG 1/29/01
	               IF(IR-2 .LE. UBOUND(IWORKP1,1)) THEN
	                  IWORKP1(IR - 2) = ALLDAT(K) - ALLDAT(IC3)
	               ENDIF
	               IF (IR .EQ. NRP1) THEN
	                  CALL GRMCV(IWP1P1,IWP1P1,IWORKP1,UNO,P1)
	                  DO L = 1,P1
	                     IWRKP1(L) = IWRKP1(L) +(YDEV * IWORKP1(L))
	                  END DO
	               ENDIF
	            END IF
	         
	         ELSE
	            ICW = (K+(IDNI(II*2)*NRP1) -IR +1) -((KI-2)*NRP1)
	            ALLDAT(IC4) = ALLDAT(IC4) + ALLDAT(ICW) * &
	                                ((ALLDAT(K) - ALLDAT(IC3))**2)
	            IF (IR .EQ. 1) YDEV = ALLDAT(ICW) * (ALLDAT(K)-ALLDAT(IC3))
	
	            IF (IR .GT. 2) THEN
	               ! Don't allow IR2 to run off the end of the
	               ! work vector TG 1/29/01
	               IF(IR-2 .LE. UBOUND(IWORKP1,1)) THEN
	                  IWORKP1(IR - 2) = ALLDAT(K) - ALLDAT(IC3)
	               ENDIF
	               IF (IR .EQ. NRP1) THEN
	                  CALL GRMCV(IWP1P1,IWP1P1,IWORKP1,ALLDAT(ICW),P1)
	                  DO L = 1,P1
	                     IWRKP1(L) = IWRKP1(L) + (YDEV * IWORKP1(L))
	                  END DO
	               ENDIF
	            END IF
	         ENDIF
	      ENDIF
	   END DO
	
	   ! calculate the standard deviations
	
	   DO K  = 1,NRP1
	      IC4    = NALL + NRP1 + NRP1 + NRP1 + K
	      IF (IWT .EQ. 0) THEN
	         ALLDAT(IC4) = DSQRT(ALLDAT(IC4) / DBLE(NTOT-1))
	      ELSE
	         ALLDAT(IC4) = DSQRT(ALLDAT(IC4) / (WSUM - 1.0D0))
	      ENDIF
	
	      !  see if the SD of the 1 random effect is 0
	      !  for intraclass correlation calculation later
	
	      IF (K .EQ. 2 .AND. R .EQ. 1 & 
	         .AND. ALLDAT(IC4) .GT. -0.001D0  &
	         .AND. ALLDAT(IC4) .LT. 0.001D0) &
	          ICCY = ICCY + 1
	
	   END DO
	
	   ! write out descriptives
	
	   WRITE(2,554) 
	   554 FORMAT(1x,'MIXNO - The program for mixed-effects nominal logistic regression analysis',/)
	
	   WRITE(2,"(1x,15A4)") HEAD
	
	   IF (r .ge. 1 .and. iunif .eq. 0) then
	      WRITE(2,"(/,1x,'Random-effects distribution: normal')") 
	   ELSEIF (r .ge. 1 .and. iunif .eq. 1) then
	      WRITE(2,"(/,1x,'Random-effects distribution: rectangular')") 
	   ENDIF
	
	   WRITE(2,255)
	   255 FORMAT(//,1x,'Numbers of observations',/,1x,'-----------------------',/)
	   IF (IWT .EQ. 0) THEN
	      WRITE(2,"(1x,'Level 1 observations = ',i6)")NTOT
	      IF (R .GE. 1) THEN
	         WRITE(2,"(1x,'Level 2 observations = ',i6)")N
	         WRITE(2,5610)
	         WRITE(2,"(1x,19I4)")(IDNI(I),I=2,N*2,2)
	         5610 FORMAT(//,1x,'The number of level 1 observations per level 2 unit are:',/)
	      ENDIF     
	   ELSE
	      WRITE(2,506)WSUM,NTOT
	      506 FORMAT(1x,'Level 1 observations = ',f8.2,/,1x, &
	         'Level 1 patterns =     ',I8)
	      IF (R .GE. 1) THEN
	         WRITE(2,516)SUM2,N
	         516 FORMAT(/,1x,'Level 2 observations = ',f8.2,/,1x, &
	            'Level 2 patterns =     ',I8)
	         WRITE(2,5070)
	         WRITE(2,"(1x,19I4)")(IDNI(I),I=2,N*2,2)
	         5070 FORMAT(//,1x, 'The number of level 1 patterns per level 2 pattern are:',/)
	      ENDIF     
	   ENDIF
	
	   WRITE(2,257)
	   257 FORMAT(//,1x,'Descriptive statistics for all variables',/,1x, &
	      '----------------------------------------',/)
	   WRITE(2,357)
	   357 FORMAT(1x,'Variable',5X,'     Minimum',5x,'     Maximum',5x, &
	             '        Mean',5x,' Stand. Dev.',/)
	  
	   IC  = NALL + 1
	   IC2 = NALL + NRP1 + 1
	   IC3 = NALL + NRP1 + NRP1 + 1
	   IC4 = NALL + NRP1 + NRP1 + NRP1 + 1
	   WRITE(2,377) YLABEL,ALLDAT(IC),ALLDAT(IC2),ALLDAT(IC3),ALLDAT(IC4)
	   
	   IF (R .GE. 1) THEN
	      DO H = 1,R
	         IC  = IC  + 1
	         IC2 = IC2 + 1
	         IC3 = IC3 + 1
	         IC4 = IC4 + 1
	         WRITE(2,377)IBLABEL(H),ALLDAT(IC),ALLDAT(IC2),ALLDAT(IC3),ALLDAT(IC4)
	      END DO
	   ENDIF
	   
	   IF (P .GE. 1) THEN
	      DO L = 1,P
	         IC  = IC  + 1
	         IC2 = IC2 + 1
	         IC3 = IC3 + 1
	         IC4 = IC4 + 1
	         WRITE(2,377)IALABEL(L),ALLDAT(IC),ALLDAT(IC2),ALLDAT(IC3),ALLDAT(IC4)
	      END DO
	   ENDIF
	
	   377 FORMAT(1x,A8,4(5X,F12.5))
	
	   WRITE(2,258)YLabel
	   258 FORMAT(//,1x,'Categories of the response variable ',A8,/,1x, &
	      '--------------------------------------------',/)
	   WRITE(2,457)
	   457 FORMAT(1X,'Category',5X,'   Frequency',5x,'  Proportion',/)
	
	   DO J = 1,MAXJ
	      IF (IWT .EQ. 0) THEN
	         DENUM = DBLE(NTOT)
	      ELSE
	         DENUM = WSUM
	      ENDIF
	      IWORKCAT(J) = ICATFQ(J) / DENUM
	      WRITE(2,477)ICODE(J),ICATFQ(J),IWORKCAT(J)
	      477 FORMAT(1X,F8.2,5X,F12.2,5x,f12.5)
	   END DO
	
	   IF (ICATYX .EQ. 1) THEN
	      WRITE(2,358)XLabel,YLabel
	      358 FORMAT(//,1x,'Crosstabulation of variable ',A8, &
	                 ' by the response variable ',A8,/,1x, &
	   '----------------------------------------------------------------------',/)
	      WRITE(2,372)YLabel,XLabel,(ICODE(J),J=1,MAXJ)
	      372 FORMAT(1X,8X,A8,/,1X,8X,'--------',/,1X,A8,8F8.2)
	      LABTEMP = '--------'
	      
	      DO J = 1,MAXJ
	         L2 = J*8
	         LABTEMP = LABTEMP(1:L2)//'--------'
	      END DO
	      
	      L2 = L2 + 8
	      LABTEMP = LABTEMP(1:L2)//'   Total'
	      WRITE(2,"(1X,A80)") LABTEMP
	
	      TEMPFORM = '(11X,'
	      L2 = 5
	      DO J = 1,MAXJ
	         TEMPFORM = TEMPFORM(1:L2)//"' (',F4.2,') ',"
	         L2 = L2 + 15
	      END DO
	      TEMPFORM = TEMPFORM(1:L2)//')'
	
	      DO J = 1,MAXXJ
	         TEMPSUM = 0.0D0
	         DO J2= 1,MAXJ
	            KIND = ((J2-1)*MAXXJ) + J
	            TEMPSUM = TEMPSUM + ICATFQX(KIND) 
	         END DO
	         DO J2= 1,MAXJ
	            KIND = ((J2-1)*MAXXJ) + J
	            IWORKCAT(J2) = ICATFQX(KIND) / TEMPSUM
	         END DO
	         WRITE(2,"(/,1X,F8.2,8F8.1)") CODEX(J), &
	            (ICATFQX(J3),J3=J,MAXJXJ,MAXXJ),TEMPSUM
	         WRITE(2,TEMPFORM)(IWORKCAT(J3),J3=1,MAXJ)
	      END DO
	      WRITE(2,"(/,1X,'Total   ',8F8.1)")(ICATFQ(J),J=1,MAXJ),DENUM
	   ENDIF
	
	   ! done writing out descriptives, get starting values
	
	   IF (ISTART .NE. 1) THEN
	
	      ! calculate the starting values for the regression coefficients
	
	      CALL INVS(IWP1P1,P1,DET,IWORKP1,IER)
	      CALL MPYM(IWP1P1,IWRKP1,IWORKP1,P1,P1,1,0,1)
	
	      WA = 0.0D0
	      IC3 = NALL + NRP1 + NRP1 + 1
	      IC4 = NALL + NRP1 + NRP1 + NRP1 + 1
	      
	      IF (P1 .EQ. P) THEN
	         ICA = 0
	         DO J = 2,MAXJ
	            DO L = 1,P
	               ICA = ICA+1
	               IALPHA(ICA) = (1.0D0+(J-2)*.1)*(0-IWORKP1(L)) /ALLDAT(IC4)
	               WA = WA + (IALPHA(ICA) * ALLDAT(IC3+R+L))
	            END DO
	         END DO
	         
	      ELSEIF (P1 .GT. P) THEN
	         IF (IPRIOR .EQ. 0 .AND. NOMU .EQ. 0) THEN
	            ICA = 0
	            DO J = 2,MAXJ
	               DO H = 2,R
	                  ICA = ICA+1
	                  IMU(ICA) = (1.0D0+(J-2)*.1)*(0-IWORKP1(H-1))/ALLDAT(IC4)
	                  WA = WA + (IMU(ICA) * ALLDAT(IC3+H))
	               END DO
	            END DO
	         ENDIF
	         ICA = 0
	         DO J = 2,MAXJ
	            DO L = 1,P
	               ICA = ICA+1 
	               L2 = L+1
	               IALPHA(ICA) = (1.0D0+(J-2)*.1) * &
	                               (0-IWORKP1(L2)) / ALLDAT(IC4)
	               WA = WA + (IALPHA(ICA) * ALLDAT(IC3+R+L))
	            END DO
	         END DO
	      ELSE
	      ENDIF
	
	      ! call starting value routine for intercepts, 
	      ! and random-effects variance-covariance matrix
	
	      IF(R .GT. 0) THEN
	         CALL STARTV2(IMU,ISIGMA,R,MAXJ, &
	           ICATFQ,DENUM,WA,NFN,IPRIOR,IDIAG,IVSEP,IVCAT)
	      ENDIF
	
	   ENDIF
	   ! write out AND save the starting values
	
	   WRITE(2,"(//,1x,'Starting values',/,1x,'---------------',/)")
	   IF (r .GE. 1 .AND. NOMU .EQ. 0) then
	       WRITE(2,"(1x,'mean       ',10F7.3)") (IMU(H),H=1,R*(MAXJ-1))
	       ALLOCATE(IMU0(RJ))
	       CALL RELOC(IMU,IMU0,RJ,1,0)
	   ENDIF
	   
	   IF (P .GT. 0) THEN
	       IF(R.EQ.0) IALPHA(1) = IMU(1)
	       WRITE(2,"(1x,'covariates ',10F7.3)") (IALPHA(L),L=1,P*(MAXJ-1))
	       ALLOCATE(IALPHA0(PJ))
	       CALL RELOC(IALPHA,IALPHA0,PJ,1,0)
	   ENDIF
	   
	   IF (r .ge. 1) then
	       WRITE(2,"(1x,'var. terms ',10F7.3)") (ISIGMA(HR),HR=1,RRJ)
	       ALLOCATE(ISIGMA0(RRJ))
	       CALL RELOC(ISIGMA,ISIGMA0,RRJ,1,0)
	   ENDIF
	   
	   ! set up for main loop
	   
	   ITLAST = 0 
	   RIDGE  = 0.0D0
	   IFIN   = 1
	   ! NORI   = 0
	   IUN    = 6
	   ND     = 0
	   rloglp = -999999999999999.0D0
	
	   IT    = 0
	
	   ! isig   = 1 when INVS has problems (0 otherwise)
	   ! ier    = 1 when INVS has problems (0 otherwise)
	   ! nonpos = 0 when CHSKY detects a non pos-def matrix (1 when pos-def)
	
	   ISIG   = 0
	   IER    = 0
	   NONPOS = 1
	   IRBAD  = 0
	   IGO = 1
	   INITIALIZED = 1
	   
	   DEALLOCATE(IXMISS, IWMISS, IAQ1,IAQ2, IWP1P1, ICATFQ, &
	      ICATFQX,IWORKCAT,IWORKP1,IWRKP1, IXIND, IWIND )
	   IF(ALLOCATED(CODEX)) DEALLOCATE(CODEX)
	
	END SUBROUTINE INIT
	
	   
	! ************************************************
	!               **********************                         
	!               * FUNCTION MainLoop()
	!               **********************                         
	!                                                                       
	! Perform the work of iterating through the dataset to converge
	! on a stable solution for the regression.  This is called repeatedly
	! by either the DLL caller, or the EXE shell until it returns a 0,
	! signalling that the calculation is complete.
	!
	! MAINLOOP()
	!                                                                       
	! The SAVE attribute is used for the local variables so that they
	! will persist between invokations in the DLL.  Otherwise, the counters
	! would reset on each call.  This means that there can only be one
	! session (thread) active in the dll at any given time.
	! ************************************************
	INTEGER FUNCTION MainLoop()
	USE MIXLIB
	IMPLICIT NONE
	   
	   SAVE
	   
	   ! **********************************************************
	   !  Perform a single iteration of the main loop
	   ! **********************************************************
	   
	   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	   ! LOCAL VARS ONLY - SHARED VARS ARE DECLARED
	   ! BELOW THE 'MODULE' STATEMENT AT THE TOP OF THE 
	   ! FILE
	   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	   
	   CHARACTER (LEN=8), ALLOCATABLE :: TEMPLABEL(:)
	  
	   INTEGER NQRR, L, J, K, IC, IC2, &
	         I, ICOUNT, IDI, IH, IND, IND2, INDD, INDD2, &
	         IW, IX, KK, KKK, NII, NSAME, NSAMES, &
	         H,HR,RNII,PNII,Q,R2,H2, TEMP_I
	   
	   REAL*8 AQALL, BIGCOR, DEN, DERIV, DET, &
	         HPROB, PRJL, PSUM, RLOGDIFF, QMULT, &
	         RSAMES, SCAL, SIGN, STEP, SUMW, &
	         WA, WTSUBI, XMU, XTB, ZJ
	         
	   REAL (KIND = 8),POINTER ::IA(:)
	   
	   REAL (KIND = 8), ALLOCATABLE :: IWRKR(:),IWORKR2(:), &
	         IAQ1(:),IAQ2(:),ICHWRK(:), ICOREC(:), IDERP(:),IDERQ(:), &
	         IYI(:),IXI(:),IWI(:),ILIK(:)
	      
	   LOGICAL CATJ
	   
	   IF(INITIALIZED .EQ. 0) THEN
	      CALL POST_ERROR("Mainloop() routine called before Init() " // PROGNAME)
	      RETURN
	   ELSEIF(INITIALIZED .EQ. 2) THEN
	      CALL POST_ERROR("Mainloop() routine called after completion. " // PROGNAME)
	      RETURN
	   ENDIF
	   
	   IF (IGO .EQ. 1 .OR. ITLAST .EQ. 1) THEN
	  
	      IT  = IT + 1
	
	      ! do step-halving for first 0 iterations
	
	      IF (IT .LE. 0) THEN
	         STEP = 0.5D0
	      ELSE
	         STEP = 1.0D0
	      ENDIF 
	
	      ! CHECK TO SEE IF LESS RANDOM-EFFECTS ARE TO
	      ! BE ESTIMATED AND MODIFY IALPHA POINTER
	      ! R, P, AND NPAR ACCORDINGLY
	
	      IF (IRBAD .GE. 1) THEN
	         
	         IRBAD   = 0
	         ISIG    = 0
	         IER     = 0
	         NONPOS  = 1
	         ! NORI    = 0
	         ! rloglp  = -999999999999999.0D0
	         
	         ! Must insert an element into first position 
	         ! will be filled in below from IALPHA0
	         
	         TEMP_I = UBOUND(IALPHA,1)
	         DEALLOCATE(IALPHA)
	         ALLOCATE(IALPHA(TEMP_I+1))
	         
	         ! Shift the contents of IALABEL up one and fill
	         ! in the first location with IBLABEL.
	         
	         TEMP_I = UBOUND(IALABEL,1)
	         IF (ALLOCATED(TEMPLABEL)) DEALLOCATE(TEMPLABEL)
	         ALLOCATE(TEMPLABEL(TEMP_I+1))
	         TEMPLABEL(2:) = IALABEL   ! Array copy to temp
	         DEALLOCATE(IALABEL)
	         ALLOCATE(IALABEL(TEMP_I+1))
	         IALABEL = TEMPLABEL       ! Array copy
	         IALABEL(1) = IBLABEL(UBOUND(IBLABEL,1))
	         
	         R  = R-1
	         RR = (R * (R+1)) / 2
	         P  = P+1      
	         IF (IPRIOR .EQ. 0 .AND. NOMU .EQ. 0) THEN
	            npar = (p+r)*(MAXJ-1)+rr*(MAXJ-1)
	         ELSE
	            npar = p*(MAXJ-1)+rr*(MAXJ-1)
	         ENDIF
	
	         ! get NEW quadrature nodes & weights  
	
	         NQRR = R*(NQ1**R)
	         IF (ALLOCATED(IAQ1)) THEN
	            IF(UBOUND(IAQ1,1) .NE. NQ1) THEN
	               DEALLOCATE(IAQ1)
	               ALLOCATE(IAQ1(NQ1))
	            ENDIF
	            IF(UBOUND(IAQ2,1) .NE. NQRR) THEN
	               DEALLOCATE(IAQ2)
	               ALLOCATE(IAQ2(NQRR))
	            ENDIF
	         ELSE
	            ALLOCATE(IAQ1(NQ1))
	            ALLOCATE(IAQ2(NQRR))
	         ENDIF
	
	         IF (IVSEP .EQ. 0) THEN
	            CALL QUADP(IBQ,IBQ1,IAQ,NQ1,NQ,R,IUNIF,IAQ1,IAQ2)
	         ELSE
	            CALL QUADP(IBQ,IBQ1,IAQ,NQ1,NQ,1,IUNIF,IAQ1,IAQ2)
	         ENDIF
	
	         ! get ORIGINAL starting values
	
	         IF (r .ge. 1) then
	             IF (NOMU.EQ.0) CALL RELOC(IMU0,IMU,R*(MAXJ-1),1,0)
	             CALL RELOC(ISIGMA0,ISIGMA,RRJ,1,0)
	         ENDIF
	         
	         IF (P .GT. 0) THEN
	            ! note that R now equals the old R minus 1
	            ! so that imu0+R equals the old imu0+R-1 
	            ! which should be the last element of mu
	             IALPHA(1) = IMU0(R*(MAXJ-1)+1)
	             CALL RELOC(IALPHA0,IALPHA(2),(P-1)*(MAXJ-1),1,0)
	         ENDIF
	
	      ENDIF
	
	      ! get the exponential of the diagonal of ISIGMA
	
	      IF (IFIN .EQ. 0) THEN
	         IH = 0
	         DO H=1,R
	            DO H2=1,H
	               IH = IH + 1
	               IF (H .EQ. H2) ISIGTAU(H) = DLOG(ISIGMA(IH))
	            END DO
	         END DO
	      ENDIF
	
	      ! ***
	      ! LIKELIHOOD INITIALIZATION
	      ! ***
	      !  calculate and save the conditional likelihood - ILIK()
	      !                     the marginal likelihood    - hprob
	      !                     the log-likelihood value   - rlogl
	
	      IDER2 = 0.0D0     ! Array assignment gets every element
	      IDER = 0.0D0      ! Array assignment gets every element
	      rlogl = 0.0d0
	
	      IF (IPRIOR .EQ. 1) THEN
	         AQALL    = 0.0D0
	         DO Q = 1,NQ
	            IAQEMP(Q) = 0.0D0
	         END DO
	      ENDIF
	
	      ! *************************************************************
	      ! GO OVER ALL SUBJECTS
	      !    IUN = UNIT NUMBER FOR OUTPUT (0 = SCREEN)
	      ! *************************************************************
	      ICOUNT = 0
	
	      IF (ALLOCATED(IDERP)) THEN
	         IF(UBOUND(IDERP,1) .NE. NPAR) THEN
	            DEALLOCATE(IDERP)
	            ALLOCATE(IDERP(NPAR))
	         ENDIF
	      ELSE
	         ALLOCATE(IDERP(NPAR))
	      ENDIF
	
	      DO I = 1,N
	
	         hprob = 0.0d0
	         IDERP = 0.0D0    ! Array assignment gets every element at once
	   
	         IC2  = 2*(I-1)+1
	         IDI  = IDNI(IC2)
	         IC2  = IC2+1
	         NII  = IDNI(IC2)
	         RNII = NII * R
	         PNII = NII * P
	   
	         ! THE YI(K) VECTOR  K = 1 ... NI(I)   DEPENDENT VARIABLE VECTOR
	   
	         IF (ALLOCATED(IYI)) THEN
	            IF(UBOUND(IYI,1) .NE. NII) THEN
	               DEALLOCATE(IYI)
	               ALLOCATE(IYI(NII))
	            ENDIF
	         ELSE
	            ALLOCATE(IYI(NII))
	         ENDIF
	         
	         DO K = 1,NII
	            IC2         = ICOUNT + (NRP1 * (K-1)) + 1 
	            IYI(K) = ALLDAT(IC2)
	         END DO
	   
	         !  calculate the number of level-2 units with no variance 
	         !  in terms of level-1 units Y values
	   
	         IF (IT .EQ. 1) THEN
	            IF (I  .EQ. 1) NSAMES = 0
	            CALL YSAME(IYI,NII,NSAME)
	            NSAMES = NSAMES  + NSAME
	            
	            IF (I  .EQ. N .AND. R .GT. 0) THEN
	               RSAMES = 100.0D0 * (DBLE(NSAMES) / DBLE(N))
	               IF (IWT.EQ.0) THEN
	                  WRITE(2,508)NSAMES,RSAMES
	                  508 FORMAT(//,1x,'==> The number of level 2 observations with non-varying responses', &
	              /,5x,'= ',I6,' ( ',F6.2,' percent )')  
	               ELSE
	                  WRITE(2,509)NSAMES,RSAMES
	                  509 FORMAT(//,1x,'==> The number of level 2 patterns with non-varying responses', &
	             /,5x,'= ',I6,' ( ',F6.2,' percent )')  
	               ENDIF
	            ENDIF
	         ENDIF
	   
	         ! THE X(K, H) MATRIX  K = 1 .. NI(I)  H = 1 .. R
	   
	         IF (R .GE. 1) THEN
	            
	            IF (ALLOCATED(IXI)) THEN
	               IF(UBOUND(IXI,1) .NE. RNII) THEN
	                  DEALLOCATE(IXI)
	                  ALLOCATE(IXI(RNII))
	               ENDIF
	            ELSE
	               ALLOCATE(IXI(RNII))
	            ENDIF
	            
	            IC = 0
	            DO H = 1,R
	               DO K = 1,NII
	                  IC  = IC + 1
	                  IC2 = ICOUNT + (NRP1 * (K-1) + H+1) 
	                  IXI(IC) = ALLDAT(IC2)
	               END DO
	            END DO
	            ELSE
	         ENDIF
	   
	         ! THE W(K, L) MATRIX  K = 1 .. NI(I)  L = 1 .. P
	   
	         IF (P .GE. 1) THEN
	            
	            IF (ALLOCATED(IWI)) THEN
	               IF(UBOUND(IWI,1) .NE. PNII) THEN
	                  DEALLOCATE(IWI)
	                  ALLOCATE(IWI(PNII))
	               ENDIF
	            ELSE
	               ALLOCATE(IWI(PNII))
	            ENDIF
	            
	            IC = 0
	            DO L  = 1,P   
	               DO K  = 1,NII
	                  IC  = IC + 1
	                  IC2 = ICOUNT + (NRP1 * (K-1) + L+R+1)
	                  IWI(IC) = ALLDAT(IC2)
	               END DO
	            END DO
	         ENDIF
	   
	         ICOUNT = ICOUNT + NII + RNII + PNII 
	   
	         ! THE WEIGHT FOR EACH LEVEL-2 UNIT
	   
	         IF (IWT .EQ. 1) THEN
	            ICOUNT = ICOUNT + 1
	            WTSUBI = ALLDAT(ICOUNT)
	         ELSE
	            WTSUBI = 1.0D0
	         ENDIF
	   
	         ! PRINT THE DATA OF THE 1ST NPR SUBJECTS
	         ! TO THE UNIT NUMBER IUN ON THE FIRST ITERATION ONLY
	   
	         IF (I .LE. NPR .AND. IT .EQ. 1) THEN
	            WRITE(IUN,11)N
	            11 FORMAT(1X,'TOTAL NUMBER OF SUBJECTS = ',I6,/)
	            WRITE(IUN,12) IDI,NII
	            12 FORMAT(1X,'DATA FOR SUBJECT ',I10,' WHO HAS',I6,' OBSERVATIONS')
	            CALL PRNT(IUN,IYI,NII,1,0,IXLAB,IXLAB,ND,HEAD,1,80,5,1,1, &
	         60H Y VECTOR                                                   )
	            
	            IF (R .GT. 0) THEN
	               CALL PRNT(IUN,IXI,NII,R,0,IXLAB,IXLAB,ND,HEAD,1,80,5,1,1, &
	         60H X MATRIX                                                   )
	            ENDIF
	            
	            IF (P .GT. 0) THEN
	               CALL PRNT(IUN,IWI,NII,P,0,IXLAB,IXLAB,ND,HEAD,1,80,5,1,1, &
	         60H W MATRIX                                                   )
	            ENDIF
	         
	            IF (IWT .EQ. 1) THEN
	               WRITE(IUN,42) IDI,WTSUBI
	               42 FORMAT(1X,'WEIGHT FOR SUBJECT ',I10,' IS',F12.5)
	            ENDIF
	         ENDIF
	   
	   
	         ! GO OVER THE QUADRATURE POINTS
	   
	         IF (ALLOCATED(IDERQ)) THEN
	            IF(UBOUND(IDERQ,1) .NE. NPAR) THEN
	               DEALLOCATE(IDERQ)
	               ALLOCATE(IDERQ(NPAR))
	            ENDIF
	         ELSE
	            ALLOCATE(IDERQ(NPAR))
	         ENDIF
	         
	         IF (ALLOCATED(ILIK)) THEN
	            IF(UBOUND(ILIK,1) .NE. NQ) THEN
	               DEALLOCATE(ILIK)
	               ALLOCATE(ILIK(NQ))
	            ENDIF
	         ELSE
	            ALLOCATE(ILIK(NQ))
	         ENDIF
	         
	         DO q=1,nq 
	   
	            psum     = 0.0d0
	            ILIK(Q)  = 0.0d0
	            
	            IDERQ = 0.0D0        ! Array assignment gets every element at once
	      
	            DO k=1,nii 
	      
	               IPRJ  = 0.0D0        ! Array assignment gets every element at once
	               INUMJ = 0.0D0        ! Array assignment gets every element at once
	         
	               ! get the MAXJ numerators NUMJ for the current level-1 unit
	               !     and denominator DEN 
	         
	               DEN       = 1.0D0
	               INUMJ(1) = 1.0D0
	         
	               DO J=2,MAXJ
	         
	                  wa  = 0.0D0
	                  IF (P .GE. 1) THEN
	                     DO l=1,p
	                        iw = k  + (l-1)*nii
	                        IA => IALPHA(P*(J-2)+L :)
	                        wa = wa + IA(1)*IWI(iw)
	                     END DO
	                  ENDIF
	            
	                  xmu = 0.0D0
	                  xtb = 0.0D0
	                  IF (R .GE. 1) THEN
	                     DO H = 1,R
	                        IF (IVSEP .EQ. 0) THEN
	                           H2 = Q + (h-1)*NQ
	                        ELSE
	                           H2 = Q
	                        ENDIF
	                        IWORKR(H) = IBQ(H2)
	                     END DO
	                     
	                     IF (IDIAG .EQ. 0 .AND. IVSEP .EQ. 0) THEN
	                        IF (IVCAT .EQ. 0) THEN
	                           IA =>ISIGMA
	                        ELSE
	                           IA =>ISIGMA(RR*(J-2)+1 :)
	                        ENDIF
	                        CALL MPYM(IA,IWORKR,IWORKR,R,R,3,0,1)
	                     ELSE
	                        IF (IVCAT .EQ. 0) THEN
	                           IA =>ISIGMA
	                        ELSE
	                           IA => ISIGMA(R*(J-2)+1:)
	                        ENDIF
	                        CALL MPYM(IA,IWORKR,IWORKR,R,R,2,0,1)
	                     ENDIF
	                     
	                     DO h=1,r
	                        ix  = k   + (h-1)    * nii
	                        IF (NOMU .EQ. 0) THEN
	                           IA => IMU(R*(J-2)+H:) 
	                           xmu = xmu + IA(1)    * IXI(ix)
	                        ELSE
	                           XMU = 0.0D0
	                        ENDIF
	                        xtb = xtb + IWORKR(H) * IXI(ix)
	                     END DO
	                  ELSE
	                  ENDIF
	            
	                  ZJ = DEXP(WA + XMU + XTB)
	                  INUMJ(J) = ZJ
	                  DEN = DEN + ZJ
	               END DO
	         
	               catj = .TRUE.
	               DO J=1,MAXJ
	                  PRJL = DLOG(INUMJ(J)) - DLOG(DEN)
	                  IPRJ(J) = DEXP(PRJL)
	                  IF (FP_EQUAL(IYI(K) , ICODE(J))) THEN 
	                     PSUM = PSUM + PRJL
	                     catj = .FALSE.
	                  ENDIF
	               END DO
	               
	               IF (catj) THEN
	                  write(6,*)'nominal response outside range for subject ',idi
	               ENDIF
	         
	               IND = 0
	         
	               IF (R .GE. 1 .and. IPRIOR .EQ. 0 .AND. NOMU .EQ. 0) THEN
	                  DO J=2,MAXJ
	                     IF (FP_EQUAL(IYI(K) , ICODE(J))) THEN 
	                        DERIV = 1.0D0 - IPRJ(J)
	                     ELSE
	                        DERIV = 0.0D0 - IPRJ(J)
	                     ENDIF
	                     DO h=1,r
	                        IND = IND+1
	                        ix  = k + (h-1) * nii
	                        IDERQ(IND) = IDERQ(IND) + deriv*IXI(ix)
	                     END DO
	                  END DO
	               ENDIF
	      
	               IF (P .GE. 1) THEN
	                  DO J=2,MAXJ
	                     IF (FP_EQUAL(IYI(K) , ICODE(J))) THEN 
	                        DERIV = 1.0D0 - IPRJ(J)
	                     ELSE
	                        DERIV = 0.0D0 - IPRJ(J)
	                     ENDIF
	                     DO l = 1,p
	                        IND= IND + 1
	                        iw = k + (l-1)*nii
	                        IDERQ(IND) = IDERQ(IND) + deriv*IWI(iw)
	                     END DO
	                  END DO
	               ENDIF
	      
	               IF (R .GE. 1) THEN
	      
	                  ! ALLOCATE  r2 = r * r 
	                  R2 = R * R
	                  
	                  IF (ALLOCATED(IWRKR)) THEN
	                     IF(UBOUND(IWRKR,1) .NE. R) THEN
	                        DEALLOCATE(IWRKR)
	                        ALLOCATE(IWRKR(R))
	                     ENDIF
	                  ELSE
	                     ALLOCATE(IWRKR(R))
	                  ENDIF
	                  
	                  IF (ALLOCATED(IWORKR2)) THEN
	                     IF(UBOUND(IWORKR2,1) .NE. R2) THEN
	                        DEALLOCATE(IWORKR2)
	                        ALLOCATE(IWORKR2(R2))
	                     ENDIF
	                  ELSE
	                     ALLOCATE(IWORKR2(R2))
	                  ENDIF
	                  
	                  DO H = 1,R
	                     ix  = k + (h-1)*nii
	                     IF (IVSEP .EQ. 0) THEN
	                        H2 = Q + (h-1)*NQ
	                     ELSE
	                        H2 = Q
	                     ENDIF
	                     IWORKR(H) = IBQ(H2)
	                     IWRKR(H) = IXI(ix)
	                  END DO
	        
	                  CALL KMPY(IWRKR,IWORKR,IWORKR2,R,1,0,1,R)   
	                  IF (IDIAG .EQ. 0 .AND. IVSEP .EQ. 0) THEN
	                     CALL CHAMS(IWORKR2,IWORKR2,R,0,3)
	                  ELSE
	                     CALL CHAMS(IWORKR2,IWORKR2,R,0,2)
	                  ENDIF
	      
	                  MJ = 2
	                  IF (IVCAT .EQ. 1) MJ=MAXJ
	                  DO J=2,MJ
	                     IF (IVCAT .EQ. 0) THEN 
	                        IF (FP_EQUAL(IYI(K) , ICODE(1))) THEN 
	                           DERIV = 0.0D0 - (1.0D0 - IPRJ(1))
	                        ELSE
	                           DERIV = IPRJ(1)
	                        ENDIF
	                     ELSE
	                        IF (FP_EQUAL(IYI(K) , ICODE(J))) THEN 
	                           DERIV = 1.0D0 - IPRJ(J)
	                        ELSE
	                           DERIV = 0.0D0 - IPRJ(J)
	                        ENDIF
	                     ENDIF
	                     IH = 0
	                     DO h  = 1,r
	                        IF (IDIAG .EQ. 1 .OR. IVSEP .EQ. 1) THEN 
	                           IH  = IH + 1
	                           IND  = IND + 1
	                           IDERQ(IND) = IDERQ(IND) + deriv * IWORKR2(IH)
	                        ELSE
	                           DO h2 = 1,h
	                              IH  = IH  + 1
	                              iC  = iC  + 1
	                              IF (IFIN .NE. 0 .OR. H .NE. H2) THEN
	                                 ! estimate sigma 
	                                  IND = IND + 1
	                                  IDERQ(IND) = IDERQ(IND) +deriv*IWORKR2(IH)
	                              ELSE
	                              ! estimate tau = ln(diagonal elements of sigma)
	                              IF (DERIV .LT. 0.0D0 .AND. IWORKR2(IH) &
	                                  .LT. 0.0D0)SIGN =  1.0D0
	                              IF (DERIV .LT. 0.0D0 .AND. IWORKR2(IH) &
	                                  .GE. 0.0D0)SIGN = -1.0D0
	                              IF (DERIV .GE. 0.0D0 .AND. IWORKR2(IH) &
	                                  .LT. 0.0D0)SIGN = -1.0D0
	                              IF (DERIV .GE. 0.0D0 .AND. IWORKR2(IH) &
	                                  .GE. 0.0D0)SIGN =  1.0D0
	                                  IND = IND + 1
	                                  IDERQ(IND) = IDERQ(IND) + SIGN *  &
	                                                   DEXP(DLOG(DABS(deriv)) &
	                                                  +DLOG(ISIGMA(IH)) &
	                                                  +DLOG(DABS(IWORKR2(IH))))
	                              ENDIF
	                           END DO
	                        ENDIF
	                     END DO
	                  END DO
	               ENDIF
	      
	            END DO  !  DO k=1,nii 
	      
	            ILIK(Q) = DEXP(psum)
	            hprob  = hprob + DEXP ( psum + DLOG(IAQ(Q)))
	            IF (IPRIOR .EQ. 1) IAQEMP(Q) = IAQEMP(Q) +  &
	                 DEXP( (DLOG(WTSUBI)) + (DLOG(ILIK(Q))) + (DLOG(IAQ(Q))))
	      
	            DO L  = 1,npar
	               IDERP(L) = IDERP(L) +IDERQ(L)*ILIK(Q)*IAQ(Q)  
	            END DO
	       
	         END DO   ! DO q=1,nq 
	   
	         IF (HPROB .LT. .1D-305) HPROB = .1D-305
	         
	         rlogl = rlogl + WTSUBI*DLOG(hprob)
	         IF (IPRIOR .EQ. 1)  &
	            AQALL = AQALL + DEXP( (DLOG(WTSUBI)) + (DLOG(HPROB)))
	   
	         scal       = DEXP( 0.0d0 - DLOG(hprob))
	         DO l = 1,npar
	            IDER(L)  = IDER(L) + IDERP(L)* DEXP(DLOG(WTSUBI) - DLOG(hprob))
	            IDERP(L) = IDERP(L)*scal
	         END DO
	   
	         CALL grmcv(IDER2,IDER2,IDERP,wtsubi,npar)
	   
	         ! ***************************************************************************
	         ! write out the BAYES estimates at final iteration (IFIN=2or3) ONLY IF IRES=1 
	   
	         IF (IFIN .GE. 2 .AND. IRES .EQ. 1) THEN
	            IH = 0
	            
	            ITHETA  = 0.0D0    ! Array assignment gets every element
	            ITHETAV = 0.0D0    ! Array assignment gets every element
	            
	            DO Q = 1,NQ
	               QMULT = ILIK(Q) * IAQ(Q)
	               DO H = 1,IRT
	                  H2 = Q + (h-1)*NQ
	                  ITHETA(H) = ITHETA(H) + IBQ(H2)*QMULT
	               END DO
	            END DO
	            
	            DO H = 1,IRT
	               ITHETA(H) = ITHETA(H)/HPROB
	            END DO
	   
	            DO Q=1,NQ
	               QMULT = ILIK(Q) * IAQ(Q)
	               DO H = 1,IRT
	                  H2 = Q + (h-1)*NQ
	                  IWORKR(H) = IBQ(H2) - ITHETA(H)
	               END DO
	               CALL GRMCV(ITHETAV,ITHETAV,IWORKR,QMULT,IRT)
	            END DO
	   
	            IH = 0
	            DO H = 1,IRT
	               DO H2 = 1,H
	                  IH = IH+1
	                  ITHETAV(IH)=ITHETAV(IH)/HPROB
	               END DO
	            END DO
	            
	            IF (IWT .EQ. 0) THEN
	               WRITE(5,"(2I15)")IDI,NII
	            ELSEIF (IWT .EQ. 1) THEN
	               WRITE(5,"(I15,F12.5,I15)")IDI,WTSUBI,NII
	            ENDIF
	            WRITE(5,"(5F15.6)")(ITHETA(H),H=1,IRT)
	            WRITE(5,"(5F15.6)")(ITHETAV(HR),HR=1,IRTT)
	   
	         ENDIF
	   
	      END DO    !   DO I = 1,N
	
	      ! ****************************************************
	      !  done with subjects                                *
	      ! ****************************************************
	
	
	      ! CALL PRNT(IUN,IDER,NPAR,1,0,IXLAB,IXLAB,ND,HEAD,1,
	      !    +80,5,1,1,
	      !    +60H derivatives                                                )
	
	      IF (ALLOCATED(ICOREC)) THEN
	         IF(UBOUND(ICOREC,1) .NE. NPAR) THEN
	            DEALLOCATE(ICOREC)
	            ALLOCATE(ICOREC(NPAR))
	         ENDIF
	      ELSE
	         ALLOCATE(ICOREC(NPAR))
	      ENDIF
	
	      !   ADD RIDGE TO DIAGONAL ELEMENTS OF DER2 PRIOR TO INVERSION
	      !   (unless at final iteration itlast=1)
	
	      IF(ITLAST.EQ.0)THEN
	         DO KKK=1,NPAR
	            KK=(KKK*(KKK+1))/2
	            IDER2(KK)=IDER2(KK) +RIDGE*IDER2(KK)
	         ENDDO
	      ENDIF
	
	      ! check to see if the matrix of second derivatives is
	      ! positive-definite
	      ! (nonpos = 1) or not positive definite (nonpos = 0)
	
	      IF (ALLOCATED(ICHWRK)) THEN
	         IF(UBOUND(ICHWRK,1) .NE. NPARR) THEN
	            DEALLOCATE(ICHWRK)
	            ALLOCATE(ICHWRK(NPARR))
	         ENDIF
	      ELSE
	         ALLOCATE(ICHWRK(NPARR))
	      ENDIF
	      
	      CALL CHSKY(IDER2,ICHWRK,npar,NONPOS)
	      CALL INVS(IDER2,npar,det,ICOREC,IER)
	      CALL MPYM(IDER2,IDER,ICOREC,npar,npar,1,0,1)
	
	      ! DO STEP-HALVING FOR FIRST 0 ITERATIONS
	
	      DO L = 1,NPAR
	         ICOREC(L) = STEP*ICOREC(L)
	      END DO
	
	      BIGCOR = 0.0D0
	      IND = 0
	
	      IF (r .ge. 1 .and. IPRIOR .EQ. 0 .AND. NOMU .EQ. 0) then
	         DO J=2,MAXJ
	            DO h = 1,r
	               IND = IND+1
	               IMU(IND)    = IMU(IND) + ICOREC(IND)
	               IF (DABS(ICOREC(IND)) .GT. BIGCOR)  &
	                   BIGCOR = DABS(ICOREC(IND))
	            END DO
	         END DO
	      ENDIF   
	
	      IF (p .ge. 1) then
	         IND2 = 0
	         DO J = 2,MAXJ
	            DO l = 1,p
	               IND = IND+1
	               IND2 = IND2+1
	               IALPHA(IND2) = IALPHA(IND2) +ICOREC(IND)
	               IF (DABS(ICOREC(IND)) .GT. BIGCOR)  &
	                   BIGCOR = DABS(ICOREC(IND))
	            END DO
	         END DO
	      ENDIF
	
	      IF (r .ge. 1) then
	         IND2 = 0
	         INDD = 1
	         INDD2= 0
	         MJ = 2
	         IF (IVCAT .EQ. 1) MJ=MAXJ
	         DO J = 2,MJ
	            DO h = 1,rr
	               IND = IND+1
	               IND2 = IND2+1
	               
	               IF (H .EQ. INDD) THEN
	                  ! on the diagonal
	                  IF (IDIAG .EQ. 0 .AND. IVSEP .EQ. 0) THEN
	                     INDD2 = INDD2 + 1
	                     INDD  = INDD  + (INDD2 + 1)
	                  ELSE
	                     INDD  = INDD  + 1
	                  ENDIF
	         
	                  ! make sure that the first variance term is positive 
	                  ! (COMMENT OUT)
	                  
	                  IF (IFIN .NE. 0) THEN
	                     IF (H.EQ.1 .AND. RR .EQ. 1 .AND.  &
	                        (0.0D0-ICOREC(IND) .GT. ISIGMA(IND2))) THEN
	                        ISIGMA(IND2) = 0.10D0* DABS(ISIGMA(IND2) &
	                            + ICOREC(IND))
	                     ELSE
	                        ISIGMA(IND2) = ISIGMA(IND2) + ICOREC(IND)
	                     ENDIF
	   
	                     ! DON'T shift to estimation of tau if sigma gets small 
	                     ! INCREASE THE RIDGE ONCE INSTEAD
	                     ! IF(ISIGMA(H) .LE. .1D0 .AND. NORI .EQ. 0) THEN
	                     !    IFIN=0
	                     !    NORI=NORI+1
	                     !    RIDGE = RIDGE +.1D0
	                     !    WRITE(2,869)
	                     !    869 FORMAT(///,1X,'==> Reparameterized Estimation Occurred')
	                     ! ENDIF
	                  
	                  ELSE
	                     ISIGTAU(INDD2) = ISIGTAU(INDD2) + ICOREC(IND)
	                     ISIGMA(IND2)   = DEXP(ISIGTAU(INDD2))
	   
	                     ! reduce number of random effects if sigma get too small
	                     IF(ISIGMA(IND2) .LE. .000000001D0) ISIG=1
	                  ENDIF
	               ELSE
	                  ! off the diagonal
	                  ISIGMA(IND2) = ISIGMA(IND2) + ICOREC(IND)
	               ENDIF 
	               
	               IF (DABS(ICOREC(IND)) .GT. BIGCOR)  &
	                   BIGCOR = DABS(ICOREC(IND))
	            END DO
	         END DO
	      ENDIF
	
	      ! PRINT OUT RESULTS FOR THIS ITERATION
	      ! unless iteration was for computation of information matrix only
	      ! (ITLAST = 1)
	
	      IF (ITLAST .EQ. 1) THEN
	         MainLoop = 0
	         INITIALIZED = 2
	         RETURN     ! escape from the loop
	      ENDIF
	
	      WRITE(IUN,65)IT
	      65 FORMAT(/,1X,'**************',/,1X,'ITERATION ',I4,/, &
	                1X,'**************',/)
	
	      WRITE(IUN,66)rlogl
	      66 FORMAT(1X,'Log Likelihood   = ',F12.3,//)
	
	      ! save current rlogl to check at next iteration
	
	      RLOGDIFF = RLOGL - RLOGLP
	      RLOGLP   = RLOGL
	
	      ! calculate versions of RLOGL 
	      AIC  = RLOGL - RNPAR
	      SBC  = RLOGL - 0.5 * RNPAR * DLOG(DBLE(N))
	      DEV  = -2*RLOGL
	      AICD = -2*(RLOGL-RNPAR)
	      SBCD = -2*(RLOGL-0.5*RNPAR * DLOG(DBLE(N)))
	
	      IF (r .ge. 1 .AND. NOMU .EQ. 0) then
	         WRITE(IUN,903) (IMU(H),H=1,R*(MAXJ-1))
	         903 FORMAT(1X,'mu         ',7F10.6)
	      ENDIF
	      
	      IF (P .GT. 0) THEN
	         WRITE(IUN,904) (IALPHA(L),L=1,P*(MAXJ-1))
	         904 FORMAT(1X,'alpha      ',7F10.6)
	      ENDIF
	      
	      IF (r .ge. 1) then
	         WRITE(IUN,905) (ISIGMA(HR),HR=1,RRJ)
	         905 FORMAT(1X,'sigma      ',7F10.6)
	      ENDIF
	      
	      WRITE(IUN,907) (ICOREC(L),L=1,NPAR)
	      907 FORMAT(1X,'corrections',7F10.6)
	
	      ! check to see if the log-likelihood doesn't increase
	      ! and if so then increase the ridge by 0.1 each time 
	      ! check to see if there are numerical problems and if so
	      ! reduce the number of random-effects by 1
	
	      IF(RLOGDIFF .LT. 0.0D0) RIDGE = RIDGE + 0.1D0
	
	      IF(IER .EQ. 1 .OR. NONPOS .EQ. 0 .OR. ISIG .EQ. 1)IRBAD = IRBAD+1
	      
	      IF(IRBAD .EQ. 1) THEN
	         WRITE(6,878)IT
	         WRITE(2,879)IT
	         878 FORMAT(///,1X,'==> WARNING!  Estimation Difficulties at Iteration', &
	               I4,/,1x,'==> will proceed with one less random effect')
	         879 FORMAT(///,1X,'==> WARNING!  Estimation Difficulties Occurred at Iteration', &
	               I4,/,1x,'==> The model was fit with one less random effect than was requested')
	      ENDIF      
	
	      ! change weights for empirical prior
	
	      IF (IPRIOR .EQ. 1) THEN
	         SUMW=0.0D0
	         ! AQBQ  = 0.0D0
	         ! AQBQ2 = 0.0D0 
	         
	         DO Q = 1,NQ
	            IAQ(Q) = IAQEMP(Q) / AQALL
	            SUMW  = SUMW   + IAQ(Q)
	            ! write(2,*)IT,q,IAQ(Q),IAQEMP(Q),aqall
	         END DO
	
	         DO Q = 1,NQ
	            IAQ(Q) = IAQ(Q)/SUMW
	            ! AQBQ  = AQBQ  + (IAQ(Q) * IBQ(Q))
	         END DO
	
	         ! AQMEAN    = AQBQ / DBLE(NQ)
	         ! DO Q = 1,NQ
	         !    write(2,*)it,q,IAQ(Q),IBQ(Q)
	         !    IBQ(Q)  = IBQ(Q) - AQBQ
	         !    AQBQ2  = AQBQ2 + (IAQ(Q)*IBQ(Q)*IBQ(Q))
	         !    write(2,*)IT,q,IAQ(Q),aqmean,AQBQ
	         ! END DO
	         ! AQSD  = DSQRT(AQBQ2)
	         ! AQBQ  = 0.0D0
	         ! AQBQ2 = 0.0D0 
	         ! DO Q = 1,NQ
	         !    IBQ(Q)  = IBQ(Q) / AQSD
	         !    AQBQ   = AQBQ  + (IAQ(Q) * IBQ(Q))
	         !    AQBQ2  = AQBQ2 + (IAQ(Q)*IBQ(Q)*IBQ(Q))
	         !    write(2,*)IT,q,IBQ(Q)
	         ! END DO
	         ! AQSD  = DSQRT(AQBQ2)
	         ! write(2,*)AQBQ,AQBQ2,AQSD
	         
	      ENDIF
	
	      ! CHECK IF CONVERGENCE HAS BEEN REACHED
	      ! IF NOT KEEP ITERATING 
	      ! IF IT HAS PROCEED TO STANDARD ERROR CALCULATIONS
	      ! AND PRINT OUT FINAL RESULTS
	
	      IGO = 0
	      IF (BIGCOR .GE. CONV) THEN
	         IGO = 1
	      ELSE
	         IFIN = IFIN + 2
	         IF (RIDGE .GT. 0.0) ITLAST = 1
	         IF (IFIN .EQ. 2)    ITLAST = 1
	         IF (IFIN .GE. 3 .AND. IRES .EQ. 1) ITLAST = 1
	         IF (IFIN .GE. 3 .AND. IRES .EQ. 0) IGO = 0
	      ENDIF
	      
	      ! Signal caller that there's more to do
	      MainLoop = 1
	   ELSE
	      
	      ! Signal caller that we're done iterations
	      MainLoop = 0
	      INITIALIZED = 2
	   ENDIF
	end function mainloop   
	
	! ************************************************
	!         **********************                         
	!         *  SUBROUTINE CLEANUP
	!         **********************                         
	!                                                                       
	!  This routine is called to output the results of 
	!  processing generated by successive calls to MainLoop.
	!  Once that is done, all global and any local arrays 
	!  that were allocated are freed in preparation for starting
	!  over again with new files.
	!                                                                       
	!  CALL CLEANUP()
	!                                                                       
	! ************************************************
	SUBROUTINE CLEANUP
	USE MIXLIB
	IMPLICIT NONE
	   
	   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	   ! LOCAL VARS ONLY - SHARED VARS ARE DECLARED
	   ! BELOW THE 'MODULE' STATEMENT AT THE TOP OF THE 
	   ! FILE
	   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	   
	   CHARACTER*4  BLANK(40), SIDEL(2), SIDET
	   CHARACTER*10 Blab
	   CHARACTER*2  CH(99)
	   
	   CHARACTER (LEN=8), ALLOCATABLE :: IPARMLAB(:)
	  
	   DATA SIDEL/' (1)',' (2)'/
	   DATA CH/'1','2','3','4','5','6','7','8','9','10', &
	   '11','12','13','14','15','16','17','18','19','20', &
	   '21','22','23','24','25','26','27','28','29','30', &
	   '31','32','33','34','35','36','37','38','39','40', &
	   '41','42','43','44','45','46','47','48','49','50', &
	   '51','52','53','54','55','56','57','58','59','60', &
	   '61','62','63','64','65','66','67','68','69','70', &
	   '71','72','73','74','75','76','77','78','79','80', &
	   '81','82','83','84','85','86','87','88','89','90', &
	   '91','92','93','94','95','96','97','98','99'/
	
	   INTEGER L, J, IR, IC, I, I2, IHI, ILO, IND, IND2, &
	         INP, INP0,J2,J3, J4, L1, NCONN, H,HR,Q
	   
	   REAL*8 PI, CORR, COV, SIDE, SIGSD, SIGVAR, VAR1, VAR2
	         
	   REAL (KIND = 8),POINTER ::ICS(:)
	   REAL (KIND=8), ALLOCATABLE :: IDER2SQ(:), &
	      ISE(:),IVAL(:),IZVAL(:),IPVAL(:), ICNTMAT(:), &
	      ICVAL(:),ICSE(:),ICZVAL(:),ICPVAL(:),ICDER2(:), IHVAL(:), IHDER2(:)
	  
	   IF(INITIALIZED .NE. 2) THEN
	      CALL POST_ERROR("Cleanup() routine called out of sequence in " // PROGNAME)
	      RETURN
	   ENDIF
	   
	   pi = 3.141592654D0
	   
	   ALLOCATE(ISE(NPAR))
	   ALLOCATE(IVAL(NPAR))
	   ALLOCATE(IZVAL(NPAR))
	   ALLOCATE(IPVAL(NPAR))
	   
	   IF (ICNT .GT. 0) THEN
	       ALLOCATE(ICNTMAT(NPAR*NPAR))
	       ALLOCATE(IHVAL(NPAR))
	       ALLOCATE(IHDER2(NPARR))
	   ENDIF
	   IF (NCON .GT. 0) THEN
	       ALLOCATE(ICVAL(NCON))
	       ALLOCATE(ICSE(NCON))
	       ALLOCATE(ICZVAL(NCON))
	       ALLOCATE(ICPVAL(NCON))
	       NCONN = (NCON * (NCON+1)) / 2
	       ALLOCATE(ICDER2(NCONN))
	   ENDIF
	
	   ALLOCATE(IPARMLAB(NPAR))
	
	    open(7,file="mixno.lik")
	    write(7,*) rlogl, npar
	    close(7)
	   WRITE(2,455)
	   455 FORMAT(///,1x,'---------------------------------------------------------', &
	         /,1x,'* Final Results - Maximum Marginal Likelihood Estimates *', &
	         /,1x,'---------------------------------------------------------',/)
	   
	   IF (R .EQ. 0) THEN
	      WRITE(2,563)IT,RIDGE,RLOGL,AIC,SBC,DEV,AICD,SBCD
	      563 FORMAT(1X,'Total Iterations  =',I4, &
	               /,1X,'Ridge             =',F8.3,//, &
	         1X,'Log Likelihood                 = ',F12.3,/, &
	         1X,"Akaike's Information Criterion = ",F12.3,/, &
	         1X,"Schwarz's Bayesian Criterion   = ",F12.3,//,&
	         1X,"==> multiplied by -2             ",      /  &
	         1X,'Log Likelihood                 = ',F12.3,/, &
	         1X,"Akaike's Information Criterion = ",F12.3,/, &
	         1X,"Schwarz's Bayesian Criterion   = ",F12.3,/)
	   ELSE
	      WRITE(2,564)IT,NQ1,RIDGE,RLOGL,AIC,SBC,DEV,AICD,SBCD
	      564 FORMAT(1X,'Total Iterations  =',I4, &
	               /,1X,'Quad Pts per Dim  =',I4, &
	               /,1X,'Ridge             =',F8.3,//, &
	         1X,'Log Likelihood                 = ',F12.3,/, &
	         1X,"Akaike's Information Criterion = ",F12.3,/, &
	         1X,"Schwarz's Bayesian Criterion   = ",F12.3,//,&
	         1X,"==> multiplied by -2             ",      /  &
	         1X,'Log Likelihood                 = ',F12.3,/, &
	         1X,"Akaike's Information Criterion = ",F12.3,/, &
	         1X,"Schwarz's Bayesian Criterion   = ",F12.3,/)
	   ENDIF
	   IF ((ICNT.EQ.0).OR.(ICNT.GT.0 .AND. (R.GT.0 .AND. IVCAT.EQ.0))) THEN
	      WRITE(2,567)
	      567 FORMAT(1X,'Response category contrasts = Reference-cell',/)   
	   ELSE IF ((ICNT.EQ.1) .AND. (R.EQ.0 .OR. (R.GT.0 .AND. IVCAT.EQ.1))) THEN
	      WRITE(2,568)
	      568 FORMAT(1X,'Response category contrasts = Adjacent-cell',/)   
	   ELSE IF ((ICNT.EQ.2) .AND. (R.EQ.0 .OR. (R.GT.0 .AND. IVCAT.EQ.1))) THEN
	      WRITE(2,569)
	      569 FORMAT(1X,'Response category contrasts = Helmert',/)   
	   ENDIF
	   ! Find the z-statistics for all terms  ZFIX and ZVAR
	   !  and their associated p-values        PFIX and PVAR
	
	   IR = 0
	   IF (R .GT. 0 .AND. IPRIOR .EQ. 0 .AND. NOMU .EQ. 0) THEN
	      DO J = 2,MAXJ
	         DO H = 1,R
	            IR = IR+1
	            IVAL(IR)  = IMU(IR)
	         END DO
	      END DO
	   ENDIF
	   
	   IF (P .GT. 0) THEN
	      L1 = 0
	      DO J = 2,MAXJ
	         DO L = 1,P
	            IR = IR + 1
	            L1 = L1 + 1
	            IVAL(IR) = IALPHA(L1)
	         END DO
	      END DO
	   ENDIF
	   
	   IF (R .GT. 0) THEN
	      L1 = 0
	      MJ = 2
	      IF (IVCAT .EQ. 1) MJ=MAXJ
	      DO J = 2,MJ
	         DO HR= 1,RR
	            ir  = ir + 1
	            L1  = L1 + 1
	            IVAL(IR) = ISIGMA(L1)
	         END DO
	      END DO
	   ENDIF
	!
	! see what kind of contrasts to get
	! ICNT = 1 : repeated contrasts
	! ICNT = 2 : helmert contrasts
	!
	   IF (ICNT .GT. 0 .AND. (R.GT.0 .AND. IVCAT .EQ. 0)) THEN
	      ICNT = 0
	      WRITE(2,565)
	      565 FORMAT(1X,'==> Non reference cell contrasts are not allowed', &
	               /,1x,'    unless the random-effect (co)variance(s) are', & 
	               /,1x,'    allowed to vary across the categories',        &
	               /,1x,'==> Reference cell contrasts are used instead')
	   ELSE IF (ICNT.GT.0 .AND. (R.EQ.0 .OR. (R.GT.0 .AND. IVCAT .EQ. 1))) THEN
	      CALL CONTRAST(ICNT,MAXJ-1,R,P,NOMU,IPRIOR,NPAR,ICNTMAT)
	      CALL MPYM(ICNTMAT,IVAL,IHVAL,NPAR,NPAR,0,0,1)
	      CALL GRAMM(ICNTMAT,IDER2,IHDER2,NPAR,NPAR,1,IZVAL)
	      CALL RELOC(IHVAL,IVAL,NPAR,1,0)
	      CALL RELOC(IHDER2,IDER2,NPAR,NPAR,1)
	      DEALLOCATE(IHVAL ,IHDER2)
	   END IF
	
	   IC = 0
	   DO I = 1,Npar
	      DO J = 1,I
	         IC = IC + 1
	         IF (I .EQ. J) THEN
	            ISE(I)   = DSQRT(IDER2(IC))
	            IZVAL(I) = IVAL(I) / ISE(I)
	            IPVAL(I) = 2.0D0 * (1.0D0 - PHIFN(DABS(IZVAL(I)),0))
	         ENDIF
	      END DO
	   END DO
	   
	   WRITE(2,57)
	   57 FORMAT(/,1X,'--------',5x,'------------',5X, &
	      '------------',5X,'------------',5X,'------------',/, &
	      1X,'Variable',5x,'    Estimate',5X,'Stand. Error',5x,    &
	      '           Z',5X,'     p-value',/,1X,'--------',5x, &
	      '------------',5X, '------------',5X,'------------', &
	      5X,'------------')
	
	   DO J=2,MAXJ
	
	      IF (ICNT .EQ. 0) THEN      
	         WRITE(2,681)ICODE(J),ICODE(1)
	      ELSE IF (ICNT .EQ. 1) THEN 
	         WRITE(2,681)ICODE(J),ICODE(J-1)
	      ELSE IF (ICNT .EQ. 2) THEN 
	         WRITE(2,6681)J-1, ICODE(J), ICODE(MAXJ), ICODE(J-1)
	      END IF
	      681 FORMAT(/,1X,'RESPONSE CODE',F4.0,' vs CODE',F4.0,/,1X, &
	         '-----------------------------',/,1x,'fixed effects')
	     6681 FORMAT(/,1X,'HELMERT CONTRAST',I3,': ', & 
	         '( code',F4.0,' to ','code',F4.0,' ) vs code',F4.0, &
	         /,1X,'--------------------',/,1x,'fixed effects')
	      IF (R .GT. 0 .AND. NOMU .EQ. 0) THEN
	         DO H = 1,R
	            INP = (J-2)*R + H
	            IF (IPRIOR .EQ. 0) THEN
	               WRITE(2,803)IBLABEL(H),IVAL(INP),ISE(INP), &
	                  IZVAL(INP),IPVAL(INP),SIDEL(2)
	            ELSE
	               WRITE(2,803)IBLABEL(H),IMU(INP)
	            ENDIF
	            803 FORMAT(1X,A8,3(5x,F12.5),F12.5,1X,A4)
	         END DO
	      ENDIF
	   
	      IF (P .GT. 0) THEN
	         DO L = 1,P
	            IF (NOMU .EQ. 0) THEN
	               INP  = (MAXJ-1)*R + (J-2)*P + L 
	            ELSE
	               INP  = (J-2)*P + L 
	            ENDIF
	            WRITE(2,804)IALABEL(L),IVAL(INP),ISE(INP), &
	                  IZVAL(INP),IPVAL(INP),SIDEL(2)
	            804 FORMAT(1X,a8,3(5x,F12.5),F12.5,1X,A4)
	         END DO
	      ENDIF
	
	      IF (r .gt. 0 .and. (IVCAT .EQ. 1 .OR.  &
	         (IVCAT .EQ. 0 .and. J.EQ.MAXJ))) then
	         IF (r .eq. 1) then
	            WRITE(2,581)
	            581 FORMAT(/,1X,'random effect variance term: expressed as a standard deviation')
	         ELSE 
	            WRITE(2,481)
	            481 FORMAT(/,1X,'random effect variance & covariance terms: cholesky of var-covariance matrix')
	         ENDIF
	         
	         IND = 1
	         IND2= 0
	         DO HR= 1,RR
	            IF (HR .EQ. IND) then
	               IF (IDIAG .EQ. 0 .AND. IVSEP .EQ. 0) THEN
	                  ind2 = ind2 + 1
	                  ind  = ind  + (ind2 + 1)
	               ELSE
	                  ind2 = ind2 + 1
	                  ind  = ind  + 1
	               ENDIF
	               
	               Blab = IBLABEL(IND2) // '  '
	               IF (IVSEP .EQ. 0) THEN
	                  side = 0.5d0
	                  SIDET= SIDEL(1)
	               ELSE
	                  side = 1.0d0
	                  SIDET= SIDEL(2)
	               ENDIF
	            ELSE
	               Blab = 'covariance'
	               side = 1.0d0
	               SIDET= SIDEL(2)
	            ENDIF
	            
	            IF (NOMU .EQ. 0 .AND. IVCAT .EQ. 0) THEN
	               INP  = (MAXJ-1)*(R+P) + HR 
	               INP0 = HR 
	            ELSEIF (NOMU .EQ. 1 .AND. IVCAT .EQ. 0) THEN
	               INP  = (MAXJ-1)*P +  HR 
	               INP0 = HR 
	            ELSEIF (NOMU .EQ. 0 .AND. IVCAT .EQ. 1) THEN
	               INP  = (MAXJ-1)*(R+P) + (J-2)*RR + HR 
	               INP0 = (J-2)*RR + HR 
	            ELSEIF (NOMU .EQ. 1 .AND. IVCAT .EQ. 1) THEN
	               INP  = (MAXJ-1)*P + (J-2)*RR + HR 
	               INP0 = (J-2)*RR + HR 
	            ENDIF
	            
	            !    IF (IPRIOR .EQ. 0) THEN
	            IPVAL(INP) = IPVAL(INP)*side
	            WRITE(2,805)Blab,IVAL(INP),ISE(INP), &
	                        IZVAL(INP),IPVAL(INP),SIDET
	            ! put transformed (or not) IVAL back into ISIGMA
	            ISIGMA(INP0) = IVAL(INP)
	            !    ELSE
	            !       WRITE(2,805)BLab,ISIGMA(HR)
	            !    ENDIF
	            805 FORMAT(1X,a10,3x,f12.5,2(5x,F12.5),F12.5,1X,A4)
	         END DO
	      ENDIF
	
	   END DO
	
	   WRITE(2,587)
	   587 FORMAT(//,1x,'note: (1) = 1-tailed p-value',/, &
	               '       (2) = 2-tailed p-value')
	
	   MJ = 2
	   IF (IVCAT .EQ. 1) MJ=MAXJ
	   
	   DO J=2,MJ
	      IF (R .EQ. 1 .AND. ICCY .EQ. 2 .AND. IUNIF .EQ. 0) THEN
	         IF (J.EQ.2) WRITE(2,811)
	         811 FORMAT(///,1X,'Calculation of the intracluster correlation' &
	             ,/,1x,'-------------------------------------------')
	         sigsd  = ISIGMA(J-1)
	         sigvar = ISIGMA(J-1)*ISIGMA(J-1)
	
	         IF (nfn .eq. 0) then
	            CORR   = sigvar / (sigvar + 1.0d0)
	            IF (J.EQ.2) WRITE(2,711)
	            711 FORMAT(1x,'residual variance = 1 (assumed)')
	            WRITE(2,611)j-1,sigsd,sigsd,sigvar,sigvar,sigvar,CORR
	            611 FORMAT(1x,i4,' cluster  variance = (',f5.3,' * ',f5.3, &
	                  ') = ',f6.3,/,6x,'intracluster correlation = ', &
	                  f6.3,' / (',f6.3,' + 1.000)',' = ',F5.3)
	
	         ELSEIF (nfn .eq. 1) then
	            CORR   = sigvar/(sigvar + (pi*pi/3.0D0))
	            IF (J.EQ.2) WRITE(2,712)
	            712 FORMAT(1x,'residual variance = pi*pi / 3 (assumed)')
	            WRITE(2,612)j-1,sigsd,sigsd,sigvar,sigvar,sigvar,CORR
	            612 FORMAT(1x,i4,' cluster  variance = (',f5.3,' * ',f5.3, &
	                  ') = ',f6.3,/,6x,'intracluster correlation = ', &
	                  f6.3,' / (',f6.3,' + (pi*pi/3))',' = ',F5.3)
	
	         ELSEIF (nfn .ge. 2) then
	            CORR   = sigvar/(sigvar + (pi*pi/6.0D0))
	            IF (J.EQ.2) WRITE(2,713)
	            713 FORMAT(1x,'residual variance = pi*pi / 6 (assumed)')
	            WRITE(2,613)j-1,sigsd,sigsd,sigvar,sigvar,sigvar,CORR
	            613 FORMAT(1x,i4,' cluster  variance = (',f5.3,' * ',f5.3, &
	                  ') = ',f6.3,/,6x,'intracluster correlation = ', &
	                  f6.3,' / (',f6.3,' + (pi*pi/6))',' = ',F5.3)
	         ENDIF
	
	      ELSEIF (R .EQ. 2 .AND. IDIAG .EQ. 0 .AND. IVSEP .EQ. 0) THEN
	         IF (J .EQ. 2) WRITE(2,812)
	         812 FORMAT(///,1X, &
	            'Calculation of the random effects variance-covariance matrix',/, &
	            1x,'------------------------------------------------------------')
	         ICS  =>ISIGMA((J-2)*RR+1:)
	         VAR1 = ICS(1)*ICS(1)
	         COV  = ICS(1)*ICS(2)
	         VAR2 = ICS(2)*ICS(2) + ICS(3)*ICS(3)
	         CORR = COV / (DSQRT(VAR1) * DSQRT(VAR2))
	         WRITE(2,813)J-1,IBLABEL(1),ICS(1),ICS(1),var1,ICS(1),ICS(2) &
	            ,cov,IBLABEL(2),ICS(2),ICS(2),ICS(3),ICS(3),var2,corr
	         813 FORMAT(1x,I4,1X,A8,' variance = (',f5.3,' * ',f5.3,') = ', &
	            f5.3,/,6x,7x,'covariance = (',f5.3,' * ',f5.3,') = ', &
	            f5.3,/,6x,a8,' variance = (',f5.3,' * ',f5.3,') + (', &
	            f5.3,' * ',f5.3,') = ',f5.3,//,6x, &
	            'Covariance expressed as a correlation = ',f5.3,//)
	      ENDIF
	   END DO
	
	   I = 0
	   IF (R .NE. 0 .AND. NOMU .NE. 1) THEN
	      DO J=2,MAXJ
	         DO H=1,R
	            I = I+1
	            ISE(I) = 1.0D0 / ISE(I)
	            IPARMLAB(I) = IBLABEL(H)
	         END DO
	      END DO
	   ENDIF
	
	   IF (P .NE. 0) THEN
	      DO J=2,MAXJ
	         DO L=1,P
	            I = I+1
	            ISE(I) = 1.0D0 / ISE(I)
	            IPARMLAB(I) = IALABEL(L)
	         END DO
	      END DO
	   ENDIF
	
	   IF (R .NE. 0) THEN
	      MJ = 2
	      IF (IVCAT .EQ. 1) MJ=MAXJ
	      DO J=2,MJ
	         DO H=1,RR
	            I = I+1
	            ISE(I) = 1.0D0 / ISE(I)
	            IPARMLAB(I) = 'VarCov' // CH(H)
	         END DO
	      END DO
	   ENDIF
	
	   BLANK(1:30) = '    '   ! array assignment gets all 30 elements
	   ND = 0
	
	   ! write out the parameter estimates  
	   ! the estimated variance-covariance estimates of these
	   ! estimates
	
	   DO I = 1,Npar
	      IZVAL(I) = IZVAL(I) / ISE(I)
	      WRITE(3,"(1X,a8,1x,F13.7)")IPARMLAB(I),IZVAL(I)
	   END DO
	   
	   ALLOCATE(IDER2SQ(NPAR*NPAR))
	   CALL CHAMS(IDER2,IDER2SQ,NPAR,1,0)
	   DO I = 1,Npar
	      ILO = I                  
	      IHI = NPAR*NPAR            
	      WRITE(4,"(1X,6F13.7)")(IDER2SQ(I2),I2=ILO,IHI,NPAR)
	   END DO
	   CLOSE(3)
	   CLOSE(4)
	
	   ! compute the contrasts of the parameter estimates
	   ! uses IZVAL and IPVAL as a work vector
	
	   IF (NCON .GE. 1) THEN
	      CALL MPYTR(IVAL,ICON,ICVAL,NPAR,1,0,NCON)
	      CALL GRMMT(ICON,IDER2,ICDER2,NPAR,NCON,1,IZVAL)
	      IC = 0
	      DO I = 1,NCON
	         DO J = 1,I
	            IC = IC + 1
	            IF (I .EQ. J) THEN
	               ICSE(I)   = DSQRT(ICDER2(IC))
	               ICZVAL(I) = ICVAL(I) / ICSE(I)
	               ICPVAL(I) = 2.0D0 *(1.0D0-PHIFN(DABS(ICZVAL(I)),0))
	            ENDIF
	         END DO
	      END DO
	   ENDIF
	
	   CALL MPDSD(ISE,IDER2,IDER2,NPAR)
	   CALL PRNT(2,IDER2,NPAR,NPAR,1,IPARMLAB,IPARMLAB,ND, &
	   BLANK,1,78,5,2,2, &
	   60HCorrelation of the Maximum Marginal Likelihood Estimates    )
	
	   ! print out the requested contrasts of the parameter estimates
	
	   IF (NCON .GE. 1) THEN
	
	      WRITE(2,655)
	      655 FORMAT(///,1x,'-------------------------------------',/, &
	                     1x,'* Transforms of parameter estimates *',/, &
	                     1x,'-------------------------------------')
	
	      CALL PRNT(2,ICON,NPAR,NCON,0,IPARMLAB,IXLAB,ND, &
	         BLANK,1,78,5,2,1, &
	         60HTranspose of the Transform Matrix (parameters by transforms))
	      WRITE(2,757)
	      757 FORMAT(/,1X,'Transform',5x,'    Estimate',5X,'Stand. Error',5x, &
	         '          Z',5X,'     p-value',/,1X,'---------',5x, &
	         '------------',5X,'------------',5X,'------------',5X, &
	         '------------')
	
	      DO L = 1,NCON
	         WRITE(2,840)L,ICVAL(L),ICSE(L),ICZVAL(L),ICPVAL(L)
	         840 FORMAT(1X,I8,4(5x,F12.5))
	         ICSE(L) = 1.0D0 / ICSE(L)
	      END DO
	
	      write(2,"(/,1x,'note: p-values are 2-tailed')")
	
	      CALL MPDSD(ICSE,ICDER2,ICDER2,NCON)
	      CALL PRNT(2,ICDER2,NCON,NCON,1,IXLAB,IXLAB,ND, &
	         BLANK,1,78,5,1,1, &
	         60HCorrelation of the MML Transformed Estimates                )
	
	   ENDIF
	
	
	   ! print out the empirical prior
	
	   IF (IPRIOR .EQ. 1) THEN
	
	      WRITE(2,656)
	      656 FORMAT(//,1x,'Empirical Prior Distribution',/,1x, &
	             '----------------------------',/)
	      IF (r .eq. 1) then
	         WRITE(2,657)
	         657 FORMAT(/,1X,' #',5X,'    Node',5x,'   Weight',/,1x,'--', &
	               6x,'-------',6x,'--------')
	         DO Q = 1,NQ
	            WRITE(2,"(1x,I2,5X,F8.5,5X,F8.5)") Q,IBQ(Q),IAQ(Q)
	         END DO
	      
	      ELSEIF (r .eq. 2) then 
	         WRITE(2,"(/,1X,'   Node',10F7.4,/)")(IBQ1(Q),Q=1,NQ1)
	         DO J  = 1,NQ1
	            J2 = (J-1)*NQ1 +  1
	            J3 = (J-1)*NQ1 + 10
	            WRITE(2,644)IBQ1(J),(IAQ(J4),J4=J2,J3)
	            WRITE(2,644)IBQ1(J),(IBQ(J4),J4=J2,J3)
	            644 FORMAT(1X,11F7.4)
	         END DO
	      ENDIF
	   ENDIF
	
	   CLOSE(2)
	   CLOSE(5)
	   
	   ! Deallocate local arrays now that we're finished
	   DEALLOCATE(ISE, IVAL ,IZVAL ,IPVAL ,IPARMLAB ,IDER2SQ)
	   
	   IF (ALLOCATED(ICVAL)) THEN
	      DEALLOCATE(ICVAL ,ICSE ,ICZVAL ,ICPVAL ,ICDER2)
	   ENDIF
	   
	   ! Now deallocate globals in prep for starting again
	   DEALLOCATE( ALLDAT, IDNI )
	   DEALLOCATE( IXLAB, IBLABEL, IALABEL )
	   DEALLOCATE( IWORKR,ISIGTAU,IAQEMP, ICODE, IBQ, IBQ1, IAQ,  &
	      IPRJ, INUMJ, IDER, IDER2, ICON, ISIGMA,IMU,IALPHA )
	
	   IF(ALLOCATED(ISIGMA0)) DEALLOCATE(ISIGMA0)
	   IF(ALLOCATED(IALPHA0)) DEALLOCATE(IALPHA0)
	   IF(ALLOCATED(IMU0)) DEALLOCATE(IMU0)
	   
	   INITIALIZED = 0
	
	END SUBROUTINE CLEANUP
	
	! *************************************************************
	!  SUBROUTINTE STARTV2 (MU1,SIGMA,R,MAXJ,CATFQ,          
	!              DENUM,WA,NFN,IPRIOR,IDIAG,IVSEP,IVCAT)     
	!                                                         
	!  GET STARTING VALUES FOR INTERCEPTS AND                 
	!  VARIANCE-COVARIANCE MATRIX OF RANDOM EFFECTS           
	!                                                         
	!  Parameters Sent                                        
	!  R      = Number of Random Effects                      
	!  MAXJ   = Number of Categories                          
	!  CATFQ  = MAXJ Category Frequencies                     
	!  DENUM  = Total Number of Level-1 observations          
	!  WA     = Predicted value of Y minus the intercept      
	!  NFN    = 0 for Probit                                  
	!         = 1 for Logistic                                
	!         = 2 for Complementary Log-Log                   
	!         = 3 for Log-Log                                 
	!                                                         
	!  IPRIOR = 0 for Specified Prior for Random Effects      
	!         = 1 for Empirical Prior for Random Effects      
	!  IDIAG  = 0 for correlated random effects               
	!         = 1 for independent random effects              
	!  IVSEP  = 0 R random effects aren't indicator variables 
	!         = 1 R random effects are indicator variables    
	!           (between-subjects indicator variables)        
	!  IVCAT  = 0 R random effects don't vary across MAXJ-1   
	!         = 1 R random effects do vary across MAXJ-1      
	!             category contrasts (simple contrasts to     
	!             the first category)                         
	!                                                         
	!  Parameters Returned                                    
	!  MU1    = Starting value for intercepts (maxj-1)        
	!  SIGMA  = (R*(R+1))/2 * [(maxj-1) or 1]                 
	!           Starting values for Variance                  
	!           Covariance matrix of random effects           
	!                                                         
	! *************************************************************
	
	SUBROUTINE STARTV2(MU1,SIGMA,R,MAXJ,CATFQ,DENUM,WA, &
	                   NFN,IPRIOR,IDIAG,IVSEP,IVCAT)
	   IMPLICIT REAL*8(A-H,O-Z),INTEGER*4(I-N)
	   DOUBLE PRECISION MU1,SIGMA,DENUM,WA,CATFQ
	   INTEGER R
	   DIMENSION SIGMA(1),MU1(1),CATFQ(1)
	   DATA PI/3.141592654D0/
	
	   cumsum = 0.0D0
	   DO j=1,maxj-1
	      cumsum   = cumsum + CATFQ(J)
	      cumprob  = cumsum / DENUM
	      cumodds  = cumsum / (DENUM - cumsum)
	      Rlncumod = DLOG(cumodds)
	      IF (nfn .eq. 0) then
	          tempparm = 0.625 * (Rlncumod + wa)
	      ELSEIF (nfn .eq. 1) then
	          tempparm = Rlncumod + wa
	      ELSEIF (nfn .eq. 2) then
	          tempparm = DLOG(0.0D0 - DLOG(1.0d0 - cumprob)) + wa
	      ELSEIF (nfn .eq. 3) then
	          tempparm = DLOG(0.0D0 - DLOG(cumprob)) + wa
	      ELSE
	      ENDIF
	      IF (j .EQ. 1) THEN
	         IF (IPRIOR .EQ. 0) THEN
	             MU1(1) = tempparm
	         ENDIF
	      ELSE
	         MU1(j) = tempparm 
	      ENDIF
	   END DO
	
	   IR = 0
	   MJ = 1
	   IF (IVCAT .EQ. 1) MJ = MAXJ-1
	   DO j = 1,MJ
	      IF (R .EQ. 1) THEN
	         IR = IR+1
	         SIGMA(IR) = .31623d0 + (J-1)*.1*.31623D0
	         IF (nfn .EQ. 1) & 
	            SIGMA(IR) = SIGMA(IR) * DSQRT(pi*pi/3.0D0)
	         IF (nfn .GE. 2) &
	            SIGMA(IR) =  SIGMA(IR) * DSQRT(pi*pi/6.0D0)
	      
	      ELSEIF (R .GT. 1 .AND. IDIAG .EQ. 0 .AND. IVSEP .EQ. 0) THEN
	         DO I  = 1,R
	            DO I2 = 1,I
	               IR = IR + 1
	               IF (I2 .EQ. I) THEN
	                  IF (I .EQ. 1) THEN
	                     SIGMA(IR) = 1.0d0 + (J-1)*.1
	                  ELSE
	                     SIGMA(IR) = 0.5d0 + (J-1)*.1
	                  ENDIF
	               ELSE
	                       SIGMA(IR) = 0.0d0
	               ENDIF
	               IF (nfn .EQ. 1) &
	                  SIGMA(IR) = SIGMA(IR)*DSQRT(pi*pi/3.0D0)
	               IF (nfn .GE. 2) &
	                  SIGMA(IR) = SIGMA(IR)*DSQRT(pi*pi/6.0D0)
	            END DO
	         END DO
	      
	      ELSEIF (R .GT. 1 .AND. (IDIAG .EQ. 1 .OR. IVSEP .EQ. 1)) THEN
	         DO I  = 1,R
	            IR = IR + 1
	            IF (I .EQ. 1) THEN
	               SIGMA(IR) = 1.0d0 + (J-1)*.1
	            ELSE
	               SIGMA(IR) = 0.5d0 + (J-1)*.1
	            ENDIF
	            IF (nfn .EQ. 1) SIGMA(IR) = SIGMA(IR)*DSQRT(pi*pi/3.0D0)
	            IF (nfn .GE. 2) SIGMA(IR) = SIGMA(IR)*DSQRT(pi*pi/6.0D0)
	         END DO
	      ENDIF
	   END DO
	
	   RETURN
	END SUBROUTINE STARTV2

END program RRM_MODULE