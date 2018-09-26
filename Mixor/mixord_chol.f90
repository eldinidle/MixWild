!  **************************************************************
!  MIXOR        
!
!  modified on 
!  10/09/02 to correct AIC & BIC statistics
!
!  modified on 8/5/04 to correct BIC statistics - use N not NTOT
!
!  10/24/02 to put in scaling terms as option for KG
!
!
! *************************************************************
!  ordinal RRM MAIN PROGRAM                               
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
! CODEX dimension = maxxj = max categories for crosstab X variable
! IDNI            = 2 * maxn (number of level-2 units)
! IDNI is allocated in the data read routine.
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
!
! this program was modified on 10/2/95 
! a) to put in right-censoring for survival analysis
!
! this program was modified on 3/23/96 
! a) to reverse the sign of the covariate coefficients
!
! this program was modified on 3/26/96 
! a) to put in interactions with the thresholds
! (assumes the first KG of the p fixed effects interact with the thresholds)
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
! a) write out estimates to MIXOR.EST (with labels)
! b) write out asymptotic variance-covariance matrix (in unpacked form)
!  to MIXOR.VAR (without labels - in the order specified in MIXOR.EST)
!
! this version was updated 10/9/96 to include
!
! IVSEP  =  0 R random effects don't indicate separate groups of subjects
!           1 R random effects do    indicate separate groups of subjects
!            (or in IRT case, X variables for R random effects change
!             but there is only 1 random effect beta)
! 5/3/99
! modified code to estimate empirical prior weights (IPRIOR=1)
!
! 2013
! program modified to calculate sigma (random effect matrix) from estimated cholesky
! adaptive quadrature option added
! option added to display either sigma or cholesky
! .dat .def .out filenames displayed in output

MODULE RRM_MODULE ! must have same name in all 4 programs
IMPLICIT NONE

   INTEGER ,POINTER ::IDNI(:)
   INTEGER H, H2, HR, ICCY, ICEN, IDIAG, IER, IFIN, IGO, &
      IPRIOR, IRBAD, IRES, IRT, IRTT, ISIG, IT, ITLAST, &
      IUN, IUNIF, IVSEP, IWT, KG, MAXJ, N, NCON, ND, NFN, &
      NGAM, NGAM1, NOMU, NONPOS, NPAR, NPARR, NPR, NQ, &
      NQ1, NRP1, P, Q, R, RR, KS, NPROP, AQUAD, CHOL, maxits, fail

   INTEGER :: INITIALIZED = 0

   REAL (KIND=8) CONV, PI, RADD, RIDGE, RLOGL, RLOGLP, AIC, SBC, &
                 DEV, AICD, SBCD, RNTOT, RNPAR, DZT, RIDGEMAX,sum2,sbcn
   REAL (KIND = 8),ALLOCATABLE,TARGET ::IALPHA(:)
   REAL (KIND = 8),POINTER ::ALLDAT(:)
   REAL (KIND=8), ALLOCATABLE :: IALPHA0(:), IAQ0(:), iaq(:), IAQ1(:), IAQ2(:), &
      IAQEMP(:), IBQ0(:), ibq(:), IBQ1(:), ICODE(:), ICON(:), IDER(:), &
      IDER2(:), IGAM(:), IGAM0(:), IMU(:), IMU0(:), ISIGMA(:), &
      ISIGMA0(:), ISIGTAU(:), ITHETA(:), ITHETAV(:), IWG(:), &
      IWORKR(:), ITAU(:), ITAU0(:), IDER2S(:), ithetas(:,:), ithetavs(:,:)
      
   CHARACTER (LEN=16), ALLOCATABLE :: IBLABEL(:),IALABEL(:)
   CHARACTER (LEN=4), ALLOCATABLE :: IXLAB(:)
   CHARACTER*4  HEAD(30)
   CHARACTER, PARAMETER:: PROGNAME*5 = 'MIXOR'


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
SUBROUTINE INIT(FILENAMES)
USE MIXLIB
IMPLICIT NONE

   CHARACTER*40 , INTENT(IN OUT):: FILENAMES(4)
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! LOCAL VARS ONLY - SHARED VARS ARE DECLARED
   ! BELOW THE 'MODULE' STATEMENT AT THE TOP OF THE 
   ! FILE
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   CHARACTER*16  YLabel,XLabel,TLABEL
   CHARACTER*80 FILEDAT,FILEOUT,filedef
   CHARACTER*160 LABTEMP,TEMPFORM
   
   INTEGER DEATHI, I, IADD, IC, IC2, IC3, IC4, ICAT, ICATYX, &
      ICW, ICX, ICXW, IDIND, II, IR, ISTART, J, J2, J3, &
      K, K2, K3, K4, KI, KIND, L, L2, MAXCOL, MAXJXJ, MAXK, &
      MAXXJ, MISS, MOD, NALL, NCATYX, NCPAR, NQR, NQRR, NTOT, &
      NXLAB, P1, P1P1, WTIND, XCIND, YIND, kg0
   INTEGER ,ALLOCATABLE :: IXIND(:), IWIND(:)
   
   REAL (KIND=8) CATADD, DENUM, DET, TEMPSUM, UNO, WA, WSUM, &
      YDEV, YMISS

   REAL (KIND=8), ALLOCATABLE :: CODEX(:), ICATFQ(:), ICATFQX(:), &
      IWMISS(:), IWORKCAT(:), IWORKP1(:), IWP1P1(:), IWRKP1(:), IXMISS(:)
      
   LOGICAL KEEP
   
   ! Start by checking if we are ready to run yet.  You 
   ! can't call init twice in a row, and must call 
   ! \ before you can start another session.
   
   IF(INITIALIZED .NE. 0) THEN
      CALL POST_ERROR("Init() routine called out of sequence in " // PROGNAME)
      RETURN
   ENDIF
   
   ! ********************************************************
   ! OPEN AND READ THE DEFINITION (.DEF) FILE
   ! ********************************************************
   ! READ IN TWO LINES (60 CHARS EACH) FOR THE TITLE 
   ! AND THEN PARAMETERS AND STARTING VALUES FROM MIXOR.DEF

!   OPEN(1,ACTION='READ,DENYWRITE', FILE=FILENAMES(1))
   OPEN(1,ACTION='READ', FILE=FILENAMES(1))
   READ(1,"(15A4)") HEAD
   READ(1,"(A80)")FILEDAT
   READ(1,"(A80)")FILEOUT
   READ(1,"(A80)")filedef
   IF(LEN_TRIM(FILENAMES(2)) .EQ. 0) FILENAMES(2) = FILEDAT 
   IF(LEN_TRIM(FILENAMES(3)) .EQ. 0) FILENAMES(3) = FILEOUT 
   
   OPEN(2, FILE= FILENAMES(3))  ! OUT file

   OPEN(3, FILE= TRIM(ADJUSTL(FILENAMES(4))) // '.EST')
   OPEN(4, FILE= TRIM(ADJUSTL(FILENAMES(4))) // '.VAR')
   OPEN(5, FILE= TRIM(ADJUSTL(FILENAMES(4))) // '.lik')

   READ(1,*) NPR,MAXCOL,R,P,CONV,MAXJ,MISS,ISTART,IWT, &
             ICATYX,NQ1,AQUAD,NFN,ICEN,KG,IADD,NCON,IDIAG,NOMU,IVSEP,NPROP,CHOL,maxits

   IF (NQ1 .LE. 1) THEN
      BACKSPACE 1
      READ(1,*) NPR,MAXCOL,R,P,CONV,MAXJ,MISS,ISTART, &
             IWT,ICATYX,IPRIOR, &
             IUNIF,NQ1,AQUAD,NFN,ICEN,KG,IADD,NCON,IDIAG,NOMU,IVSEP,NPROP,CHOL,maxits
   ENDIF
    
   IF (IADD .NE. 1 .AND. IADD .NE. -1) IADD = -1

   ! IRES = 0 DON'T PRINT OUT INDIVIDUAL PARAMETERS          
   !      = 1 PRINT INDIVIDUAL PARAMETERS TO FILERES         

   IRES=0
   IF(R.GE.1) IRES=1
   IF (IRES .EQ. 1) THEN
      ! Use either the module name, or the last name on the cmdline
      OPEN(9, FILE= TRIM(ADJUSTL(FILENAMES(4))) // '.RES')
   ENDIF          

   ! ALLOCATE  R 

   ALLOCATE (IXMISS(R))
   ALLOCATE (IXIND(R))
   ALLOCATE (IWORKR(R))
   ALLOCATE(ISIGTAU(R))

   ! IRT = dimension of theta

    IRT = R
    if(ivsep .eq. 1) irt = 1
   ALLOCATE(ITHETA(IRT))

   ! ALLOCATE  rr = (r * (r+1) / 2) 

   RR = (R * (R+1)) / 2
   if(ivsep .eq. 1 .or. idiag .eq. 1) RR = R
   ALLOCATE(ISIGMA(RR))
   
   IRTT = (IRT*(IRT+1))/2
   ALLOCATE(ITHETAV(IRTT))

   ! ALLOCATE  nquad points per dimension (nq1)

   ALLOCATE(IBQ1(NQ1))

   ! ALLOCATE  nq1**r

   NQR = NQ1**IRT

   ALLOCATE(IAQ0(NQR))
   ALLOCATE(IAQ(NQR))
   ALLOCATE(IAQEMP(NQR))

   ! ALLOCATE  R * (nq1**r)

   NQRR = IRT*NQR
   ALLOCATE(IBQ0(NQRR))
   ALLOCATE(IBQ(NQRR))

   ! ALLOCATE  these vectors consectutively 
   !           so that random-effects can 
   !           be reassigned as fixed-effects
   !           if numerical problems develop 

   ALLOCATE(IMU(MAX(1,R)))   ! Lahey doesn't like 0 length allocs
   ALLOCATE(IALPHA(P))
   ALLOCATE(IBLABEL(MAX(1,R)))
   ALLOCATE(IALABEL(P))

   ! ALLOCATE  p 

   ALLOCATE(IWMISS(P))
   ALLOCATE(IWIND(P))

   ! ALLOCATE  maxj

   ALLOCATE(ICODE(MAXJ))

   ! ALLOCATE    ngam                       
   !    INCREASE NGAM BY 1 IF THERE IS CENSORING
kg0 = kg
   IF (ICEN.EQ.1) THEN
      NGAM1 = MAXJ-1
   ELSE 
      NGAM1 = MAXJ-2
   ENDIF

   IF (NPROP.EQ.1) THEN
      KS = 0
      KG = KG
   ELSE
      KS = KG
      KG = 0
   ENDIF
   NGAM = NGAM1*(KG+1)

   ALLOCATE(IWG(NGAM1))
   ALLOCATE(IGAM(NGAM))

   ! ALLOCATE KS for scaling terms 

   ALLOCATE(ITAU(KS))

   ! ALLOCATE    (npar * (npar+1)) / 2
   !    where    npar = r + rr + p + ngam + ks

   npar = p+r+rr+ngam+ks
   IF (NOMU .EQ. 1) npar = p+rr+ngam+ks

   RNPAR = DBLE(NPAR)

   ALLOCATE(IDER(NPAR))

   NCPAR = NPAR*NCON
   NPARR = (NPAR * (NPAR+1)) / 2
   ALLOCATE(ICON(NCPAR))
   ALLOCATE(IDER2(NPARR))
   ALLOCATE(IDER2S(NPAR*NPAR))

   ! MAXCOL is the dimension of the data matrix to read in     
   ! IDIND (integer) indicates which column contains the level 2 ID
   ! YIND (integer) indicates which column contains the dependent var
   ! R is the dimension of the X matrix to be used in the analysis
   ! IXIND() (integer) indicates which columns are to be used for X           
   ! P is the dimension of the W matrix to be used in the analysis
   ! IWIND() (integer) indicates which columns are to be used for W           
   ! WTIND (integer) indicates which column to be used for the 
   !       weighting of level-2 units           
   !   
   ! DEATHI (integer) indicates which column to be used for the 
   !       event indicator (1=event & 0=censored) IF ICEN=1
   !   
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
   ! ICEN   = 0 no right-censoring (or censored observations are coded
   !                                EXPLICITLY as last category)
   !          1 right-censoring
   !
   ! KG     = number of covariates to interact with threshold parameters
   !          (the first KG variables in W)                        
   !
   ! IADD   =  1      add covariates & mean of random effects to   thresholds 
   !        = -1 subtract covariates & mean of random effects from thresholds 
   !
   ! NCON   =  number of transforms of estimated parameters 
   !           (linear re-expressions)                          
   !
   !    IDIAG  =  0 correlated random effects                       
   !          1 independent random effects                              
   !
   !    NOMU   =  0 estimate the mean of the random effects         
   !          1 don't estimate the mean of the random effects 
   !            (assumed to be zero)
   !
   !    IVSEP  =  0 R random effects don't indicate separate groups of subjects
   !          1 R random effects do    indicate separate groups of subjects
   !
   !
   !  IFIN = 0 ESTIMATE TAU (ln of diagonaol elements of cholesky)
   !       = 1 ESTIMATE SIGMA (cholesky of variance-covariance
   !                           matrix of random effects)
   !
   !   CHOL = 0 Display sigma and its var/covar
   !        = 1 Display estimated Cholesky
   !
   
   READ(1,*) IDIND,YIND

   IF (R .GE. 1) READ(1,*) (IXIND(H), H=1,R)
   IF (P .GE. 1) READ(1,*) (IWIND(L), L=1,P)
   IF (IWT .EQ. 1) READ(1,*) WTIND
   IF (ICEN.EQ. 1) READ(1,*) DEATHI

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
         ELSE
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

   IF (MISS .EQ. 1 .or. miss .eq. 2) THEN
      READ(1,*) YMISS
      IF (R .GE. 1) READ(1,*) (IXMISS(H), H=1,R)
      IF (P .GE. 1) READ(1,*) (IWMISS(L), L=1,P)
   ENDIF

   READ(1,*) YLabel
   IF (YLABEL .EQ. '        ') READ(1,8) YLABEL

   IF (R .NE. 0) THEN
      READ(1,*) TLABEL
      IF (TLABEL .EQ. '        ') THEN
         READ(1,*) (IBLABEL(H), H=1,R)
      ELSE
         BACKSPACE 1
         READ(1,*) (IBLABEL(H), H=1,R)
      ENDIF
      IF (NOMU .EQ. 0 .AND. ISTART .EQ. 1) READ(1,*)(IMU(H), H=1,R)
   ENDIF

   IF (P .NE. 0) THEN
!      READ(1,*) TLABEL
!      IF (TLABEL .EQ. '        ') THEN
         READ(1,*) (IALABEL(L), L=1,P)
!      ELSE
!         BACKSPACE 1
!         READ(1,*) (IALABEL(L), L=1,P)
!      ENDIF
      IF (ISTART .EQ. 1) READ(1,*) (IALPHA(L), L=1,P)
      8 FORMAT(10A16)

      IF (ICATYX .EQ. 1) THEN
         IF (ICXW .EQ. 1) XLabel = IBLABEL(ICX)
         IF (ICXW .EQ. 2) XLabel = IALABEL(ICX)
      ELSE
      ENDIF
   ENDIF

   IF(ISTART .EQ. 1) THEN
      IF (R .NE. 0) READ(1,*) (ISIGMA(HR), HR=1,RR)
      IF (NGAM .NE. 0) READ(1,*) (IGAM(L), L=1,ngam)
      IF (KS .NE. 0) READ(1,*) (ITAU(L), L=1,ks)
   ENDIF

   IF (NCON .NE. 0) READ(1,*) (ICON(L), L=1,NCPAR)

   CLOSE(1)   ! Finished reading the DEF file

OPEN(1,FILE=FILEDEF)
    WRITE(1,"(15A4)") HEAD
    WRITE(1,"(A80)")FILEDAT
    WRITE(1,"(A80)")FILEOUT
    WRITE(1,"(A80)")FILEDEF
    WRITE(1,*) NPR,MAXCOL,R,P,CONV,MAXJ,MISS,ISTART,IWT,ICATYX,IPRIOR, &
             IUNIF,NQ1,AQUAD,NFN,ICEN,KG0,IADD,NCON,IDIAG,NOMU,IVSEP,NPROP,CHOL,maxits
             
    WRITE(1,*) IDIND, YIND
    IF (R .GE. 1) THEN
        WRITE(1,*) (IXIND(I), I=1,R)
    END IF
    IF (P .GE. 1) THEN
        WRITE(1,*) (IWIND(I), I=1,P)
    END IF
   IF (IWT .EQ. 1) WRITE(1,*) WTIND
   IF (ICEN.EQ. 1) WRITE(1,*) DEATHI

   WRITE(1,*)(ICODE(J), J = 1,MAXJ)
      
   IF (ICATYX .EQ. 1) THEN
      WRITE(1,*) XCIND,MAXXJ,(CODEX(J),J=1,MAXXJ)
   ENDIF

   IF (MISS .EQ. 1 .or. miss .eq. 2) THEN
      WRITE(1,*) YMISS
      IF (R .GE. 1) WRITE(1,*) (IXMISS(H), H=1,R)
      IF (P .GE. 1) WRITE(1,*) (IWMISS(L), L=1,P)
   ENDIF

   WRITE(1,8) YLabel
   
   IF (R .NE. 0) THEN
      WRITE(1,8) (IBLABEL(H), H=1,R)
      IF (NOMU .EQ. 0 .AND. ISTART .EQ. 1) WRITE(1,*)(IMU(H), H=1,R)
   ENDIF

   IF (P .NE. 0) THEN
      WRITE(1,8) (IALABEL(L), L=1,P)
      IF (ISTART .EQ. 1) WRITE(1,*) (IALPHA(L), L=1,P)
   ENDIF

   IF(ISTART .EQ. 1) THEN
      IF (R .NE. 0) WRITE(1,*) (ISIGMA(HR), HR=1,RR)
      IF (NGAM .NE. 0) WRITE(1,*) (IGAM(L), L=1,ngam)
      IF (KS .NE. 0) WRITE(1,*) (ITAU(L), L=1,ks)
   ENDIF

   IF (NCON .NE. 0) WRITE(1,*) (ICON(L), L=1,NCPAR)
    CLOSE(1)
   ! parameters         

   pi = 3.141592654D0

   NRP1 = 1 + R + P
   IF (ICEN .EQ. 1) NRP1 = NRP1 + 1
   UNO  = 1.0D0

   ! get quadrature nodes & weights  

   ALLOCATE(IAQ1(NQ1))
   ALLOCATE(IAQ2(NQRR))

   CALL QUADP(IBQ0,IBQ1,IAQ0,NQ1,NQ,IRT,IUNIF,IAQ1,IAQ2)

   ! ********************************************************
   ! GET THE DATA 
   ! ********************************************************
   ! N    (total number of level 2 units)
   ! NTOT (total number of observations)
   ! 

   CALL READAT(FILEDAT,N,NTOT,MAXK,MAXCOL,R,P,ALLDAT, &
         IDNI,IDIND,YIND,IXIND,IWIND, &
         MISS,YMISS,IXMISS,IWMISS, &
         IWT,WTIND,ICEN,DEATHI)

   ! NALL is the number of elements for the data read in from filedat 
   ! that ALLDAT() contains 

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
      
      IF(KEEP) THEN
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
            ALLDAT(IC3) = ALLDAT(IC3) + ALLDAT(K)
         ELSE
            ICW = (K + (IDNI(II*2)*NRP1) -IR+1) - ((KI-2)*NRP1)
            ALLDAT(IC3) = ALLDAT(IC3) + ALLDAT(ICW)*ALLDAT(K)
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
               ENDIF
            END DO
         ENDIF

         ! for one of the Xs or Ws get the crosstab by Y
         ! and add the weight if appropriate

         IF (ICATYX .EQ. 1 .AND. IR .EQ. NCATYX) THEN
            DO K4 = 1,MAXXJ
               IF (ALLDAT(K) .GE. (CODEX(K4)-0.0001D0)  &
                  .AND. ALLDAT(K) .LE. (CODEX(K4)+0.0001D0)) THEN
                   KIND = ((ICAT-1)*MAXXJ) + K4
                   ICATFQX(KIND) = ICATFQX(KIND) + CATADD
               ELSE
               ENDIF
            END DO
         ELSE
         ENDIF
      ENDIF
   END DO

   !    calculate the means

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

      IF (K .EQ. 2 .AND. R .EQ. 1 .AND. ALLDAT(IC3)  &
         .GT. 0.999D0 .AND. ALLDAT(IC3) .LT. 1.001D0) &
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
            ALLDAT(IC4) = ALLDAT(IC4) + ((ALLDAT(K) - ALLDAT(IC3))**2)
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
                     IWRKP1(L) = IWRKP1(L) + (YDEV * IWORKP1(L))
                  END DO
               ENDIF
            ENDIF
         ELSE
            ICW = (K + (IDNI(II*2) * NRP1) - IR+1) - ((KI-2)*NRP1)
            ALLDAT(IC4) = ALLDAT(IC4) + ALLDAT(ICW) * &
                   ((ALLDAT(K) - ALLDAT(IC3))**2)
            
            IF (IR .EQ. 1) YDEV = ALLDAT(ICW)*(ALLDAT(K) - ALLDAT(IC3))

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
            ENDIF
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

      IF (K .EQ. 2 .AND. R .EQ. 1  &
         .AND. ALLDAT(IC4) .GT. -0.001D0  &
         .AND. ALLDAT(IC4) .LT. 0.001D0)  &
          ICCY = ICCY + 1
   END DO

   ! write out descriptives

   WRITE(2,554) 
   554 FORMAT(1x,'MIXOR - The program for mixed-effects ordinal regression analysis', &
             /,1x,'         (version 2)',/)

   WRITE(2,"(1x,15A4)") HEAD

   WRITE(2,"(/,1X,4(1x,A40))") (FILENAMES(2))
   WRITE(2,"(1X,4(1x,A40))") (FILENAMES(3))
   WRITE(2,"(1X,4(1x,A40))") (filedef)

   IF (nfn .eq. 0) THEN
      WRITE(2,"(/,1x,'Response function: normal',/)") 
   ELSEIF (nfn .eq. 1) THEN
      WRITE(2,"(/,1x,'Response function: logistic',/)") 
   ELSEIF (nfn .eq. 2) THEN
      WRITE(2,"(/,1x,'Response function: complementary log-log',/)") 
   ELSEIF (nfn .eq. 3) THEN
      WRITE(2,"(/,1x,'Response function: log-log',/)") 
   ENDIF

   IF (r .ge. 1 .and. iunif .eq. 0 .and. iprior .eq. 0) THEN
      WRITE(2,"(1x,'Random-effects distribution: normal',/)") 
   ELSEIF (r .ge. 1 .and. iunif .eq. 1 .and. iprior .eq. 0) THEN
      WRITE(2,"(1x,'Random-effects distribution: rectangular',/)") 
   ELSEIF (r .ge. 1 .and. iprior .eq. 1) THEN
      WRITE(2,"(1x,'Random-effects distribution: empirical',/)")
   ENDIF

   IF (IADD .EQ. -1) THEN
      RADD = -1.0D0
      WRITE(2,5574) 
      5574 FORMAT(1x,'Covariate(s) and random-effect(s) mean subtracted from thresholds', &
         /,1x,'==> positive coefficient = positive association between regressor', &
         /,1x,'    and ordinal outcome',/)
   ELSEIF (IADD .EQ. 1) THEN
      RADD =  1.0D0
      WRITE(2,5575) 
      5575 FORMAT(1x,'Covariate(s) and random-effect(s) mean added to thresholds', &
         /,1x,'==> positive coefficient = negative association between regressor', &
         /,1x,'    and ordinal outcome',/)
   ENDIF

   WRITE(2,255)
   255 FORMAT(//,1x,'Numbers of observations',/,1x,'-----------------------',/)
   
   IF (IWT .EQ. 0) THEN
      WRITE(2,"(1x,'Level 1 observations = ',i6)") NTOT
      IF (R .GE. 1) THEN
         WRITE(2,"(1x,'Level 2 observations = ',i6)") N
         WRITE(2,5610)
         5610 FORMAT(//,1x,'The number of level 1 observations per level 2 unit are:',/)
         WRITE(2,"(1x,19I4)")(IDNI(I),I=2,N*2,2)
      ENDIF     
   ELSE
      WRITE(2,506) WSUM,NTOT
      506 FORMAT(1x,'Level 1 observations = ',f8.2, &
               /,1x,'Level 1 patterns =     ',I8)
      IF (R .GE. 1) THEN
         WRITE(2,516)SUM2,N
         516 FORMAT(/,1x,'Level 2 observations = ',f8.2, &
                   /, 1x,'Level 2 patterns =     ',I8)
         WRITE(2,5070)
         5070 FORMAT(//,1x,'The number of level 1 patterns per level 2 pattern are:',/)
         WRITE(2,"(1x,19I4)") (IDNI(I),I=2,N*2,2)

      ENDIF     
   ENDIF

   WRITE(2,257)
   257 FORMAT(//,1x,'Descriptive statistics for all variables', &
               /,1x,'----------------------------------------',/)
   WRITE(2,357)
   357 FORMAT(1x,'Variable',5X,'             Minimum',5x,'     Maximum', &
           5x,'        Mean',5x,' Stand. Dev.',/)
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
         WRITE(2,377)IBLABEL(H),ALLDAT(IC),ALLDAT(IC2), &
                  ALLDAT(IC3),ALLDAT(IC4)
      END DO
   ENDIF
   
   IF (P .GE. 1) THEN
      DO L = 1,P
         IC  = IC  + 1
         IC2 = IC2 + 1
         IC3 = IC3 + 1
         IC4 = IC4 + 1
         WRITE(2,377)IALABEL(L),ALLDAT(IC), ALLDAT(IC2), &
                  ALLDAT(IC3),ALLDAT(IC4)
      END DO
   ENDIF
   
   IF (ICEN .GE. 1) THEN
      IC  = IC  + 1
      IC2 = IC2 + 1
      IC3 = IC3 + 1
      IC4 = IC4 + 1
      WRITE(2,377)'Event   ',ALLDAT(IC),ALLDAT(IC2) &
                  ,ALLDAT(IC3),ALLDAT(IC4)
   ENDIF

   377 FORMAT(1x,A16,4(5X,F12.5))

   WRITE(2,258) YLabel
   258 FORMAT(//,1x,'Categories of the response variable ',A8, &
               /,1x,'--------------------------------------------',/)
   WRITE(2,457)
   457 FORMAT(1X,'Category',5X,'   Frequency',5x,'  Proportion',/)

   DO J = 1,MAXJ
      IF (IWT .EQ. 0) THEN
         DENUM = DBLE(NTOT)
      ELSE
         DENUM = WSUM
      ENDIF
      IWORKCAT(J) = ICATFQ(J) / DENUM
      WRITE(2,"(1X,F8.2,5X,F12.2,5x,f12.5)") ICODE(J),ICATFQ(J),IWORKCAT(J)
   END DO

   IF (ICATYX .EQ. 1) THEN
      WRITE(2,358)XLabel,YLabel
      358 FORMAT(//,1x,'Crosstabulation of variable ',A8, &
                 ' by the response variable ',A8, &
                 /,1x,'----------------------------------------------------------------------',/)
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
         
         WRITE(2,"(/,1X,F8.2,8F8.1)")  &
            CODEX(J),(ICATFQX(J3),J3=J,MAXJXJ,MAXXJ),TEMPSUM
         WRITE(2,TEMPFORM)(IWORKCAT(J3),J3=1,MAXJ)
      END DO
      WRITE(2,"(/,1X,'Total   ',8F8.1)")(ICATFQ(J),J=1,MAXJ),DENUM
   ELSE
   ENDIF

   ! done writing out descriptives, get starting values

   IF (ISTART .NE. 1) THEN

      ! calculate the starting values for the regression coefficients
      IF (P1 .GT. 0) THEN
        CALL INVS(IWP1P1,P1,DET,IWORKP1,IER,.FALSE.)
        CALL MPYM(IWP1P1,IWRKP1,IWORKP1,P1,P1,1,0,1)

        WA = 0.0D0
        IC3 = NALL + NRP1 + NRP1 + 1
        IC4 = NALL + NRP1 + NRP1 + NRP1 + 1

!       Changed the next IF since it seems like if nomu=1 then p1 = p-1 and not p 
!       IF (P1 .EQ. P) THEN
        IF (P1 .EQ. (P-1)) THEN
           DO L = 1,P
              IALPHA(L) = (0-IWORKP1(L)) / ALLDAT(IC4)
              WA = WA + (IALPHA(L) * ALLDAT(IC3+R+L))
           END DO
        ELSEIF (P1 .GT. P) THEN
           IF (NOMU .EQ. 0) THEN
              DO H = 2,R
                 IMU(H) = (0-IWORKP1(H-1)) / ALLDAT(IC4)
                 WA = WA + (IMU(H) * ALLDAT(IC3+H))
              END DO
           ENDIF
         
           DO L = 1,P
              L2 = L+1
              IALPHA(L) = (0-IWORKP1(L2)) / ALLDAT(IC4)
              WA = WA + (IALPHA(L) * ALLDAT(IC3+R+L))
           END DO
        ENDIF
      ELSE
        WA = 0.0D0
      ENDIF
      imu(1) = 0
      ! call starting value routine for intercept, 
      ! thresholds, and random-effects variance-covariance matrix
      IF (R .GT. 0 .OR. NGAM .GT. 0) THEN
        CALL STARTV2(IMU(1),ISIGMA,IGAM,R,MAXJ, &
             NGAM,ICATFQ,DENUM,WA,NFN,ICEN,KG,RADD,IDIAG,IVSEP)
      ENDIF

      ! starting values for the TAU scaling terms
      IF (KS .GT. 0) THEN
         ITAU(:) = 0.1D0
      ENDIF

   ENDIF     ! ISTART .NE. 1
   ! write out AND save the starting values

   WRITE(2,"(//,1x,'Starting values',/,1x,'---------------',/)")
   IF (r .ge. 1 .AND. NOMU .EQ. 0) THEN
      WRITE(2,"(1x,'mean       ',10F7.3)") (IMU(H),H=1,R)
      ALLOCATE(IMU0(R))
      CALL RELOC(IMU,IMU0,R,1,0)
   ENDIF
   
   IF (P .GT. 0) THEN
      IF(R.EQ.0) IALPHA(1) = IMU(1)
      WRITE(2,"(1x,'covariates ',10F7.3)") (IALPHA(L),L=1,P)
      ALLOCATE(IALPHA0(P))
      CALL RELOC(IALPHA,IALPHA0,P,1,0)
   ENDIF
   
   IF (r .ge. 1) THEN
      WRITE(2,"(1x,'var. terms ',10F7.3)") (ISIGMA(HR),HR=1,RR)
      ALLOCATE(ISIGMA0(RR))
      CALL RELOC(ISIGMA,ISIGMA0,RR,1,0)
   ENDIF
   
   IF (NGAM .GT. 0) THEN
      WRITE(2,"(1x,'thresholds ',10F7.3)") (IGAM(L),l=1,ngam)
      ALLOCATE(IGAM0(NGAM))
      CALL RELOC(IGAM,IGAM0,NGAM,1,0)
   ENDIF

   IF (KS .GT. 0) THEN
      WRITE(2,"(1x,'scale terms',10F7.3)") (ITAU(L),l=1,ks)
      ALLOCATE(ITAU0(KS))
      CALL RELOC(ITAU,ITAU0,KS,1,0)
   ENDIF
    write(2,*)
    write(2,*)

   ! *****************************************************
   !  start iterations 
   ! *****************************************************
   fail = 0   
   ITLAST = 0 
   RIDGE  = 0.0D0
   RIDGEMAX = 0.0D0
   IFIN   = 1
   !  NORI   = 0
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
   
   allocate(ithetas(n,irt))
   allocate(ithetavs(n,irtt))
   
   ! free up arrays only used as temporary in the init process
   IF (ALLOCATED(CODEX)) DEALLOCATE(CODEX)
   DEALLOCATE( IXMISS, IWMISS, ICATFQ,ICATFQX,IWORKCAT, &
               IWORKP1,IWRKP1,IWP1P1, IXIND, IWIND)
   
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
   ! The SAVE statement preserves the values in local 
   ! arrays between calls to the MainLoop function
   
   ! **********************************************************
   !  Perform a single iteration of the main loop
   ! **********************************************************
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! LOCAL VARS ONLY - SHARED VARS ARE DECLARED
   ! BELOW THE 'MODULE' STATEMENT AT THE TOP OF THE 
   ! FILE
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   INTEGER IC, IC0, IC2, ICOUNT, IDI, IH, IND, INDC, INDD, INDD2, &
      IW, J, K, KK, KKK, L, L0, L1, L2, LM, NII, NQR, NQRR, NSAME, NSAMES, &
      I, IX, PNII, R2, RNII, TEMP_I, RIDGEIT, Radj, RRadj, place, place0
      
   CHARACTER (LEN=8), ALLOCATABLE :: TEMPLABEL(:)
   REAL (KIND=8), ALLOCATABLE :: ICHWRK(:), ICOREC(:), IDERP(:), &
      IDERQ(:), IEV(:), ILIK(:), IWI(:), IWORKR2(:), IWRKR(:), &
      IXI(:), IYI(:)
   REAL (KIND=8) BIGCOR, CSTAT, DERIV, DERIV0, DERIV1, DET, HPROB, PROB, &
      PROBP0, PROBP1, PSUM, QMULT, RLOGDIFF, RSAMES, SCAL, SIGN, STEP, &
      SUMW, WA, WTSUBI, XMU, XTB, Z, Z0, Z1, WTAU, temp, sdev, phiRatio
   
   LOGICAL CATJ
   
   IF(INITIALIZED .EQ. 0) THEN
      CALL POST_ERROR("Mainloop() routine called before Init() " // PROGNAME)
      RETURN
   ELSEIF(INITIALIZED .EQ. 2) THEN
      CALL POST_ERROR("Mainloop() routine called after completion. " // PROGNAME)
      RETURN
   ENDIF
open(23, file='mixor.its')
   
             Radj = IRT
             RRadj = Radj*(Radj+1)/2
   
   IF (IGO .EQ. 1 .OR. ITLAST .EQ. 1) THEN
      
      IT  = IT + 1

      ! do step-halving for first 10 iterations

      IF (IT .LE. 10) THEN
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
         fail = 1
      
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
         IF (NOMU .EQ. 0) THEN
            npar = p+r+rr+ngam+ks
         ELSE
            npar = p+rr+ngam+ks
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

         CALL QUADP(IBQ0,IBQ1,IAQ0,NQ1,NQ,IRT,IUNIF,IAQ1,IAQ2)

         ! get ORIGINAL starting values

         IF (r .ge. 1) THEN
            IF (NOMU .EQ. 0) CALL RELOC(IMU0,IMU,R,1,0)
            CALL RELOC(ISIGMA0,ISIGMA,RR,1,0)
         ENDIF
         
         IF (P .GT. 0) THEN
         
            ! note that R now equals the old R minus 1
            ! so that imu0(R+1) equals the old imu0(R) 
            ! which should be the last element of mu
            
            IALPHA(1) = IMU0(R+1)
            CALL RELOC(IALPHA0,IALPHA(2),P-1,1,0)
         ENDIF
         
         IF (NGAM .GT. 0) THEN
            CALL RELOC(IGAM0,IGAM,NGAM,1,0)
         ENDIF

         IF (KS .GT. 0) THEN
            CALL RELOC(ITAU0,ITAU,KS,1,0)
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
      !                                              
      !  calculate and save the conditional likelihood - ILIK() 
      !                     the marginal likelihood    - hprob 
      !                     the log-likelihood value   - rlogl

      CALL gen(IDER2,0.0d0,npar,npar,1)
      CALL gen(IDER,0.0d0,npar,1,0)
      rlogl = 0.0d0

      IF (IPRIOR .EQ. 1) THEN
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

      ILOOP:DO I = 1,N
   
         hprob = 0.0d0
         CALL gen(IDERP,0.0d0,npar,1,0)
   
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
                  IC           = IC + 1
                  IC2          = ICOUNT + (NRP1 * (K-1) + H+1) 
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
                  IC           = IC + 1
                  IC2          = ICOUNT + (NRP1 * (K-1) + L+R+1)
                  IWI(IC) = ALLDAT(IC2)
               END DO
            END DO
         ENDIF
   
         ! THE EVENT VARIABLE variable EV(K) vector  K = 1 .. NI(I) 
   
         IF (ICEN .GE. 1) THEN
            
            IF (ALLOCATED(IEV)) THEN
               IF(UBOUND(IEV,1) .NE. NII) THEN
                  DEALLOCATE(IEV)
                  ALLOCATE(IEV(NII))
               ENDIF
            ELSE
               ALLOCATE(IEV(NII))
            ENDIF
            
            IC = 0
            DO K  = 1,NII
               IC           = IC + 1
               IC2          = ICOUNT + (NRP1 * (K-1) + P+R+2)
               IEV(IC) = ALLDAT(IC2)
            END DO
         ELSE
         ENDIF
   
         IF (ICEN .EQ. 0) THEN
            ICOUNT = ICOUNT + NII + RNII + PNII 
         ELSE
            ICOUNT = ICOUNT + NII + RNII + PNII + NII
         ENDIF
   
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
            WRITE(IUN,"(1X,'TOTAL NUMBER OF SUBJECTS = ',I6,/)") N
            WRITE(IUN,12) IDI,NII
            12 FORMAT(1X,'DATA FOR SUBJECT ',I10,' WHO HAS',I6,' OBSERVATIONS')
            CALL PRNT(IUN,IYI,NII,1,0,IXLAB,IXLAB,ND,HEAD,1,80,5,1,1, &
                'Y VECTOR                                                   ')
            IF (R .GT. 0) THEN
               CALL PRNT(IUN,IXI,NII,R,0,IXLAB,IXLAB,ND,HEAD,1,80,5,1,1, &
                'X MATRIX                                                   ')
            ENDIF
            
            IF (P .GT. 0) THEN
               CALL PRNT(IUN,IWI,NII,P,0,IXLAB,IXLAB,ND,HEAD,1,80,5,1,1, &
                'W MATRIX                                                   ')
            ENDIF
            
            IF (ICEN .GT. 0) THEN
               CALL PRNT(IUN,IEV,NII,1,0,IXLAB,IXLAB,ND,HEAD,1,80,5,1,1, &
                'EVENT VARIABLE - 1=EVENT 0=CENSOR                          ')
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
         
!Rough Adapted quadrature
      !NQRR = Radj*NQ
      !ibq is stored with all points for 1st dimension followed by all points for second

        ibq = ibq0
        iaq = iaq0
        if(it >= 5 .and. aquad .ne. 0) then
             place = 1
             do k=1,Radj
                 call getSDev(k,Radj,ithetavs(i,:), sdev)
                 do q=1,nq
                     ibq(place) = ithetas(i,k) + sdev*ibq0(place)
                     call get_phi_ratio(ibq(place),ibq0(place), phiRatio)
                     !write(9,*) ibq(place), ibq0(place), phiRatio
                     place0 = mod(place-1,nq)+1
                     iaq(place0) = sdev*phiRatio*iaq(place0)
                     place = place + 1
                     !if(it == 5) write(9,*) k, place0, sdev, phiRatio, iaq(place0)
                 end do
             end do

        end if
                            
         QLOOP:DO q=1,nq    
      
            psum          = 0.0d0
            ILIK(Q)  = 0.0d0
            CALL gen(IDERQ,0.0d0,npar,1,0)
      
            KLOOP:DO k=1,nii   
         
               IF (NGAM .GE. 1) THEN
                  DO L=1,NGAM1
                     IWG(L) = IGAM(L)
                  END DO
               ENDIF
         
               wa  = 0.0D0
               wtau = 0.0D0
               IF (P .GE. 1) THEN
                  DO l=1,p
                     iw = k  + (l-1)*nii
                     wa = wa + IALPHA(L)*IWI(iw)
            
                     IF (KG .GT. 0 .AND. L .LE. KG) THEN
                        DO L0=1,NGAM1
                           LM = L*NGAM1+L0
                           IWG(L0)=IWG(L0)+IGAM(LM)*IWI(IW)*RADD
                     !       IWG(L0) = IWG(L0) + IGAM(LM)*IWI(IW)
                        END DO
                     ENDIF
                     IF (KS .GT. 0 .AND. L .LE. KS) THEN
                        wtau = wtau + ITAU(L)*IWI(iw)
                     ENDIF
                  END DO
               ENDIF
         
               xmu = 0.0D0
               xtb = 0.0D0
               IF (R .GE. 1) THEN
                  DO H = 1,R
                      !ivsep
                     IF (IVSEP .EQ. 0) THEN
                        H2 = Q + (h-1)*NQ
                     ELSE
                        H2 = Q
                     ENDIF
                     IWORKR(H) = IBQ(H2)
                  END DO
                  IF (IDIAG .EQ. 0 .AND. IVSEP .EQ. 0) THEN
                     CALL MPYM(ISIGMA,IWORKR,IWORKR,R,R,3,0,1)
                  ELSE
                     CALL MPYM(ISIGMA,IWORKR,IWORKR,R,R,2,0,1)
                  ENDIF
                  DO h=1,r
                     ix  = k   + (h-1)    * nii
                     IF (NOMU .EQ. 0) THEN
                        xmu = xmu + IMU(H)    * IXI(ix)
                     ELSE
                        XMU = 0.0D0
                     ENDIF
                     xtb = xtb + IWORKR(H) * IXI(ix)
                  END DO
               ENDIF
         
               !   probp0 is j-1 part of numerator  pdf function in derivatives for all 
               !   probp1 is jth part of numerator  pdf function in derivatives for all
               !   prob  is denominator cdf function in derivatives for all parameters
         
               catj = .TRUE.
               j    =  0
               
               DO while (catj)
                  j    = j + 1 
                  IF (j .gt. MAXJ) THEN
                     write(6,*)'response outside range for subject ',idi
                     catj = .FALSE.
                  ELSEIF (FP_EQUAL(IYI(K) , ICODE(J))) THEN 
                     z    = xmu + wa + xtb
         
                     IF (ICEN .EQ. 0) THEN
                        INDC = 1
                     ELSEIF (ICEN .EQ. 1) THEN
                        IF (IEV(K) .EQ. 0) INDC = 0
                        IF (IEV(K) .NE. 0) INDC = 1
                     ENDIF
         
                     IF (j .eq. 1) THEN 
                        z1     =  RADD*z / DEXP(WTAU)
                        IF (INDC.EQ.1) THEN
                           prob   = phifn(z1,nfn)
                           probp0 = 0.0d0
                           probp1 = phiy(z1,nfn)
                        ELSE
                           prob   = 1.0D0 - phifn(z1,nfn)                 
                           probp0 = 0.0d0
                           probp1 = 0.0d0 - phiy(z1,nfn)
                        ENDIF
                        dzt = (0.0d0 - probp1)/prob * z1
                     ELSEIF (j .gt. 1 .and. j .lt. maxj) THEN
                        IF (j .eq. 2) THEN 
                           Z0   = radd*Z  / DEXP(WTAU)
                        ELSE
                           z0   = (IWG(j-2) + RADD*z) / DEXP(WTAU)
                        ENDIF
                        z1   = (IWG(j-1) + RADD*z)  / DEXP(WTAU)
                        IF (INDC.EQ.1) THEN
                           prob   = phifn(z1,nfn) - phifn(z0,nfn)
                           probp0 = 0.0d0 - phiy(z0,nfn)
                           probp1 = phiy(z1,nfn)
                        ELSE
                           prob = 1.0D0 - phifn(z1,nfn)                 
                           probp0 = 0.0d0
                           probp1 = 0.0d0 - phiy(z1,nfn)
                        ENDIF
                     dzt = ((0.0d0 - probp0)*z0 - probp1*z1)/prob
                     ELSEIF (j .eq. maxj) THEN
                        IF (j .eq. 2) THEN 
                           z0   = RADD*z  / DEXP(WTAU)
                        ELSE
                           z0   = (IWG(j-2) + RADD*z) / DEXP(WTAU)
                        ENDIF
                        IF (ICEN.EQ.0) THEN
                           prob   = 1.0d0 - phifn(z0,nfn) 
                           probp0 = 0.0d0 - phiy(z0,nfn)
                           probp1 = 0.0d0
                           dzt = (0.0d0 - probp0)/prob * z0
                        ELSEIF(ICEN.EQ.1) THEN
                           z1     = (IWG(j-1) + RADD*z) / DEXP(WTAU)
                           IF(FP_EQUAL(IEV(K) , 1.0D0)) THEN
                              prob   = phifn(z1,nfn) - phifn(z0,nfn)
                              probp0 = 0.0d0 - phiy(z0,nfn)
                              probp1 = phiy(z1,nfn)
                              dzt = (0.0d0 - probp0)/prob * z0
                           ELSE
                              prob = 1.0D0 - phifn(z1,nfn)                 
                              probp0 = 0.0d0
                              probp1 = 0.0d0 - phiy(z1,nfn)
                              dzt = (0.0d0 - probp1)/prob * z0
                           ENDIF
                        ENDIF
                     ENDIF
            
                     IF (PROB .LT. .1D-305) PROB = .1D-305
                     psum = psum + DLOG(prob)
            
                     !    write(6,*)' before deriv',i,q,k
                     IF (prob  .le. 0.0D0) THEN
                        deriv = 0.0d0
                     ELSE
                        deriv = RADD*(probp0+probp1)/prob
                     ENDIF
            
                     !    write(6,*)' after deriv'
            
                     IF (R .GE. 1 .AND. NOMU .EQ. 0) THEN
                        DO h=1,r
                           ix  = k + (h-1) * nii
                           IDERQ(H) = IDERQ(H)+deriv*IXI(ix)/DEXP(WTAU)
                        END DO
                     ELSE
                     ENDIF
            
                     IF (P .GE. 1) THEN
                        IF (NOMU .EQ. 0) THEN
                           ind = r
                        ELSE
                           ind = 0
                        ENDIF
                        do l = 1,p
                           ind= ind + 1
                           iw = k + (l-1)*nii
                           IDERQ(IND) = IDERQ(IND) + deriv*IWI(iw)/DEXP(WTAU)
                        END DO
                     ELSE
                     ENDIF

                     IF (KS .GT. 0) THEN
                        IF (NOMU .EQ. 0) THEN
                           ind = r + p
                        ELSE
                           ind = p
                        ENDIF
                        DO L=1,KS
                           ind= ind + 1
                           iw = k + (l-1)*nii
                           IDERQ(IND)   = IDERQ(IND) + DZT*IWI(iw)
                        END DO
                     ENDIF
            
                     IF (R .GE. 1) THEN
                        IF (NOMU .EQ. 0) THEN
                           ind = r + p + ks
                        ELSE
                           ind = p + ks
                        ENDIF
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
            
                        IH = 0
                        DO h  = 1,r
                           IF (IDIAG .EQ. 1 .OR. IVSEP .EQ. 1) THEN 
                              IH  = IH  + 1
                              ind = ind + 1
                              IDERQ(IND) = IDERQ(IND) + deriv*IWORKR2(IH) & 
                                           /DEXP(WTAU)
                           ELSE
                           
                              DO h2 = 1,h
                                 IH  = IH  + 1
                                 ind = ind + 1
                                 ! estimate sigma 
                                 IF (IFIN .NE. 0 .OR. H .NE. H2) THEN
                                    IDERQ(IND) = IDERQ(IND)+deriv*IWORKR2(IH) &
                                                 /DEXP(WTAU)
                                 ELSE
                                    ! estimate tau = ln(diagonal elements of sigma)
                                    IF (DERIV .LT. 0.0D0 .AND. IWORKR2(IH)  &
                                        .LT. 0.0D0) SIGN =  1.0D0
                                    IF (DERIV .LT. 0.0D0 .AND. IWORKR2(IH)  &
                                        .GE. 0.0D0) SIGN = -1.0D0
                                    IF (DERIV .GE. 0.0D0 .AND. IWORKR2(IH)  &
                                        .LT. 0.0D0) SIGN = -1.0D0
                                    IF (DERIV .GE. 0.0D0 .AND. IWORKR2(IH)  &
                                        .GE. 0.0D0) SIGN =  1.0D0
                                    IDERQ(IND) = IDERQ(IND) + SIGN *  &
                                        DEXP(DLOG(DABS(deriv)) +DLOG(ISIGMA(IH)) &
                                        +DLOG(DABS(IWORKR2(IH)))) / DEXP(WTAU)
                                 ENDIF
                              END DO
                           ENDIF
                        END DO ! h = 1,r
                     ELSE
                     ENDIF
         
                     !    WRITE(6,*)' before deriv0 and deriv1'
         
                     IF (prob .le. 0.0D0) THEN
                        deriv0 = 0.0d0
                        deriv1 = 0.0d0
                     ELSE
                        deriv0 =         (probp0/prob)  / DEXP(WTAU)
                        deriv1 =         (probp1/prob)  / DEXP(WTAU)
                     ENDIF
            
                     !    WRITE(6,*)' after deriv0 and deriv1'
            
                     IF (NGAM .gt. 0) THEN
                        IF (NOMU .EQ. 0) THEN
                           ic = p+r+rr+ks
                        ELSE
                           ic = p+rr+ks
                        ENDIF
            
                        CSTAT = 1.0D0
                        IF (ICEN .EQ. 1 .AND. IEV(K) .EQ. 0) CSTAT = 0.0D0
            
                        IF (J.GT.2) IDERQ(IC+J-2)=IDERQ(IC+J-2)+deriv0*CSTAT 
                        
                        IF(J.GE.2.AND.J.LE.NGAM1+1) &
                            IDERQ(IC+J-1)=IDERQ(IC+J-1)+deriv1 
         
                        IF(KG.GT.0) THEN
                           IF (NOMU .EQ. 0) THEN
                              ic0= p+r+rr+ngam1
                           ELSE
                              ic0= p+rr+ngam1
                           ENDIF
                           DO L=1,KG
                              ic = ic0+ (l-1)*ngam1
                              iw = k + (l-1)*nii
                              IF (J.GT.2) IDERQ(IC+J-2) =  &
                                   IDERQ(IC+J-2)+deriv0*CSTAT*IWI(iw)*RADD
                              IF (J.GE.2 .AND. J.LE.NGAM1+1) IDERQ(IC+J-1) =  &
                                   IDERQ(IC+J-1)+deriv1*IWI(iw)*RADD
                           END DO
                        ENDIF
                     ENDIF
         
                     catj = .FALSE.
         
                  ENDIF
               END DO !while
         
            END DO KLOOP   ! DO k=1,nii 
      
            ILIK(Q) = DEXP(psum)
            hprob  = hprob + DEXP ( psum + DLOG(IAQ(Q)))
      
            DO L  = 1,npar
               IDERP(L) = IDERP(L) + IDERQ(L)*ILIK(Q)*IAQ(Q)  
            END DO
      
         END DO QLOOP    ! DO q=1,nq 
   
         ! ** done with quadrature points
   
         IF (HPROB .LT. .1D-305) HPROB = .1D-305
         rlogl = rlogl + WTSUBI*DLOG(hprob)
   
         scal       = DEXP( 0.0d0 - DLOG(hprob))
         
         DO l   = 1,npar
            IDER(L)  = IDER(L) + IDERP(L)* DEXP(DLOG(WTSUBI) -DLOG(hprob))
            IDERP(L) = IDERP(L)*scal
         END DO
         
         IF (IPRIOR .EQ. 1) THEN
            DO q   = 1,nq
               IAQEMP(Q) = IAQEMP(Q) + DEXP( (DLOG(WTSUBI)) + &
                  (DLOG(ILIK(Q))) + (DLOG(IAQ(Q))) - DLOG(hprob))
            END DO
         ENDIF
   
         ! scal2 = (1.0d0 + hprob)*wtsubi
         ! CALL grmcv(IDER2,IDER2,IDERP,scal2,npar)
         CALL grmcv(IDER2,IDER2,IDERP,wtsubi,npar)
   
         ! ***************************************************************************
         ! write out the BAYES estimates at final iteration (IFIN=2or3) ONLY IF IRES=1 
   
         IF ((it >= 4 .and. aquad .ne. 0) .or. (IFIN .GE. 2 .AND. IRES .EQ. 1)) THEN
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
            
            ithetas(i,:) = itheta
            ithetavs(i,:) = ithetav
            if (IFIN .GE. 2 .AND. IRES .EQ. 1) then
                IF (IWT .EQ. 0) THEN
                   WRITE(9,"(2I10)")IDI,NII
                ELSEIF (IWT .EQ. 1) THEN
                   WRITE(9,"(I10,F12.5,I15)")IDI,WTSUBI,NII
                ENDIF
                WRITE(9,"(5F15.6)")(ITHETAs(i,H),H=1,IRT)
                WRITE(9,"(5F15.6)")(ITHETAVs(i,HR),HR=1,IRTT)                
            end if
   
         ENDIF
   
      END DO ILOOP  ! end of subject loop:  DO I = 1,N
   
      ! ************************************************
      !  done with subjects                            
      ! ************************************************

      ! CALL PRNT(IUN,IDER,NPAR,1,0,IXLAB,IXLAB,ND,HEAD,1,80,5,1,1,
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
            IDER2(KK)=IDER2(KK)+RIDGE*IDER2(KK)
         END DO
      ENDIF

      ! check to see if the matrix of second derivatives is positive-definite
      ! (nonpos = 1) or not positive definite (nonpos = 0)

      IF (ALLOCATED(ICHWRK)) THEN
         IF(UBOUND(ICHWRK,1) .NE. NPARR) THEN
            DEALLOCATE(ICHWRK)
            ALLOCATE(ICHWRK(NPARR))
         ENDIF
      ELSE
         ALLOCATE(ICHWRK(NPARR))
      ENDIF
! do the inversion of the information matrix via the Cholesky      
      CALL CHSKY(IDER2,ICHWRK,npar,NONPOS)
      CALL INVT(ICHWRK,NPAR,DET)
      CALL CHAMS(ICHWRK,IDER2S,NPAR,3,0)
      CALL GRAMT(IDER2S,IDER2,NPAR,NPAR)
!     CALL INVS(IDER2,npar,det,ICOREC,IER,.FALSE.)
      CALL MPYM(IDER2,IDER,ICOREC,npar,npar,1,0,1)

      ! DO STEP-HALVING FOR FIRST 10 ITERATIONS

      DO L = 1,NPAR
         ICOREC(L) = STEP*ICOREC(L)
      END DO

      BIGCOR = 0.0D0

      IF (r .ge. 1 .AND. NOMU .EQ. 0) THEN
         DO h = 1,r
            IMU(H)    = IMU(H) + ICOREC(H)
            IF (DABS(ICOREC(H)) .GT. BIGCOR) BIGCOR = DABS(ICOREC(H))
         END DO
      ENDIF   

      IF (p .ge. 1) THEN
         DO l = 1,p
            IF (NOMU .EQ. 0) THEN
               ind = l + r
            ELSE
               ind = l
            ENDIF
            IALPHA(L) = IALPHA(L) + ICOREC(IND)
            IF (DABS(ICOREC(IND)) .GT. BIGCOR) BIGCOR = DABS(ICOREC(IND))
         END DO
      ENDIF

! keep the scaling terms at their starting values for the first 10 iterations
      IF (it .gt. 10 .and. ks .ge. 1) THEN
         DO l = 1,ks
            IF (NOMU .EQ. 0) THEN
               ind = l + r + p
            ELSE
               ind = l + p
            ENDIF
            ITAU(L) = ITAU(L) + ICOREC(IND)
            IF (DABS(ICOREC(IND)) .GT. BIGCOR) BIGCOR = DABS(ICOREC(IND))
         END DO
      ENDIF

! keep the variance-covariance parameters at their starting values for the first 5 iterations
      IF (it .gt. 5 .and. r .ge. 1) THEN
         INDD = 1
         INDD2= 0
         
         DO h = 1,rr
            IF (NOMU .EQ. 0) THEN
               ind = r + p + ks + h 
            ELSE
               ind = p + ks + h
            ENDIF
            
            IF (H .EQ. INDD) THEN
               ! on the diagonal
               IF (IDIAG .EQ. 0 .AND. IVSEP .EQ. 0) THEN
                  INDD2 = INDD2 + 1
                  INDD  = INDD  + (INDD2 + 1)
               ELSE
                  INDD  = INDD  + 1
               ENDIF
               
               ! make sure that the first variance term is positive
               IF (IFIN .NE. 0) THEN
                  IF (H.EQ.1 .AND.  &
                     (0.0D0-ICOREC(IND) .GT. ISIGMA(H)))THEN
                     ISIGMA(H) = 0.10D0* DABS(ISIGMA(H) + ICOREC(IND))
                  ELSE
                     ISIGMA(H) = ISIGMA(H) + ICOREC(IND)
                  ENDIF

                  ! DON'T shift to estimation of tau if sigma gets small 
                  ! INCREASE THE RIDGE ONCE INSTEAD
                  !  IF(ISIGMA(H) .LE. .1D0 .AND. NORI .EQ. 0) THEN
                  !     IFIN=0
                  !     NORI=NORI+1
                  !     RIDGE = RIDGE +.1D0
                  !     WRITE(2,"(///,1X,'==> Reparameterized Estimation Occurred')")
                  !  ENDIF
      
               ELSE
                  ISIGTAU(INDD2) = ISIGTAU(INDD2) + ICOREC(IND)
                  ISIGMA(H)  = DEXP(ISIGTAU(INDD2))

                  ! reduce number of random effects if sigma get too small
                  IF(ISIGMA(H) .LE. .000000001D0) ISIG=1
               ENDIF
            ELSE
               ! off the diagonal
               ISIGMA(H)  = ISIGMA(H) + ICOREC(IND)
            ENDIF 
            
            IF (DABS(ICOREC(IND)) .GT. BIGCOR) BIGCOR = DABS(ICOREC(IND))
         END DO ! h = 1,rr
      ENDIF

      IF (ngam .gt. 0) THEN
         L = 0
         DO l1 = 1,KG+1 
            DO l2 = 1,ngam1 
               L = L + 1
               IF (NOMU .EQ. 0) THEN
                  ind = r + p + rr + ks + l
               ELSE
                  ind = p + rr + ks + l
               ENDIF
               
               ! make sure that gamma is increasing 
               ! (NOT FOR THE INTERACTIONS WITH THRESHOLDS)
               
               IGAM(L) = IGAM(L) + ICOREC(IND)
               IF (L1 .EQ. 1 .AND. L2 .GT. 1 .AND. IGAM(L) .LE.  &
                   IGAM(L-1)) IGAM(L) =  IGAM(L-1) + .1D0
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
         CALL FREE_LOCALS()
         INITIALIZED = 2
         RETURN     ! escape from the loop
      ENDIF

      WRITE(23,65)IT
      65 FORMAT(/,1X,'**************',/,1X,'ITERATION ',I4,/, &
                  1X,'**************',/)

      WRITE(23,"(1X,'Log Likelihood   = ',F12.3,//)") rlogl

      ! save current rlogl to check at next iteration

      RLOGDIFF = RLOGL - RLOGLP
      RLOGLP   = RLOGL

      ! calculate versions of RLOGL 
      IF (IWT .EQ. 0) THEN
          SBCN = DBLE(N)
      ELSE
          SBCN = SUM2
      END IF
      AIC  = RLOGL - RNPAR
      SBC  = RLOGL - 0.5 * RNPAR * DLOG(DBLE(sbcn))
      DEV  = -2*RLOGL
      AICD = -2*(RLOGL-RNPAR)
      SBCD = -2*(RLOGL-0.5*RNPAR*DLOG(DBLE(sbcn)))

      IF (r .ge. 1 .AND. NOMU .EQ. 0) THEN
         WRITE(23,"(1X,'mu         ',7F11.6)") (IMU(H),H=1,R)
      ENDIF
      
      IF (P .GT. 0) THEN
         WRITE(23,"(1X,'alpha      ',7F11.6)") (IALPHA(L),L=1,P)
      ENDIF
      
      IF (KS .GT. 0) THEN
         WRITE(23,"(1X,'tau        ',7F11.6)") (ITAU(L),L=1,KS)
      ENDIF
      
      IF (r .ge. 1) THEN
         WRITE(23,"(1X,'sigma      ',7F11.6)") (ISIGMA(HR),HR=1,RR)
      ENDIF
      
      IF (ngam .gt. 0) THEN
         WRITE(23,"(1X,'gamma      ',7F11.6)") (IGAM(L),l=1,ngam)
      ENDIF
      
      WRITE(23,"(1X,'corrections',20F11.6)") (ICOREC(L),L=1,NPAR)
      WRITE(23,"(1X,'derivatives',20g11.3)") (Ider(L),L=1,NPAR)

      ! check to see if the log-likelihood doesn't increase
      ! and if so THEN increase the ridge by 0.1 each time 
!     IF(RLOGDIFF .LT. 0.0D0) RIDGE = RIDGE + 0.1D0

     ! determine if an iteration is bad and up the ridge
     ! take the ridge off after 20 good iterations
     IF (RLOGDIFF/RLOGLP > .000001) THEN
        RIDGEIT = 0
        RIDGE = RIDGE + .1D0
        IF (RIDGE .GT. RIDGEMAX) RIDGEMAX = RIDGE
     END IF
     IF (RLOGDIFF/RLOGLP <= .000001 .AND. RIDGEIT < 20) THEN
        RIDGEIT = RIDGEIT+1
     ELSE IF (RLOGDIFF/RLOGLP <= .000001 .AND. RIDGEIT >= 20) THEN
        RIDGE = 0.0D0
     END IF

      ! check to see if there are numerical problems and if so
      ! reduce the number of random-effects by 1
      IF(IER .EQ. 1 .OR. NONPOS .EQ. 0 .OR. ISIG .EQ. 1) IRBAD =IRBAD+1
      IF(IRBAD .EQ. 1) THEN
         WRITE(6,878)IT, ier, nonpos, isig
         WRITE(2,879)IT, ier, nonpos, isig
         878 FORMAT(///,1X,'==> WARNING!  Estimation Difficulties at Iteration', &
                   I4,/,1x,'==> will proceed with one less random effect', &
                      /,' ier = ',i4,' nonpos = ',i4,' isig = ',i4)
         879 FORMAT(///,1X,'==> WARNING!  Estimation Difficulties Occurred at Iteration', &
                   I4,/,1x,'==> The model was fit with one less random effect than was requested', &
                      /,' ier = ',i4,' nonpos = ',i4,' isig = ',i4)
         fail = 1
      ENDIF      

      ! change weights and nodes for empirical prior

      IF (IPRIOR .EQ. 1) THEN
         SUMW=0.0D0
         
         DO Q = 1,NQ
            IAQ(Q) = IAQEMP(Q) / DBLE(N)
            SUMW  = SUMW   + IAQ(Q)
         END DO
         DO Q = 1,NQ
            IAQ(Q) = IAQ(Q)/SUMW
         END DO
         WRITE(IUN,"(1X,'quad weight',7F10.6)") (IAQ(L),L=1,NQ)

      ENDIF

      ! CHECK IF CONVERGENCE HAS BEEN REACHED
      ! IF NOT RETURN TO ITERATIONS 
      ! IF IT HAS PROCEED TO STANDARD ERROR CALCULATIONS
      ! AND PRINT OUT FINAL RESULTS

      IGO = 0
      if(bigcor < conv) then
         IFIN = IFIN + 2
         IF (RIDGE .GT. 0.0) ITLAST = 1
         IF (IFIN .EQ. 2)    ITLAST = 1
         IF (IFIN .GE. 3 .AND. IRES .EQ. 1) ITLAST = 1
         IF (IFIN .GE. 3 .AND. IRES .EQ. 0) IGO = 0
    else if(it .ge. maxits) then
            WRITE(2,'("NOTE: CONVERGENCE CRITERION WAS NOT ACHIEVED")')
            WRITE(2,'("FINAL FIRST DERIVATIVE AND CORRECTION VALUES")')
            DO L=1,NPAR
              WRITE(2,*)iDER(L),iCORec(L)
            END DO
    else
         IGO = 1
      ENDIF
      
      ! Signal caller that there's more to do
      MainLoop = 1
   ELSE
      
      ! Signal caller that we're done iterations
      MainLoop = 0
      INITIALIZED = 2
      ! free up arrays only used in the main loop, now that 
      ! we are finished with them.  Note that the SAVE statement
      ! preserves the values in these arrays between calls to
      ! the MainLoop function
      CALL FREE_LOCALS()
   ENDIF
   
CONTAINS

! ************************************************
!               **********************                         
!               *  SUBROUTINE FREE_LOCALS()
!               **********************                         
!                                                                       
! Deallocate locally used arrays prior to exiting the routine
! on the final iteration.  It is in a routine because there
! are two exit points from the routine and this must be called
! from each.  Note that the locals MUST be explicitly deallocated
! because otherwise the next session will find they are already
! allocated and will generate an error when it attempts to 
! allocate them.  Tony Gray 3/7/01
!
! CALL FREE_LOCALS()
!                                                                       
! ************************************************
SUBROUTINE FREE_LOCALS()
   IF(ALLOCATED(TEMPLABEL)) DEALLOCATE(TEMPLABEL)
   IF(ALLOCATED(IWRKR)) THEN
      DEALLOCATE(IWRKR)
      DEALLOCATE(IWORKR2)
   ENDIF
   IF (ALLOCATED(IDERP)) DEALLOCATE(IDERP)
   IF (ALLOCATED(IDERQ)) THEN
      DEALLOCATE(IDERQ)
      DEALLOCATE(ILIK)
   ENDIF
   IF (ALLOCATED(IYI)) DEALLOCATE(IYI)
   IF (ALLOCATED(IXI)) DEALLOCATE(IXI)
   IF (ALLOCATED(IWI)) DEALLOCATE(IWI)
   IF (ALLOCATED(IEV)) DEALLOCATE(IEV)
   IF (ALLOCATED(ICOREC)) DEALLOCATE(ICOREC)
   IF (ALLOCATED(ICHWRK)) DEALLOCATE(ICHWRK)
END SUBROUTINE FREE_LOCALS
   
END FUNCTION MainLoop
   
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
   
   !Cleanup vars
   CHARACTER*4  BLANK(40)
   CHARACTER*10 Blab
   CHARACTER*4  SIDEL(2),SIDET
   
   DATA SIDEL/' (1)',' (2)'/
   CHARACTER*2  CH(99)
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
  
   CHARACTER (LEN=8), ALLOCATABLE :: IPARMLAB(:)
   
   REAL (KIND=8), ALLOCATABLE :: ICDER2(:), ICONT(:), ICPVAL(:), &
      ICSE(:), ICVAL(:), ICZVAL(:), IDER2SQ(:), IPVAL(:), &
      ISE(:), ITHETAV2(:), IVAL(:), IZVAL(:), &
      sstar(:), sstar2(:), work(:), sstar2t(:), sigma(:), &
      asstar2(:), adjVar(:), asstar2t(:), varSigma(:)
   REAL (KIND=8) CORR, COV, QMULT, SIDE, SIGVAR, VAR1, VAR2
   
   INTEGER I, I1, I2, I3, I4, IC, ICT, IH, IHI, ILO, IND, IND2, &
      IR, IS, IS2, J, L, L1, L2, NCONK, NCONN, NCT, NCTPAR, pabove, pbelow, &
      asCount, sCount, parCount
   
   IF(INITIALIZED .NE. 2) THEN
      CALL POST_ERROR("Cleanup() routine called out of sequence in " // PROGNAME)
      RETURN
   ENDIF
   
    !adjustment back to sigma
        allocate(sstar(RR*RR))
        allocate(sstar2(RR*RR))
        allocate(sstar2t(RR*RR))
        allocate(work(npar*npar))
        allocate(adjVar(npar*npar))
        allocate(asstar2(npar*npar))
        allocate(asstar2t(npar*npar))
            allocate(sigma(RR))

    if(CHOL .ne. 1) THEN
            if(IVSEP .ne. 1 .and. IDIAG .ne. 1) then
                call getSStar(isigma,r,RR,sstar)
                call scm(sstar,2.0D0,sstar2,RR,RR,0)
                call trp(sstar2,sstar2t,RR,RR)
                call mpym(sstar,isigma,sigma,RR,RR,0,0,1)
            else
                sstar2 = 0
                do i=1,r
                    sstar2(i+(i-1)*r) = isigma(i)*2
                    sigma(i) = isigma(i)*isigma(i)
                end do
                sstar2t = sstar2
            end if

            pabove = p+ks+(1-nomu)*r
            pbelow = ngam
            asstar2 = 0
            asstar2t = 0
            asCount = 1
            sCount = 1
            do i=1,npar
                do j=1,npar
                    if(i > pabove .and. i <= npar - pbelow) then
                        if(j > pabove .and. j <= npar - pbelow) then
                            asstar2(asCount) = sstar2(sCount)
                            asstar2t(asCount) = sstar2t(sCount)
                            sCount = sCount + 1
                        end if
                    else if(i .eq. j) then
                        asstar2(asCount) = 1
                        asstar2t(asCount) = 1
                    end if
                    asCount = asCount + 1
                end do
            end do
                
            call mpym(asstar2,ider2,work,npar,npar,0,1,npar)
            call mpym(work,asstar2t,adjVar,npar,npar,0,0,npar)

    else
        call chams(ider2, adjvar, npar, 1, 0)
        sigma = isigma
    end if

   ALLOCATE(ISE(NPAR))
   ALLOCATE(IVAL(NPAR))
   ALLOCATE(IZVAL(NPAR))
   ALLOCATE(IPVAL(NPAR))
   ALLOCATE(ITHETAV2(RR))
   NCT    = KG*NGAM1
   
   IF (NCON .GT. 0 .OR. KG .GT. 0) THEN
      IF (NCON .GT. NCT) THEN
         NCONK = NCON
      ELSE
         NCONK = NCT
      ENDIF
      ALLOCATE(ICVAL(NCONK))
      ALLOCATE(ICSE(NCONK))
      ALLOCATE(ICZVAL(NCONK))
      ALLOCATE(ICPVAL(NCONK))
      NCONN = (NCONK * (NCONK+1)) / 2
      ALLOCATE(ICDER2(NCONN))
   ENDIF

   ALLOCATE(IPARMLAB(NPAR))
   if(fail==0) then
      write(5,*) rlogl, it
   else
      write(5,*) rlogl, maxits+10
      write(*,*) "Bad run of mixor!"
   end if
   close(5)
   
   WRITE(2,455)
   455 FORMAT(///,1x,'---------------------------------------------------------', &
                /,1x,'* Final Results - Maximum Marginal Likelihood Estimates *', &
                /,1x,'---------------------------------------------------------',/)
   IF (R .EQ. 0) THEN
      WRITE(2,563)IT,RIDGEMAX,RLOGL,AIC,SBC,DEV,AICD,SBCD
      563 FORMAT(1X,'Total Iterations  =',I4, &
               /,1X,'Maximum Ridge     =',F8.3,//, &
         1X,'Log Likelihood                 = ',F12.3,/, &
         1X,"Akaike's Information Criterion = ",F12.3,/, &
         1X,"Schwarz's Bayesian Criterion   = ",F12.3,//,&
         1X,"==> multiplied by -2             ",      /  &
         1X,'Log Likelihood                 = ',F12.3,/, &
         1X,"Akaike's Information Criterion = ",F12.3,/, &
         1X,"Schwarz's Bayesian Criterion   = ",F12.3,/)
   ELSE if(aquad .eq. 0) then
      WRITE(2,564)IT,NQ1,RIDGEMAX,RLOGL,AIC,SBC,DEV,AICD,SBCD
      564 FORMAT(1X,'Total Iterations  =',I4, &
               /,1X,'Quad Pts per Dim  =',I4, &
               1X," (non-adaptive)", &
               /,1X,'Maximum Ridge     =',F8.3,//, &
         1X,'Log Likelihood                 = ',F12.3,/, &
         1X,"Akaike's Information Criterion = ",F12.3,/, &
         1X,"Schwarz's Bayesian Criterion   = ",F12.3,//,&
         1X,"==> multiplied by -2             ",      /  &
         1X,'Log Likelihood                 = ',F12.3,/, &
         1X,"Akaike's Information Criterion = ",F12.3,/, &
         1X,"Schwarz's Bayesian Criterion   = ",F12.3,/)
    else
      WRITE(2,565)IT,NQ1,RIDGEMAX,RLOGL,AIC,SBC,DEV,AICD,SBCD
      565 FORMAT(1X,'Total Iterations  =',I4, &
               /,1X,'Quad Pts per Dim  =',I4, &
               1X," (adaptive)", &
               /,1X,'Maximum Ridge     =',F8.3,//, &
         1X,'Log Likelihood                 = ',F12.3,/, &
         1X,"Akaike's Information Criterion = ",F12.3,/, &
         1X,"Schwarz's Bayesian Criterion   = ",F12.3,//,&
         1X,"==> multiplied by -2             ",      /  &
         1X,'Log Likelihood                 = ',F12.3,/, &
         1X,"Akaike's Information Criterion = ",F12.3,/, &
         1X,"Schwarz's Bayesian Criterion   = ",F12.3,/)
   ENDIF
   
   WRITE(2,57)
   57 FORMAT(/,1X,'Variable',5x,'            Estimate',5X,'Stand. Error', &
        5x,'           Z',5X,'     p-value',/,1X,'----------------',5x, &
        '------------',5X,'------------',5X,'------------',5X,'------------')
   ! *
   ! * Find the z-statistics for all terms  ZFIX and ZVAR
   ! * and their associated p-values        PFIX and PVAR
   ! *
    DO I = 1,Npar
        ISE(I)   = DSQRT(adjvar((i-1)*npar+i))
    end do

 parCount = 0
   IF (R .GT. 0 .AND. NOMU .EQ. 0) THEN
      DO H = 1,R
         IVAL(H)  = IMU(H)
      END DO
      parCount = parCount + R
   ENDIF
   
   IF (P .GT. 0) THEN
      DO L = 1,P
         IVAL(L+parCount) = IALPHA(L)
      END DO
      parCount = parCount + P
   ENDIF
   
   
   IF (KS .GT. 0) THEN
      DO L = 1,KS
         IVAL(L+parCount) = ITAU(L)
      END DO
      parCount = parCount + ks
   ENDIF
   
   IF (R .GT. 0) THEN
       hr = 1
       if(idiag .ne. 1 .and. ivsep .ne. 1) then
             allocate(varSigma(RR))
             varSigma = ise(parCount+1:parCount+RR)
            DO i=1,r
                do j=1,i
                    iVAL(parCount+HR) = sIGMA(r*(j-1)-j*(j-1)/2+i)
                    !Have to shuffle order of standard errors and parameters
                    !Incorrect ordering of packed matrix originally
                    ise(parCount+HR) = varSigma(r*(j-1)-j*(j-1)/2+i)
                    hr = hr + 1
                end do
            END DO
        else
            do i=1,RR
                ival(parcount+HR) = sigma(hr)
                hr = hr + 1
            end do
        end if
      parCount = ParCount+RR
   ENDIF
   
   IF (ngam .gt. 0) THEN
      DO L = 1,ngam
         IVAL(parCount+L) = IGAM(L)
      END DO
   ENDIF


   IC = 0
   DO I = 1,Npar
        IZVAL(I) = IVAL(I) / ISE(I)
        IPVAL(I) = 2.0D0 * (1.0D0 - PHIFN(DABS(IZVAL(I)),0))
   END DO

    parCount = 0

   IF (R .GT. 0 .AND. NOMU .EQ. 0) THEN
      DO H = 1,R
         WRITE(2,"(1X,A16,3(5x,F12.5),F12.5,1X,A4)") &
            IBLABEL(H),IMU(H),ISE(H),IZVAL(H),IPVAL(H),SIDEL(2)
      END DO
      parCount = parCount + R
   ENDIF
   

   IF (P .GT. 0) THEN
      DO L = 1,P
         L1 = L + parCount
         WRITE(2,"(1X,a16,3(5x,F12.5),F12.5,1X,A4)") &
            IALABEL(L),IALPHA(L),ISE(L1),IZVAL(L1),IPVAL(L1),SIDEL(2)
      END DO
      parCount = ParCount + P
   ENDIF

   IF (KS .GT. 0) THEN
      WRITE(2,681)
      681 FORMAT(/,1X,'Scaling terms')
      DO L = 1,KS
         L1 = L + parCount
         WRITE(2,"(1X,a16,3(5x,F12.5),F12.5,1X,A4)") &
            IALABEL(L),ITAU(L),ISE(L1),IZVAL(L1),IPVAL(L1),SIDEL(2)
      END DO
      parCount = parCount + KS
   ENDIF

   IF (r .gt. 0) THEN
      IF (r .eq. 1) THEN
         WRITE(2,581)
         581 FORMAT(/,1X,'Random effect variance term (standard deviation)')
      ELSE if (CHOL .ne. 1) then
         WRITE(2,481)
         481 FORMAT(/,1X,'Random effect variance & covariance terms')
      else
         WRITE(2,4810)
         4810 FORMAT(/,1X,'Random effect variance & covariance terms (Cholesky of var-covariance matrix)')
      ENDIF
      
      IND = 1
      IND2= 0
      DO HR= 1,RR
          ih = parCount + hr
         IF (HR .EQ. IND) THEN
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
         IPVAL(IH) = IPVAL(IH)*side
         WRITE(2,"(1X,a10,3x,f12.5,2(5x,F12.5),F12.5,1X,A4)") &
            Blab,ival(iH),ISE(IH),IZVAL(IH),IPVAL(IH),SIDET
      END DO
      parCount = parCount + RR
   ENDIF
   
   IF (ngam .GT. 0) THEN
      IS = 0
      DO I1 = 1,KG+1
         IF (I1 .EQ. 1) THEN 
            WRITE(2,584)
            584 FORMAT(/,1X,'Thresholds (for identification: threshold 1 = 0)')
         ELSE
            ! remember interaction terms start with i1=2
            WRITE(2,583)IALABEL(I1-1)
            583 FORMAT(/,1X,'Interaction of Thresholds by',A10)
         ENDIF
         
         DO I2 = 1,ngam1
            IS2   = parCount + I2 
            IF (I1 .EQ. 1) THEN
               IPVAL(IS2) = IPVAL(IS2)*0.5D0
               SIDET = SIDEL(1)
            ELSE
               SIDET = SIDEL(2)
            ENDIF
            WRITE(2,"(1X,4X,I2,2x,3(5x,F12.5),F12.5,1X,A4)") &
               I2+1,IGAM(I2),ISE(IS2),IZVAL(IS2),IPVAL(IS2),SIDET
         END DO
      END DO
   ENDIF
   
   write(2,587)
   587 format(//,1x,'note: (1) = 1-tailed p-value',/, &
                   '       (2) = 2-tailed p-value')
   IF (R .EQ. 1 .AND. ICCY .EQ. 2 .AND. IUNIF .EQ. 0  &
      .AND. IPRIOR .EQ. 0) THEN
      sigvar = SIGMA(1)*SIGMA(1)

      IF (nfn .eq. 0) THEN
         CORR   = sigvar / (sigvar + 1.0d0)
         WRITE(2,811)SIGMA(1),SIGMA(1),sigvar,sigvar,sigvar,CORR
         811 FORMAT(///,1X,'Calculation of the intracluster correlation',/, &
               1x,'-------------------------------------------',/, &
               1x,'residual variance = 1 (assumed)',/, &
               1x,'cluster  variance = (',f5.3,' * ',f5.3,') = ',f6.3,//, &
               1x,'intracluster correlation = ',f6.3,' / (',f6.3,' + 1.000)',' = ',F5.3)
  
      ELSEIF (nfn .eq. 1) THEN
         CORR   = sigvar/(sigvar + (pi*pi/3.0D0))
         WRITE(2,8111)ISIGMA(1),ISIGMA(1),sigvar,sigvar,sigvar,CORR
         8111 FORMAT(///,1X,'Calculation of the intracluster correlation',/, &
               1x,'-------------------------------------------',/,  &
               1x,'residual variance = pi*pi / 3 (assumed)',/, &
               1x,'cluster  variance = (',f5.3,' * ', f5.3,') = ',f6.3,//, &
               1x,'intracluster correlation = ',f6.3,' / (',f6.3,' + (pi*pi/3))',' = ',F5.3)
      ELSEIF (nfn .ge. 2) THEN
         CORR   = sigvar/(sigvar + (pi*pi/6.0D0))
         WRITE(2,8110)ISIGMA(1),ISIGMA(1),sigvar,sigvar,sigvar,CORR
         8110 FORMAT(///,1X,'Calculation of the intracluster correlation',/, &
               1x,'-------------------------------------------',/, &
               1x,'residual variance = pi*pi / 6 (assumed)',/, &
               1x,'cluster  variance = (',f5.3,' * ',f5.3,') = ',f6.3,//, &
               1x,'intracluster correlation = ',f6.3,' / (',f6.3,' + (pi*pi/6))',' = ',F5.3)
      ENDIF

   ELSEIF (R .EQ. 2 .AND. IDIAG .EQ. 0 .AND. IVSEP .EQ. 0) THEN
    if(CHOL .eq. 1) then
      VAR1 = ISIGMA(1)*ISIGMA(1)
      COV  = ISIGMA(1)*ISIGMA(2)
      VAR2 = ISIGMA(2)*ISIGMA(2) + ISIGMA(3)*ISIGMA(3)
      CORR = COV / (DSQRT(VAR1) * DSQRT(VAR2))
    else
      COV  = sigma(2)
      CORR = COV / (DSQRT(sigma(1)) * DSQRT(sigma(3)))
    end if
      WRITE(2,812) corr
      812 FORMAT(///,1x,'Covariance expressed as a correlation = ',f6.3)
   ENDIF

   ! GIVE LABELS TO THE PARAMETERS and create the contrast matrix for
   ! the threshold interaction terms
   !
   NCTPAR = NPAR*NCT
   ALLOCATE(ICONT(NCTPAR))
   ICT = 0

   IF (NOMU .EQ. 0)THEN
      DO I = 1,Npar
         IF (I .LE. R) THEN
            IPARMLAB(I) = IBLABEL(I)
            DO I3 = 1,NCT
               ICT = ICT + 1
               ICONT(ICT) = 0.0D0
            END DO

         ELSEIF (I .GT. R .AND. I .LE. R+P) THEN
            I2 = I - R
            IPARMLAB(I) = IALABEL(I2)
            DO I3 = 1,KG
               DO I4 = 1,NGAM1
                  IF (I2 .EQ. I3) THEN 
                     ICT = ICT + 1
                     ICONT(ICT) = 1.0D0
                  ELSE
                     ICT = ICT + 1
                     ICONT(ICT) = 0.0D0
                  ENDIF
               END DO
            END DO

         ELSEIF (I .GT. R+P .AND. I .LE. R+P+KS) THEN
            IPARMLAB(I) = IALABEL(I-R-P)
            DO I3 = 1,NCT
               ICT = ICT + 1
               ICONT(ICT) = 0.0D0
            END DO

         ELSEIF (I .GT. R+P+KS .AND. I .LE. R+P+KS+RR) THEN
            IPARMLAB(I) = 'VarCov' // CH(I-R-P-KS)
            DO I3 = 1,NCT
               ICT = ICT + 1
               ICONT(ICT) = 0.0D0
            END DO

         ELSEIF (I .GT. R+P+KS+RR .AND. I .LE. R+P+KS+RR+NGAM1) THEN
            IPARMLAB(I) = 'Thresh' // CH(I-R-P-KS-RR+1)
            DO I3 = 1,NCT
               ICT = ICT + 1
               ICONT(ICT) = 0.0D0
            END DO

         ELSEIF (I .GT. R+P+KS+RR+NGAM1) THEN
            I2 = I-R-P-KS-RR-NGAM1
            IPARMLAB(I) = 'ThrInt' // CH(I2)
            DO I3 = 1,NCT  
               IF (I2 .EQ. I3) THEN 
                  ICT = ICT + 1
                  ICONT(ICT) = 1.0D0
               ELSE
                  ICT = ICT + 1
                  ICONT(ICT) = 0.0D0
               ENDIF
            END DO
         ENDIF
         ISE(I) = 1.0D0 / ISE(I)
      END DO

   ELSE
      DO I = 1,Npar
         IF (I .LE. P) THEN
            IPARMLAB(I) = IALABEL(I)
            DO I3 = 1,KG
               DO I4 = 1,NGAM1
                  IF (I .EQ. I3) THEN 
                     ICT = ICT + 1
                     ICONT(ICT) = 1.0D0
                  ELSE
                     ICT = ICT + 1
                     ICONT(ICT) = 0.0D0
                  ENDIF
               END DO
            END DO

         ELSEIF (I .GT. P .AND. I .LE. P+KS) THEN
            IPARMLAB(I) = IALABEL(I-P)
            DO I3 = 1,NCT
               ICT = ICT + 1
               ICONT(ICT) = 0.0D0
            END DO

         ELSEIF (I .GT. P+KS .AND. I .LE. P+KS+RR) THEN
            IPARMLAB(I) = 'VarCov' // CH(I-P-KS)
            DO I3 = 1,NCT
               ICT = ICT + 1
               ICONT(ICT) = 0.0D0
            END DO

         ELSEIF (I .GT. P+KS+RR .AND. I .LE. P+KS+RR+NGAM1) THEN
            IPARMLAB(I) = 'Thresh' // CH(I-P-KS-RR+1)
            DO I3 = 1,NCT
               ICT = ICT + 1
               ICONT(ICT) = 0.0D0
            END DO

         ELSEIF (I .GT. P+KS+RR+NGAM1) THEN
            I2=I-P-KS-RR-NGAM1
            IPARMLAB(I) = 'ThrInt' // CH(I2)
            DO I3 = 1,NCT  
               IF (I2 .EQ. I3) THEN 
                  ICT = ICT + 1
                  ICONT(ICT) = 1.0D0
               ELSE
                  ICT = ICT + 1
                  ICONT(ICT) = 0.0D0
               ENDIF
            END DO
         ENDIF
         ISE(I) = 1.0D0 / ISE(I)
      END DO
   ENDIF

   DO I = 1,30
      BLANK(I) = '    '
   END DO
   ND = 0


   ! write out the parameter estimates  
   ! the estimated variance-covariance estimates of these estimates

   DO I = 1,Npar
      !IZVAL(I) = IZVAL(I) / ISE(I)
      WRITE(3,"(1X,a8,1x,F13.7)")IPARMLAB(I),IVAL(I)
   END DO
   
   DO I = 1,Npar
      WRITE(4,"(1X,16F13.7)")(adjvar(I2),I2=(i-1)*npar+1,i*npar)
   END DO
   CLOSE(3)
   CLOSE(4)

   ! CALL MPDSD(ISE,IDER2,IDER2,NPAR)
   ! CALL PRNT(2,IDER2,NPAR,NPAR,1,IPARMLAB,IPARMLAB,ND,
   !    +BLANK,1,78,5,2,2,
   !    +60HCorrelation of the Maximum Marginal Likelihood Estimates    )


   ! compute the EFFECTS OF THE COVARIATES ON THE THRESHOLDS 
   ! (besides the first threshold)               
   ! uses IZVAL and IPVAL as a work vector

   IF (KG .GE. 1) THEN
      CALL MPYM(ICONT,IVAL,ICVAL,NCT,NPAR,0,0,1)
      CALL GRAMM(ICONT,IDER2,ICDER2,NCT,NPAR,1,IZVAL)
      IC = 0
      DO I = 1,NCT 
         DO J = 1,I
            IC = IC + 1
            IF (I .EQ. J) THEN
               ICSE(I)   = DSQRT(ICDER2(IC))
               ICZVAL(I) = ICVAL(I) / ICSE(I)
               ICPVAL(I) = 2.0D0*(1.0D0 - PHIFN(DABS(ICZVAL(I)),0))
            ENDIF
         END DO
      END DO
   ENDIF

   ! print out the requested contrasts of the parameter estimates

   IF (KG .GE. 1) THEN

      WRITE(2,1655)
      1655 FORMAT(///,1x,'------------------------------------',/, &
                      1x,'*  Variable effects on Thresholds  *',/, &
                      1x,'------------------------------------')

      WRITE(2,1757)
      1757 FORMAT(/,1X,'Thresh',2x, 'Variable ', &
         4x,'    Estimate',4X,'Stand. Error',4x,'         Z',    &
         4X,'   p-value',/,1X,'------',2x,'---------',4x,'------------', &
         4X,'------------',4X,'----------',4X,'----------')

      L = 0
      DO L1 = 1,KG
         DO L2 = 1,NGAM1
         L = L+1
            WRITE(2,"(1X,I4,4X,A8,1X,2(4x,F12.5),2(4x,F10.5))") &
               L2+1,IALABEL(L1),ICVAL(L),ICSE(L),ICZVAL(L),ICPVAL(L)
            ICSE(L) = 1.0D0 / ICSE(L)
         END DO
      END DO

      write(2,"(/,1x,'note: p-values are 2-tailed')")

      ! CALL MPDSD(ICSE,ICDER2,ICDER2,NCT)
      ! CALL PRNT(2,ICDER2,NCT,NCT,1,IXLAB,IXLAB,ND,BLANK,1,78,5,1,1,
      !    +60HCorrelation of the MML Transformed Estimates                )

   ENDIF


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
               ICPVAL(I) = 2.0D0 *(1.0D0 - PHIFN(DABS(ICZVAL(I)),0))
            ELSE
            ENDIF
         END DO
      END DO
   ENDIF

   ! print out the requested contrasts of the parameter estimates

   IF (NCON .GE. 1) THEN

      WRITE(2,655)
      655 FORMAT(///,1x,'-------------------------------------',/, &
                     1x,'* Transforms of parameter estimates *',/, &
                     1x,'-------------------------------------')

      CALL PRNT(2,ICON,NPAR,NCON,0,IPARMLAB,IXLAB,ND,BLANK,1,78,5,2,1, &
              'Transpose of the Transform Matrix (parameters by transforms)')
      WRITE(2,757)
      757 FORMAT(/,1X,'Transform',5x,'    Estimate',5X,'Stand. Error',5x, &
         '          Z',5X,'     p-value',/,1X,'---------',5x,'------------',5X, &
         '------------',5X,'------------',5X,'------------')

      DO L = 1,NCON
         WRITE(2,840)L,ICVAL(L),ICSE(L),ICZVAL(L),ICPVAL(L)
         840 FORMAT(1X,I8,4(5x,F12.5))
         ICSE(L) = 1.0D0 / ICSE(L)
      END DO

      write(2,"(/,1x,'note: p-values are 2-tailed')")

      ! CALL MPDSD(ICSE,ICDER2,ICDER2,NCON)
      ! CALL PRNT(2,ICDER2,NCON,NCON,1,IXLAB,IXLAB,ND,
      !    +BLANK,1,78,5,1,1,
      !    +60HCorrelation of the MML Transformed Estimates                )

   ENDIF

   ! print out the empirical prior

   IF (IPRIOR .EQ. 1) THEN

      WRITE(2,656)
      656 FORMAT(///,1x,'--------------------------------', &
                   /,1X,'* Empirical Prior Distribution *', &
                   /,1x,'--------------------------------')

      IF (r .eq. 1) THEN
         WRITE(2,657)
         657 FORMAT(/,1X,' #',5X,'  Weight',4x,'     Node', &
                    /,1x,'--',5x,'--------',4x,'---------')
      ELSEIF (r .gt. 1) THEN
         WRITE(2,6657)
         6657 FORMAT(/,1X,' #',5X,'  Weight',4x,'     Node (r = 1, .. , R)', &
                     /,1x,'--',5x,'--------',4x,'-------------------------')
      ENDIF

      ITHETAV  = 0.0D0    ! Array assignment gets every element
      ITHETAV2 = 0.0D0    ! Array assignment gets every element
      
      DO Q = 1,NQ
         DO H = 1,IRT
            H2 = Q + (h-1)*NQ
            ITHETA(H) = ITHETA(H) + IBQ(H2)*IAQ(Q)
         END DO
      END DO
      
      CALL GRMCV(ITHETAV2,ITHETAV2,ITHETA,DBLE(N),IRT)

      DO Q=1,NQ
         QMULT = DBLE(N) * IAQ(Q)
         DO H = 1,IRT
            H2 = Q + (h-1)*NQ
            IWORKR(H) = IBQ(H2) 
         END DO
         WRITE(2,6659) Q,IAQ(Q),(IWORKR(H2),H2=1,IRT)
         6659 FORMAT(1x,I2,5X,F8.5,5X,8F9.5)
         CALL GRMCV(ITHETAV,ITHETAV,IWORKR,QMULT,IRT)
      END DO

      IH = 0
      DO H = 1,IRT
         DO H2 = 1,H
            IH = IH+1
            ITHETAV(IH)=(ITHETAV(IH)-ITHETAV2(IH)) /DBLE(N-1)
            IF (H.EQ.H2) THEN
               IWORKR(H) = 1.0D0 / DSQRT(ITHETAV(IH))
            ENDIF
         END DO
      END DO

      WRITE(2,"(/,1X,'Empirical Prior Mean Vector')")
      WRITE(2,3388)(ITHETA(H),H=1,IRT)
      WRITE(2,667)
      667 FORMAT(/,1X,'Empirical Prior Variance-Covariance Matrix (in packed form)')
      WRITE(2,3388)(ITHETAV(HR),HR=1,IRTT)
      IF (IRT .GT.1) THEN
         CALL MPDSD(IWORKR,ITHETAV,ITHETAV,IRT)
         WRITE(2,668)
         668 FORMAT(/,1X,'Empirical Prior Correlation matrix (in packed form)')
         WRITE(2,3388)(ITHETAV(HR),HR=1,IRTT)
      ENDIF
       
      3388 FORMAT(5F15.6)

   ENDIF

   !Deallocate local arrays
   DEALLOCATE(ISE, IVAL, IZVAL, IPVAL, ITHETAV2, ICONT, IPARMLAB)
   IF(ALLOCATED(ICVAL)) THEN
      DEALLOCATE( ICVAL, ICSE, ICZVAL, ICPVAL, ICDER2 )
   ENDIF
   
   !Deallocate global arrays
   DEALLOCATE( ALLDAT, IDNI )
   IF (ALLOCATED(IGAM0)) DEALLOCATE(IGAM0)
   IF (ALLOCATED(ISIGMA0)) DEALLOCATE(ISIGMA0)
   IF (ALLOCATED(IALPHA0)) DEALLOCATE(IALPHA0)
   IF (ALLOCATED(ITAU0)) DEALLOCATE(ITAU0)
   IF (ALLOCATED(IMU0)) DEALLOCATE(IMU0)
      
   DEALLOCATE ( IALPHA, IAQ, IAQ1, IAQ2, IAQEMP, IBQ, IBQ1, ICODE,  &
      ICON, IDER, IDER2, IGAM, IMU, ISIGMA, ISIGTAU, ITHETA, ITHETAV, &
      IWG, IWORKR, ITAU, IDER2S )
      deallocate(iaq0, ibq0, ithetas, ithetavs)
   DEALLOCATE( IBLABEL, IALABEL, IXLAB)
    deallocate(sstar)
    deallocate(sstar2)
    deallocate(sstar2t)
    deallocate(sigma)
    deallocate(work)
    deallocate(adjVar)
    deallocate(asstar2)
    deallocate(asstar2t)
    
   CLOSE(2)
   CLOSE(9)
   INITIALIZED = 0
   
END SUBROUTINE CLEANUP

! *************************************************************
!  SUBROUTINTE STARTV2 (MU1,SIGMA,GAMMA,R,MAXJ,NGAM,CATFQ, 
!              DENUM,WA,NFN,ICEN,KG,RADD,IDIAG,IVSEP)     
!                                                         
!  GET STARTING VALUES FOR INTERCEPT, THRESHOLDS, AND     
!  VARIANCE-COVARIANCE MATRIX OF RANDOM EFFECTS           
!                                                         
!  Parameters Sent                                        
!  R      = Number of Random Effects                      
!  MAXJ   = Number of Ordered Categories                  
!  NGAM   = Number of thresholds to estimate              
!  CATFQ  = MAXJ Category Frequencies                     
!  DENUM  = Total Number of Level-1 observations          
!  WA     = Predicted value of Y minus the intercept      
!  NFN    = 0 for Probit                                  
!         = 1 for Logistic                                
!         = 2 for Complementary Log-Log                   
!         = 3 for Log-Log                                 
!  RADD   = 1 add alpha and mu to thresholds              
!           (covariates and mean of random effects)       
!         =-1 subtract beta and mu from thresholds        
!  IDIAG  = 0 for correlated random effects               
!         = 1 for independent random effects              
!  IVSEP  = 0 R random effects aren't indicator variables 
!         = 1 R random effects are indicator variables    
!           (between-subjects indicator variables)        
!                                                         
!  Parameters Returned                                    
!  MU1    = Starting value for intercept                  
!  GAMMA  = NGAM Starting values for Thresholds           
!  SIGMA  = (R*(R+1))/2 Starting values for Variance      
!           Covariance matrix of random effects           
!                                                         
! *************************************************************

SUBROUTINE STARTV2(MU1,SIGMA,GAMMA,R,MAXJ,NGAM,CATFQ,DENUM,WA, &
                      NFN,ICEN,KG,RADD,IDIAG,IVSEP)
   IMPLICIT REAL*8(A-H,O-Z),INTEGER*4(I-N)
   DOUBLE PRECISION MU1,SIGMA,GAMMA,DENUM,WA,CATFQ,pi
   INTEGER R
   DIMENSION SIGMA(1),GAMMA(1),CATFQ(1)

   PI=3.141592654D0

   cumsum = 0.0D0
   DO j=1,maxj-1
      cumsum   = cumsum + CATFQ(J)
      cumprob  = cumsum / DENUM
      cumodds  = cumsum / (DENUM - cumsum)
      Rlncumod = DLOG(cumodds)
      IF (nfn .eq. 0) THEN
         tempparm = 0.625 * (Rlncumod + wa)
      ELSEIF (nfn .eq. 1) THEN
         tempparm = Rlncumod + wa
      ELSEIF (nfn .eq. 2) THEN
         tempparm = DLOG(0.0D0 - DLOG(1.0d0 - cumprob)) + wa
      ELSEIF (nfn .eq. 3) THEN
         tempparm = DLOG(0.0D0 - DLOG(cumprob)) + wa
      ELSE
      ENDIF
      IF (j .EQ. 1) THEN
         MU1 = RADD*tempparm
      ELSE
         GAMMA(j-1) = tempparm - RADD*MU1
      ENDIF
   END DO

   ! FOR RIGHT CENSORING AN EXTRA GAMMA NEEDS TO BE ESTIMATED
   !  arbitrarily assume half of the last category is censored 

   NGAM1 = MAXJ-2
   IF (ICEN .EQ. 1) THEN
      NGAM1 = MAXJ-1
      cumprob  = cumprob + (0.5d0 * (1.0d0 - cumprob))
      cumodds  = cumprob / (1.0d0 - cumprob)
      Rlncumod = DLOG(cumodds)
      IF (nfn .eq. 0) THEN
         tempparm = 0.625 * (Rlncumod + wa)
      ELSEIF (nfn .eq. 1) THEN
         tempparm = Rlncumod + wa
      ELSEIF (nfn .eq. 2) THEN
         tempparm = DLOG(0.0D0 - DLOG(1.0d0 - cumprob)) + wa
      ELSEIF (nfn .eq. 3) THEN
         tempparm = DLOG(0.0D0 - DLOG(cumprob)) + wa
      ELSE
      ENDIF
      GAMMA(MAXJ-1) = tempparm - RADD*MU1
   ELSE 
   ENDIF

   ! FOR interactions with GAMMA - assign the same GAMMAs for all

   IF (KG .GT. 0) THEN
      DO KGN = NGAM1+1,NGAM
         KGM = MOD(KGN,NGAM1)  
         IF (KGM .EQ. 0) KGM=NGAM1
         GAMMA(KGN) = GAMMA(KGM)*0
      END DO
   ENDIF


   IF (R .EQ. 1) THEN
      SIGMA(1) = .31623d0
      IF (nfn .EQ. 1) SIGMA(1) = SIGMA(1) * DSQRT(pi*pi/3.0D0)
      IF (nfn .GE. 2) SIGMA(1) = SIGMA(1) * DSQRT(pi*pi/6.0D0)
   ELSEIF (R .GT. 1 .AND. IDIAG .EQ. 0 .AND. IVSEP .EQ. 0) THEN
      IR = 0
      DO I  = 1,R
         DO I2 = 1,I
            IR = IR + 1
            IF (I2 .EQ. I) THEN
               IF (I .EQ. 1) THEN
                  SIGMA(IR) = 1.0d0
               ELSE
                  SIGMA(IR) = 0.5d0
               ENDIF
            ELSE
               SIGMA(IR) = 0.0d0
            ENDIF
            IF (nfn .EQ. 1) SIGMA(IR) = SIGMA(IR)*DSQRT(pi*pi/3.0D0)
            IF (nfn .GE. 2) SIGMA(IR) = SIGMA(IR)*DSQRT(pi*pi/6.0D0)
         END DO
      END DO
   ELSEIF (R .GT. 1 .AND. (IDIAG .EQ. 1 .OR. IVSEP .EQ. 1)) THEN
      IR = 0
      DO I  = 1,R
         IR = IR + 1
         IF (I .EQ. 1) THEN
            SIGMA(IR) = 1.0d0
         ELSE
            SIGMA(IR) = 1.0d0
         ENDIF
         IF (nfn .EQ. 1) SIGMA(IR) = SIGMA(IR)*DSQRT(pi*pi/3.0D0)
         IF (nfn .GE. 2) SIGMA(IR) = SIGMA(IR)*DSQRT(pi*pi/6.0D0)
      END DO
   ENDIF

   RETURN
END SUBROUTINE STARTV2

END MODULE RRM_MODULE

! ************************************************
! The DLL wrapper and EXE stub routines are identical for all 
! four modules, but need to be included in each file rather than
! compiled as as module.  So the INCLUDE directive is used to 
! avoid duplication - a single file is used in all modules.
! ************************************************
INCLUDE "DLLSTUB.F90"
