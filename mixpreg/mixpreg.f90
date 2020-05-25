!  **************************************************************
!  MIXPREG        
!
!  modified on 
!  10/04/02 to compute EB estimates for multiple random effects
!
!  10/08/02 to correct AIC & BIC statistics
!
!  modified on 8/5/04 to correct BIC statistics - use N not NTOT
!
!
! *************************************************************
!  POISSON RRM MAIN PROGRAM                              
!                                                        
!  Model                                                 
!  Y = X BETA + W ALPHA + ERROR                          
!                                                        
!  Indices                                               
!     N  = TOTAL NUMBER OF SUBECTS         I = 1 .. N    
!     R  = DIMENSION OF RANDOM BETA SPACE  H = 1 .. R    
!     P  = DIMENSION OF FIXED ALPHA SPACE  L = 1 .. P    
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
! *****************************************************************
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
! this program was modified on 9/27/96
! a) allow the variance covariance matrix of the random effects
!    to be a diagonal matrix (uncorrelated random effects)
!    (IDIAG=1)
! b) allow the mean of the random effects
!    to be zero 
!    (NOMU=1)
!
! this version was updated 9/27/96 to 
! a) write out estimates to MIXPREG.EST (with labels)
! b) write out asymptotic variance-covariance matrix (in unpacked form)
!    to MIXPREG.VAR (without labels - in the order specified in MIXPREG.EST)
!
! this version was updated 10/9/96 to include
!
! IVSEP  =  0 R random effects don't indicate separate groups of subjects
!           1 R random effects do    indicate separate groups of subjects
!            (or in IRT case, X variables for R random effects change
!             but there is only 1 random effect beta)
!

! Note that the module that INCLUDE's this must define a character
! constant called PROGNAME that contains the program name without
! the 'b' on the end, like MIXOR,MIXNO,MIXREG or MIXPREG.
PROGRAM RRM_EXE
       INTEGER ,POINTER ::IDNI(:)
   INTEGER I2, IC, IC2, IDIAG, IER, IFIN, IGO, IND, IND2,  &
      IOFS, IPRIOR, IR, IRBAD, IRES, ISIG, IT, ITLAST, IUN, IUNIF, IVSEP, &
      IWT, J, J3, K, L, L1, N, NCON, ND, NOMU, NONPOS, NPAR, NPARR, &
      NPR, NQ, NQ1, NQRR, NRP1, IRT, IRTT 
   INTEGER H, HR, P, R, RR, ICCY 
   INTEGER :: INITIALIZED = 0
   
   CHARACTER (LEN=4), ALLOCATABLE :: IXLAB(:)
   CHARACTER*4 HEAD(30)
   CHARACTER (LEN=16), ALLOCATABLE :: IALABEL(:), IBLABEL(:)
   
   REAL (KIND=8) CONV, DET, GMEAN, RIDGE, RLOGL, RLOGLP, AIC, SBC, &
                 DEV, AICD, SBCD, RNTOT, RNPAR
   REAL (KIND = 8),POINTER :: ALLDAT(:)   
   
   REAL (KIND=8), ALLOCATABLE :: IALPHA(:), IALPHA0(:), IAQ(:), IAQ1(:), &
      IAQ2(:), IAQEMP(:), IBQ(:), IBQ1(:), ICON(:), IDER(:), &
      IDER2(:), IMEAN(:), IMU(:), IMU0(:), ISIGTAU(:), IWORKR(:), &
      ITHETA(:), ITHETAV(:)      
   
   REAL (KIND=8), ALLOCATABLE,TARGET :: ISIGMA(:), ISIGMA0(:)
   REAL (KIND=8), ALLOCATABLE :: ICHWRK(:), ICOREC(:), IDERP(:), &
      IDERP2(:), IDERQ(:), IDERQ2(:), IDZ(:), ILIK(:), &
      IOV(:), IWI(:), IWORKR2(:), IWRKR(:), IXI(:), IYI(:)
     
CHARACTER, PARAMETER:: PROGNAME*7 = 'MIXNO'
   CHARACTER (LEN=8), ALLOCATABLE :: TEMPLABEL(:) 
	

   CALL INIT()  ! supply command line params to INIT to override
                       ! filenames for DEF, DAT, etc.
   DO WHILE(MainLoop() .NE. 0) 
   ENDDO
   CALL CLEANUP()

contains   
     
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
   
   INTEGER I,IC3, IC4, ICL, ICM, ICX, ICXW, IDIND, IHI, &
      II, IR2, ISTART, K2, K3, K4, KIN1, KIN2, KIN3, &
      KIND, MAXJXJ, MAXJXJ2, MAXK, MAXXJ, MEANYX, MISS, &
      NALL, NCATYX, NCPAR, NF, NFF, NQR, NTOT, NXLAB 
   INTEGER OFIND, WTIND, XCIND, MAXCOL, YIND
   
   REAL (KIND=8) SUM2, UNO, WSUM, YMISS, YVAL, YVAL2 
   
   CHARACTER*8 TLABEL, XLABEL, YLABEL
   
   CHARACTER*80 FILEDAT, FILEOUT
   
   INTEGER ,ALLOCATABLE::IWIND(:), IXIND(:)
   
   LOGICAL CATJ  
   REAL (KIND=8), ALLOCATABLE ::CODEX(:), ICATFQX(:), &
      IWMISS(:), IWNFNF(:), IWORKNF(:), IWRKNF(:), IXMISS(:) 

   ! Start by checking if we are ready to run yet.  You 
   ! can't call init twice in a row, and must call 
   ! cleanup before you can start another session.
   
   IF(INITIALIZED .NE. 0) THEN
      CALL POST_ERROR("Init() routine called out of sequence in " // PROGNAME)
      RETURN
   ENDIF
   
   !  Start by reading in the DEF file that describes the dataset
   !  and the desired statistical analysis
   !  Since the sizes of the arrays read in later in the DEF file
   !  depend on the numerical values read earlier in the file, some
   !  Allocation must be done during the file reading.  


   ! MAXCOL is the dimension of the data matrix to read in     
   !   
   ! IDIND (integer) indicates which column contains the level 2 ID
   !   
   ! YIND (integer) indicates which column contains the dependent var
   !   
   ! R is the dimension of the X matrix to be used in the analysis
   !   
   ! IXIND (integer) indicates which columns are to be used for X           
   !   
   ! P is the dimension of the W matrix to be used in the analysis
   !   
   ! IWIND (integer) indicates which columns are to be used for W           
   !   
   ! WTIND (integer) indicates which column to be used for the 
   !       weighting of level-2 units           
   !   
   ! OFIND (integer) indicates which column contains the offset variable
   !   
   !   
   ! ISTART = 0 FOR AUTOMATIC STARTING VALUES
   !          1 FOR STARTING VALUES READ IN
   !
   ! IWT    = 0 EACH PERSON (2ND LEVEL UNIT) IS WEIGHTED EQUALLY
   !          1 FOR DIFFERENTIAL WEIGHTING FOR EACH PERSON (2ND LEVEL)
   !
   !
   ! MEANYX = 0 NO MEANS FOR ANY (X or W) BY Y                        
   !          1 GET MEANS FOR A SPECIFIED (X or W) BY Y   
   !          (the x or w variable is specified below 
   !           as being in column XCIND)
   !
   ! IPRIOR = 0 specified PRIOR FOR RANDOM EFFECTS                 
   !          1 EMPIRICAL PRIOR FOR RANDOM EFFECTS                     
   !
   ! IUNIF  = 0 NORMAL PRIOR FOR RANDOM EFFECTS                 
   !          1 rectagular PRIOR FOR RANDOM EFFECTS                     
   !          2 log gamma  PRIOR FOR RANDOM EFFECTS                     
   !
   ! NQ1    = number of quadrature nodes per dimension          
   !
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
   !    IOFS   =  0 no offset variable                                         
   !          1 read in an offset variable                                  
   !
   !
   !  IFIN = 0 ESTIMATE TAU (ln of diagonaol elements of cholesky)
   !       = 1 ESTIMATE SIGMA (cholesky of variance-covariance
   !                           matrix of random effects)
   !

   !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   ! READ THE DEF FILE THAT DESCRIBES THE DATASET AND THE 
   ! ANALYSIS OPTIONS.
   !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


   ! READ IN TWO LINES (60 CHARS EACH) FOR THE TITLE 
   ! AND THEN PARAMETERS AND STARTING VALUES FROM MIXPREG.DEF

!  OPEN(1,ACTION='READ,DENYWRITE', FILE=FILENAMES(1))
   OPEN(1,ACTION='READ', FILE="mixpreg.def")
   READ(1,"(15A4)") HEAD
   READ(1,"(A80)")FILEDAT
   READ(1,"(A80)")FILEOUT
   READ(1,"(80X)")         ! Skip over the DEF file name
   
   OPEN(2, FILE= fileout)
   OPEN(3, FILE= 'mixpreg.EST')
   OPEN(4, FILE= 'mixpreg.VAR')

   READ(1,*) NPR,MAXCOL,R,P,CONV,MISS,ISTART,IWT,MEANYX,IPRIOR,IUNIF, &
             NQ1,NCON,IDIAG,NOMU,IVSEP,IOFS

   ! IRES = 0 DON'T PRINT OUT INDIVIDUAL PARAMETERS          
   !      = 1 PRINT INDIVIDUAL PARAMETERS TO FILERES         

   write(*,*) NPR,MAXCOL,R,P,CONV,MISS,ISTART,IWT,MEANYX,IPRIOR,IUNIF, &
             NQ1,NCON,IDIAG,NOMU,IVSEP,IOFS
   IRES=0
   IF(R.GE.1) IRES=1
   IF (IRES .EQ. 1) OPEN(5, FILE='mixpreg.RES')

   WRITE(6,*) R,' random terms'
   WRITE(6,*) P,' fixed  terms'

   ! ALLOCATE  IWORKR,IXMISS,ISIGTAU,IXIND,   length: R 
   ALLOCATE (IXMISS(R))
   ALLOCATE (IXIND(R))
   ALLOCATE (IWORKR(R))
   ALLOCATE (ISIGTAU(R))
   
   ! IRT = dimension of theta

   IF (IVSEP .EQ. 0) THEN
      IRT = R
   ELSEIF (IVSEP .EQ. 1) THEN
      IRT = 1
   ENDIF          
   ALLOCATE(ITHETA(IRT))

   ! ALLOCATE  ISIGMA, length: rr = (r * (r+1) / 2) 

   IF (IDIAG .EQ. 0 .AND. IVSEP .EQ. 0) THEN
      RR = (R * (R+1)) / 2
   ELSE
      RR = R
   ENDIF
   ALLOCATE(ISIGMA(MAX(RR,1)))
   
   IRTT = (IRT*(IRT+1))/2
   ALLOCATE(ITHETAV(IRTT))

   ! ALLOCATE  IBQ1, length: NQ1 (nquad points per dimension)

   ALLOCATE(IBQ1(NQ1))

   ! ALLOCATE  IAQ, IAQEMP, length: nq1**r or nq1 if IVSEP<> 0

   IF (IVSEP .EQ. 0) THEN
      NQR = NQ1**R
   ELSE
      NQR = NQ1
   ENDIF
   ALLOCATE(IAQ(NQR))
   ALLOCATE(IAQEMP(NQR))

   ! ALLOCATE  IBQ, length: NQRR = R * (nq1**r) or NQ1 if IVSEP <> 0

   IF (IVSEP .EQ. 0) THEN
      NQRR = R*NQR
   ELSE
      NQRR = NQ1
   ENDIF
   ALLOCATE(IBQ(NQRR))

   ! ALLOCATE  IMU, length: R    IALPHA, length: P
   ! ALLOCATE  these vectors consectutively 
   !           so that random-effects can 
   !           be reassigned as fixed-effects
   !           if numerical problems develop 

   ALLOCATE(IMU(R))
   ALLOCATE(IALPHA(P))

   ! ALLOCATE  IBLABEL, length: R for random effect labels
   ! ALLOCATE  IALABEL, length: P for covariate labels

   ALLOCATE(IBLABEL(R))
   ALLOCATE(IALABEL(P))

   ! ALLOCATE  IWMISS,IWIND length: p 
   ALLOCATE(IWMISS(P))
   ALLOCATE(IWIND(P))
   
   ! ALLOCATE    IDER2, length: NPARR = (npar * (npar+1)) / 2
   !   where    npar = r + rr + p 

   IF (IPRIOR .EQ. 0 .AND. NOMU .EQ. 0) THEN
      NF   = P+R
   ELSEIF (IPRIOR .EQ. 1 .OR. NOMU .EQ. 1) THEN
      NF   = P
   ENDIF
   NPAR = NF + RR
   RNPAR = DBLE(NPAR)

   NCPAR = NPAR*NCON
   NPARR = (NPAR * (NPAR+1)) / 2
   ALLOCATE(ICON(NCPAR))
   ALLOCATE(IDER2(NPARR))

   ! ALLOCATE  IMEAN, length: NF  (either r+p  or  p)

   NFF = NF
   IF (IOFS .EQ. 1) NFF = NFF+1
   ALLOCATE(IMEAN(NFF))

   READ(1,*) IDIND,YIND

   IF (R .GE. 1) READ(1,*) (IXIND(H), H=1,R)
   IF (P .GE. 1) READ(1,*) (IWIND(L), L=1,P)
   IF (IWT .EQ. 1) READ(1,*) WTIND
   IF (IOFS .EQ. 1) READ(1,*) OFIND

   IF (MEANYX .EQ. 1) THEN
      
      READ(1,*) XCIND,MAXXJ
      ALLOCATE(CODEX(MAXXJ))
      BACKSPACE(1)
      READ(1,*) XCIND,MAXXJ,(CODEX(J),J=1,MAXXJ)

      CATJ = .TRUE.
      IC = 0
      DO WHILE (CATJ)
         IC = IC + 1
         IF (IC .LE. R) THEN
            IF (XCIND .EQ. IXIND(IC)) THEN
               ICXW = 1
               ICX  = IC
               NCATYX = IC + 1
               CATJ = .FALSE.
            ELSE
            ENDIF
         ELSEIF (IC .GT. R .AND. IC .LE. R+P) THEN
            IC2 = IC - R
            IF (XCIND .EQ. IWIND(IC2)) THEN
               ICXW = 2
               ICX  = IC2
               NCATYX = IC2 + R + 1
               CATJ = .FALSE.
            ELSE
            ENDIF
         ELSE
            WRITE(6,*)'XCIND does not match XIND or WIND'
            CATJ = .FALSE.
            MEANYX = 0
         ENDIF
      END DO
   ENDIF

   IF (MISS .EQ. 1) THEN
      READ(1,*) YMISS
      IF (R .GE. 1) READ(1,*) (IXMISS(H), H=1,R)
      IF (P .GE. 1) READ(1,*) (IWMISS(L), L=1,P)
   ELSE
   ENDIF

   READ(1,*) YLabel

   IF (R .NE. 0) THEN
      READ(1,*) (IBLABEL(H), H=1,R)
      IF (NOMU .EQ. 0 .AND. (ISTART .EQ. 1 .or. IPRIOR .EQ. 1)) &
          READ(1,*) (IMU(H), H=1,R)
   ENDIF

   IF (P .NE. 0) THEN
      READ(1,*) (IALABEL(L), L=1,P)
      
      IF (ISTART.EQ.1)READ(1,*)(IALPHA(L), L=1,P)
      
      IF (MEANYX .EQ. 1) THEN
         IF (ICXW .EQ. 1) XLabel = IBLABEL(ICX)
         IF (ICXW .EQ. 2) XLabel = IALABEL(ICX)
      ENDIF
   ENDIF

   IF (R .NE. 0 .AND. ISTART .EQ. 1) READ(1,*) (ISIGMA(HR), HR=1,RR)
   
   IF (NCON .NE. 0) READ(1,*) (ICON(L), L=1,NCPAR)

   CLOSE(1)
   
   !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   ! END of DEF file read
   !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

   ! get quadrature nodes & weights  

   ALLOCATE(IAQ1(NQ1))
   ALLOCATE(IAQ2(NQRR))
   
   IF (IVSEP .EQ. 0) THEN
      CALL QUADP(IBQ,IBQ1,IAQ,NQ1,NQ,R,IUNIF,IAQ1,IAQ2,ISIGMA(1))
   ELSE
      CALL QUADP(IBQ,IBQ1,IAQ,NQ1,NQ,1,IUNIF,IAQ1,IAQ2,ISIGMA(1))
   ENDIF


   ! parameters         

   NRP1 = 1 + R + P
   IF (IOFS .EQ. 1) NRP1 = NRP1 + 1
   UNO  = 1.0D0

   ! GET THE DATA 
   ! N    (total number of level 2 units)
   ! NTOT (total number of observations)

   CALL READAT(FILEDAT,N,NTOT,MAXK,MAXCOL,R,P,ALLDAT,IDNI,IDIND, &
               YIND,IXIND,IWIND,MISS,YMISS,IXMISS,IWMISS, &
               IWT,WTIND,IOFS,OFIND)
  
   !WRITE(6,"(' Sizeof(Alldat): ',I6,' UBOUND(IDNI): ',I6)")  &
   !      UBOUND(ALLDAT),UBOUND(IDNI)
   !WRITE(6,"(' first 10 of alldat: ',10F15.6)")(ALLDAT(J),J=1,10)
   !WRITE(6,"(' first 10 of idni: ',10I6)")(IDNI(J),J=1,10)

   ! NALL is the number of elements for the data read in from filedat 
   ! that ALLDAT contains - all other data put in BM must follow
   ! this point

   RNTOT = DBLE(NTOT)

   IF (IWT .EQ. 0) THEN
      NALL = NRP1*NTOT  
   ELSEIF (IWT .EQ. 1) THEN
      NALL = NRP1*NTOT + N
   ENDIF
   !WRITE(6,"(' NALL(',I4,') should equal size of alldat')") NALL*2+1

   ! XLAB dimension is largest of maxni or npar (= p + r + rr)
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
   
   MAXJXJ = 3*MAXXJ
   ALLOCATE(ICATFQX(MAXJXJ))
   
   IF (MEANYX .EQ. 1) THEN
      DO K3 = 1,MAXXJ
         DO K  = 1,3
            KIND = ((K-1)*MAXXJ) + K3
            ICATFQX(KIND) = 0.0D0
         END DO
      END DO
   ELSE
   ENDIF

   II = 1
   WSUM = 0.0D0
   SUM2 = 0.0D0
   
   DO K  = 1,NALL
      K2 = K - (II - 1)
      IR     = MOD(K2,NRP1)
      IF (IR .EQ. 0) IR = NRP1
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

      ALLDAT(IC3) = ALLDAT(IC3) + ALLDAT(K)

      ! for one of the Xs or Ws get the sum & n  by Y

      IF (MEANYX .EQ. 1) THEN
         IF (IR .EQ. 1) THEN
            YVAL      = ALLDAT(K)
         ELSEIF (IR .EQ. NCATYX) THEN
            CATJ = .true.
            K4    =  0
            DO WHILE (CATJ)
               K4    = K4 + 1 
               ! Test for equality with FP-safe comparison  TG
               IF (FP_EQUAL(ALLDAT(K) , CODEX(K4))) THEN
                  KIND = K4
                  ICATFQX(KIND) = ICATFQX(KIND) + YVAL
                  KIND = (2*MAXXJ) + K4
                  ICATFQX(KIND) = ICATFQX(KIND) + 1
                  CATJ = .FALSE.
               ELSE
               ENDIF
            END DO
         ENDIF
      ENDIF

   END DO


   ! calculate the means

   ICCY = 0
   ICM  = 0
   
   DO K  = 1,NRP1
      IC3    = NALL + NRP1 + NRP1 + K
      ALLDAT(IC3) = ALLDAT(IC3) / DBLE(NTOT)

      ! see if the mean of the 1 random effect is 1
      ! for intraclass correlation calculation later

      IF (K .EQ. 2 .AND. R .EQ. 1 .AND. ALLDAT(IC3) .GT.  &
          0.999D0 .AND. ALLDAT(IC3) .LT. 1.001D0)  &
          ICCY = ICCY + 1

      ! save the means of the fixed effects variables AND OFFSET 

      IF ((K.GE.2 .AND. NF.GT.P) .OR. (K.GE.R+2 .AND. NF.EQ.P))THEN
         ICM = ICM+1
         IMEAN(ICM) = ALLDAT(IC3)
      ENDIF
   END DO

   IF (MEANYX .EQ. 1) THEN
      DO K4 = 1,MAXXJ
         KIN3 = (2*MAXXJ) + K4
         KIN1 = K4
         ICATFQX(KIN1) = ICATFQX(KIN1) / ICATFQX(KIN3)
      END DO
   ENDIF  

   ! get starting values for the p regression coefficient plus
   ! the means of the random effects                 
   !
   ! use IWORKNF() and IWRKNF() as work vectors (nf = p + r)
   !     IWNFNF() as a work vector (nf * (nf+1) / 2)

   IF (ISTART .EQ. 0) THEN
      ALLOCATE(IWORKNF(NF))
      ALLOCATE(IWRKNF(NF))
      ALLOCATE(IWNFNF(NF * (NF+1)/2))
      IC = 0       
      DO L = 1,NF
         IWRKNF(L) = 0.0D0
         IWORKNF(L)  = 0.0D0
         DO L1 = 1,L
            IC = IC + 1
            IWNFNF(IC) = 0.0D0
         END DO
      END DO
   ELSE
   ENDIF
       
   ! get the sums of squared deviations about the means
   ! and work vectors for the regression coefficients

   YVAL2 = 0.0D0
   II = 1
   IHI = 1
   IF (NOMU.EQ.1) IHI=R+1
   
   DO K  = 1,NALL
      K2 = K - (II - 1)
      IR     = MOD(K2,NRP1)
      IF (IR .EQ. 0) IR = NRP1
      IC3    = NALL + NRP1 + NRP1 + IR
      IC4    = NALL + NRP1 + NRP1 + NRP1 + IR
      ALLDAT(IC4) = ALLDAT(IC4) +((ALLDAT(K) - ALLDAT(IC3))**2)
      IF (IR .EQ. 1) THEN
         YVAL  = ALLDAT(K) 
         IF (ISTART .EQ. 0) YVAL2 = YVAL2 + YVAL*YVAL
      ELSE
      ENDIF

      IF (IR .GT. IHI .AND. ISTART .EQ. 0) THEN
         IR2 = IR-IHI
         IF(IR2 .LE. NF) IWORKNF(IR2) = ALLDAT(K) 
         IF (IR .EQ. NRP1) THEN
            CALL GRMCV(IWNFNF,IWNFNF,IWORKNF,UNO,NF)
            DO L = 1,NF
               IWRKNF(L) = IWRKNF(L) + (DLOG(YVAL+1) * IWORKNF(L))
            END DO
         ELSE
         ENDIF
      ELSE 
      ENDIF

      ! for one of the Xs or Ws get the sum of square deviations by Y

      IF (MEANYX .EQ. 1 .AND. IR .EQ. NCATYX) THEN
          CATJ = .true.
          K4    =  0
          DO WHILE (CATJ)
             K4    = K4 + 1 
               ! Test for equality with FP-safe comparison  TG
             IF (FP_EQUAL(ALLDAT(K) , CODEX(K4))) THEN
                KIN1 = K4
                KIN2 = MAXXJ + K4
                ICATFQX(KIN2)=ICATFQX(KIN2)+ &
                                   ((YVAL-ICATFQX(KIN1))**2)
                CATJ = .FALSE.
             ELSE
             ENDIF
          END DO
      ELSE
      ENDIF

   END DO

   !  calculate the standard deviations

   DO K  = 1,NRP1
      IC4    = NALL + NRP1 + NRP1 + NRP1 + K
      ALLDAT(IC4) = DSQRT(ALLDAT(IC4) / DBLE(NTOT-1))

      !  see if the SD of the 1 random effect is 0
      !  for intraclass correlation calculation later

      IF (K .EQ. 2 .AND. R .EQ. 1 .AND. ALLDAT(IC4) .GT.  &
          -0.001D0 .AND. ALLDAT(IC4) .LT. 0.001D0)  &
          ICCY = ICCY + 1

   END DO
   
   IF (MEANYX .EQ. 1) THEN
      DO K4 = 1,MAXXJ
         KIN2 = MAXXJ + K4
         KIN3 = (2*MAXXJ) + K4
         ICATFQX(KIN2) = DSQRT(ICATFQX(KIN2)/ &
                         (ICATFQX(KIN3) - 1.0D0))
      END DO
   ELSE
   ENDIF

   ! write out descriptives

   WRITE(2,554) 
   554 FORMAT(1x,'MIXPREG - The program for mixed-effects poisson regression analysis',/)
   WRITE(2,55) HEAD
   55 FORMAT(1x,15A4)

   IF (r .ge. 1 .and. iunif .eq. 0) then
      WRITE(2,5564) 
      5564 FORMAT(/,1x,'Random-effects distribution: normal')
   ELSEIF (r .ge. 1 .and. iunif .eq. 1) then
      WRITE(2,5565) 
      5565 FORMAT(/,1x,'Random-effects distribution: rectangular')
   ELSEIF (r .ge. 1 .and. iunif .eq. 2) then
      WRITE(2,5566) 
      5566 FORMAT(/,1x,'Random-effects distribution: log gamma')
   ENDIF

   WRITE(2,255)
   255 FORMAT(//,1x,'Numbers of observations',/,1x,'-----------------------',/)
   IF (IWT .EQ. 0) THEN
      WRITE(2,56)NTOT
      56 FORMAT(1x,'Level 1 observations = ',i6)
      IF (R .GE. 1) THEN
      WRITE(2,562)N
      562 FORMAT(1x,'Level 2 observations = ',i6)
         WRITE(2,5610)
         WRITE(2,561)(IDNI(I),I=2,N*2,2)
         5610 FORMAT(//,1x,'The number of level 1 observations per level 2 unit are:',/)
         561 FORMAT(1x,19I4)
      ELSE
      ENDIF     
   ELSE
      WRITE(2,506)WSUM,NTOT
      506 FORMAT(1x,'Level 1 observations = ',f8.2,/,1x,'Level 1 patterns =     ',I8)
      IF (R .GE. 1) THEN
      WRITE(2,516)SUM2,N
      516 FORMAT(/,1x,'Level 2 observations = ',f8.2,/,1x,'Level 2 patterns =     ',I8)
         WRITE(2,5070)
         WRITE(2,507)(IDNI(I),I=2,N*2,2)
         5070 FORMAT(//,1x,'The number of level 1 patterns per level 2 pattern are:',/)
         507 FORMAT(1x,19I4)
      ELSE
      ENDIF     
   ENDIF


   WRITE(2,257)
   257 FORMAT(//,1x,'Descriptive statistics for all variables',/, &
                 1x,'----------------------------------------',/)
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
         WRITE(2,377)IBLABEL(H),ALLDAT(IC),ALLDAT(IC2) &
                     ,ALLDAT(IC3),ALLDAT(IC4)
      END DO
   ELSE 
   ENDIF
   
   IF (P .GE. 1) THEN
      DO L = 1,P
         IC  = IC  + 1
         IC2 = IC2 + 1
         IC3 = IC3 + 1
         IC4 = IC4 + 1
         WRITE(2,377)IALABEL(L),ALLDAT(IC),ALLDAT(IC2) &
                     ,ALLDAT(IC3),ALLDAT(IC4)
      END DO
   ELSE 
   ENDIF
   
   IF (IOFS .GE. 1) THEN
      IC  = IC  + 1
      IC2 = IC2 + 1
      IC3 = IC3 + 1
      IC4 = IC4 + 1
      WRITE(2,377)'Offset  ',ALLDAT(IC),ALLDAT(IC2) &
                  ,ALLDAT(IC3),ALLDAT(IC4)
   ELSE 
   ENDIF

   377 FORMAT(1x,A16,4(5X,F12.5))

   IF (MEANYX .EQ. 1) THEN
      WRITE(2,358)YLabel,XLabel
      358 FORMAT(//,1x,'Descriptives for response variable ',A16,' by the variable ',A16,/, &
                    1x,'--------------------------------------------------------------------',/)
      WRITE(2,372)YLabel,XLabel
      372 FORMAT(1X,8X,A16,/,1X,8X,'--------',/,1X,A16,13X,'Mean',6x,'Stand. Dev.',10x,'N',/, &
                 1x,'-----------------------------------------------------')

      MAXJXJ2 = 2*MAXXJ
      DO J = 1,MAXXJ
         WRITE(2,373)CODEX(J),(ICATFQX(J3),J3=J,MAXJXJ2,MAXXJ), &
             INT(ICATFQX(J+MAXJXJ2))
         373 FORMAT(/,1X,F8.3,2(5x,F12.5),5x,i6)
      END DO
   ELSE
   ENDIF


   ! done writing out descriptives, get starting values

   IF (ISTART .NE. 1) THEN

      ! calculate the starting values for the regression coefficients
 
      IER    = 0
      CALL INVS(IWNFNF,NF,DET,IWORKNF,IER)
      CALL MPYM(IWNFNF,IWRKNF,IWORKNF,NF,NF,1,0,1)
 
      IC3 = NALL + NRP1 + NRP1 + 1
      IC4 = NALL + NRP1 + NRP1 + NRP1 + 1
 
      ICL = 0
      GMEAN = 0.0D0
      IF (R .GE. 1 .and. IPRIOR .EQ. 0 .AND. NOMU .EQ. 0) THEN
         DO H = 1,R
            ICL = ICL+1
            IMU(H) = IWORKNF(ICL) 
            GMEAN = GMEAN + IMU(H)*IMEAN(ICL)
         END DO
      ELSE 
      ENDIF
      
      IF (P .GT. 0) THEN
         DO L = 1,P
            ICL = ICL+1
            IALPHA(L) = IWORKNF(ICL) 
            GMEAN = GMEAN + IALPHA(L)*IMEAN(ICL)
         END DO
      ELSE 
      ENDIF
 
      IF (R .EQ. 1) THEN
         ISIGMA(1) = .3d0*DSQRT(GMEAN)
      
      ELSEIF (R .GT. 1 .AND. IDIAG .EQ. 0 .AND. IVSEP .EQ. 0) THEN
         IR = 0
         DO I  = 1,R
            DO I2 = 1,I
               IR = IR + 1
               IF (I2 .EQ. I) THEN
                  IF (I .EQ. 1) THEN
                     ISIGMA(IR) = 1.0d0*DSQRT(GMEAN) 
                  ELSE
                     ISIGMA(IR) = 0.5d0*DSQRT(GMEAN) 
                  ENDIF
               ELSE
                       ISIGMA(IR) = 0.0d0
               ENDIF
            END DO
         END DO
      
      ELSEIF (R .GT. 1 .AND. (IDIAG .EQ. 1 .OR. IVSEP .EQ. 1)) THEN
         IR = 0
         DO I  = 1,R
            IR = IR + 1
            IF (I .EQ. 1) THEN
               ISIGMA(IR) = 1.0d0*DSQRT(GMEAN) 
            ELSE
               ISIGMA(IR) = 0.5d0*DSQRT(GMEAN) 
            ENDIF
         END DO
      ELSE
      ENDIF
   ENDIF

   ! write out the starting values

   WRITE(2,256)
   256 FORMAT(//,1x,'Starting values',/,1x,'---------------',/)

   IF (R .GT. 0 .AND. NOMU .EQ. 0) THEN
       WRITE(2,603) (IMU(H),H=1,R)
       603 FORMAT(1x,'mean      ',7F10.4)
       ALLOCATE(IMU0(R))
       CALL RELOC(IMU,IMU0,R,1,0)
   ELSE 
   ENDIF
   
   IF (P .GT. 0) THEN
       WRITE(2,604) (IALPHA(L),L=1,P)
       604 FORMAT(1x,'covariates',7F10.4)
       ALLOCATE(IALPHA0(P))
       CALL RELOC(IALPHA,IALPHA0,P,1,0)
   ELSE
   ENDIF
   
   IF (R .GT. 0) THEN
       WRITE(2,605) (ISIGMA(HR),HR=1,RR)
       605 FORMAT(1x,'var. terms',7F10.4)
       ALLOCATE(ISIGMA0(RR))
       CALL RELOC(ISIGMA,ISIGMA0,RR,1,0)
   ELSE
   ENDIF

   ! **********************************************************
   !               
   !  start iterations 
   !              
   ! **********************************************************

   ITLAST = 0 
   RIDGE  = 0.0D0
   IFIN   = 1
   ! NORI   = 0
   IUN    = 6
   ND     = 0
   RLOGLP = -999999999999999.0D0

   ALLOCATE(IDER(NPAR))

   IT    = 0

   ! isig   = 1 when INVS has problems (0 otherwise)
   ! ier    = 1 when INVS has problems (0 otherwise)
   ! nonpos = 0 when CHSKY detects a non pos-def matrix (1 when pos-def)

   ISIG   = 0
   IER    = 0
   NONPOS = 1
   IRBAD  = 0
   IGO = 1

   DEALLOCATE(IWIND, IXIND)
   IF(ALLOCATED(CODEX)) DEALLOCATE(CODEX)
   DEALLOCATE(ICATFQX, IWMISS, IXMISS) 
   IF(ALLOCATED(IWNFNF)) DEALLOCATE( IWNFNF, IWORKNF, IWRKNF) 
   INITIALIZED = 1
   
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
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! LOCAL VARS ONLY - SHARED VARS ARE DECLARED
   ! BELOW THE 'MODULE' STATEMENT AT THE TOP OF THE 
   ! FILE
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   INTEGER I, ICOUNT, IDI, IH, INDD, INDD2, IW, IX, &
      KK, KKK, NII, NSAME, NSAMES  
   INTEGER H2, PNII, Q, R2, RNII, TEMP_I
   
   REAL (KIND=8) AQALL, BIGCOR, DDER, DERIV, HPROB, LAM, &
      LPROB, OFFS, PRA, PSUM, QMULT, RLOGDIFF, RSAMES, &
      SCAL, SIGN, STEP, SUMW, UBOUND, WA, WTSUBI, &
      XMU, XTB, Z
   
   
   
   ! **********************************************************
   !  Perform a single iteration of the main loop
   ! **********************************************************
   
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
      ! (this assumes IALPHA follows IMU in BM),
      ! R, P, AND NPAR ACCORDINGLY
 
      IF (IRBAD .GE. 1) THEN
         IRBAD   = 0
         ISIG    = 0
         IER     = 0
         NONPOS  = 1
         ! NORI    = 0
         ! RLOGLP  = -999999999999999.0D0
         
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
            npar = p+r+rr
         ELSE
            npar = p+rr
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
            CALL QUADP(IBQ,IBQ1,IAQ,NQ1,NQ,R,IUNIF,IAQ1,IAQ2,ISIGMA(1))
         ELSE
            CALL QUADP(IBQ,IBQ1,IAQ,NQ1,NQ,1,IUNIF,IAQ1,IAQ2,ISIGMA(1))
         ENDIF
 
         ! get ORIGINAL starting values
 
         IF (r .ge. 1) then
            IF (NOMU.EQ.0) CALL RELOC(IMU0,IMU,R,1,0)
            CALL RELOC(ISIGMA0,ISIGMA,RR,1,0)
         ELSE 
         ENDIF
         
         IF (P .GT. 0) THEN

            ! note that R now equals the old R minus 1
            ! so that imu0+R equals the old imu0+R-1 
            ! which should be the last element of mu
            
            IALPHA(1) = IMU0(R+1)
            CALL RELOC(IALPHA0,IALPHA(2),P-1,1,0)
         ELSE
         ENDIF
 
      ELSE
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
      ELSE
      ENDIF
 
      ! ***
      ! LIKELIHOOD INITIALIZATION
      ! ***
      !                                                    */
      !  calculate and save the conditional likelihood - ILIK()   */
      !                     the marginal likelihood    - hprob */
      !                     the log-likelihood value   - RLOGL */
      !                                                        */ 
      CALL gen(IDER2,0.0d0,npar,npar,1)
      CALL gen(IDER,0.0d0,npar,1,0)
      RLOGL = 0.0d0
 
      IF (IPRIOR .EQ. 1) THEN
         AQALL    = 0.0D0
         DO Q = 1,NQ
            IAQEMP(Q) = 0.0D0
         END DO
      ELSE
      ENDIF
 
      ! *****************************************************************
      ! GO OVER ALL SUBJECTS
      !    IUN = UNIT NUMBER FOR OUTPUT (0 = SCREEN)
      ! *****************************************************************
      ICOUNT = 0
 
      IF (ALLOCATED(IDERP)) THEN
         IF(UBOUND(IDERP,1) .NE. NPAR) THEN
            DEALLOCATE(IDERP)
            ALLOCATE(IDERP(NPAR))
         ENDIF
      ELSE
         ALLOCATE(IDERP(NPAR))
      ENDIF

      IF (ALLOCATED(IDERP2)) THEN
         IF(UBOUND(IDERP2,1) .NE. NPARR) THEN
            DEALLOCATE(IDERP2)
            ALLOCATE(IDERP2(NPARR))
         ENDIF
      ELSE
         ALLOCATE(IDERP2(NPARR))
      ENDIF

      DO I = 1,N   ! loop once per subject
    
         hprob = 0.0d0
         CALL gen(IDERP,0.0d0,npar,1,0)
         CALL gen(IDERP2,0.0d0,npar,npar,1)
    
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
     
         ! calculate the number of level-2 units with no variance 
         ! in terms of level-1 units Y values
    
         IF (IT .EQ. 1) THEN
            IF (I  .EQ. 1) NSAMES = 0
            CALL YSAME(IYI,NII,NSAME)
            NSAMES = NSAMES  + NSAME
            IF (I  .EQ. N .AND. R .GT. 0) THEN
               RSAMES = 100.0D0 * (DBLE(NSAMES) / DBLE(N))
               IF (IWT.EQ.0) THEN
               WRITE(2,508)NSAMES,RSAMES
               508 FORMAT(//,1x,'==> The number of level 2 observations with non-varying responses',/, &
                  5x,'= ',I6,' ( ',F6.2,' percent )')  
               ELSE
               WRITE(2,509)NSAMES,RSAMES
               509 FORMAT(//,1x,'==> The number of level 2 patterns with non-varying responses',/, &
                  5x,'= ',I6,' ( ',F6.2,' percent )')  
               ENDIF
            ELSE
            ENDIF
         ELSE
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
         ELSE
         ENDIF
    
   
         ! THE OFFSET VARIABLE variable OV(K) vector  K = 1 .. NI(I) 
   
         IF (IOFS .GE. 1) THEN
            IF (ALLOCATED(IOV)) THEN
               IF(UBOUND(IOV,1) .NE. NII) THEN
                  DEALLOCATE(IOV)
                  ALLOCATE(IOV(NII))
               ENDIF
            ELSE
               ALLOCATE(IOV(NII))
            ENDIF
            
            IC = 0
            DO K  = 1,NII
               IC           = IC + 1
               IC2          = ICOUNT + (NRP1 * (K-1) + P+R+2)
               IOV(IC) = ALLDAT(IC2)
            END DO
         ELSE
         ENDIF
    
         IF (IOFS .EQ. 0) THEN
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
            WRITE(IUN,11)N
            11 FORMAT(1X,'TOTAL NUMBER OF SUBJECTS = ',I6,/)
            WRITE(IUN,12) IDI,NII
            12 FORMAT(1X,'DATA FOR SUBJECT ',I10,' WHO HAS',I6, &
                   ' OBSERVATIONS')
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
            
            IF (IOFS .GT. 0) THEN
               CALL PRNT(IUN,IOV,NII,1,0,IXLAB,IXLAB,ND,HEAD,1,80,5,1,1, &
                  60H OFFSET VARIABLE                                            )
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
               DEALLOCATE(IDZ)
               ALLOCATE(IDERQ(NPAR))
               ALLOCATE(IDZ(NPAR))
            ENDIF
         ELSE
            ALLOCATE(IDERQ(NPAR))
            ALLOCATE(IDZ(NPAR))
         ENDIF
         
         IF (ALLOCATED(ILIK)) THEN
            IF(UBOUND(ILIK,1) .NE. NQ) THEN
               DEALLOCATE(ILIK)
               ALLOCATE(ILIK(NQ))
            ENDIF
         ELSE
            ALLOCATE(ILIK(NQ))
         ENDIF
         
         IF (ALLOCATED(IDERQ2)) THEN
            IF(UBOUND(IDERQ2,1) .NE. NPARR) THEN
               DEALLOCATE(IDERQ2)
               ALLOCATE(IDERQ2(NPARR))
            ENDIF
         ELSE
            ALLOCATE(IDERQ2(NPARR))
         ENDIF
         
         DO q=1,nq 
    
            ! for log gamma get new weights - they depend on sigma
       
            IF (IUNIF .EQ. 2) THEN
               CALL QLGAM(IBQ,IAQ,NQ,ISIGMA(1))
            ENDIF
       
            PSUM          = 0.0D0
            ILIK(Q)  = 0.0d0
            CALL gen(IDERQ,0.0d0,npar,1,0)
            CALL gen(IDERQ2,0.0d0,npar,npar,1)
       
            DO K=1,NII 
       
               WA = 0.0D0
               IF (P .GE. 1) THEN
                  DO l=1,p
                     iw = k  + (l-1)*nii
                     wa = wa + IALPHA(L)*IWI(iw)
                  END DO
               ELSE
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
                     CALL MPYM(ISIGMA,IWORKR,IWORKR,R,R,3,0,1)
                  ELSE
                     CALL MPYM(ISIGMA,IWORKR,IWORKR,R,R,2,0,1)
                  ENDIF
                  
                  DO h=1,r
                     ix  = k   + (h-1)    * nii
                     IF (NOMU .EQ. 0) THEN
                        xmu = xmu + IMU(H)    *IXI(ix)
                     ELSE
                        XMU = 0.0D0
                     ENDIF
                     xtb = xtb + IWORKR(H) * IXI(ix)
                  END DO
               ELSE
               ENDIF
          
               IF (IOFS .EQ. 0) THEN
                  OFFS = 0.0D0
               ELSEIF (IOFS .EQ. 1) THEN
                  OFFS = DLOG(IOV(K))
               ENDIF
       
               Z     = XMU + WA + XTB + OFFS
               LAM   = DEXP(Z)
               LPROB = IYI(K)*Z - LAM - DLOGFACT(IDINT(IYI(K)))
               PSUM  = PSUM + LPROB
       
               IND = 0
            
               DERIV = IYI(K) - LAM
       
               IF (R .GE. 1 .and. IPRIOR .EQ. 0 .AND. NOMU .EQ. 0) THEN
                  DO h=1,r
                     IND = IND+1
                     ix  = k + (h-1) * nii
                     IDZ(IND)   = IXI(ix)
                     IDERQ(IND) = IDERQ(IND)+deriv*IXI(ix)
                  END DO
               ENDIF
       
               IF (P .GE. 1) THEN
                  DO l = 1,p
                     IND= IND + 1
                     iw = k + (l-1)*nii
                     IDZ(IND)   = IWI(iw)
                     IDERQ(IND) = IDERQ(IND) + deriv*IWI(iw)
                  END DO
               ENDIF
       
               IF (R .GE. 1) THEN
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
                        IH  = IH + 1
                        IND  = IND + 1
                        IDERQ(IND) = IDERQ(IND) + deriv * IWORKR2(IH)
                     ELSE
                     
                        DO h2 = 1,h
                           IH  = IH  + 1
                           IND = IND + 1
                           
                           ! estimate sigma 
                           IF (IFIN .NE. 0 .OR. H .NE. H2) THEN
                               IDZ(IND)   = IWORKR2(IH)
                               IDERQ(IND) = IDERQ(IND) + deriv * IWORKR2(IH)
                           ELSE
                           ! estimate tau = ln(diagonal elements of sigma)
                           IF (DERIV .LT. 0.0D0 .AND. IWORKR2(IH) .LT. 0.0D0) &
                               SIGN =  1.0D0
                           IF (DERIV .LT. 0.0D0 .AND. IWORKR2(IH) .GE. 0.0D0) &
                               SIGN = -1.0D0
                           IF (DERIV .GE. 0.0D0 .AND. IWORKR2(IH) .LT. 0.0D0) &
                               SIGN = -1.0D0
                           IF (DERIV .GE. 0.0D0 .AND. IWORKR2(IH) .GE. 0.0D0) &
                               SIGN =  1.0D0
                               IDZ(IND)=IWORKR2(IH)*SIGN*ISIGMA(IH)
                               IDERQ(IND) = IDERQ(IND) + SIGN *  &
                                          DEXP(DLOG(DABS(deriv)) &
                                         +DLOG(ISIGMA(IH)) &
                                         +DLOG(DABS(IWORKR2(IH))))
                           ENDIF
                        END DO
                     ENDIF
                  END DO
               ELSE
               ENDIF
       
               !  CALL PRNT(IUN,IDZ,NPAR,1,0,IXLAB,IXLAB,ND,
               ! +HEAD,1,80,5,1,1,
               ! +60H DZ  derivatives                                            )
      
               SCAL = 0.0D0 + LAM
               CALL GRMCV(IDERQ2,IDERQ2,IDZ,SCAL,NPAR)
       
            END DO   ! K=1,NII
       
            ILIK(Q) = DEXP(PSUM)
            PRA =  DEXP ( PSUM + DLOG(IAQ(Q)))
            hprob  = hprob + PRA
            IF (IPRIOR .EQ. 1) IAQEMP(Q) = IAQEMP(Q) + DEXP( (DLOG(WTSUBI)) &
                               + (DLOG(ILIK(Q))) + (DLOG(IAQ(Q))))
       
            DO L  = 1,npar
               IDERP(L) = IDERP(L) + IDERQ(L)*PRA
            END DO
       
            ! CALL PRNT(IUN,IDERQ2,NPAR,NPAR,1,IXLAB,IXLAB,ND,
            !+HEAD,1,80,5,1,1,
            !+60H 2ND derivatives before adding 2nd part                     )
            
            CALL GRMCV(IDERQ2,IDERQ2,IDERQ,1.0D0,NPAR)
            
            ! CALL PRNT(IUN,IDERQ2,NPAR,NPAR,1,IXLAB,IXLAB,ND,
            !+HEAD,1,80,5,1,1,
            !+60H 2ND derivatives after  adding 2nd part                     )
            
            CALL SCM(IDERQ2,PRA,IDERQ2,NPAR,NPAR,1)
            CALL ADDM(IDERP2,IDERQ2,IDERP2,NPAR,NPAR,1)      
       
         END DO ! do q=1,nq 
    
         IF (HPROB .LT. .1D-305) HPROB = .1D-305
         RLOGL = RLOGL + WTSUBI*DLOG(hprob)
         IF (IPRIOR .EQ. 1) AQALL = AQALL + DEXP( (DLOG(WTSUBI)) +  &
                                       (DLOG(HPROB)))
    
         scal = DEXP( 0.0d0 - DLOG(hprob))
         
         DO l = 1,npar
            IDER(L) = IDER(L) + IDERP(L)* DEXP(DLOG(WTSUBI)-DLOG(hprob))
            sign=1.0d0
            IF (IDERP(L) .LT. 0) SIGN=-1.0D0
            DDER = DABS(IDERP(L)) 
            IF (dder .LT. .1D-305)  dder  = .1D-305
            IDERP(L) = SIGN*(DEXP(DLOG(DDER) - DLOG(hprob)))
            ! WRITE(6,*)' L = ',L,'D = ',IDERP(L)
         END DO
    
         ! CALL PRNT(IUN,IDERP2,NPAR,NPAR,1,IXLAB,IXLAB,ND,
         !+HEAD,1,80,5,1,1,
         !+60H 2ND derivatives before SCAL                                )
   
         CALL SCM(IDERP2,SCAL,IDERP2,NPAR,NPAR,1)
         
         ! CALL PRNT(IUN,IDERP2,NPAR,NPAR,1,IXLAB,IXLAB,ND,
         !+HEAD,1,80,5,1,1,
         !+60H 2ND derivatives after SCAL                                 )
         ! CALL SCM(IDERP2,HPROB,IDERP2,NPAR,NPAR,1)
   
         CALL grmcv(IDERP2,IDERP2,IDERP,-1.0D0,npar)
         
         ! scal       = DEXP( 0.0d0 - DLOG(hprob) - DLOG(hprob))
         ! CALL SCM(IDERP2,SCAL,IDERP2,NPAR,NPAR,1)
         ! CALL grmcv(IDERP2,IDERP2,IDERP,-1.0D0,npar)
   
         CALL ADDM(IDER2,IDERP2,IDER2,NPAR,NPAR,1)      
    
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
    
      END DO   ! loop once per subject   DO I = 1,N  
  
      ! ****************************************************
      !  done with subjects                                *
      ! ****************************************************
 
      !  CALL PRNT(IUN,IDER,NPAR,1,0,IXLAB,IXLAB,ND,HEAD,1,
      ! +80,5,1,1,
      ! +60H derivatives                                                )
      !  CALL PRNT(IUN,IDER2,NPAR,NPAR,1,IXLAB,IXLAB,ND,
      ! +HEAD,1,80,7,1,1,
      ! +60H 2ND derivatives                                            )
 
 
      IF (ALLOCATED(ICOREC)) THEN
         IF(UBOUND(ICOREC,1) .NE. NPAR) THEN
            DEALLOCATE(ICOREC)
            ALLOCATE(ICOREC(NPAR))
         ENDIF
      ELSE
         ALLOCATE(ICOREC(NPAR))
      ENDIF

      ! ADD RIDGE TO DIAGONAL ELEMENTS OF DER2 PRIOR TO INVERSION
      ! (unless at final iteration itlast=1)
                
      IF(ITLAST.EQ.0)THEN
         DO KKK=1,NPAR
            KK=(KKK*(KKK+1))/2
            IDER2(KK)=IDER2(KK)+RIDGE*IDER2(KK)
         ENDDO
      ENDIF
 
      ! check to see if the matrix of second derivatives is 
      ! positive-definite (nonpos = 1) or not positive definite 
      ! (nonpos = 0)
 
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
         DO h = 1,r
            IND = IND+1
            IMU(IND)    = IMU(IND) + ICOREC(IND)
            IF (DABS(ICOREC(IND)) .GT. BIGCOR) BIGCOR =  &
                DABS(ICOREC(IND))
         END DO
      ELSE
      ENDIF   
 
      IF (p .ge. 1) then
         IND2 = 0
         DO l = 1,p
            IND = IND+1
            IND2 = IND2+1
            IALPHA(IND2) = IALPHA(IND2) + ICOREC(IND)
            IF (DABS(ICOREC(IND)) .GT. BIGCOR)  &
                BIGCOR = DABS(ICOREC(IND))
         END DO
      ENDIF
 
      IF (r .ge. 1) then
         IND2 = 0
         INDD = 1
         INDD2= 0
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

               ! make sure that the first variance term is 
               ! positive (COMMENT OUT)
               IF (IFIN .NE. 0) THEN
                  IF (H.EQ.1 .AND. RR .EQ. 1 .AND.  &
                  (0.0D0-ICOREC(IND) .GT. ISIGMA(IND2)))THEN
                  ISIGMA(IND2) = 0.10D0* DABS(ISIGMA(IND2) &
                                 + ICOREC(IND))
                  ELSE
                  ISIGMA(IND2) = ISIGMA(IND2) &
                                 + ICOREC(IND)
                  ENDIF

                  ! DON'T shift to estimation of tau if sigma gets small 
                  ! INCREASE THE RIDGE ONCE INSTEAD
                  
                  ! IF(ISIGMA(H) .LE. .1D0 .AND. NORI .EQ. 0) THEN
                  !    IFIN=0
                  !    NORI=NORI+1
                  !    RIDGE = RIDGE +.1D0
                  !    WRITE(2,869)
                  ! 869   FORMAT(///,1X,'==> Reparameterized Estimation Occurred')
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
            IF (DABS(ICOREC(IND)) .GT. BIGCOR) BIGCOR =  &
                DABS(ICOREC(IND))
         END DO  ! DO h = 1,rr
      ELSE 
      ENDIF
 
      ! ******
      ! PRINT OUT RESULTS FOR THIS ITERATION
      ! unless iteration was for computation of information matrix only
      ! (ITLAST = 1)
      ! ******
 
      IF (ITLAST .EQ. 1) THEN
         MainLoop = 0
         CALL FREE_LOCALS()
         INITIALIZED = 2
         RETURN     ! escape from the loop
      ENDIF
 
      WRITE(IUN,65)IT
      65 FORMAT(/,1X,'**************',/,1X,'ITERATION ',I4,/, &
                1X,'**************',/)
 
      WRITE(IUN,66)RLOGL
      66 FORMAT(1X,'Log Likelihood   = ',F12.3,//)
   
      ! save current RLOGL to check at next iteration
      RLOGDIFF = RLOGL - RLOGLP
      RLOGLP   = RLOGL

      ! calculate versions of RLOGL 
      AIC  = RLOGL - RNPAR
      SBC  = RLOGL - 0.5 * RNPAR * DLOG(DBLE(N))
      DEV  = -2*RLOGL
      AICD = -2*(RLOGL-RNPAR)
      SBCD = -2*(RLOGL-0.5*RNPAR * DLOG(DBLE(N)))

      IF (r .ge. 1 .AND. NOMU .EQ. 0) then
         WRITE(IUN,903) (IMU(H),H=1,R)
         903 FORMAT(1X,'mu         ',7F10.6)
      ELSE 
      ENDIF
      
      IF (P .GT. 0) THEN
         WRITE(IUN,904) (IALPHA(L),L=1,P)
         904 FORMAT(1X,'alpha      ',7F10.6)
      ELSE
      ENDIF
      
      IF (r .ge. 1) then
         WRITE(IUN,905) (ISIGMA(HR),HR=1,RR)
         905 FORMAT(1X,'sigma      ',7F10.6)
      ELSE 
      ENDIF
      
      WRITE(IUN,907) (ICOREC(L),L=1,NPAR)
      907 FORMAT(1X,'corrections',7F10.6)
 
      ! check to see if the log-likelihood doesn't increase
      ! and if so then increase the ridge by 0.1 each time 
      ! check to see if there are numerical problems and if so
      ! reduce the number of random-effects by 1
 
      IF(RLOGDIFF .LT. 0.0D0) RIDGE = RIDGE + 0.1D0
 
      IF(IER .EQ. 1 .OR. NONPOS .EQ. 0 .OR. ISIG .EQ. 1) IRBAD = IRBAD+1
      
      IF(IRBAD .EQ. 1) THEN
         WRITE(6,878)IT
         WRITE(2,879)IT
         878 FORMAT(///,1X,'==> WARNING!  Estimation Difficulties at Iteration',I4,/, &
                        1x,'==> will proceed with one less random effect')
         879 FORMAT(///,1X,'==> WARNING!  Estimation Difficulties Occurred at Iteration',I4,/, &
                        1x,'==> The model was fit with one less random effect than was requested')
      ELSE
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

      ELSE
      ENDIF
 

      ! CHECK IF CONVERGENCE HAS BEEN REACHED
      ! IF NOT RETURN TO ITERATIONS 
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
      ! free up arrays only used in the main loop, now that 
      ! we are finished with them.  Note that the SAVE statement
      ! preserves the values in these arrays between calls to
      ! the MainLoop function
      CALL FREE_LOCALS()
      INITIALIZED = 2
   ENDIF
END FUNCTION MainLoop
   
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
   
   DEALLOCATE(ICHWRK, ICOREC, IDERP, IDERP2, &
      IDERQ, IDERQ2, IDZ, ILIK)
   
   IF(ALLOCATED(IOV)) DEALLOCATE(IOV)
   IF(ALLOCATED(IWI)) DEALLOCATE(IWI)
   IF(ALLOCATED(IWRKR)) DEALLOCATE(IWRKR,IWORKR2)
   IF(ALLOCATED(IXI)) DEALLOCATE(IXI)
   IF(ALLOCATED(IYI)) DEALLOCATE(IYI)

END SUBROUTINE FREE_LOCALS
   
   
   
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
   
   REAL (KIND=8) CORR, COV, SIDE, VAR1, VAR2
   
   INTEGER I, INP, INP0, J2, J4, NCONN, Q
   
   CHARACTER (LEN=8), ALLOCATABLE :: IPARMLAB(:)
   CHARACTER*10 BLAB   
   CHARACTER*4 BLANK(40)

   REAL (KIND = 8),POINTER ::ICS(:) 
   
   REAL (KIND=8), ALLOCATABLE :: ICDER2(:), ICPVAL(:), ICSE(:), ICVAL(:), &
      ICZVAL(:), IDER2SQ(:), IPVAL(:), ISE(:), IVAL(:), IZVAL(:) 

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
  
   IF(INITIALIZED .NE. 2) THEN
      CALL POST_ERROR("Cleanup() routine called out of sequence in " // PROGNAME)
      RETURN
   ENDIF
   
   ALLOCATE(ISE(NPAR))
   ALLOCATE(IVAL(NPAR))
   ALLOCATE(IZVAL(NPAR))
   ALLOCATE(IPVAL(NPAR))
   IF (NCON .GT. 0) THEN
      ALLOCATE(ICVAL(NCON))
      ALLOCATE(ICSE(NCON))
      ALLOCATE(ICZVAL(NCON))
      ALLOCATE(ICPVAL(NCON))
      NCONN = (NCON * (NCON+1)) / 2
      ALLOCATE(ICDER2(NCONN))
   ENDIF

   ALLOCATE(IPARMLAB(NPAR))

      open(29,file="mixpreg.lik")
      write(29,*) rlogl,npar
      close(29)
      
   WRITE(2,455)
   455 FORMAT(///,1x,'---------------------------------------------------------',/, &
                  1x,'* Final Results - Maximum Marginal Likelihood Estimates *',/, &
                  1x,'---------------------------------------------------------',/)
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
   
   WRITE(2,57)
   57 FORMAT(/,1X,'Variable',5x,'    Estimate',5X,'Stand. Error',5x,    &
      '           Z',5X,'     p-value',/,1X,'--------',5x, &
      '------------',5X, '------------',5X,'------------',5X,'------------')
   
   ! * Find the z-statistics for all terms  ZFIX and ZVAR
   ! * and their associated p-values        PFIX and PVAR

   IR = 0
   GMEAN = 0.0D0
   IF (R .GT. 0 .AND. IPRIOR .EQ. 0 .AND. NOMU .EQ. 0) THEN
      DO H = 1,R
         IR = IR+1
         GMEAN = GMEAN  + IMU(IR)*IMEAN(IR)
         IVAL(IR)  = IMU(IR)
      END DO
   ELSE
   ENDIF
   
   IF (P .GT. 0) THEN
      L1 = 0
      DO L = 1,P
         IR = IR + 1
         L1 = L1 + 1
         GMEAN = GMEAN  + IALPHA(L1)*IMEAN(IR)
         IVAL(IR) = IALPHA(L1)
      END DO
   ELSE
   ENDIF
   
   ! add in the mean of the offset variable
   IF (IOFS .EQ. 1) GMEAN = GMEAN + DLOG(IMEAN(IR+1))

   IF (R .GT. 0) THEN
      L1 = 0
      DO HR= 1,RR
         ir  = ir + 1
         L1  = L1 + 1
         IVAL(IR) = ISIGMA(L1)
      END DO
   ENDIF

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

   IF (R .GT. 0 .AND. NOMU .EQ. 0) THEN
      DO H = 1,R
         INP = H
         IF (IPRIOR .EQ. 0) THEN
            WRITE(2,803)IBLABEL(H),IMU(INP),ISE(INP), &
                        IZVAL(INP),IPVAL(INP),SIDEL(2)
         ELSE
            WRITE(2,803)IBLABEL(H),IMU(INP)
         ENDIF
         803 FORMAT(1X,A16,3(5x,F12.5),F12.5,1X,A4)
      END DO
   ELSE 
   ENDIF
   
   IF (P .GT. 0) THEN
      DO L = 1,P
         IF (NOMU .EQ. 0) THEN
            INP  = R + L 
         ELSE
            INP  = L 
         ENDIF
         INP0 = L 
         WRITE(2,804)IALABEL(L),IALPHA(INP0),ISE(INP), &
                      IZVAL(INP),IPVAL(INP),SIDEL(2)
         804 FORMAT(1X,a16,3(5x,F12.5),F12.5,1X,A4)
      END DO
   ENDIF

   IF (r .gt. 0) then
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
            SIDET = SIDEL(2)
         ENDIF
         
         IF (NOMU .EQ. 0) THEN
            INP  = R+P + HR 
            INP0 = HR 
         ELSEIF (NOMU .EQ. 1) THEN
            INP  = P +  HR 
            INP0 = HR 
         ENDIF
         ! IF (IPRIOR .EQ. 0) THEN
            IPVAL(INP) = IPVAL(INP)*side
            WRITE(2,805)Blab,ISIGMA(INP0),ISE(INP), &
                        IZVAL(INP),IPVAL(INP),SIDET
         ! ELSE
         !    WRITE(2,805)BLab,ISIGMA(HR)
         ! ENDIF
         805 FORMAT(1X,a10,3x,f12.5,2(5x,F12.5),F12.5,1X,A4)
      END DO
   ELSE
   ENDIF


   WRITE(2,587)
   587 FORMAT(//,1x,'note: (1) = 1-tailed p-value',/, &
               '       (2) = 2-tailed p-value')

   IF (R .EQ. 1 .AND. ICCY .EQ. 2 .AND. IUNIF .EQ. 0) THEN
      !        WRITE(2,811)
      ! 811    FORMAT(///,1X,'Calculation of the intracluster correlation',/,
      !    +1x,'-------------------------------------------')
      !        sigsd  = ISIGMA(1)
      !        sigvar = ISIGMA(1)*ISIGMA(1)
      !
      !           CORR   = sigvar / (sigvar + DEXP(GMEAN))
      !        WRITE(2,711)DEXP(GMEAN)
      ! 711    FORMAT(1x,'residual variance = ',F6.3,' (exp of grand mean)')
      !        WRITE(2,611)sigsd,sigsd,sigvar,sigvar,sigvar,DEXP(gmean),CORR
      ! 611    FORMAT(1x,'cluster  variance = (',f6.3,' * ',f6.3,
      !    +' ) = ',f6.3,/,1x,'intracluster correlation = ',f6.3,' / (',f6.3,
      !    +' + ',f6.3,' ) = ',F5.3)

   ELSEIF (R .EQ. 2 .AND. IDIAG .EQ. 0 .AND. IVSEP .EQ. 0) THEN
      WRITE(2,812) 
      812 FORMAT(///,1X,'Calculation of the random effects variance-covariance matrix',/, &
                     1x,'------------------------------------------------------------')
      ICS  => ISIGMA 
      VAR1 = ICS(1)*ICS(1)
      COV  = ICS(1)*ICS(2)
      VAR2 = ICS(2)*ICS(2) + ICS(3)*ICS(3)
      CORR = COV / (DSQRT(VAR1) * DSQRT(VAR2))
      WRITE(2,813)IBLABEL(1),ICS(1),ICS(1),VAR1,ICS(1),ICS(2),COV,IBLABEL(2), &
             ICS(2),ICS(2),ICS(3),ICS(3),var2,corr
      813 FORMAT(1x,A16,' variance = ( ',f6.3,' * ',f6.3,' ) = ',f6.3,/, &
         1x,7x,'covariance = ( ',f6.3,' * ',f6.3,' ) = ',f6.3,/, &
         1x,a8,' variance = ( ',f6.3,' * ',f6.3,' ) + ( ',f6.3,' * ',f6.3, &
         ' ) = ',f6.3,//,1x,'Covariance expressed as a correlation = ' &
         ,f6.3,//)
   ENDIF

   I = 0
   IF (R .NE. 0 .AND. NOMU .NE. 1) THEN
      DO H=1,R
         I = I+1
         ISE(I) = 1.0D0 / ISE(I)
         IPARMLAB(I) = IBLABEL(H)
      END DO
   ENDIF
   
   IF (P .NE. 0) THEN
      DO L=1,P
         I = I+1
         ISE(I) = 1.0D0 / ISE(I)
         IPARMLAB(I) = IALABEL(L)
      END DO
   ENDIF

   IF (R .NE. 0) THEN
      DO H=1,RR
         I = I+1
         ISE(I) = 1.0D0 / ISE(I)
         IPARMLAB(I) = 'VarCov' // CH(H)
      END DO
   ENDIF
  
   DO I = 1,30
      BLANK(I) = '    '
   END DO
   
   ND = 0

   ! write out the parameter estimates  
   ! the estimated variance-covariance estimates of these estimates

   DO I = 1,Npar
      IZVAL(I) = IZVAL(I) / ISE(I)
      WRITE(3,704)IPARMLAB(I),IZVAL(I)
      704 FORMAT(1X,a16,1x,F13.7)
   END DO
   write(*,*)
   
   ALLOCATE(IDER2SQ(NPAR*NPAR))
   CALL CHAMS(IDER2,IDER2SQ,NPAR,1,0)
   
   DO I = 1,Npar
      WRITE(4,705)(IDER2SQ(I2),I2=I,NPAR*NPAR,NPAR)
      705 FORMAT(1X,6F13.7)
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
               ICPVAL(I) = 2.0D0 *(1.0D0 - PHIFN(DABS(ICZVAL(I)),0))
            ELSE
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
         '          Z',5X,'     p-value',/,1X,'---------',5x,'------------',5X &
         ,'------------',5X,'------------',5X,'------------')

      DO L = 1,NCON
         WRITE(2,840)L,ICVAL(L),ICSE(L), ICZVAL(L),ICPVAL(L)
         840 FORMAT(1X,I8,4(5x,F12.5))
         ICSE(L) = 1.0D0 / ICSE(L)
      END DO

      WRITE(2,589)
      589 FORMAT(/,1x,'note: p-values are 2-tailed')

      CALL MPDSD(ICSE,ICDER2,ICDER2,NCON)
      CALL PRNT(2,ICDER2,NCON,NCON,1,IXLAB,IXLAB,ND,BLANK,1,78,5,1,1, &
            60HCorrelation of the MML Transformed Estimates                )

   ENDIF


   ! print out the empirical prior

   IF (IPRIOR .EQ. 1) THEN

      WRITE(2,656)
      656 FORMAT(//,1x,'Empirical Prior Distribution',/, &
                    1x,'----------------------------',/)
      IF (r .eq. 1) THEN
         WRITE(2,657)
         657 FORMAT(/,1X,' #',5X,'    Node',5x,'   Weight',/,1x,'--', &
                6x,'-------',6x,'--------')
         DO Q = 1,NQ
            WRITE(2,659) Q,IBQ(Q),IAQ(Q)
            659 FORMAT(1x,I2,5X,F8.5,5X,F8.5)
         END DO
         
      ELSEIF (r .eq. 2) THEN 
         WRITE(2,637)(IBQ1(Q),Q=1,NQ1)
         637 FORMAT(/,1X,'   Node',10F7.4,/)
         DO J  = 1,NQ1
            J2 = (J-1)*NQ1 +  1
            J3 = (J-1)*NQ1 + 10
            WRITE(2,644)IBQ1(J),(IAQ(J4),J4=J2,J3)
            WRITE(2,644)IBQ1(J),(IBQ(J4),J4=J2,J3)
            644 FORMAT(1X,11F7.4)
         END DO
      ENDIF

   ELSE 
   ENDIF

   CLOSE(2)
   CLOSE(5)

   ! Deallocate local arrays
   DEALLOCATE(ISE, IVAL, IZVAL, IPVAL, IPARMLAB, IDER2SQ)
   IF(ALLOCATED(ICVAL)) DEALLOCATE( ICVAL, ICSE, ICZVAL, ICPVAL, ICDER2)

   ! Deallocated global arrays in preparation for starting over
   DEALLOCATE( ALLDAT, IDNI )
   DEALLOCATE( IBLABEL, IALABEL, IXLAB)
   
   DEALLOCATE( IWORKR, ISIGTAU, ISIGMA, IBQ1, ICON, &
      IBQ, IMU, IALPHA, IMEAN )
   
   DEALLOCATE(IAQ, IAQEMP, IAQ1, IAQ2, IDER, IDER2)
   IF(ALLOCATED(IALPHA0)) DEALLOCATE(IALPHA0)
   IF(ALLOCATED(IMU0)) DEALLOCATE(IMU0)
   IF(ALLOCATED(ISIGMA0)) DEALLOCATE(ISIGMA0)
   INITIALIZED = 0
   close(3)
   close(4)
   
END SUBROUTINE CLEANUP


! *****************************************************************
! LOG FACTORIAL FUNCTION                                 
! *****************************************************************
REAL (KIND=8) FUNCTION DLOGFACT(N)
   INTEGER I,N

   DLOGFACT = 0.0D0
   DO I = 2,N
      DLOGFACT = DLOGFACT + DLOG(DBLE(I))
   END DO

   RETURN
END FUNCTION DLOGFACT


! *************************************************************
! SUBROUTINE QLGAM                         
!      PURPOSE: CALCULATE QUADRATURE WEIGHTS FOR LOG-GAMMA
! ************************************************************
SUBROUTINE QLGAM(X,A,NN,SCA)
   IMPLICIT REAL*8(A-H,O-Z),INTEGER*4(I-N)
   DIMENSION X(*),A(*)

   GSCA = DEXP(0.0D0-SCA)*DEXP((SCA*0.5D0)*DLOG(SCA))* &
          DSQRT(2.0D0*3.141592654D0)
   SUM = 0.0D0
   DO I = 1,NN
      A(I) = (DEXP(SCA*X(I)) * DEXP(0.0D0 - DEXP(X(I))) ) / GSCA
      SUM = SUM + A(I)
   ENDDO 

   DO I = 1,NN
      A(I) = A(I) / SUM                              
   END DO 
   RETURN
END SUBROUTINE QLGAM

end program