!  **************************************************************
!   MIXREG        
!
!  modified on 7/10/00 to add NOMU NOCOV options
!
!  modified on 7/24/00 to add NCON option             
!
!  modified on 9/10/02 to correct AIC & BIC statistics
! *************************************************************
!  RRM MAIN PROGRAM                                       
!                                                         
!  Model                                                  
!  Y = X BETA + W ALPHA + OMEGA ERROR                     
!                                                         
!  Indices                                                
!     N  = TOTAL NUMBER OF SUBECTS         I = 1 .. N     
!     R  = DIMENSION OF RANDOM BETA SPACE  H = 1 .. R     
!     P  = DIMENSION OF FIXED ALPHA SPACE  L = 1 .. P     
!  NI(I) = NUMBER OF OBSERVED TIMEPTS      K = 1 .. NI(I) 
!          WITH MAXNI = MAXIMUM OF NI(I)                  
!     S  = NUMBER OF AUTOCORR TERMS       IS = 1 .. S     
!          WITH MAXIMUM = MAXNI                           
!   RNI  = R x MAXNI - MAX SIZE OF X VECTOR               
!   PNI  = P x MAXNI - MAX SIZE OF W VECTOR               
!  RNII  = SIZE OF EACH SUBJ' X VECTOR    HR = 1 .. RNII  
!  PNII  = SIZE OF EACH SUBJ' W VECTOR    LP = 1 .. PNII  
!   NINI = SIZE OF EACH SUBJ OMEGA VECTOR IE = 1 .. NINI  
!          NINI = (NI(I) * NI(I)+1) / 2                   
!          WITH MIMI  = MAXIMUM OF NINI                   
!                                                         
!  Parameters                                             
!  NI  = VECTOR (N) WITH NUMBER OF OBS FOR EACH SUBJECT   
!  TIME= VECTOR (MAXNI) WITH TOTAL POSSIBLE TIMEPOINTS    
!  ID  = VECTOR (N) WITH EACH SUBJECT'S ID                
!  Y   = MATRIX (N, MAXNI) DEPENDENT VARIABLE             
!  X   = MATRIX (N, RNI) DESIGN FOR RANDOM EFFECTS        
! BETA = MATRIX (N, R) OF RANDOM EFFECTS                  
!  W   = MATRIX (N, PNI) FIXED COVARIATES                 
! XTI  = VECTOR (NI(I)) OF TIMEPOINTS (LINEAR) FOR SUBJ I 
! ALPHA= VECTOR (P) OF FIXED EFFECTS                      
! OMEGA= MATRIX (N, MIMI) DESIGN FOR ERROR                
! ERROR= SCALAR ERROR                                     
!  MU  = VECTOR (R) OF MEAN BETA VALUES                   
! VARCO= VECTOR (RR) OF BETA VARIANCE COVARIANCE          
!        WHERE HR = 1 .. RR    RR = (R * (R+1)) / 2       
! VARPO= VECTOR (RR) OF POSTERIOR BETA VARIANCE COVAR     
!        WHERE HR = 1 .. RR    RR = (R * (R+1)) / 2       
!                                                         
! *************************************************************

MODULE RRM_MODULE ! must have same name in all 4 programs
IMPLICIT NONE

   REAL (KIND=8), POINTER:: ALLDAT(:)
   INTEGER ,POINTER ::IDNI(:)
   
   CHARACTER (LEN=4), ALLOCATABLE :: IXLAB(:)
   CHARACTER(LEN=4) :: HEAD(36)
   
   CHARACTER (LEN=16), ALLOCATABLE :: IALABEL(:), IBLABEL(:)
   
   INTEGER  H, HR, I, IC, IC2, ICCY, IER, IFIN, IFS, IGO
   INTEGER  IR, IR2R2, IRBAD, IRBADC
   INTEGER  IRBADN, IRES, IRRRR, IS, IT, ITPOS, ITVEC, IUN, J, K, L
   INTEGER  L1, MAXCOL, MAXNI, N, NAUTO, NCON, ND, NEC, NEM, NERR, NERR2, NF
   INTEGER  NFNF, NOCOV, NOMU, NPAR, NPR, NRP1, NS, NTOT, NV, NVNV, P
   INTEGER  PP, R, R2, R2R2, RP, RR, RR2, RRA, RRR2, RRU
   INTEGER  RRUN, RUN, RUN2, S, SRR, SRUN, TIMEIND
   
   INTEGER ,ALLOCATABLE:: IXIND(:)
   
   INTEGER :: INITIALIZED = 0
   
   REAL (KIND = 8), POINTER :: IALGR(:), IALIN(:), IAUTOGR(:), IAUTOIN(:), &
     IWOMX(:), IWORKSR(:)
   
   REAL*8  WORK1(1),ERREM(1),ERRGR(1)
   
   REAL (KIND=8)  AVECOR, AVECORP, BIGCOR, CONV, ERROR, ERROR0, PLIK, &
      REALN, RLOGL, RORIG, RNINV, UNO, UNON 
      
   REAL (KIND=8), ALLOCATABLE ::  IALPHA(:), IALPHA0(:), IAUERIN(:),  &
      IAUTLIN(:), IAUTO(:), IAUTO0(:), IAUV1IN(:), ICON(:), &
      ICORFIX(:), ICORVAR(:), IERRINF(:), IERV1IN(:), IERVAIN(:), &
      IGRAFIX(:), &
      IGRAVAR(:), IINFFIX(:), IINFVAR(:), IMU(:), IMU0(:), IMU1IN(:), &
      IMUEMGR(:), IMUIN(:), ITRANG(:), ITRANH(:), ITRANJ(:), ITRANN(:), &
      IVARCEM(:), IVARCIN(:), IVARCK(:), IVARCO(:), IVARCO0(:), IVARCOI(:), &
      IVARCON(:), IVARCSQ(:), IVARD(:), IVARGR(:), IVARIN(:), IVARL(:), &
      IVARLD(:), IVARLSQ(:), IVARPI(:), IW3(:), IW3G(:), IWKEC(:), &
      IWKER2(:), IWORKP(:), IWORKPP(:), IWORKR(:), IWORKR2(:), IWORKR3(:), &
      IWORKRR(:), TIME(:)
   
   REAL (KIND=8), ALLOCATABLE, TARGET :: IAL1EM(:), IAL2EM(:), &
      IALM1IN(:), IALMUIN(:), IAUTOEM(:), IAUTOEM2(:), &
      IAUTVIN(:) 
   
   CHARACTER, PARAMETER:: PROGNAME*6 = 'MIXREG'
     
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
USE RRMSET
IMPLICIT NONE

   CHARACTER*40 , INTENT(IN OUT):: FILENAMES(4)
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! LOCAL VARS ONLY - SHARED VARS ARE DECLARED
   ! BELOW THE 'MODULE' STATEMENT AT THE TOP OF THE 
   ! FILE
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   CHARACTER*20   AUTOLAB
   CHARACTER*16    XLABEL
   CHARACTER*16    YLABEL
   CHARACTER*80   FILEDAT, FILEOUT, FILERES, LINE7
   
   INTEGER  H2, IC3, IC4, ICX, ICXW, IDIND, II
   INTEGER  IR2, ISTART, K2, K3, K4, KIN1, KIN2, KIN3
   INTEGER  KIND, L2, MAXJXJ, MAXJXJ2, MAXK, MAXKK
   INTEGER  MAXXJ, MEANYX, MISS, NALL, NCATYX, NFPAR
   INTEGER  NXLAB, SS, XCIND, YIND
   
   INTEGER ,ALLOCATABLE:: IWIND(:)
   
   LOGICAL  CATJ
   
   REAL (KIND=8)  DET, YMISS, YVAL, YVAL2
   
   REAL (KIND=8), ALLOCATABLE ::  CODEX(:), IXMISS(:), ICATFQX(:), &
      IWMISS(:), IWNFNF(:), IWORKNF(:), IWRKNF(:)
   
   ! Start by checking if we are ready to run yet.  You 
   ! can't call init twice in a row, and must call 
   ! cleanup before you can start another session.
   
   IF(INITIALIZED .NE. 0) THEN
      CALL POST_ERROR("Init() routine called out of sequence in " // PROGNAME)
      RETURN
   ENDIF
   
   ! ********************************************************
   ! OPEN AND READ THE DEFINITION (.DEF) FILE
   ! ********************************************************
   ! READ IN TWO LINES (60 CHARS EACH) FOR THE TITLE 
   ! AND THEN PARAMETERS AND STARTING VALUES FROM RRM.DEF

   OPEN(1,ACTION='READ', FILE=FILENAMES(1))
    READ(1,9) HEAD
 9  FORMAT(18A4)
   READ(1,"(A80)")FILEDAT
   READ(1,"(A80)")FILEOUT
   READ(1,"(80X)")         ! Skip over the DEF file name
   READ(1,"(A80)")FILERES  
   IF(LEN_TRIM(FILENAMES(2)) .EQ. 0) FILENAMES(2) = FILEDAT 
   IF(LEN_TRIM(FILENAMES(3)) .EQ. 0) FILENAMES(3) = FILEOUT 
   
   OPEN(2, FILE= FILENAMES(3))
   OPEN(3, FILE= TRIM(ADJUSTL(FILENAMES(4))) // '.EST')
   OPEN(4, FILE= TRIM(ADJUSTL(FILENAMES(4))) // '.VAR')
      
   ! Old style DEF files have 12 params on line 7, 
   ! new style have 15.  Count the number of params
   ! on the line and read them, so we can accept either style.
   READ(1,"(A80)") Line7
   BACKSPACE(1)
   IF (NumericInputCount(Line7) > 12) THEN
      READ(1,*) NPR,NEM,MAXCOL,R,P,CONV, &
      MISS,ISTART,MEANYX,IFIN,IRES,NAUTO, &
      NOMU,NOCOV,NCON
   ELSE
      READ(1,*) NPR,NEM,MAXCOL,R,P,CONV, &
      MISS,ISTART,MEANYX,IFIN,IRES,NAUTO
   END IF
   
   ! for the time being use the unparameterized solution if
   ! only a fixed-effects model is being fit
   IF (R .LE. 0) IFIN = 1

   !      CONV = CONVERGENCE CRITERION                            
   !      IFIN = 0 ESTIMATE PHI TAU & PI IN FISHER                
   !           = 1 ESTIMATE AUTO ERROR & VARCO IN FISHER          
   !      IRES = 0 DON'T PRINT OUT INDIVIDUAL BAYES PARAMETERS          
   !           = 1 PRINT INDIVIDUAL BAYES PARAMETERS TO FILERES         
   !    ISTART = 0 AUTOMATIC STARTING VALUES
   !           = 1 USER DEFINED STARTING VALUES
   !    MEANYX = 0 DON'T GENERATE TABLE
   !           = 1 GENERATE TABLE OF MEAN OF THE DEPENDENT VAR BROKEN DOWN
   !               BY VALUES OF ANOTHER VARIABLE
   !      MISS = 0 NO MISSING VALUES ARE PRESENT
   !           = 1 MISSING VALUES ARE PRESENT
   !       NEM = NUMBER OF EM CYCLES                              
   !       NFS = NUMBER OF FISHER SCORING CYCLES TO FOLLOW EM     
   !       NPR = NUMBER OF SUBJECTS TO PRINT DATA OF              
   !                                                              
   !     NAUTO = 0  FIX AUTOCORR AT ZERO & ESTIMATE ALL ELSE
   !           = 1  FIX AUTOCORR NONZERO & ESTIMATE ALL ELSE
   !           = 2  ESTIMATE EVERYTHING
   !      NOMU = 0 IF MEAN OF RANDOM EFFECTS IS TO BE ESTIMATED
   !           = 1 IF IT IS SET TO 0
   !     NOCOV = 0 IF RANDOM EFFECT COVARIANCE TERMS ARE TO BE ESTIMATED
   !           = 1 IF THEY ARE SET TO 0
   !      NCON =  number of transforms of estimated FIXED-EFFECTS parameters 
   !              (linear re-expressions)                          
   !         P = NUMBER OF FIXED EFFECTS
   !         R = NUMBER OF RANDOM EFFECTS
   !


   ! ALLOCATE R VECTORS

   ALLOCATE (IWORKR(R))
   ALLOCATE (IXMISS(R))
   ALLOCATE (IXIND(R))

   ! ALLOCATE rr = (r * (r+1)) / 2  (unless nocov=1)

   IF (NOCOV .EQ. 0) THEN
      RR = (R * (R+1)) / 2
      RRA = RR
   ELSEIF (NOCOV .EQ. 1) THEN
      RR = R
      RRA = (R * (R+1)) / 2
   ENDIF
   ALLOCATE(IVARCO(RR))

   ! random-effects can 
   ! be reassigned as fixed-effects
   ! if numerical problems develop 

   ! set RORIG as the original value of R if nomu=1 so that if 
   ! computational troubles develop and R is decreased
   ! the W matrix is still defined correctly
   
   IF (NOMU .EQ. 1) RORIG = R
   IF (NOMU .EQ. 0) ALLOCATE(IMU(R))
   ALLOCATE(IALPHA(P))

   ALLOCATE(IBLABEL(R))
   ALLOCATE(IALABEL(P))

   ! ALLOCATE  p

   ALLOCATE(IWORKP(P))
   ALLOCATE(IWMISS(P))
   ALLOCATE(IWIND(P))

   ! ALLOCATE  NFNF = (r+p * (r+p+1)) / 2   (unless nomu=1)

   NF   = R*(1-NOMU) + P
   NFNF = (NF* (NF+1)) / 2
   ALLOCATE(IINFFIX(NFNF))

   WRITE(6,*) 'Perform ',NEM,' EM iterations'
   WRITE(6,*) R,' random terms'
   WRITE(6,*) P,' fixed  terms'
   IF(IRES.EQ.1) WRITE(6,*)'Individual parameters will be printed'


   ! MAXCOL =  dimension of the data matrix to read in     
   ! IDIND  =  column containing the level 2 ID
   ! YIND   =  column containing the dependent var
   ! R      =  dimension of the X matrix to be used in the analysis
   !           which is the number of random effects
   ! IXIND()   indicates which columns are to be used for X 
   ! P      =  dimension of the W matrix to be used in the analysis
   !           which is the number of fixed effects
   ! IWIND()   indicates which columns are to be used for W
   !   
   ! ISTART = 0 FOR AUTOMATIC STARTING VALUES
   !          1 FOR STARTING VALUES READ IN
   !
   ! MEANYX = 0 NO MEANS FOR ANY (X or W) BY Y                        
   !          1 GET MEANS FOR A SPECIFIED (X or W) BY Y   
   !          (the x or w variable is specified below 
   !           as being in column XCIND)

   READ(1,*) IDIND,YIND
   IF (R .GE. 1) READ(1,*) (IXIND(H), H=1,R)
   IF (P .GE. 1) READ(1,*) (IWIND(L), L=1,P)

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
            ENDIF
         ELSE
            WRITE(6,*)'XCIND does not match XIND or WIND'
            CATJ = .FALSE.
            MEANYX = 0
         ENDIF
      END DO

   ENDIF

   IF (MISS .EQ. 1 .or. miss .eq. 2) THEN
      READ(1,*) YMISS
      IF (R .GE. 1) READ(1,*) (IXMISS(H), H=1,R)
      IF (P .GE. 1) READ(1,*) (IWMISS(L), L=1,P)
   ENDIF

   READ(1,*) YLabel
   IF (YLABEL .EQ. '        ') READ(1,*) YLABEL

   IF (R .NE. 0) THEN
      READ(1,*) (IBLABEL(H), H=1,R)
      IF (IBLABEL(1) .EQ. '                ') READ(1,*) (IBLABEL(H), H=1,R)
      IF (ISTART .EQ. 1 .AND. NOMU .EQ. 0) READ(1,*) (IMU(H), H=1,R)
   ENDIF
   
   IF (p .NE. 0) THEN
      ! reads 8-character labels for covariates from 
      ! a maximum of 2 lines
      READ(1,*) (IALABEL(L), L=1,P)
      IF (IALABEL(1) .EQ. '                ') READ(1,*) (IALABEL(L), L=1,P)
      IF (ISTART .EQ. 1) READ(1,*) (IALPHA(L), L=1,P)
      8 FORMAT(10A8)

   ENDIF
   
   IF (MEANYX .EQ. 1) THEN
      IF (ICXW .EQ. 1) XLabel = IBLABEL(ICX)
      IF (ICXW .EQ. 2) XLabel = IALABEL(ICX)
   ENDIF

   IF (R .NE. 0) THEN
      IF (ISTART .EQ. 1) READ(1,*) (IVARCO(HR), HR=1,RR)
   ENDIF
   
   IF (ISTART .EQ. 1) READ(1,*) ERROR
   
   IF (NAUTO  .GT. 0) THEN
      READ(1,*) TIMEIND,MAXNI
      ALLOCATE(TIME(MAXNI))
      BACKSPACE(1)
      READ(1,*) TIMEIND,MAXNI,(TIME(K),K=1,MAXNI),NS,S

      ! ALLOCATE  IAUTO to size S
      ALLOCATE(IAUTO(S))

      IF (ISTART .EQ. 1 .OR. NAUTO .EQ. 1) READ(1,*) (IAUTO(IS), IS=1,S)

      ! NS    = 0      STATIONARY AR1 ERROR STRUCTURE
      !       = 1  NON-STATIONARY AR1 ERROR STRUCTURE
      !       = 2      STATIONARY MA1 ERROR STRUCTURE
      !       = 3      STATIONARY ARMA(1,1) ERROR STRUCTURE
      !       = 4      SYMMETRIC TOEPLITZ (s) ERROR STRUCTURE
      !       = 5      TOEPLITZ SMOOTHED VIA SPECTRAL ANALYSIS

      IF (NS .EQ. 0) THEN
         AUTOLAB = 'AR(1)'
      ELSEIF (NS .EQ. 1) THEN
         AUTOLAB = 'Non-Stationary AR(1)'
      ELSEIF (NS .EQ. 2) THEN
         AUTOLAB = 'MA(1)' 
      ELSEIF (NS .EQ. 3) THEN
         AUTOLAB = 'ARMA(1,1)'
      ELSEIF (NS .EQ. 4 .OR. NS .EQ. 5) THEN
        write(AUTOLAB, '(a10, i1, a1)') 'Toeplitz(', S, ')' 
      ENDIF

      ! itvec indicates which matrix contains the time vector
      !    0 = X   1 = W
      ! itpos indicates which column of this matrix is the time vector
 
      CATJ = .TRUE.
      IC = 0
      DO WHILE (CATJ)
         IC = IC + 1
         IF (IC .LE. R) THEN
            IF (TIMEIND .EQ. IXIND(IC)) THEN
               ITVEC = 0
               ITPOS = IC
               CATJ = .FALSE.
            ENDIF
         ELSEIF (IC .GT. R .AND. IC .LE. R+P) THEN
            IC2 = IC - R
            IF (TIMEIND .EQ. IWIND(IC2)) THEN
               ITVEC = 1
               ITPOS = IC2
               CATJ = .FALSE.
            ENDIF
         ELSE
            WRITE(6,*)'TIMEIND does not match XIND or WIND'
            CATJ = .FALSE.
         ENDIF
      END DO

   ELSE
      S  = 0
      NS = 0
   ENDIF

   ! dimension = ((nerr + rr) * (nerr + rr + 1)) / 2

   NERR = S + 1 
   NV = RR + NERR
   NVNV  = (NV* (NV+1)) / 2
   ALLOCATE(IINFVAR(NVNV))

   NPAR  = NF + NV
   NFPAR = NF*NCON
   
   ALLOCATE(ICON(NFPAR))

   IF (NCON .NE. 0) THEN
      READ(1,*) (ICON(L), L=1,NFPAR)
   ENDIF

   CLOSE(1)  ! FINISHED READING DEF FILE
   
   IF (IRES .EQ. 1) THEN
      IF(TRIM(FILENAMES(4)) .EQ. PROGNAME ) THEN
         FILENAMES(4) = FILERES
      ELSE
         FILENAMES(4) = PROGNAME // '.RES'
      ENDIF
      OPEN(5, FILE=FILENAMES(4))
   ENDIF

   ! COMPUTE A FEW INDICES
   ! IUN = UNIT NUMBER FOR OUTPUT (6 = SCREEN)

   NRP1 = 1 + R + P
   UNO  =  1.0D0
   UNON = -1.0D0
   IUN = 6

   ! ********************************************************
   ! GET THE DATA 
   ! ********************************************************
   ! N    (total number of level 2 units)
   ! NTOT (total number of observations)
   
   ! Note that ALLDAT and IDNI are allocated in READAT
   CALL READAT(FILEDAT,N,NTOT,MAXK,MAXCOL,R,P,ALLDAT,IDNI,IDIND, &
               YIND,IXIND,IWIND,MISS,YMISS,IXMISS,IWMISS)

   ! NALL is the number of elements for the data read in from filedat 
   ! that ALLDAT contains 

   NALL = NRP1*NTOT  

   REALN = DBLE(N)
   RNINV = 1.0D0 / REALN

   ! XLAB dimension is largest of maxKK (maxk * maxk-1 /2) 
   ! or npar (= nf + nv)
   ! NOTE - CAN'T CALL PRNT BEFORE THIS POINT!!

   MAXKK = (MAXK * (MAXK-1)) / 2
   IF (MAXKK .GE. NPAR) THEN
      NXLAB = MAXKK
   ELSE
      NXLAB = NPAR
   ENDIF
   
   IF (NXLAB .LT. NCON) NXLAB = NCON
   ALLOCATE(IXLAB(NXLAB))
   DO K = 1,NXLAB
      IXLAB(K) = '    '
   END DO

   PLIK = DBLE(NTOT) * DLOG(2.D0 * 3.1415926D0)

   ! calculate means, standard deviations, and proportions 

   MAXJXJ = 3*MAXXJ
   ALLOCATE(ICATFQX(MAXJXJ))
   
   IF (MEANYX .EQ. 1) THEN
      DO K3 = 1,MAXXJ
         DO K  = 1,3
            KIND = ((K-1)*MAXXJ) + K3
            ICATFQX(KIND) = 0.0D0
         END DO
      END DO
   ENDIF

   II = 1
   
   DO K  = 1,NALL
      K2 = K - (II - 1)
      IR     = MOD(K2,NRP1)
      IF (IR .EQ. 0) IR = NRP1
      IC     = NALL + IR
      IC2    = NALL + NRP1 + IR
      IC3    = NALL + NRP1 + NRP1 + IR
       
      !     get the minimums
      IF (K .LE. NRP1) THEN
         ALLDAT(IC) = ALLDAT(K)
      ELSE
         IF (ALLDAT(K) .LT. ALLDAT(IC)) ALLDAT(IC) = ALLDAT(K)
      ENDIF
       
      !     get the maximums
      IF (K2 .LE. NRP1) THEN
         ALLDAT(IC2) = ALLDAT(K)
      ELSE
         IF (ALLDAT(K) .GT. ALLDAT(IC2)) ALLDAT(IC2) = ALLDAT(K)
      ENDIF
       
      !     get the sums
      ALLDAT(IC3) = ALLDAT(IC3) + ALLDAT(K)
       
      !     for one of the Xs or Ws get the sum & n  by Y
      IF (MEANYX .EQ. 1) THEN
         IF (IR .EQ. 1) THEN
            YVAL      = ALLDAT(K)
         ELSEIF (IR .EQ. NCATYX) THEN
            CATJ = .true.
            K4    =  0
            DO WHILE (CATJ)
               K4    = K4 + 1 
               IF (FP_EQUAL(ALLDAT(K), CODEX(K4))) THEN
                  KIND = K4
                  ICATFQX(KIND) = ICATFQX(KIND) + YVAL
                  KIND = (2*MAXXJ) + K4
                  ICATFQX(KIND) = ICATFQX(KIND) + 1
                  CATJ = .FALSE.
               ENDIF
            END DO
         ENDIF
      ENDIF

   END DO
     
   !    calculate the means
    
   ICCY = 0
   DO K  = 1,NRP1
      IC3    = NALL + NRP1 + NRP1 + K
      ALLDAT(IC3) = ALLDAT(IC3) / DBLE(NTOT)

      !  see if the mean of the 1 random effect is 1
      !  for intraclass correlation calculation later

      IF (K .EQ. 2 .AND. R .EQ. 1 .AND. ALLDAT(IC3) .GT.  &
          0.999D0 .AND. ALLDAT(IC3) .LT. 1.001D0)  &
          ICCY = ICCY + 1
   END DO

   IF (MEANYX .EQ. 1) THEN
      DO K4 = 1,MAXXJ
         KIN3 = (2*MAXXJ) + K4
         KIN1 = K4
         ICATFQX(KIN1) = ICATFQX(KIN1) / ICATFQX(KIN3)
      END DO
   ENDIF  

   ! use IWORKNF() and IWRKNF() as work vectors (nf = p + r)
   !     IWNFNF() as a work vector (nf * (nf+1) / 2)

   IF (ISTART .EQ. 0) THEN
      ALLOCATE(IWORKNF(NF))
      ALLOCATE(IWRKNF(NF))
      ALLOCATE(IWNFNF(NFNF))
      IC = 0       
      DO L = 1,NF
         IWRKNF(L) = 0.0D0
         IWORKNF(L)  = 0.0D0
         DO L1 = 1,L
            IC = IC + 1
            IWNFNF(IC) = 0.0D0
         END DO
      END DO
   ENDIF
       
   ! get the sums of squared deviations about the means
   ! and work vectors for the regression coefficients

   YVAL2 = 0.0D0
   II = 1
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
      ENDIF

      ! IF (IR .GT. 1 .AND. ISTART .EQ. 0) THEN
      !    IR2 = IR-1
      !    IWORKNF(IR2) = ALLDAT(K) 

      IF ((ISTART .EQ. 0 .AND. NOMU .EQ. 0 .AND. IR .GT. 1) .OR. &
          (ISTART .EQ. 0 .AND. NOMU .EQ. 1 .AND. IR .GT. R+1)) THEN
         IR2 = (IR-1) - (1 + R*NOMU)
         IWORKNF(IR2+1) = ALLDAT(K) 
         
         IF (IR .EQ. NRP1) THEN
            CALL GRMCV(IWNFNF,IWNFNF,IWORKNF,UNO,NF)
            DO L = 1,NF
               IWRKNF(L) = IWRKNF(L) +(YVAL * IWORKNF(L))
            END DO
         ENDIF
      ENDIF

      ! for one of the Xs or Ws get the sum of square deviations by Y

      IF (MEANYX .EQ. 1 .AND. IR .EQ. NCATYX) THEN
         CATJ = .true.
         K4    =  0
         
         DO WHILE (CATJ)
            K4    = K4 + 1 
            
            IF (FP_EQUAL(ALLDAT(K), CODEX(K4))) THEN
               KIN1 = K4
               KIN2 = MAXXJ + K4
               ICATFQX(KIN2)=ICATFQX(KIN2)+ &
                                  ((YVAL-ICATFQX(KIN1))**2)
               CATJ = .FALSE.
            ENDIF
         END DO
      ENDIF

   END DO

   ! calculate the standard deviations

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
         ICATFQX(KIN2) = DSQRT(ICATFQX(KIN2)/ (ICATFQX(KIN3) - 1.0D0))
      END DO
   ENDIF

   ! write out descriptive statistics

   WRITE(2,554) 
   554 FORMAT(1x,'MIXREG - The program for mixed-effects linear regression analysis',/)
   
   IF (NAUTO .EQ. 2) THEN
      WRITE(2,5554)AUTOLAB
      5554 FORMAT(1X,'Autocorrelated error structure: ',A20,/)
   ELSEIF (NAUTO .EQ. 1) THEN
      WRITE(2,5556)AUTOLAB,(IAUTO(IS),IS=1,S)
      5556 FORMAT(1X,'Autocorrelated error structure: ',A20,/,1X,'==> specified value(s) = ',10F6.3)
      WRITE(2,*)
   ENDIF
    
   WRITE(2,9) HEAD

   WRITE(2,255)
   255 FORMAT(//,1x,'Numbers of observations',/,1x,'-----------------------',/)
   WRITE(2,56)N,NTOT
   56 FORMAT(1x,'Level 2 observations = ',I6,/,1x,'Level 1 observations = ',i6)
   
   IF (R .GT. 0) THEN
      WRITE(2,5610)
      WRITE(2,561)(IDNI(I),I=2,N*2,2)
      5610 FORMAT(//,1x,'The number of level 1 observations per level 2 unit  are:',/)
      561 FORMAT(1x,19I4)
   ENDIF     

   WRITE(2,257)
   257 FORMAT(//,1x,'Descriptive statistics for all variables',/, &
                 1x,'----------------------------------------',/)
   WRITE(2,357)
   357 FORMAT(1X,'Variable',13X,'     Minimum',5x,'     Maximum',5x, &
             '        Mean',5x,' Stand. Dev.',/)
   IC  = NALL + 1
   IC2 = NALL + NRP1 + 1
   IC3 = NALL + NRP1 + NRP1 + 1
   IC4 = NALL + NRP1 + NRP1 + NRP1 + 1
   
   WRITE(2,3377)
   3377 FORMAT(1X,'dependent variable',/,1X,'------------------')
   WRITE(2,377) YLABEL,ALLDAT(IC),ALLDAT(IC2),ALLDAT(IC3),ALLDAT(IC4)
   WRITE(2,*) '    '
   
   IF (R .GE. 1) THEN
      WRITE(2,3378)
      3378 FORMAT(1X,'random-effect variable(s)',/,1X,'-------------------------')
      DO H = 1,R
         IC  = IC  + 1
         IC2 = IC2 + 1
         IC3 = IC3 + 1
         IC4 = IC4 + 1
         WRITE(2,377) IBLABEL(H),ALLDAT(IC), &
            ALLDAT(IC2),ALLDAT(IC3),ALLDAT(IC4)
      END DO
      WRITE(2,*) '    '
   ENDIF
   
   IF (P .GE. 1) THEN
      WRITE(2,3379)
      3379 FORMAT(1X,'fixed regressor(s)',/,1X,'------------------')
      DO L = 1,P
         IC  = IC  + 1
         IC2 = IC2 + 1
         IC3 = IC3 + 1
         IC4 = IC4 + 1
         WRITE(2,377)IALABEL(L),ALLDAT(IC), &
            ALLDAT(IC2),ALLDAT(IC3),ALLDAT(IC4)
      END DO
   ENDIF

   377 FORMAT(1x,A16,4(5X,F12.5))

   IF (MEANYX .EQ. 1) THEN
      WRITE(2,358)YLabel,XLabel
      358 FORMAT(//,1x,'Descriptives for response variable ',A8, &
                 ' by the variable ',A8,/,1x, &
                 '--------------------------------------------------------------------',/)
      WRITE(2,372)YLabel,XLabel
      372 FORMAT(1X,8X,A8,/,1X,8X,'--------',/,1X,A8,13X,'Mean',6x, &
                 'Stand. Dev.',10x,'N',/,1x, &
                 '-----------------------------------------------------')

      MAXJXJ2 = 2*MAXXJ
      DO J = 1,MAXXJ
         WRITE(2,373)CODEX(J), &
         (ICATFQX(H),H=J,MAXJXJ2,MAXXJ),INT(ICATFQX(J+MAXJXJ2))
         373 FORMAT(/,1X,F8.3,2(5x,F12.5),5x,i6)
      END DO
   ENDIF

   ! done writing out descriptives, get starting values

   IF (ISTART .NE. 1) THEN

      ! calculate the starting values for the 
      ! regression coefficients

      CALL INVS(IWNFNF,NF,DET,IWORKNF)
      CALL MPYM(IWNFNF,IWRKNF,IWORKNF,NF,NF,1,0,1)

      IC3 = NALL + NRP1 + NRP1 + 1
      IC4 = NALL + NRP1 + NRP1 + NRP1 + 1

      IF (R .GT. 0 .AND. NOMU .EQ. 0) THEN
         DO H = 1,R
            IMU(H) = IWORKNF(H) 
         END DO
      ENDIF
      
      IF (P .GT. 0) THEN
         DO L = 1,P
            L2 = L+R*(1-NOMU)
            IALPHA(L) =  IWORKNF(L2) 
         END DO
      ENDIF
 
      CALL INVS(IWNFNF,NF,DET,IWRKNF)
      CALL GRMMT(IWORKNF,IWNFNF,WORK1,NF,1,1,IWRKNF)      
      ERROR = (YVAL2 - WORK1(1)) / (DBLE(NTOT - NF))
      ERROR0 = ERROR
 
      IF (NAUTO .GT. 1) THEN
         DO L = 1,S
            IAUTO(L) = 0.2D0**L
         END DO
      ENDIF
 
      IF (R .EQ. 1) THEN
         IVARCO(1) = .2d0 * ERROR
      ELSEIF (R .GT. 1) THEN
         HR = 0
         IF (NOCOV .EQ. 0) THEN 
            DO H  = 1,R
               DO H2 = 1,H
                  HR = HR + 1
                  IF (H2 .EQ. H) THEN
                     IF (H .EQ. 1) THEN
                        IVARCO(HR) = 1.0d0 * ERROR
                     ELSE
                        IVARCO(HR) = 0.5d0 * ERROR
                     ENDIF
                  ELSE
                     IVARCO(HR) = 0.0d0 * ERROR 
                  ENDIF
               END DO
            END DO
         ELSEIF (NOCOV .EQ. 1) THEN
            DO H = 1,R
               HR = HR+1
               IF (H .EQ. 1) THEN
                  IVARCO(HR) = 1.0d0 * ERROR
               ELSE
                  IVARCO(HR) = 0.5d0 * ERROR
               ENDIF
            END DO
         ENDIF
      ENDIF
   ENDIF

   ! write out the starting values

   WRITE(2,256)
   256 FORMAT(//,1x,'Starting values',/,1x,'---------------',/)
   
   IF (R .GT. 0 .AND. NOMU .EQ. 0) THEN
      WRITE(2,603) (IMU(H),H=1,R)
      603 FORMAT(1x,'mean      ',7F16.4)
      ALLOCATE(IMU0(R))
      CALL RELOC(IMU,IMU0,R,1,0)
   ENDIF
   
   IF (P .GT. 0) THEN
      WRITE(2,604) (IALPHA(L),L=1,P)
      604 FORMAT(1x,'covariates',7F16.4)
      ALLOCATE(IALPHA0(P))
      CALL RELOC(IALPHA,IALPHA0,P,1,0)
   ENDIF
   
   IF (R .GT. 0) THEN
      WRITE(2,605) (IVARCO(HR),HR=1,RR)
      605 FORMAT(1x,'var. terms',7F16.4)
      ALLOCATE(IVARCO0(RR))
      CALL RELOC(IVARCO,IVARCO0,RR,1,0)
   ENDIF
   
   WRITE(2,606) ERROR
   606 FORMAT(1X,'residual  ',F16.4)
   
   IF (NAUTO .GT. 1) THEN
      WRITE(2,607) (IAUTO(HR),HR=1,S)
      607 FORMAT(1x,'auto terms',7F16.4)
      ALLOCATE(IAUTO0(S))
      CALL RELOC(IAUTO,IAUTO0,S,1,0)
   ENDIF

   ! ******************************************
   ! START ITERATIONS
   !  IT = CURRENT ITERATION
   !  ND = DISPLAY NUMBER FOR MATCAL OUTPUT
   ! ******************************************
   IFS = NEM + 1
   ND  = 0

   ! ALLOCATE  IVARPI and IVARD of size R
   
   IF (NOMU .EQ. 0) ALLOCATE(IMUEMGR(R))
   ALLOCATE(IVARPI(R))
   ALLOCATE(IVARD(R))
   
   ! ALLOCATE  RR

   ALLOCATE(IWORKRR(RRA))
   ALLOCATE(IVARCEM(RR))
   IF (NOMU .EQ. 0) ALLOCATE(IMUIN(RRA))
   IF (NOMU .EQ. 0) ALLOCATE(IMU1IN(RRA))
   ALLOCATE(IVARL(RR))
   ALLOCATE(IERV1IN(RR))
   ALLOCATE(IVARCON(RR))
   ALLOCATE(IVARCOI(RR))

   ! ALLOCATE R2 = r * r

   R2   = R * R
   ALLOCATE(IVARCSQ(R2))
   ALLOCATE(IVARLSQ(R2))
   ALLOCATE(IVARGR(R2))
   ALLOCATE(IWORKR2(R2))

   ! ALLOCATE  (rr * (rr+1)) / 2   where  rr = (r * (r+1)) / 2

   RRU  = (RR * (RR+1)) / 2
   ALLOCATE(IVARCIN(RRU))

   ! ALLOCATE    r2 * r2

   R2R2 = R2 * R2
   IR2R2 = (R2 *(R2+1)) / 2
   ALLOCATE(IVARCK(R2R2))
   ALLOCATE(IVARIN(IR2R2))
   ALLOCATE(IW3(R2R2))
   ALLOCATE(ITRANN(R2R2))

   ! ALLOCATE    run * r2     where  r2 = r * r
   !                           and  run = (r * (r-1)) / 2  

   RUN  = (R * (R-1)) / 2
   RUN2 = R2 * RUN
   ALLOCATE(ITRANJ(RUN2))

   ! ALLOCATE    run * r

   RRUN = R * RUN
   ALLOCATE(IVARLD(RRUN))

   ! ALLOCATE    r2 * r

   RR2  = R2 * R
   ALLOCATE(ITRANH(RR2))
   ALLOCATE(IWORKR3(RR2))

   ! ALLOCATE    r2 * rr

   RRR2 = R2 * RR
   ALLOCATE(ITRANG(RRR2))

   ! ALLOCATE  rr * rr            where  rr = (r * (r+1)) / 2

   IRRRR = RR * RR
   ALLOCATE(IW3G(IRRRR))

   ! ALLOCATE    rr * nerr      nerr = number of error variance terms
   !                             max = s + 1

   IF (NERR .GT. 0) THEN
      NEC  = RR * NERR
      ALLOCATE(IERVAIN(NEC))
      ALLOCATE(IWKEC(NEC))
   
      ! ALLOCATE    rr + nerr
   
      ALLOCATE(IGRAVAR(NV))
      ALLOCATE(ICORVAR(NV))
   
      ! ALLOCATE  (nerr * (nerr+1)) / 2
   
      NERR2= (NERR * (NERR+1)) / 2
      ALLOCATE(IERRINF(NERR2))
      ALLOCATE(IWKER2(NERR2))
   ENDIF

   IF(NAUTO .GT. 1) THEN

      ! ALLOCATE  S

      ALLOCATE(IAUTOEM(S))
      IAUTOGR => IAUTOEM
      ALLOCATE(IAUERIN(S))

      ! ALLOCATE  SS = S * (S+1) / 2

      SS = S * (S+1) / 2
      ALLOCATE(IAUTOEM2(SS))
      IAUTOIN => IAUTOEM2

      ! ALLOCATE  s  *  rr [ = (r * (r+1) / 2)]

      ! SRR = S * RR        RR = (R * (R+1)) / 2
      SRR = S * (R * (R+1)) / 2
      ALLOCATE(IAUTVIN(SRR))
      ALLOCATE(IAUV1IN(SRR))
      IWORKSR => IAUTVIN

      ! ALLOCATE  s  *  run [ = (r * (r-1)) / 2]

      SRUN = S * RUN
      ALLOCATE(IAUTLIN(SRUN))

   ENDIF

   ! ALLOCATE  P

   ALLOCATE(IAL2EM(P))
   IALGR => IAL2EM

   ! ALLOCATE pp = (p * (p+1)) / 2

   PP   = (P * (P+1)) / 2
   ALLOCATE(IAL1EM(PP))
   ALLOCATE(IWORKPP(PP))
   IALIN => IAL1EM

   ! ALLOCATE rp = r * p

   RP   = R * P
   ALLOCATE(IALMUIN(RP))
   ALLOCATE(IALM1IN(RP))
   IWOMX => IALMUIN

   ! ALLOCATE  nf = r + p

   ALLOCATE(IGRAFIX(NF))
   ALLOCATE(ICORFIX(NF))

   ! ******************************************
   ! GET COMMUTATION MATRIX N AND
   ! TRANSFORMATION MATRICES FOR DIAGONAL,
   ! STRICTLY LOWER TRIANGULAR, & SYMMETRIC
   !
   ! CALL TRANNN(ITRANN(),R,R2,R2R2,WORKRR,WORKR2,IW3)
   ! gives N_(rr) in ITRANN() (r2r2)
   ! needs work space WORKRR(r) WORKR2(r2) IW3(r2r2)
   !
   ! CALL TRANDI(ITRANH(),R,RR2)
   ! gives PSI_(r) in ITRANH() (r * r*r)
   ! for diagonal transformation
   !
   ! CALL TRANLT(ITRANJ(),R,RUN2)
   ! gives tilde L_(r) in ITRANJ() (run * r*r)
   ! for strictly lower triangular trans
   !
   ! CALL TRANGG(ITRANG(),R,RRR2,R2R2,ITRANN())
   ! gives G_(r) in ITRANG() (r*r * rr)
   ! for symmetric transformation
   !
   ! r2   = r  * r
   ! rr   = (r * (r+1))/2
   ! r2r2 = r2 * r2
   ! rr2  = r2 * r
   ! rrr2 = r2 * rr
   ! run  = (r * (r-1)) / 2
   ! run2 = r2 * run
   ! ******************************************

   IF (R .GT. 0) THEN
      CALL TRANLT(ITRANJ,R,RUN2)
      CALL TRANDI(ITRANH,R,RR2)
      CALL TRANNN(ITRANN,R,R2,R2R2,IWORKRR,IWORKR2,IW3)
      CALL TRANGG(ITRANG,R,RRR2,R2R2,ITRANN)
   ENDIF

   ! Set up scalars for the big iteration loop
   BIGCOR  = 999.0D0
   AVECORP = 999.0D0
   AVECOR  = 999.0D0

   IT      = 0
   IRBAD   = 0
   IRBADC  = 0
   IRBADN  = 0
   IER     = 0
   IGO = 1
   
   ! Get rid of init-only allocated arrays
   IF (ALLOCATED(CODEX)) DEALLOCATE(CODEX)
   DEALLOCATE( IXMISS, ICATFQX, IWMISS, IWIND) 
   IF (ALLOCATED(IWNFNF)) DEALLOCATE(IWNFNF, IWORKNF, IWRKNF)
   
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
USE RRMSET
IMPLICIT NONE
   
   SAVE
   ! The SAVE statement preserves the values in local 
   ! arrays between calls to the MainLoop function
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! LOCAL VARS ONLY - SHARED VARS ARE DECLARED
   ! BELOW THE 'MODULE' STATEMENT AT THE TOP OF THE 
   ! FILE
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   CHARACTER (LEN=8), ALLOCATABLE :: TEMPLABEL(:)
   
   INTEGER  ICI2, ICID, ICII, ICOUNT, ICSUM, IDI, IH
   INTEGER  IJ, IL, IN, IP, IRR, IRRR, IS2, LP, MODID
   INTEGER  NII, NII2, NINI, NINS, NNS, NOC, NSPEC, NSPECS
   INTEGER  PNII, RC, RNII, RNIIO, TEMP_I
   
   REAL (KIND = 8), POINTER :: IAMAT(:), IXVARP(:)
   
   REAL (KIND=8)  AUTOOLD, AVECORD, CORSUM, DET, DETB, DETE
   REAL (KIND=8)  DETP, DETR, DETV, DPHI, ERLIK, XNII
   REAL (KIND=8)  ERRINV, PHI, QUERR, TAU, VALIK
   
   REAL (KIND=8), ALLOCATABLE ::  IBETA(:), IDEV(:), IEI(:), &
      IOMEGA(:), IOMEGAD(:), IOMEGAD0(:), IOMEGSQ(:), IUU(:), &
      IUVEC(:), IVARPO(:), IWAL(:), IWI(:), IWOM(:), &
      IWORKE(:), IWORKNI2(:), IWORKRNI(:), IWORKRRR(:), &
      IWRKNI2(:), IWRKNINI(:), IWRKNINS(:), IXI(:), &
      IXMU(:), IXTI(:), IYI(:), IYXB(:)
   
   REAL (KIND=8), ALLOCATABLE, TARGET :: IBE1(:), IXVARX(:)
     
   IF(INITIALIZED .EQ. 0) THEN
      CALL POST_ERROR("Mainloop() routine called before Init() " // PROGNAME)
      RETURN
   ELSEIF(INITIALIZED .EQ. 2) THEN
      CALL POST_ERROR("Mainloop() routine called after completion. " // PROGNAME)
      RETURN
   ENDIF
   
   ! **********************************************************
   !  Perform a single iteration of the main loop
   ! **********************************************************
   
   IF (IGO .EQ. 1 ) THEN
      IT      = IT + 1

      ! LIKELIHOOD INITIALIZATION

      DETE  = 0.0D0
      DETB  = 0.0D0
      ERLIK = 0.0D0
      VALIK = 0.0D0
 
      ! CHECK TO SEE IF FEWER RANDOM-EFFECTS ARE TO
      ! BE ESTIMATED AND MODIFY IALPHA ARRAY
      ! AND ALL INDICES ACCORDINGLY
 
      IF (IRBADC .GE. 1) THEN
         IT      = 1
         IRBAD   = 0
         IRBADC  = 0
         IER     = 0
         BIGCOR  = 999.0D0
         AVECORP = 999.0D0
         AVECOR  = 999.0D0
         
         IF (NOMU .EQ. 0) THEN
            ! Must insert an element into first position 
            ! will be filled in below from IALPHA0
            TEMP_I = UBOUND(IALPHA,1)
            DEALLOCATE(IALPHA)
            ALLOCATE(IALPHA(TEMP_I+1))
            
            TEMP_I = UBOUND(IALABEL,1)
            IF (ALLOCATED(TEMPLABEL)) DEALLOCATE(TEMPLABEL)
            ALLOCATE(TEMPLABEL(TEMP_I+1))
            TEMPLABEL(2:) = IALABEL   ! Array copy to temp
            DEALLOCATE(IALABEL)
            ALLOCATE(IALABEL(TEMP_I+1))
            IALABEL = TEMPLABEL       ! Array copy
            IALABEL(1) = IBLABEL(UBOUND(IBLABEL,1))
         ENDIF
 
         ! Decrement the count of random effects, and all 
         ! of the dependent counters
         R  = R-1
         IF (NOCOV .EQ. 0) THEN
            RR = (R * (R+1)) / 2
            RRA = RR
         ELSEIF (NOCOV .EQ. 1) THEN
            RR = R
            RRA = (R * (R+1)) / 2
         ENDIF
         
         R2   = R * R
         RRU  = (RR * (RR+1)) / 2
         R2R2 = R2 * R2
         IR2R2 = (R2 *(R2+1)) / 2
         RUN  = (R * (R-1)) / 2
         RUN2 = R2 * RUN
         RRUN = R * RUN
         RR2  = R2 * R
         RRR2 = R2 * RR
         IRRRR = RR * RR
         NEC  = RR * NERR
         SRR = S * RR
         SRUN = S * RUN
         NV = RR + NERR
         NVNV  = (NV* (NV+1)) / 2
 
         ! Increase the number of fixed effects if NOMU .EQ. 0
         
         IF (NOMU .EQ. 0) P  = P+1      
         PP   = (P * (P+1)) / 2
         NF   = R*(1-NOMU) + P
         NFNF = (NF* (NF+1)) / 2
 
         IF (NOMU .EQ. 0) THEN
            NRP1 = 1 + R + P
         ELSEIF (NOMU .EQ. 1) THEN
            NRP1 = 1 + R + P + 1
            ! this is to preserve NRP1 to the original R 
         ENDIF
         
         RP   = R * P
         NPAR = NF + NV
 
         IF (R .GT. 0) THEN
            CALL TRANLT(ITRANJ,R,RUN2)
            CALL TRANDI(ITRANH,R,RR2)
            CALL TRANNN(ITRANN,R,R2,R2R2,IWORKRR,IWORKR2,IW3)
            CALL TRANGG(ITRANG,R,RRR2,R2R2,ITRANN)
         ELSE
            IFIN = 1
         ENDIF
 
         ! RE ALLOCATE  P

         IF (ALLOCATED(IAL2EM)) DEALLOCATE(IAL2EM)
         ALLOCATE(IAL2EM(P))
         
         IF (ALLOCATED(IWORKP)) DEALLOCATE(IWORKP)
         ALLOCATE(IWORKP(P))
         IALGR => IAL2EM

         ! RE ALLOCATE pp = (p * (p+1)) / 2

         PP   = (P * (P+1)) / 2
         IF (ALLOCATED(IAL1EM)) DEALLOCATE(IAL1EM)
         ALLOCATE(IAL1EM(PP))
         IF (ALLOCATED(IWORKPP)) DEALLOCATE(IWORKPP)
         ALLOCATE(IWORKPP(PP))
         IALIN => IAL1EM

         ! RE ALLOCATE rp = r * p

         RP   = R * P
         IF (R .GT. 0 .AND. NOMU .EQ. 0) THEN
            IF (ALLOCATED(IALMUIN)) DEALLOCATE(IALMUIN)
            ALLOCATE(IALMUIN(RP))
            IF (ALLOCATED(IALM1IN)) DEALLOCATE(IALM1IN)
            ALLOCATE(IALM1IN(RP))
            IWOMX => IALMUIN
          ENDIF
         
         ! get ORIGINAL starting values
         
         IF (R .GT. 0 .AND. NOMU .EQ. 0) THEN
            CALL RELOC(IMU0,IMU,R,1,0)
            CALL RELOC(IVARCO0,IVARCO,RR,1,0)
         ELSEIF (R .GT. 0 .AND. NOMU .EQ. 1) THEN
            CALL RELOC(IVARCO0,IVARCO,RR,1,0)
         ENDIF
         
         IF (P .GT. 0 .AND. NOMU .EQ. 0) THEN
            ! note that R now equals the old R minus 1
            ! so that imu0+R equals the old imu0+R-1 
            ! which should be the last element of mu
            DO IL = 1,IRBADN
               IALPHA(IL) = IMU0(R+IL)
            END DO
            ! copy back the rest of the array from IALPHA0
            CALL RELOC(IALPHA0,IALPHA(IRBADN+1),P-IRBADN,1,0)
         ELSEIF (P .GT. 0 .AND. NOMU .EQ. 1) THEN
            CALL RELOC(IALPHA0,IALPHA,P,1,0)
         ENDIF
         
         ERROR = ERROR0
         IF (NAUTO .GT. 1) CALL RELOC(IAUTO0,IAUTO,S,1,0)

         ! reassign TIMEIND if appropriate

         IF (ITVEC .EQ. 0 .AND. TIMEIND .EQ. IXIND(R)) THEN
            ITVEC = 1
            ITPOS = 1
         ENDIF
 
      ENDIF !  IRBADC .GE. 1
 
      IF (IT .LE. NEM) THEN

         ! EM INITIALIZATION
         IF (NOMU .EQ. 0) THEN
            DO H = 1,R
               IMUEMGR(H)  = 0.0D0
            END DO
         ENDIF
         DO L = 1,PP
            IAL1EM(L)   = 0.0D0
         END DO
         DO L = 1,P
            IAL2EM(L)   = 0.0D0
         END DO
         DO H = 1,RR
            IVARCEM(H)  = 0.0D0
         END DO
         ERREM(1) = 0.0D0
         IC = 0
         DO IS = 1,S
            IAUTOEM(IS)= 0.0D0
            DO IS2= 1,IS
               IC = IC + 1
               IAUTOEM2(IC) = 0.0D0
            END DO
         END DO
      ELSE
         
         ! FISHER SCORING INITIALIZATION
         IF (NOMU .EQ. 0) THEN
            DO H = 1,R
               IMUEMGR(H) = 0.0D0
            END DO
            DO L = 1,RP
               IALM1IN(L) = 0.0D0
            END DO
         ENDIF
         
         DO L = 1,P
            IALGR(L)   = 0.0D0
         END DO
         
         DO L = 1,PP
            IALIN(L)   = 0.0D0
         END DO
         
         IF (NOMU .EQ. 0) THEN 
            DO H = 1,RRA
               IMU1IN(H)  = 0.0D0
            END DO
         ENDIF
         
         DO H = 1,UBOUND(IVARGR,1)
            IVARGR(H)  = 0.0D0
         END DO
         DO H = 1,UBOUND(IVARCIN,1)
            IVARCIN(H) = 0.0D0
         END DO
         
         IC = 0
         IC2= 0
         
         DO IS= 1,S
            IAUTOGR(IS)= 0.0D0
            IAUERIN(IS)= 0.0D0
            DO IS2= 1,IS
               IC = IC + 1
               IAUTOIN(IC) = 0.0D0
            END DO
            DO H = 1,RR
               IC2= IC2+ 1
               IAUV1IN(IC2)= 0.0D0
            END DO
         END DO
         
         ERRGR(1) = 0.0D0
         IERRINF(1) = 0.0D0
         DO H = 1,RR
            IERV1IN(H) = 0.0D0
         END DO
         
         DO H = 1, UBOUND(IVARIN,1)
            IVARIN(H) = 0.0D0
         END DO
      ENDIF
      

      ! INVERT VARIANCES FOR CALCULATIONS 

      IF (R .GT. 0 .AND. NOCOV .EQ. 0) THEN
         CALL RELOC(IVARCO,IVARCOI,R,R,1)
         CALL INVS(IVARCOI,R,DETV,IWORKRR)
         CALL CHAMS(IVARCOI,IVARCSQ,R,1,0)
      ELSEIF (R .GT. 0 .AND. NOCOV .EQ. 1) THEN
         CALL RELOC(IVARCO,IVARCOI,R,R,2)
         CALL INVD(IVARCOI,R,DETV)
         CALL CHAMS(IVARCOI,IVARCSQ,R,2,1)
      ENDIF
 

      ! get the LDL' factorization for the first Fisher scoring iteration
      ! if r=0  then IFIN=1 so don't do this 

      IF (R .GT. 0 .AND. IT .EQ. IFS) THEN
         IF (NOCOV .EQ. 0) THEN
            CALL GFACT(IVARCO,IVARL,IVARD,R,RR)
            DO H=1,R
               IVARPI(H) = DLOG(IVARD(H))
            END DO
         ELSEIF (NOCOV .EQ. 1) THEN
            DO H=1,R
               IVARPI(H) = DLOG(IVARCO(H))
            END DO
         ENDIF
         TAU = DLOG(ERROR)
         IF (NAUTO .GT. 1 .AND. NS .LE. 2) THEN
            PHI  = IAUTO(1) / (DSQRT(1.0D0 - IAUTO(1)*IAUTO(1)))
            DPHI = 1.0D0 / ((1.0D0 + PHI*PHI) * (DSQRT(1.0D0 + PHI*PHI)))
         ENDIF
      ENDIF
 
      ERRINV = 1.0D0 / ERROR

      ! ***********************************************************
      ! GO OVER ALL SUBJECTS
      !    IUN = UNIT NUMBER FOR OUTPUT (0 = SCREEN)
      ! ***********************************************************

      ICOUNT = 0
 
      IF (R .GT. 0) THEN
         
         IF (ALLOCATED(IBETA)) THEN
            IF(UBOUND(IBETA,1) .NE. R) THEN
               DEALLOCATE(IBETA)
               DEALLOCATE(IDEV)
               ALLOCATE(IBETA(R))
               ALLOCATE(IDEV(R))
            ENDIF
         ELSE
            ALLOCATE(IBETA(R))
            ALLOCATE(IDEV(R))
         ENDIF
   
         ! must reallocate in case R has changed.
         IF (ALLOCATED(IVARPO)) THEN
            IF(UBOUND(IVARPO,1) .NE. RRA) THEN
               DEALLOCATE(IVARPO)
               ALLOCATE(IVARPO(RRA))
            ENDIF
         ELSE
            ALLOCATE(IVARPO(RRA))
         ENDIF
   
      ENDIF
 
      !WRITE(6,"('first 10 of idni: ',10I6)")(IDNI(J),J=1,10)
      
      DO I = 1,N   ! loop through each subject. 
 
         IC2  = 2*(I-1)+1
         IDI  = IDNI(IC2)
         IC2  = IC2+1
         NII  = IDNI(IC2)
         NII2 = NII * NII
         RNII = NII * R
         IF (NOMU .EQ. 1) RNIIO = NII * RORIG
         PNII = NII * P
         NINI = (NII * (NII + 1)) / 2
         NINS = NINI  * S
         XNII = DBLE(NII)
         
         !WRITE(6,"('Subject ',I4,' has id ',I4,' and ', I4,' obs.')")
         !+      I,IDI,NII
         
         ! ALLOCATE maxni
   
         IF (ALLOCATED(IXTI)) THEN
            IF(UBOUND(IXTI,1) .NE. NII) THEN
               DEALLOCATE(IXTI)
               DEALLOCATE(IWORKE)
               DEALLOCATE(IXMU)
               DEALLOCATE(IWAL)
               DEALLOCATE(IYXB)
               DEALLOCATE(IUVEC)
               ALLOCATE(IXTI(NII))
               ALLOCATE(IWORKE(NII))
               ALLOCATE(IXMU(NII))
               ALLOCATE(IWAL(NII))
               ALLOCATE(IYXB(NII))
               ALLOCATE(IUVEC(NII))
            ENDIF
         ELSE
            ALLOCATE(IXTI(NII))
            ALLOCATE(IWORKE(NII))
            ALLOCATE(IXMU(NII))
            ALLOCATE(IWAL(NII))
            ALLOCATE(IYXB(NII))
            ALLOCATE(IUVEC(NII))
         ENDIF
   
         ! ALLOCATE  NINI = (nii * (nii+1)) / 2
   
         IF (NAUTO .GT. 0) THEN
            IF (ALLOCATED(IOMEGA)) THEN
               IF(UBOUND(IOMEGA,1) .NE. NINI) THEN
                  DEALLOCATE(IOMEGA)
                  ALLOCATE(IOMEGA(NINI))
               ENDIF
            ELSE
               ALLOCATE(IOMEGA(NINI))
            ENDIF
         ENDIF
      
         IF (ALLOCATED(IXVARX)) THEN
            IF(UBOUND(IXVARX,1) .NE. NINI) THEN
               DEALLOCATE(IXVARX)
               DEALLOCATE(IWRKNINI)
               ALLOCATE(IXVARX(NINI))
               ALLOCATE(IWRKNINI(NINI))
            ENDIF
         ELSE
            ALLOCATE(IXVARX(NINI))
            ALLOCATE(IWRKNINI(NINI))
         ENDIF
   
         IAMAT => IXVARX
   
         ! ALLOCATE   NII2 =  nii * nii
   
         IF (ALLOCATED(IEI)) THEN
            IF(UBOUND(IEI,1) .NE. NII2) THEN
               DEALLOCATE(IEI)
               DEALLOCATE(IUU)
               ALLOCATE(IEI(NII2))
               ALLOCATE(IUU(NII2))
            ENDIF
         ELSE
            ALLOCATE(IEI(NII2))
            ALLOCATE(IUU(NII2))
         ENDIF

         IF (NAUTO .GT. 0) THEN
            IF (ALLOCATED(IOMEGSQ)) THEN
               IF(UBOUND(IOMEGSQ,1) .NE. NII2) THEN
                  DEALLOCATE(IOMEGSQ)
                  DEALLOCATE(IWORKNI2)
                  DEALLOCATE(IWRKNI2)
                  ALLOCATE(IOMEGSQ(NII2))
                  ALLOCATE(IWORKNI2(NII2))
                  ALLOCATE(IWRKNI2(NII2))
               ENDIF
            ELSE
               ALLOCATE(IOMEGSQ(NII2))
               ALLOCATE(IWORKNI2(NII2))
               ALLOCATE(IWRKNI2(NII2))
            ENDIF
         ENDIF
   
         ! ALLOCATE   NINS =  s * NINI 
   
         IF (NAUTO .EQ. 1) THEN
            IF (ALLOCATED(IOMEGAD)) THEN
               IF(UBOUND(IOMEGAD,1) .NE. NINS) THEN
                  DEALLOCATE(IOMEGAD)
                  ALLOCATE(IOMEGAD(NINS))
               ENDIF
            ELSE
               ALLOCATE(IOMEGAD(NINS))
            ENDIF
         ELSEIF (NAUTO .EQ. 2) THEN
            IF (ALLOCATED(IOMEGAD)) THEN
               IF(UBOUND(IOMEGAD,1) .NE. NINS) THEN
                  DEALLOCATE(IOMEGAD)
                  DEALLOCATE(IOMEGAD0)
                  DEALLOCATE(IWRKNINS)
                  ALLOCATE(IOMEGAD(NINS))
                  ALLOCATE(IOMEGAD0(NINS))
                  ALLOCATE(IWRKNINS(NINS))
               ENDIF
            ELSE
               ALLOCATE(IOMEGAD(NINS))
               ALLOCATE(IOMEGAD0(NINS))
               ALLOCATE(IWRKNINS(NINS))
            ENDIF
         ENDIF
   
         ! ALLOCATE   RNII = r * maxni
   
         IF (ALLOCATED(IBE1)) THEN
            IF(UBOUND(IBE1,1) .NE. RNII) THEN
               DEALLOCATE(IBE1)
               DEALLOCATE(IWORKRNI)
               ALLOCATE(IBE1(RNII))
               ALLOCATE(IWORKRNI(RNII))
            ENDIF
         ELSE
            ALLOCATE(IBE1(RNII))
            ALLOCATE(IWORKRNI(RNII))
         ENDIF
   
         IXVARP => IBE1
   
         ! ALLOCATE   PNII = p * maxni 
   
         IF (ALLOCATED(IWOM)) THEN
            IF(UBOUND(IWOM,1) .NE. PNII) THEN
               DEALLOCATE(IWOM)
               ALLOCATE(IWOM(PNII))
            ENDIF
         ELSE
            ALLOCATE(IWOM(PNII))
         ENDIF
   
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
   
   
         ! THE X(K, H) MATRIX  K = 1 .. NI(I)  H = 1 .. R
         ! THE IXTI(K) VECTOR  K = 1 ... NI(I)   LINEAR TERM - ACTUAL TIMEPTS
         ! which is read from X if R > 1
         !                    W if R = 1
   
         IF (R .GT. 0) THEN
            IF (ALLOCATED(IXI)) THEN
               IF(UBOUND(IXI,1) .NE. RNII) THEN
                  DEALLOCATE(IXI)
                  ALLOCATE(IXI(RNII))
               ENDIF
            ELSE
               ALLOCATE(IXI(RNII))
            ENDIF
            
            HR = 0
            IC = 0
            
            DO H = 1, R
               DO K = 1, NII
                  IC  = IC + 1
                  IC2 = ICOUNT + (NRP1 * (K-1) + H+1) 
                  IXI(IC) = ALLDAT(IC2)
                  HR  = HR + 1
                  IF (ITVEC .EQ. 0) THEN
                     IF ((HR .GT. (ITPOS-1)*NII)  &
                         .AND. (HR .LE. ITPOS*NII)) THEN
                         IXTI(K) = ALLDAT(IC2)
                      ENDIF
                  ENDIF
               END DO
            END DO
         ENDIF
   
         ! THE W(K, L) MATRIX  K = 1 .. NI(I)  L = 1 .. P
   
         IF (P .GT. 0) THEN
            IF (ALLOCATED(IWI)) THEN
               IF(UBOUND(IWI,1) .NE. PNII) THEN
                  DEALLOCATE(IWI)
                  ALLOCATE(IWI(PNII))
               ENDIF
            ELSE
               ALLOCATE(IWI(PNII))
            ENDIF
            
            LP = 0
            IC = 0
            DO L = 1,P
               DO K = 1,NII
                  IC         = IC + 1
                  IF (NOMU .EQ. 0) THEN
                     IC2        = ICOUNT + (NRP1 * (K-1) + L+R+1)
                  ELSEIF (NOMU .EQ. 1) THEN
                     IC2        = ICOUNT + (NRP1 * (K-1) + L+RORIG+1)
                  ENDIF
                  
                  IWI(IC) = ALLDAT(IC2)
                  LP = LP + 1
                  IF (ITVEC .EQ. 1) THEN
                     IF ((LP .GT. (ITPOS-1)*NII)  &
                        .AND. (LP .LE. ITPOS*NII)) THEN
                          IXTI(K) = ALLDAT(IC2)
                     ENDIF
                  ENDIF
               END DO
            END DO
         ENDIF
    
         IF (NOMU .EQ. 0) THEN
            ICOUNT = ICOUNT + NII + RNII + PNII 
         ELSEIF (NOMU .EQ. 1) THEN
            ICOUNT = ICOUNT + NII + RNIIO + PNII 
         ENDIF
    
         ! THE IOMEGA(IE) MATRIX IE = 1 ... NINI  ERROR MATRIX 
          
         ! NOTE: IWRKNI2() AND IWORKNI2() ARE USED HERE FOR WORK VECTORS
         !       OF SIZE NII2 = NII * NII
          
         ! for EM iterations use AR1 case for MA1 and NS AR1 
   
         NNS = NS
         IF (IT .LE. NEM .AND. (NS .EQ. 1 .OR. NS .EQ. 2)) NNS = 0
    
         ! used to set the IOMEGA matrix to an identity if NAUTO .EQ. 0
         ! but this was commented out.  TG 1/10/00
         
         IF (NAUTO .NE. 0) THEN
   
            CALL OMEGAI(NNS,MAXNI,TIME,IAUTO,S,NII,NII2,NINI, &
                        IXTI,RC,IOMEGA,IOMEGAD,IWRKNI2,IWORKNI2)
       
            IF (RC .LT. 0) THEN
               WRITE(IUN,331) IDI
               331 FORMAT(1X,'TOO MANY MISSING TIMEPOINTS FOR SUBJECT ',I8)
               EXIT ! skip over this subject by going back to loop start
      
            ELSEIF (RC .GT. 0) THEN
               WRITE(IUN,333) IDI
               333 FORMAT(1X,'TOO FEW MISSING TIMEPOINTS FOR SUBJECT ',I8)
               EXIT ! skip over this subject by going back to loop start
      
            ELSE    ! RC .EQ. 0
               IF (NAUTO .GT. 1) THEN
                  IF (NS .LE. 2) THEN
                     !
                     ! ===> leave omegad as symmetric
                     !
                     !  CALL CHAMS(IOMEGAD,IOMEGAD,NII,1,0)
                     IF (NAUTO .GT. 1 .AND. IFIN .EQ. 0 .AND. IT .GT. NEM) &
                        CALL SCM(IOMEGAD,DPHI,IOMEGAD,NII,NII,1)
                        ! CALL SCM(IOMEGAD,DPHI,IOMEGAD,NII,NII,0)
      
                     ! **************************************************************
                     ! at this point IOMEGAD() (for ns >2) is a matrix (in vec form)
                     ! indicating the lag (or diagonal) index for each element
                     ! thus this matrix needs to be passed to routines ARMADER and
                     ! TRANGA to get the real derivative matrices
                     ! **************************************************************
                  ELSEIF (NS .EQ. 3) THEN   
                     
                     !  ARMA(1,1) TERMS (ns= 3)
                     CALL RELOC(IOMEGAD,IWRKNI2,NII,NII,1)
                     CALL ARMADER(IOMEGAD,S,NINI,NINS,IWRKNI2,IAUTO)
   
                  ELSE  
                     
                     ! TOEPLITZ TERMS (ns = 4 or 5)
                     CALL RELOC(IOMEGAD,IWRKNI2,NII,NII,1)
                     CALL TRANGA(IOMEGAD,S,NINI,NINS,IWRKNI2)
                  
                  ENDIF  ! NS .LE. 2/.EQ. 3
                  
               ENDIF  ! NAUTO .GT. 1
            ENDIF  ! RC .LT./.GT./.EQ. 0
         ENDIF ! NAUTO .NE. 0
   
         ! PRINT THE DATA OF THE 1ST NPR SUBJECTS
         ! TO THE UNIT NUMBER IUN ON THE FIRST ITERATION ONLY
   
         IF (IRBADN .EQ. 0 .AND. IT .EQ. 1 .AND. I .LE. NPR) THEN
         
            WRITE(IUN,22)N
            22 FORMAT(1X,'Total Number of Level-2 Units = ',I6,/)
            WRITE(IUN,12) IDI,NII
            12 FORMAT(1X,'Data for Level-2 Unit ',I10,' which has', &
               I6,' observations nested within')
            
            CALL PRNT(IUN,IYI,NII,1,0,IXLAB,IXLAB,ND,HEAD,1,80,5,1,1, &
                      "Dependent Variable Vector")
         
            IF (NAUTO .GT. 0) THEN
               CALL PRNT(IUN,IXTI,NII,1,0,IXLAB,IXLAB,ND,HEAD,1,80,5,1,1, &
                      "Vector of Timepoints")
            ENDIF
         
            IF (R .GT. 0) THEN
               CALL PRNT(IUN,IXI,NII,R,0,IXLAB,IXLAB,ND,HEAD,1,80,5,1,1, &
                      "Random-effect Design Matrix")
            ENDIF
         
            IF (P .GT. 0) THEN
               CALL PRNT(IUN,IWI,NII,P,0,IXLAB,IXLAB,ND,HEAD,1,80,5,1,1, &
                      "Covariate Matrix")
            ENDIF
         
            IF (NAUTO .GT. 0) THEN
               CALL PRNT(IUN,IOMEGA,NII,NII,1,IXLAB,IXLAB,ND,HEAD,1,80,7,1,1, &
                      "Autocorrelation Matrix")
            ENDIF
         
            IF (NAUTO .GT. 1) THEN
               IF (NS .LE. 2) THEN
                  CALL PRNT(IUN,IOMEGAD,NII,NII,1,IXLAB,IXLAB,ND,HEAD,1,80,7,1,1, &
                      "OMEGA DERIVATIVE MATRIX w/ respect to AUTO")
               ELSE
                  CALL PRNT(IUN,IOMEGAD,NINI,S,0,IXLAB,IXLAB,ND,HEAD,1,80,7,1,1, &
                      "OMEGA DERIVATIVE MATRIX w/ respect to AUTO")
               ENDIF
            ENDIF
         ENDIF
    
         ! print out the IOMEGA() matrices for the first Fisher iteration
         ! since they were set to the AR1 case for NS AR1 and MA1 in EM
    
         IF (IRBADN .EQ. 0 .AND. IT .EQ. IFS .AND. I .LE. NPR .AND.  &
              NAUTO .GT. 1 .AND. (NS .EQ. 1 .OR. NS .EQ. 2)) THEN
            CALL PRNT(IUN,IOMEGA,NII,NII,1,IXLAB,IXLAB,ND,HEAD,1,80,7,1,1, &
                   "OMEGA MATRIX")
            CALL PRNT(IUN,IOMEGAD,NII,NII,1,IXLAB,IXLAB,ND,HEAD,1,80,7,1,1, &
                   "OMEGA DERIVATIVE MATRIX w/ respect to AUTO           ")
         ENDIF
   
   
         ! *********************************************************
         ! START - COMPUTE POSTERIOR BETA VARIANCE COVARIANCE
         ! *********************************************************
   
         IF (R .GT. 0) THEN  ! R is count of random effects
   
            IF (NAUTO .EQ. 0) THEN
               DETE = DETE + XNII*DLOG(ERROR)
               CALL GRAMT(IXI,IWORKRR,NII,R)
            ELSE
               CALL INVS(IOMEGA,NII,DET,IWORKE,IER)
               IF (DET .LE. 0.000000000000001D0) DET = 0.000000000000001D0
               DETE = DETE + (DLOG(DET) + XNII*DLOG(ERROR))
               CALL GRMMT(IXI,IOMEGA,IWORKRR,NII,R,1,IWORKE)
            ENDIF
   
            CALL SCM(IWORKRR,ERRINV,IWORKRR,R,R,1)
            IF (NOCOV .EQ. 0) THEN
               CALL ADDM(IWORKRR,IVARCOI,IVARPO,R,R,1)
            ELSEIF (NOCOV .EQ. 1) THEN
               CALL ADDM(IWORKRR,IVARCSQ,IVARPO,R,R,1)
            ENDIF
            CALL INVS(IVARPO,R,DET,IWORKRR,IER)
    
            ! SUM THE DETERMINANTS OF VARPO  (which equal 1 / DET)
            ! since DET is the determinant of the inverse of VARPO
            
            IF (DET .LE. 0.000000000000001D0) DET = 0.000000000000001D0
            DETB = DETB + DLOG(1.0D0) - DLOG(DET)
   
            ! NOW CALCULATE THE BETA VECTOR 
   
            CALL RELOC(IYI,IWORKE,NII,1,0)
            IF (NOMU.EQ.0) CALL MPYM(IXI,IMU,IXMU,NII,R,0,0,1)
            
            IF (P .EQ. 0) THEN
               IF (NOMU .EQ. 0) CALL SUBM(IWORKE,IXMU,IWORKE,NII,1,0)
            ELSE
               CALL MPYM(IWI,IALPHA,IWAL,NII,P,0,0,1)
               CALL SUBM(IWORKE,IWAL,IWORKE,NII,1,0)
               IF (NOMU .EQ. 0) CALL SUBM(IWORKE,IXMU,IWORKE,NII,1,0)
            ENDIF
   
            IF (NAUTO .EQ. 0) THEN
               CALL TRP(IXI,IWORKRNI,NII,R)
            ELSE
               CALL MPYTR(IXI,IOMEGA,IWORKRNI,NII,R,1,NII)
            ENDIF
   
            CALL SCM(IWORKRNI,ERRINV,IWORKRNI,R,NII,0)
            CALL MPYM(IVARPO,IWORKRNI,IBE1,R,R,1,0,NII)
            CALL MPYM(IBE1,IWORKE,IDEV,R,NII,0,0,1)
            IF (NOMU .EQ. 0) THEN
               CALL ADDM(IDEV,IMU,IBETA,R,1,0)
            ELSEIF (NOMU .EQ. 1) THEN
               CALL RELOC(IDEV,IBETA,R,1,0)
            ENDIF
            
            ! If final iteration (IFIN = 2 or 3) then print out the final
            ! individual statistics to FILERES
   
            IF (IFIN .GE. 2 .AND. IRES .EQ. 1) THEN
               WRITE(5,38)IDI,NII
               WRITE(5,338)(IBETA(H),H=1,R)
               WRITE(5,338)(IVARPO(HR),HR=1,RRA)
               38 FORMAT(2I15)
               338 FORMAT(5F15.6)
            ENDIF
   
            ! *******************************************************
            ! OBTAIN SOME NECESSARY VECTORS AND MATRICES
            ! USED BY LIKELIHOOD AND BOTH EM & FISHER
            ! *******************************************************
              
            !   VECTOR U - ESTIMATED DEVIATIONS
    
            CALL MPYM(IXI,IBETA,IWORKE,NII,R,0,0,1)
            IF (P .EQ. 0) THEN
               CALL SUBM(IYI,IWORKE,IUVEC,NII,1,0)
            ELSE         
               CALL SUBM(IYI,IWORKE,IYXB,NII,1,0)
               CALL SUBM(IYXB,IWAL,IUVEC,NII,1,0)
            ENDIF
         
         ELSE    ! R .LE. 0, NO random effects
            CALL MPYM(IWI,IALPHA,IWAL,NII,P,0,0,1)
            CALL SUBM(IYI,IWAL,IUVEC,NII,1,0)
            IF (NAUTO .EQ. 0) THEN
               DETE = DETE + XNII*DLOG(ERROR)
            ELSE
               CALL INVS(IOMEGA,NII,DET,IWORKE,IER)
               IF (DET .LE. 0.000000000000001D0) DET = 0.000000000000001D0
               DETE = DETE + (DLOG(DET) + XNII*DLOG(ERROR))
            ENDIF
            
         ENDIF
     
         ! *********************************************************
         ! END - COMPUTE POSTERIOR BETA VARIANCE COVARIANCE
         ! *********************************************************
    
         ! ERROR TERMS
    
         IF (R .GT. 0) THEN
            CALL GRAMM(IXI,IVARPO,IXVARX,NII,R,1,IWORKRR)
            CALL SCM(IXVARX,ERRINV,IEI,NII,NII,1)
            IF (NAUTO .EQ. 0) THEN
               ICID = 0
               DO ICII = 1,NII
                  DO ICI2 = 1,ICII
                     ICID = ICID + 1
                     IF (ICI2 .EQ. ICII) THEN
                        IEI(ICID) = 1.0D0 - IEI(ICID)
                     ELSE
                        IEI(ICID) = 0.0D0 - IEI(ICID)
                     ENDIF
                  END DO
               END DO
            ELSE
               CALL INVS(IOMEGA,NII,DET,IWORKE)
               CALL SUBM(IOMEGA,IEI,IEI,NII,NII,1)
            ENDIF
         ELSE
            IF (NAUTO .EQ. 0) THEN
    
               ! note that this next part builds an identity matrix as E
               ! which can be changed at some point to be more efficient
               ! if the equations are worked out without E for the
               ! case of r=0
    
               ICID = 0
               DO ICII = 1,NII
                  DO ICI2 = 1,ICII
                     ICID = ICID + 1
                     IF (ICI2 .EQ. ICII) THEN
                        IEI(ICID) = 1.0D0 
                     ELSE
                        IEI(ICID) = 0.0D0 
                     ENDIF
                  END DO
               END DO
            ELSE
               CALL INVS(IOMEGA,NII,DET,IWORKE)
               CALL RELOC(IOMEGA,IEI,NII,NII,1)
            ENDIF
         ENDIF
   
   
         CALL SCM(IEI,ERROR,IEI,NII,NII,1)
         CALL CHAMS(IEI,IEI,NII,1,0)
         IF (NAUTO .GT. 0) THEN
            CALL INVS(IOMEGA,NII,DET,IWORKE)
            CALL CHAMS(IOMEGA,IOMEGSQ,NII,1,0)
         ENDIF
         CALL GRAM(IUVEC,IUU,NII,1)
         CALL CHAMS(IUU,IUU,NII,1,0)
     
         ! W TRANSPOSED x ERROR INVERSE MATRIX - FOR ALPHA
    
         IF (P .GT. 0) THEN
            IF (NAUTO .EQ. 0) THEN
               CALL TRP(IWI,IWOM,NII,P)
            ELSE
               CALL MPYTR(IWI,IOMEGA,IWOM,NII,P,1,NII)
            ENDIF
         ENDIF
    
         ! LIKELIHOOD TERMS
    
         IF (NAUTO .EQ. 0) THEN 
            CALL TRAM(IUU,WORK1,NII,0)
         ELSE
            CALL MPYM(IOMEGSQ,IUU,WORK1,1,NII2,0,0,1)
         ENDIF
         
         ERLIK = ERLIK + WORK1(1)
         IF (R .GT. 0) THEN
            IF (NOCOV .EQ. 0) THEN
               CALL GRMMT(IDEV,IVARCOI,WORK1,R,1,1,IWORKRR)
            ELSEIF (NOCOV .EQ. 1) THEN
               CALL GRMMT(IDEV,IVARCOI,WORK1,R,1,2,IWORKRR)
            ENDIF
            VALIK = VALIK + WORK1(1)
         ENDIF
    
         ! AUTOCORRELATION TERMS
    
         IF (NAUTO .GT. 1) THEN
            CALL RELOC(IOMEGAD,IOMEGAD0,NINI,S,0)                
            CALL SYMTGG(IOMEGAD,NII,S)                        
         ENDIF
    
         ! *******************************************************
         ! For the first NEM iterations, do EM solutions, after that, Fischer
         ! *******************************************************
         
         IF (IT .LE. NEM) THEN
            
            ! ****************************************************************
            ! EM ALGORITHM
            ! 
            !  IMUEMGR()   = VECTOR (R)  SUM OF BETA OVER SUBJECTS
            !  IAL1EM()    = VECTOR (PP) SUM OF FIRST PART OF ALPHA EQUATION
            !  IAL2EM()    = VECTOR (P)  SUM OF SECOND PART OF ALPHA EQUATION
            !  VARCEM      = VECTOR (RR) SUM OF BETA SQUARED PLUS POSTERIOR VAR
            !  ERREM       = SCALAR      SUM FOR ERROR VARIANCE
            !  IAUTOEM()   = VECTOR (S)  SUM FOR AUTOCORRELATION PARAMETERS
            ! ****************************************************************
    
            IF (R .GT. 0) THEN
               IF (NOMU .EQ. 0) CALL ADDM(IMUEMGR,IBETA,IMUEMGR,R,1,0)
               CALL GRMCV(IVARPO,IWORKRR,IBETA,UNO,R)
               IF (NOCOV .EQ. 0) THEN
                  CALL ADDM(IVARCEM,IWORKRR,IVARCEM,R,R,1)
               ELSEIF (NOCOV .EQ. 1) THEN
                  CALL CHAMS(IWORKRR,IWORKRR,R,1,2)
                  CALL ADDM(IVARCEM,IWORKRR,IVARCEM,R,R,2)
               ENDIF
            ENDIF
    
            IF (P .GT. 0) THEN
               IF (NAUTO .EQ. 0) THEN
                  CALL GRAMT(IWI,IWORKPP,NII,P)
               ELSE
                  CALL GRMMT(IWI,IOMEGA,IWORKPP,NII,P,1,IWORKE)
               ENDIF
               CALL ADDM(IAL1EM,IWORKPP,IAL1EM,P,P,1)
               IF (R .GT. 0) THEN
                  CALL MPYM(IWOM,IYXB,IWORKPP,P,NII,0,0,1)
               ELSE
                  CALL MPYM(IWOM,IYI,IWORKPP,P,NII,0,0,1)
               ENDIF
               CALL ADDM(IAL2EM,IWORKPP,IAL2EM,P,1,0)
            ENDIF
    
            IF (NAUTO .EQ. 0) THEN
               WORK1(1) = 0.0D0
               IF (R .GT. 0) CALL TRAM(IXVARX,WORK1,NII,1)
               CALL TRAM(IUU,IWORKE,NII,0)
               IWORKE(1) = IWORKE(1) + WORK1(1)
            ELSE
               IF (R .GT. 0) THEN
                  CALL CHAMS(IXVARX,IWRKNI2,NII,1,0)
                  CALL ADDM(IUU,IWRKNI2,IWORKNI2,NII,NII,0)
               ELSE
                  CALL RELOC(IUU,IWORKNI2,NII,NII,0)
               ENDIF
               CALL MPYM(IOMEGSQ,IWORKNI2,IWORKE,1,NII2,0,0,1)
            ENDIF
            CALL ADDM(ERREM,IWORKE,ERREM,1,1,0)
    
    
            IF (NAUTO .GT. 1) THEN
    
               CALL CHAMS(IWORKNI2,IWORKNI2,NII,0,1)
               CALL GRAMM(IOMEGSQ,IWORKNI2,IWRKNINI,NII,NII,1,IWORKE)
               CALL MPYTR(IOMEGAD,IWRKNINI,IWORKE,NINI,S,0,1)
               CALL ADDM(IAUTOEM,IWORKE,IAUTOEM,S,1,0)
    
               NOC = 0
               CALL AVECHA(IOMEGSQ,IOMEGAD0,IWRKNINS,NII,NII,S,IWORKE,NOC,NOCOV)
               CALL MPYTR(IOMEGAD,IWRKNINS,IWORKNI2,NINI,S,0,S)
               CALL CHAMS(IWORKNI2,IWORKNI2,S,0,1)
               CALL ADDM(IAUTOEM2,IWORKNI2,IAUTOEM2,S,S,1)
            ENDIF
    
         ELSE ! IT .GT. NEM
   
            ! **********************************************************
            ! FISHER SCORING 
            ! **********************************************************
    
            ! ********
            ! GRADIENT
            ! 
            !  IMUEMGR() = VECTOR (R) SUM OF V OVER SUBJECTS
            !  IALGR()   = VECTOR (P) SUM OF W x ERR INV x U FOR ALPHA
            !  VARGR     = VECTOR (RR) SUM OF BETA VARIANCE COMPONENTS
            !  ERRGR     = SCALAR SUM FOR ERROR VARIANCE
            !  IAUTOGR() = VECTOR (S) SUM FOR AUTOCORRELATION
            ! ********
    
            ! WRITE(6,*)'line 1712 A-OK'
            IF (R .GT. 0) THEN
               IF (NOMU .EQ. 0) CALL ADDM(IMUEMGR,IDEV,IMUEMGR,R,1,0)
               CALL GRMCV(IVARPO,IWORKRR,IDEV,UNO,R)
               IF (NOCOV .EQ. 0) THEN
                  CALL ADDM(IVARGR,IWORKRR,IVARGR,R,R,1)
               ELSEIF (NOCOV .EQ. 1) THEN
                  CALL CHAMS(IWORKRR,IWORKRR,R,1,2)
                  CALL ADDM(IVARGR,IWORKRR,IVARGR,R,R,2)
               ENDIF
            ENDIF
     
            IF (P .GT. 0) THEN
               CALL MPYM(IWOM,IUVEC,IWORKPP,P,NII,0,0,1)
               CALL ADDM(IALGR,IWORKPP,IALGR,P,1,0)
            ENDIF
    
    
            CALL SUBM(IUU,IEI,IUU,NII,NII,0)
            IF (NAUTO .EQ. 0) THEN 
               CALL TRAM(IUU,IWORKE,NII,0)
            ELSE
               CALL MPYM(IOMEGSQ,IUU,IWORKE,1,NII2,0,0,1)
            ENDIF
            CALL ADDM(ERRGR,IWORKE,ERRGR,1,1,0)
   
            ! autocorrelation equations
   
            IF (NAUTO .GT. 1) THEN
               CALL CHAMS(IUU,IUU,NII,0,1)
               CALL GRAMM(IOMEGSQ,IUU,IWRKNINI,NII,NII,1,IWORKE)
               CALL MPYTR(IOMEGAD,IWRKNINI,IWORKE,NINI,S,0,1)
               CALL ADDM(IAUTOGR,IWORKE,IAUTOGR,S,1,0)
            ENDIF
   
            ! ********
            ! INFORMATION MATRIX
            ! 
            !   IMU1IN()  = MATRIX (RRA)   SUM OF POSTERIOR VARIANCE COVARIANCE
            !   IALIN()   = MATRIX (PP)    SUM FOR ALPHA INFORMATION
            !   IALM1IN() = MATRIX (PR)    SUM FOR ALPHA x MU INFORMATION
            !   IVARCIN() = MATRIX (R x R) SUM OF HADAMARD PRODUCTS
            !   IERRINF() = SCALAR         SUM FOR ERROR INFORMATION
            !   IERV1IN() = MATRIX (RR) OF ERROR x BETA VARIANCE INFORMATION
            !   IAUTOIN() = MATRIX (SS)    SUM FOR AUTOCORRELATION INFORMATION
            !   IAUV1IN() = MATRIX (S x RR) OF AUTOCORR x BETA VAR INFORMATION
            !   IAUERIN() = VECTOR (S)     SUM FOR AUTOCORR x ERROR INFORMATION
            ! ********
    
            ! WRITE(6,*)'line 1757 A-OK'
            IF (R .GT. 0) THEN
               IF (NOMU .EQ. 0) CALL ADDM(IMU1IN,IVARPO,IMU1IN,R,R,1)
               CALL MPYM(IXI,IVARPO,IXVARP,NII,R,0,1,R)
            ENDIF
     
            IF (P .GT. 0) THEN
               CALL CHAMS(IEI,IEI,NII,0,1)
               CALL GRAMM(IWOM,IEI,IWORKPP,P,NII,1,IWORKE)
               CALL CHAMS(IEI,IEI,NII,1,0)
               CALL ADDM(IALIN,IWORKPP,IALIN,P,P,1)
   
               IF (R .GT. 0 .AND. NOMU .EQ. 0) THEN
                  CALL MPYM(IWOM,IXVARP,IWOMX,P,NII,0,0,R)
                  CALL ADDM(IALM1IN,IWOMX,IALM1IN,P,R,0)
               ENDIF
            ENDIF
    
            IF (R .GT. 0) THEN
               IF (NOCOV .EQ. 0) THEN
                  CALL RELOC(IVARCO,IWORKRR,R,R,1)
                  CALL SUBM(IWORKRR,IVARPO,IWORKRR,R,R,1)
                  CALL KMPY(IWORKRR,IWORKRR,IW3,R,R,1,R,R)
                  CALL ADDM(IVARIN,IW3,IVARIN,R2,R2,1)
               ELSEIF (NOCOV .EQ. 1) THEN
                  CALL CHAMS(IVARCO,IWORKRR,R,2,1)
                  CALL SUBM(IWORKRR,IVARPO,IWORKRR,R,R,1)
                  CALL HMPY(IWORKRR,IWORKRR,IW3,R,R,1)
                  CALL ADDM(IVARIN,IW3,IVARIN,R,R,1)
               ENDIF
                
               !  CALL PRNT(IUN,IVARCO,R,R,1,IXLAB,IXLAB,ND,HEAD,1,80,4,1,1, &
               ! "BETA VARIANCE COVARIANCE MATRIX                          ")
               !  CALL PRNT(IUN,IVARPO,R,R,1,IXLAB,IXLAB,ND,HEAD,1,80,4,1,1, &
               ! "POSTERIOR BETA VARIANCE COVARIANCE MATRIX                ")
    
            ENDIF
    
            IF (NAUTO .EQ. 0) THEN
                CALL GRAM(IEI,IWRKNINI,NII,NII)
                CALL TRAM(IWRKNINI,WORK1,NII,1)
            ELSE
                CALL CHAMS(IEI,IEI,NII,0,1)
                CALL GRAMM(IOMEGSQ,IEI,IAMAT,NII,NII,1,IWORKE)
                ! Used to use IWRKNINI as the work vector, but this
                ! was not big enough to hold the output!
                ! TG altered 1/29/01 to use IWORKNI2, which is big
                ! enough.
                CALL MPYM(IEI,IAMAT,IWORKNI2,NII,NII,1,1,NII)
                CALL TRAM(IWORKNI2,WORK1,NII,0)
            ENDIF
            
            CALL ADDM(IERRINF,WORK1,IERRINF,1,1,1)
    
            IF (R .GT. 0) THEN
               IF (NAUTO .EQ. 0) THEN
                   CALL GRAMT(IXVARP,IWORKR2,NII,R)
               ELSE
                   CALL GRMMT(IXVARP,IOMEGA,IWORKR2,NII,R,1,IWORKE)
               ENDIF
               IF (NOCOV .EQ. 1) CALL CHAMS(IWORKR2,IWORKR2,R,1,2)
               CALL ADDM(IERV1IN,IWORKR2,IERV1IN,1,RR,0)
            ENDIF
    
            !WRITE(6,*)'line 1804 A-OK'
            IF (NAUTO .GT. 1) THEN
               CALL CHAMS(IAMAT,IWORKNI2,NII,1,0)
               NOC = 0
               CALL AVECHA(IWORKNI2,IOMEGAD0,IWRKNINS,NII,NII, &
                           S,IWORKE,NOC,NOCOV)
               CALL MPYTR(IOMEGAD,IWRKNINS,IWRKNI2,NINI,S,0,S)
               CALL CHAMS(IWRKNI2,IWRKNI2,S,0,1)
               CALL ADDM(IAUTOIN,IWRKNI2,IAUTOIN,S,S,1)
    
               CALL INVS(IOMEGA,NII,DET,IWORKE)
               CALL GRAMM(IWORKNI2,IOMEGA,IWRKNINI,NII,NII,1,IWORKE)
               CALL MPYTR(IOMEGAD,IWRKNINI,IWORKE,NINI,S,0,1) 
               CALL ADDM(IAUERIN,IWORKE,IAUERIN,S,1,0)
    
               IF (R .GT. 0) THEN
                  CALL INVS(IOMEGA,NII,DET,IWORKE)
                  CALL MPYTR(IXVARP,IOMEGA,IWORKRNI,NII,R,1,NII)
                  NOC = NOCOV
                  CALL AVECHA(IWORKRNI,IOMEGAD0,IWORKSR,R,NII, &
                              S,IWORKE,NOC,NOCOV)
                  CALL ADDM(IAUV1IN,IWORKSR,IAUV1IN,RR,S,0)
               ENDIF
            ENDIF
         ENDIF  ! EM vs. Fischer estimation
         
      END DO ! end of single subject  DO I = 1,N   
 
      ! **********************************************************
      ! DONE WITH SUBJECTS  COMPUTE VARCEM AND FINISH EM SOLUTIONS
      ! **********************************************************

      ! WRITE(6,*)'line 1832 A-OK'
      IF (IT .LE. NEM) THEN
         
         ! *******************************************************
         ! First NEM iterations, do EM solutions, after that, Fischer
         ! *******************************************************
 
         IF (P .GT. 0) THEN
            CALL INVS(IAL1EM,P,DET,IWORKPP)
            CALL MPYM(IAL1EM,IAL2EM,IALPHA,P,P,1,0,1)
         ENDIF
 
         IF (R .GT. 0) THEN
            
            IF (NOCOV .EQ. 0) THEN
               CALL SCM(IVARCEM,RNINV,IVARCEM,R,R,1)
            ELSEIF (NOCOV .EQ. 1) THEN
               CALL SCM(IVARCEM,RNINV,IVARCEM,R,R,2)
            ENDIF
            
            IF (NOMU .EQ. 0) THEN
               CALL SCM(IMUEMGR,RNINV,IMU,R,1,0)
               IF (NOCOV .EQ. 0) THEN
                  CALL GRMCV(IVARCEM,IVARCO,IMU,UNON,R)
               ELSEIF (NOCOV .EQ. 1) THEN
                  CALL GRAMD(IMU,IWORKRR,R,1)
                  CALL SUBM(IVARCEM,IWORKRR,IVARCO,R,R,2)
               ENDIF
            ELSEIF (NOMU .EQ. 1) THEN
               IF (NOCOV .EQ. 0) THEN
                  CALL RELOC(IVARCEM,IVARCO,R,R,1)
               ELSEIF (NOCOV .EQ. 1) THEN
                  CALL RELOC(IVARCEM,IVARCO,R,R,2)
               ENDIF
            ENDIF
         ENDIF
 
         ERROR = ERREM(1) / DBLE(NTOT)
 
         IF (NAUTO .GT. 1) THEN
            IF (NS .LE. 2) AUTOOLD = IAUTO(1)
            CALL INVS(IAUTOEM2,S,DET,IWORKE)
            CALL MPYM(IAUTOEM2,IAUTOEM,IAUTO,S,S,1,0,1)
            CALL SCM(IAUTO,ERRINV,IAUTO,S,1,0)
            IF (NS .LE. 2) THEN
               IAUTO(1) = IAUTO(1) / (1.0D0 - AUTOOLD*AUTOOLD)
               IAUTO(1) = (DSQRT(1.0D0 + 4.0D0 * IAUTO(1)*IAUTO(1)) - &
                          1.0D0) / (2.0D0 * IAUTO(1))
            ENDIF
         ENDIF

      ELSE ! IT .GT. NEM
 
         ! **********************************************************
         ! FINISH FISHER SCORING RESULTS
         ! **********************************************************

         ! WRITE(6,*)'line 1866 A-OK'
         IF (R .GT. 0) THEN
            IF (NOCOV.EQ.0) CALL SCM(IVARCO,REALN,IVARCON,R,R,1)
            IF (NOCOV.EQ.1) CALL SCM(IVARCO,REALN,IVARCON,R,R,2)
 
            IF (NOMU .EQ. 0 .AND. NOCOV .EQ. 0) THEN
               CALL MPYM(IVARCOI,IMUEMGR,IGRAFIX,R,R,1,0,1)
               CALL SUBM(IVARCON,IMU1IN,IMU1IN,R,R,1)
               CALL GRMMT(IVARCSQ,IMU1IN,IMUIN,R,R,1, IWORKRR)
            ELSEIF (NOMU .EQ. 0 .AND. NOCOV .EQ. 1) THEN
               CALL MPYM(IVARCOI,IMUEMGR,IGRAFIX,R,R,2,0,1)
               CALL CHAMS(IVARCON,IWORKRR,R,2,1)
               CALL SUBM(IWORKRR,IMU1IN,IMU1IN,R,R,1)
               CALL MPDSD(IVARCOI,IMU1IN,IMUIN,R)
            ENDIF
         ENDIF
 
         IF (P .GT. 0 .AND. R .GT. 0 .AND. NOMU .EQ. 0) THEN         
            CALL SCM(IALGR,ERRINV,IALGR,P,1,0)
            IF (NOCOV .EQ. 0) CALL MPYM(IALM1IN,IVARCOI,IALMUIN,P,R,0,1,R)
            IF (NOCOV .EQ. 1) CALL MPYM(IALM1IN,IVARCOI,IALMUIN,P,R,0,2,R)
            CALL SCM(IALMUIN,ERRINV,IALMUIN,P,R,0)
            QUERR = ERRINV * ERRINV
            CALL SCM(IALIN,QUERR,IALIN,P,P,1)

            CALL ADJR(IGRAFIX,IALGR,IGRAFIX,R,1,P)

            CALL PARTIV(IMUIN,R,RRA,IALMUIN,RP,IALIN,P,PP, &
                 IINFFIX,NFNF,IWORKRR,IALM1IN,IWORKPP,IER)

         ELSEIF (P .LE. 0 .AND. R .GT. 0 .AND. NOMU .EQ. 0) THEN
            CALL INVS(IMUIN,R,DETR,IWORKR) 
            CALL RELOC(IMUIN,IINFFIX,R,R,1)
         
         ELSEIF (P .GT. 0 .AND. (R .LE. 0 .OR. NOMU .EQ. 1)) THEN
            CALL SCM(IALGR,ERRINV,IGRAFIX,P,1,0)
            QUERR = ERRINV * ERRINV
            CALL SCM(IALIN,QUERR,IALIN,P,P,1)
            CALL INVS(IALIN,P,DETP,IWORKP)
            CALL RELOC(IALIN,IINFFIX,P,P,1)
         ENDIF
 
         CALL MPYM(IINFFIX,IGRAFIX,ICORFIX,NF,NF,1,0,1)
          
         ! CALL PRNT(2,IGRAFIX,NF,1,0,IXLAB,IXLAB,ND,HEAD,1,80,6,1,1, &
         !    "Gradient of fixed terms                                  ")
         ! CALL PRNT(2,IINFFIX,NF,NF,1,IXLAB,IXLAB,ND,HEAD,1,80,6,1,1, &
         !    "Information matrix of fixed terms                        ")
 
         IF (IFIN .LT. 2) THEN
 
            ! WRITE(6,*)'line 1902 A-OK'
            BIGCOR = 0.0D0
            CORSUM = 0.0D0
            ICSUM  = 0
             
            IF (R .GT. 0 .AND. NOMU .EQ. 0) THEN
               DO H = 1,R
                  IF (DABS(ICORFIX(H)) .GT. BIGCOR) BIGCOR = DABS(ICORFIX(H))
                  CORSUM = CORSUM + DABS(ICORFIX(H))
                  ICSUM  = ICSUM  + 1
                  IMU(H) = IMU(H) + ICORFIX(H)
               END DO
            ENDIF
    
            IF (P .NE. 0) THEN
               DO L = 1,P
                  IP = (1-NOMU)*R + L
                  IF (DABS(ICORFIX(IP)) .GT. BIGCOR) BIGCOR = DABS(ICORFIX(IP))
                  CORSUM = CORSUM + DABS(ICORFIX(IP))
                  ICSUM  = ICSUM  + 1
                  IALPHA(L) = IALPHA(L) + ICORFIX(IP)
               END DO
            ENDIF
         ENDIF
 
         IF (R .GT. 0) THEN
            IF (NOCOV .EQ. 0) THEN
               CALL SUBM(IVARGR,IVARCON,IWORKR2,R,R,1)
               CALL CHAMS(IVARGR,IVARGR,R,1,0)
               CALL CHAMS(IWORKR2,IWORKR2,R,1,0)
            ELSEIF (NOCOV .EQ. 1) THEN
               CALL SUBM(IVARGR,IVARCON,IWORKR2,R,R,2)
            ENDIF
         ENDIF
 
         ! finish information part for either reparameterized (ifin=0)
         ! or original solution (ifin > 0)
 
         IF (IFIN .NE. 0) THEN
            IF (R .GT. 0) THEN
               IF (NOCOV .EQ. 0) THEN          
                  CALL KMPY(IVARCOI,IVARCOI,IVARCK,R,R,1,R,R)
                  CALL CHAMS(IVARCK,IVARCK,R2,1,0)
                  CALL MPYTR(ITRANG,IVARCK,IW3,R2,RR,0,R2)
                  !
                  ! the next two lines were added 5/93
                  !
                  CALL MPYM(IW3,ITRANG,IW3G,RR,R2,0,0,RR)
                  CALL CHAMS(IW3G,IW3G,RR,0,1)
  
                  CALL MPYM(IW3,IWORKR2,IGRAVAR,RR,R2,0,0,1)
                  CALL SCM(IGRAVAR,0.5D0,IGRAVAR,RR,1,0)
 
                  CALL GRAMM(IW3,IVARIN,IVARCIN,RR,R2,1,IWORKR2)
                  CALL SCM(IVARCIN,0.5D0,IVARCIN,RR,RR,1)
 
               ELSEIF (NOCOV .EQ. 1) THEN
                  CALL HMPY(IVARCOI,IWORKR2,IWORKR2,R,R,2)
                  CALL HMPY(IWORKR2,IVARCOI,IGRAVAR,R,R,2)
                  CALL SCM(IGRAVAR,0.5D0,IGRAVAR,R,R,2)
                  CALL HMPY(IVARCOI,IVARCOI,IWORKR2,R,R,2)
                  CALL MPDSD(IWORKR2,IVARIN,IVARCIN,R)
                  CALL SCM(IVARCIN,0.5D0,IVARCIN,R,R,1)
               ENDIF
            ENDIF
 
            QUERR = 0.5D0 * ERRINV * ERRINV
            ERRGR(1)  =   ERRGR(1) * QUERR
 
            IF (R .GT. 0) THEN
               CALL ADJR(IGRAVAR,ERRGR,IGRAVAR,RR,1,1)
 
               IF (NOCOV .EQ. 0) THEN
                  CALL MPYM(IERV1IN,IW3G,IERVAIN,1,RR,0,1,RR)
               ELSEIF (NOCOV .EQ. 1) THEN
                  CALL MPYM(IERV1IN,IWORKR2,IERVAIN,1,R,0,2,R)
                  ! remember r = rr if NOCOV=1
               ENDIF
               CALL SCM(IERVAIN,QUERR,IERVAIN,1,RR,0)
            ELSEIF (R .EQ. 0) THEN
               IGRAVAR(1) = ERRGR(1)
            ENDIF
 
            QUERR = QUERR * ERRINV * ERRINV
            IERRINF(1) =  IERRINF(1) * QUERR
 
            ! WRITE(6,*)'line 1947 A-OK'
            IF (NAUTO .GT. 1) THEN
               QUERR = 0.5D0 * ERRINV 
               CALL SCM(IAUTOGR,QUERR,IAUTOGR,S,1,0)
               IR = RR+1 
               CALL ADJR(IGRAVAR,IAUTOGR,IGRAVAR,IR,1,S)
               QUERR = QUERR * ERRINV 
               CALL SCM(IAUTOIN,QUERR,IAUTOIN,S,S,1)
               QUERR = QUERR * ERRINV
               CALL SCM(IAUERIN,QUERR,IAUERIN,S,1,0)
               CALL ADJRC(IERRINF,IAUERIN,IAUTOIN,IERRINF,S,1)
 
               IF (R .GT. 0) THEN
                  QUERR = 0.5D0 * ERRINV
                  ! multiply by the transpose of IAUV1IN() 
                  IF (NOCOV .EQ. 0) THEN
                     CALL MPYTR(IAUV1IN,IW3G,IAUTVIN,RR,S,1,RR)
                  ELSEIF (NOCOV .EQ. 1) THEN
                     CALL MPYTR(IAUV1IN,IWORKR2,IAUTVIN,R,S,2,R)
                  ENDIF
                  CALL SCM(IAUTVIN,QUERR,IAUTVIN,S,RR,0)
                  CALL ADJR(IERVAIN,IAUTVIN,IERVAIN,1,RR,S)
               ENDIF
 
            ENDIF
 
         ! the following is for the reparameterized solution (ifin = 0)
 
         ELSE   ! IFIN .EQ. 0
         
            IF (NOCOV .EQ. 0) THEN
               CALL RELOC(IVARL,IVARLSQ,R,R,3)
               CALL INVT(IVARLSQ,R,DET)
               CALL CHAMS(IVARLSQ,IVARLSQ,R,3,0)
 
               CALL INVD(IVARD,R,DET)
               CALL KMPY(IVARLSQ,IVARLSQ,IVARCK,R,R,0,R,R)
               CALL MPYM(IVARD,ITRANH,IWORKR3,R,R,2,0,R2)
               CALL MPYM(IWORKR3,IVARCK,IW3,R,R2,0,0,R2)
 
               CALL MPYM(IW3,IWORKR2,IGRAVAR,R,R2,0,0,1)
               CALL SCM(IGRAVAR,0.5D0,IGRAVAR,R,1,0)
            ELSEIF (NOCOV .EQ. 1) THEN
               CALL INVD(IVARD,R,DET)
               CALL HMPY(IVARD,IWORKR2,IGRAVAR,R,R,2)
               CALL SCM(IGRAVAR,0.5D0,IGRAVAR,R,R,2)
            ENDIF


            IF (NOCOV .EQ. 0) THEN
               ! ALLOCATE  r * rr            where  rr = (r * (r+1)) / 2
               IRRR  = R * RR
               
               IF (ALLOCATED(IWORKRRR)) THEN
                  IF(UBOUND(IWORKRRR,1) .NE. IRRR) THEN
                     DEALLOCATE(IWORKRRR)
                     ALLOCATE(IWORKRRR(IRRR))
                  ENDIF
               ELSE
                  ALLOCATE(IWORKRRR(IRRR))
               ENDIF
         
               CALL MPYM(IW3,ITRANG,IWORKRRR,R,R2,0,0,RR) 
               CALL TRP(IWORKRRR,IW3G,R,RR)
            
               CALL MPYRT(IVARIN,IW3,IWORKR3,R2,R2,1,R)
 
               CALL GRAMM(IW3,IVARIN,IVARCIN,R,R2,1,IWORKR2)
               CALL SCM(IVARCIN,0.5D0,IVARCIN,R,R,1)
            ELSEIF (NOCOV .EQ. 1) THEN
               CALL MPDSD(IVARD,IVARIN,IVARCIN,R)
               CALL SCM(IVARCIN,0.5D0,IVARCIN,R,R,1)
            ENDIF
 
            QUERR = 0.5D0 * ERRINV 
            IF (NOCOV .EQ. 0) THEN
               CALL MPYM(IERV1IN,IW3G,IERVAIN,1,RR,0,0,R)
            ELSEIF (NOCOV .EQ. 1) THEN
               CALL MPYM(IERV1IN,IVARD,IERVAIN,1,R,0,2,R)
               ! remember r = rr if NOCOV=1
            ENDIF
            
            CALL SCM(IERVAIN,QUERR,IERVAIN,1,R,0)
            ERRGR(1)  =  ERRGR(1) * QUERR
            QUERR = QUERR * ERRINV 
            IERRINF(1) = IERRINF(1) * QUERR
 
            !     WRITE(6,*)'line 2008 A-OK'
            IF (NAUTO .GT. 1) THEN
               QUERR = 0.5D0 * ERRINV 
               CALL SCM(IAUTOGR,QUERR,IAUTOGR,S,1,0)
               QUERR = QUERR * ERRINV 
               CALL SCM(IAUTOIN,QUERR,IAUTOIN,S,S,1)
               CALL SCM(IAUERIN,QUERR,IAUERIN,S,1,0)
               CALL ADJRC(IERRINF,IAUERIN,IAUTOIN,IERRINF,S,1)
  
               QUERR = 0.5D0 * ERRINV
               IF (NOCOV .EQ. 0) THEN
                  CALL MPYTR(IAUV1IN,IW3G,IAUTVIN,RR,S,0,R)
               ELSEIF (NOCOV .EQ. 1) THEN
                  CALL MPYTR(IAUV1IN,IVARD,IAUTVIN,R,S,2,R)
               ENDIF
               CALL SCM(IAUTVIN,QUERR,IAUTVIN,S,R,0)
  
            ENDIF
 
            ! Finished with IW3() as VARD*ITRANH()*IVARCK() 
            ! for diagonal PI terms.  
            ! Now put ITRANJ()*IVARCK() into IW3()
 
            IF (NOCOV .EQ. 0) THEN
               CALL KMPY(IVARLSQ,IVARCSQ,IVARCK,R,R,0,R,R)
               CALL MPYM(ITRANJ,IVARCK,IW3,RUN,R2,0,0,R2)
 
               ! the next TWO lines were added 5/93

               ! ALLOCATE  run * rr , where  rr = (r * (r+1)) / 2

               IRRR  = RUN * RR
               IF (ALLOCATED(IWORKRRR)) THEN
                  IF(UBOUND(IWORKRRR,1) .NE. IRRR) THEN
                     DEALLOCATE(IWORKRRR)
                     ALLOCATE(IWORKRRR(IRRR))
                  ENDIF
               ELSE
                  ALLOCATE(IWORKRRR(IRRR))
               ENDIF
         
               CALL MPYM(IW3,ITRANG,IWORKRRR,RUN,R2,0,0,RR)           
               CALL TRP(IWORKRRR,IW3G,RUN,RR)
  
               CALL MPYM(IW3,IVARGR,IWORKRR,RUN,R2,0,0,1)
               CALL ADJR(IGRAVAR,IWORKRR,IGRAVAR,R,1,RUN)
            ENDIF 
            
            CALL ADJR(IGRAVAR,ERRGR,IGRAVAR,RR,1,1)

            IF (NOCOV .EQ. 0) THEN
               CALL MPYM(IW3,IWORKR3,IVARLD,RUN,R2,0,0,R)
  
               CALL MPYM(IVARIN,ITRANN,IVARCK,R2,R2,1,0,R2)
               CALL MPYM(IW3,IVARCK,IVARIN,RUN,R2,0,0,R2)
               CALL MPYRT(IVARIN,IW3,IVARCK,RUN,R2,0,RUN)
               CALL CHAMS(IVARCK,IVARCK,RUN,0,1)
               CALL SCM(IVARCK,2.0D0,IVARCK,RUN,RUN,1)
               CALL ADJRC(IVARCIN,IVARLD,IVARCK,IVARCIN, RUN,R)

               ! WRITE(6,*)'line 2056 A-OK'
               CALL MPYM(IERV1IN,IW3G,IWORKRR,1,RR,0,0,RUN)
               CALL SCM(IWORKRR,ERRINV,IWORKRR,1,RUN,0)
               CALL ADJC(IERVAIN,IWORKRR,IERVAIN,1,R,RUN)
            ENDIF

            IF (NAUTO .GT. 1) THEN
               IF (NOCOV .EQ. 0) THEN
                  CALL MPYTR(IAUV1IN,IW3G,IAUTLIN,RR,S,0,RUN)
                  CALL SCM(IAUTLIN,ERRINV,IAUTLIN,S,RUN,0)
               ENDIF
 
               IRR  = RR + 1
               CALL ADJR(IGRAVAR,IAUTOGR,IGRAVAR,IRR,1,S)
               IF (NOCOV .EQ. 0) CALL ADJC(IAUTVIN,IAUTLIN,IAUTVIN,S,R,RUN)
               CALL ADJR(IERVAIN,IAUTVIN,IERVAIN,1,RR,S)
 
            ENDIF
         ENDIF

         !  CALL PRNT(IUN,IVARCIN,RR,RR,1,IXLAB,IXLAB,ND,HEAD,1,80,6,1,1,
         ! +"INFORMATION MATRIX OF BETA VARIANCE                      ")
         !  CALL PRNT(IUN,IERVAIN,NERR,RR,0,IXLAB,IXLAB,ND,HEAD,1,80,6,1,1,
         ! +"INFORMATION MATRIX OF ERROR * BETA VARIANCE              ")
         !  CALL PRNT(IUN,IERRINF,NERR,NERR,1,IXLAB,IXLAB,ND,HEAD,1,80,6,1,1,
         ! +"INFORMATION MATRIX OF ERROR VARIANCE (tau & auto)        ")
         !  WRITE(6,*)'line 2085 A-OK'

         IF (R .GT. 0) THEN
            CALL PARTIV(IVARCIN,RR,RRU,IERVAIN,NEC,IERRINF, &
                 NERR,NERR2,IINFVAR,NVNV,IW3,IWKEC,IWKER2,IER)
         ELSEIF (R .EQ. 0) THEN
            ! note using iwker2 as work vector but only need 
            ! work vector of size nerr
            CALL INVS(IERRINF,NERR,DETP,IWKER2)
            CALL RELOC(IERRINF,IINFVAR,NERR,NERR,1)
         ENDIF

         CALL MPYM(IINFVAR,IGRAVAR,ICORVAR,NV,NV,1,0,1)

         !  CALL PRNT(IUN,IGRAVAR,NV,1,0,IXLAB,IXLAB,ND,HEAD,1,80,6,1,1,
         ! +"Gradient of variance terms                               ")
         !  CALL PRNT(IUN,IINFVAR,NV,NV,1,IXLAB,IXLAB,ND,HEAD,1,80,6,1,1,
         ! +"Information matrix of variance terms                     ")
 
         IF(IFIN .NE. 0 .AND. IFIN .LT. 2) THEN
 
            IF (R .GT. 0) THEN
               DO HR = 1,RR
                  IF (DABS(ICORVAR(HR)) .GT. BIGCOR) BIGCOR = DABS(ICORVAR(HR))
                  CORSUM = CORSUM + DABS(ICORVAR(HR))
                  ICSUM  = ICSUM  + 1
                  IVARCO(HR) = IVARCO(HR) + ICORVAR(HR)
               END DO
            ENDIF
  
            IR = RR + 1
            IF (DABS(ICORVAR(IR)) .GT. BIGCOR) BIGCOR = DABS(ICORVAR(IR))
            CORSUM = CORSUM + DABS(ICORVAR(IR))
            ICSUM  = ICSUM  + 1
            ERROR = ERROR + ICORVAR(IR)
            IF (ERROR .LE. 0.000000000000001D0) ERROR = 0.000000000000001D0
     
  
            IF (NAUTO .GT. 1) THEN
               DO IS = 1,S
                  IR  = RR + 1 + 1
                  IF (DABS(ICORVAR(IR)) .GT. BIGCOR) BIGCOR = DABS(ICORVAR(IR))
                  CORSUM = CORSUM + DABS(ICORVAR(IR))
                  ICSUM  = ICSUM  + 1
                  IAUTO(IS)= IAUTO(IS) + ICORVAR(IR)
                  IF (DABS(IAUTO(IS)) .GE. 1.0D0) IAUTO(IS) = 0.95
               END DO
            ENDIF
            ! WRITE(6,*)'line 2120 A-OK'
         
         ELSEIF (IFIN .EQ. 0) THEN
            DO H = 1,R
               IF (DABS(ICORVAR(H)) .GT. BIGCOR) BIGCOR = DABS(ICORVAR(H))
               CORSUM = CORSUM + DABS(ICORVAR(H))
               ICSUM  = ICSUM  + 1
               IVARPI(H) = IVARPI(H) + ICORVAR(H)
            END DO
            DO H = 1,R
               IVARD(H) = DEXP(IVARPI(H))
            END DO
 
            IF (NOCOV .EQ. 0) THEN           
 
               !  CALL PRNT(IUN,IVARL,R,R,3,IXLAB,IXLAB,ND,HEAD,1,80,6,1,1,
               ! +"L MATRIX  (before improvement)                           ")

               L = 0
               DO IH = 1,R
                  IR = R
                  DO IJ = 1,IH
                     L = L+1
                     IF (IJ .LT. IH) THEN
                        IF (IJ .EQ. 1) THEN
                           IR = IR + (IH - IJ)
                        ELSE
                           IR = IR + (R  - IJ)
                        ENDIF
                        IF (DABS(ICORVAR(IR)) .GT. BIGCOR) &
                             BIGCOR = DABS(ICORVAR(IR))
                        CORSUM = CORSUM + DABS(ICORVAR(IR))
                        ICSUM  = ICSUM  + 1
                        IVARL(L) = IVARL(L) + ICORVAR(IR)
                     ENDIF
                  END DO
               END DO
            
               !  CALL PRNT(IUN,IVARL,R,R,3,IXLAB,IXLAB,ND,HEAD,1,80,6,1,1,
               ! +"L MATRIX  (after  improvement)                           ")
 
            ENDIF
 
            IR = RR + 1
            IF (DABS(ICORVAR(IR)) .GT. BIGCOR) BIGCOR = DABS(ICORVAR(IR))
            CORSUM = CORSUM + DABS(ICORVAR(IR))
            ICSUM  = ICSUM  + 1
            TAU   = TAU   + ICORVAR(IR)
            ERROR = DEXP(TAU)
 
            ! WRITE(6,*)'line 2159 A-OK'
            IF (NAUTO .GT. 1) THEN
               IF (NS .LE. 2) THEN
                  IR  = RR + 1 + 1
                  IF (DABS(ICORVAR(IR)) .GT. BIGCOR) BIGCOR = DABS(ICORVAR(IR))
                  CORSUM = CORSUM + DABS(ICORVAR(IR))
                  ICSUM  = ICSUM  + 1
                  PHI = PHI + ICORVAR(IR)
                  IAUTO(1) = PHI / (DSQRT(1.0D0 + PHI*PHI))
                  DPHI= 1.0D0 / ((1.0D0 + PHI*PHI) * (DSQRT(1.0D0 + PHI*PHI)))
               ELSE
                  IR = RR + 1 
                  DO IS = 1,S
                     IR  = IR + 1
                     IF (DABS(ICORVAR(IR)) .GT. BIGCOR) BIGCOR = DABS(ICORVAR(IR))
                     CORSUM = CORSUM + DABS(ICORVAR(IR))
                     ICSUM  = ICSUM  + 1
                     IAUTO(IS)= IAUTO(IS) + ICORVAR(IR)
                    
                     !WRITE(6,*)'line 2215 A-OK'
                     IF (DABS(IAUTO(IS)) .GE. 1.0D0) IAUTO(IS) = 0.95
                  END DO
               ENDIF
            ENDIF
           
            ! WRITE(6,*)'line 2221 A-OK'
            IF (NOCOV .EQ. 0) THEN
               CALL CHAMS(IVARL,IVARLSQ,R,3,0)
               CALL GRAMM(IVARLSQ,IVARD,IVARCO,R,R,2,IWORKRR)
            ELSEIF (NOCOV .EQ. 1) THEN
               CALL RELOC(IVARD,IVARCO,R,R,2) 
            ENDIF
            
            !  WRITE(6,*)'line 2223 A-OK'
            !  CALL PRNT(IUN,IVARL,R,R,3,IXLAB,IXLAB,ND,HEAD,1,80,6,1,1,
            ! +"L MATRIX  (after  improvement)                           ")
            !  CALL PRNT(IUN,IVARLSQ,R,R,0,IXLAB,IXLAB,ND,HEAD,1,80,6,1,1,
            ! +"L MATRIX  (square)"  ")
            !  CALL PRNT(IUN,IVARD,R,R,2,IXLAB,IXLAB,ND,HEAD,1,80,6,1,1,
            ! +"D MATRIX"            ")
            !  WRITE(6,*)'line 2226 A-OK'
         
         ENDIF  !  IFIN .NE. 0 test
      ENDIF  ! IT .LE. NEM test
      
      ! WRITE(6,*)'line 2186 A-OK'
      
      ! COMPUTE LIKELIHOOD FOR THIS ITERATION
 
      IF (R .GT. 0) THEN
         IF (DETV .LE. 0.000000000000001D0) DETV = 0.000000000000001D0
         RLOGL = -0.5D0 * (PLIK + DETE + REALN*DLOG(DETV) + ERLIK/ERROR &
                + VALIK - DETB)
      ELSE
         RLOGL = -0.5D0 * (PLIK + DETE + ERLIK/ERROR)
      ENDIF
      
      ! PRINT OUT RESULTS FOR THIS ITERATION
      ! unless iteration was for computation of information matrix only
      ! (IFIN = 2)
 
      IF (IFIN .GE. 2) THEN  ! EXIT the main iteration loop
         MainLoop = 0
         CALL FREE_LOCALS()
         INITIALIZED = 2
         RETURN     ! escape from the loop
      ENDIF
 
      WRITE(IUN,65)IT
      65 FORMAT(/,1X,'**************',/,&
                  1X,'ITERATION ',I4,/, &
                  1X,'**************',/)
      WRITE(IUN,66) RLOGL,ERROR,(IAUTO(IS),IS=1,S)
      66 FORMAT(1X,'Log Likelihood  = ',F12.3,//, &
                1X,'Error Variance  = ',F12.6,/, &
                1X,'Autocorrelation = ',F12.6,6F10.6)
 
      IF (R .GT. 0 .AND. NOMU .EQ. 0) THEN
         WRITE(IUN,903) (IMU(H),H=1,R)
         903 FORMAT(1X,'MU        ',7F10.6)
      ENDIF
      
      IF (P .GT. 0) THEN
          WRITE(IUN,904) (IALPHA(L),L=1,P)
          904 FORMAT(1X,'ALPHA     ',7F10.6)
      ENDIF
      
      IF (R .GT. 0) THEN
         WRITE(IUN,905) (IVARCO(HR),HR=1,RR)
         905 FORMAT(1X,'BETA VAR  ',7F10.6)
      ENDIF

      ! check for positive definiteness of Toeplitz errors (MODID = 1)

      IF (NS .EQ. 5) THEN
            
         MODID = 0
         DO WHILE (MODID .EQ. 0)
            MODID = 1
            CALL TOECHECK(IAUTO,S,MAXNI,IWORKNI2,IWRKNI2, &
                         MODID,iun,nd,IXLAB,head)
            IF (MODID .EQ. 0) THEN
               CALL AUTOSPEC(IAUTO,IWORKNI2,IWRKNI2,S,NSPEC,NSPECS)
               WRITE(IUN,959) (IAUTO(IS),IS=1,S)
               
               959 FORMAT(/,1X,'===> The toeplitz autocorrelation matrix is not positive definite',/, &
                            1X,'===> it has been smoothed via spectral methods', &
                         //,1X,'     The smoothed autocorrelations are',/,1X,5F12.5)
            ELSE
               WRITE(IUN,960)
               960 FORMAT(/,1X,'===> The toeplitz autocorrelation matrix is positive definite')
     
            ENDIF
         END DO
      ENDIF
 
      ! WRITE(6,*)'line 2245 A-OK'
      ! if still on EM iterations go back for more 
 
      IF (IT .LE. NEM) THEN
         MainLoop = 1
         RETURN     ! escape from the loop
      ENDIF
 
      ! For Fisher Scoring iterations
      ! write out the corrections and then
      ! check to see if there are numerical problems and if so
      ! reduce the number of random-effects by 1

      ! WRITE(IUN,906) (IVARD(H),H=1,R)
      ! 906 FORMAT(1X,'BVAR DIAG ',7F10.6)
      WRITE(IUN,901) (ICORFIX(IN),IN=1,NF)
      901 FORMAT(/,1X,'CORR FIX  ',7F10.6)
      WRITE(IUN,902) (ICORVAR(IN),IN=1,NV)
      902 FORMAT(1X,'CORR VAR  ',7F10.6)
 
      ! save current avecor to check at next iteration
 
      AVECOR   = DABS(CORSUM / DBLE(ICSUM))
      AVECORD  = (AVECOR - AVECORP) / AVECORP
      AVECORP  = AVECOR

      WRITE(IUN,900) AVECOR,AVECORD
      900 FORMAT(/,1X,'AVE CORR  ',F10.6,/,1X,'AVECOR DIF',F10.6)
 
      IF(AVECORD .GT. 0.4D0 .OR. IER .EQ. 1) IRBAD = IRBAD + 1
      IF(IRBAD .GE. 1) THEN
         IRBADC = 1
         IRBADN = IRBADN + 1
         WRITE(6,878)IT
         WRITE(2,879)IT
         878 FORMAT(///,1X,'==> WARNING!  Estimation Difficulties at Iteration',I4,/, &
                        1x,'==> will restart iterations with one less random effect')
         879 FORMAT(///,1X,'==> WARNING!  Estimation Difficulties Occured at Iteration',I4,/, &
                        1x,'==> The model was refit with one less random effect than was requested',/, &
                        1x,'==> and the iterations were restarted after this point')
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
         IF (IFIN .EQ. 2) IGO = 1
         IF (IFIN .EQ. 3 .AND. IRES .EQ. 1) IGO = 1
         IF (IFIN .EQ. 3 .AND. IRES .EQ. 0) IGO = 0
      ENDIF
      
      ! Signal caller that there's more to do
      MainLoop = 1
   ELSE
      ! Signal caller that we're done iterations
      MainLoop = 0
      CALL FREE_LOCALS()
      ! free up arrays only used in the main loop, now that 
      ! we are finished with them.  Note that the SAVE statement
      ! preserves the values in these arrays between calls to
      ! the MainLoop function
      INITIALIZED = 2
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
   
   IF(ALLOCATED(IBETA)) DEALLOCATE(IBETA, IDEV, IVARPO)
   
   IF(ALLOCATED(IEI)) DEALLOCATE(IEI, IUU) 
   
   IF(ALLOCATED(IOMEGA)) DEALLOCATE(IOMEGA)
   
   IF (ALLOCATED(IOMEGAD)) THEN
      DEALLOCATE(IOMEGAD)
      DEALLOCATE(IOMEGAD0)
      DEALLOCATE(IWRKNINS)
   ENDIF
   IF (ALLOCATED(IOMEGSQ)) THEN
      DEALLOCATE(IOMEGSQ)
      DEALLOCATE(IWORKNI2)
      DEALLOCATE(IWRKNI2)
   ENDIF
   IF (ALLOCATED(IXTI)) THEN
      DEALLOCATE(IXTI)
      DEALLOCATE(IWORKE)
      DEALLOCATE(IXMU)
      DEALLOCATE(IWAL)
      DEALLOCATE(IYXB)
      DEALLOCATE(IUVEC)
   ENDIF
   IF (ALLOCATED(IWI)) DEALLOCATE(IWI)
   IF (ALLOCATED(IWOM)) DEALLOCATE(IWOM)
   IF (ALLOCATED(IBE1)) THEN
      DEALLOCATE(IBE1)
      DEALLOCATE(IWORKRNI)
   ENDIF
   IF (ALLOCATED(IWORKRRR)) DEALLOCATE(IWORKRRR)
   IF (ALLOCATED(IXVARX)) THEN
      DEALLOCATE(IXVARX)
      DEALLOCATE(IWRKNINI)
   ENDIF
   IF (ALLOCATED(IXI)) DEALLOCATE(IXI)
   IF (ALLOCATED(IYI)) DEALLOCATE(IYI)
      
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
   
   CHARACTER (LEN=8), ALLOCATABLE :: IPARMLAB(:), IPARMLB(:)
   CHARACTER*16 BLAB
   CHARACTER*4  BLANK(40)
   
   INTEGER  I2, IHI, ILO, IND, IND2, NCONN
   REAL (KIND=8)  CORR, SIDE
   
   REAL (KIND=8), ALLOCATABLE ::  ICDER2(:), ICPVAL(:), ICSE(:), &
      ICVAL(:), ICZVAL(:), &
      IINFSQ(:), IPFIX(:), IPVAR(:), ISEFIX(:), &
      ISEVAR(:), IVAL(:), IZFIX(:), IZVAR(:)
   
   IF(INITIALIZED .NE. 2) THEN
      CALL POST_ERROR("Cleanup() routine called out of sequence in " // PROGNAME)
      RETURN
   ENDIF
   
   WRITE(2,455)
   455 FORMAT(///,1x,'---------------------------------------------------------------',/, &
                  1x,'* Final Results - Maximum Marginal Likelihood (MML) Estimates *',/, &
                  1x,'---------------------------------------------------------------',/)
   WRITE(2,562)NEM,IT-NEM,IT,RLOGL,RLOGL-NPAR, &
   RLOGL-0.5D0*DBLE(NPAR)*DLOG(DBLE(NTOT)),0-2*RLOGL,0-2*(RLOGL-NPAR), &
   0-2*(RLOGL-0.5D0*DBLE(NPAR)*DLOG(DBLE(NTOT)))
   562 FORMAT(1X,'EM     Iterations =',I4,/,1X,'Fisher Iterations =',I4,/, &
         1X,'Total  Iterations =',I4,//, &
         1X,'Log Likelihood                 = ',F12.3,/, &
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

    open(100, file="mixreg.lik")
    write(100,*) rlogl, npar
    close(100)
    
   ! Find the z-statistics for all terms  IZFIX() and IZVAR()
   ! and their associated p-values        IPFIX() and IPVAR()
     
   ! save the estimates of the fixed effects in IVAL
   
   ! Start allocating back at the base, overwriting vars used
   ! during the iterations.  TG 
   ALLOCATE(IVAL(NF))
   ALLOCATE(ISEFIX(NF))
   ALLOCATE(IZFIX(NF))
   ALLOCATE(IPFIX(NF))
   
   ALLOCATE(ISEVAR(NV))
   ALLOCATE(IZVAR(NV))
   ALLOCATE(IPVAR(NV))

   IF (NCON .GT. 0) THEN
      ALLOCATE(ICVAL(NCON))
      ALLOCATE(ICSE(NCON))
      ALLOCATE(ICZVAL(NCON))
      ALLOCATE(ICPVAL(NCON))
      NCONN = (NCON * (NCON+1)) / 2
      ALLOCATE(ICDER2(NCONN))
   ENDIF

   IF (R .GT. 0 .AND. NOMU .EQ. 0) THEN
      DO H = 1,R
         IVAL(H)  = IMU(H)
      END DO
   ENDIF
   
   IF (P .GT. 0) THEN
      DO L = 1,P
         L1 =   L + R*(1-NOMU)
         IVAL(L1) = IALPHA(L)
      END DO
   ENDIF
   
   IF (R .GT. 0) THEN
      DO HR= 1,RR
         IZVAR(HR) = IVARCO(HR)
      END DO
   ENDIF
   
   IR  = RR + 1
   IZVAR(IR) = ERROR
   IF (NAUTO .GT. 1) THEN
      IR = RR + 1 
      DO IS = 1,S
         IR   = IR + 1
         IZVAR(IR) = IAUTO(IS)
      END DO
   ENDIF

   IC = 0
   DO I = 1,NF
      DO J = 1,I
         IC = IC + 1
         IF (I .EQ. J) THEN
            ISEFIX(I) = DSQRT(IINFFIX(IC))
            IZFIX(I)  = IVAL(I) / ISEFIX(I)
            IPFIX(I) = 2.0D0 *(1.0D0 - PHIFN(DABS(IZFIX(I)),0))
         ENDIF
      END DO
   END DO

   IC = 0
   DO I = 1,NV
      DO J = 1,I
         IC = IC + 1
         IF (I .EQ. J) THEN
            ISEVAR(I) = DSQRT(IINFVAR(IC))
            IZVAR(I)  = IZVAR(I) / ISEVAR(I)
            IPVAR(I) = 2.0D0 *(1.0D0 - PHIFN(DABS(IZVAR(I)),0))
         ENDIF
      END DO
   END DO

    804 FORMAT(A16,4(4x,F12.5))

   IF (R .GT. 0 .AND. NOMU .EQ. 0) THEN
      DO H = 1,R
         WRITE(2,804) IBLABEL(H),IMU(H),ISEFIX(H),IZFIX(H),IPFIX(H)
      END DO
   ENDIF
   
   IF (P .GT. 0) THEN
      DO L = 1,P
         L1 = L + R*(1-NOMU)
         WRITE(2,804)IALABEL(L),IALPHA(L),ISEFIX(L1),IZFIX(L1),IPFIX(L1)
      END DO
   ENDIF

   IF (R .GT. 0) THEN
      WRITE(2,581)
      581 FORMAT(/,'Random-effect variance & covariance term(s)')
      IF (NOCOV .EQ. 0) THEN
         IND = 1
         IND2= 0
         DO HR= 1,RR
            IF (HR .EQ. IND) THEN
               IND2 = IND2 + 1
               IND  = IND  + (IND2 + 1)
               Blab = IBLABEL(IND2)
               IWORKR(IND2) = 1.0D0 / DSQRT(IVARCO(HR))
               SIDE = 0.5D0
            ELSE
               Blab = 'covariance      '
               SIDE = 1.0D0
            ENDIF
            IPVAR(HR) = IPVAR(HR)*SIDE
            WRITE(2,804)Blab,IVARCO(HR),ISEVAR(HR),IZVAR(HR),IPVAR(HR)
            805 FORMAT(A16,2x,f12.5,3(4x,F12.5))
         END DO
      ELSEIF (NOCOV .EQ. 1) THEN 
         DO HR= 1,RR
            Blab = IBLABEL(HR) // '  '
            SIDE = 0.5D0
            IPVAR(HR) = IPVAR(HR)*SIDE
            WRITE(2,804)Blab,IVARCO(HR),ISEVAR(HR),IZVAR(HR),IPVAR(HR)
            815 FORMAT(A16,2x,f12.5,3(4x,F12.5))
         END DO
      ENDIF
   ENDIF

   WRITE(2,583)
   583 FORMAT(/,'Residual variance')
   IR = RR+ 1
   IPVAR(IR) = IPVAR(IR)*0.5D0
   WRITE(2,806)ERROR,ISEVAR(IR),IZVAR(IR),IPVAR(IR)
   806 FORMAT(16x,4(4x,F12.5))
   
   IF (NAUTO .GT. 1) THEN
      WRITE(2,584)
      584 FORMAT(/,'Autocorrelation term(s)')
      IR = RR + 1 
      DO IS = 1,S
         IR   = IR + 1
         WRITE(2,810)IAUTO(IS),ISEVAR(IR),IZVAR(IR),IPVAR(IR)
         810 FORMAT(16x,4(4x,F12.5))
      END DO
   ENDIF
   
   WRITE(2,587)
   587 FORMAT(//,'note: p-values are 2-tailed except for variances which are 1-tailed')

   IF (R .EQ. 1 .AND. ICCY .EQ. 2) THEN
      CORR = IVARCO(1) / (IVARCO(1) + ERROR)
      WRITE(2,811)ERROR,IVARCO(1),IVARCO(1),IVARCO(1),ERROR,CORR
      811 FORMAT(///,1X,'Calculation of the intracluster correlation',/, &
              1x,'-------------------------------------------',/, & 
              1x,'residual variance = ',F8.3,/, &
              1x,'cluster  variance = ',f8.3,//, &
              1x,'intracluster correlation = ',f8.3,' / (',f8.3,' + ',F8.3,')',' = ',F8.3)
   ELSEIF (R .GT. 1 .AND. NOCOV .EQ. 0) THEN
      DO I = 1,30
         BLANK(I) = '    '
      END DO
      ND = 0

      CALL MPDSD(IWORKR,IVARCO,IVARCO,R)
      CALL PRNT(2,IVARCO,R,R,1,IBLABEL,IBLABEL,ND,BLANK,1,78,5,2,2, &
         "Random-effect covariances expressed as correlations      ")
   ENDIF

   ! compute the contrasts of the parameter estimates
   ! uses IZFIX as a work vector

   IF (NCON .GE. 1) THEN
      CALL MPYTR(IVAL,ICON,ICVAL,NF,1,0,NCON)
      CALL GRMMT(ICON,IINFFIX,ICDER2, NF,NCON,1,IZFIX)
      IC = 0
      DO I = 1,NCON
         DO J = 1,I
            IC = IC + 1
            IF (I .EQ. J) THEN
               ICSE(I)   = DSQRT(ICDER2(IC))
               ICZVAL(I) = ICVAL(I) / ICSE(I)
               ICPVAL(I) = 2.0D0 *(1.0D0 -PHIFN(DABS(ICZVAL(I)),0))
            ENDIF
         END DO
      END DO
   ENDIF

   ! write out the parameter estimates and 
   ! the estimated variance-covariance estimates of these estimates

   ! write out the estimates of the fixed effects 

   IF (NF .GT. 1) THEN
      ALLOCATE(IPARMLAB(NF))

      DO I = 1,NF
         IF (I .LE. R .AND. NOMU .EQ. 0) THEN
            IPARMLAB(I) = IBLABEL(I)
         ELSEIF ((NOMU .EQ. 0 .AND. I .GT. R) .OR. NOMU .EQ. 1) THEN
            I2 = I - R*(1-NOMU)
            IPARMLAB(I) = IALABEL(I2)
         ENDIF
         ISEFIX(I) = 1.0D0 / ISEFIX(I)
         WRITE(3,704)IPARMLAB(I),IVAL(I)
         704 FORMAT(1X,a16,1x,F13.7)
      END DO
if(r .eq. 0) write(3,704) "Residual",error
      ALLOCATE(IINFSQ(NF*NF))
      CALL CHAMS(IINFFIX,IINFSQ,NF,1,0)

      DO I = 1,NF
         ILO = I                  
         IHI = ILO-1  + NF*NF            
         WRITE(4,705)(IINFSQ(I2),I2=ILO,IHI,NF)
         705 FORMAT(1X,6F17.7)
      END DO

if(r .eq. 0) write(4,705) iinfvar((rr+1)*(rr+2)/2)

      CALL MPDSD(ISEFIX,IINFFIX,IINFFIX,NF)
!      CALL PRNT(2,IINFFIX,NF,NF,1,IPARMLAB,IPARMLAB,ND,BLANK,1,78,5,2,2, &
!               "Correlation of the MML estimates of the fixed terms      ")
   ENDIF

   IF (NV .GT. 1) THEN
      ALLOCATE(IPARMLB(NV))

      DO I = 1,NV
         IF (I .LE. RR) THEN
             write(iparmlb(i),'(a6,i1)') 'VarCov',I
         ELSEIF (I .EQ. RR+1) THEN
            IPARMLB(I) = 'Residual'
         ELSEIF (I .GT. RR+1) THEN
            write(iparmlb(i),'(a7,i1)') 'AutoCor', I-RR-1
         ENDIF
         ISEVAR(I) = 1.0D0 / ISEVAR(I)
      END DO


      DO I = 1,NV
         IZVAR(I) = IZVAR(I) / ISEVAR(I)
         WRITE(3,706)IPARMLB(I),IZVAR(I)
         706 FORMAT(1X,a16,1x,F13.7)
      END DO

      IF(ALLOCATED(IINFSQ)) DEALLOCATE(IINFSQ)
      ALLOCATE(IINFSQ(NV*NV))
      CALL CHAMS(IINFVAR,IINFSQ,NV,1,0)

      DO I = 1,NV
         ILO = I                  
         IHI = ILO-1  + NV*NV            
         WRITE(4,707)(IINFSQ(I2),I2=ILO,IHI,NV)
         707 FORMAT(1X,6F13.7)
      END DO

      CALL MPDSD(ISEVAR,IINFVAR,IINFVAR,NV)
      CALL PRNT(2,IINFVAR,NV,NV,1,IPARMLB,IPARMLB,ND,BLANK,1,78,5,2,2, &
          "Correlation of the MML estimates of variance-related terms")
   ENDIF

   CLOSE(3)
   CLOSE(4)

   ! print out the requested contrasts of the parameter estimates

   IF (NCON .GE. 1) THEN

      WRITE(2,655)
      655 FORMAT(///,1x,'----------------------------------------------------------',/, &
              1x,'* Transforms of parameter estimates of the fixed effects *',/, &
              1x,'----------------------------------------------------------')

      CALL PRNT(2,ICON,NF,NCON,0,IPARMLAB,IXLAB,ND,BLANK,1,78,5,2,1, &
         "Transpose of the Transform Matrix (parameters by transforms)")

      WRITE(2,757)
      757 FORMAT(/,1X,'Transform',5x,'    Estimate',5X,'Stand. Error',5x, &
          '          Z',5X,'     p-value',/,1X,'---------',5x,'------------',5X, &
          '------------',5X,'------------',5X,'------------')

      DO L = 1,NCON
         WRITE(2,840)L,ICVAL(L),ICSE(L),ICZVAL(L),ICPVAL(L)
         840 FORMAT(1X,I8,4(5x,F12.5))
         ICSE(L) = 1.0D0 / ICSE(L)
      END DO

      WRITE(2,589)
      589 FORMAT(/,1x,'note: p-values are 2-tailed')

      CALL MPDSD(ICSE,ICDER2,ICDER2,NCON)
      CALL PRNT(2,ICDER2,NCON,NCON,1,IXLAB,IXLAB,ND,BLANK,1,78,5,1,1, &
            "Correlation of the MML Transformed Estimates              ")

   ENDIF

   CLOSE(2)
   IF (IRES .EQ. 1) CLOSE(5)
   
   ! Deallocate local arrays
      
   IF(ALLOCATED(IPARMLAB)) DEALLOCATE(IPARMLAB)
   IF(ALLOCATED(IPARMLB)) DEALLOCATE(IPARMLB)
   DEALLOCATE(IVAL, ISEFIX, IZFIX, IPFIX, ISEVAR, IZVAR, IPVAR)
   
   IF(ALLOCATED(IINFSQ)) DEALLOCATE(IINFSQ)
   IF(ALLOCATED(ICVAL)) DEALLOCATE( ICVAL, ICSE, ICZVAL, ICPVAL, ICDER2)
      
   ! Deallocate shared arrays in preparation for starting again.
   ! This is needed when we are operating as a DLL, where the 
   ! program never really 'ends'. 
   
!   DEALLOCATE(ALLDAT, IDNI, IXLAB, IALABEL, IBLABEL, IXIND, IWORKP, IINFFIX)
   
   DEALLOCATE( IALPHA, ICON, IINFVAR, IVARCO, IGRAFIX, ICORFIX)
   
   IF(ALLOCATED(IALPHA0)) DEALLOCATE(IALPHA0)
   
   IF(ALLOCATED(IAUERIN)) DEALLOCATE(IAUERIN, IAUTOEM, &
      IAUTOEM2, IAUTVIN, IAUV1IN, IAUTLIN)
   
   IF(ALLOCATED(IAUTO)) DEALLOCATE(TIME,IAUTO)
   IF(ALLOCATED(IAUTO0)) DEALLOCATE(IAUTO0) 
   
   IF (ALLOCATED(IMUEMGR)) DEALLOCATE(IMUEMGR)
   IF (ALLOCATED(IMUIN)) DEALLOCATE(IMUIN)
   IF (ALLOCATED(IMU1IN)) DEALLOCATE(IMU1IN)
   
!   DEALLOCATE(IVARPI, IVARD, IWORKR, IWORKRR, IVARCEM , IVARL, &
!      IERV1IN , IVARCON , IVARCOI , IVARCSQ , IVARLSQ, &
!      IVARGR , IWORKR2 , IVARCIN , IVARCK , IVARIN , IW3, &
!     ITRANN , ITRANJ , IVARLD , ITRANH , IWORKR3 , ITRANG, IW3G , IWORKPP )
   
   IF (ALLOCATED(IERVAIN)) DEALLOCATE( IERVAIN, IWKEC, IGRAVAR,  &
            ICORVAR, IERRINF, IWKER2) 
      
   IF (ALLOCATED(IMU)) DEALLOCATE(IMU)
   IF (ALLOCATED(IMU0)) DEALLOCATE(IMU0)
      
   IF (ALLOCATED(IVARCO0)) DEALLOCATE(IVARCO0)
   IF (ALLOCATED(IAL1EM)) DEALLOCATE(IAL1EM)
   IF (ALLOCATED(IAL2EM)) DEALLOCATE(IAL2EM)
   IF (ALLOCATED(IALM1IN)) DEALLOCATE(IALM1IN)
   IF (ALLOCATED(IALMUIN)) DEALLOCATE(IALMUIN)
   
   INITIALIZED = 0

END SUBROUTINE CLEANUP

      
      
! ************************************************
!         **********************                         
!         *  FUNCTION NumericInputCount(InString)
!         **********************                         
!                                                                       
!     Counts the number of contiguous nonblank character sequences
!     in the string passed in, to determine how many numeric values
!     are present on the given line.  No attempt to validate the 
!     numeric character of the tokens is made.  Used to distinguish
!     between versions of the DEF file, which differ in the number
!     of parameters supplied on some lines.
!                                                                       
! InString ... The string with the tokens to count
!                                                                       
! ************************************************
INTEGER FUNCTION NumericInputCount( InString )
   CHARACTER (len=*),INTENT(IN):: InString
   INTEGER TokenCount,InStringPos
   LOGICAL InToken
   
   InToken = .false.
   TokenCount = 0
   
   DO InStringPos = 1, LEN(InString)
      IF (InString(InStringPos:InStringPos) == ' ') THEN
         InToken = .false.
      ELSE
         IF (.NOT. InToken) THEN
            InToken = .true.
            TokenCount = TokenCount + 1
         END IF
      END IF
   END DO
   
   NumericInputCount = TokenCount
   
END FUNCTION

END MODULE RRM_MODULE

! ************************************************
! The DLL wrapper and EXE stub routines are identical for all 
! four modules, but need to be included in each file rather than
! compiled as as module.  So the INCLUDE directive is used to 
! avoid duplication - a single file is used in all modules.
! ************************************************
INCLUDE "DLLSTUB.F90"
