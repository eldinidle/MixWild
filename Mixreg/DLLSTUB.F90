! Note that the module that INCLUDE's this must define a character
! constant called PROGNAME that contains the program name without
! the 'b' on the end, like MIXOR,MIXNO,MIXREG or MIXPREG.
PROGRAM RRM_EXE

   USE RRM_MODULE
   USE MIXLIB, ONLY: SET_FILENAMES,POST_ERROR
   CHARACTER*256 CMDLINE
   CHARACTER*40 FILENAMES(4)
   
   ! CALL GETCL(CMDLINE)
   CMDLINE = " " 
   ! Check to see if any of the default names for the 
   ! input or output files have been overridden by the 
   ! caller.
   ! DEFFILE,DATFILE,OUTFILE,BASENAME
   CALL SET_FILENAMES(FILENAMES, CMDLINE, PROGNAME )
   
   CALL INIT(FILENAMES)  ! supply command line params to INIT to override
                       ! filenames for DEF, DAT, etc.
   DO WHILE(MainLoop() .NE. 0) 
   ENDDO
   CALL CLEANUP()
END PROGRAM RRM_EXE
   
SUBROUTINE CALL_INIT(CMDLINE)
   USE RRM_MODULE
   USE MIXLIB, ONLY: SET_FILENAMES,POST_ERROR
   CHARACTER (LEN=*), INTENT(IN OUT):: CMDLINE
   CHARACTER (LEN=40) :: FILENAMES(4)=" "
!   DLL_EXPORT CALL_INIT
   
   ! Check to see if any of the default names for the 
   ! input or output files have been overridden by the 
   ! caller.
   ! DEFFILE,DATFILE,OUTFILE,BASENAME
   CALL SET_FILENAMES(FILENAMES, CMDLINE, PROGNAME )
   
   OPEN(6, FILE= TRIM(ADJUSTL(FILENAMES(4))) // '.TXT')
   CALL INIT(FILENAMES)
END SUBROUTINE
   
INTEGER FUNCTION CALL_MainLoop()
   USE RRM_MODULE
!   DLL_EXPORT CALL_MainLoop
   CALL_MainLoop = MainLoop()
END FUNCTION
   
SUBROUTINE CALL_CLEANUP()
   USE RRM_MODULE
!   DLL_EXPORT CALL_CLEANUP
   CALL CLEANUP()
   CLOSE(6)  ! Make sure output file is updated
END SUBROUTINE

INTEGER FUNCTION CALL_RETRIEVE_ERROR(ErrorString)
   USE MIXLIB, ONLY: RETRIEVE_ERROR
!   DLL_EXPORT CALL_RETRIEVE_ERROR
   CHARACTER (len=*), INTENT(IN OUT):: ErrorString
   CALL_RETRIEVE_ERROR = RETRIEVE_ERROR(ErrorString)
   RETURN
END FUNCTION CALL_RETRIEVE_ERROR
