#include "PACKAGES_CONFIG.h"
#include "CPP_OPTIONS.h"

CBOP
C     !ROUTINE: CPL_ISSM
C     !INTERFACE:
      SUBROUTINE CPL_ISSM( myTime, myIter, myThid )

C     !DESCRIPTION: \bv
C     *==================================================================
C     | SUBROUTINE cpl_issm
C     | o Couple MITgcm ocean model with ISSM ice sheet model
C     *==================================================================
C     \ev

C     !USES:
      IMPLICIT NONE
C     == Global variables ==
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "DYNVARS.h"
#include "GRID.h"
#include "FFIELDS.h"
#include "SHELFICE_OPTIONS.h"
#include "SHELFICE.h"
#ifdef ALLOW_EXF
# include "EXF_OPTIONS.h"
# include "EXF_FIELDS.h"
#endif

      LOGICAL  DIFFERENT_MULTIPLE
      EXTERNAL DIFFERENT_MULTIPLE

C     !LOCAL VARIABLES:
C     mytime - time counter for this thread (seconds)
C     myiter - iteration counter for this thread
C     mythid - thread number for this instance of the routine.
      _RL     mytime
      INTEGER myiter, mythid 
CEOP

#ifdef ALLOW_CPL_ISSM
#include "EESUPPORT.h"
      COMMON /CPL_MPI_ID/ mpiMyWid, toissmcomm
      INTEGER mpiMyWid, toissmcomm, mpiRC
      INTEGER mpistatus(MPI_STATUS_SIZE)
      INTEGER i, j, bi, bj, buffsize
      COMMON /CPL_ISSM_TIME/ CouplingTime
      _R8 CouplingTime, IceModelTime
      _R8 xfer_array(Nx,Ny)
      _R8 local(1:sNx,1:sNy,nSx,nSy)
      CHARACTER*(MAX_LEN_MBUF) suff

# ifdef ALLOW_DEBUG
      IF (debugMode) CALL DEBUG_ENTER('CPL_ISSM',myThid)
# endif

      IF ( ( myTime .NE. startTime )
     &     .AND. ( myTime .NE. endtime )
     &     .AND. ( MOD(myTime,CouplingTime) .LT. deltaT/2. )
     &   ) THEN

C Recurring steps C3, C4, C5, C1, and C2:

C   C3b. ISSM sends ice model time IceTimeTag
         IF( myProcId .EQ. 0 ) THEN
            _BEGIN_MASTER( myThid )
            call MPI_Recv(IceModelTime,1,MPI_DOUBLE,0,10001001,
     &           toissmcomm,mpistatus,mpiRC)
#ifdef ALLOW_DEBUG
            IF ( debugLevel .GE. debLevD ) THEN
               print*, 'Ocean received IceModelTime: ', IceModelTime
            ENDIF
#endif
            _END_MASTER( myThid )
         ENDIF

C   C3c. MITgcm sends ocean model time OceanTimeTag
         IF( myProcId .EQ. 0 ) THEN
            _BEGIN_MASTER( myThid )
            call MPI_Send(myTime,1,MPI_DOUBLE,0,10001002,
     &           toissmcomm,mpistatus)
            _END_MASTER( myThid )
         ENDIF

C  C4. MITgcm sends
C      (N-1)*CouplingTime <= OceanModelTime < N*CouplingTime
C      time-mean melt rate to ISSM
         DO bj=1,nSy
            DO bi=1,nSx
               DO j=1,sNy
                  DO i=1,sNx
                     local(i,j,bi,bj)=shelficeFreshWaterFlux(i,j,bi,bj)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         CALL BAR2( myThid ) 
         CALL GATHER_2D_R8( xfer_array, local, Nx, Ny,
     &        .FALSE., .FALSE., myThid )
         IF( myProcId .EQ. 0 ) THEN
            _BEGIN_MASTER( myThid )
            buffsize = Nx*Ny
            CALL MPI_SEND(xfer_array,buffsize,MPI_DOUBLE_PRECISION,
     &           0,10001007,toissmcomm,mpistatus)
            _END_MASTER( myThid )
         ENDIF
         CALL BAR2( myThid )
#ifdef ALLOW_DEBUG
         IF ( debugLevel .GE. debLevD ) THEN
            print*,'Done Sending shelficeFreshWaterFlux array.'
         ENDIF
#endif

C   C2a. ISSM sends (dynamic) IceMass to MITgcm
         IF( myProcId .EQ. 0 ) THEN
            _BEGIN_MASTER( myThid )         
            call MPI_Recv(xfer_array,buffsize,MPI_DOUBLE_PRECISION,
     &           0,10001008,toissmcomm,mpistatus,mpiRC)
            _END_MASTER( myThid )
         ENDIF
         CALL BAR2( myThid ) 
         CALL SCATTER_2D_R8( xfer_array, local, Nx, Ny,
     &        .FALSE., .FALSE., myThid )
         DO bj = myByLo(myThid), myByHi(myThid)
            DO bi = myBxLo(myThid), myBxHi(myThid)
               DO j=1,sNy
                  DO i=1,sNx
                     shelficeMassInit(i,j,bi,bj ) =
     &                    shelficeMass(i,j,bi,bj)
                    IF( local(i,j,bi,bj).LT.9998 ) THEN
                        shelficeMass(i,j,bi,bj) = local(i,j,bi,bj)
                     ELSE
                        shelficeMass(i,j,bi,bj) = 0. _d 0
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO
         ENDDO

C   C2b. MITgcm computes DIceMassDT (dynamic ice mass change + correction)
         DO bj = myByLo(myThid), myByHi(myThid)
            DO bi = myBxLo(myThid), myBxHi(myThid)
               DO j=1-OLy,sNy+OLy
                  DO i=1-OLx,sNx+OLy
                     shelfIceMassDynTendency(i,j,bi,bj) =
     &                    ( shelficeMass(i,j,bi,bj) -
     &                    shelficeMassInit(i,j,bi,bj ) ) / CouplingTime
                     shelficeMassInit(i,j,bi,bj ) =
     &                    shelficeMass(i,j,bi,bj)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
C End of recurring steps C3, C4, C5, C1, and C2

      ELSEIF ( myTime .EQ. endtime ) THEN
C     Termination steps C3, C4, and C5

C   C3b. ISSM sends ice model time IceTimeTag
         IF( myProcId .EQ. 0 ) THEN
            _BEGIN_MASTER( myThid )
            call MPI_Recv(IceModelTime,1,MPI_DOUBLE,0,10001001,
     &           toissmcomm,mpistatus,mpiRC)
#ifdef ALLOW_DEBUG
            IF ( debugLevel .GE. debLevD ) THEN
               print*, 'Ocean received IceModelTime: ', IceModelTime
            ENDIF
#endif
            _END_MASTER( myThid )
         ENDIF

C   C3c. MITgcm sends ocean model time OceanTimeTag
         IF( myProcId .EQ. 0 ) THEN
            _BEGIN_MASTER( myThid )
            call MPI_Send(myTime,1,MPI_DOUBLE,0,10001002,
     &           toissmcomm,mpistatus)
            _END_MASTER( myThid )
         ENDIF

C  C4. MITgcm sends
C      (N-1)*CouplingTime <= OceanModelTime < N*CouplingTime
C      time-mean melt rate to ISSM
         DO bj=1,nSy
            DO bi=1,nSx
               DO j=1,sNy
                  DO i=1,sNx
                     local(i,j,bi,bj)=shelficeFreshWaterFlux(i,j,bi,bj)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         CALL BAR2( myThid ) 
         CALL GATHER_2D_R8( xfer_array, local, Nx, Ny,
     &        .FALSE., .FALSE., myThid )
         IF( myProcId .EQ. 0 ) THEN
            _BEGIN_MASTER( myThid )
            buffsize = Nx*Ny
            CALL MPI_SEND(xfer_array,buffsize,MPI_DOUBLE_PRECISION,
     &           0,10001007,toissmcomm,mpistatus)
            _END_MASTER( myThid )
         ENDIF
         CALL BAR2( myThid )
#ifdef ALLOW_DEBUG
         IF ( debugLevel .GE. debLevD ) THEN
            print*,'Done Sending shelficeFreshWaterFlux array.'
         ENDIF
#endif
      ENDIF
C     End of termination steps C3, C4, and C5

#endif /* ALLOW_CPL_ISSM */

      RETURN
      END
