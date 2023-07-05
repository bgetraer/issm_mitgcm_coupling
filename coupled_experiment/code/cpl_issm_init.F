#include "PACKAGES_CONFIG.h"
#include "CPP_OPTIONS.h"

CBOP
C     !ROUTINE: CPL_ISSM_INIT
C     !INTERFACE:
      SUBROUTINE CPL_ISSM_INIT( myThid )

C     !DESCRIPTION: \bv
C     *==================================================================
C     | SUBROUTINE cpl_issm_init
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
C     mythid - thread number for this instance of the routine.
      INTEGER mythid 
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
      IF (debugMode) CALL DEBUG_ENTER('CPL_ISSM_INIT',myThid)
# endif

C Initialization steps I0, I1, and I2:

# ifdef ALLOW_DEBUG
      IF (debugMode) CALL DEBUG_MSG('Initialization',myThid)
# endif

C   I0a. ISSM sends CouplingTime, the interval at which we couple
      IF( myProcId .EQ. 0 ) THEN
         _BEGIN_MASTER( myThid )
         call MPI_Recv(CouplingTime,1,MPI_DOUBLE,0,10001000,
     &        toissmcomm,mpistatus,mpiRC)
         _END_MASTER( myThid )
      ENDIF
#ifdef ALLOW_DEBUG
      IF (debugMode)
     &     CALL DEBUG_MSG('MIT waiting for CouplingTime',myThid)
#endif
      _BEGIN_MASTER( myThid )
      CALL MPI_BCAST(CouplingTime,1,MPI_DOUBLE,0,
     &     MPI_COMM_MODEL,mpiRC)
      _END_MASTER( myThid )
#ifdef ALLOW_DEBUG
      IF ( debugLevel .GE. debLevD ) THEN
         print*, 'Ocean received CouplingTime: ', CouplingTime
      ENDIF
#endif

C   I0b. MITgcm sends grid size (NX and NY)
      IF( myProcId .EQ. 0 ) THEN
         _BEGIN_MASTER( myThid )
         call MPI_Send(Nx,1,MPI_INT,0,10001003,
     &        toissmcomm,mpistatus)
         call MPI_Send(Ny,1,MPI_INT,0,10001004,
     &        toissmcomm,mpistatus)
         _END_MASTER( myThid )
      ENDIF
#ifdef ALLOW_DEBUG
      IF ( debugLevel .GE. debLevD ) THEN
         print*, 'grid size NX: ', NX
         print*, 'grid size NY: ', NY
      ENDIF
#endif

C   I0c. MITgcm sends grid coordinates of center of cells
C       (longitude -180 <= XC < 180 and latitude YC)
C     Send longitude East of center of cell
      DO bj=1,nSy
         DO bi=1,nSx
            DO j=1,sNy
               DO i=1,sNx
                  local(i,j,bi,bj) = xC(i,j,bi,bj)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      CALL BAR2( myThid ) 
      CALL GATHER_2D_R8( xfer_array, local, Nx, Ny,
     &     .FALSE., .FALSE., myThid )
      IF( myProcId .EQ. 0 ) THEN
         _BEGIN_MASTER( myThid )
         buffsize = Nx*Ny
         CALL MPI_SEND(xfer_array,buffsize,MPI_DOUBLE_PRECISION,
     &        0,10001005,toissmcomm,mpistatus)
         _END_MASTER( myThid )
      ENDIF
      CALL BAR2( myThid )
#ifdef ALLOW_DEBUG
      IF (debugMode) CALL DEBUG_STATS_RL(1,XC,'XC',myThid)
#endif         
C     Send latitude North of center of cell
      DO bj=1,nSy
         DO bi=1,nSx
            DO j=1,sNy
               DO i=1,sNx
                  local(i,j,bi,bj) = yC(i,j,bi,bj)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      CALL BAR2( myThid ) 
      CALL GATHER_2D_R8( xfer_array, local, Nx, Ny,
     &     .FALSE., .FALSE., myThid )
      IF( myProcId .EQ. 0 ) THEN
         _BEGIN_MASTER( myThid )
         buffsize = Nx*Ny
         CALL MPI_SEND(xfer_array,buffsize,MPI_DOUBLE_PRECISION,
     &        0,10001006,toissmcomm,mpistatus)
         _END_MASTER( myThid )
      ENDIF
      CALL BAR2( myThid )
#ifdef ALLOW_DEBUG
      IF (debugMode) CALL DEBUG_STATS_RL(1,YC,'YC',myThid)
#endif         

C   I0d. ISSM sends IceModelTime=0 IceMass to MITgcm
      IF( myProcId .EQ. 0 ) THEN
         _BEGIN_MASTER( myThid )         
         call MPI_Recv(xfer_array,buffsize,MPI_DOUBLE_PRECISION,
     &        0,10001008,toissmcomm,mpistatus,mpiRC)
         _END_MASTER( myThid )
      ENDIF
      CALL BAR2( myThid ) 
      CALL SCATTER_2D_R8( xfer_array, local, Nx, Ny,
     &     .FALSE., .FALSE., myThid )
      DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
            DO j=1,sNy
               DO i=1,sNx
                  IF( local(i,j,bi,bj).GT.0 ) THEN
                     shelficeMassInit(i,j,bi,bj) = local(i,j,bi,bj)
                     R_shelfIce(i,j,bi,bj) =
     &                    - local(i,j,bi,bj) / rhoShelfIce
                  ELSE
                     shelficeMassInit(i,j,bi,bj) = 0. _d 0
                     R_shelfIce(i,j,bi,bj) = 0. _d 0
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDDO

C   I2a. ISSM sends (dynamic) IceMass to MITgcm
      IF( myProcId .EQ. 0 ) THEN
         _BEGIN_MASTER( myThid )         
         call MPI_Recv(xfer_array,buffsize,MPI_DOUBLE_PRECISION,
     &        0,10001008,toissmcomm,mpistatus,mpiRC)
         _END_MASTER( myThid )
      ENDIF
      CALL BAR2( myThid ) 
      CALL SCATTER_2D_R8( xfer_array, local, Nx, Ny,
     &     .FALSE., .FALSE., myThid )
      DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
            DO j=1,sNy
               DO i=1,sNx
                  IF( local(i,j,bi,bj).LT.9998 ) THEN
                     shelficeMass(i,j,bi,bj) = local(i,j,bi,bj)
                  ELSE
                     shelficeMass(i,j,bi,bj) = 0. _d 0
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDDO

C   I2b. MITgcm computes DIceMassDT (dynamic ice mass change + correction)
      DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
            DO j=1-OLy,sNy+OLy
               DO i=1-OLx,sNx+OLy
                  shelfIceMassDynTendency(i,j,bi,bj) =
     &                 ( shelficeMass(i,j,bi,bj) -
     &                 shelficeMassInit(i,j,bi,bj ) ) / CouplingTime
               ENDDO
            ENDDO
         ENDDO
      ENDDO

C End initialization steps I0, I1, and I2.

#endif /* ALLOW_CPL_ISSM */

      RETURN
      END
