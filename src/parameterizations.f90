! subroutine para3d(X,n,m,l,ndim,U,V,W,uu,vv,ww)
!   !***DESCRIPTION
!   !
!   !     Written by Gaetan Kenway
!   !
!   !     Abstract para3d calculates the parametric locations for a 3d block
!   !
!   !     Description of Arguments
!   !     Input
!   !     X       - Real,size(n,m,l,ndim): Coordiantes
!   !     Output
!   !     U,V,W   - Real,size(n,m,l): The u,v,w parametric positions
!   !     u       - Real, size(n): The averaged u parameters
!   !     v       - Real, size(m): The averaged v parameters
!   !     w       - Real, size(l): The averaged w parameters

!   implicit none
  
!   ! Input
!   double precision , intent(in)   :: X(n,m,l,ndim)
!   integer          , intent(in)   :: n,m,l,ndim

!   ! Output
!   double precision , intent(out)  :: U(n,m,l),V(n,m,l),W(n,m,l)
!   double precision , intent(out)  :: u(n),v(m),w(l)

!   ! Working 
!   integer                         :: i,j,k


!   ! Zero all entries
!   u(:) = 0.0
!   v(:) = 0.0
!   w(:) = 0.0
!   U(:,:,:) = 0.0
!   V(:,:,:) = 0.0
!   W(:,:,:) = 0.0



   SUBROUTINE PARAM3DM (IL, JL, KL, XYZ, S)

!     ******************************************************************
!     *   PARAM3DM parameterizes the volume of one block of a multi-   *
!     *   block grid structure by setting up the normalized arc-length *
!     *   increments in all three index directions.                    *
!     *                                                                *
!     *   11/29/95  D.Saunders  Adaptation of PARAMXYZ for specialized *
!     *                         WARP-BLK used by FLO107-MB.            *
!     *   06/19/96      "       Allow for degenerate edges.            *
!     *                                                                *
!     *   David Saunders/James Reuther, NASA Ames Research Center, CA. *
!     ******************************************************************

      IMPLICIT REAL*8 (A-H,O-Z) ! Take out when all compilers have a switch

!     Arguments.

      INTEGER   IL, JL, KL                  ! I Grid array dimensions.
      double precision XYZ(3,0:IL+1,0:JL+1,0:KL+1) ! I Grid coordinates
      DIMENSION S(3,0:IL+1,0:JL+1,0:KL+1)   ! O Normalized arc-lengths:
                                            !   S(1,1,J,K) = 0.,
                                            !   S(2,I,1,K) = 0.,
                                            !   S(3,I,J,1) = 0.,
                                            !   S(1,IL,J,K) = 1.,etc.
!     Local constants.

!     REAL      ONE, ZERO
      PARAMETER (ONE = 1.E+0, ZERO = 0.E+0)

!    Local variables.

      INTEGER   I, J, K

!     Local functions.

!     REAL      DELI, DELJ, DELK

      DELI(I,J,K) = SQRT ((XYZ(1,I,J,K) - XYZ(1,I-1,J,K)) ** 2 + &
           (XYZ(2,I,J,K) - XYZ(2,I-1,J,K)) ** 2 + &
           (XYZ(3,I,J,K) - XYZ(3,I-1,J,K)) ** 2)

      DELJ(I,J,K) = SQRT ((XYZ(1,I,J,K) - XYZ(1,I,J-1,K)) ** 2 + &
           (XYZ(2,I,J,K) - XYZ(2,I,J-1,K)) ** 2 + &
           (XYZ(3,I,J,K) - XYZ(3,I,J-1,K)) ** 2)

      DELK(I,J,K) = SQRT ((XYZ(1,I,J,K) - XYZ(1,I,J,K-1)) ** 2 + &
           (XYZ(2,I,J,K) - XYZ(2,I,J,K-1)) ** 2 + &
           (XYZ(3,I,J,K) - XYZ(3,I,J,K-1)) ** 2)

!     Execution.
!     ----------

!     Zero the three low-end faces (or edges if one plane is specified).

      DO K = 1, KL
         DO J = 1, JL
            S(1,1,J,K) = ZERO
         END DO

         DO I = 1, IL
            S(2,I,1,K) = ZERO
         END DO
      END DO

      DO J = 1, JL
         DO I = 1, IL
            S(3,I,J,1) = ZERO
         END DO
      END DO

!     Set up the low-end edge lines because they are missed by the
!     following loops over most of the low-end faces:

      DO I = 2, IL
         S(1,I,1,1) = S(1,I-1,1,1) + DELI(I,1,1)
      END DO

      DO J = 2, JL
         S(2,1,J,1) = S(2,1,J-1,1) + DELJ(1,J,1)
      END DO

      DO K = 2, KL
         S(3,1,1,K) = S(3,1,1,K-1) + DELK(1,1,K)
      END DO

!     Set up the rest of the low-end face lines because they are
!     missed by the the main loop over most of the volume.

      DO K = 2, KL
         DO J = 2, JL
            S(2,1,J,K) = S(2,1,J-1,K) + DELJ(1,J,K)
            S(3,1,J,K) = S(3,1,J,K-1) + DELK(1,J,K)
         END DO
         DO I = 2, IL
            S(1,I,1,K) = S(1,I-1,1,K) + DELI(I,1,K)
            S(3,I,1,K) = S(3,I,1,K-1) + DELK(I,1,K)
         END DO
      END DO

      DO J = 2, JL
         DO I = 2, IL
            S(1,I,J,1) = S(1,I-1,J,1) + DELI(I,J,1)
            S(2,I,J,1) = S(2,I,J-1,1) + DELJ(I,J,1)
         END DO
      END DO

!     Traverse the block just once for all lines except those within
!     the low-end faces.

      DO K = 2, KL
         DO J = 2, JL
            DO I = 2, IL
               S(1,I,J,K) = S(1,I-1,J,K) + DELI(I,J,K)
               S(2,I,J,K) = S(2,I,J-1,K) + DELJ(I,J,K)
               S(3,I,J,K) = S(3,I,J,K-1) + DELK(I,J,K)
            END DO
         END DO
      END DO

!     Normalizing requires another pass through the volume.
!     Handle lines of zero length first by inserting uniform
!     distributions.  Then the standard normalization can be
!     applied safely everywhere.

      DO K = 1, KL

!        Zero-length lines in the I direction?

         DO J = 1, JL
            IF (S(1,IL,J,K) .EQ. ZERO) THEN
               DO I = 2, IL
                  S(1,I,J,K) = I - 1
               END DO
            END IF
         END DO

!        Zero-length lines in the J direction?

         DO I = 1, IL
            IF (S(2,I,JL,K) .EQ. ZERO) THEN
               DO J = 2, JL
                  S(2,I,J,K) = J - 1
               END DO
            END IF
         END DO
      END DO

!     Zero-length lines in the K direction?

      DO J = 1, JL
         DO I = 1, IL
            IF (S(3,I,J,KL) .EQ. ZERO) THEN
               DO K = 2, KL
                  S(3,I,J,K) = K - 1
               END DO
            END IF
         END DO
      END DO

!     Normalize:

      DO K = 1, KL
         DO J = 1, JL
            DO I = 1, IL
               S(1,I,J,K) = S(1,I,J,K) / S(1,IL,J,K)
               S(2,I,J,K) = S(2,I,J,K) / S(2,I,JL,K)
               S(3,I,J,K) = S(3,I,J,K) / S(3,I,J,KL)
            END DO
         END DO
      END DO

!     Finally, precise 1s for the three high-end faces:

      DO K = 1, KL
         DO J = 1, JL
            S(1,IL,J,K) = ONE
         END DO

         DO I = 1, IL
            S(2,I,JL,K) = ONE
         END DO
      END DO

      DO J = 1, JL
         DO I = 1, IL
            S(3,I,J,KL) = ONE
         END DO
      END DO
      
    END SUBROUTINE PARAM3DM
