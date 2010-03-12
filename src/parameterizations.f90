subroutine para3d(X,n,m,l,ndim,S,u,v,w)
  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract para3d calculates the parametric locations for a 3d block
  !
  !     Description of Arguments
  !     Input
  !     X       - Real,size(n,m,l,ndim): Coordiantes
  !     Output
  !     S       - Real,size(n,m,l,3): The u,v,w parametric positions
  !     u       - Real, size(n): The averaged u parameters
  !     v       - Real, size(m): The averaged v parameters
  !     w       - Real, size(l): The averaged w parameters

  implicit none

  ! Input
  integer          , intent(in)   :: n,m,l,ndim
  double precision , intent(in)   :: X(n,m,l,ndim)
  

  ! Output
  double precision , intent(out)  :: S(n,m,l,ndim)
  double precision , intent(out)  :: u(n),v(m),w(l)

  ! Working 
  integer                         :: i,j,k

  double precision DELI,DELJ,DELK

  DELI(I,J,K) = SQRT ((X(I,J,K,1) - X(I-1,J,K,1)) ** 2 + &
       (X(I,J,K,2) - X(I-1,J,K,2)) ** 2 + &
       (X(I,J,K,3) - X(I-1,J,K,3)) ** 2)

  DELJ(I,J,K) = SQRT ((X(I,J,K,1) - X(I,J-1,K,1)) ** 2 + &
       (X(I,J,K,2) - X(I,J-1,K,2)) ** 2 + &
       (X(I,J,K,3) - X(I,J-1,K,3)) ** 2)

  DELK(I,J,K) = SQRT ((X(I,J,K,1) - X(I,J,K-1,1)) ** 2 + &
       (X(I,J,K,2) - X(I,J,K-1,2)) ** 2 + &
       (X(I,J,K,3) - X(I,J,K-1,3)) ** 2)


  u(:) = 0.0
  v(:) = 0.0
  w(:) = 0.0

  !     Zero the three low-end faces (or edges if one plane is specified).
  
  DO K = 1, l
     DO J = 1, m
        S(1,J,K,1) = 0.0
     END DO
     
     DO I = 1, n
        S(I,1,K,2) = 0.0
     END DO
  END DO
  
  DO J = 1, m
     DO I = 1, n
        S(I,J,1,3) = 0.0
     END DO
  END DO
  
  !     Set up the low-end edge lines because they are missed by the
  !     following loops over most of the low-end faces:

  DO I = 2, n
     S(I,1,1,1) = S(I-1,1,1,1) + DELI(I,1,1)
  END DO

  DO J = 2, m
     S(1,J,1,2) = S(1,J-1,1,2) + DELJ(1,J,1)
  END DO

  DO K = 2, l
     S(1,1,K,3) = S(1,1,K-1,3) + DELK(1,1,K)
  END DO

  !     Set up the rest of the low-end face lines because they are
  !     missed by the the main loop over most of the volume.

  DO K = 2, l
     DO J = 2, m
        S(1,J,K,2) = S(1,J-1,K,2) + DELJ(1,J,K)
        S(1,J,K,3) = S(1,J,K-1,3) + DELK(1,J,K)
     END DO
     DO I = 2, n
        S(I,1,K,1) = S(I-1,1,K,1) + DELI(I,1,K)
        S(I,1,K,3) = S(I,1,K-1,3) + DELK(I,1,K)
     END DO

  END DO

  DO J = 2, m
     DO I = 2, n
        S(I,J,1,1) = S(I-1,J,1,1) + DELI(I,J,1)
        S(I,J,1,2) = S(I,J-1,1,2) + DELJ(I,J,1)
     END DO
  END DO

  !     Traverse the block just once for all lines except those within
  !     the low-end faces.

  DO K = 2, l
     DO J = 2, m
        DO I = 2, n
           S(I,J,K,1) = S(I-1,J,K,1) + DELI(I,J,K)
           S(I,J,K,2) = S(I,J-1,K,2) + DELJ(I,J,K)
           S(I,J,K,3) = S(I,J,K-1,3) + DELK(I,J,K)
        END DO
     END DO
  END DO

  !     Normalizing requires another pass through the volume.
  !     Handle lines of zero length first by inserting uniform
  !     distributions.  Then the standard normalization can be
  !     applied safely everywhere.

  DO K = 1, l

     !        Zero-length lines in the I direction?

     DO J = 1, m
        IF (S(n,J,K,1) == 0.0) THEN
           DO I = 2, n
              S(I,J,K,1) = I - 1
           END DO
        END IF
     END DO

     !        Zero-length lines in the J direction?

     DO I = 1, n
        IF (S(I,m,K,2) == 0.0) THEN
           DO J = 2, m
              S(I,J,K,2) = J - 1
           END DO
        END IF
     END DO
  END DO

  !     Zero-length lines in the K direction?

  DO J = 1, m
     DO I = 1, n
        IF (S(I,J,l,3) == 0.0) THEN
           DO K = 2, l
              S(I,J,K,3) = K - 1
           END DO
        END IF
     END DO
  END DO

  !     Normalize:

  DO K = 1, l
     DO J = 1, m
        DO I = 1, n
           S(I,J,K,1) = S(I,J,K,1) / S(n,J,K,1)
           S(I,J,K,2) = S(I,J,K,2) / S(I,m,K,2)
           S(I,J,K,3) = S(I,J,K,3) / S(I,J,l,3)
        END DO
     END DO
  END DO

  !     Finally, precise 1s for the three high-end faces:

  DO K = 1, l
     DO J = 1, m
        S(n,J,K,1) = 1.0
     END DO
     
     DO I = 1, n
        S(I,m,K,2) = 1.0
     END DO
  END DO
  
  DO J = 1, m
     DO I = 1, n
        S(I,J,l,3) = 1.0
     END DO
  END DO

  ! Now get an average u,v,w

  ! Average u
  do j=1,m
     do k=1,l
        u = u + S(:,j,k,1)
     end do
  end do
  u = u/(l*m)

  ! Average v
  do i=1,n
     do k=1,l
        v = v + S(i,:,k,2)
     end do
  end do
  v = v/(n*l)

  ! Average w
  do i=1,n
     do j=1,m
        w = w + S(i,j,:,3)
     end do
  end do
  w = w/(n*m)


end subroutine para3d

