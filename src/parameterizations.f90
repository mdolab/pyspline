subroutine para3d(X,n,m,l,ndim,S,u,v,w)
  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract para3d calculates the parametric locations for a 3d block
  !
  !     Description of Arguments
  !     Input
  !     X       - Real,size(ndim,l,m,n): Coordiantes
  !     Output
  !     S       - Real,size(n,m,l,3): The u,v,w parametric positions
  !     u       - Real, size(n): The averaged u parameters
  !     v       - Real, size(m): The averaged v parameters
  !     w       - Real, size(l): The averaged w parameters

  implicit none

  ! Input
  integer          , intent(in)   :: n,m,l,ndim
  double precision , intent(in)   :: X(ndim,l,m,n)
  

  ! Output
  double precision , intent(out)  :: S(ndim,l,m,n)
  double precision , intent(out)  :: u(n),v(m),w(l)

  ! Working 
  integer                         :: i,j,k

  double precision DELI,DELJ,DELK

  DELI(K,J,I) = SQRT ((X(1,K,J,I) - X(1,K,J,I-1)) ** 2 + &
                      (X(2,K,J,I) - X(2,K,J,I-1)) ** 2 + &
                      (X(3,K,J,I) - X(3,K,J,I-1)) ** 2)

  DELJ(K,J,I) = SQRT ((X(1,K,J,I) - X(1,K,J-1,I)) ** 2 + &
                      (X(2,K,J,I) - X(2,K,J-1,I)) ** 2 + &
                      (X(3,K,J,I) - X(3,K,J-1,I)) ** 2)

  DELK(K,J,I) = SQRT ((X(1,K,J,I) - X(1,K-1,J,I)) ** 2 + &
                      (X(2,K,J,I) - X(2,K-1,J,I)) ** 2 + &
                      (X(3,K,J,I) - X(3,K-1,J,I)) ** 2)


  u(:) = 0.0
  v(:) = 0.0
  w(:) = 0.0

  !     Zero the three low-end faces (or edges if one plane is specified).
  
  DO J = 1, m
     DO K = 1, l
        S(1,K,J,1) = 0.0
     END DO
     
     DO I = 1, n
        S(2,K,1,I) = 0.0
     END DO
  END DO
  

  DO I = 1, n
     DO J = 1, m
        S(3,1,J,I) = 0.0
     END DO
  END DO
  
  !     Set up the low-end edge lines because they are missed by the
  !     following loops over most of the low-end faces:

  DO I = 2, n
     S(1,1,1,I) = S(1,1,1,I-1) + DELI(1,1,I)
  END DO

  DO J = 2, m
     S(2,1,J,1) = S(2,1,J-1,1) + DELJ(1,J,1)
  END DO

  DO K = 2, l
     S(3,K,1,1) = S(3,k-1,1,1) + DELK(K,1,1)
  END DO

  !     Set up the rest of the low-end face lines because they are
  !     missed by the the main loop over most of the volume.

  DO K = 2, l
     DO J = 2, m
        S(2,K,J,1) = S(2,K,J-1,1) + DELJ(K,J,1)
        S(3,K,J,1) = S(3,K-1,J,1) + DELK(K,J,1)
     END DO
  end DO

  DO I = 2, n
     DO K = 2, l
        S(1,K,1,I) = S(1,K,1,I-1) + DELI(K,1,I)
        S(3,K,1,I) = S(3,K-1,1,I) + DELK(K,1,I)
     END DO

  END DO

  DO I = 2, n
     DO J = 2, m
        S(1,1,J,I) = S(1,1,J,I-1) + DELI(1,J,I)
        S(2,1,J,I) = S(2,1,J-1,I) + DELJ(1,J,I)
     END DO
  END DO

  !     Traverse the block just once for all lines except those within
  !     the low-end faces.

  
  DO I = 2, n
     DO J = 2, m
        DO K = 2, l
           S(1,K,J,I) = S(1,K,J,I-1) + DELI(K,J,I)
           S(2,K,J,I) = S(2,K,J-1,I) + DELJ(K,J,I)
           S(3,K,J,I) = S(3,K-1,J,I) + DELK(K,J,I)
        END DO
     END DO
  END DO

  !     Normalizing requires another pass through the volume.
  !     Handle lines of zero length first by inserting uniform
  !     distributions.  Then the standard normalization can be
  !     applied safely everywhere.

  DO J = 1, m
     !        Zero-length lines in the I direction?

     DO K = 1, l
        IF (S(1,K,J,n) == 0.0) THEN
           DO I = 2, n
              S(1,K,J,I) = I - 1
           END DO
        END IF
     END DO
  end DO
!      !        Zero-length lines in the J direction?

  DO I = 1, n
     DO K = 1, l
        IF (S(2,K,m,I) == 0.0) THEN
           DO J = 2, m
              S(2,K,J,I) = J - 1
           END DO
        END IF
     END DO
  END DO

!   !     Zero-length lines in the K direction?

  DO I = 1, n
     DO J = 1, m
        IF (S(3,l,J,I) == 0.0) THEN
           DO K = 2, l
              S(3,K,J,I) = K - 1
           END DO
        END IF
     END DO
  END DO

  !     Normalize:

  DO I = 1, n
     DO J = 1, m
        DO K = 1, l
           S(1,K,J,I) = S(1,K,J,I) / S(1,K,J,N)
           S(2,K,J,I) = S(2,K,J,I) / S(2,K,M,I)
           S(3,K,J,I) = S(3,K,J,I) / S(3,L,J,I)
        END DO
     END DO
  END DO

  !     Finally, precise 1s for the three high-end faces:

  DO J = 1, m
     DO K = 1, l
        S(1,K,J,n) = 1.0
     END DO
  end DO

  DO I = 1, n
     DO K = 1, l
        S(2,K,M,I) = 1.0
     END DO
  END DO
  
  DO I = 1, n
     DO J = 1, m
        S(3,L,J,I) = 1.0
     END DO
  END DO

  ! Now get an average u,v,w

  ! Average u

  do j=1,m
     do k=1,l
        u = u + S(1,k,j,:)
     end do
  end do
  u = u/(l*m)

  ! Average v
  do i=1,n
     do k=1,l
        v = v + S(2,k,:,i)
     end do
  end do
  v = v/(n*l)

  ! Average w
  do i=1,n
     do j=1,m
        w = w + S(3,:,j,i)
     end do
  end do
  w = w/(n*m)


end subroutine para3d

