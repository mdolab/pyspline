      SUBROUTINE BTPCF_mod(X,NDATA,FCN,LDF,NF,T,BCOEF,WORK,gk)
C***BEGIN PROLOGUE  BTPCF
C***REFER TO  B2INK,B3INK
C***ROUTINES CALLED  BINTK,BNSLV
C***END PROLOGUE  BTPCFdf
C
C  -----------------------------------------------------------------
C  BTPCF COMPUTES B-SPLINE INTERPOLATION COEFFICIENTS FOR NF SETS
C  OF DATA STORED IN THE COLUMNS OF THE ARRAY FCN. THE B-SPLINE
C  COEFFICIENTS ARE STORED IN THE ROWS OF BCOEF HOWEVER.
C  EACH INTERPOLATION IS BASED ON THE N ABCISSA STORED IN THE
C  ARRAY X, AND THE N+K KNOTS STORED IN THE ARRAY T. THE ORDER
C  OF EACH INTERPOLATION IS K. THE WORK ARRAY MUST BE OF LENGTH
C  AT LEAST 2*K*(N+1).
C  -----------------------------------------------------------------
C
C  ------------
C  DECLARATIONS
C  ------------
C
C  PARAMETERS
C
      INTEGER
     *        N, LDF, K
      REAL
     *     X(NDATA), FCN(LDF,NF), T(1),BCOEF(NF+2,NDATA+2),WORK(100000)
!     *     WORK2(5*(NDATA+2))
C
      logical gk
C  LOCAL VARIABLES
C
      INTEGER
     *        I, J, K1, K2, IQ, IW
C
C  ---------------------------------------------
C  CHECK FOR NULL INPUT AND PARTITION WORK ARRAY
C  ---------------------------------------------
C
C***FIRST EXECUTABLE STATEMENT
      IF (NF .LE. 0)  GO TO 500
      k=4
      K1 = K - 1
      K2 = K1 + K
      IQ = 1 + NDATA+2
      IW = IQ + K2*(NDATA+2)+1
C  -----------------------------
C  COMPUTE B-SPLINE COEFFICIENTS
C  -----------------------------
C
C     
C   FIRST DATA SET
C
      print *,'in btpcf_mod.f'
C  ALL REMAINING DATA SETS BY BACK-SUBSTITUTION
C     


      DO 100 J=1,NF
         !print *,'j:',j
         !print *,'fcn:',fcn(:,j)
         CALL BINT4(X,FCN(:,j),NDATA,2,2,0.0,0.0,1,T,WORK,
     *     N,K,WORK(IW))
         !print *,'coefs:'
         DO 60 I=1,N
            BCOEF(J,I) = WORK(I)
          !  print *,work(i)
   60    CONTINUE
  100 CONTINUE
C
      print *,'nf,N:',nf,N
c      do i=1,N
c         bcoef(1,i) = bcoef(2,i)
c         bcoef(NF+2,i) = bcoef(NF+1,i)
c      end do
C     
      print *,'bceof:'
      do i=1,nf+2
         do j=1,ndata+2
            print *,'i,j,bcoef:',i,j,bcoef(i,j)
         end do
      end do
C  ----
C  EXIT
C  ----
C
  500 CONTINUE
      RETURN
      END
