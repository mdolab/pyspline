subroutine fit_foil(M,N,A,X1,X2,Y1,Y2,COEFS)
!-----------------------------------------------------------------------
! Purpose: Fit_foil performs a constrained Least Mean Squares Fit on a
! set of x-y points in 2-space defined by a parametric value s. The
! main functionality is the enformenet of 1st order continutity at the
! join of the curves. This is useful when two splines are required to
! define the upper and lower surfaces of an airfoil since we want to
! ensure at least C1 continutity. 
!
!
!
! Inputs:
!
! A: This is the matrix of partial derivaties of size M by N
!
! M: This is the number of Datapoints we are fitting
!
! N: This is the number of control points plus derivative constraints
! we have
!
! 
  
integer, intent(in):: M
integer, intent(in):: N

double precision, intent(in) :: A(M,N)
double precision, intent(in) :: X1(M)
double precision, intent(in) :: Y1(M)
double precision, intent(in) :: X2(M)
double precision, intent(in) :: Y2(M)
double precision, intent(out):: COEFS(N,4)
integer P
double precision :: A_temp(4*M,4*N)
double precision :: X(4*M)
double precision :: B(10,4*N)
double precision :: D(10)
double precision :: temp(4*N)
!-- variables needed for QR decomposition
integer :: lwork, info
double precision :: work(1+4*M+4*N+10)
double precision err
double precision xa,xb,xc,ya,yb,yc,tangent
character trans

external dgglse

! Code for doing constrained least squares
P = 10
lwork = 1+4*M+4*N+10

B(:,:) = 0.0
B(1,    1) = 1
B(2,  N  ) = 1
B(3,  N+1) = 1
B(4,2*N  ) = 1
B(5,2*N+1) = 1
B(6,3*N  ) = 1
B(7,3*N+1) = 1
B(8,4*N  ) = 1
B(9,2)     = 1
B(10,N+2) = 1
D(1) = x1(1)
D(2) = x1(M)
D(3) = x2(1)
D(4) = x2(M)
D(5) = y1(1)
D(6) = y1(M)
D(7) = y2(1)
D(8) = y2(M)
D(9) = 0
D(10) = 0
A_temp(:,:) = 0
A_temp(    1:  M,    1:  N) = A
A_temp(  M+1:2*M,  N+1:2*N) = A
A_temp(2*M+1:3*M,2*N+1:3*N) = A
A_temp(3*M+1:4*M,3*N+1:4*N) = A

X(    1:  M) = X1
X(  M+1:2*M) = X2
X(2*M+1:3*M) = Y1
X(3*M+1:4*M) = Y2

! print *,'before the call:'
! print *,'M:',M
! print *,'N:',N
! print *,'P:',P
! print *,'A:',A
! print *,'A_temp:',A_temp
! print *,'B:',B
! print *,'X:',X
! print *,'D:',D

!call dgglse(M, N, P, A, LDA, B, LDB, C , D, X,    WORK, LWORK, INFO )
call  dgglse(4*M, 4*N, P, A_temp, 4*M,   B, P, X, D, temp, work, lwork, info )

print *,'info1:',info
coefs(:,1) = temp(    1:  N)
coefs(:,3) = temp(  N+1:2*N)
coefs(:,2) = temp(2*N+1:3*N)
coefs(:,4) = temp(3*N+1:4*N)

! Now calculate the le constraint

xa = coefs(2,2)
xb = coefs(1,1)
xc = coefs(1,2)

ya = coefs(3,2)
yb = coefs(2,1)
yc = coefs(2,2)

tangent = xa*(yb-yc)+xb*(yc-ya)+xc*(ya-yb)
print *,'tangent:',tangent

end subroutine fit_foil
