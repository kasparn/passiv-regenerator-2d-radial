module matrix_inversion
IMPLICIT NONE
contains

subroutine invMatrixLUD( A, Ainv, B, X, n )
real,dimension(n,n),intent(in) :: A
real,dimension(n,n),intent(inout) :: Ainv
real,dimension(n),intent(in) :: B
real,dimension(n),intent(inout) :: X
integer,intent(in) :: n
integer :: i,j
integer :: info, LDA
integer,allocatable,dimension(:) :: IPIV
integer :: DeAllocateStatus
real :: t1,t2,t3

external DGETRF
external DGETRS
call cpu_time( t1 )
LDA = n


ALLOCATE (IPIV(n))
ipiv = 0
!::
!::DGETRF computes an LU factorization of a general M-by-N matrix A
!::using partial pivoting with row interchanges.

LDA=N

!::Store A in Ainv to prevent it from being overwritten by LAPACK

Ainv = A
CALL DGETRF( n, n, Ainv, n, IPIV, INFO )

IF(INFO.EQ.0)THEN
 PRINT '(" LU decomposition successful ")'
ENDIF
IF(INFO.LT.0)THEN
 PRINT '(" LU decomposition:  illegal value ")'
 STOP
ENDIF
IF(INFO.GT.0)THEN
 WRITE(*,*) INFO
!35       FORMAT( 'LU decomposition: U(',I4,',',I4,') = 0 ')
ENDIF

!write(*,*) ipiv

!::DGETRI computes the inverse of a matrix using the LU factorization
!::computed by DGETRF.
!call cpu_time( t2 )
write(*,*) 'after LU dec'
!write(*,*) t2-t1
!CALL DGETRI(N, Ainv, N, IPIV, WORK, LWORK, INFO)
X = B
call dgetrs( 'N', n, 1, Ainv, n, IPIV, X, n, info )
write(*,*) 'info',info
call cpu_time( t3 )
write(*,*) 'after dgetri'
write(*,*) t3-t1

PRINT '(" ")'
IF (info.NE.0) THEN
 stop 'Matrix inversion failed!'
ELSE
 PRINT '(" Inverse Successful ")'
ENDIF
stop
DEALLOCATE (IPIV, STAT = DeAllocateStatus)

PRINT '(" Done ")'

END subroutine invMatrixLUD

end module matrix_inversion
