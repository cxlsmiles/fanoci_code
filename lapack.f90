module lapack
    implicit none
    contains

!************************************************
!  A LAPACK subroutine for diagonalization
!  of real symmetric matrices N x N
!************************************************
    subroutine RealSymm(mat, N, eigval)
        integer                                 :: N
        real(kind=8), dimension(N, N)           :: mat
        real(kind=8), dimension(N)              :: eigval
        integer                                 :: INFO, LWORK
        real(kind=8), allocatable, dimension(:) :: WORK

        LWORK = 5*N
        allocate(WORK(LWORK))

        call dsyev('V', 'U', N, mat, N, eigval, WORK, LWORK, INFO)
                !JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO
!*  JOBZ    (input) CHARACTER*1
!*          = 'N':  Compute eigenvalues only;
!*          = 'V':  Compute eigenvalues and eigenvectors.
!*
!*  UPLO    (input) CHARACTER*1
!*          = 'U':  Upper triangle of A is stored;
!*          = 'L':  Lower triangle of A is stored.
!*
!*  N       (input) INTEGER
!*          The order of the matrix A.  N >= 0.
!*
!*  A       (input/output) DOUBLE PRECISION array, dimension (LDA, N)
!*          On entry, the symmetric matrix A.  If UPLO = 'U', the
!*          leading N-by-N upper triangular part of A contains the
!*          upper triangular part of the matrix A.  If UPLO = 'L',
!*          the leading N-by-N lower triangular part of A contains
!*          the lower triangular part of the matrix A.
!*          On exit, if JOBZ = 'V', then if INFO = 0, A contains the
!*          orthonormal eigenvectors of the matrix A.
!*          If JOBZ = 'N', then on exit the lower triangle (if UPLO='L')
!*          or the upper triangle (if UPLO='U') of A, including the
!*          diagonal, is destroyed.
!*
!*  LDA     (input) INTEGER
!*          The leading dimension of the array A.  LDA >= max(1,N).
!*
!*  W       (output) DOUBLE PRECISION array, dimension (N)
!*          If INFO = 0, the eigenvalues in ascending order.
!*
!*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
!*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!*
!*  LWORK   (input) INTEGER
!*          The length of the array WORK.  LWORK >= max(1,3*N-1).
!*          For optimal efficiency, LWORK >= (NB+2)*N,
!*          where NB is the blocksize for DSYTRD returned by ILAENV.
!*
!*          If LWORK = -1, then a workspace query is assumed; the routine
!*          only calculates the optimal size of the WORK array, returns
!*          this value as the first entry of the WORK array, and no error
!*          message related to LWORK is issued by XERBLA.
!*
!*  INFO    (output) INTEGER
!*          = 0:  successful exit
!*          < 0:  if INFO = -i, the i-th argument had an illegal value
!*          > 0:  if INFO = i, the algorithm failed to converge; i
!*                off-diagonal elements of an intermediate tridiagonal
!*                form did not converge to zero.
    end subroutine RealSymm

end module lapack
