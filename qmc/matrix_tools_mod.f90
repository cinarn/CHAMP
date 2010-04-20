module matrix_tools_mod

  use basic_tools_mod
  use strings_tools_mod
  use variables_mod

  contains

! ==============================================================================
  subroutine inverse_by_svd (matrix, matrix_inv, dim, threshold)
! ------------------------------------------------------------------------------
! Description   : Calculate inverse of square matrix by SVD with a threshold on singular values
!
! Created       : J. Toulouse, 04 Nov 2005
! ------------------------------------------------------------------------------
  implicit none

! input
  real(dp), intent (in)      :: matrix (:,:)
  integer,          intent (in)      :: dim
  real(dp), intent (in)      :: threshold

! output
  real(dp), intent (out)      :: matrix_inv (:,:)

! local
  character(len=max_string_len_rout), save :: lhere = 'inverse_by_svd'
  integer i, j, k
  real(dp), allocatable :: mat_u(:,:)
  real(dp), allocatable :: mat_v(:,:)
  real(dp), allocatable :: mat_w(:)
  real(dp), allocatable :: mat_w_inv(:)
  integer w_kept_nb
  real(dp), allocatable :: mat_a(:,:)
  real(dp), allocatable :: mat_vt(:,:)
  real(dp), allocatable :: work (:)
  integer lwork
  integer info

! begin
  if (dim == 0) return

! temporary arrays for SVD
  call alloc ('mat_u', mat_u, dim, dim)
  call alloc ('mat_v', mat_v, dim, dim)
  call alloc ('mat_w', mat_w, dim)
  call alloc ('mat_w_inv', mat_w_inv, dim)

  lwork = 10 * dim
  call alloc ('mat_a', mat_a, dim, dim)
  call alloc ('mat_vt', mat_vt, dim, dim)
  call alloc ('mat_w', mat_w, dim)
  call alloc ('work', work, lwork)

! SVD from numerical recipes
!  mat_u = matrix ! matrix to be inverted
!  call svdcmp (mat_u, dim, dim, dim, dim, mat_w, mat_v)
!
!  write(6,*) trim(lhere),': SVD from numerical recipes:'
!  write(6,*) trim(lhere),': mat_u=',mat_u
!  write(6,*) trim(lhere),': mat_v=',mat_v
!  write(6,*) trim(lhere),': mat_w=',mat_w

! SVD from Lapack
  mat_a = matrix ! matrix to be inverted
!  write(6,*) trim(lhere),': mat_a=',mat_a
!  write(6,*) trim(lhere),': before dgesvd'
  call dgesvd( 'A', 'A', dim, dim, mat_a, dim, mat_w, mat_u, dim, mat_vt, dim, work, lwork, info)
!  write(6,*) trim(lhere),': after dgesvd'
  if (info /= 0) then
   call die (lhere, 'problem in dgesvd')
  endif

 mat_v = transpose(mat_vt)

!  write(6,*) trim(lhere),': SVD from Lapack:'
!  write(6,*) trim(lhere),': mat_u=',mat_u
!  write(6,*) trim(lhere),': mat_v=',mat_v
!  write(6,*) trim(lhere),': mat_w=',mat_w

! Singular values
!JT  do i = 1, dim
!JT   write(6,*) trim(lhere), ': i=',i,' mat_w=',mat_w(i)
!JT  enddo

! Inverse singular values  (drop small singular values)
!JT  write(6,*) trim(lhere), ': threshold on singular values: ',threshold
  w_kept_nb = dim
  do i = 1, dim
   if (mat_w (i) < threshold) then
     mat_w_inv (i) = 0.d0
     w_kept_nb = w_kept_nb - 1
   else
     mat_w_inv (i) = 1.d0 / mat_w (i)
   endif
  enddo
!JT  write(6,*) trim(lhere), ': number of singular values dropped:', dim - w_kept_nb

! calculate inverse matrix
  matrix_inv = 0.d0

   do i = 1, dim
     do j = 1, dim
       do k = 1, dim
         matrix_inv (i, j) = matrix_inv (i, j) + mat_v (i, k) * mat_w_inv (k) * mat_u (j, k)
       enddo
     enddo
   enddo

!  write(6,*) trim(lhere), ': matrix_inv=',matrix_inv

! release arrays for SVD
  call release ('mat_u', mat_u)
  call release ('mat_v', mat_v)
  call release ('mat_w', mat_w)
  call release ('mat_w_inv', mat_w_inv)
  call release ('mat_a', mat_a)
  call release ('mat_vt', mat_vt)
  call release ('work', work)

!  write(6,*) trim(lhere), ': exiting'
  end subroutine inverse_by_svd

! ==============================================================================
  subroutine eigensystem (matrix, eigenvectors, eigenvalues, dim)
! ------------------------------------------------------------------------------
! Description   : compute eigenvectors and eigenvalues of a real symmetrix matrix
! Description   : using Lapack routine
!
! Created       : J. Toulouse, 11 Jan 2006
! ------------------------------------------------------------------------------
  implicit none

! input
  real(dp), intent(in)   :: matrix (:,:)
  integer,  intent(in)   :: dim

! output
  real(dp), intent(out)  :: eigenvectors (:,:)
  real(dp), intent(out)  :: eigenvalues (:)

! local
  character(len=max_string_len_rout), save :: lhere = 'eigensystem'
  real(dp), allocatable :: mat_a(:,:)
  real(dp), allocatable :: work (:)
  real(dp) :: matrix_check
  integer lwork, info, i, j, k
  integer, save :: warnings_nb = 0
  integer       :: warnings_nb_max = 20

! begin
  if (dim == 0) return

  lwork = 10 * dim
  call alloc ('mat_a', mat_a, dim, dim)
  call alloc ('work', work, lwork)

  mat_a (:,:) = matrix (:,:)

  call dsyev ('V','U',dim, mat_a, dim, eigenvalues, work, lwork, info)
  if (info /= 0) then
   call die (lhere, 'exiting dsyev with info='+info+' /= 0.')
  endif

  eigenvectors (:,:) = mat_a (:,:)

! checkings
!  call is_a_number_or_die ('eigenvalues', eigenvalues)

! check recovery of original matrix after diagonalization
  do i = 1, dim
   do j = 1, dim
     matrix_check = 0.d0
     do k = 1, dim
      matrix_check = matrix_check + eigenvectors (i, k) * eigenvalues (k) * eigenvectors (j, k)
     enddo ! k
     if(warnings_nb < warnings_nb_max .and. abs(matrix_check-matrix(i,j)) > 1.d-7) then
       write(6,'(''Warning: low accuracy in diagonalization; the error on a matrix element'',2i4,'' is'',d12.4)') &
         i,j,matrix_check-matrix(i,j)
       l_warning = .true.
       warnings_nb = warnings_nb + 1
       if (warnings_nb == warnings_nb_max) then
        write(6,'(a)') 'all further similar warnings will be suppressed'
       endif
     endif
! JT: Warning: comment out this stop for now
!     if(abs(matrix_check-matrix(i,j)) > 1.d-2) then
!       call die (lhere, 'low accuracy in diagonalization; the error on a matrix element is '+abs(matrix_check-matrix(i,j)))
!     endif
   enddo ! j
  enddo ! i

! release arrays
  call release ('mat_a', mat_a)
  call release ('work', work)

  end subroutine eigensystem

! ==============================================================================
  subroutine to_the_power (matrix, dim, power_n, matrix_out)
! ------------------------------------------------------------------------------
! Description   : returns the matrix to the power n
! Description   : valid for a real symmetric matrix (uses eigensystem)
!
! Created       : B. Mussard, 09 Mar 2010
! ------------------------------------------------------------------------------
  implicit none

! input
  real(dp), intent(in)   :: matrix (:,:)
  integer,  intent(in)   :: dim
  real(dp), intent(in)   :: power_n

! output
  real(dp), intent(out)  :: matrix_out (:,:)

! local
  character(len=max_string_len_rout), save :: lhere = 'to_the_power'
  integer eigen_i, bas_i, bas_j, bas_k
  real(dp), allocatable  :: eigenvectors (:,:)
  real(dp), allocatable  :: eigenvalues (:)

! begin
  call alloc('eigenvectors',eigenvectors,dim,dim)
  call alloc('eigenvalues',eigenvalues,dim)

  call eigensystem (matrix, eigenvectors, eigenvalues, dim)
 
  if (power_n < 1) then
     do eigen_i = 1, dim
        if (eigenvalues (eigen_i) < (-10**-5)) then
            write(6,'(a)') 'one eigenvalue is too negative to be equalled to zero'
        else if (eigenvalues (eigen_i) < 0) then
            eigenvalues (eigen_i) = 0
        endif
     enddo
  endif
            

  do bas_i = 1, dim
    do bas_j = 1, dim
       matrix_out(bas_i, bas_j) = 0.d0
       do bas_k = 1, dim
          matrix_out(bas_i, bas_j) = matrix_out(bas_i, bas_j) + eigenvectors (bas_i, bas_k) * eigenvalues(bas_k)**power_n * eigenvectors (bas_j, bas_k)
       enddo ! bas_k
    enddo ! bas_j
  enddo ! bas_i

  end subroutine to_the_power

! ==============================================================================
  real(dp) function matrix_determinant (matrix, n)
! ------------------------------------------------------------------------------
! Description   : returns determinant of a matrix
! WARNING: not checked!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Function to find the determinant of a square matrix
! Author : Louisda16th a.k.a Ashwith J. Rego
! Description: The subroutine is based on two key points:
! 1] A determinant is unaltered when row operations are performed: Hence, using this principle,
! row operations (column operations would work as well) are used
! to convert the matrix into upper traingular form
! 2]The determinant of a triangular matrix is obtained by finding the product of the diagonal elements
!
! Created       : J. Toulouse, 29 Oct 2009
! ------------------------------------------------------------------------------
  IMPLICIT NONE
  REAL(dp), DIMENSION(n,n) :: matrix
  INTEGER, INTENT(IN) :: n
  REAL(dp) :: m, temp
  INTEGER :: i, j, k, l
  LOGICAL :: DetExists = .TRUE.

  l = 1
  !Convert to upper triangular form
  DO k = 1, n-1
          IF (matrix(k,k) == 0) THEN
                  DetExists = .FALSE.
                  DO i = k+1, n
                          IF (matrix(i,k) /= 0) THEN
                                  DO j = 1, n
                                          temp = matrix(i,j)
                                          matrix(i,j)= matrix(k,j)
                                          matrix(k,j) = temp
                                  END DO
                                  DetExists = .TRUE.
                                  l=-l
                                  EXIT
                          ENDIF
                  END DO
                  IF (DetExists .EQV. .FALSE.) THEN
                          matrix_determinant = 0
                          return
                  END IF
          ENDIF
          DO j = k+1, n
                  m = matrix(j,k)/matrix(k,k)
                  DO i = k+1, n
                          matrix(j,i) = matrix(j,i) - m*matrix(k,i)
                  END DO
          END DO
  END DO
  
  !Calculate determinant by finding product of diagonal elements
  matrix_determinant = l
  DO i = 1, n
          matrix_determinant = matrix_determinant * matrix(i,i)
  END DO
             
  end function matrix_determinant

! ==============================================================================
  real(dp) function trace (matrix)
! ------------------------------------------------------------------------------
! Description   : returns trace of a square matrix
!
! Created       : J. Toulouse, 20 Jan 2006
! ------------------------------------------------------------------------------
  implicit none

! input
  real(dp)             , intent (in)   :: matrix (:,:)

! local
  character(len=max_string_len_rout), save :: lhere = 'trace'
  integer dim1, dim2
  integer i

! begin
  dim1 = size(matrix,1)
  dim2 = size(matrix,2)

  if (dim1 /= dim2) then
   write(6,*) trim(lhere),': in routine', trim(here),' trace requested of a non-square matrix!'
   write(6,*) trim(lhere),': dim1=',dim1,' /= dim2=',dim2
   call die (lhere)
  endif

  trace = 0.d0
  do i = 1, dim1
   trace = trace + matrix (i,i)
  enddo

  end function trace

! ==============================================================================
  integer function kronecker_delta (i, j)
! ------------------------------------------------------------------------------
! Description   : kronecker_delta
!
! Created       : J. Toulouse, 19 Mar 2006
! ------------------------------------------------------------------------------
  implicit none

! input
  integer, intent (in)   :: i, j

! begin
  if (i == j) then
   kronecker_delta = 1
  else
   kronecker_delta = 0
  endif

  end function kronecker_delta

end module matrix_tools_mod
