module csr


implicit none

! real(kind=8) matrix
type csrp_t
   real(kind=8), allocatable :: mat(:) ! (nnz)
   integer, allocatable :: col_ind(:) ! (nnz)
   integer, allocatable :: row_ptr(:) ! (nrow+1)
end type csrp_t

contains
  
  subroutine init_csrp(spm, nrow, nnz)
    type(csrp_t), intent(out) :: spm
    integer, intent(in) :: nrow, nnz
    integer :: ierr
    
    allocate(spm%mat(nnz), stat=ierr)
    allocate(spm%col_ind(nnz), stat=ierr)
    allocate(spm%row_ptr(nrow+1), stat=ierr)
    spm%row_ptr(nrow+1) = nnz+1

    spm%mat=0
    spm%col_ind=0
    spm%row_ptr=0
  end subroutine init_csrp
  
  subroutine end_csrp(spm)
    type(csrp_t), intent(inout) :: spm
    integer :: ierr
    deallocate(spm%mat, stat=ierr)
    deallocate(spm%col_ind, stat=ierr)
    deallocate(spm%row_ptr, stat=ierr)
  end subroutine end_csrp
  
  subroutine csrpgemv(apply_transpose, zero_input_vector, a, spm, x, y)
    
    ! Perform y <- a*m*x + y or y <- a*m*x, where a is a constant, m is a
    ! sparse matrix and x and y are dense vectors.
    ! In:
    !   apply_transpose: If true then perform the matrix multiplication
    !      using the transpose of the matrix stored in spm.
    !   zero_input_vector: If true then y <- a*m*x is performed. Otherwise,
    !      y <- a*m*x + y is performed.
    !   a: constant factor which is included in the matrix multiplication.
    !   spm: sparse matrix (real(kind=8) in csr format. See module-level notes
    !      about storage format.
    !   x: dense vector.  If apply_transpose is false then the number of
    !      elements must be at least the number of columns in spm, else it
    !      must be at least the number of rows in spm.
    ! In/Out:
    !   y: dense vector.  Holds a*m*x or a*m*x + y on exit, depending on
    !      the input value of zero_input_vector. If apply_transpose is false
    !      then the number of elements must be at least the number of rows
    !      in spm, with all additional elements set to 0, else it must be
    !      at least the number of columns in spm, with all additional
    !      elements set to 0.
        
    logical, intent(in) :: apply_transpose, zero_input_vector
    real(kind=8), intent(in) :: a
    type(csrp_t), intent(in) :: spm
    real(kind=8), intent(in) :: x(:)
    real(kind=8), intent(inout) :: y(:)
    
    integer :: irow, icol, iz
    
    if (zero_input_vector) y = 0.0_8
    
    if (apply_transpose) then
       do irow = 1, size(spm%row_ptr)-1
          do iz = spm%row_ptr(irow), spm%row_ptr(irow+1)-1
             icol = spm%col_ind(iz)
             y(icol) = y(icol) + a*spm%mat(iz)*x(irow)
          end do
       end do
    else
       do irow = 1, size(spm%row_ptr)-1
          do iz = spm%row_ptr(irow), spm%row_ptr(irow+1)-1
             icol = spm%col_ind(iz)
             y(irow) = y(irow) + a*spm%mat(iz)*x(icol)
          end do
       end do
    end if
  end subroutine csrpgemv
  
  subroutine csrpgemv_single_row(spm, x, irow, y_irow)

    ! Calculate a single value in the vector y = m*x, where m is a sparse
    ! matrix and x and y are dense vectors.
    
    ! In:
    !   spm: sparse matrix (real(kind=8) in csr format. See module-level notes
    !        about storage format.
    !   x: dense vector.  Number of elements must be at least the number of
    !      columns in spm.
    !   irow: The index of the row of the Hamiltonian to multiply with.
    ! Out:
    !   y_irow: Holds \sum_j m_{irow,j}*x_j on exit.
    
    type(csrp_t), intent(in) :: spm
    real(kind=8), intent(in) :: x(:)
    integer, intent(in) :: irow
    real(kind=8), intent(out) :: y_irow
    
    integer :: icol, iz
    
    y_irow = 0.0_8
    do iz = spm%row_ptr(irow), spm%row_ptr(irow+1)-1
       icol = spm%col_ind(iz)
       y_irow = y_irow + spm%mat(iz)*x(icol)
    end do
    
  end subroutine csrpgemv_single_row
  
end module csr
