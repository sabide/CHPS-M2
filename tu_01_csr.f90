program main
  use csr
  implicit none
  type(csrp_t) :: a
  integer i,j,nrows,nnz,n,p
  real(kind=8),dimension(:),allocatable :: x,y,b
  
  
  n=100
  allocate( x(n*n) )
  allocate( y(n*n) )
  nnz=0
  nrows=0
  do j=1,n
     do i=1,n
        nrows = nrows + 1
        nnz = nnz + 1
        if (i.ne.1) nnz = nnz+1
        if (i.ne.n) nnz = nnz+1
        if (j.ne.1) nnz = nnz+1
        if (j.ne.n) nnz = nnz+1
     end do
  end do
  call init_csrp(spm=a, nrow=nrows, nnz=nnz)
  
  nnz=0
  nrows=0
  a%row_ptr(1) = 1
  do j=1,n
     do i=1,n
        
        nrows = nrows + 1
        nnz = nnz + 1
        a%mat(nnz) = 4
        a%col_ind(nnz) = nrows
        a%row_ptr(nrows+1) = a%row_ptr(nrows) + 1
        if (i.ne.1) then

        end if
        if (i.ne.n) then

        end if
        
        if (j.ne.1) then

        end if
        if (j.ne.n) then
           
        end if
     end do
  end do
  
  !print * , a%mat
  x=1
  call csrpgemv(apply_transpose=.false., zero_input_vector=.true., a=1d0, spm=a, x=x, y=y)

  p=0
  do j=1,n
     do i=1,n
        p=p+1
        print'(3(e15.8,1x))',real(i),real(j),y(p)
     end do
  end do
end program main
