program main
  use csr
  implicit none
  type(csrp_t) :: a
  integer i,j,nrows,nnz,n,pp
  real(kind=8),dimension(:),allocatable :: x,y,b
  real(kind=8) :: alf,bet,r_dot_r
  real(kind=8),dimension(:),allocatable :: r,p,q
  integer :: iter,nn_diag
  !>================================
  n=200
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
        nn_diag = nnz
        a%mat(nn_diag) = -4
        a%col_ind(nnz) = nrows
        a%row_ptr(nrows+1) = a%row_ptr(nrows) + 1
        if (i.ne.1) then
           nnz = nnz+1
           a%mat(nnz) = +1
!           a%mat(nn_diag) = a%mat(nn_diag) +1
           a%col_ind(nnz) = nrows-1
           a%row_ptr(nrows+1) = a%row_ptr(nrows+1) + 1
        end if
        if (i.ne.n) then
           nnz = nnz+1
           a%mat(nnz) = +1
 !          a%mat(nn_diag) = a%mat(nn_diag) +1
           a%col_ind(nnz) = nrows+1
           a%row_ptr(nrows+1) = a%row_ptr(nrows+1) + 1
        end if
        
        if (j.ne.1) then
           nnz = nnz+1
           a%mat(nnz) = +1
  !         a%mat(nn_diag) = a%mat(nn_diag) +1
           a%col_ind(nnz) = nrows-n
           a%row_ptr(nrows+1) = a%row_ptr(nrows+1) + 1
        end if
        
        if (j.ne.n) then
           nnz = nnz+1
           a%mat(nnz) = +1
  !         a%mat(nn_diag) = a%mat(nn_diag) +1
           a%col_ind(nnz) = nrows+n
           a%row_ptr(nrows+1) = a%row_ptr(nrows+1) + 1
        end if
     end do
  end do
  
  !print * , a%mat
  !x=1
  !call csrpgemv(apply_transpose=.false., zero_input_vector=.true., a=1d0, spm=a, x=x, y=y)
  


  x=0
  allocate( p(n*n) )
  allocate( q(n*n) )
  allocate( r(n*n) )
  allocate( b(n*n) )
  !> r = b-ax 
  b=-real(1)/real(n*n)
  r = b 
  call csrpgemv(apply_transpose=.false., zero_input_vector=.false., a=-1.0, spm=a, x=x, y=r)
  p = r
  
  
  do iter=1,100000
     !> q = A r  
     call csrpgemv(apply_transpose=.false., zero_input_vector=.true., a=+1.0, spm=a, x=p, y=q)
     r_dot_r = dot_product(r,r)
     alf = r_dot_r/dot_product(p,q)
     x = x + alf*p
     r = r - alf*q
     bet = dot_product(r,r)/r_dot_r
     p = r + bet*p
     if (sqrt(dot_product(r,r)).lt.1e-20) exit
     !print*,sqrt(dot_product(r,r))
  end do



  pp=0
  do j=1,n
     do i=1,n
        pp=pp+1
        print'(3(e15.8,1x))',real(i),real(j),x(pp)
     end do
  end do
end program main
