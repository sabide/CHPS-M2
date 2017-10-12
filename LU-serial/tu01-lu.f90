program main
  use module_lu_chps
  implicit none
  real(kind=8),dimension(:,:),allocatable :: A,L,U
  real(kind=8),dimension(:),allocatable :: b,x,y
  integer :: i,j,p
  integer :: nx,ny

  allocate(A(1:4,1:4))
  allocate(L(1:4,1:4))
  allocate(U(1:4,1:4))
  allocate(b(1:4),x(1:4),y(1:4))
  
  A(1,1:4) = (/2,1,1,0/)
  A(2,1:4) = (/4,3,3,1/)
  A(3,1:4) = (/8,7,9,5/)
  A(4,1:4) = (/6,7,9,8/)

  call random_number(x)
  b = matmul(A,X)

  do i=1,4
     print*,x(i),b(i)
  end do
  
  call LU_FACTORIZE(A)
  
  L=0
  L(1,1)=1
  DO I=2,4
     L(I,I)=1
     L(I,1:I-1)=A(I,1:I-1)
  END DO
  U=0
  DO I=1,4
     U(I,I:4)=A(I,I:4)
  END DO
  A=matmul(L,U)
  
  
  call L_SOLVE(A,Y,B)
  call U_SOLVE(A,X,Y)
  
  do i=1,4
     print*,A(i,1:4),x(i)
  end do




!!$  nx=5
!!$  ny=5
!!$  call build_laplacian(A,nx,ny)
!!$  do i=1,nx*ny
!!$     print'(16(e15.8,1x))',A(i,1:nx*ny)
!!$  end do
!!$  deallocate(b,x,y)
!!$  allocate(b(nx*ny),x(nx*ny),y(nx*ny))
!!$  b=1
!!$  call LU_FACTORIZE(A)
!!$  call L_SOLVE(A,Y,B)
!!$  call U_SOLVE(A,X,Y)
!!$
!!$  do i=1,nx
!!$     do j=1,ny
!!$        p = i + (j-1)*nx
!!$        print*,i,j,x(p)
!!$        write(20,*),i,j,x(p)
!!$     end do
!!$     write(20,*)
!!$  end do

end program main
