program main
  use module_lu_chps
  implicit none
  real(kind=8),dimension(:,:),allocatable :: A,L,U
  real(kind=8),dimension(:),allocatable :: b,x,dg
  integer :: i,j,p
  integer :: n

  integer :: nsys
  real(kind=8) :: h
  real(kind=8),dimension(:),allocatable :: xf

  n = 100
  allocate(xf(0:n))
  h=1./n
  do i=0,n
     xf(i) = i*h
  end do
  do i=0,n
     print*,i,xf(i)
  end do

  nsys =  N-1

  allocate(A(1:nsys,1:nsys))
  allocate(b(1:nsys))
  allocate(x(1:nsys))
  allocate(dg(1:nsys))
  A=0
  do i=1,nsys
     A(i,i) = -2
     b(i) =h**2
     b(i) =                        1*h**2
     if (i.ne.1)    A(i,i-1) = +1
     if (i.ne.nsys) A(i,i+1) = +1
  end do
  do i=1,nsys
     print'(10(f15.8,1x))',A(i,:),b(i)
  end do

  call SOLVE(A,X,B,DG)
  do i=1,nsys
     print'(2(f15.8,1x))',xf(i),x(i)
     write(33,'(2(e15.8,1x))'),xf(i),x(i)
  end do
  !call LU_FACTORIZE(A)
  !call L_SOLVE(A,Y,B)
  !call U_SOLVE(A,X,Y)
  
end program main
