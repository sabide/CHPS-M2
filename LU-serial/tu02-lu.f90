program main
  use module_lu_chps
  implicit none
  real(kind=8),dimension(:,:),allocatable :: As,A,L,U
  real(kind=8),dimension(:),allocatable :: b,x,y,xs
  real(kind=8) :: err_x,err_lu,t1,t2,t_lu
  integer :: i,j,p
  integer :: nx,ny,n


  ! INITIALISATION
  n=1000
  allocate(A (1:n,1:n))
  allocate(L (1:n,1:n))
  allocate(U (1:n,1:n))
  allocate(As(1:n,1:n))
  allocate(b(1:n),x(1:n),y(1:n),xs(1:n))


  ! ALEATOIRE
  call random_number(A)
  DO p=1,1
     A= matmul(A,A)
  end DO
  
  call random_number(x)
  b = matmul(A,X)

  
  
  !> factorisation et sauvegarde de A
  As = A
  call cpu_time(t1)
  call LU_FACTORIZE(A)
  call cpu_time(t2)
  t_lu = t2-t1
  !> résolution et sauvegarde de x
  XS=X
  call L_SOLVE(A,Y,B)
  call U_SOLVE(A,X,Y)

  !> erreur sur x 
  err_x = maxval(abs(XS-X))


  !> calcul de l'erreur de la reconstruction
  L=0
  L(1,1)=1
  DO I=2,N
     L(I,I)=1
     L(I,1:I-1)=A(I,1:I-1)
  END DO
  
  U=0
  DO I=1,N
     U(I,I:N)=A(I,I:N)
  END DO

  A=matmul(L,U)
  err_lu = maxval(abs(AS-A))



  print*,n,err_x,err_lu,t_lu
  

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
