module module_lu_chps

contains

  !> Abide => bug sur la numerotation non diagonale i \in {2,5} j\in {1,4} par exemple 
  subroutine LU_FACTORIZE(A)
    implicit none
    real(kind=8),dimension(:,:),allocatable :: A
    integer :: i,j,k,n
    integer :: is(2),ie(2)
    
    is = lbound(A)
    ie = ubound(A)
   
    print*,is,ie

    do k=is(2),ie(2)-1
       do i=k+1,ie(1)
          a(i,k) = a(i,k)/a(k,k)
       end do
       do j=k+1,ie(2)
          do i=k+1,ie(1)
             a(i,j) = a(i,j)-a(i,k)*a(k,j)
          end do
       end do
    end do
    
  end subroutine LU_FACTORIZE
  !>
  subroutine L_SOLVE(L,Y,B)
    implicit none
    real(kind=8),dimension(:,:),allocatable :: L
    real(kind=8),dimension(:)  ,allocatable :: Y,B
    integer :: i,j,k,n
    integer :: is(2),ie(2)
   
    
    is = lbound(L)
    ie = ubound(L)

    do i=is(1),ie(1)
       do j=is(2),i-1
          b(i) = b(i)-l(i,j)*y(j)
       end do
       y(i)=b(i)
    end do
    
    
  end subroutine L_SOLVE

  subroutine SOLVE(A,X,B,DG)
    implicit none
    real(kind=8),dimension(:,:),allocatable :: A !> matrice de départ
    real(kind=8),dimension(:)  ,allocatable :: X,B,DG
    integer :: i,j,k,n
    integer :: is(2),ie(2)
   
    call LU_FACTORIZE(A)
    call L_SOLVE(A,DG,B)
    call U_SOLVE(A,X,DG)
    
  end subroutine SOLVE


  !>
  subroutine U_SOLVE(U,X,Y)
    implicit none
    real(kind=8),dimension(:,:),allocatable :: U
    real(kind=8),dimension(:)  ,allocatable :: X,Y
    integer :: i,j,k,n
    integer :: is(2),ie(2)
    
    is = lbound(U)
    ie = ubound(U)
   
    do i=ie(1),is(1),-1
       do j=i+1,ie(2)
          y(i) = y(i)-u(i,j)*x(j)
       end do
       x(i) = y(i)/u(i,i)
    end do
    
  end subroutine U_SOLVE


  subroutine L_INV(L,Y,B)
    implicit none
    real(kind=8),dimension(:,:),allocatable :: L
    real(kind=8),dimension(:,:),allocatable :: Y,B
    integer :: i,j,k,n
    integer :: is(2),ie(2)
    
    is = lbound(L)
    ie = ubound(L)
    
  end subroutine L_INV

  
  subroutine U_INV(U,X,Y)
    implicit none
    real(kind=8),dimension(:,:),allocatable :: U
    real(kind=8),dimension(:,:)  ,allocatable :: X,Y
    integer :: i,j,k,n
    integer :: is(2),ie(2)
    
    is = lbound(U)
    ie = ubound(U)
    
  
  end subroutine U_INV

end module module_lu_chps

