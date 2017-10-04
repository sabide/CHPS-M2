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
    
  end subroutine L_SOLVE
  !>
  subroutine U_SOLVE(U,X,Y)
    implicit none
    real(kind=8),dimension(:,:),allocatable :: U
    real(kind=8),dimension(:)  ,allocatable :: X,Y
    integer :: i,j,k,n
    integer :: is(2),ie(2)
    
    is = lbound(U)
    ie = ubound(U)
    
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

