module matrix_mod
    use kinds
    implicit none
    ! @var : [amata] lhs matrix in Ax = B
    ! @var : [bmata] rhs matrix in Ax = B
    ! @var : [cmata] used for generating bmata
    real(rk),allocatable :: amata(:,:,:),bmata(:,:)
    real(rk),allocatable :: cmata(:,:,:)
    ! @var : [indx] store pivot info in LU decomposision
    integer, allocatable :: indx(:,:)
    !
    !
    !real(rk),allocatable ::  lhs(:,:,:),rhs(:,:),rhs_c(:,:,:)
contains
   subroutine init_matrix_var(nsys,nnode,nnrml)
       implicit none
       integer,intent(in) :: nsys,nnode,nnrml
       if (allocated(amata)) deallocate(amata)
       if (allocated(bmata)) deallocate(bmata)
       if (allocated(cmata)) deallocate(cmata)
       if (allocated(indx)) deallocate(indx)

       allocate(amata(nnode,nnode,nsys),&
           &    bmata(nnode,nsys), indx(nnode,nsys))
       allocate(cmata(nnode,nnrml,nsys))
       amata=0.0d0
       bmata=0.0d0
       cmata = 0.0d0

   end subroutine

end module
