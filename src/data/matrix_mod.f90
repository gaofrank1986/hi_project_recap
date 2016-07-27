module matrix_mod
    use kinds
    implicit none
    integer, allocatable :: indx(:,:)
    real(rk),allocatable :: amata(:,:,:),bmata(:,:)
    real(rk),allocatable :: cmata(:,:,:)
    !real(rk),allocatable ::  lhs(:,:,:),rhs(:,:),rhs_c(:,:,:)
    ! new name lhs,rhs,rhs constructor
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
