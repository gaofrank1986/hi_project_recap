module misc_var
    use kinds
    implicit none

    real*8,allocatable:: angle(:)
    real(8),allocatable ::fra3(:,:),frc31(:,:)
    real(8),allocatable :: frc32(:,:),frc33(:,:)
    integer,allocatable :: nrml_2_node(:)
contains
    subroutine init_misc_var(nsys,nnode,nnrml)
        implicit none
        integer,intent(in) ::nsys,nnode,nnrml
        if(allocated(angle)) deallocate(angle) 
        ! maybe a error if frc31-33 not allocated
        if(allocated(fra3)) deallocate(fra3,frc31,frc32,frc33) 
        if(allocated(angle)) deallocate(angle) 
        if(allocated(nrml_2_node)) deallocate(nrml_2_node) 

        allocate(angle(nnode))
        allocate(fra3(nnode,nsys),frc31(nnode,nsys),frc32(nnode,nsys),frc33(nnode,nsys))

        allocate(nrml_2_node(nnrml))
     end subroutine
end module
