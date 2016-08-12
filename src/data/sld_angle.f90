module info_mod
    use kinds
    implicit none

    type  info
        real,allocatable :: angle(:)
    contains
        procedure :: read_solid_angle
        procedure :: write_solid_angle
    end type
contains
    subroutine read_solid_angle(this,nnode)
        implicit none
        class(info) :: this
        integer :: nnode,tmp,i

        allocate(this%angle(nnode))
        open(99,file='./input/solid_angle.txt')
        do i=1,nnode
            read(99,*) tmp,this%angle(i)
        end do
        close(99)
    end subroutine

    subroutine write_solid_angle(this,nnode)
        implicit none
        class(info) :: this
        integer :: nnode,tmp,i

        open(99,file='./test_solid_angle.txt')
        do i=1,nnode
            write(99,'(i5,f14.8)') i,this%angle(i)
        end do
        close(99)
    end subroutine
end module

