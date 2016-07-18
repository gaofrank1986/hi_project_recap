!-------------------------------------------------------------------------------
! DUTWAV
!-------------------------------------------------------------------------------
!  MODULE: proj_cnst
!
!> @brief
!! <Some Constant Values used in the project>
!!
!! @author
!! Song Gao,DUT 
!!
!! @date
!! 07 Nov 2013
!! 
!! @note <A note here.>
!! <Or starting here...>
!
! REVISION HISTORY:
!
! 07 Mar 2015 - 
!
!-------------------------------------------------------------------------------
module proj_cnst
    use kinds
    implicit none   
    real(rk),parameter :: rsn(4,4) = reshape((/1.,  1.,  1.,  1., &
        &            1., -1.,  1., -1., 1.,  1., -1., -1., &
        &            1., -1., -1.,  1./),(/4,4/)) 
    !
    real(rk),parameter :: ex(4) =  (/  1.0d0,  1.0d0, -1.0d0, -1.0d0/)
    real(rk),parameter :: ey(4) =  (/  1.0d0, -1.0d0, -1.0d0,  1.0d0/)
    real(rk),parameter :: xiqsi(8) = (/-1.0d0, 0.0d0, 1.0d0, 1.0d0, 1.0d0, 0.0d0,-1.0d0,-1.0d0/)
    real(rk),parameter :: xiqet(8) = (/-1.0d0,-1.0d0,-1.0d0, 0.0d0, 1.0d0, 1.0d0, 1.0d0, 0.0d0/)    
    !used in sg0_1
    real(rk),parameter :: new8(8) = (/1, 5, 9, 13, 17, 21, 25, 29/)
    real(rk),parameter :: pi4 = 12.56637061435917295385d0
    ! @var : transform matrix for hi_intgl order
    integer,parameter  :: trf(8,8) = reshape(  &
                                    [   1,0,0,0,0,0,0,0, &
                                        0,0,0,0,1,0,0,0, &
                                        0,1,0,0,0,0,0,0, &
                                        0,0,0,0,0,1,0,0, &
                                        0,0,1,0,0,0,0,0, &
                                        0,0,0,0,0,0,1,0, &
                                        0,0,0,1,0,0,0,0, &
                                        0,0,0,0,0,0,0,1],[8,8])
contains 
   ! subroutine init_trf()
        !trf(:,:)=0
        !trf(1,1)=1
        !trf(2,3)=1
        !trf(3,5)=1
        !trf(4,7)=1
        !trf(5,2)=1
        !trf(6,4)=1
        !trf(7,6)=1
        !trf(8,8)=1
    !end subroutine

    subroutine print_rsn()
        implicit none
        integer i,j
        do i=1,4
            write (*,1000) (rsn(i,j),j=1,4)

        end do
        1000 format(4f10.6)

    end subroutine

    function cross_product(a,b) result(c)
        implicit none

        real(rk),intent(in) :: a(3),b(3)
        real(rk),dimension(3) :: c

        !c = spread(a,dim=2,ncopies=size(b)) * spread(b,dim=1,ncopies=size(a))
        c(1) = a(2) * b(3) - a(3) * b(2)
        c(2) = a(3) * b(1) - a(1) * b(3)
        c(3) = a(1) * b(2) - a(2) * b(1)
    end function

end module

