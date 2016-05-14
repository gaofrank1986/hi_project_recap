!-------------------------------------------------------------------------------
! DUTWAV
!-------------------------------------------------------------------------------
!  MODULE: Green Function
!
!> @brief
!! <BriefDescription>
!!
!! @author
!! Song Gao,DUT 
!!
!! @date
!! 07 Mar 2016
!! 
!! @note <A note here.>
!! <Or starting here...>
!
! REVISION HISTORY:
!
! 07 Mar 2015 - Added gcombo to wrap mirrored source or sink. 
!
!-------------------------------------------------------------------------------
module green_funcs 
    use kinds
    implicit none
    real(rk),parameter :: pi = 3.14159265358979
contains

    !<  Green Function
    function GFunc(p,p0) result(ans)
        implicit none
        real(rk),dimension(3) :: p,p0
        real(rk) :: r,ans
        r = norm2(p-p0)
        ans = -1/(4*pi)*1/r
    end function

    !<  dG/dxi 
    !!
    function DGFunc(p,p0) result(ans)
        implicit none
        real(rk),dimension(3),intent(in) :: p,p0
        real(rk),dimension(3) :: ans

        real(rk) :: r
        r = norm2(p-p0)
        ans = 1/(4*pi)*(p-p0)/r**3
    end function

    !<  DG/dy3
    !!
    function Dy3GFunc(p,p0) result(ans)
        implicit none
        real(rk),dimension(3),intent(in) :: p,p0
        real(rk) :: ans,tmp(3)

        tmp = DGFunc(p0,p)
        ans=tmp(3)
    end function

    !< d^2{G}/{dy3 dxi}
    function Dy3DGFunc(p,p0) result(ans)
        implicit none
        real(rk),dimension(3),intent(in) :: p,p0
        real(rk),dimension(3) :: ans,dp

        real(rk) :: r

        r = norm2(p-p0)
        dp = p-p0
        ans(1:2) = 1/(4*pi)*1/r**5*(3*dp(1:2)*dp(3))
        ans(3) = 1/(4*pi)*1/r**3*(3*dp(3)**2/r**2-1)
    end function

    !> get mirror point position
    function mirror(h,p) 
        implicit none
        real(rk),intent(in) :: h,p(3)
        real(rk),dimension(3) :: mirror 
        mirror=p
        mirror(3) = -(2*h+p(3))
    end function

    !> ------------------------------------------
    !>             Combo Function
    !> ------------------------------------------
    !> used for normal boudary integral equation 
    !> -------------------------------------------
    subroutine Gcombo0(h,p,p0,gxf)
        real(rk),intent(in) :: h,p(3),p0(3)
        real(rk),intent(out) :: gxf(4)

        !added mirrored src point
        gxf(1) = GFunc(p,p0)+GFunc(p,mirror(h,p0))
        gxf(2:4) = DGFunc(p,p0)+DGFunc(p,mirror(h,p0))

    end subroutine
    
    !> -------------------------------------------------
    !> used for hypersingular boudary integral equation
    !> -------------------------------------------------
    subroutine Gcombo1(h,p,p0,gxf)
        real(rk),intent(in) :: h,p(3),p0(3)
        real(rk),intent(out) :: gxf(4)

        !add mirrored sink point
        gxf(1) = Dy3GFunc(p,p0)-Dy3GFunc(p,mirror(h,p0))
        gxf(2:4) = Dy3DGFunc(p,p0)-Dy3DGFunc(p,mirror(h,p0))

    end subroutine

    !> hypersingular BIE with only src on boudary surface
    subroutine Gcombo1_1(h,p,p0,gxf)
        real(rk),intent(in) :: h,p(3),p0(3)
        real(rk),intent(out) :: gxf(4)

        gxf(1) = Dy3GFunc(p,p0)
        gxf(2:4) = Dy3DGFunc(p,p0)

    end subroutine

    !> hypersingular BIE with only mirrored sink
    subroutine Gcombo1_2(h,p,p0,gxf)
        real(rk),intent(in) :: h,p(3),p0(3)
        real(rk),intent(out) :: gxf(4)

        gxf(1) = -Dy3GFunc(p,mirror(h,p0))
        gxf(2:4) = -Dy3DGFunc(p,mirror(h,p0))

    end subroutine
end module               




