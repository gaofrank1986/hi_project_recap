    !> @file
    !! The green_funcs module provides several routines to calculate Green Functions 
    !! and their derivatives.
    !<
    !> @defgroup GreenFunc Lib
    !!
    module green_funcs 
        implicit none
        real(8),parameter :: pi = 3.14159265358979
    contains
        function GFunc(p,p0) result(ans)
            implicit none
            real(8),dimension(3) :: p,p0
            real(8) :: r,ans
            r = norm2(p-p0)
            ans = -1/(4*pi)*1/r
        end function

        function DGFunc(p,p0) result(ans)
            implicit none
            real(8),dimension(3),intent(in) :: p,p0
            real(8),dimension(3) :: ans

            real(8) :: r
            r = norm2(p-p0)
            ans = 1/(4*pi)*(p-p0)/r**3
        end function

        function Dy3GFunc(p,p0) result(ans)
            implicit none
            real(8),dimension(3),intent(in) :: p,p0
            real(8) :: ans,tmp(3)

            tmp = DGFunc(p0,p)
            ans=tmp(3)
        end function

        function Dy3DGFunc(p,p0) result(ans)
            implicit none
            real(8),dimension(3),intent(in) :: p,p0
            real(8),dimension(3) :: ans,dp

            real(8) :: r

            r = norm2(p-p0)
            dp = p-p0
            ans(1:2) = 1/(4*pi)*1/r**5*(3*dp(1:2)*dp(3))
            ans(3) = 1/(4*pi)*1/r**3*(3*dp(3)**2/r**2-1)
        end function

        !> get mirror point position
        function mirror(h,p) 
            implicit none
            real(8),intent(in) :: h,p(3)
            real(8),dimension(3) :: mirror 
            mirror=p
            mirror(3) = -(2*h+p(3))
        end function

        !> used for normal boudary integral equation 
        subroutine Gcombo0(h,p,p0,gxf)
            real(8),intent(in) :: h,p(3),p0(3)
            real(8),intent(out) :: gxf(4)

            !added mirrored src point
            gxf(1) = GFunc(p,p0)+GFunc(p,mirror(h,p0))
            gxf(2:4) = DGFunc(p,p0)+DGFunc(p,mirror(h,p0))

        end subroutine

        !> used for hypersingular boudary integral equation
        subroutine Gcombo1(h,p,p0,gxf)
            real(8),intent(in) :: h,p(3),p0(3)
            real(8),intent(out) :: gxf(4)

            !add mirrored sink point
            gxf(1) = Dy3GFunc(p,p0)-Dy3GFunc(p,mirror(h,p0))
            gxf(2:4) = Dy3DGFunc(p,p0)-Dy3DGFunc(p,mirror(h,p0))

        end subroutine

    end module               




