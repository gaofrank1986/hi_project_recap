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

      function mirror(h,p) !result(pz)
         implicit none
         real(8),intent(in) :: h,p(3)
         real(8),dimension(3) :: mirror 
         mirror=p
         mirror(3) = -(2*h+p(3))
      end function
         
      subroutine test_square()
           implicit none 
           integer :: p(3)

           p =(/1,2,3/)
           print *,p**2
      end subroutine
end module               

               


