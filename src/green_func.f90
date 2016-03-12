module green_funcs 
      implicit none
contains
      function GFunc(p,p0) result(ans)
         implicit none
         real(8),dimension(3) :: p,p0
         real(8),parameter :: pi=3.1415926
         real(8) :: r,ans
         r = norm2(p-p0)
         ans = -1/(4*pi)*1/r
      end function

      function DGFunc(p,p0) result(ans)
         implicit none
         real(8),dimension(3),intent(in) :: p,p0
         real(8),dimension(3) :: ans

         real(8),parameter :: pi=3.1415926
         real(8) :: r
         r = norm2(p-p0)
         ans = 1/(4*pi)*(p-p0)/r**3
      end function
      
      function DyDGFunc(p,p0) result(ans)
         implicit none
         real(8),dimension(3),intent(in) :: p,p0
         real(8),dimension(3) :: ans

         real(8),parameter :: pi=3.1415926
         real(8) :: r

         r = norm2(p-p0)
         ans = 1/(4*pi)*1/r**3*(3*(p-p0)**2/r**2-1)
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

               


