module proj_cnst
     implicit none   
     real(8),parameter :: rsn(4,4) = reshape((/1.,  1.,  1.,  1., &
     &            1., -1.,  1., -1., 1.,  1., -1., -1., &
     &            1., -1., -1.,  1./),(/4,4/)) 
!
     real(8),parameter :: ex(4) =  (/  1.0d0,  1.0d0, -1.0d0, -1.0d0/)
     real(8),parameter :: ey(4) =  (/  1.0d0, -1.0d0, -1.0d0,  1.0d0/)
     real(8),parameter :: xiqsi(8) = (/-1.0d0, 0.0d0, 1.0d0, 1.0d0, 1.0d0, 0.0d0,-1.0d0,-1.0d0/)
     real(8),parameter :: xiqet(8) = (/-1.0d0,-1.0d0,-1.0d0, 0.0d0, 1.0d0, 1.0d0, 1.0d0, 0.0d0/)    
     real(8),parameter :: new8(8) = (/1, 5, 9, 13, 17, 21, 25, 29/)
     real(8),parameter :: pi4 = 12.56637061435917295385d0
contains 
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

         real(8),intent(in) :: a(3),b(3)
         real(8),dimension(3) :: c

         !c = spread(a,dim=2,ncopies=size(b)) * spread(b,dim=1,ncopies=size(a))
         c(1) = a(2) * b(3) - a(3) * b(2)
         c(2) = a(3) * b(1) - a(1) * b(3)
         c(3) = a(1) * b(2) - a(2) * b(1)
     end function

end module

