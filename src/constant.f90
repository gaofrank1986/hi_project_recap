module proj_cnst
     implicit none   
     real(8),parameter :: rsn(4,4) = reshape((/1.,  1.,  1.,  1., &
     &            1., -1.,  1., -1., 1.,  1., -1., -1., &
     &            1., -1., -1.,  1./),(/4,4/)) 
!
     real(8),parameter :: ex(4) =  (/  1.0d0,  1.0d0, -1.0d0, -1.0d0/)
     real(8),parameter :: ey(4) =  (/  1.0d0, -1.0d0, -1.0d0,  1.0d0/)
contains 
     subroutine print_rsn()
          implicit none
          integer i,j
          do i=1,4
                write (*,1000) (rsn(i,j),j=1,4)

          end do
          1000 format(4f10.6)

     end subroutine
     
end module

