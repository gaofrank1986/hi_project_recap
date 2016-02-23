module free_term
     real(8),allocatable,protected :: fterm(:,:,:)
     integer,private :: nsys,nnode
contains
     subroutine init_ft(a,b)
          integer,intent(in) :: a,b
          nsys = a
          nnode = b
     end subroutine
     
     subroutine get_free_term()
          implicit none
          integer :: m,i,j,k
          !nsys = 1
          !nnode = 5
          print *,'nsys=',nsys,'   ,  nnode=',nnode
          allocate(fterm(nnode,nsys,1:4))
          open(1,file = './fterm_final.txt',status='old')
          do j = 1,nnode
               do i=1,nsys
                    !read (1,*) (fterm(j,i,k) k=1,4)
                    !FIXME
                    read(1,*) m,(fterm(m,i,k),k=1,4)

          enddo;enddo
          close(1)
      end subroutine

      subroutine output_fterms()
           implicit none
           integer :: j,i,k
           open(1,file='./ft_reout.txt',status='unknown')
           do j =1,nnode
               do i=1,nsys
               write (1,1000) (fterm(j,i,k),k=1,4)
           enddo;enddo
           1000  format(4f14.8)
           close(1)
      end subroutine

           
end module
          
