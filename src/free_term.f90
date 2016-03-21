module free_term
     real(8),allocatable,public,save :: fterm(:,:,:)
     integer,private :: nsys,nnode
contains
     subroutine init_ft(a,b)
          integer,intent(in) :: a,b
          nsys = a
          nnode = b
          allocate(fterm(nnode,nsys,1:4))
          fterm=0.0d0
     end subroutine
     
     subroutine get_free_term()
          implicit none
          integer :: m,i,j,k
          real(8)::tmp
          !nsys = 1
          !nnode = 5
          !print *,'nsys=',nsys,'   ,  nnode=',nnode
          open(1001,file = './INPUT/fterm_final.txt',status='old')
          do j = 1,nnode
               do i=1,nsys
                    !read (1,*) (fterm(j,i,k) k=1,4)
                    !FIXME
                    read(1001,*) (fterm(j,i,k),k=1,4),tmp
                    !print *,(fterm(j,i,k),k=1,4)
                    !pause

          enddo;enddo
          close(1001)
      end subroutine

      subroutine output_fterms()
           implicit none
           integer :: j,i,k
           open(1000,file='./ft_reout.txt',status='unknown')
           do j =1,nnode
               do i=1,nsys
               write (1000,1000) (fterm(j,i,k),k=1,4)
           enddo;enddo
           1000  format(4f14.8)
           close(1000)
      end subroutine

           
end module
          
