module free_term
     real(8),allocatable,public,save :: fterm(:,:,:)
     integer,private,save :: nsys,nnf
contains
     subroutine init_ft(a,b)
          integer,intent(in) :: a,b
          nsys = a
          nnf = b
          allocate(fterm(nnf,nsys,1:4))
          fterm=0.0d0
     end subroutine

     subroutine get_free_term()
          implicit none
          integer :: m,i,j,k
          real(8)::tmp
          open(1001,file = './INPUT/fterm_final.txt',status='old')
          do j = 1,nnf
              do i=1,nsys
                  read(1001,*) (fterm(j,i,k),k=1,4),tmp
              enddo
          enddo
          close(1001)
      end subroutine

      subroutine output_fterms()
           implicit none
           integer :: j,i,k
           open(1000,file='./ft_reout.txt',status='unknown')
           do j =1,nnf
               do i=1,nsys
                   write (1000,'(4f14.8)') (fterm(j,i,k),k=1,4)
               enddo
           enddo
           close(1000)
      end subroutine

      subroutine calc_fterms()
          use mesh,only:xyz
          implicit none
          integer i,j

          forall(i=1:nnf,j=1:nsys)
              fterm(i,j,2:4)=fterm(i,j,2:4)-fterm(i,j,1)*xyz(1:3,i)
          end forall

          do i=1,nnf
              write(*,'(i6,4f14.6)') i,fterm(i,1,1:4)
          enddo
        end subroutine
end module

