module gaussian_info
      implicit none
      real(8),allocatable :: samb(:,:,:),sambxy(:,:,:),dsamb(:,:,:)
contains
      include './include/bodyfd.f90'  
end module

