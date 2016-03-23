!C  TASSB0
!C *********************************************************
!C *                                                       *
!C * Calculate the element contribution and assembly the   *
!C * coefficients of the corresponding system of equationn *
!C *                                                       *
!C *********************************************************
    include './add_on/common_block.f90'
    subroutine comp_link(ielem,inode,ii) 
        use mvar_mod
        ! return how many nodes in this element 
        implicit none
        integer,intent(in) :: inode,ielem
        integer,intent(out) :: ii

        integer :: i

        ii = 0 
        do i=1, nodnoe(inode)!!list of linked element by inode (how manys times inode appear in element)
        if(ielem .eq. nodele(inode,i)) then!
          ii=ii+1
        endif
        enddo
    end subroutine 

    logical function is_connected(ielem,inode) 
        use mvar_mod
        ! return how many nodes in this element 
        implicit none
        integer,intent(in) :: inode,ielem
        !logical :: is_connected

        integer :: i,ii

        ii = 0 
        do i=1, nodnoe(inode)!!list of linked element by inode (how manys times inode appear in element)
        if(ielem .eq. nodele(inode,i)) then!
          ii=ii+1
        endif
        enddo
        if (ii==0) then
            is_connected=.false.
        else
            is_connected=.true.
        endif
    end function 

    subroutine topology_analysis()
        use mvar_mod
        implicit none
        integer :: inode,ielem,j,l
        do 50 inode=1, nnode 
        l=0
        do 40 ielem=1,  nelem
        do 30 j=1,      ncn(ielem)
        if(inode.eq.ncon(ielem,j)) then
        l=l+1
        nodele(inode,l)=ielem! elem num linked to inode
        nodelj(inode,l)=j!in node-linked-element, inode appear as j-th node 
        endif
30      continue
40      continue
              nodnoe(inode)=l !total number of links
!        below related to symmetry
        nodqua(inode)=0
        if( nsys .ge. 2) then
          if( abs(xyz(2,inode)).lt.1.0e-06 ) then
          nodqua(inode)=2
          end if
        end if
!
        if( nsys .eq. 4) then
          if( abs(xyz(1,inode)).lt.1.0e-06.and.&
     &        abs(xyz(2,inode)).lt.1.0e-06) then
           nodqua(inode)=5
          else if( abs(xyz(1,inode)).lt.1.0e-06 ) then
           nodqua(inode)=4
          endif
        end if
!
50      continue
        PRINT *,"topology analysis finished"
    end subroutine

    subroutine tassb0
        use mvar_mod
        use pvar_mod
        use body_property
        use free_term,only:fterm,output_fterms,calc_fterms
        use mfunc_mod
        use proj_cnst,only:rsn

        !use sebsm_mod

        implicit   none  
        integer  inode,ielem,j,jnode,ind,indd,ip
        integer ::    i,ii,is,l
        real(8)  xp,yp,zp,dpox,dpoy,dpoz,phi2
        real(8)  bmatrix(4,8),amatrix(4,8),bmat(4)

        real(8)  s_angle
        !real(8) :: fterm_coef(0:3,4)
        real(8) :: dsign

        interface 
            logical function is_connected(ielem,inode)
                integer,intent(in)::ielem,inode
            end function
        end interface

        call topology_analysis()

        !
!`
!  ----------------------------------------------------
          WRITE(10, *)   ' IN TASSB0 '
         
          !DSDT(:)=0.0
!                 


        amata = 0.0d0
        cmata = 0.0d0

! =======================================================================
        !call output_fterms()
        pause
        !print *,"finished fterm output"
        do     inode=1,  nnf
                !print *,inode
            xp=xyz(1,inode)
            yp=xyz(2,inode)
            zp=xyz(3,inode)

            !fterm_coef = 0
            call solidangle(inode,nnode,nelem,ncn,ncon,nodqua,&
             &                        h,xyz,dxyze,s_angle)    

            angle(inode)=1.0d0- s_angle
            !cmata(inode,inode,1:nsys)= -1.0d0+s_angle!angle(inode)
            amata(inode,inode,1:nsys)= 1.0d0-s_angle!angle(inode)
            !  ---------------------------
            !  Integration on the free surface
                
            do   ielem=1,  nelemf
                if (not(is_connected(ielem,inode)))   then! if src node not on element 
                    call norm_elem_wrapper(ielem,xp,yp,zp,amatrix,bmatrix,2)

                else 
                    call sing_elem_wrapper(inode,ielem,nodqua(inode),xp,yp,zp,&
                        &                   amatrix,bmatrix,2)
                end if 
                call common_block(0,0,ielem,inode,amatrix,bmatrix)!,fterm_coef)
            end do
            !  Integration on the body surface
            do    ielem=nelemf+1,  nelem

                if (not(is_connected(ielem,inode)))   then 
                    call norm_elem_wrapper(ielem,xp,yp,zp,amatrix,bmatrix,2)
                else 
                    call sing_elem_wrapper(inode,ielem,nodqua(inode),xp,yp,zp,&
                     &                   amatrix,bmatrix,2)
                end if
                call common_block(1,0,ielem,inode,amatrix,bmatrix)

            end do
            
            !do ip = 1,nsys
                    !fra3(inode,ip) = fterm(inode,ip,1)!
                    !frc31(inode,ip)=fterm(inode,ip,2)-fterm(inode,ip,1)*xp!
                    !frc32(inode,ip)=fterm(inode,ip,3)-fterm(inode,ip,1)*yp!
                    !frc33(inode,ip)=fterm(inode,ip,4)-fterm(inode,ip,1)*zp!
            !end do
            !do ip = 1,nsys
                !fterm(inode,ip,2:4) = fterm(inode,ip,2:4)-fterm(inode,ip,1)*xyz(1:3,inode)
            !end do

            !write(*,'(i6,4f14.6)') inode,fterm(inode,1,1:4)
            !write(404,5001) xp,yp,fterm(inode,1,1:4),angle(inode)
            !write(405,5000) inode,fra3(inode,1),frc31(inode,1),frc32(inode,1)
            5000 format(I6,3f14.6)
            5001 format(7f14.8)
            


        enddo
        

        !forall(i=1:nnf,j=1:nsys)
            !fterm(i,j,2:4)=fterm(i,j,2:4)-fterm(i,j,1)*xyz(1:3,i)
        !end forall

        !do i=1,nnf
            !write(*,'(i6,4f14.6)') i,fterm(i,1,1:4)
        !enddo
        call calc_fterms()
        ! =======================================================================
        !    Source point is on the body surface
        !
        do inode=nnf+1, nnode   

            xp=xyz(1,inode)
            yp=xyz(2,inode)
            zp=xyz(3,inode) 

            !fterm_coef = 0
            call solidangle(inode,nnode,nelem,ncn,ncon,nodqua,&
                &                    h,xyz,dxyze,s_angle) 

            angle(inode)=1.0d0 - s_angle

            amata(inode,inode,1:nsys)= 1.0d0-s_angle! angle(inode)

            do   ielem=1,  nelemf

                ii=0   
                call norm_elem_wrapper(ielem,xp,yp,zp,amatrix,bmatrix,1)
                call common_block(0,1,ielem,inode,amatrix,bmatrix)

            end do

            do  ielem=1+nelemf, nelem

                if (not(is_connected(ielem,inode)))   then 
                    call norm_elem_wrapper(ielem,xp,yp,zp,amatrix,bmatrix,1)
                else 
                    call sing_elem_wrapper(inode,ielem,nodqua(inode),xp,yp,zp,&
                        &                   amatrix,bmatrix,1)
                end if

                call common_block(1,1,ielem,inode,amatrix,bmatrix)
            end do
            !write(404,5001) xp,yp,angle(inode)
        enddo
!
! =============================================

!        if( nsys .eq. 2) then
            !do inode=1, nnf
            !if(nodqua(inode) .eq. 2) then
                !amata(inode,inode,2)=1.0e20         
            !endif
            !enddo
        !else if( nsys .eq. 4) then
            !do inode=1, nnf
            !if(nodqua(inode) .eq. 2) then
                !amata(inode,inode,2)=1.0e20
                !amata(inode,inode,4)=1.0e20            
            !else if(nodqua(inode) .eq. 4) then
                !amata(inode,inode,3)=1.0e20
                !amata(inode,inode,2)=1.0e20
            !else if(nodqua(inode) .eq. 5) then
                !amata(inode,inode,2)=1.0e20
                !amata(inode,inode,3)=1.0e20            
                !amata(inode,inode,4)=1.0e20            
            !endif
            !enddo
        !endif

!! =============================================
        print *,"I am here!"
!
        !write(102, *) '  =========== before rludcmp =============='

        do i = 1,nnode
            do j = 1,nnode
                write(400,*) amata(i,j,1:nsys)
        end do;end do

        do i = 1,nnode
            do j = 1,nnoded
                write(402,*) cmata(i,j,1:nsys)
        end do;end do
        stop

  !      do i = 1,nnode
                !write(401,*) amata(i,i,1:nsys)
        !end do
        write(*,*) "Begin inversing LHS matrix............"
        do ip=1, nsys
            call rludcmp(ip,amata,nnode,nnode,nsys,indx,dsign)  
        enddo
        write(*,*) "Finished inversing LHS matrix............"

                write(102, *) 
        write(102, *)
        write(102, *) '  =========== after rludcmp =============='

      return
      end

