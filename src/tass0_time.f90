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

    subroutine tassb0
        use mvar_mod
        use pvar_mod
                !use body_property
                use free_term,only:fterm,output_fterms
                use mfunc_mod
                use sebsm_mod

                implicit   none  
                integer  inode,ielem,j,jnode,ind,indd,ip
                integer ::    i,ii,is,l
                real(8)  xp,yp,zp,dpox,dpoy,dpoz,phi2
                real(8)  rsn(4,4)
                real(8)  bmatrix(4,8),amatrix(4,8),bmat(4)

                real(8)  s_angle
                !real(8) :: fterm_coef(0:3,4)
                real(8) :: dsign

                DATA RSN /1.,  1.,  1.,  1., &
             &            1., -1.,  1., -1., &
             &            1.,  1., -1., -1., &
             &            1., -1., -1.,  1./ 
        !
!`
!  ----------------------------------------------------
          WRITE(10, *)   ' IN TASSB0 '
         
          !DSDT(:)=0.0
!                 
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

        amata = 0.0d0
        cmata = 0.0d0

! =======================================================================
        call output_fterms()
        pause
        print *,"finished fterm output"
        do  500   inode=1,  nnf
                print *,inode
            xp=xyz(1,inode)
            yp=xyz(2,inode)
            zp=xyz(3,inode)

            !fterm_coef = 0
            call solidangle(inode,nnode,nelem,ncn,ncon,nodqua,&
             &                        h,xyz,dxyze,s_angle)    

            angle(inode)= s_angle
            amata(inode,inode,1:nsys)= s_angle!angle(inode)
            !  ---------------------------
            !  Integration on the free surface
                
            do   ielem=1,  nelemf
                call comp_link(ielem,inode,ii)
                if (ii .eq. 0)   then! if src node not on element 
                    call norm_ele1(ielem,xp,yp,zp,amatrix,bmatrix)

                else if (ii .ne. 0)   then 
                    call sing_ele1(inode,ielem,nodqua(inode),xp,yp,zp,&
                        &                   amatrix,bmatrix)
                end if 
                call common_block(0,0,ielem,inode,amatrix,bmatrix)!,fterm_coef)
            end do
            !  Integration on the body surface
            do    ielem=nelemf+1,  nelem

                call comp_link(ielem,inode,ii)
                if (ii .eq. 0)   then 
                    call norm_ele1(ielem,xp,yp,zp,amatrix,bmatrix)
                else if (ii .ne. 0)   then 
                    call sing_ele1(inode,ielem,nodqua(inode),xp,yp,zp,&
                     &                   amatrix,bmatrix)
                end if
                call common_block(1,0,ielem,inode,amatrix,bmatrix)

            end do
            
            do ip = 1,nsys
                    fra3(inode,ip) = fterm(inode,ip,1)!
                    frc31(inode,ip)=fterm(inode,ip,2)-fterm(inode,ip,1)*xp!
                    frc32(inode,ip)=fterm(inode,ip,3)-fterm(inode,ip,1)*yp!
                    !frc33(inode,ip)=fterm(inode,ip,4)-fterm(inode,ip,0)*zp!
            end do

            !TODO only works for ip=1
            write(404,5001) xp,yp,fterm(inode,1,1:4),angle(inode)
            write(405,5000) inode,fra3(inode,1),frc31(inode,1),frc32(inode,1)
            5000 format(I6,3f14.6)
            5001 format(7f14.8)
            


        500     continue
!
        
! =======================================================================
!    Source point is on the body surface
!
        do  1000   inode=nnf+1, nnode   

            xp=xyz(1,inode)
            yp=xyz(2,inode)
            zp=xyz(3,inode) 
            
            !fterm_coef = 0
            call solidangle(inode,nnode,nelem,ncn,ncon,nodqua,&
             &                    h,xyz,dxyze,s_angle) 

            angle(inode)=s_angle

            amata(inode,inode,1:nsys)= s_angle! angle(inode)

            do   ielem=1,  nelemf

                ii=0   
                call norm_ele0(ielem,xp,yp,zp,amatrix,bmatrix)
                call common_block(0,1,ielem,inode,amatrix,bmatrix)

            end do

            do  ielem=1+nelemf, nelem

                call comp_link(ielem,inode,ii)!
                if (ii .eq. 0)   then 
                    call norm_ele0(ielem,xp,yp,zp,amatrix,bmatrix)
                else if (ii .ne. 0)   then 
                    call sing_ele0(inode,ielem,nodqua(inode),xp,yp,zp,&
                     &                   amatrix,bmatrix)
                end if

                call common_block(1,1,ielem,inode,amatrix,bmatrix)
            end do
            write(404,5001) xp,yp,angle(inode)
1000     continue
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
        write(102, *) '  =========== before rludcmp =============='

  !      do i = 1,nnode
            !do j = 1,nnode
                !write(400,*) amata(i,j,1:nsys)
        !end do;end do
        !stop
        write(*,*) "Beginning LHS matrix inversion............"
        do ip=1, nsys
            call rludcmp(ip,amata,nnode,nnode,nsys,indx,dsign)  
        enddo
        write(*,*) "LHS matrix inversion  finished............"
        
        write(102, *) 
        write(102, *)
        write(102, *) '  =========== after rludcmp =============='

      return
      end

