    !C  TASSB0
    !C *********************************************************
    !C *                                                       *
    !C * Calculate the element contribution and assembly the   *
    !C * coefficients of the corresponding system of equationn *
    !C *                                                       *
    !C *********************************************************

    subroutine tassb0
        use free_term,only:fterm,output_fterms,calc_fterms
        use mesh,only:is_connected,topology_analysis
        use mesh,only:xyz,nnf,nnode,nelem,ncn,ncon,nodqua,nsys,nelemf,dxyze
        use wave,only:h
        use matrix_mod
        use misc_var,only:angle
        use linalg,only:rludcmp
        use proj_cnst,only:rsn
        use io

        implicit   none  
        integer  inode,ielem,ip
        integer ::    ii,hi
        real(8)  xp,yp,zp
        real(8)  bmatrix(4,8),amatrix(4,8)

        real(8)  s_angle
        !real(8) :: fterm_coef(0:3,4)
        real(8) :: dsign
        type(Ostream) :: fstream 
        fstream = Ostream("tass0",6)



        call fstream%fout('lhs matrix computing start......')
        

        !DSDT(:)=0.0
        !                 


        amata = 0.0d0
        cmata = 0.0d0

        ! =======================================================================
        !call output_fterms()
        pause

        do     inode=1, nnf
            !do     inode=874,877
            !print *,inode
            call fstream%fout(fstream%toString(inode))
            if (inode <= nnf) then
                hi = 2
            end if

            xp=xyz(1,inode)
            yp=xyz(2,inode)
            zp=xyz(3,inode)

            !fterm_coef = 0
            call solidangle(inode,nnode,nelem,ncn,ncon,nodqua,&
                &                        h,xyz,dxyze,s_angle)    

            angle(inode)=1.0d0- s_angle
            !cmata(inode,inode,1:nsys)= -1.0d0+s_angle!angle(inode)
            !---up---for solid angle term in rhs

            amata(inode,inode,1:nsys)= 1.0d0-s_angle!angle(inode)
            !  ---------------------------
            !  Integration on the free surface

            do   ielem=1,  nelemf
                if (not(is_connected(ielem,inode)))   then! if src node not on element 
                    call norm_elem_wrapper(ielem,xp,yp,zp,amatrix,bmatrix,hi)

                else 
                    call sing_elem_wrapper(inode,ielem,nodqua(inode),xp,yp,zp,&
                        &                   amatrix,bmatrix,hi)
                end if 
                call common_block(0,hi,ielem,inode,amatrix,bmatrix)!,fterm_coef)
            end do
            !  Integration on the body surface
            do    ielem=nelemf+1,  nelem

                if (not(is_connected(ielem,inode)))   then 
                    call norm_elem_wrapper(ielem,xp,yp,zp,amatrix,bmatrix,hi)
                else 
                    call sing_elem_wrapper(inode,ielem,nodqua(inode),xp,yp,zp,&
                        &                   amatrix,bmatrix,hi)
                end if
                call common_block(1,hi,ielem,inode,amatrix,bmatrix)

            end do

            !do ip = 1,nsys
            !fra3(inode,ip) = fterm(inode,ip,1)!
            !frc31(inode,ip)=fterm(inode,ip,2)-fterm(inode,ip,1)*xp!
            !frc32(inode,ip)=fterm(inode,ip,3)-fterm(inode,ip,1)*yp!
            !frc33(inode,ip)=fterm(inode,ip,4)-fterm(inode,ip,1)*zp!
            !end do
            do ip = 1,nsys
                fterm(inode,ip,2:4) = fterm(inode,ip,2:4)-fterm(inode,ip,1)*xyz(1:3,inode)
            end do

            !write(*,'(i6,5f14.6)') inode,fterm(inode,1,1:4),angle(inode)
            !write(398,'(i6,5f14.6)') inode,fterm(inode,1,1:4),angle(inode)
            !write(404,5001) xp,yp,fterm(inode,1,1:4),angle(inode)
            !write(405,5000) inode,fra3(inode,1),frc31(inode,1),frc32(inode,1)
            5000 format(I6,3f14.6)
            5001 format(7f14.8)



        enddo



        !do i=1,nnf
        !write(*,'(i6,4f14.6)') i,fterm(i,1,1:4)
        !enddo
        !call calc_fterms()

        ! =======================================================================
        !    Source point is on the body surface
        !
        do inode=nnf+1, nnode   

            hi=1

            xp=xyz(1,inode)
            yp=xyz(2,inode)
            zp=xyz(3,inode) 

            call solidangle(inode,nnode,nelem,ncn,ncon,nodqua,&
                &                    h,xyz,dxyze,s_angle) 

            angle(inode)=1.0d0 - s_angle

            amata(inode,inode,1:nsys)= 1.0d0-s_angle! angle(inode)

            do   ielem=1,  nelemf

                ii=0   
                call norm_elem_wrapper(ielem,xp,yp,zp,amatrix,bmatrix,hi)
                call common_block(0,hi,ielem,inode,amatrix,bmatrix)

            end do
            do  ielem=1+nelemf, nelem

                if (not(is_connected(ielem,inode)))   then 
                    call norm_elem_wrapper(ielem,xp,yp,zp,amatrix,bmatrix,hi)
                else 
                    call sing_elem_wrapper(inode,ielem,nodqua(inode),xp,yp,zp,&
                        &                   amatrix,bmatrix,hi)
                end if

                call common_block(1,hi,ielem,inode,amatrix,bmatrix)
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
        !print *,"I am here!"
        !
        !write(102, *) '  =========== before rludcmp =============='

        !do i = 1,nnode
        !do j = 1,nnode
        !write(400,*) amata(i,j,1:nsys)
        !end do;end do

        !do i = 1,nnode
        !do j = 1,nnoded
        !write(402,*) cmata(i,j,1:nsys)
        !end do;end do
        !stop

        !      do i = 1,nnode
        !write(401,*) amata(i,i,1:nsys)
        !end do

        call fstream%fout('Begin inversing LHS matrix............')
        do ip=1, nsys
            call rludcmp(ip,amata,nnode,nnode,nsys,indx,dsign)  
        enddo
        call fstream%fout('Finished inversing LHS matrix............')

        !write(102, *) 
        !write(102, *)
        !write(102, *) '  =========== after rludcmp =============='

    end subroutine

    subroutine common_block(flag1,hi,ielem,inode,amatrix,bmatrix)!,fterm_coef)

        use matrix_mod
        use mesh,only:ncn,ncon,ncond,nsys,xyz,dxyz,nnf
        use free_term,only:fterm
        use proj_cnst,only :rsn,ex,ey
        use wave_funcs_simple,only:dinp1

        implicit none
        integer,intent(in) :: ielem,inode,flag1,hi
        real(8),intent(in) :: amatrix(4,8),bmatrix(4,8)

        integer :: jncon,jnrml,j,ip,is,i
        real(8) :: dpdn,phi
        real(8) :: jp(3),jnp(3),prefix(3),dpdx(3)


        !if (flag1.eq.0) then
        !write (399,'(2i5,8f10.6)') inode,ielem,bmatrix(1,:)
        !end if

        do     j=1,  ncn(ielem) 
            jncon=ncon(ielem,j)!jth node in ielem
            jnrml=ncond(ielem,j)  !jth nrml in ielem
            do     ip=1, nsys          

                prefix=(/ex(ip),ey(ip),1.0d0/)
                jp = prefix*xyz(1:3,jncon)
                jnp = prefix*dxyz(1:3,jnrml)

                if (flag1.eq.0) then
                    do   is=1, nsys    
                        amata(inode,jncon,ip)=amata(inode,jncon,ip)+&
                            &              rsn(is,ip)*bmatrix(is,j)

                        cmata(inode,jnrml,ip)=cmata(inode,jnrml,ip)+rsn(is,ip)&
                            &           *amatrix(is,j)
                    enddo
                else 
                    do  is=1, nsys    
                        if(jncon .gt. nnf)  then
                            amata(inode,jncon,ip)=amata(inode,jncon,ip)-&
                                &                   rsn(is,ip)*amatrix(is,j)
                        else
                            cmata(inode,jncon,ip)=cmata(inode,jncon,ip)+rsn(is,ip)*amatrix(is,j)!*phi2
                        endif

                        cmata(inode,jnrml,ip)=cmata(inode,jnrml,ip)-rsn(is,ip)*bmatrix(is,j)!*dpdn
                    enddo
                end if
                
                ! indirect method for c3x,a3
                if (hi.eq.2) then!
                    do i = 1,4
                        call dinp1(i,jp,phi,dpdx)       
                        dpdn=dot_product(jnp,dpdx)
                        do    is=1, nsys             

                            fterm(inode,ip,i)=fterm(inode,ip,i)-rsn(is,ip)*bmatrix(is,j)*dpdn        
                            fterm(inode,ip,i) =fterm(inode,ip,i)+rsn(is,ip)*amatrix(is,j)*phi
                        enddo
                    end do

                end if
                !end if!fterm
            end do;end do !j,ip
        end subroutine
