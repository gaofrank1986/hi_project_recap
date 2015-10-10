    subroutine common_block(flag1,flag2,ielem,inode,amatrix,bmatrix,fterm_coef)
        use mvar_mod
        use mfunc_mod 

        implicit none
        integer,intent(in) :: ielem,inode,flag1,flag2
        real(8),intent(in) :: amatrix(4,8),bmatrix(4,8)
        real(8),intent(inout) :: fterm_coef(0:3,4)

        integer :: cur_node,cur_nrml,j,ip,is,i
        real(8) :: xsb,ysb,zsb,nx,ny,nz,dpdn,phi
        real(8) :: dpox,dpoy,dpoz
        real*8 rsn(4,4),ex(4),ey(4)

        data rsn /1.,  1.,  1.,  1., &
     &            1., -1.,  1., -1.,&
     &            1.,  1., -1., -1., &
     &            1., -1., -1.,  1./ 
!
        data ex /  1.0d0,  1.0d0, -1.0d0, -1.0d0/
        data ey /  1.0d0, -1.0d0, -1.0d0,  1.0d0/

        do     j=1,  ncn(ielem) 
            cur_node=ncon(ielem,j)
            cur_nrml=ncond(ielem,j)  
            do     ip=1, nsys          
                xsb=ex(ip)*xyz(1,cur_node)
                ysb=ey(ip)*xyz(2,cur_node)
                zsb=       xyz(3,cur_node)
                nx=ex(ip)*dxyz(1,cur_nrml)
                ny=ey(ip)*dxyz(2,cur_nrml)
                nz=       dxyz(3,cur_nrml)

                call dinp(xsb,ysb,zsb,dpox,dpoy,dpoz)       
                dpdn=dpox*nx+dpoy*ny+dpoz*nz

                if (flag1.eq.0) then
                    do   is=1, nsys    
                    amata(inode,cur_node,ip)=amata(inode,cur_node,ip)+&
                        &              rsn(is,ip)*bmatrix(is,j)               

                    bmata(inode,ip)=bmata(inode,ip)+rsn(is,ip)&
                        &           *amatrix(is,j)*poxy(xsb,ysb,zsb)
                    enddo
                else 
                    do  is=1, nsys    
                    if(cur_node .gt. nnf)  then
                        amata(inode,cur_node,ip)=amata(inode,cur_node,ip)-&
                            &                   rsn(is,ip)*amatrix(is,j)
                    else
                        phi=poxy(xsb,ysb,zsb)
                        bmata(inode,ip)=bmata(inode,ip)+rsn(is,ip)*amatrix(is,j)*phi   !  * ******
                    endif

                    bmata(inode,ip)=bmata(inode,ip)-rsn(is,ip)*bmatrix(is,j)*dpdn  !  * ******
                    enddo
                end if
                if (flag2.eq.0) then!
                    do i = 0,3
                    call dinp0(i,xsb,ysb,zsb,phi,dpox,dpoy,dpoz)       
                    dpdn=dpox*nx+dpoy*ny+dpoz*nz
                    do    is=1, nsys             
                    fterm_coef(i,ip)=fterm_coef(i,ip)-rsn(is,ip)*bmatrix(is,j)*dpdn        
                    fterm_coef(i,ip)=fterm_coef(i,ip)+rsn(is,ip)*amatrix(is,j)*phi
                    enddo
                    end do
                end if
                end do;end do !j,ip
     end subroutine
