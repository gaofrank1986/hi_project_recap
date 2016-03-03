    subroutine common_block(flag1,flag2,ielem,inode,amatrix,bmatrix,fterm_coef)
        use mvar_mod
        use mfunc_mod 

        implicit none
        integer,intent(in) :: ielem,inode,flag1,flag2
        real(8),intent(in) :: amatrix(4,8),bmatrix(4,8)
        real(8),intent(inout) :: fterm_coef(0:3,4)

        integer :: jncon,cur_nrml,j,ip,is,i
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
            jncon=ncon(ielem,j)
            cur_nrml=ncond(ielem,j)  
            do     ip=1, nsys          
                xsb=ex(ip)*xyz(1,jncon)
                ysb=ey(ip)*xyz(2,jncon)
                zsb=       xyz(3,jncon)
                nx=ex(ip)*dxyz(1,cur_nrml)
                ny=ey(ip)*dxyz(2,cur_nrml)
                nz=       dxyz(3,cur_nrml)

                call dinp(xsb,ysb,zsb,dpox,dpoy,dpoz)   !get initial condition    
                dpdn=dpox*nx+dpoy*ny+dpoz*nz !get initial condition

                if (flag1.eq.0) then
                    do   is=1, nsys    
                    amata(inode,jncon,ip)=amata(inode,jncon,ip)+&
                        &              rsn(is,ip)*bmatrix(is,j)

                    bmata(inode,ip)=bmata(inode,ip)+rsn(is,ip)&
                        &           *amatrix(is,j)*poxy(xsb,ysb,zsb)
                    enddo
                else 
                    do  is=1, nsys    
                    !basically do the same job as above
                    !and the trick part is for fs node in body elem
                    !whcih should
                    if(jncon .gt. nnf)  then
                        amata(inode,jncon,ip)=amata(inode,jncon,ip)-&
                            &                   rsn(is,ip)*amatrix(is,j)
                    else
                        phi=poxy(xsb,ysb,zsb)! get intial conditon
                        bmata(inode,ip)=bmata(inode,ip)+rsn(is,ip)*amatrix(is,j)*phi
                    endif

                    bmata(inode,ip)=bmata(inode,ip)-rsn(is,ip)*bmatrix(is,j)*dpdn
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
                else 
                    !call dinp0(0,xsb,ysb,zsb,phi,dpox,dpoy,dpoz)      
                    phi = 1.0d0
                    do    is=1, nsys             
                    !fterm_coef(0,ip) = fterm_coef(0,ip)-rsn(is,ip)*bmatrix(is,j)*dpdn
                    fterm_coef(0,ip) = fterm_coef(0,ip)+rsn(is,ip)*amatrix(is,j)*phi
                    end do
 
                end if
                end do;end do !j,ip
     end subroutine
