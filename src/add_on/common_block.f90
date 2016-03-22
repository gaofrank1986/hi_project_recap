    subroutine common_block(flag1,ielem,inode,amatrix,bmatrix,is_ft)!,fterm_coef)
        use mvar_mod,only :ncn,ncon,ncond,nsys,nnf,amata,cmata,xyz,dxyz
        use mfunc_mod,only: dinp0
        use free_term,only:fterm
        use proj_cnst,only:rsn,ex,ey

        implicit none
        integer,intent(in) :: ielem,inode,flag1
        real(8),intent(in) :: amatrix(4,8),bmatrix(4,8)
        integer :: is_ft

        integer :: jncon,jnrml,j,ip,is,i
        real(8) :: dpdn,phi,prefix(3),np_j(3),p_j(3)
        real(8) :: dpdx(3)
        
        do     j=1,  ncn(ielem) 
            jncon=ncon(ielem,j)!jth node in ielem
            jnrml=ncond(ielem,j)  !jth nrml in ielem
            do     ip=1, nsys          
                prefix=(/ex(ip),ey(ip),1.0d0/)
                p_j=dot_product(prefix,xyz(:,jncon))
                np_j=dot_product(prefix,dxyz(1:3,jnrml))

                if (flag1.eq.0) then
                    do   is=1, nsys    
                        amata(inode,jncon,ip)=amata(inode,jncon,ip)+rsn(is,ip)*bmatrix(is,j)

                        cmata(inode,jnrml,ip)=cmata(inode,jnrml,ip)+rsn(is,ip)*amatrix(is,j)
                    enddo
                else 
                    do  is=1, nsys    
                        !waterline node need special treatment
                        if(jncon .gt. nnf)  then
                            amata(inode,jncon,ip)=amata(inode,jncon,ip)-rsn(is,ip)*amatrix(is,j)
                        else
                            cmata(inode,jncon,ip)=cmata(inode,jncon,ip)+rsn(is,ip)*amatrix(is,j)!*phi2
                        endif
                        cmata(inode,jnrml,ip)=cmata(inode,jnrml,ip)-rsn(is,ip)*bmatrix(is,j)!*dpdn
                    enddo
                end if

                if (is_ft.eq.1) then!
                    do i = 1,4
                        !call dinp0(i,xsb,ysb,zsb,phi,dpox,dpoy,dpoz)       
                        !dpdn=dpox*nx+dpoy*ny+dpoz*nz
                        call dinp0(i,p_j,phi,dpdx)
                        dpdn=dot_product(dpdx,np_j)
                        do    is=1, nsys             
                            !fterm_coef(i,ip)=fterm_coef(i,ip)-rsn(is,ip)*bmatrix(is,j)*dpdn        
                            !fterm_coef(i,ip)=fterm_coef(i,ip)+rsn(is,ip)*amatrix(is,j)*phi
                            fterm(inode,ip,i)=fterm(inode,ip,i)-rsn(is,ip)*bmatrix(is,j)*dpdn        
                            fterm(inode,ip,i)=fterm(inode,ip,i)+rsn(is,ip)*amatrix(is,j)*phi
                        enddo
                    end do
                else 
                    !call dinp0(0,xsb,ysb,zsb,phi,dpox,dpoy,dpoz)      
                    !phi = 1.0d0
                    !do    is=1, nsys             
                    !!fterm_coef(0,ip) = fterm_coef(0,ip)-rsn(is,ip)*bmatrix(is,j)*dpdn
                    !fterm_coef(0,ip) = fterm_coef(0,ip)+rsn(is,ip)*amatrix(is,j)*phi
                    !end do
 
                end if
                end do;end do !j,ip

     end subroutine
    !subroutine common_block(flag1,flag2,ielem,inode,amatrix,bmatrix)!,fterm_coef)
        !use mvar_mod
        !use mfunc_mod 

        !implicit none
        !integer,intent(in) :: ielem,inode,flag1,flag2
        !real(8),intent(in) :: amatrix(4,8),bmatrix(4,8)

        !integer :: jncon,jnrml,j,ip,is,i
        !real(8) :: xsb,ysb,zsb,nx,ny,nz,dpdn,phi
        !real(8) :: dpox,dpoy,dpoz
        !real*8 rsn(4,4),ex(4),ey(4)
        
        !data rsn /1.,  1.,  1.,  1., &
     !&            1., -1.,  1., -1.,&
     !&            1.,  1., -1., -1., &
     !&            1., -1., -1.,  1./ 
!!
        !data ex /  1.0d0,  1.0d0, -1.0d0, -1.0d0/
        !data ey /  1.0d0, -1.0d0, -1.0d0,  1.0d0/
        
        
        !do     j=1,  ncn(ielem) 
            !jncon=ncon(ielem,j)!jth node in ielem
            !jnrml=ncond(ielem,j)  !jth nrml in ielem
            !do     ip=1, nsys          
                !!xsb=ex(ip)*xyz(1,jncon)
                !!ysb=ey(ip)*xyz(2,jncon)
                !!zsb=       xyz(3,jncon)
                !!nx=ex(ip)*dxyz(1,jnrml)
                !!ny=ey(ip)*dxyz(2,jnrml)
                !!nz=       dxyz(3,jnrml)
                !!TODO ??? change dinp for time domain??
                !!call dinp(xsb,ysb,zsb,dpox,dpoy,dpoz)   !get initial condition    
                !!dpdn=dpox*nx+dpoy*ny+dpoz*nz !get initial condition

                !if (flag1.eq.0) then
                    !do   is=1, nsys    
                    !amata(inode,jncon,ip)=amata(inode,jncon,ip)+&
                        !&              rsn(is,ip)*bmatrix(is,j)

                    !cmata(inode,jnrml,ip)=cmata(inode,jnrml,ip)+rsn(is,ip)&
                        !&           *amatrix(is,j)
                    !enddo
                !else 
                    !do  is=1, nsys    
                    !if(jncon .gt. nnf)  then
                        !amata(inode,jncon,ip)=amata(inode,jncon,ip)-&
                            !&                   rsn(is,ip)*amatrix(is,j)
                        !!cmata(inode,jnrml,ip)=cmata(inode,jnrml,ip)-rsn(is,ip)*bmatrix(is,j)!*dpdn
                    !else
                        !cmata(inode,jncon,ip)=cmata(inode,jncon,ip)+rsn(is,ip)*amatrix(is,j)!*phi2
                    !endif

                    !cmata(inode,jnrml,ip)=cmata(inode,jnrml,ip)-rsn(is,ip)*bmatrix(is,j)!*dpdn
                    !enddo
                !end if

                !!if (present(fterm_coef)) then
                        !!print *,"run with optional"
                !!if (flag2.eq.0) then!
                    !!do i = 0,3
                    !!call dinp0(i,xsb,ysb,zsb,phi,dpox,dpoy,dpoz)       
                    !!dpdn=dpox*nx+dpoy*ny+dpoz*nz
                    !!do    is=1, nsys             
                    !!fterm_coef(i,ip)=fterm_coef(i,ip)-rsn(is,ip)*bmatrix(is,j)*dpdn        
                    !!fterm_coef(i,ip)=fterm_coef(i,ip)+rsn(is,ip)*amatrix(is,j)*phi
                    !!enddo
                    !!end do
                !!else 
                    !!!call dinp0(0,xsb,ysb,zsb,phi,dpox,dpoy,dpoz)      
                    !!phi = 1.0d0
                    !!do    is=1, nsys             
                    !!!fterm_coef(0,ip) = fterm_coef(0,ip)-rsn(is,ip)*bmatrix(is,j)*dpdn
                    !!fterm_coef(0,ip) = fterm_coef(0,ip)+rsn(is,ip)*amatrix(is,j)*phi
                    !!end do
 
                !!end if
                !!end if!fterm
                !end do;end do !j,ip
     !end subroutine
