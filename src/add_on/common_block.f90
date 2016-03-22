    subroutine common_block(flag1,flag2,ielem,inode,amatrix,bmatrix)!,fterm_coef)
        use mvar_mod
        use mfunc_mod 
        use free_term,only:fterm

        implicit none
        integer,intent(in) :: ielem,inode,flag1,flag2
        real(8),intent(in) :: amatrix(4,8),bmatrix(4,8)

        integer :: jncon,jnrml,j,ip,is,i
        real(8) :: xsb,ysb,zsb,nx,ny,nz,dpdn,phi
        real(8) :: dpox,dpoy,dpoz,p(3),np(3),prefix(3),dpdx(3)
        real*8 rsn(4,4),ex(4),ey(4)
        
        data rsn /1.,  1.,  1.,  1., &
     &            1., -1.,  1., -1.,&
     &            1.,  1., -1., -1., &
     &            1., -1., -1.,  1./ 
!
        data ex /  1.0d0,  1.0d0, -1.0d0, -1.0d0/
        data ey /  1.0d0, -1.0d0, -1.0d0,  1.0d0/
        
        
        do     j=1,  ncn(ielem) 
            jncon=ncon(ielem,j)!jth node in ielem
            jnrml=ncond(ielem,j)  !jth nrml in ielem
            do     ip=1, nsys          
                xsb=ex(ip)*xyz(1,jncon)
                ysb=ey(ip)*xyz(2,jncon)
                zsb=       xyz(3,jncon)
                nx=ex(ip)*dxyz(1,jnrml)
                ny=ey(ip)*dxyz(2,jnrml)
                nz=       dxyz(3,jnrml)
                prefix=(/ex(ip),ey(ip),1.0d0/)
                !print *,shape(prefix)
                !print *,shape(xyz(1:3,jncon))
                !print *,shape(p)
                !print *,shape(prefix*xyz(1:3,jncon))
                !stop
                p = prefix*xyz(1:3,jncon)
                np = prefix*dxyz(1:3,jnrml)
                !print *,xsb,ysb,zsb
                !print *,p
                !print *,nx,ny,nz
                !print *,np
                !pause
                !TODO ??? change dinp for time domain??
                !call dinp(xsb,ysb,zsb,dpox,dpoy,dpoz)   !get initial condition    
                !dpdn=dpox*nx+dpoy*ny+dpoz*nz !get initial condition

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
                        !cmata(inode,jnrml,ip)=cmata(inode,jnrml,ip)-rsn(is,ip)*bmatrix(is,j)!*dpdn
                    else
                        cmata(inode,jncon,ip)=cmata(inode,jncon,ip)+rsn(is,ip)*amatrix(is,j)!*phi2
                    endif

                    cmata(inode,jnrml,ip)=cmata(inode,jnrml,ip)-rsn(is,ip)*bmatrix(is,j)!*dpdn
                    enddo
                end if

                !if (present(fterm_coef)) then
                        !print *,"run with optional"
                if (flag2.eq.0) then!
                    !do i = 0,3
                    do i = 1,4
                    !call dinp0(i,xsb,ysb,zsb,phi,dpox,dpoy,dpoz)       
                    call dinp1(i,p,phi,dpdx)       
                    !dpdn=dpox*nx+dpoy*ny+dpoz*nz
                    dpdn=dot_product(np,dpdx)
                    do    is=1, nsys             
                    !fterm(inode,ip,i+1)=fterm(inode,ip,i+1)-rsn(is,ip)*bmatrix(is,j)*dpdn        
                    !fterm(inode,ip,i+1) =fterm(inode,ip,i+1)+rsn(is,ip)*amatrix(is,j)*phi
                    fterm(inode,ip,i)=fterm(inode,ip,i)-rsn(is,ip)*bmatrix(is,j)*dpdn        
                    fterm(inode,ip,i) =fterm(inode,ip,i)+rsn(is,ip)*amatrix(is,j)*phi
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
                !end if!fterm
                end do;end do !j,ip
     end subroutine
