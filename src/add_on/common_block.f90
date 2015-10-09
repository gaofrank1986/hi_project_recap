    subroutine common_block(flag1,ielem,inode,amatrix,bmatrix,fterm_coef)
        use mvar_mod
        use mfunc_mod 

        implicit none
        integer,intent(in) :: ielem,inode,flag1
        real(8),intent(in) :: amatrix(4,8),bmatrix(4,8)
        real(8),intent(inout) :: fterm_coef(0:3,4)

        integer :: jncon,kncon,j,ip,is,i
        real(8) :: xsb,ysb,zsb,nx,ny,nz,dpdn,phi
        real(8) :: dpox,dpoy,dpoz
          REAL*8 RSN(4,4),EX(4),EY(4)

        DATA RSN /1.,  1.,  1.,  1., &
     &            1., -1.,  1., -1.,&
     &            1.,  1., -1., -1., &
     &            1., -1., -1.,  1./ 
!
          DATA EX /  1.0d0,  1.0d0, -1.0d0, -1.0d0/                       
        DATA EY /  1.0d0, -1.0d0, -1.0d0,  1.0d0/

        do     j=1,  ncn(ielem) 
            jncon=ncon(ielem,j)
            kncon=ncond(ielem,j)  
        do     ip=1, nsys          
            xsb=ex(ip)*xyz(1,jncon)
            ysb=ey(ip)*xyz(2,jncon)
            zsb=       xyz(3,jncon)
            nx=ex(ip)*dxyz(1,kncon)
            ny=ey(ip)*dxyz(2,kncon)
            nz=       dxyz(3,kncon)
            call dinp(xsb,ysb,zsb,dpox,dpoy,dpoz)       
            dpdn=dpox*nx+dpoy*ny+dpoz*nz
        if (flag1.eq.0) then

            do   is=1, nsys    
                amata(inode,jncon,ip)=amata(inode,jncon,ip)+&
                    &              rsn(is,ip)*bmatrix(is,j)               

                bmata(inode,ip)=bmata(inode,ip)+rsn(is,ip)&
                            &*amatrix(is,j)*poxy(xsb,ysb,zsb)
            enddo
        else 
              do  is=1, nsys    
             if(jncon .gt. nnf)  then
                 amata(inode,jncon,ip)=amata(inode,jncon,ip)-&
         &                          rsn(is,ip)*amatrix(is,j)
             else
                 phi=poxy(xsb,ysb,zsb)
         bmata(inode,ip)=bmata(inode,ip)+rsn(is,ip)*amatrix(is,j)*phi   !  * ******
             endif

         bmata(inode,ip)=bmata(inode,ip)-rsn(is,ip)*bmatrix(is,j)*dpdn  !  * ******
             enddo
        end if
!
        do i = 0,3
            call dinp0(i,xsb,ysb,zsb,phi,dpox,dpoy,dpoz)       
            dpdn=dpox*nx+dpoy*ny+dpoz*nz
            do    is=1, nsys             
                fterm_coef(i,ip)=fterm_coef(i,ip)-rsn(is,ip)*bmatrix(is,j)*dpdn        
                fterm_coef(i,ip)=fterm_coef(i,ip)+rsn(is,ip)*amatrix(is,j)*phi
            enddo
        end do

        end do;end do !j,ip
     end subroutine
