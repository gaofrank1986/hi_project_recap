!
!  Ax = B problem
!  this program generating B matrix and solving for x
!
!C *********************************************************
!C *                                                       *
!C * Calculate the element contribution and assembly the   *
!C * coefficients of the corresponding system of equationn *
!C *                                                       *
!C *********************************************************
!C
    subroutine tassbt
  
      use data_all
      use motion
      use linalg,only:rlubksb
      use wave_funcs_simple,only:dinp,poxy
      use gradient,only:eval_gradient
      use free_term,only:fterm
      use proj_cnst,only:rsn,ex,ey

      implicit none  

      integer  inode,ind,is,ip,j,i,n1,n2
      real(8) xp,yp,zp,dpox,dpoy,dpoz,dpdn
      real(8) bmat(nsys),cmat(nnoded,nsys)
      real(8) :: dpoxyz_save(2,4,nnoded)

      print *,"Entering tassbt",timerk


      ! establish nrml to node relationship
      do i=1,nelem
          do j=1,ncn(i)
              n1 = ncon(i,j)
              n2 = ncond(i,j)
              nrml_2_node(n2) = n1
          enddo
      enddo

      dpoxyz_save = 0.0

      ip=1!fixme for initialize xp,yp,zp,symmetry not considered
      !// assign  dpoxyz with analytical vlaue
      do inode =1,nnf
          xp=ex(ip)*xyz(1,inode)
          yp=ey(ip)*xyz(2,inode)
          zp=       xyz(3,inode)

          !bkn(inode,1) =poxy(xp,yp,zp) 
          !fixme dpoxy use analytic value
          call dinp(xp,yp,zp,dpox,dpoy,dpoz)
          dpoxyz_save(1,1,inode) = dpox 
          dpoxyz_save(2,1,inode) = dpoy 
      end do

      !//get computed dpoxyz_save
      do  inode=1,  nnf  
!          do is=1, nsys; do ip=1, nsys
              !cmat(inode,is)=cmat(inode,is)+rsn(is,ip)*bkn(inode,ip)
          !enddo;enddo

          !fixme commented dpoxyz with computed value
          !
          !      do j=1,ncn(nodele(inode,1))
          !tmp(j) = bkn(ncon(nodele(inode,1),j),1)!get surface value
          !end do
          !call eval_gradient(inode,tmp,tmp2)
          !dpoxyz_save(1,1,inode) = tmp2(1,1)
          !dpoxyz_save(2,1,inode) = tmp2(2,1)
      end do

      !//======================================
      !  ----generate cmat-----------------
      !  ------------------------------

      cmat(:,:)=0.0d0

      !Surface boudary assignment | first part of cmat
      ! bk is the potential from time stepping----on surface node
      do  inode=1,  nnf  
          do is=1, nsys; do ip=1, nsys
              cmat(inode,is)=cmat(inode,is)+rsn(is,ip)*bkn(inode,ip)
          enddo;enddo
      enddo

      ! body surface BC assignment|  second part of cmat
      ! assign boudary value has dp/dn ---------------:
      do 40 inode=nnf+1, nnoded
          do 40 ip=1, nsys 

              !FIXME inode is wrong|corrected
              n2 = nrml_2_node(inode)
              xp=ex(ip)*xyz(1,n2)
              yp=ey(ip)*xyz(2,n2)
              zp=       xyz(3,n2)

              !timerk is read in in modlue wave-func-simple
              call dinp(xp,yp,zp,dpox,dpoy,dpoz)   !get initial condition    

              !================================================================
              !timerk at each rugga kutta time step
              !call dpoxyz(h,g,ampn,phi_w,beta,wkn,freq,timerk,rampf,xp,yp,zp,&
              !&              nfreq,nwave,iorder,dpox,dpoy,dpoz)

              dpdn=(dpox*ex(ip)*dxyz(1,inode)+&
                  &           dpoy*ey(ip)*dxyz(2,inode)+&
                  &           dpoz       *dxyz(3,inode) )
              !       
              ! consider body velocity contribution to dpdn,commented for
              !dpdn=dpdn-dsdt(1)*ex(ip)*dxyz(1,inode)
              !dpdn=dpdn-dsdt(2)*ey(ip)*dxyz(2,inode)
              !dpdn=dpdn-dsdt(3)*       dxyz(3,inode)
              !dpdn=dpdn-dsdt(4)*ey(ip)*dxyz(4,inode)
              !dpdn=dpdn-dsdt(5)*ex(ip)*dxyz(5,inode)
              !dpdn=dpdn-dsdt(6)*ex(ip)*ey(ip)*dxyz(6,inode)

              ! cmat[nff+1:nnoded] store dpdn
              do  is=1, nsys 
                  cmat(inode,is)=cmat(inode,is)+rsn(is,ip)*dpdn
              enddo

        40      continue   


        ! ==============================================
        !
        ! Generate bmata matrix---------
        !|||----------------------------------

        bmata(:,:)=0.0d0 !//initialize bmata

        do is=1, nsys; do ind=1,   nnode
            !---------------loop-body--------------------------
            bmata(ind,is)=dot_product(cmata(ind,:,is),cmat(:,is))
            ! .. cmat d\phi/dp   ... phi
            ! .. cmata
            if (ind<=nnf) then!potential only
                bmata(ind,is) = bmata(ind,is)-fterm(ind,is,1)*cmat(ind,is)&
                    &-fterm(ind,is,2)*dpoxyz_save(1,is,ind) &
                    &-fterm(ind,is,3)*dpoxyz_save(2,is,ind)
            end if
            !---------------------------------------------------
        enddo;enddo 

        !Solving Ax=B      // B is time varying
        ! result saved in bmata
        do  is=1, nsys   
            call rlubksb(is,amata,nnode,nnode, 1,nsys, 1,indx,bmata)
        enddo


        !do i = 1,nnode
        !write(401,*) bmata(i,1:nsys)
        !end do
        !do i=1,nnoded
        !!print *,cmat(i,1)
        !write(403,*) cmat(i,1)
        !enddo

        ! applied symmetry to get bmat
        do ip=1, nsys; do ind=1, nnode 
            bmat(ip)=(0.0d0, 0.0d0)
            do  is=1, nsys
                bmat(ip)=bmat(ip)+bmata(ind,is)*rsn(ip,is)
            enddo
            bmat(ip)=bmat(ip)/nsys
            unkn(ind,ip)=bmat(ip)        
        enddo;enddo

    end subroutine
