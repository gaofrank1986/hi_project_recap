!  TASSBT+ TINBODT+ NSBDT
!
!C *********************************************************
!C *                                                       *
!C * Calculate the element contribution and assembly the   *
!C * coefficients of the corresponding system of equationn *
!C *                                                       *
!C *********************************************************
!C
        SUBROUTINE TASSBT
  
      USE MVAR_MOD
      USE PVAR_MOD
      !use wave2,only:dpoxyz
      use mfunc_mod,only:rlubksb,dinp,poxy,dinp0
      use gradient,only:eval_gradient

      implicit none  

      integer  inode,ind,is,ip,j,i,n1,n2
      real(8) xp,yp,zp,dpox,dpoy,dpoz,dpdn,nx,ny,nz
      real(8) rsn(4,4),ex(4),ey(4)
      real(8) bmat(nsys),cmat(nnoded,nsys)
      real(8) :: dpoxyz_save(2,4,nnoded),tmp(8),tmp2(2,1)

      data rsn /1.,  1.,  1.,  1.,& 
     &            1., -1.,  1., -1.,&
     &            1.,  1., -1., -1.,& 
     &            1., -1., -1.,  1./ 
      data ex / 1.,  1., -1., -1./                                     
      data ey / 1., -1., -1.,  1./

       print *,"Entering tassbt"
!
!        WRITE(102,*)
!        WRITE(102,*)   'T=', TimeRK                 

      ! establish nrml to node relationship
      do i=1,nelem
          do j=1,ncn(i)
              n1 = ncon(i,j)
              n2 = ncond(i,j)
              nrml_2_node(n2) = n1
      enddo;enddo

!C ***************************************
!C
      cmat(:,:)=0.0
      dpoxyz_save = 0.0
      !bkn=0
      !do inode =1,nnf
          !xp=ex(ip)*xyz(1,inode)
          !yp=ey(ip)*xyz(2,inode)
          !zp=       xyz(3,inode)

          !bkn(inode,1) =poxy(xp,yp,zp) 
          !call dinp(xp,yp,zp,dpox,dpoy,dpoz)
          !dpoxyz_save(1,1,inode) = dpox 
          !dpoxyz_save(2,1,inode) = dpoy 
      !end do

        ! assign boundary value has potentials--------
        ! bkn is the potential from time stepping----on surface node
        do  20   inode=1,  nnf  
            do is=1, nsys
                do ip=1, nsys
                  cmat(inode,is)=cmat(inode,is)+rsn(is,ip)*bkn(inode,ip)
            enddo;enddo
            
            do j=1,ncn(nodele(inode,1))
                 tmp(j) = bkn(ncon(nodele(inode,1),j),1)!get surface value
            end do
            call eval_gradient(inode,tmp,tmp2)
            dpoxyz_save(1,1,inode) = tmp2(1,1)
            dpoxyz_save(2,1,inode) = tmp2(2,1)
        
 
        20   continue
    1021 format(8f10.6)
    !print *,"cmat",cmat
!
!        WRITE(102,*)  'INODE,    CMAT(INODE,1) --20'
!        DO    INODE=1,  NNF
!        WRITE(102, 620)  INODE,CMAT(INODE,1)
!        ENDDO


! ==============================================
!
!         WRITE(102,*) 'DSDT=',DSDT

        ! assign boudary value has dp/dn ---------------:
        do 40 inode=nnf+1, nnoded
        do 40 ip=1, nsys 
        !FIXME inode is wrong
        n2 = nrml_2_node(inode)
        xp=ex(ip)*xyz(1,n2)
        yp=ey(ip)*xyz(2,n2)
        zp=       xyz(3,n2)
        call dinp(xp,yp,zp,dpox,dpoy,dpoz)   !get initial condition    
        !================================================================
        !timerk at each rugga kutta time step
        !call dpoxyz(h,g,ampn,phi_w,beta,wkn,freq,timerk,rampf,xp,yp,zp,&
     !&              nfreq,nwave,iorder,dpox,dpoy,dpoz)

        dpdn=(dpox*ex(ip)*dxyz(1,inode)+&
    &           dpoy*ey(ip)*dxyz(2,inode)+&
    &           dpoz       *dxyz(3,inode) )
!       
!         WRITE(102,*) 'DPDN=',DPDN,DPOX,DPOY,DPOZ
        ! consider body velocity contribution to dpdn,commented for
        ! testing 
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

             !this is only incident wave, missing the diffraction wave
40      continue   
!

!        WRITE(102,*)  'INODE,    CMAT(INODE,1) --40'
!        DO    INODE=NNF+1,  NNODED
!        WRITE(102, 620)  INODE,CMAT(INODE,1)
!        ENDDO

! ==============================================
!
        ! Generate bmata matrix---------
        !|||----------------------------------

        bmata(:,:)=0.0d0
        do 200 is=1, nsys
            do 200 ind=1,   nnode
                !do inode=1, nnoded
        !---------------loop-body--------------------------
        bmata(ind,is)=dot_product(cmata(ind,:,is),cmat(:,is))
                !end do
                if (ind<=nnf) then!potential only
                        bmata(ind,is) = bmata(ind,is)-fra3(ind,is)*cmat(ind,is)&
                         &-frc31(ind,is)*dpoxyz_save(1,is,ind) &
                         &-frc32(ind,is)*dpoxyz_save(2,is,ind)
                end if
        !---------------------------------------------------
         200     continue   
        

          !do 300 is=1, nsys   
            !call rlubksb(is,amata,nnode,nnode, 1,nsys, 1,indx,bmata)
!300       continue


        do i = 1,nnode
            write(401,*) bmata(i,1:nsys)
        end do
        do i=1,nnoded
            !print *,cmat(i,1)
            write(403,*) cmat(i,1)
        enddo
!C                 
!C ** output the results, compute unkn[1:NNF] is dpdn,unkn[NNF+1:nnode]
! is potential
!C      
        do 500 ip=1, nsys 
        do 360 ind=1, nnode 
        bmat(ip)=(0.0d0, 0.0d0)
        do 350 is=1, nsys
350     bmat(ip)=bmat(ip)+bmata(ind,is)*rsn(ip,is)
        bmat(ip)=bmat(ip)/nsys
        unkn(ind,ip)=bmat(ip)        
360     continue
500     continue

!        WRITE(102,*) ' INODE,XP,YP,ZP,BMATA(INODE,1),UNKN(INODE,1' 
!        DO    INODE=1,  NNODE
!        XP=XYZ(1,INODE)
!        YP=XYZ(2,INODE)
!        ZP=XYZ(3,INODE) 
!        WRITE(102, 620)  INODE,XP,YP,ZP,BMATA(INODE,1),UNKN(INODE,1)
!        ENDDO
  
 620   FORMAT(1X,I4,5(1X,E13.6))                            

            
       RETURN
       END
