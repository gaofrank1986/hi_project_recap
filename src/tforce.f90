!C 
!C ****************************************************** 
!C *                                                    * 
!C *  Evaluate the wave force on a 3-D body             * 
!C *                                                    * 
!C ****************************************************** 
!C 
SUBROUTINE TFORCE 

    use kinds
    use gaussian_info
    use motion
    use potential_mod,only:dpdt
    use time_mod
    use wave_func,only:dpot2
    use proj_cnst,only:ex,ey

    implicit   none 

    integer ip,ielem,nsamb,n,k 
    real(rk)  :: p(3),xp,yp,zp     
    real(rk)  ddum 

    forcew=0.0d0 

    do ielem=nelemf+1,  nelem  
        nsamb=16 
        if(ncn(ielem) .eq.6 ) nsamb=4 

        do  ip=1,  nsys    
            do  n =1, nsamb   

                ddum=0.0d0 
                do k=1,  ncn(ielem)  
                    ddum=ddum+ dpdt(ncon(ielem,k),ip)*samb(ielem,n,k)  
                end do

                p = [ex(ip),ey(ip),0.0d0]*sambxy(ielem,n,1:3)
                !xp=ex(ip)*sambxy(ielem,n,1) 
                !yp=ey(ip)*sambxy(ielem,n,2)   
                !zp=          sambxy(ielem,n,3)       

                ddum=ddum+dpot2(h,g,ampn,phi_w,beta,wkn,freq,timerk,rampf, &
                    p(1),p(2),p(3),nfreq,nwave,iorder)*samb(ielem,n,0)  

                forcew(1:3)= forcew(1:3) +ddum*(/ex(ip),ey(ip),1.0d0/)*dsamb(ielem,n,1:3)
                forcew(4:6) = forcew(4:6) +ddum*(/ey(ip),ex(ip),ex(ip)*ey(ip)/)*dsamb(ielem,n,4:6)
                !forcew(1)=forcew(1)+ddum* 
                !1                        exy(ip,1)*dsamb(ielem,n,1) 
                !forcew(2)=forcew(2)+ ddum* 
                !1                        exy(ip,2)*dsamb(ielem,n,2) 

                !forcew(3)=forcew(3) + ddum* dsamb(ielem,n,3)        
                !c 
                !forcew(4)=forcew(4)+ddum* 
                !2               exy(ip,2)* dsamb(ielem,n,4)  
                !forcew(5)=forcew(5)+ddum* 
                !2               exy(ip,1)* dsamb(ielem,n,5) 
                !forcew(6)=forcew(6)+ddum* 
                !2               exy(ip,2)* exy(ip,1)*dsamb(ielem,n,6) 


            end do
        end do
    end do
    forcew(:)=rho*forcew(:) 

end subroutine


