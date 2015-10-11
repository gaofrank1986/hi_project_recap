
C  
C **************************************************************** 
C *                                                              * 
C *    Test program of a novel BEM method basded on hybrid       * 
C *   integral equations for internal problem.                   * 
C *							         *
C *  							         *
C * 	 December 1, 2014                                        *
C *                                                              * 
C **************************************************************** 
C 
	PROGRAM Hyper_Integral_Equation

	 USE MVAR_MOD
	  USE PVAR_MOD
	  use INTVar_mod
	  USE WMAT_MOD
      
          !use hi_SIEPPEM
          use hi_intg

	  IMPLICIT  NONE  
	  INTEGER I,IFWKO,NDRFA,IS,IND,M,IP,IPOL,K,NDMAX
	  INTEGER NTnum,IFLAG_T

	  real(8) ::  AMFJ,R,EX(4),Alpha,WAVES  
	  real(8) ::  FAMPR(6),FAMPI(6),FORCER(6),FORCEI(6)
	  real(8) ::  PL_AMP(6),FORAMP
	  real(8) ::  FCD_AMR,  FCD_AMI
C !
C ! IFLAG_T,FCD_AMR, FCD_AMI: not used in this program
C !
C         DATA EX /  1.,  1., -1., -1./ 
C !
C ! FAMPR, FAMPI:   real and imaginary parts of body motion 
C !                 from frequency domain calculation
C ! FORCER, FORCEI: real and imaginary parts of wave force 
C !                 from frequency domain calculation
C ! PL_AMP:         amplitude of plotting variable
C
C ----------------------------------------
C Input data files

C         OPEN(8, FILE='INPUT/SOLIDANGLE.txt',    STATUS='OLD') 
C         open(67, FILE='INPUT/M.txt',    STATUS='OLD')
C         read(67,*) wl
C         stop


C ! -----------------------------------
C !  Output data files
C !
        OPEN(6,  FILE='OUTPUT/OUTScreen.txt',    STATUS='UNKNOWN')
        OPEN(9,  FILE='OUTPUT/OUTPUT1.txt',    STATUS='UNKNOWN')
        OPEN(10, FILE='OUTPUT/OUTPUT.txt' ,    STATUS='UNKNOWN')
        OPEN(11, FILE='OUTPUT/mesh_track.txt' ,    STATUS='UNKNOWN')
        OPEN(12, FILE='OUTPUT/output_amatrix.txt' ,    STATUS='UNKNOWN')


        OPEN(101,FILE='OUTPUT/OUTAMT.txt' ,    STATUS='UNKNOWN')	  
        OPEN(102,FILE='OUTPUT/OUTBMT.txt' ,    STATUS='UNKNOWN')
        OPEN(103,FILE='OUTPUT/OUTCMT.txt' ,    STATUS='UNKNOWN')
        OPEN(109,FILE='OUTPUT/OUTPUT9.txt' ,    STATUS='UNKNOWN')
        OPEN(108,FILE='OUTPUT/OUTPUT7.txt' ,    STATUS='UNKNOWN')
        OPEN(110,FILE='OUTPUT/OUTPUT10.txt' ,    STATUS='UNKNOWN')    

C            OPEN(108,FILE='OUTPUT/sing1_debug_info.txt',STATUS='UNKNOWN')


C
C ----------------------------------------------- 
C
C
	MFREE = 1
	NBETA=0
	 
        WRITE(10, *) ' Test on huper-singular integration'
        call read_wav_data()
        call output_wav_data()

        XC=0
        YC=0
        ZC=0

       call read_mesh()
!
	  write(11,*) ' isys=',isys,' nsys=',nsys
	  write(11,*) ' nelemb=',nelemb,' nelemf=',nelemf
	  write(11,*) ' nelem=',nelem,'  iopl=',ipol

       allocate(samb(nelem,16,0:8),sambxy(nelem,16,3))
       allocate(DSAMB(NELEM,16,6))


        call pre_mesh_2()!
!  --------------------------------------------
        ALLOCATE(ANGLE(NNODE),FrA3(NNODE),
     1	        FrC31(NNODE),FrC32(NNODE),FrC33(NNODE))

       ALLOCATE(AMATA(NNODE,NNODE,NSYS),
	1		  BMATA(NNODE,NSYS), INDX(NNODE,NSYS))
!       
       ALLOCATE(UNKN(NNODE,NSYS),  BKN(NNODED,NSYS),
	1		  UNKN_O(NNODE,NSYS),BKN_O(NNODED,NSYS),
	1	      DPDT(NNODE,NSYS))
	 ALLOCATE(DH(4,NNF,NSYS),DP(4,NNF,NSYS),Dposi(4,6))
!
        CALL BODYFD                  

        print *,"initialise hi data"
        call init_hi_var() 

         CALL TASSB0   
C 	  WRITE(10,*) '    AFTER TASSB0' 
! 
! =================
! 
C         do i = 1,NELEMF
C           write (200,999) i,ncon(i,1:8)
C           x = 0.5*(XYZ(1,NCON(2))+XYZ(1,NCON(6)))
C           Y = 0.5*(XYZ(2,NCON(2))+XYZ(2,NCON(6)))
C         end DO
        print *,"=================== main program ends ==============="
        end  program      

