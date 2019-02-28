C     File name: FEM1DUNSYM                             Date: 22 Jan 2018                
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *          
C     *                                                                 *        
C     *                  Program  FEM1DUNSM by JN Reddy                 *               
C     *  Nonlinear finite element analysis program for one-dimensional  *
C     *  equations of the form (Chapter 3 of the nonlinear FEM book):   *
C     *  																  *
C     * -(d/dx)[a(x,u)(du/dx)]+b(x,u){du/dx)+c(x,u)u=f(x),  0 < x < AL  *
C     *                                                                 *
C     *  and based on its UNSYMMETRIC weak form.                        *
C     *                                                                 *               
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *              
C                                                                               
C     DESCRIPTION OF VARIABLES USED IN THE PROGRAM
C
C     MXELM    Maximum number of ELEMENTS the program is compiled to run
C     MXNEQ    Maximum number of EQUATIONS the program is compiled to run
C     MXFBW    Maximum allowable FULL BAND WIDTH of array GLK 
C     MXEBC    Maximum number of specified ESSENTIAL boundary conditions
C     MXBBC    Maximum number of speci. NATURAL or NEWTON boundary cond.
C     MXNOD    Maximum number of NODES the program is compiled to run 
C     *****    Values of all of the parameters defined above can be set by
C     *****    the user in PARAMETER statement of the program.
C
C     GLK      Global coefficient matrix, [K]; Note that the global
C              source vector is stored in the last column of [GLK]
C     GLX      Vector of global coordinates of the global nodes, {x}
C     GLU      Vector of global nodal values of the primary variable, {U}
C     GP1,GP2  Vectors used to define a weighted average of GLU for 
C              accelerated convergence [GP2 (r-2)nd iteration solution;
C              and GP1 (r-1)st iteration, also the initial guess vector]
C     NOD      Connectivity array (relates global nodes to element nodes)
C     DP       Vector of fractions of the source magnitude
C     ISPV     ISPV(I,1) node number of the Ith B.C. on PRIMARY variable
C              ISPV(I,2) Degree of freedom at the node of the Ith B.C.
C     VSPV     Specified value of the PRIMARY variable in the Ith B.C.
C     ISSV     ISPV(I,1) node number of the Ith B.C. on SECONDARY variable
C              ISPV(I,2) Degree of freedom at the node of the Ith B.C. 
C     VSSV     Specified value of the SECONDARY variable in the Ith B.C.
C     ISMB     ISMB(I,1) node number of the Ith MIXED boundary condition
C              ISMB(I,2) Degree of freedom at the node of Ith MIXED B.C.
C     UREF     Variable defined in the MIXED B.C. of the form 
C                             Q+(BETA0+BETAU*U+...)*(U-UREF)=0
C     BETA0,BETAU,... See the above line for the meaning of these variables
C
C     ELK      Element coefficient matrix, [Ke]
C     ELF      element source vector, {fe}
C     ELX      Vector of global coordinates of the element nodes
C     ELU      Vector of the primary variable at the element nodes
C
C     AX0,etc  Coefficients in the definition of a(x,u)=AX:
C              AX=AX0+AX1*x+AU1*u+AU2*u*u+AUX1*(du/dx)+AUX2*(du/dx)*(du/dx)
C     BX0,etc  Coefficients in the definition of c(x,u)=BX:
C              BX=BC0+BX1*x+BU1*u+BU2*u*u+BUX1*(du/dx)+BUX2*(du/dx)*(du/dx)
C     CX0,etc  Coefficients in the definition of c(x,u)=CX:
C              CX=AC0+CX1*x+CU1*u+CU2*u*u+CUX1*(du/dx)+CUX2*(du/dx)*(du/dx)
C     FX0,etc  Coefficients in the definition of f(x)=FX:
C              FX=FX0+FX1*x+FX2*x*x
C 
C     NONLIN   Flag for type of analysis:                  
C               = 0, LINEAR analysis                                              
C               = 1, NONLINEAR analysis with Direct iteration 
C               > 1, NONLINEAR analysis with Newton's iteration 
C 
C     NPRNT    Flag for print out of processor computations:                  
C               = 0, No printing of element matrices                                             
C               = 1, Print element matrices only
C               = 2, Print global matrices before BCs only 
C               > 2, Print global matrices after BCs only           
C   _________________________________________________________________________          
C                                                                               
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*20 FIN,FOUT
	PARAMETER(MXELM=50,MXNEQ=200,MXFBW=25,MXEBC=10,MXNBC=10,MXNOD=100)        
      DIMENSION GLK(MXNEQ,MXFBW),GLX(MXNOD),GLU(MXNEQ),DP(20),TITLE(20),
     *          GP1(MXNEQ),GP2(MXNEQ),NOD(MXELM,3),ISPV(MXEBC,2),
     *          ISSV(MXNBC,2),ISMB(MXNBC,2),VSPV(MXEBC),VSSV(MXNBC),
     *          UREF(MXNBC),BETA0(MXNBC),BETAU(MXNBC)
      COMMON/STF/ELK(3,3),ELF(3),ELX(3),ELU(3),W0(3),W1(3),W2(3),            
     1           AX0,AX1,AU1,AU2,AUX1,AUX2,CX0,CX1,CU1,CU2,CUX1,CUX2,
     2           BX0,BX1,BU1,BU2,BUX1,BUX2,FX0,FX1,FX2                           
      DATA IN,IT,NDF/5,6,1/                                        
C     
      WRITE(*,*)('INPUT FILE NAME = ')
      READ(*,1100)FIN
      WRITE(*,*)('OUTPUT FILE NAME = ')
      READ(*,1100)FOUT
      
      IN=5
      IT=6
      OPEN(IN, FILE=FIN, STATUS='OLD')
      OPEN(IT, FILE=FOUT, STATUS='UNKNOWN')
1100  FORMAT(20A)    
      CALL ECHO (IN,IT)                                                         
      READ(IN,470)TITLE                                                         
      READ(IN,*)AX0,AX1,AU1,AU2,AUX1,AUX2
	READ(IN,*)BX0,BX1,BU1,BU2,BUX1,BUX2
	READ(IN,*)CX0,CX1,CU1,CU2,CUX1,CUX2
      READ(IN,*)FX0,FX1,FX2
      READ(IN,*)X0,AL        ! origin of the coordinate system; domain length                              
      READ(IN,*)IEL,NEM     
      NDF=1                                                       
      NPE=IEL+1
      NNM=IEL*NEM+1                                                             
      NEQ=NNM*NDF                                                               
      READ(IN,*)NSPV,NSSV,NSMB ! no. of BCs for PVs, SVs, and mixed BCs    
	IF(NSPV.GT.0)THEN
         DO 10 NB=1,NSPV                                                      
         READ(IN,*) (ISPV(NB,J),J=1,2),VSPV(NB)                            
   10    CONTINUE                                                              
      ENDIF    
      IF (NSSV.GT.0)THEN                                                         
          DO 20 NB=1,NSSV                                                       
          READ(IN,*) (ISSV(NB,J),J=1,2),VSSV(NB)                            
   20     CONTINUE                                                              
      ENDIF
      IF(NSMB.NE.0)THEN                                                         
          DO 30 I=1, NSMB                                                       
   30     READ(IN,*) (ISMB(I,J),J=1,2),UREF(I),BETA0(I),BETAU(I)                          
      ENDIF
	READ(IN,*)NONLIN,NPRNT                                                                          
      IF(NONLIN.GT.0)THEN                                                        
         READ(IN,*)NLS,ITMAX,EPS,GAMA                                              
         READ(IN,*)(DP(I),I=1,NLS)                                                 
         READ(IN,*)(GP1(I),I=1,NNM)
      ENDIF
C                                                                               
C     ---------------------------------------------------------------           
C              E N D    O F   T H E   I N P U T    D A T A                      
C     ---------------------------------------------------------------           
C                                                                               
      DL=AL/NEM/IEL                                                        
      DO I=1,NNM                                                             
         GLX(I)=X0+DL*(I-1)                                                          
	ENDDO
      NHBW=(IEL+1)*NDF
	NBW=2*NHBW                                                          
      NN=NPE*NDF                                                                
      DO I=1,NPE
         NOD(1,I)=I                                                                
      ENDDO   
      DO 40 N=2,NEM                                                             
      DO 40 I=1,NPE                                                             
   40 NOD(N,I)=NOD(N-1,I)+IEL
C                                                   
      WRITE(IT,480)                                                             
      WRITE(IT,470)TITLE                                                        
      WRITE(IT,365)
      WRITE(IT,370)X0,AL,AX0,AX1,AU1,AU2,AUX1,AUX2,
     *                   BX0,BX1,BU1,BU2,BUX1,BUX2,
     *			       CX0,CX1,CU1,CU2,CUX1,CUX2,FX0,FX1,FX2
      WRITE(IT,410)NEM,IEL,NNM                                      
      WRITE(IT,360)(GLX(I),I=1,NNM)
C
	IF(NSPV.NE.0)THEN                                                         
         WRITE(IT,380)NSPV                                                         
         DO 50 IB=1,NSPV                                                      
         WRITE(IT,500)(ISPV(IB,J),J=1,2),VSPV(IB)                         
   50    CONTINUE                                                              
      ENDIF                                                                     
C                                                                               
      IF(NSSV.NE.0)THEN                                                         
         WRITE(IT,385)NSSV                                                         
         DO 60 IB=1,NSSV                                                      
   60    WRITE(IT,500)(ISSV(IB,J),J=1,2),VSSV(IB)                             
      ENDIF 
	IF(NSMB.NE.0)THEN                                                         
         WRITE(IT,390)NSMB                                                         
         DO 70 IB=1,NSMB                                                      
   70    WRITE(IT,500)(ISMB(IB,J),J=1,2),BETA0(IB),BETAU(IB),UREF(IB)                             
      ENDIF                                                                     
C                                               
      IF(NONLIN.GT.0)THEN
	   IF(NONLIN.EQ.1)THEN
	      WRITE(IT,396)
	   ELSE 
	      WRITE(IT,397)
	   ENDIF
         WRITE(IT,375)NLS,ITMAX,EPS,GAMA      
	   WRITE(IT,376)(GP1(I),I=1,NEQ)                    
      ELSE
         NLS=1
         ITMAX=1
	   WRITE(IT,398)
      ENDIF
      DO I=1,NEQ                                                             
         GP2(I)=0.0
      ENDDO                                                                
      NCOUNT=0
C     
      F0=0.0                                                             
      DO 250 NL=1,NLS                                                                     
      F0=F0+DP(NL)                                                                           
      ITER=0                                                                    
   80 ITER=ITER+1                                                               
      NCOUNT=NCOUNT+1                                                           
      IF(ITER.GT.ITMAX)THEN                                                
         WRITE(IT,330)                                                             
         STOP
      ENDIF                                                                     
C     
C     Initialize the global coefficient matrix (full bandwidth)                                    
C                                                                               
      DO 90 I=1,NEQ                                                            
          DO 90 J=1,NBW                                                           
   90        GLK(I,J)=0.0                                                            
C                                                                               
C     Element calculations and assembly (in full bandwidth) begins here                             
C                                                                               
      DO 140 N=1,NEM                                                            
         DO 100 I=1,NPE                                                             
         NI=NOD(N,I)                                                               
         ELU(I)=GAMA*GP2(NI)+(1.0-GAMA)*GP1(NI)                                   
  100    ELX(I)=GLX(NI)                                                              
         CALL ELMATRCS1D(IEL,NPE,NONLIN,F0)              
         IF(NPRNT.GT.0)THEN                                                        
            IF(N.EQ.1 .AND. NCOUNT.EQ.1)THEN                                        
               WRITE(IT,450)                                                             
               DO 110 I=1,NN                                                              
  110          WRITE(IT,350)(ELK(I,J),J=1,NN)                                           
               WRITE(IT,460)                                                             
               WRITE(IT,350)(ELF(I),I=1,NN)                                              
            ENDIF                                                                  
         ENDIF
         DO 130 I=1,NPE                                                            
            NR=(NOD(N,I)-1)*NDF                                                       
            DO 130 II=1,NDF                                                           
               NR=NR+1                                                                   
               L=(I-1)*NDF+II                                                            
               GLK(NR,NBW)=GLK(NR,NBW)+ELF(L)                                                    
               DO 120 J=1,NPE                                                            
                  NCL=(NOD(N,J)-1)*NDF                                                      
                  DO 120 JJ=1,NDF                                                           
                     M=(J-1)*NDF+JJ                                                            
                     NC=NCL-NR+JJ+NHBW                                                            
                     IF(NC.GT.0)THEN                                                       
                        GLK(NR,NC)=GLK(NR,NC)+ELK(L,M)
	               ENDIF
  120          CONTINUE                                                                  
  130    CONTINUE                                                                  
  140 CONTINUE                                                                  
C
      IF(NPRNT.EQ.2)THEN                                                    
         WRITE(IT,320)                                                             
         DO 150 I=1,NEQ                                                            
  150    WRITE(IT,350)(GLK(I,J),J=1,NBW-1)
	   WRITE(IT,325)
	   WRITE(IT,350)(GLK(I,NBW),I=1,NEQ)                                        
      ENDIF   
      IF(NONLIN.GT.0)THEN
	   DO 160 I=1,NEQ
  160	      GLU(I)=GAMA*GP2(I)+(1.0-GAMA)*GP1(I)                                   
	ENDIF    
C                                                                               
C     Modify the GLK to include specified BOUNDARY CONDITIONS             
C     
C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
	IF(NONLIN.GT.1)THEN
      IF(ITER.GT.0)THEN
	  DO NB=1,NSPV
	    VSPV(NB)=0.0
	  ENDDO
      ENDIF
	ENDIF
C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
      
	CALL BNDRYUNS1D(NONLIN,MXNEQ,MXFBW,MXEBC,MXNBC,NDF,NHBW,GLK,GLU,
	1     NSPV,NSSV,NSMB,ISPV,ISSV,ISMB,VSPV,VSSV,BETA0,BETAU,UREF)
C
      IF(NPRNT.GT.2)THEN                                                    
         WRITE(IT,320)                                                             
         DO 180 I=1,NEQ                                                            
  180    WRITE(IT,350)(GLK(I,J),J=1,NBW-1)
	   WRITE(IT,325)
	   WRITE(IT,350)(GLK(I,NBW),I=1,NEQ)                                        
      ENDIF       
C                    
C     Solve the assembled equations for the velocities using a banded           
C     equation solver. The solution is stored in 'GSTIF(I,NBW)'                 
C                                                                               
      CALL SOLVRUNSYM(GLK,MXNEQ,MXFBW,NEQ,NHBW)                                   
      IF(NONLIN.EQ.0 .OR. NCOUNT.EQ.1)THEN
         WRITE(IT,395)F0
         WRITE(IT,350)(GLK(I,NBW),I=1,NEQ)                                 
         IF(NONLIN.EQ.0)THEN
            STOP
         ENDIF
      ENDIF
C                                                                               
C     Save the previous iteration solution and update the current solution      
C                                                                               
C	 ******** write your statements ********
C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
      DO 200 I=1,NEQ
         GP2(I) = GP1(I)
         IF(NONLIN.LE.1) THEN
            GP1(I) = GLK(I,NBW)
         ELSE
            GP1(I) = GP1(I)+GLK(I,NBW)
         ENDIF
  200    CONTINUE      
C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
      
C                                                                               
C     Compute error and check for the convergence of the solution                                  
C                                                                               
      
C	 ******** write your statements ********
C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
      DNORM = 0.0
      DINORM = 0.0
      DO 220 IE = 1,NEQ
      DNORM = DNORM+GP1(IE)*GP1(IE)
  220 DINORM = DINORM+(GP1(IE)-GP2(IE))*(GP1(IE)-GP2(IE))
      TOLR = DSQRT(DINORM/DNORM)
      IF(TOLR.GT.EPS) THEN
         WRITE(IT,440) ITER,TOLR
         WRITE(IT,350)(GP1(I),I=1,NEQ)
         GOTO 80
      ELSE
         WRITE(IT,400) NL,F0
         WRITE(IT,420) ITER,TOLR
         WRITE(IT,350)(GP1(I),I=1,NEQ)
      ENDIF 
C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
      
  250 CONTINUE                                                                  
      PAUSE                                                                    
C                                                                               
C                    F   O   R   M   A   T    S                                 
C                                                                              
320   FORMAT(/,8X,'Global COEFFICIENT MATRIX',/)     
325   FORMAT(/,8X,'Global SOURCE VECTOR',/)                                                          
330   FORMAT(/,5X,'***** Convergence criterion is NOT satisfied *****')        
335   FORMAT(/,5X,'***** Number of time steps specified exceeded ****')        
350   FORMAT(5X,4E15.5)                                                        
360   FORMAT(5X,5E13.4)                                                        
365   FORMAT(/,5X,'A ONE-D SECOND-ORDER problem is analyzed')                
370   FORMAT(/,3X,'P R O B L E M   D A T A:',//,        
     1    5X,'Coordinate of node 1 ....................... =',E13.4,/,          
     2    5X,'Length of the domain modeled ............... =',E13.4,/,          
     3    5X,'Coefficient of 1 in AX, AX0 ................ =',E13.4,/,          
     4    5X,'Coefficient of x in AX, AX1 ................ =',E13.4,/,     
     5    5X,'Coefficient of u in AX, AU1 ................ =',E13.4,/,                 
     6    5X,'Coefficient of u*u in AX, AU2 .............. =',E13.4,/,
	7    5X,'Coefficient of du/dx in AX, AUX1 ........... =',E13.4,/,
	8    5X,'Coefficient of (du/dx)*(du/dx) in AX, AUX2.. =',E13.4,/,
	3    5X,'Coefficient of 1 in BX, BX0 ................ =',E13.4,/,          
     4    5X,'Coefficient of x in BX, BX1 ................ =',E13.4,/,     
     5    5X,'Coefficient of u in BX, BU1 ................ =',E13.4,/,                 
     6    5X,'Coefficient of u*u in BX, BU2 .............. =',E13.4,/,
	7    5X,'Coefficient of du/dx in BX, BUX1 ........... =',E13.4,/,
	8    5X,'Coefficient of (du/dx)*(du/dx) in BX, BUX2.. =',E13.4,/,
     3    5X,'Coefficient of 1 in CX, CX0 ................ =',E13.4,/,          
     4    5X,'Coefficient of x in CX, CX1 ................ =',E13.4,/,     
     5    5X,'Coefficient of u in CX, CU1 ................ =',E13.4,/,                 
     6    5X,'Coefficient of u*u in CX, CU2 .............. =',E13.4,/,
	7    5X,'Coefficient of du/dx in CX, CUX1 ........... =',E13.4,/,
	8    5X,'Coefficient of (du/dx)*(du/dx) in CX, CUX2.. =',E13.4,/,
     3    5X,'Coefficient of 1 in FX, FX0 ................ =',E13.4,/,     
     4    5X,'Coefficient of x in FX, FX1................. =',E13.4,/,                 
     5    5X,'Coefficient of x*x in FX, FX2 .............. =',E13.4,/)      
375   FORMAT(/,3X,'NONLINEAR  A N A L Y S I S   D A T A:',//,                
     1    5X,'Number of load steps used, NLS ............. =',I5,/,             
     2    5X,'Maximum allowable number of iterations ..... =',I5,/,             
     3    5X,'Allowable convergence tolerance (EPS) ...... =',E13.4,/,          
     4    5X,'Acceleration parameter:GAMA*GP2+(1-GAMA)*GP1 =',E13.4)
376   FORMAT(/,3X,'INITIAL (GUESS) VECTOR:',//,3X,8E13.5,8E13.5,8E13.5)            
380   FORMAT(/,3X,'B O U N D A R Y   D A T A:',//,               
     2    5X,'Number of specified PRIMARY variables ...... =',I5,/,             
     3    5X,'Specified PRIMARY degrees of freedom ....... :',/)                
385   FORMAT(/,5X,'Number of specified SECONDARY variables =',I5,         
     2  /,5X,'Specified SECONDARY variables .............. :',/)    
390   FORMAT(/,5X,'Number of specified MIXED B.C. ........ =',I5,         
     2  /,5X,'Specified MIXED variables (BETA0+BETAU)*(U-UREF).. :',/)                           
395   FORMAT(/,5X,'First iteration solution at source value ... =',
     1   E12.4,/)  
396   FORMAT(/,5X,'*** NONLINEAR Analysis with DIRECT ITERATION ***',/)  
397   FORMAT(/,5X,'*** NONLINEAR Analysis with NEWTON ITERATION ***',/)        
398   FORMAT(/,5X,'*******  L I N E A R    A N A L Y S I S  *******',/)        
400   FORMAT(//,3X,'STEP NUMBER =',I3,3X,'SOURCE VALUE =',E12.4)
410   FORMAT(/,3X,'M E S H   D A T A:',//,                       
     1    5X,'Number of elements in the mesh, NEM......... =',I5,/,             
     2    5X,'Type of element (1: linear; 2: quadratic) .. =',I5,/,             
     3    5X,'Number of nodes in the mesh, NNM ........... =',I5,/,             
     4    5X,'GLOBAL COORDINATES OF THE GLOBAL NODES ..... =',/)                
420   FORMAT(/,5X,'Number of iterations taken for convergence.. =',I3,/,          
     2    5X,'The Euclidean norm of the solution error ... =',E12.4,/,         
     3    5X,'The CONVERGED solution is:',/)
430   FORMAT (2X,63('_'))                                                       
440   FORMAT(/,5X,'Iteration Number =',I3,5X,'Error=',E12.4,/)
450   FORMAT (/,8X,'Element COEFFICENT matrix:',/)                               
460   FORMAT (//,8X,'Element SOURCE vector:',/)                                  
470   FORMAT (20A4)                                                             
480   FORMAT (80('_'),/,8X,'Output from computer program  FEM1DIUNSYM            
     *  by J. N. Reddy',/,80('_'),//)                                      
490   FORMAT (10X,14I5)                                                         
500   FORMAT (10X,2I5,3E13.4)                                         
      END
	                                                                       

      SUBROUTINE ELMATRCS1D(IEL,NPE,NONLIN,F0)
C     __________________________________________________________________
C
C     SFL(I)      Shape (interpolation) functions of the Lagrange type
C     GDSFL(I)    Global derivative (w.r.t. x) of SFL(I)
C     GJ          Jacobian of the transformation
C     GAUSPT(I,J) Ith Gauss point in the J-point Gauss rule
C     GAUSWT(I,J) Ith Gauss weight in the J-point Gauss rule
C     NGP         Number of Gauss points used to evaluate integrals
C     EL          Length of the element
C     XI          Normalized local coordinate such that -1 < xi < 1
C     U           Primary variable (function of x)
C     DU          Global derivative of U (dU/dx)
C     __________________________________________________________________
C    	        
      IMPLICIT REAL*8(A-H,O-Z)                                                  
      DIMENSION GAUSPT(5,5),GAUSWT(5,5),TANG(3,3)                        
      COMMON/SHP/SFL(4),GDSFL(4)                       
      COMMON/STF/ELK(3,3),ELF(3),ELX(3),ELU(3),W0(3),W1(3),W2(3),            
     1           AX0,AX1,AU1,AU2,AUX1,AUX2,CX0,CX1,CU1,CU2,CUX1,CUX2,
     2           BX0,BX1,BU1,BU2,BUX1,BUX2,FX0,FX1,FX2                           
C                                                                               
      DATA GAUSPT/5*0.0D0,-.57735027D0,.57735027D0,3*0.0D0,-.77459667D0,
     1 0.0D0,.77459667D0,2*0.0D0,-.86113631D0,-0.33998104D0,.33998104D0,
     2 0.86113631D0,.0D0,-.906180D0,-.538469D0,.0D0,.538469D0,.906180D0/
      DATA GAUSWT/2.0D0,4*0.0D0,2*1.0D0,3*0.0D0,.55555555D0,.88888888D0, 
     1 0.55555555D0,2*0.0D0,.34785485D0,2*.65214515D0,.34785485D0,0.0D0,
     2 0.236927D0,0.478629D0,0.568889D0,0.478629D0,0.236927D0/         
C    
C     Initialize coefficient matrices and column vectors
C     
C	 ******** write your statements ********
C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
      NGP = IEL+1
      EL = ELX(NPE)-ELX(1)
      DO 10 I=1,NPE
          ELF(I) = 0.0
          DO 10 J=1,NPE
              IF(NONLIN.GT.1)THEN
                  TANG(I,J) = 0.0
              ENDIF
10            ELK(I,J) = 0.0
      
      DO 50 NI=1,NGP
          XI = GAUSPT(NI,NGP)
      CALL INTERPLN1D(ELX,GJ,IEL,NPE,XI)
      CNST = GJ*GAUSWT(NI,NGP)
      

C      X = ELX(1)+0.5*(1.0+XI)*EL
C      AX = AX0+AX1*X
C      BX = BX0+BX1*X
C      CX = CX0+CX1*X
C      FX = FX0+FX1*X+FX2*X*X
      
      IF(NONLIN.GT.0) THEN
          U = 0.0
          DU = 0.0
          X = 0.0
      DO 20 I=1,NPE
          U = U+SFL(I)*ELU(I)
          X = X+SFL(I)*ELX(I)
          DU = DU+GDSFL(I)*ELU(I)
      
20    CONTINUE
      
           
C    
C     Compute the ELK(I,J), and ELKT(I,J)=ELK(I,J)+TANG(I,J) 
C                using Gauss quadrature
C     
C	 ******** write your statements ********
      
      AX = AX0+AX1*X+AU1*U+AUX1*DU+AU2*U*U+AUX2*DU*DU
      BX = BX0+BX1*X+BU1*U+BUX1*DU+BU2*U*U+BUX2*DU*DU
      CX = CX0+CX1*X+CU1*U+CUX1*DU+CU2*U*U+CUX2*DU*DU
      FX = FX0+FX1*X+FX2*X*X
      
      IF(NONLIN.GT.1) THEN
              
              HA1 = (AU1+2.0*AU2*U)*DU
              HA2 = (AUX1+2.0*AUX2*DU)*DU
              HB1 = (BU1+2.0*BU2*U)*DU
              HB2 = (BUX1+2.0*BUX2*DU)*DU
              HC1 = (CU1+2.0*CU2*U)*U
              HC2 = (CUX1+2.0*CUX2*DU)*U
      ENDIF
      ENDIF
      
      DO 40 I=1,NPE
          ELF(I) = ELF(I)+FX*SFL(I)*CNST
          DO 40 J=1,NPE
              S11 = GDSFL(I)*GDSFL(J)*CNST
              S22 = SFL(I)*GDSFL(J)*CNST
              S33 = SFL(I)*SFL(J)*CNST
              S44 = GDSFL(I)*SFL(J)*CNST
              
              ELK(I,J) = ELK(I,J)+AX*S11+BX*S22+CX*S33
          IF(NONLIN.GT.1) THEN
              HA = HA1*S44+HA2*S11
              HB = HB1*S33+HB2*S22
              HC = HC1*S33+HC2*S22
              TANG(I,J)=TANG(I,J)+HA+HB+HC
          ENDIF          
40    CONTINUE
50    CONTINUE
      
      
C
C     Compute the residual vector and tangent coefficient matrix
C     
C	 ******** write your statements ********
      
      IF(NONLIN.GT.1) THEN
      DO 100 I=1,NPE
          DO 100 J=1,NPE
100           ELF(I) = ELF(I)-ELK(I,J)*ELU(J)
      DO 120 I=1,NPE
          DO 120 J=1,NPE
120           ELK(I,J) = ELK(I,J)+TANG(I,J)
      ENDIF   
      
C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
      
      RETURN                                                                    
      END
	                                                                       
      
      SUBROUTINE INTERPLN1D(ELX,GJ,IEL,NPE,XI)
C     __________________________________________________________________
C
C            DSFL(I)    Derivative of SFL(I) with respect to xi
C     __________________________________________________________________
C
      IMPLICIT REAL*8(A-H,O-Z)                                                  
      COMMON /SHP/SFL(4),GDSFL(4)                       
      DIMENSION DSFL(4),ELX(3)                                   
C                                                                               
C     Interpolations functions and their derivatives                                   
C                                                                               
C     LINEAR shape (interpolation) functions
C                                                                               
      IF(IEL.EQ.1)THEN
         SFL(1)=0.5*(1.0-XI)                                                       
         SFL(2)=0.5*(1+XI)                                                         
         DSFL(1)=-0.5                                                              
         DSFL(2)=0.5                                                               
      ENDIF                                                                     
C                                                                               
C     QUADRATIC shape (interpolations) functions                                                
C                                                                               
      IF(IEL.EQ.2)THEN
         SFL(1)=0.5*(XI-1.0)*XI                                                    
         SFL(2)=1.0-XI*XI                                                          
         SFL(3)=0.5*(XI+1.0)*XI                                                    
         DSFL(1)=XI-0.5                                                            
         DSFL(2)=-2.0*XI                                                           
         DSFL(3)=XI+0.5                                                            
      ENDIF
      GJ=0.0                                                                    
      DO 10 I=1,NPE                                                             
10    GJ=GJ+DSFL(I)*ELX(I)                                                      
      DO 30 I=1,NPE                                                             
      GDSFL(I)=DSFL(I)/GJ                                                       
30    CONTINUE                                                                  
      RETURN                                                                    
      END     
                                                             
      
      SUBROUTINE SOLVRUNSYM(A,NRMAX,NCMAX,N,ITERM)                                 
C     _________________________________________________________________         
C                                                                               
C       Solver for BANDED UNSYMMETRIC system of algebraic equations         
C     _________________________________________________________________         
C                                                                               
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
      DIMENSION A(NRMAX,NCMAX)                                                  
      CERO=1.0D-15                                                              
      PARE=CERO**2                                                              
      NBND=2*ITERM                                                              
      NBM=NBND-1                                                                
C                                                                               
C     Begin elimination of the lower left
C                                                                               
      DO 80 I=1,N                                                               
      IF (DABS(A(I,ITERM)).LT.CERO) GO TO 10                                    
      GO TO 20                                                                  
   10 IF (DABS(A(I,ITERM)).LT.PARE) GO TO 110                                   
   20 JLAST=MIN0(I+ITERM-1,N)                                                   
      L=ITERM+1                                                                 
      DO 40 J=I,JLAST                                                           
      L=L-1                                                                     
      IF (DABS(A(J,L)).LT.PARE) GO TO 40                                        
      B=A(J,L)                                                                  
      DO 30 K=L,NBND                                                            
   30 A(J,K)=A(J,K)/B                                                           
      IF (I.EQ.N) GO TO 90                                                      
   40 CONTINUE                                                                  
      L=0                                                                       
      JFIRST=I+1                                                                
      IF (JLAST.LE.I) GO TO 80                                                  
      DO 70 J=JFIRST,JLAST                                                      
      L=L+1                                                                     
      IF (DABS(A(J,ITERM-L)).LT.PARE) GO TO 70                                  
      DO 50 K=ITERM,NBM                                                         
   50 A(J,K-L)=A(J-L,K)-A(J,K-L)                                                
      A(J,NBND)=A(J-L,NBND)-A(J,NBND)                                           
      IF (I.GE.N-ITERM+1) GO TO 70                                              
      DO 60 K=1,L                                                               
   60 A(J,NBND-K)=-A(J,NBND-K)                                                  
   70 CONTINUE                                                                  
   80 CONTINUE                                                                  
   90 L=ITERM-1                                                                 
      DO 100 I=2,N                                                              
      DO 100 J=1,L                                                              
      IF (N+1-I+J.GT.N) GO TO 100                                               
      A(N+1-I,NBND)=A(N+1-I,NBND)-A(N+1-I+J,NBND)*A(N+1-I,ITERM+J)              
  100 CONTINUE                                                                  
      RETURN                                                                    
  110 WRITE (6,140) I,A(I,ITERM)                                                
      STOP                                                                      
C                                                                               
  140 FORMAT (/,2X,'Computation stopped in SLVUNSYM because zero appears        
     * on the main diagonal *** Eqn no. and value:',I5,E12.4)                   
      END
	
           
      SUBROUTINE ECHO (IN,IT)                                                   
C                                                                               
      DIMENSION AA(20)                                                          
      WRITE(IT,40)                                                              
C
      DO                                                                       
         READ(IN,30,END=20) AA                                                     
         WRITE(IT,30) AA                                                           
      ENDDO                                                                  
C             
   20 CONTINUE                                                      
      REWIND(IN)                                                                
      WRITE(IT,50)                                                              
      RETURN                                                                    
   30 FORMAT(20A4)                                                              
   40 FORMAT(5X,'*** ECHO OF THE INPUT DATA STARTS ***',/)                      
   50 FORMAT(/,5X,'**** ECHO OF THE INPUT DATA ENDS ****',/)                      
      END                                                                       


	SUBROUTINE BNDRYUNS1D(NONLIN,MXNEQ,MXFBW,MXEBC,MXNBC,NDF,NHBW,GLK,
	1    GLU,NSPV,NSSV,NSMB,ISPV,ISSV,ISMB,VSPV,VSSV,BETA0,BETAU,UREF)    
C     __________________________________________________________________        
C                                                                               
C     The subroutine is used to implement specified boundary conditions         
C           on BANDED UNSYMMETRIC system of finite element equations               
C     __________________________________________________________________        
C                                                                               
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
      DIMENSION ISPV(MXEBC,2),ISSV(MXNBC,2),ISMB(MXNBC,2),VSPV(MXEBC),
     1          VSSV(MXEBC),UREF(MXNBC),BETA0(MXNBC),BETAU(MXNBC),
     2          GLU(MXNEQ),GLK(MXNEQ,MXFBW)
C
      NBW=2*NHBW
C                                                                               
C     Include specified PRIMARY degrees of freedom                              
C                                                                               
 	IF(NSPV.NE.0)THEN                                                                          
         DO 120 NP=1,NSPV                                                            
            NB=(ISPV(NP,1)-1)*NDF+ISPV(NP,2)                                  
            DO 110 J=1,NBW                                                            
  110          GLK(NB,J)=0.0D0                                                         
            GLK(NB,NHBW)=1.0D0                                                      
  120       GLK(NB,NBW)=VSPV(NP)
      ENDIF
C                                                                               
C     Modify the source vector to include specified nonzero SECONDARY VARIABLES             
C                                                                               
      IF(NSSV.NE.0)THEN                                                        
         DO 130 NS=1,NSSV
	   NB=(ISSV(NS,1)-1)*NDF+ISSV(NS,2)                                  
  130    GLK(NB,NBW)=GLK(NB,NBW)+VSSV(NS)                                         
      ENDIF                                                                     
C                                                                               
C     Implement the specified MIXED BOUNDARY CONDITIONS                       
C         
	IF(NSMB.NE.0)THEN                                                                      
         DO 150 MB=1,NSMB                                                            
         NB=(ISMB(MB,1)-1)*NDF+ISMB(MB,2)
	   IF(NONLIN.LE.1)THEN 
         GLK(NB,NHBW)=GLK(NB,NHBW)+BETA0(MB)+BETAU(MB)*GLU(NB)
         GLK(NB,NBW)=GLK(NB,NBW)+UREF(MB)*(BETA0(MB)+BETAU(MB)*GLU(NB))                             
         ELSE
         GLK(NB,NHBW)=GLK(NB,NHBW)+BETA0(MB)+2.0*BETAU(MB)*GLU(NB)
     *	                                  -UREF(MB)*BETAU(MB)
         GLK(NB,NBW)=GLK(NB,NBW)+UREF(MB)*(BETA0(MB)+BETAU(MB)*GLU(NB))
     *                          -(BETA0(MB)+BETAU(MB)*GLU(NB))*GLU(NB)
         ENDIF                               
  150	   CONTINUE
      ENDIF
	RETURN
	END