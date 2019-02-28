C     Date: April 8, 2010 (revised program FEM2DUNS.FOR)
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C     *                     Program  FEM2DUNS                       *
C     *   (AN IN-CORE FINITE ELEMENT ANALYSIS COMPUTER PROGRAM)     *
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C ----------------------------------------------------------------------
C    . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
C    .                                                               .
C    .   This is a nonlinear finite element computer program for the .
C    . analysis of heat transfer  and  other field problems governed .
C    . by general Poisson's equation in  two dimensions.  Convection .
C    . type boundary conditions are included.                        .
C    .                                                               .
C    .     The program  uses linear and quadratic RECTANGULAR (only) .
C    . elements  with  isoparametric  formulations (quadrilaterals). .
C    . Only  one  type  of element are allowed for a problem  (i.e., .
C    . two different types of elements cannot be used in a problem). .
C    . Only linear problems can be solved with this program.         .
C    .                                                               .
C    .     The main objective of this  program  is to illustrate how .
C    . finite element formulations presented in Chapters  8, 9 and 13.
C    . can be implemented on a computer  and used in the analysis of .
C    . engineering problems.
C    .                                                               .
C    . This is a modified version of the computer program  FEM2D  in .
C    . AN INTRODUCTION TO THE FINITE ELEMENT METHOD  by J. N. Reddy, .
C    . McGraw-Hill, New York, 2006 (Third Edition).                  .
C    . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
C    ___________________________________________________________________
C
C           DESCRIPTION OF SOME KEY VARIABLES USED IN THE PROGRAM
C
C      {ELF}     Vector of element nodal source vector
C      [ELK]     Element coefficient matrix
C      {ELU}     Vector of element nodal values of primary variables
C      {ELXY}    Vector of elemental global coordinates:
C                ELXY(I,1)=x-coordinate;  ELXY(I,2)=y-coordinate
C      [GLK]     Global coefficient matrix (in general, unsymmetric)
C      {GPU}     Global solution vector from the Previous iteration
C      {GLU}     Global solution vector from the Current iteration
C
C      NDF       Number of degrees of freedom per node (NDF=1 here)
C      NEQ       Total number of equations in the problem (=NNM*NDF)
C      NHBW      Half band width of the global coefficient matrix, GLK
C      NBW       =2 times NHBW (full bandwidth + 1)
C      NN        Total number of degrees of freedom per element (NN=NPE)
C      ---------------------------------------------------------------
C
C          DESCRIPTION OF PARAMETERS USED TO DIMENSION THE ARRAYS
C
C      MAXCNV... Maximum number of elements with convection B.C.
C      MAXELM... Maximum number of elements allowed in the program
C      MAXNOD... Maximum number of nodes allowed in the program
C      MAXNX.... Maximum number of allowed subdivisions DX(I) along x
C      MAXNY.... Maximum number of allowed subdivisions DY(I) along y
C      MAXSPV... Maximum number of specified primary variables
C      MAXSSV... Maximum number of specified secondary variables
C      NCMAX.... Actual column dimension of: [GLK],[GLM],{GLU}, and {GLF}
C      NRMAX.... Actual row dimension of: [GLK] and [GLM]
C
C      NOTES:___ The actual row dimension of the assembled coefficient
C                matrix  should be greater than or equal to  the total
C                number of algebraic equations in the FE model.
C
C                The actual column-dimension of  assembled coefficient
C                matrix  should be greater than or equal to  the  half
C                bandwidth for static analysis.
C
C                The values of MAXNEQ, MAXFBW, MAXELM, MAXNOD, MAXCNV,
C                MAXSSV and MAXSPV in the 'PARAMETER' statement should
C                modified as required by the size of the problem.
C ______________________________________________________________________
C
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (MAXNEQ=1000, MAXFBW=1000, MAXELM=1000, MAXNOD=2000,
     1           MAXSPV=1000, MAXSSV=500, MAXCNV=20, MAXNX=50, MAXNY=50)
C
      DIMENSION ISPV(MAXSPV,2),VSPV(MAXSPV),ISSV(MAXSSV,2),VSSV(MAXSSV),
     1          VSSV0(MAXSSV)
	DIMENSION IBN(MAXCNV),INOD(MAXCNV,3),BETA(MAXCNV),UREF(MAXCNV)
      DIMENSION GLF(MAXNEQ),TITLE(20),IBS(3),IBP(3),GLK(MAXNEQ,MAXFBW)
      DIMENSION NOD(MAXELM,9),GLXY(MAXNOD,2),DX(MAXNX),DY(MAXNY),
     1          GLXYP(MAXNOD,2)
      DIMENSION GLU(MAXNEQ),GPU(MAXNEQ),DP(100),GCU(MAXNEQ)
C      DIMENSION GLV(MAXNEQ),GLA(MAXNEQ),GLM(MAXNEQ,MAXFBW),GLP(MAXNEQ)
C
      COMMON/STF/ELF(18),ELK(18,18),ELXY(9,2),ELU(18),ELXYP(9,2)
C      COMMON/SOURCE/FX0,FY0,FXX,FXY,FYX,FYY
      COMMON/CNTC/C(3,3)

      
      COMMON/IO/IN,IT,IP,IE,IS
C
C     ****************************************************************
C     *                                                              *
C     *               P R E P R O C E S S O R   U N I T              *
C     *                                                              *
C     ****************************************************************
C
      IN=5
      IT=6
      IP = 7
      IE = 8
      IS = 9
      open(in, file='EX945.inp')
      open(IT, file='EX945.out')
      OPEN(IP, FILE='U.dat')
      OPEN(IE, FILE='STRESS.dat')
      OPEN(IS, FILE='PKSTRESS.dat')
      CALL ECHO(IN,IT)
      ICONV=0
C      NSSV=0
C
C     R E A D I N   T H E   I N P U T    D A T A    H E R E
C
      READ(IN,400) TITLE
      READ(IN,*) NPE,NGPF,NGPR,MESH,NPRNT,IGRAD,NONLIN,NEWTON
      
      IF(NPE.LE.4) THEN
         IEL=1
      ELSE
         IEL=2
      ENDIF
C
C     If mesh CANNOT be generated by the program, read the mesh data in
C     the next three statements
C
      IF(MESH.EQ.0) THEN
         READ(IN,*) NEM,NNM
         DO 10 N=1,NEM
   10    READ(IN,*) (NOD(N,I),I=1,NPE)
         READ(IN,*) ((GLXY(I,J),J=1,2),I=1,NNM)
      ELSE
C
C     When mesh is to be generated for rectangular domains, call program
C     MSH2DRectangular, which requires the following data:
C
         READ(IN,*) NX,NY
         READ(IN,*) X0,(DX(I),I=1,NX)
         READ(IN,*) Y0,(DY(I),I=1,NY)
         

         CALL MSH2DR (IEL,NX,NY,NPE,NNM,NEM,NOD,DX,DY,X0,Y0,
     *                GLXY,MAXELM,MAXNOD,MAXNX,MAXNY)
      ENDIF
      
            
      DO I=1,NNM
          GLXYP(I,1) = GLXY(I,1)
          GLXYP(I,2) = GLXY(I,2)
      ENDDO
      
      IF (NPE.EQ.8)THEN
      INDEX = 17
      ELSEIF(NPE.EQ.9)THEN
      INDEX = 22
      ENDIF
      
           
C
      NDF = 2
      NEQ=NNM*NDF
      NN=NPE*NDF
C
C     Compute the half bandwidth of the global coefficient matrix
C
      NHBW=0
      DO 20 N=1,NEM
      DO 20 I=1,NPE
      DO 20 J=1,NPE
      NW=(IABS(NOD(N,I)-NOD(N,J))+1)*NDF
   20 IF (NHBW.LT.NW) NHBW=NW
      NBW=2*NHBW
C
C     Read specified primary and secondary degrees of freedom: node
C     number, local degree of freedom number, and specified value.
C
      READ(IN,*) NSPV
      IF(NSPV.NE.0) THEN
         READ(IN,*) ((ISPV(I,J),J=1,2),I=1,NSPV)
         READ(IN,*) (VSPV(I),I=1,NSPV)
      ENDIF
      READ(IN,*) NSSV
      IF(NSSV.NE.0) THEN
         READ(IN,*) ((ISSV(I,J),J=1,2),I=1,NSSV)
         READ(IN,*) (VSSV(I),I=1,NSSV)
      ENDIF
C
      WRITE(IT,400) TITLE
      WRITE(IT,910)
      WRITE(IT,890)
      WRITE(IT,910)
C
C     Read the coefficients of the differential equation modeled
C
C     A11 = A10 + A1X*X + A1Y*Y; A22 = A20 + A2X*X + A2Y*Y; A00=CONST.
C
C      READ(IN,*)A10,A1X,A1Y
C      READ(IN,*)A20,A2X,A2Y
C      READ(IN,*)A00,A0X,A0Y
C      WRITE(IT,420) A10,A1X,A1Y,A20,A2X,A2Y,A00,A0X,A0Y
C      READ(IN,*)ICONV
      IF(ICONV.NE.0) THEN
         READ(IN,*) NBE
         READ(IN,*) (IBN(I),BETA(I),UREF(I),I=1,NBE)
         READ(IN,*) ((INOD(I,J),J=1,2),I=1,NBE)
         WRITE(IT,440) NBE
         DO 30 I=1,NBE
   30    WRITE(IT,860) IBN(I),(INOD(I,J),J=1,2),BETA(I),UREF(I)
      ENDIF
C      READ(IN,*) F0,FX,FY
C      WRITE(IT,430) F0,FX,FY
C
C     Input data for nonlinear analysis
C
C      IF(NONLIN.GT.0)THEN
C         READ(IN,*)A1U,A1UX,A1UY
C         READ(IN,*)A2U,A2UX,A2UY
C         WRITE(IT,980)A1U,A1UX,A1UY,A2U,A2UX,A2UY
C      ENDIF
C      
      READ(IN,*)LFORM
      READ(IN,*)AL,BL,THIKNS

      READ(IN,*)E,ANU
      E1 = E
	E2 = E
      ANU12 = ANU
	G12 = E/2.0/(1.0 + ANU)
      
      IF(NONLIN.GT.0)THEN
         READ(IN,*)NLS,ITMAX,EPS,REX
         READ(IN,*)(DP(I),I=1,NLS)
         WRITE(IT,990)NLS,ITMAX,NONLIN,EPS,REX
      ENDIF
C
C     *****     E N D   O F   T H E    D A T A   I N P U T       *****
C
      WRITE(IT,800)
      WRITE(IT,560) IEL,NPE,NDF,NEM,NNM,NEQ,NHBW
      IF(MESH.EQ.1) WRITE(IT,570) NX,NY
      WRITE(IT,700) NSPV
      IF(NSSV.NE.0) THEN
         WRITE(IT,710) NSSV
         WRITE(IT,720)
         DO 80 IB=1,NSSV
   80    WRITE(IT,960) (ISSV(IB,JB),JB=1,2),VSSV(IB)
      ENDIF
C
      IF(NPRNT.EQ.1) THEN
         WRITE(IT,690)
         DO 100 I=1,NEM
  100    WRITE(IT,900) I,(NOD(I,J),J=1,NPE)
      ENDIF
C
      WRITE(IT,910)
      WRITE(IT,580)
      WRITE(IT,910)
      DO 150 IM=1,NNM
      DO 110 K=1,NDF
      IBP(K)=0
  110 IBS(K)=0
      IF(NSPV.NE.0) THEN
         DO 120 JP=1,NSPV
         NODE=ISPV(JP,1)
         NDOF=ISPV(JP,2)
         IF(NODE.EQ.IM) THEN
            IBP(NDOF)=NDOF
         ENDIF
  120    CONTINUE
      ENDIF
C
      IF(NSSV.NE.0) THEN
         DO 140 JS=1,NSSV
         NODE=ISSV(JS,1)
         NDOF=ISSV(JS,2)
         IF(NODE.EQ.IM) THEN
            IBS(NDOF)=NDOF
         ENDIF
  140    CONTINUE
      ENDIF
C
      IF(NDF.EQ.1) THEN
         WRITE(IT,870)IM,(GLXY(IM,J),J=1,2),(IBP(K),K=1,NDF),
     *                   (IBS(K),K=1,NDF)
      ELSE
         WRITE(IT,920)IM,(GLXY(IM,J),J=1,2),(IBP(K),K=1,NDF),
     *                   (IBS(K),K=1,NDF)
      ENDIF
  150 CONTINUE
      WRITE(IT,910)
      WRITE(IT,480) NGPF,NGPR
	IF(NONLIN.GT.0)THEN
	   IF(NONLIN.EQ.1)THEN
	      WRITE(IT,830)
	   ELSE
	      WRITE(IT,820)
	   ENDIF
	ENDIF
      
      H = THIKNS
      ANU21 = ANU12*E2/E1
      C11 = E1/(1.0-ANU12*ANU21)
      C22 = C11*E2/E1
      C12 = ANU12*C22
      C66 = G12
      
      DO 151 I=1,3
      DO 151 J=1,3
151   C(I,J) = 0.0
      
      C(1,1) = C11
      C(1,2) = C12
      C(2,1) = C12
      C(2,2) = C22
      C(3,3) = C66


C
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C     *                                                             *
C     *                   P R O C E S S O R   U N I T               *
C     *                                                             *
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      DO 158 I=1,NEQ
         GCU(I) = 0.0
  158	   GPU(I) = 0.0      
      
      FS = 0.0
      NCOUNT = 0
      
      DO 175 NL=1,NLS
231         FS = FS + DP(NL)
      PRINT*,DABS(FS)

      ITER=0
  160 ITER=ITER+1
      NCOUNT = NCOUNT + 1
      IF(ITER.GT.ITMAX)THEN
         WRITE(IT,1200)ITMAX
         WRITE(IT,1050)
         WRITE(IT,940)
         WRITE(IT,670)
         WRITE(IT,940)
         IF(NDF.EQ.1)THEN
         DO I=1,NNM
            WRITE(IT,950)I,(GLXY(I,J),J=1,2),GCU(I)
         ENDDO
         ELSE
         DO I=1,NNM
            I1 = (I-1)*NDF +1
            I2 = I*NDF
            WRITE(IT,950)I,(GLXY(I,J),J=1,2),(GCU(K),K=I1,I2)
         ENDDO
         ENDIF
         WRITE(IT,940)
         STOP
      ENDIF
C
C     Initialize the global coefficient matrices and vectors
C
      DO 180 I=1,NEQ
          DO 180 J=1,NBW                                                           
  180        GLK(I,J)=0.0                                                            
C
C     Do-loop on the number of ELEMENTS to compute element matrices
C     and their assembly begins here
C
      DO 250 N=1,NEM
      LL = 0
      DO 200 I=1,NPE
      NI=NOD(N,I)
      IF(NDF.EQ.1)THEN
         ELU(I) = (1.0 - REX)*GCU(NI) + REX*GPU(NI)  ! TRANSFER OF CURRENT SOLUTION
      ELSE
         ND=(NI-1)*NDF
C         L=NI*NDF-1
C         K=I*NDF-1
C         ELU(K) = (1.0 - REX)*GCU(L) + REX*GPU(L)
C         ELU(K+1) = (1.0 - REX)*GCU(L+1) + REX*GPU(L+1)
	   DO 95 J=1,NDF
		   ND=ND+1
	       LL=LL+1
   95	       ELU(LL)=(1.0 - REX)*GCU(ND) + REX*GPU(ND)
      ENDIF
      ELXY(I,1)=GLXY(NI,1)
      ELXY(I,2)=GLXY(NI,2)
  200 CONTINUE
      
      
C
C     Call subroutine ELKMFR (for Rectangular elements) to compute
C     the eLement "K", "M" and "F" matrices.
C
      F0 = 0.0
      CALL ELMATRCS2D(ICONV,IEL,NPE,NN,NGPF,NGPR,NONLIN,
     *    NOD,GLXY,NBE,IBN,INOD,BETA,UREF,MAXELM,MAXNOD,MAXCNV,
     *    MAXSPV,MAXSSV,LFORM,THIKNS,F0)
      
      
      
      IF(NPRNT.EQ.1 .OR. NPRNT.EQ.3) THEN
         IF(N.EQ.1) THEN
C
C     Print element matrices and vectors (only when NPRNT=1 or NPRNT=3)
C
            WRITE(IT,610)
            DO 210 I=1,NN
  210       WRITE(IT,930)(ELK(I,J),J=1,NN)
            WRITE(IT,630)
            WRITE(IT,930)(ELF(I),I=1,NN)
         ENDIF
      ENDIF
C
C     ASSEMBLE element matrices to obtain global matrices:______________
C
      DO 240 I=1,NPE                                                            
         NR=(NOD(N,I)-1)*NDF                                                       
         DO 240 II=1,NDF                                                           
            NR=NR+1                                                                   
            L=(I-1)*NDF+II                                                            
            GLK(NR,NBW)=GLK(NR,NBW)+ELF(L)                                                    
            DO 220 J=1,NPE                                                            
               NCL=(NOD(N,J)-1)*NDF                                                      
               DO 220 JJ=1,NDF                                                           
                  M=(J-1)*NDF+JJ                                                            
                  NC=NCL-NR+JJ+NHBW                                                            
                  IF(NC.GT.0)THEN                                                       
                     GLK(NR,NC)=GLK(NR,NC)+ELK(L,M)
	            ENDIF
  220       CONTINUE                                                                  
  240 CONTINUE                                                                  
  250 CONTINUE
      
      
C
C     Print global matrices when NPRNT > 2
C
      IF(NPRNT.GE.2) THEN
         WRITE(IT,640)
         DO 260 I=1,NEQ
  260    WRITE(IT,930) (GLK(I,J),J=1,NBW-1) ! Global coefficient matrix
         WRITE(IT,650)
         WRITE(IT,930) (GLK(I,NBW),I=1,NEQ) ! Global source vector
      ENDIF
      
      IF(ITER.EQ.1)THEN
          IF(NONLIN.GT.1)THEN ! Setting the specified values to zero after
              DO I=1,NSPV      ! the first iteration of Newton's method
	         VSPV(I)=0.0
              END DO
              DO NB=1,NSSV
               VSSV0(NB)=VSSV(NB)*FS
              ENDDO
	      ENDIF
      ENDIF
C
C     Impose BOUNDARY CONDITIONS on primary and secondary variables
C
      CALL BNDRYUNSYM(NONLIN,MAXNEQ,MAXFBW,MAXSPV,MAXSSV,NDF,NHBW,GLK,
	1                NSPV,NSSV,ISPV,ISSV,VSPV,VSSV0)    
C
C     Call subroutine SOLVRUNS to solve unsymmetric system of linear equations
C     The solution is returned in the array GLK(I,NBW)
C
C      IRES=0
      CALL SOLVRUNS(GLK,MAXNEQ,MAXFBW,NEQ,NHBW)     
      
   
C
C     Save the previous iteration solution and update the current one
C
         DO I=1,NNM
         L = (I-1)*NDF +1
C         GPU(L) = GPU(L) + GLK(L,NBW)
C         GPU(L+1) = GPU(L+1) + GLK(L+1,NBW)
C         GCU(L) = GCU(L) + GLK(L,NBW)
C         GCU(L+1) = GCU(L+1) + GLK(L+1,NBW)
         IF (LFORM.GT.1) THEN
              GLXY(I,1) = GLXY(I,1) + GLK(L,NBW)
              GLXY(I,2) = GLXY(I,2) + GLK(L+1,NBW)
         ENDIF
         ENDDO

      IF(NONLIN.GT.0)THEN
         DO 270 I=1,NEQ
              GPU(I) = GCU(I)
C            !GPU=previous iteration solution; GCU=current solution
            IF(NONLIN.EQ.1)THEN
               GCU(I)=GLK(I,NBW)! GLK(I,NBW) is the full solution
            ELSE
               GCU(I)=GCU(I)+GLK(I,NBW) ! GLK(I,NBW) is the solution increment
            ENDIF
  270    CONTINUE
         
         
         IF(NEWTON.EQ.0)THEN
         IF(ITER.EQ.1)THEN
C		  IF(NONLIN.GT.1)THEN ! Setting the specified values to zero after
C	         DO I=1,NSPV      ! the first iteration of Newton's method
C	         VSPV(I)=0.0
C	         END DO
C             DO NB=1,NSSV
C	             VSSV0(NB)=VSSV(NB)*FS
C	         ENDDO
C	      ENDIF
            WRITE(IT,1000)
            WRITE(IT,1175)FS
            WRITE(IT,940)
            WRITE(IT,670)
            WRITE(IT,940)
            IF(NDF.EQ.1)THEN
            DO I=1,NNM
               WRITE(IT,950)I,(GLXY(I,J),J=1,2),GCU(I)
            ENDDO
            ELSEIF(NDF.EQ.2)THEN
            DO I=1,NNM
            I1 = (I-1)*NDF + 1
            I2 = I*NDF
               WRITE(IT,950)I,(GLXY(I,J),J=1,2),(GCU(K),K=I1,I2)
               
              IF(I.EQ.INDEX)THEN
	             WRITE(IP,950)I,(GLXY(I,J),J=1,2),(GCU(K),K=I1,I2),FS
              ENDIF
            ENDDO           
            WRITE(IT,940)
            ENDIF
            ENDIF
            
            IF (DABS(FS).LE.NLS) GOTO 175
            ENDIF
      
         
C         PRINT*,GPU(1),GLU(1)
C
C     Test for the convergence of the solution 
C
         IF(NONLIN.GT.0)THEN
            DNORM = 0.0
            DINORM = 0.0
            DO 280 I=1,NEQ
               DNORM = DNORM + GCU(I)*GCU(I)
               IF(NONLIN.EQ.1)THEN
                  DINORM = DINORM + GLK(I,NBW)*GLK(I,NBW)
               ELSE
                  DINORM = DINORM+GLK(I,NBW)*GLK(I,NBW)
               ENDIF
  280       CONTINUE
            
            ERROR=DSQRT(DINORM/DNORM)         
            WRITE(IT,1100)ITER,ERROR
	      IF(ERROR.GT.EPS)THEN
              GOTO 160
            ENDIF
         ENDIF
        ENDIF
C
C     Print the solution (i.e., nodal values of the primary variables)
C
      IF(NDF.EQ.1)THEN
      WRITE(IT,660)
      WRITE(IT,940)
      WRITE(IT,670)
      WRITE(IT,940)
      DO 300 I=1,NNM
  300 WRITE(IT,950)I,(GLXY(I,J),J=1,2),GCU(I)
      WRITE(IT,940)
      
      ELSEIF(NDF.EQ.2)THEN
      WRITE(IT,660)
      WRITE(IT,1175)FS
      WRITE(IT,940)
      WRITE(IT,680)
      WRITE(IT,940)
      DO 301 I=1,NNM
      I1 = (I-1)*NDF+1
      I2 = I*NDF
      WRITE(IT,950)I,(GLXY(I,J),J=1,2),(GCU(K),K=I1,I2)
      IF(I.EQ.INDEX)THEN
	     WRITE(IP,950)I,(GLXY(I,J),J=1,2),(GCU(K),K=I1,I2),FS
      ENDIF
  301 CONTINUE
      WRITE(IT,940)
      ENDIF   
        
        
        
C
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C     *                                                             *
C     *             P O S T P R O C E S S O R   U N I T             *
C     *                                                             *
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
      IF(IGRAD.NE.0)THEN   
         WRITE(IT,970)
C         WRITE(IT,730)
         IF(IGRAD.EQ.1) THEN
            WRITE(6,735)
         ELSE
            WRITE(6,750)
         ENDIF
         WRITE(IT,970)
C
C     Compute the GRADIENT of the solution
C
         DO 320 N=1,NEM
         LL = 0
         DO 310 I=1,NPE
         NI=NOD(N,I)
         ELXY(I,1)=GLXY(NI,1)
         ELXY(I,2)=GLXY(NI,2)
         ELXYP(I,1)=GLXYP(NI,1)
         ELXYP(I,2)=GLXYP(NI,2)
         IF(NDF.EQ.1)THEN
         ELU(I)=GCU(NI)
         ELSE
              ND=(NI-1)*NDF
		    DO 96 J=1,NDF
		        ND = ND+1
	            LL = LL+1
   96	            ELU(LL)=GCU(ND)
         ENDIF
  310    CONTINUE
  320    CALL POSTPROC2D(ELXY,IGRAD,NDF,NPE,ELU,NGPR,NONLIN,N,ELXYP)
         WRITE(IT,970)
      ENDIF
      
      
  175 CONTINUE
      PAUSE
C
C                       F  O  R  M  A  T  S
C
  400 FORMAT(20A4)
  420 FORMAT (/,5X,'COEFFICIENTS OF THE DIFFERENTIAL EQUATION:',//,
     *        8X,'Coefficient, A10 ........................=',E12.4,/,
     *        8X,'Coefficient, A1X ........................=',E12.4,/,
     *        8X,'Coefficient, A1Y ........................=',E12.4,/,
     *        8X,'Coefficient, A20 ........................=',E12.4,/,
     *        8X,'Coefficient, A2X ........................=',E12.4,/,
     *        8X,'Coefficient, A2Y ........................=',E12.4,/,
     *        8X,'Coefficient, A00 ........................=',E12.4,/,
     *        8X,'Coefficient, A0X ........................=',E12.4,/,
     *        8X,'Coefficient, A0Y ........................=',E12.4,/)
  430 FORMAT (/,5X,'CONTINUOUS SOURCE COEFFICIENTS:',//,
     *        8X,'Coefficient, F0  ........................=',E12.4,/,
     *        8X,'Coefficient, FX  ........................=',E12.4,/,
     *        8X,'Coefficient, FY  ........................=',E12.4,/)
  440 FORMAT (/,5X,'CONVECTIVE HEAT TRANSFER DATA:',//,
     *        8X,'Number of elements with convection, NBE .=',I4,/,
     *        8X,'Elements, their LOCAL nodes and convective',/,
     *        8X,'heat transfer data:',/,
     *        8X,'Ele. No.',4X,'End Nodes',8X,'Film Coeff.',6X,
     *           'T-Infinity',/)
  480 FORMAT (/,5X,'NUMERICAL INTEGRATION DATA (Gauss Points):',//,
     *        8X,'Full quadrature (NGPF x NGPF) rule, NGPF =',I4,/,
     *        8X,'Reduced quadrature (NGPR x NGPR),   NGPR =',I4,/)
  530 FORMAT (/,5X,'PARAMETERS OF THE DYNAMIC ANALYSIS:',//,
     *8X,'Coefficient, C0 .........................=',E12.4,/,
     *8X,'Coefficient, CX .........................=',E12.4,/,
     *8X,'Coefficient, CY .........................=',E12.4)
  560 FORMAT (/,5X,'FINITE ELEMENT MESH INFORMATION:',//,
     *        8X,'Element type: 1 = LINEAR; 2=QUADRATIC)...=',I4,/,
     *        8X,'Number of nodes per element, NPE ........=',I4,/,
     *        8X,'No. of primary deg. of freedom/node, NDF =',I4,/,
     *        8X,'Number of elements in the mesh, NEM .....=',I4,/,
     *        8X,'Number of nodes in the mesh, NNM ........=',I4,/,
     *        8X,'Number of equations to be solved, NEQ ...=',I4,/,
     *        8X,'Half bandwidth of the matrix GLK, NHBW ..=',I4)
  570 FORMAT (8X,'Mesh subdivisions, NX and NY ............=',2I4,/)
  580 FORMAT (5X,'Node   x-coord.   y-coord.    Speci. primary & seconda
     *ry variables',/,38X,'(0, unspecified; >0, specified)',
     *               /,41X,'Primary DOF  Secondary DOF')
  600 FORMAT(/,1X,'TIME =',E12.4,5X,'Time step number =',I4,/) 
  610 FORMAT (/,5X,'Element coefficient matrix: ',/)
  620 FORMAT (/,5X,'Element mass matrix: ',/)
  630 FORMAT (/,5X,'Element source vector:',/)
  640 FORMAT (/,5X,'Global coefficient matrix:',/)
  650 FORMAT (/,5X,'Global source vector:',/)
  660 FORMAT (/,5X,'S O L U T I O N :',/)
  670 FORMAT (5X,'Node    x-coord.      y-coord.     SOLUTION  ')
  680 FORMAT (5X,'Node    x-coord.      y-coord.     U      V  ')
  690 FORMAT (/,5X,'Connectivity Matrix, [NOD]',/)
  700 FORMAT (8X,'No. of specified PRIMARY variables, NSPV =',I4)
  710 FORMAT (8X,'No. of speci.  SECONDARY variables, NSSV =',I4,/)
  720 FORMAT (6X,'Node  DOF     Value',/)
  730 FORMAT (4X,'The orientation of  gradient vector is measured from
     1 the positive x-axis',/)
  735 FORMAT (12X,'x-coord.     y-coord.     SigmaXX     SigmaYY',
     1        5X,'SigmaXY',/ 40X,'SXX         SYY          SXY')
  740 FORMAT (4X,'x-coord.     y-coord.   -a11(du/dx)  -a22(du/dy)',
     1        3X,'Flux Mgntd  Orientation')
  750 FORMAT (4X,'x-coord.     y-coord.    a22(du/dy)  -a11(du/dx)',
     1        3X,'Flux Mgntd  Orientation')
  800 FORMAT (/,8X,'*** A mesh of QUADRILATERALS is chosen by user ***')
  810 FORMAT (/,8X,'*** A mesh of TRIANGLES is chosen by user ***')
  820 FORMAT (/,8X,'*** The NEWTON iteration scheme is used ***')
  821 FORMAT (/,8X,'*** A TRANSIENT PROBLEM IS ANALYZED ***')
  830 FORMAT (/,8X,'*** The DIRECT iteration scheme is used ***')
  860 FORMAT (8X,I5,5X,2I5,6X,E13.5,5X,E13.5)
  870 FORMAT (5X,I3,2E12.4,8X,I9,9X,I5)
  880 FORMAT (5X,I3,2E12.4,7X,3I4,2X,3I4)
  890 FORMAT (12X,'OUTPUT FROM PROGRAM  *FEM2DUNSYM*  BY J.N. REDDY')
  900 FORMAT (10X,10I5)
  910 FORMAT (2X,70('_'),/)
  920 FORMAT (5X,I3,2E12.4,8X,2I5,4X,2I5)
  930 FORMAT (8X,5E14.5)
  940 FORMAT (2X,50('_'),/)
  950 FORMAT (5X,I3,5E14.5)
  955 FORMAT (1X,5E14.5)
  960 FORMAT (5X,I5,I4,E14.5)
  970 FORMAT (2X,77('_'),/)
  980 FORMAT (/,5X,'NONLINEAR COEFFICIENTS OF THE EQUATION:',//,
     *        8X,'Coefficient, A1U  .......................=',E12.4,/,
     *        8X,'Coefficient, A1UX .......................=',E12.4,/,
     *        8X,'Coefficient, A1UY .......................=',E12.4,/,
     *        8X,'Coefficient, A2U  .......................=',E12.4,/,
     *        8X,'Coefficient, A2UX .......................=',E12.4,/,
     *        8X,'Coefficient, A2UY .......................=',E12.4,/)
  985 FORMAT (/,5X,'VISCOUS FLOW DATA:',//,
     *        8X,'Flow viscousity, AMU ....................=',E12.4,/,
     *        8X,'Penalty parameter, GAMA .................=',E12.4,/)
  990 FORMAT (/,5X,'NONLINEAR ANALYSIS DATA:',//,
     *        8X,'Number of load steps, NLS ...............=',I4,/,
     *        8X,'Max. no. of iterations allowed, ITMAX ...=',I4,/,
     *        8X,'Type of iter. method used (1: DI; >1 NR).=',I4,/,
     *        8X,'Error tolerance, EPS ....................=',E12.4,/,
     *        8X,'Acceleration parameter, REX .............=',E12.4,/)
 1000 FORMAT (/,5X,'LINEAR SOLUTION at the first iteration:',/)
 1050 FORMAT (/,5X,'NONLINEAR SOLUTION at ITMAX:',I3,/)
 1175 FORMAT (/,5X,'LOAD, FS................................=',E12.4)      
 1100 FORMAT (/,5X,'Iteration No.',I3,
     *          5X,'Error (root-mean-square value)=',E12.4)
 1101 FORMAT (/,5X,'Time, TIME ............................ =',E12.4,/,
     *          5X,'Converged NONLINEAR SOLUTION at iteration',I3,/,  
     *          5X,'Error (root-mean-square value)..........=',E12.4)
 1200 FORMAT (/,5X,'Convergence is NOT ACHIEVED with ITMAX =',I3)
C
      END

	SUBROUTINE BNDRYUNSYM(NONLIN,MAXNEQ,MAXFBW,MAXSPV,MAXSSV,NDF,
	1           NHBW,GLK,NSPV,NSSV,ISPV,ISSV,VSPV,VSSV)    
C     __________________________________________________________________        
C                                                                               
C     The subroutine is used to implement specified boundary conditions         
C           on BANDED SYMMETRIC system of finite element equations               
C     __________________________________________________________________        
C                                                                               
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
      DIMENSION ISPV(MAXSPV,2),ISSV(MAXSSV,2),VSPV(MAXSPV),VSSV(MAXSSV),
     1          GLK(MAXNEQ,MAXFBW)
C
      NBW=2*NHBW
      
C                                                                               
C     Modify the source vector to include specified nonzero SECONDARY VARIABLES             
C                                                                               
      IF(NSSV.NE.0)THEN                                                        
         DO 130 NS=1,NSSV
	   NB=(ISSV(NS,1)-1)*NDF+ISSV(NS,2)
C         DO 125 J=1,NBW                                                            
C  125      GLK(NB,J)=0.0D0                                                         
C         GLK(NB,NHBW)=1.0D0         
  130    GLK(NB,NBW)=GLK(NB,NBW)+VSSV(NS)                                         
      ENDIF 
      
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
                                                                    
	RETURN
	END
      
      
      
      

      SUBROUTINE ELMATRCS2D(ICONV,IEL,NPE,NN,NGPF,NGPR,NONLIN,
     *    NOD,GLXY,NBE,IBN,INOD,BETA,UREF,MAXELM,MAXNOD,MAXCNV,
     *    MAXSPV,MAXSSV,LFORM,THIKNS,F0)
C     ________________________________________________________________
C
C     Element calculations based on linear and quadratic rectangular
C                elements with isoparametric formulation.
C     ________________________________________________________________
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/STF/ELF(18),ELK(18,18),ELXY(9,2),ELU(18)
C      COMMON/PST/A10,A1X,A1Y,A20,A2X,A2Y,A00,A0X,A0Y,F0,FX,FY,
C     *           A1U,A1UX,A1UY,A2U,A2UX,A2UY,C0,CX,CY,A1,A2,A3,A4,A5
      COMMON/SHP/SFL(9),GDSFL(2,9)
C      COMMON/SOURCE/FX0,FY0,FXX,FXY,FYX,FYY
	DIMENSION IBN(MAXCNV),INOD(MAXCNV,3),BETA(MAXCNV),UREF(MAXCNV)
      DIMENSION NOD(MAXELM,9),GLXY(MAXNOD,2)
      DIMENSION GAUSPT(5,5),GAUSWT(5,5),TANG(18,18)
      COMMON/CNTC/C(3,3)
      COMMON/IO/IN,IT,IP,IE,IS
C
      DATA GAUSPT/5*0.0D0, -0.57735027D0, 0.57735027D0, 3*0.0D0,
     2  -0.77459667D0, 0.0D0, 0.77459667D0, 2*0.0D0, -0.86113631D0,
     3  -0.33998104D0, 0.33998104D0, 0.86113631D0, 0.0D0, -0.90617984D0,
     4  -0.53846931D0,0.0D0,0.53846931D0,0.90617984D0/
C
      DATA GAUSWT/2.0D0, 4*0.0D0, 2*1.0D0, 3*0.0D0, 0.55555555D0,
     2   0.88888888D0, 0.55555555D0, 2*0.0D0, 0.34785485D0,
     3 2*0.65214515D0, 0.34785485D0, 0.0D0, 0.23692688D0,
     4   0.47862867D0, 0.56888888D0, 0.47862867D0, 0.23692688D0/
C
C     Initialize the arrays
C
      NDF = NN/NPE

      DO 100 I = 1,NN
         ELF(I)  = 0.0
      DO 100 J = 1,NN
         IF(NONLIN.GT.1)THEN
            TANG(I,J)=0.0
         ENDIF
  100    ELK(I,J)= 0.0
C
C     Do-loops on numerical (Gauss) integration begin here. Subroutine
C     SHPRCT (SHaPe functions for ReCTangular elements) is called here
C
      DO 200 NI = 1,NGPF
      DO 200 NJ = 1,NGPF
         XI  = GAUSPT(NI,NGPF)
         ETA = GAUSPT(NJ,NGPF)
         CALL INTERPLN2D(NPE,XI,ETA,DET,ELXY)
         CNST = DET*GAUSWT(NI,NGPF)*GAUSWT(NJ,NGPF)*THIKNS
C
C     Compute x, y, U, UX, UY needed
C
C*        ***** YOUR STATEMENTS *****
         X = 0.0
         Y = 0.0
         U = 0.0
         UX = 0.0
         UY = 0.0
         VX = 0.0
         VY = 0.0
         DO 20 I=1,NPE
              L = (I-1)*NDF + 1
              X = X + SFL(I)*ELXY(I,1)
              Y = Y + SFL(I)*ELXY(I,2)
              U = U + SFL(I)*ELU(L)
              UX = UX + GDSFL(1,I)*ELU(L)
              UY = UY + GDSFL(2,I)*ELU(L)
              VX = VX + GDSFL(1,I)*ELU(L+1)
20            VY = VY + GDSFL(2,I)*ELU(L+1)
          
          UX2 = UX*UX
          UY2 = UY*UY
          VX2 = VX*VX
          VY2 = VY*VY
          
          IF(LFORM.EQ.1)THEN
          UXP1 = 1.0 + UX
          UXP2 = UXP1*UXP1
          VYP1 = 1.0 + VY
          VYP2 = VYP1*VYP1
          
C     define the green strain and second piola-kirchoff stress
          EXX = UX + 0.5*(UX2 + VX2)
          EYY = VY + 0.5*(UY2 + VY2)
          EXY = UY + VX + UX*UY + VX*VY
          S11 = C(1,1)*EXX + C(1,2)*EYY
          S22 = C(1,2)*EXX + C(2,2)*EYY
          S12 = C(3,3)*EXY
          
          ELSE
          
C     define the euler-almansi strain and cauchy stress tensor         
          EXX = UX - 0.5*(UX2 + VX2)
          EYY = VY - 0.5*(UY2 + VY2)
          EXY = UY + VX - UX*UY - VX*VY
          S11 = C(1,1)*EXX + C(1,2)*EYY
          S22 = C(1,2)*EXX + C(2,2)*EYY
          S12 = C(3,3)*EXY
          ENDIF
          
      II = 1
      DO 50 I=1,NPE
C
C     Imbalance force coefficients for the TOTAL Lagrange formulation
C          
          IF(LFORM.EQ.1)THEN
          ELF(II) = ELF(II) + (F0*SFL(I)-UXP1*GDSFL(1,I)*S11-
     *                         UY*GDSFL(2,I)*S22-(UXP1*GDSFL(2,I)+
     *                         UY*GDSFL(1,I))*S12)*CNST
          ELF(II+1) = ELF(II+1) + (F0*SFL(I)-VX*GDSFL(1,I)*S11-
     *                         VYP1*GDSFL(2,I)*S22-(VYP1*GDSFL(1,I)+
     *                         VX*GDSFL(2,I))*S12)*CNST         
          ELSE
C
C     Imbalance force coefficients for the UPDATED Lagrange formulation
C
          ELF(II)=ELF(II)+(F0*SFL(I)-GDSFL(1,I)*S11-GDSFL(2,I)*S12)*CNST
          ELF(II+1) = ELF(II+1)+(F0*SFL(I)-GDSFL(1,I)*S12-
     *                              GDSFL(2,I)*S22)*CNST
          ENDIF
          
      JJ = 1
      DO 40 J=1,NPE
      SIG = S11*GDSFL(1,I)*GDSFL(1,J)+S22*GDSFL(2,I)*GDSFL(2,J)
     *      + S12*(GDSFL(1,I)*GDSFL(2,J)+GDSFL(2,I)*GDSFL(1,J))
C
C     Stiffness coefficients for the TOTAL Lagrange formulation
C
      IF(LFORM.EQ.1)THEN
      
      ELK(II,JJ) = ELK(II,JJ)+(C(1,1)*UXP2*GDSFL(1,I)*GDSFL(1,J)
     *  +C(2,2)*UY2*GDSFL(2,I)*GDSFL(2,J)
     *  +C(1,2)*UXP1*UY*(GDSFL(1,I)*GDSFL(2,J)+GDSFL(2,I)*GDSFL(1,J))
     *  +C(3,3)*(UXP1*GDSFL(2,I)+UY*GDSFL(1,I))*
     *          (UXP1*GDSFL(2,J)+UY*GDSFL(I,J))+SIG)*CNST
      ELK(II+1,JJ+1) = ELK(II+1,JJ+1)+(C(1,1)*VX2*GDSFL(1,I)
     *    *GDSFL(1,J)+C(2,2)*VYP2*GDSFL(2,I)*GDSFL(2,J)
     *    +C(1,2)*VYP1*VX*(GDSFL(1,I)*GDSFL(2,J)+GDSFL(2,I)*GDSFL(1,J))
     *    +C(3,3)*(VYP1*GDSFL(1,I)+VX*GDSFL(2,I))*
     *            (VYP1*GDSFL(1,J)+VX*GDSFL(2,J))+SIG)*CNST
      ELK(II,JJ+1) = ELK(II,JJ+1)+(C(1,1)*UXP1*VX*GDSFL(1,I)
     *       *GDSFL(1,J)+C(2,2)*VYP1*UY*GDSFL(2,I)*GDSFL(2,J)
     *       +C(1,2)*(UXP1*VYP1*GDSFL(1,I)*GDSFL(2,J)
     *       +UY*VX*GDSFL(2,I)*GDSFL(1,J))
     *       +C(3,3)*(UXP1*GDSFL(2,I)+UY*GDSFL(1,I))*
     *       (VYP1*GDSFL(1,J)+VX*GDSFL(2,J)))*CNST
      ELK(II+1,JJ) = ELK(II+1,JJ) +(C(1,1)*UXP1*VX*GDSFL(1,J)
     *     *GDSFL(1,I)+C(2,2)*VYP1*UY*GDSFL(2,J)*GDSFL(2,I)
     *     +C(1,2)*(UXP1*VYP1*GDSFL(1,J)*GDSFL(2,I)
     *     +UY*VX*GDSFL(2,J)*GDSFL(1,I))
     *     +C(3,3)*(UXP1*GDSFL(2,J)+UY*GDSFL(1,J))*
     *             (VYP1*GDSFL(1,I)+VX*GDSFL(2,I)))*CNST
      
      ELSE
C
C     Stiffness coefficients for the UPDATED Lagrange formulation
C     
      ELK(II,JJ) = ELK(II,JJ)+(C(1,1)*GDSFL(1,I)*GDSFL(1,J)
     *                 +C(3,3)*GDSFL(2,I)*GDSFL(2,J)+SIG)*CNST
      ELK(II+1,JJ+1) = ELK(II+1,JJ+1)+(C(3,3)*GDSFL(1,I)
     *      *GDSFL(1,J)+C(2,2)*GDSFL(2,I)*GDSFL(2,J)+SIG)*CNST
      ELK(II,JJ+1) = ELK(II,JJ+1)+(C(1,2)*GDSFL(1,I)*GDSFL(2,J)
     *      +C(3,3)*GDSFL(2,I)*GDSFL(1,J))*CNST
      ELK(II+1,JJ) = ELK(II+1,JJ)+(C(1,2)*GDSFL(2,I)*GDSFL(1,J)
     *      +C(3,3)*GDSFL(1,I)*GDSFL(2,J))*CNST
      
      ENDIF
40    JJ = NDF*J + 1
50    II = NDF*I + 1
200   CONTINUE
C
      IF(ICONV.NE.0) THEN
C
C     Add the convective terms to ELK and ELF for CONVECTION type 
C     boundary conditions (exact only for straight sided elements)
C
          DO 211 M = 1,NBE
          IF(IBN(M).EQ.N) THEN
              M1  = INOD(M,1)
              M2  = INOD(M,2)
              NM1 = NOD(N,M1)
              NM2 = NOD(N,M2)
              DL  = DSQRT((GLXY(NM2,1)-GLXY(NM1,1))**2
     *              +(GLXY(NM2,2)-GLXY(NM1,2))**2)
              BL  = BETA(M)*DL
              UF  = UREF(M)*BL
              IF(IEL.EQ.1)THEN
                  ELK(M1,M1)=ELK(M1,M1)+BL/3.0
                  ELK(M1,M2)=ELK(M1,M2)+BL/6.0
                  ELK(M2,M1)=ELK(M2,M1)+BL/6.0
                  ELK(M2,M2)=ELK(M2,M2)+BL/3.0
                  ELF(M1)=ELF(M1)+0.5*UF
                  ELF(M2)=ELF(M2)+0.5*UF
              ELSE
                  IF(NPE.GE.8) THEN
                      NPEL=4
                  ELSE
                      NPEL=3
                  ENDIF
                  M3=M1+NPEL
                  ELK(M1,M1)=ELK(M1,M1)+4.0*BL/30.0
                  ELK(M1,M3)=ELK(M1,M3)+2.0*BL/30.0
                  ELK(M1,M2)=ELK(M1,M2)-BL/30.0
                  ELK(M3,M1)=ELK(M3,M1)+2.0*BL/30.0
                  ELK(M3,M3)=ELK(M3,M3)+16.0*BL/30.0
                  ELK(M2,M3)=ELK(M2,M3)+2.0*BL/30.0
                  ELK(M2,M1)=ELK(M2,M1)-BL/30.0
                  ELK(M3,M2)=ELK(M3,M2)+2.0*BL/30.0
                  ELK(M2,M2)=ELK(M2,M2)+4.0*BL/30.0
                  ELF(M1)=ELF(M1)+UF/6.0
                  ELF(M3)=ELF(M3)+4.0*UF/6.0
                  ELF(M2)=ELF(M2)+UF/6.0
              ENDIF
          ENDIF
211       CONTINUE
      ENDIF
C
C     Compute the residual vector {R} and tangent matrix [T]            
C
C      IF(NONLIN.GT.1)THEN
C         DO 221 I=1,NN
C            DO 221 J=1,NN
C  221          ELF(I)=ELF(I)-ELK(I,J)*ELU(J)
C
C         DO 240 I=1,NN
C            DO 240 J=1,NN
C  240          ELK(I,J)=ELK(I,J)+TANG(I,J)
C      ENDIF
C  

      RETURN
      END

      SUBROUTINE MSH2DR(IEL,NX,NY,NPE,NNM,NEM,NOD,DX,DY,X0,Y0,
     1                  GLXY,MAXNEM,MAXNNM,MAXNX,MAXNY)
C     ________________________________________________________________
C
C     The subroutine generates arrays [NOD] and [GLXY] for rectangular
C     domains. The domain is divided  into  NX  subdivisions along the
C     x-direction and NY subdivisions in the y-direction. The subdivi-
C     sions define  rectangular elements  of the type required.  For a
C     triangular element mesh, the subdivision defines two linear ele-
C     ments with their common diagonal being inclined to the right.
C     ________________________________________________________________
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION NOD(MAXNEM,9),GLXY(MAXNNM,2),DX(MAXNX),DY(MAXNY)
      COMMON/IO/IN,IT
C
      NEX1 = NX+1
      NEY1 = NY+1
      NXX  = IEL*NX
      NYY  = IEL*NY
      NXX1 = NXX + 1
      NYY1 = NYY + 1
      NEM  = NX*NY
      NNM=NXX1*NYY1
      IF(NPE.EQ.8)NNM = NXX1*NYY1 - NX*NY
C
C     RECTANGULAR ELEMENTS
C
          K0 = 0
          IF (NPE .EQ. 9) K0=1
C
          NOD(1,1) = 1
          NOD(1,2) = IEL+1
          NOD(1,3) = NXX1+(IEL-1)*NEX1+IEL+1
          IF(NPE .EQ. 9) NOD(1,3)=4*NX+5
          NOD(1,4) = NOD(1,3) - IEL
          IF(NPE .GT. 4) THEN
              NOD(1,5) = 2
              NOD(1,6) = NXX1 + (NPE-6)
              NOD(1,7) = NOD(1,3) - 1
              NOD(1,8) = NXX1+1
              IF(NPE .EQ. 9) THEN
                  NOD(1,9)=NXX1+2
              ENDIF
          ENDIF
C
          IF(NY .GT. 1) THEN
              M = 1
              DO 110 N = 2,NY
              L = (N-1)*NX + 1
              DO 100 I = 1,NPE
  100         NOD(L,I) = NOD(M,I)+NXX1+(IEL-1)*NEX1+K0*NX
  110         M=L
          ENDIF
C
          IF(NX .GT .1) THEN
              DO 140 NI = 2,NX
              DO 120  I = 1,NPE
              K1 = IEL
              IF(I .EQ. 6 .OR. I .EQ. 8)K1=1+K0
  120         NOD(NI,I) = NOD(NI-1,I)+K1
              M = NI
              DO 140 NJ = 2,NY
              L = (NJ-1)*NX+NI
              DO 130 J = 1,NPE
  130         NOD(L,J) = NOD(M,J)+NXX1+(IEL-1)*NEX1+K0*NX
  140         M = L
          ENDIF
C
C     Generate the global coordinates of the nodes, [GLXY]:_____________
C
      DX(NEX1)=0.0
      DY(NEY1)=0.0
      XC=X0
      YC=Y0
      IF(NPE .EQ. 8) THEN
          DO 180 NI = 1, NEY1
          I = (NXX1+NEX1)*(NI-1)+1
          J = 2*NI-1
          GLXY(I,1) = XC
          GLXY(I,2) = YC
          DO 150 NJ = 1,NX
          DELX=0.5*DX(NJ)
          I=I+1
          GLXY(I,1) = GLXY(I-1,1)+DELX
          GLXY(I,2) = YC
          I=I+1
          GLXY(I,1) = GLXY(I-1,1)+DELX
          GLXY(I,2) = YC
  150     CONTINUE
          IF(NI.LE.NY) THEN
              I = I+1
              YC= YC+0.5*DY(NI)
              GLXY(I,1) = XC
              GLXY(I,2) = YC
              DO 160 II = 1, NX
              I = I+1
              GLXY(I,1) = GLXY(I-1,1)+DX(II)
  160         GLXY(I,2) = YC
          ENDIF
  180     YC = YC+0.5*DY(NI)
C
      ELSE
          YC=Y0
          DO 200 NI = 1, NEY1
          XC = X0
          I = NXX1*IEL*(NI-1)
          DO 190 NJ = 1, NEX1
          I=I+1
          GLXY(I,1) = XC
          GLXY(I,2) = YC
          IF(NJ.LT.NEX1) THEN
              IF(IEL.EQ.2) THEN
                  I=I+1
                  XC = XC + 0.5*DX(NJ)
                  GLXY(I,1) = XC
                  GLXY(I,2) = YC
              ENDIF
          ENDIF
  190     XC = XC + DX(NJ)/IEL
          XC = X0
          IF(IEL.EQ.2) THEN
              YC = YC + 0.5*DY(NI)
              DO 195 NJ = 1, NEX1
              I=I+1
              GLXY(I,1) = XC
              GLXY(I,2) = YC
              IF(NJ.LT.NEX1) THEN
                  I=I+1
                  XC = XC + 0.5*DX(NJ)
                  GLXY(I,1) = XC
                  GLXY(I,2) = YC
              ENDIF
  195     XC = XC + 0.5*DX(NJ)
          ENDIF
  200     YC = YC + DY(NI)/IEL
      ENDIF
      RETURN
      END


      SUBROUTINE POSTPROC2D(ELXY,IGRAD,NDF,NPE,ELU,NGPR,NONLIN,N,ELXYP)
C     ________________________________________________________________
C
C     Computes the derivatives of the solution for  heat transfer like
C     problems and stresses for fluid flow, plane elasticity and plate
C     bending problems.
C     ________________________________________________________________
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION  ELXY(9,2),ELU(18),GAUSPT(5,5),ELXYP(9,2)
C      COMMON/PST/A10,A1X,A1Y,A20,A2X,A2Y,A00,A0X,A0Y,F0,FX,FY,
C     *           A1U,A1UX,A1UY,A2U,A2UX,A2UY
      COMMON/SHP/SFL(9),GDSFL(2,9)
      COMMON/CNTC/C(3,3)
      COMMON/IO/IN,IT,IP,IE,IS
C
      DATA GAUSPT/5*0.0D0, -0.57735027D0, 0.57735027D0, 3*0.0D0,
     2  -0.77459667D0, 0.0D0, 0.77459667D0, 2*0.0D0, -0.86113631D0,
     3  -0.33998104D0, 0.33998104D0, 0.86113631D0, 0.0D0, -0.90617984D0,
     4  -0.53846931D0,0.0D0,0.53846931D0,0.90617984D0/
C
      PI=4.0D0*DATAN(1.0D0)
      CONST=180.0D0/PI
C
C     Calculation of the gradient of the solution at the reduced Gauss
C     points of RECTANGULAR ELEMENTS:_________________________________
C
      
      DO 100 NI=1,NGPR
      DO 100 NJ=1,NGPR
         XI  = GAUSPT(NI,NGPR)
         ETA = GAUSPT(NJ,NGPR)
         CALL INTERPLN2D(NPE,XI,ETA,DET,ELXY)
      XC  = 0.0
      YC  = 0.0
      DO 50 I=1,NPE
      XC  = XC + SFL(I)*ELXY(I,1)
   50 YC  = YC + SFL(I)*ELXY(I,2)
      
C
C*        ***** YOUR STATEMENTS *****
C
         U=0.0
	   UX=0.0
	   UY=0.0
	   VX=0.0
	   VY=0.0
         DO 610 I=1,NPE
		  L = (I-1)*NDF + 1
		  U = U + ELU(L)*SFL(I)
		  UX = UX + ELU(L)*GDSFL(1,I)
		  UY = UY + ELU(L)*GDSFL(2,I)
		  VX = VX + ELU(L+1)*GDSFL(1,I)
  610		  VY = VY + ELU(L+1)*GDSFL(2,I)
         UX2 = UX*UX
	   UY2 = UY*UY
	   VX2 = VX*VX
	   VY2 = VY*VY
C
C     Compute the Almansi (Eular) strain and Cauchy stresses
C
	      CEX = UX - 0.5*(UX2 + VX2)
		  CEY = VY - 0.5*(UY2 + VY2)
	      CEXY = UY + VX - UX*UY - VX*VY
		  CS11 = C(1,1)*CEX + C(1,2)*CEY
		  CS22 = C(1,2)*CEX + C(2,2)*CEY
		  CS12 = C(3,3)*CEXY
            
C
C     Compute the Green-Lagrange strains and second Piola-Kirchhoff stress
C            
      CALL INTERPLN2D(NPE,XI,ETA,DET,ELXYP)
      XCP = 0.0
      YCP = 0.0
      DO 51 I=1,NPE
      XCP = XCP + SFL(I)*ELXYP(I,1)
   51 YCP = YCP + SFL(I)*ELXYP(I,2)
      
         U=0.0
	   UX=0.0
	   UY=0.0
	   VX=0.0
	   VY=0.0
         DO 611 I=1,NPE
		  L = (I-1)*NDF + 1
		  U = U + ELU(L)*SFL(I)
		  UX = UX + ELU(L)*GDSFL(1,I)
		  UY = UY + ELU(L)*GDSFL(2,I)
		  VX = VX + ELU(L+1)*GDSFL(1,I)
  611		  VY = VY + ELU(L+1)*GDSFL(2,I)
         UX2 = UX*UX
	   UY2 = UY*UY
	   VX2 = VX*VX
	   VY2 = VY*VY

            PEX = UX + 0.5*(UX2 + VX2)
		  PEY = VY + 0.5*(UY2 + VY2)
	      PEXY = UY + VX + UX*UY + VX*VY
		  PS11 = C(1,1)*PEX + C(1,2)*PEY
		  PS22 = C(1,2)*PEX + C(2,2)*PEY
		  PS12 = C(3,3)*PEXY
            
	WRITE(IT,600) N,XC,YC,CS11,CS22,CS12,PS11,PS22,PS12

      IF(DABS(XC-0.0).LT.0.5 .AND.DABS(YC-1.0).LT.0.5)THEN
      	WRITE(IE,500) XC,YC,CS11,CS22,CS12
	ENDIF
      
      IF(DABS(XCP-0.0).LT.0.5 .AND.DABS(YCP-1.0).LT.0.5)THEN
          WRITE(IS,500) XCP,YCP,PS11,PS22,PS12
	ENDIF
      
100   CONTINUE
500   FORMAT(10E13.5)
600   FORMAT(I5,3X,5E13.5,/,34X,5E13.5)
      
      RETURN
      END

      
      SUBROUTINE STRS2D(ELXY,IGRAD,NDF,NPE,ELU,NGPR,NONLIN,GAMA,RHO,AMU)
C     ________________________________________________________________
C
C     Computes the derivatives of the solution for  heat transfer like
C     problems and stresses for fluid flow, plane elasticity and plate
C     bending problems.
C     ________________________________________________________________
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION  ELXY(9,2),ELU(18),GAUSPT(5,5)
C      COMMON/PST/A10,A1X,A1Y,A20,A2X,A2Y,A00,A0X,A0Y,F0,FX,FY,
C     *           A1U,A1UX,A1UY,A2U,A2UX,A2UY
      COMMON/SHP/SFL(9),GDSFL(2,9)
      COMMON/IO/IN,IT
C
      DATA GAUSPT/5*0.0D0, -0.57735027D0, 0.57735027D0, 3*0.0D0,
     2  -0.77459667D0, 0.0D0, 0.77459667D0, 2*0.0D0, -0.86113631D0,
     3  -0.33998104D0, 0.33998104D0, 0.86113631D0, 0.0D0, -0.90617984D0,
     4  -0.53846931D0,0.0D0,0.53846931D0,0.90617984D0/
      
C     
C     Calculation of stress and strain-rates in the element
C        
      DO 50 NI=1,NGPR
      DO 50 NJ=1,NGPR
         XI  = GAUSPT(NI,NGPR)
         ETA = GAUSPT(NJ,NGPR)
         CALL INTERPLN2D(NPE,XI,ETA,DET,ELXY)
         
         X = 0.0
         Y = 0.0
         U = 0.0
         UX = 0.0
         UY = 0.0
         VX = 0.0
         VY = 0.0
         DO 30 I=1,NPE
              X = X + ELXY(I,1)*SFL(I)
              Y = Y + ELXY(I,2)*SFL(I)
              L = (I-1)*NDF + 1
              UX = UX + ELU(L)*GDSFL(1,I)
              UY = UY + ELU(L)*GDSFL(2,I)
              VX = VX + ELU(L+1)*GDSFL(1,I)
              VY = VY + ELU(L+1)*GDSFL(2,I)
30       CONTINUE
         DIV = UX + VY
         PRS = -GAMA*DIV
         SX = 2.0*AMU*UX - PRS
         SY = 2.0*AMU*VY - PRS
         SXY = AMU*(UY + VX)
50       WRITE(IT,60) X,Y,SX,SY,SXY,PRS
60       FORMAT(6(2X,E11.4))
            
      RETURN
      END

      SUBROUTINE INTERPLN2D(NPE,XI,ETA,DET,ELXY)
C     ________________________________________________________________
C
C     The subroutine evaluates the interpolation functions (SFL(I)) and
C     their derivatives with respect to global coordinates (GDSFL(I,J))
C     for Lagrange linear & quadratic rectangular elements, using  the
C     isoparametric formulation.
C
C     SFL(I)........Interpolation function for node I of the element
C     DSFL(J,I).....Derivative of SFL(I) with respect to XI if J=1 and
C                  and ETA if J=2
C     GDSFL(J,I)....Derivative of SFL(I) with respect to  X if J=1 and
C                  and  Y  if J=2
C     XNODE(I,J)...J-TH (J=1,2) Coordinate of node I of the element
C     NP(I)........Array of element nodes (used to define SF and DSF)
C     GJ(I,J)......Jacobian matrix
C     GJINV(I,J)...Inverse of the jacobian matrix
C     ________________________________________________________________
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION ELXY(9,2),XNODE(9,2),NP(9),DSFL(2,9),GJ(2,2),GJINV(2,2)
      COMMON/SHP/SFL(9),GDSFL(2,9)
      COMMON/IO/IN,IT
      DATA XNODE/-1.0D0, 2*1.0D0, -1.0D0, 0.0D0, 1.0D0, 0.0D0, -1.0D0,
     *        0.0D0, 2*-1.0D0, 2*1.0D0, -1.0D0, 0.0D0, 1.0D0, 2*0.0D0/
      DATA NP/1,2,3,4,5,7,6,8,9/
C
      FNC(A,B) = A*B
      IF(NPE.EQ.4) THEN
C
C     LINEAR Lagrange interpolation functions for FOUR-NODE element
C
         DO 10 I = 1, NPE
         XP  = XNODE(I,1)
         YP  = XNODE(I,2)
         XI0 = 1.0+XI*XP
         ETA0=1.0+ETA*YP
         SFL(I)   = 0.25*FNC(XI0,ETA0)
         DSFL(1,I)= 0.25*FNC(XP,ETA0)
   10    DSFL(2,I)= 0.25*FNC(YP,XI0)
      ELSE
         IF(NPE.EQ.8) THEN
C
C     QUADRATIC Lagrange interpolation functions for EIGHT-NODE element
C
            DO 20 I = 1, NPE
            NI   = NP(I)
            XP   = XNODE(NI,1)
            YP   = XNODE(NI,2)
            XI0  = 1.0+XI*XP
            ETA0 = 1.0+ETA*YP
            XI1  = 1.0-XI*XI
            ETA1 = 1.0-ETA*ETA
            IF(I.LE.4) THEN
               SFL(NI)    = 0.25*FNC(XI0,ETA0)*(XI*XP+ETA*YP-1.0)
               DSFL(1,NI) = 0.25*FNC(ETA0,XP)*(2.0*XI*XP+ETA*YP)
               DSFL(2,NI) = 0.25*FNC(XI0,YP)*(2.0*ETA*YP+XI*XP)
            ELSE
               IF(I.LE.6) THEN
                  SFL(NI)    = 0.5*FNC(XI1,ETA0)
                  DSFL(1,NI) = -FNC(XI,ETA0)
                  DSFL(2,NI) = 0.5*FNC(YP,XI1)
               ELSE
                  SFL(NI)    = 0.5*FNC(ETA1,XI0)
                  DSFL(1,NI) = 0.5*FNC(XP,ETA1)
                  DSFL(2,NI) = -FNC(ETA,XI0)
               ENDIF
            ENDIF
   20       CONTINUE
         ELSE
C
C     QUADRATIC Lagrange interpolation functions for NINE-NODE element
C
            DO 30 I=1,NPE
            NI   = NP(I)
            XP   = XNODE(NI,1)
            YP   = XNODE(NI,2)
            XI0  = 1.0+XI*XP
            ETA0 = 1.0+ETA*YP
            XI1  = 1.0-XI*XI
            ETA1 = 1.0-ETA*ETA
            XI2  = XP*XI
            ETA2 = YP*ETA
            IF(I .LE. 4) THEN
               SFL(NI)   = 0.25*FNC(XI0,ETA0)*XI2*ETA2
               DSFL(1,NI)= 0.25*XP*FNC(ETA2,ETA0)*(1.0+2.0*XI2)
               DSFL(2,NI)= 0.25*YP*FNC(XI2,XI0)*(1.0+2.0*ETA2)
            ELSE
               IF(I .LE. 6) THEN
                  SFL(NI)    = 0.5*FNC(XI1,ETA0)*ETA2
                  DSFL(1,NI) = -XI*FNC(ETA2,ETA0)
                  DSFL(2,NI) = 0.5*FNC(XI1,YP)*(1.0+2.0*ETA2)
               ELSE
                  IF(I .LE. 8) THEN
                     SFL(NI)    = 0.5*FNC(ETA1,XI0)*XI2
                     DSFL(2,NI) = -ETA*FNC(XI2,XI0)
                     DSFL(1,NI) = 0.5*FNC(ETA1,XP)*(1.0+2.0*XI2)
                  ELSE
                     SFL(NI)    = FNC(XI1,ETA1)
                     DSFL(1,NI) = -2.0*XI*ETA1
                     DSFL(2,NI) = -2.0*ETA*XI1
                  ENDIF
               ENDIF
            ENDIF
   30       CONTINUE
         ENDIF
      ENDIF
C
C     Compute Jacobian matrix [GJ]=[DSFL][ELXY] and its inverse [GJINV]
C
      DO 40 I = 1,2
      DO 40 J = 1,2
      GJ(I,J)  = 0.0
      DO 40 K = 1,NPE
   40 GJ(I,J)  = GJ(I,J) + DSFL(I,K)*ELXY(K,J)
C
      DET = GJ(1,1)*GJ(2,2)-GJ(1,2)*GJ(2,1)
      GJINV(1,1) = GJ(2,2)/DET
      GJINV(2,2) = GJ(1,1)/DET
      GJINV(1,2) = -GJ(1,2)/DET
      GJINV(2,1) = -GJ(2,1)/DET
C
C     Compute the derivatives of the interpolation functions with
C     respect to the global coordinates (x,y): [GDSFL]=[GJINV][DSFL]
C
      DO 50 I  = 1,2
      DO 50 J  = 1,NPE
      GDSFL(I,J) = 0.0
      DO 50 K  = 1, 2
   50 GDSFL(I,J) = GDSFL(I,J) + GJINV(I,K)*DSFL(K,J)
      RETURN
      END


	SUBROUTINE SOLVRUNS(A,NRMAX,NCMAX,N,ITERM)                                 
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
