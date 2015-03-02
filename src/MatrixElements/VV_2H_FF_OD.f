      REAL*8 FUNCTION SVV_2H_F1F2_OD(L3,L4,JF,VMASS,F1MASS,F2MASS,
     &          SMASSES,SWIDTHS,GIN,GOUT,P1, P2, P3, P4)
C  
C FUNCTION GENERATED BY MADGRAPH
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P1,P2,P3,P4,...
C  
C FOR PROCESS : v v  -> f F~  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NEXTERNAL,   NCOMB,NCOMB2,ICOMB,FCOMB                     
      PARAMETER (NEXTERNAL=4, NCOMB= 81)
C  
C ARGUMENTS 
C  
      REAL*8 P1(0:3),P2(0:3),P3(0:3),P4(0:3)                                     
      REAL*8 GIN(2)  , GOUT(2,2)
      REAL*8 VMASS, F1MASS, F2MASS
      REAL*8 SMASSES(2), SWIDTHS(2,2)
      INTEGER L3,L4,JF
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY       
      REAL*8 T
      REAL*8 VV_2H_F1F2_OD                                                              
      INTEGER IHEL
      LOGICAL GOODHEL(NCOMB)
      DATA GOODHEL/NCOMB*.FALSE./
      DATA NTRY/0/
C INCLUDE HELICITY COMBINATIONS DEPENDING ON MASS
      IF((VMASS.EQ.0D0)) THEN
      ICOMB=4
      FCOMB=4
      CALL massless_VV_to_FF(NHEL)
      ELSE
      ICOMB=9
      FCOMB=4
      CALL massive_VV_to_FF(NHEL)
      ENDIF
C NEW NCOMB
      NCOMB2=ICOMB*FCOMB
C ----------
C BEGIN CODE
C ----------
      SVV_2H_F1F2_OD = 0d0
      NTRY=NTRY+1
      DO IHEL=1,NCOMB2
          IF (GOODHEL(IHEL) .OR. NTRY .LT. 10) THEN
C HELICITIES.
          IF(((NHEL(3,IHEL).EQ.L3).OR.(L3.EQ.9)).AND.
     &      ((NHEL(4,IHEL).EQ.L4).OR.(L4.EQ.9)))THEN
                 T=VV_2H_F1F2_OD(JF,VMASS,F1MASS,F2MASS,SMASSES,
     &                    SWIDTHS,GIN,GOUT,P1, P2, P3, P4,NHEL(1,IHEL)) 
                 SVV_2H_F1F2_OD = SVV_2H_F1F2_OD + T
          ENDIF
              IF (T .GT. 0D0 .AND. .NOT. GOODHEL(IHEL)) THEN
                  GOODHEL(IHEL)=.TRUE.
C                   WRITE(*,*) IHEL,T
              ENDIF
         ENDIF
      ENDDO
      SVV_2H_F1F2_OD = SVV_2H_F1F2_OD /  FLOAT(ICOMB)
      END
      
      REAL*8 FUNCTION VV_2H_F1F2_OD(JF,VMASS,F1MASS,F2MASS,SMASSES,
     &          SWIDTHS,GIN,GOUT,P1, P2, P3, P4,NHEL)
C  
C FUNCTION GENERATED BY MADGRAPH
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT IN PHASE SPACE P1,P2,P3,P4,...
C AND HELICITY NHEL(1),NHEL(2),....
C  
C FOR PROCESS : z z  -> f F~  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN,    NEXTERNAL       
      PARAMETER (NGRAPHS=  2,NEIGEN=  1,NEXTERNAL=4)    
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 P1(0:3),P2(0:3),P3(0:3),P4(0:3)                                    
      INTEGER NHEL(NEXTERNAL)                                                    
      REAL*8 GIN(2)  , GOUT(2,2)
      REAL*8 VMASS, F1MASS, F2MASS
      REAL*8 SMASSES(2), SWIDTHS(2,2)
      INTEGER JF
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      REAL*8 EIGEN_VAL(NEIGEN), EIGEN_VEC(NGRAPHS,NEIGEN)
      REAL*8 TEMPGIN(2,2)
      COMPLEX*16 ZTEMP
      COMPLEX*16 AMP(NGRAPHS),A_OUT(2)
      COMPLEX*16 W1(6)  , W2(6)  , W3(6)  , W4(6)  , W5(6)
      COMPLEX*16 GIN_EFF(2,2),GOUT_EFF(2,2),C_EFF(2)
      COMPLEX*16 GOUT_ONE(2),GOUT_TWO(2)
      REAL*8 M_EFF,W_EFF
      COMPLEX*16 G_LEFT(2),G_RIGHT(2)
C COLOR DATA
C  
      DATA EIGEN_VAL(1  )/       1.0000000000000000E+00 /                  
      DATA EIGEN_VEC(1  ,1  )/   1.0000000000000000E+00 /                  
      DATA EIGEN_VEC(2  ,1  )/   1.0000000000000000E+00 /  
C DUMMY CHIRAL COUPLINGS
      DATA G_LEFT(1) /(1D0,0D0)/             
      DATA G_LEFT(2) /(0D0,0D0)/             
      DATA G_RIGHT(1)/(0D0,0D0)/             
      DATA G_RIGHT(2)/(1D0,0D0)/             
C DIAGONALISE BASIS TO GET P-DEPENDENT COUPLINGS
      DO I=1,2
          TEMPGIN(1,I)=GIN(I)
          TEMPGIN(2,I)=0D0
      ENDDO
      CALL DIAG(P1,P2,SMASSES,SWIDTHS,TEMPGIN,GOUT,
     &                M_EFF,W_EFF,GIN_EFF,GOUT_EFF,C_EFF)
      DO I=1,2
          GOUT_ONE(I)=GOUT_EFF(I,1)
          GOUT_TWO(I)=GOUT_EFF(I,2)
      ENDDO
C ----------
C BEGIN CODE
C ----------
      CALL VXXXXX(P1  ,VMASS, NHEL(1  ),-1,W1  )                        
      CALL VXXXXX(P2  ,VMASS, NHEL(2  ),-1,W2  )                        
      CALL OXXXXX(P3  ,F2MASS,NHEL(3  ), 1,W3  )                       
      CALL IXXXXX(P4  ,F1MASS,NHEL(4  ),-1,W4  )    
                         
      CALL HVVXXX(W1  ,W2  ,1D0,M_EFF,W_EFF,W5  ) 
C DUMMY LEFT, RIGHT AMPS              
      CALL IOSXXX(W4  ,W3  ,W5  ,G_LEFT,A_OUT(1))    
      CALL IOSXXX(W4  ,W3  ,W5  ,G_RIGHT,A_OUT(2))    
C      CALL IOSXXX(W4  ,W3  ,W5  ,GOUT_ONE,AMP(1))    
C      CALL IOSXXX(W4  ,W3  ,W5  ,GOUT_TWO,AMP(2))     
      DO I=1,2
          AMP(I)=C_EFF(I)*GIN_EFF(1,I)*
     &          (A_OUT(1)*GOUT_EFF(1,I)+A_OUT(2)*GOUT_EFF(2,I))
      ENDDO 
C       DO I=1,2
C           AMP(I)=C_EFF(I)*GIN_EFF(1,I)*AMP(I)
C       ENDDO 
      VV_2H_F1F2_OD = 0.D0 
      DO I = 1, NEIGEN
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NGRAPHS
              ZTEMP = ZTEMP + EIGEN_VEC(J,I)*AMP(J)
          ENDDO
          VV_2H_F1F2_OD=VV_2H_F1F2_OD+ZTEMP*EIGEN_VAL(I)*CONJG(ZTEMP)
      ENDDO
      IF(JF.LE.2)VV_2H_F1F2_OD=VV_2H_F1F2_OD*3.D0

C      CALL GAUGECHECK(AMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NGRAPHS,NEIGEN)
      END