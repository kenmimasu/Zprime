      REAL*8 FUNCTION SV1V1_2H_V2V2_OD(L3,L4,V1MASS,V2MASS,SMASSES,
     &          SWIDTHS,GIN,GOUT,P1, P2, P3, P4)
C  
C FUNCTION GENERATED BY MADGRAPH
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P1,P2,P3,P4,...
C FOR PROCESS : V V  -> S1,S2 -> V V  
C 
C REAL*8 V1MASS : MASS OF I VECTOR
C REAL*8 V1MASS : MASS OF F VECTOR
C REAL*8 GIN(2) : COUPLINGS OF INITIAL STATE TO SCALARS 1,2
C REAL*8 GOUT(2): COUPLINGS OF FINAL STATE TO SCALARS 1,2
     
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NEXTERNAL,   NCOMB                     
      PARAMETER (NEXTERNAL=4, NCOMB= 81)
C  
C ARGUMENTS 
C  
      REAL*8 P1(0:3),P2(0:3),P3(0:3),P4(0:3)
      REAL*8 V1MASS,V2MASS
      REAL*8 SMASSES(2) ,SWIDTHS(2,2)
      REAL*8 GIN(2),GOUT(2) 
      INTEGER L3,L4
      
                                     
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY,NCOMB2,ICOMB,FCOMB                                         
      REAL*8 T
      REAL*8 V1V1_2H_V2V2_OD                                                               
      INTEGER IHEL
      LOGICAL GOODHEL(NCOMB)
      DATA GOODHEL/NCOMB*.FALSE./
      DATA NTRY/0/
C SET HELICITY COMBINATIONS DEPENDING ON MASSES
      IF      ((V1MASS.EQ.0D0).AND.(V2MASS.EQ.0D0)) THEN
          ICOMB=4
          FCOMB=4
          CALL massless_VV_to_massless_VV(NHEL)
      ELSE IF ((V1MASS.EQ.0D0).AND.(V2MASS.NE.0D0)) THEN
          ICOMB=4
          FCOMB=9
          CALL massless_VV_to_massive_VV(NHEL)
      ELSE IF ((V1MASS.NE.0D0).AND.(V2MASS.EQ.0D0)) THEN
          ICOMB=9
          FCOMB=4
          CALL massive_VV_to_massless_VV(NHEL)
      ELSE IF ((V1MASS.NE.0D0).AND.(V2MASS.NE.0D0)) THEN
          ICOMB=9
          FCOMB=9
          CALL massive_VV_to_massive_VV(NHEL)
      ENDIF
C NEW NCOMB
      NCOMB2=ICOMB*FCOMB
C ----------
C BEGIN CODE
C ----------
      SV1V1_2H_V2V2_OD = 0d0
      NTRY=NTRY+1
      DO IHEL=1,NCOMB2
          IF (GOODHEL(IHEL) .OR. NTRY .LT. 10) THEN
C HELICITIES.
          IF(((NHEL(3,IHEL).EQ.L3).OR.(L3.EQ.9)).AND.
     &      ((NHEL(4,IHEL).EQ.L4).OR.(L4.EQ.9)))THEN
             T=V1V1_2H_V2V2_OD(V1MASS,V2MASS,SMASSES,
     &          SWIDTHS,GIN,GOUT,P1, P2, P3, P4, NHEL(1,IHEL))
             SV1V1_2H_V2V2_OD = SV1V1_2H_V2V2_OD + T
          ENDIF
              IF (T .GT. 0D0 .AND. .NOT. GOODHEL(IHEL)) THEN
                  GOODHEL(IHEL)=.TRUE.
C                   WRITE(*,*) IHEL,T
              ENDIF
          ENDIF
      ENDDO
      SV1V1_2H_V2V2_OD = SV1V1_2H_V2V2_OD / FLOAT(ICOMB) 
      END
       
       
      REAL*8 FUNCTION V1V1_2H_V2V2_OD(V1MASS,V2MASS,SMASSES,
     &          SWIDTHS,GIN,GOUT,P1, P2, P3, P4, NHEL)
C  
C FUNCTION GENERATED BY MADGRAPH
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT IN PHASE SPACE P1,P2,P3,P4,...
C AND HELICITY NHEL(1),NHEL(2),....
C  
C FOR PROCESS : z z  -> z z  
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
      REAL*8 V1MASS,V2MASS
      REAL*8 SMASSES(2) ,SWIDTHS(2,2)
      REAL*8 GIN(2),GOUT(2) 
                                      
      INTEGER NHEL(NEXTERNAL)                                                    
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      REAL*8 EIGEN_VAL(NEIGEN), EIGEN_VEC(NGRAPHS,NEIGEN)
      REAL*8 TEMPGIN(2,2),TEMPGOUT(2,2)
      COMPLEX*16 ZTEMP
      COMPLEX*16 AMP(NGRAPHS),A_EFF
      COMPLEX*16 W1(6)  , W2(6)  , W3(6)  , W4(6)  , W5(3)  , W6(3)  
      COMPLEX*16 GIN_EFF(2,2),GOUT_EFF(2,2),C_EFF(2)
      REAL*8 M_EFF,W_EFF
C COLOR DATA
C  
      DATA EIGEN_VAL(1  )/       1.0000000000000000E+00 /                  
      DATA EIGEN_VEC(1  ,1  )/   0.70710678118655E+00 /                  
      DATA EIGEN_VEC(2  ,1  )/   0.70710678118655E+00 /                  
C DIAGONALISE BASIS TO GET P-DEPENDENT COUPLINGS
      DO I=1,2
          TEMPGIN(1,I)=GIN(I)
          TEMPGOUT(1,I)=GOUT(I)
          TEMPGIN(2,I)=0D0
          TEMPGOUT(2,I)=0D0
      ENDDO
      CALL DIAG(P1,P2,SMASSES,SWIDTHS,TEMPGIN,TEMPGOUT,
     &                M_EFF,W_EFF,GIN_EFF,GOUT_EFF,C_EFF)
C ----------
C BEGIN CODE
C ----------
      CALL VXXXXX(P1  ,V1MASS,NHEL(1  ),-1,W1  )                            
      CALL VXXXXX(P2  ,V1MASS,NHEL(2  ),-1,W2  )                            
      CALL VXXXXX(P3  ,V2MASS,NHEL(3  ), 1,W3  )                            
      CALL VXXXXX(P4  ,V2MASS,NHEL(4  ), 1,W4  )
      
      CALL HVVXXX(W1  ,W2  ,1D0,M_EFF,W_EFF,W5  )
      CALL VVSXXX(W4  ,W3  ,W5  ,1D0,A_EFF)
      
      DO I=1,2
      AMP(I)=C_EFF(I)*GIN_EFF(1,I)*A_EFF*GOUT_EFF(1,I)
      ENDDO
      V1V1_2H_V2V2_OD = 0.D0 
      DO I = 1, NEIGEN
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NGRAPHS
              ZTEMP = ZTEMP + EIGEN_VEC(J,I)*AMP(J)
          ENDDO
          V1V1_2H_V2V2_OD =V1V1_2H_V2V2_OD+
     &              ZTEMP*EIGEN_VAL(I)*CONJG(ZTEMP) 
      ENDDO
C      CALL GAUGECHECK(AMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NGRAPHS,NEIGEN)
      END