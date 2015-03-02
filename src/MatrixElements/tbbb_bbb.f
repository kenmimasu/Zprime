      REAL*8 FUNCTION TBBB_BBB(LAM1,LAM2,P1, P2, P3, P4)
C  
C FUNCTION GENERATED BY MADGRAPH
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P1,P2,P3,P4,...
C FOR PROCESS : b b~  -> b b~ (via t-channel gluon) 
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NEXTERNAL,   NCOMB                     
      PARAMETER (NEXTERNAL=4, NCOMB= 16)
C  
C ARGUMENTS 
C  
      INTEGER LAM1,LAM2
      REAL*8 P1(0:3),P2(0:3),P3(0:3),P4(0:3)
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T
      REAL*8 BBB_BBB_T
      INTEGER IHEL
      LOGICAL GOODHEL(NCOMB)
      DATA GOODHEL/NCOMB*.FALSE./
      DATA NTRY/0/
      DATA (NHEL(IHEL,  1),IHEL=1,4) / -1, -1, -1, -1/
      DATA (NHEL(IHEL,  2),IHEL=1,4) / -1, -1, -1,  1/
      DATA (NHEL(IHEL,  3),IHEL=1,4) / -1, -1,  1, -1/
      DATA (NHEL(IHEL,  4),IHEL=1,4) / -1, -1,  1,  1/
      DATA (NHEL(IHEL,  5),IHEL=1,4) / -1,  1, -1, -1/
      DATA (NHEL(IHEL,  6),IHEL=1,4) / -1,  1, -1,  1/
      DATA (NHEL(IHEL,  7),IHEL=1,4) / -1,  1,  1, -1/
      DATA (NHEL(IHEL,  8),IHEL=1,4) / -1,  1,  1,  1/
      DATA (NHEL(IHEL,  9),IHEL=1,4) /  1, -1, -1, -1/
      DATA (NHEL(IHEL, 10),IHEL=1,4) /  1, -1, -1,  1/
      DATA (NHEL(IHEL, 11),IHEL=1,4) /  1, -1,  1, -1/
      DATA (NHEL(IHEL, 12),IHEL=1,4) /  1, -1,  1,  1/
      DATA (NHEL(IHEL, 13),IHEL=1,4) /  1,  1, -1, -1/
      DATA (NHEL(IHEL, 14),IHEL=1,4) /  1,  1, -1,  1/
      DATA (NHEL(IHEL, 15),IHEL=1,4) /  1,  1,  1, -1/
      DATA (NHEL(IHEL, 16),IHEL=1,4) /  1,  1,  1,  1/
C ----------
C BEGIN CODE
C ----------
C      COMMON/HELAMP/AMPHEL(-1:1,-1:1) 
C      REAL*8 AMPHEL
C      AMPHEL(-1,-1)=0.D0
C      AMPHEL(-1,+1)=0.D0
C      AMPHEL(+1,-1)=0.D0
C      AMPHEL(+1,+1)=0.D0
      TBBB_BBB = 0d0
      NTRY=NTRY+1
      DO IHEL=1,NCOMB
C          IF (GOODHEL(IHEL) .OR. NTRY .LT. 10) THEN
             T=BBB_BBB_T(LAM1,LAM2,P1, P2, P3, P4,NHEL(1,IHEL)) 
             TBBB_BBB = TBBB_BBB + T
C              IF (T .GT. 0D0 .AND. .NOT. GOODHEL(IHEL)) THEN
C                  GOODHEL(IHEL)=.TRUE.
C                  WRITE(*,*) IHEL,T
C              ENDIF
C          ENDIF
      ENDDO
      TBBB_BBB = TBBB_BBB /  4D0 
      END
       
       
      REAL*8 FUNCTION BBB_BBB_T(LAM1,LAM2,P1, P2, P3, P4,NHEL)
C  
C FUNCTION GENERATED BY MADGRAPH
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT IN PHASE SPACE P1,P2,P3,P4,...
C AND HELICITY NHEL(1),NHEL(2),....
C  
C FOR PROCESS : q q~  -> t t~ (via gluon) 
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN,    NEXTERNAL       
      PARAMETER (NGRAPHS=  1,NEIGEN=  1,NEXTERNAL=4)    
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      INTEGER LAM1,LAM2
      REAL*8 P1(0:3),P2(0:3),P3(0:3),P4(0:3)
      INTEGER NHEL(NEXTERNAL) 
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      REAL*8 EIGEN_VAL(NEIGEN), EIGEN_VEC(NGRAPHS,NEIGEN)
      COMPLEX*16 ZTEMP
      COMPLEX*16 AMP(NGRAPHS)
      COMPLEX*16 W1(6)  , W2(6)  , W3(6)  , W4(6)  , W5(6)  
C  
C GLOBAL VARIABLES
C  
      REAL*8           GG(2), G
      COMMON /COUPQCD/ GG,    G
      REAL*8            FMASS(12), FWIDTH(12)
      COMMON /FERMIONS/ FMASS,     FWIDTH
C  
C COLOR DATA
C  
      DATA EIGEN_VAL(1  )/       2.2222222222222221E-01 /                  
      DATA EIGEN_VEC(1  ,1  )/  -1.0000000000000000E+00 /                  
C ----------
C BEGIN CODE
C ----------
c      print*,"BBB_BBB_T: TB,fmass=",TB,FMASS(13-TB )
C      COMMON/HELAMP/AMPHEL(-1:1,-1:1) 
C      REAL*8 AMPHEL
      IF((NHEL(3).EQ.LAM1).AND.(NHEL(4).EQ.LAM2))THEN
        CONTINUE
      ELSE
        BBB_BBB_T = 0.D0 
        RETURN
      END IF
c should be t-channel g exchange
      CALL IXXXXX(P1  ,FMASS(12 ),NHEL(1  ), 1,W1  )                       
      CALL OXXXXX(P2  ,FMASS(12 ),NHEL(2  ),-1,W2  )                       
      CALL OXXXXX(P3  ,FMASS(12 ),NHEL(3  ), 1,W3  )                       
      CALL IXXXXX(P4  ,FMASS(12 ),NHEL(4  ),-1,W4  )                       
      CALL JIOXXX(W3  ,W2  ,GG,ZERO,ZERO,W5  )                             
      CALL IOVXXX(W4  ,W5  ,W1  ,GG,AMP(1  ))                              
      BBB_BBB_T = 0.D0 
      DO I = 1, NEIGEN
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NGRAPHS
              ZTEMP = ZTEMP + EIGEN_VEC(J,I)*AMP(J)
          ENDDO
          BBB_BBB_T =BBB_BBB_T+ZTEMP*EIGEN_VAL(I)*CONJG(ZTEMP) 
      ENDDO
C      CALL GAUGECHECK(AMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NGRAPHS,NEIGEN)
C SAVE HELICITY AMPLITUDES.
C      IF((NHEL(3).EQ.-1).AND.(NHEL(4).EQ.-1))THEN
C          DO J = 1, NGRAPHS
C              AMPHEL(-1,-1) = AMPHEL(-1,-1) + AMP(J)
C          ENDDO
C      ELSE IF((NHEL(3).EQ.-1).AND.(NHEL(4).EQ.+1))THEN
C          DO J = 1, NGRAPHS
C              AMPHEL(-1,+1) = AMPHEL(-1,+1) + AMP(J)
C          ENDDO
C      ELSE IF((NHEL(3).EQ.+1).AND.(NHEL(4).EQ.-1))THEN
C          DO J = 1, NGRAPHS
C              AMPHEL(+1,-1) = AMPHEL(+1,-1) + AMP(J)
C          ENDDO
C      ELSE IF((NHEL(3).EQ.+1).AND.(NHEL(4).EQ.+1))THEN
C          DO J = 1, NGRAPHS
C              AMPHEL(+1,+1) = AMPHEL(+1,+1) + AMP(J)
C          ENDDO
C      END IF 
      END
