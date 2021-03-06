      REAL*8 FUNCTION 
     &SQQ_2V_FF_OD(L3,L4,IQ,JF,P1,P2,P3,P4,
     &         MASSF,VMASSES,VWIDTHS,ODWIDTHS,CONT,EW,BSM,BB_INIT)
C  ODWIDTHSS
C FUNCTION GENERATED BY MADGRAPH
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C SELECTING HELICITIES OF FINAL STATE PARTICLES L3,L4 UNLESS L3,L4=0 
C IN WHICH CASE HELICITIES ARE SUMMED OVER
C FOR THE POINT IN PHASE SPACE P1,P2,P3,P4,...
C FOR PROCESS : q q~  -> f f~ (where f is SM fermion)  
C c 
      IMPLICIT NONE
C
C CONSTANTS
C c
      INTEGER    NEXTERNAL,   NCOMB                     
      PARAMETER (NEXTERNAL=4, NCOMB= 16)
C  
C ARGUMENTS 
C c 
      REAL*8 P1(0:3),P2(0:3),P3(0:3),P4(0:3)
      INTEGER*4 L3,L4,IQ,JF,CONT,EW,BSM,BB_INIT
      REAL*8 MASSF,VMASSES(10),VWIDTHS(10),ODWIDTHS(2,2)
C  
C LOCAL VARIABLES 
C c
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T
      REAL*8 QQ_2V_FF_OD
      INTEGER IHEL
      DATA NTRY/0/
C HELICITY COMBINATIONS
      
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
      SQQ_2V_FF_OD = 0d0
      NTRY=NTRY+1
      DO IHEL=1,NCOMB
C HELICITIES.
      IF(((NHEL(3,IHEL).EQ.L3).OR.(L3.EQ.9)).AND.
     &   ((NHEL(4,IHEL).EQ.L4).OR.(L4.EQ.9)))THEN
         T=QQ_2V_FF_OD(IQ,JF,P1,P2,P3,P4,MASSF,VMASSES,VWIDTHS,
     &                 ODWIDTHS,NHEL(1,IHEL),CONT,EW,BSM,BB_INIT)
         SQQ_2V_FF_OD = SQQ_2V_FF_OD + T
      ENDIF
      ENDDO
      SQQ_2V_FF_OD = SQQ_2V_FF_OD /  4D0 
      END
       
      REAL*8 FUNCTION 
     &QQ_2V_FF_OD(IQ,JF,P1,P2,P3,P4,MASSF,
     &            VMASSES,VWIDTHS,ODWIDTHS,NHEL,CONT,EW,BSM,BB_INIT)
C  
C FUNCTION GENERATED BY MADGRAPH
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT IN PHASE SPACE P1,P2,P3,P4,...
C AND HELICITY NHEL(1),NHEL(2),....
C  
C FOR PROCESS : q q~  -> f f~ (where f is SM fermion)  
C c
      IMPLICIT NONE
C  
C CONSTANTS
C c 
      INTEGER    NGRAPHS,    NEIGEN,    NEXTERNAL       
      PARAMETER (NGRAPHS=  4,NEIGEN=  1,NEXTERNAL=4)    
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C c 
      REAL*8 P1(0:3),P2(0:3),P3(0:3),P4(0:3)
      INTEGER IQ,JF,NHEL(NEXTERNAL),CONT,EW,BSM,BB_INIT
      REAL*8 MASSF,VMASSES(10),VWIDTHS(10),ODWIDTHS(2,2)
C  
C LOCAL VARIABLES 
C c 
      INTEGER I,J,K
      REAL*8 EIGEN_VAL(NEIGEN), EIGEN_VEC(NGRAPHS,NEIGEN)
      COMPLEX*16 ZTEMP,ZTEMP3,ZTEMP4
      COMPLEX*16 AMP(NGRAPHS),A_OUT(2,2,2)
      COMPLEX*16 W1(6)  , W2(6)  , W3(6)  , W4(6)  , W5(6)        
      COMPLEX*16 W6(6)  , W7(6)  , W8(6)  , W9(6)  , W10(6)
      COMPLEX*16 GQ_EFF(2,2),GF_EFF(2,2)
      REAL*8 MASSI,M_EFF(2),W_EFF(2)
      REAL*8 G_LEFT(2),G_RIGHT(2)
      REAL*8 GAQ(2),GZQ(2),GAF(2),GZF(2),TEMPF(2),TEMPQ(2)
      REAL*8 GVQ(2,2),GVF(2,2),GZZQ(2,10),GZZF(2,10),MZZ(2)
      REAL*8 QQ_2V_FF_OD_TEMP
C  
C COMMONS
      include 'ewparams.inc'
      include 'ewcoups.inc'
C Z' COUPLINGS
      include 'LRcoups.inc'
      include 'zpparams.inc'

C
C
C COLOR DATA
C c 
      DATA EIGEN_VAL(1  )/       6.6666666666666674D-01 /
      DATA EIGEN_VEC(1  ,1  )/   7.0710678118654746D-01 /
      DATA EIGEN_VEC(2  ,1  )/   7.0710678118654746D-01 /
      DATA EIGEN_VEC(3  ,1  )/   7.0710678118654746D-01 /
      DATA EIGEN_VEC(4  ,1  )/   7.0710678118654746D-01 /
C DUMMY CHIRAL COUPLINGS
      DATA G_LEFT(1) /1D0/             
      DATA G_LEFT(2) /0D0/             
      DATA G_RIGHT(1)/0D0/             
      DATA G_RIGHT(2)/1D0/    
C ----------
C BEGIN CODE
C ----------
C RESET.
      DO I=1,NGRAPHS
      AMP(I)=0.D0
      ENDDO
      DO I=1,2
        DO J=1,2
          DO K=1,2
            A_OUT(i,j,k)=0d0
          ENDDO
        ENDDO
      ENDDO
C PHOTON,Z COUPLINGS (LEFT/- IS ARRAY 1; RIGHT/+ IS ARRAY 2).
      IF(IQ.EQ.1)THEN
C D-QUARK.
        DO I=1,2
          GAQ(I)=GAD(I)
          GZQ(I)=GZD(I)*DgZd(I)
          DO J=1,10
            GZZQ(I,J)=GZZD(I,J)
          END DO
        END DO
      ELSE IF(IQ.EQ.2)THEN
C U-QUARK.
        DO I=1,2
          GAQ(I)=GAU(I)
          GZQ(I)=GZU(I)*DgZu(I)
          DO J=1,10
            GZZQ(I,J)=GZZU(I,J)
          END DO
        END DO
      ELSE IF(IQ.EQ.3)THEN
C S-QUARK.
        DO I=1,2
          GAQ(I)=GAD(I)
          GZQ(I)=GZD(I)*DgZd(I)
          DO J=1,10
            GZZQ(I,J)=GZZS(I,J)
          END DO
        END DO
      ELSE IF(IQ.EQ.4)THEN
C C-QUARK.
        DO I=1,2
          GAQ(I)=GAU(I)
          GZQ(I)=GZU(I)*DgZu(I)
          DO J=1,10
            GZZQ(I,J)=GZZC(I,J)
          END DO
        END DO
      ELSE IF(IQ.EQ.5)THEN
C B-QUARK.
        DO I=1,2
          GAQ(I)=GAD(I)
          GZQ(I)=GZD(I)*DgZb(I)
          DO J=1,10
            GZZQ(I,J)=GZZB(I,J)
          END DO
        END DO
      END IF
      IF(JF.EQ.1)THEN
C D-QUARK.
        DO I=1,2
          GAF(I)=GAD(I)
          GZF(I)=GZD(I)*DgZd(I)
          DO J=1,10
            GZZF(I,J)=GZZD(I,J)
          END DO
        END DO
      ELSE IF(JF.EQ.2)THEN
C U-QUARK.
        DO I=1,2
          GAF(I)=GAU(I)
          GZF(I)=GZU(I)*DgZd(I)
          DO J=1,10
            GZZF(I,J)=GZZU(I,J)
          END DO
        END DO
      ELSE IF(JF.EQ.3)THEN
C LEPTON.
        DO I=1,2
          GAF(I)=GAL(I)
          GZF(I)=GZL(I)*DgZe(I)
          DO J=1,10
            GZZF(I,J)=GZZE(I,J)
          END DO
        END DO
      ELSE IF(JF.EQ.4)THEN
C NEUTRINO.
        DO I=1,2
          GAF(I)=0.D0
          GZF(I)=GZN(I)*DgZn(I)
          DO J=1,10
            GZZF(I,J)=GZZN(I,J)
          END DO
        END DO
      ELSE IF(JF.EQ.5)THEN
C BOTTOM (corrected).
        DO I=1,2
          GAF(I)=GAD(I)
          GZF(I)=GZD(I)*DgZb(I)
          DO J=1,10
            GZZF(I,J)=GZZB(I,J)
          END DO
        END DO
      ELSE IF(JF.EQ.6)THEN
C TOP (corrected).
        DO I=1,2
          GAF(I)=GAU(I)
          GZF(I)=GZU(I)*DgZt(I)
          DO J=1,10
            GZZF(I,J)=GZZT(I,J)
          END DO
        END DO
      END IF      
      DO I=1,2
        MZZ(I)=VMASSES(I)
        DO J=1,2
            GVQ(I,J)=GZZQ(I,J)
            GVF(I,J)=GZZF(I,J)
        ENDDO
      ENDDO
C DIAGONALISE BASIS TO GET P-DEPENDENT COUPLINGS
      CALL DIAG(MZZ,ODWIDTHS,GVQ,GVF,
     &                M_EFF,W_EFF,GQ_EFF,GF_EFF)
CC WAVEFUNCTIONS.
      IF(BB_INIT.eq.1)THEN
        MASSI=RMB
      ELSE
        MASSI=ZERO
      ENDIF
      CALL IXXXXX(P1  ,MASSI, NHEL(1  ), 1,W1  )                       
      CALL OXXXXX(P2  ,MASSI, NHEL(2  ),-1,W2  )                       
      CALL IXXXXX(P3  ,MASSF, NHEL(3  ),-1,W3  )                       
      CALL OXXXXX(P4  ,MASSF, NHEL(4  ), 1,W4  )  
      IF(EW.eq.1)then                     
CC PHOTON DIAGRAM.
      CALL JIOXXX(W1  ,W2  ,GAQ,AMASS,AWIDTH,W5  )  
      CALL IOVXXX(W3  ,W4  ,W5  , GAF, AMP(1  ))
C Z DIAGRAM.                             
      CALL JIOXXX(W1  ,W2  ,GZQ,ZMASS,ZWIDTH,W6  )
      CALL IOVXXX(W3  ,W4  ,W6  , GZF, AMP(2  ))
      ENDIF
C Z' DIAGRAMS.
C LEFT & RIGHT CURRENTS
      IF((CONT.ne.3).and.(BSM.eq.1))then
      CALL JIOXXX(W1  ,W2  ,G_LEFT, M_EFF(1),W_EFF(1),W7 )              
      CALL JIOXXX(W1  ,W2  ,G_RIGHT,M_EFF(1),W_EFF(1),W8 )
      CALL JIOXXX(W1  ,W2  ,G_LEFT, M_EFF(2),W_EFF(2),W9 )              
      CALL JIOXXX(W1  ,W2  ,G_RIGHT,M_EFF(2),W_EFF(2),W10)
C CONTRACT WITH ALL COMBINATIONS OF LEFT & RIGHT VERTICES 
      CALL IOVXXX(W3  ,W4  ,W7   ,G_LEFT, A_OUT(1,1,1))
      CALL IOVXXX(W3  ,W4  ,W7   ,G_RIGHT,A_OUT(1,2,1))
      CALL IOVXXX(W3  ,W4  ,W8   ,G_LEFT, A_OUT(2,1,1))
      CALL IOVXXX(W3  ,W4  ,W8   ,G_RIGHT,A_OUT(2,2,1))
      CALL IOVXXX(W3  ,W4  ,W9   ,G_LEFT, A_OUT(1,1,2))
      CALL IOVXXX(W3  ,W4  ,W9   ,G_RIGHT,A_OUT(1,2,2))
      CALL IOVXXX(W3  ,W4  ,W10  ,G_LEFT, A_OUT(2,1,2))
      CALL IOVXXX(W3  ,W4  ,W10  ,G_RIGHT,A_OUT(2,2,2))
C COMBINE CHIRAL AMPLITUDES & MULTIPLY BY EFFECTIVE COUPLINGS
      DO I=1,2
        DO J=1,2
          AMP(3)=AMP(3)+GQ_EFF(I,1)*A_OUT(I,J,1)*GF_EFF(J,1)
          AMP(4)=AMP(4)+GQ_EFF(I,2)*A_OUT(I,J,2)*GF_EFF(J,2)
        ENDDO
      ENDDO
      ENDIF
C TOTAL M2.
      QQ_2V_FF_OD = 0.D0 
      DO I = 1, NEIGEN
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NGRAPHS
              ZTEMP = ZTEMP + EIGEN_VEC(J,I)*AMP(J)
          ENDDO
          QQ_2V_FF_OD =QQ_2V_FF_OD+ZTEMP*EIGEN_VAL(I)*CONJG(ZTEMP) 
cccccccccccccccccccccccccccccccccccccc
          IF(CONT.ne.0)then
            QQ_2V_FF_OD_TEMP = QQ_2V_FF_OD
            DO J = 1, NGRAPHS
              ZTEMP =  EIGEN_VEC(J,I)*AMP(J)
              IF (J.EQ.3) ZTEMP3=ZTEMP
              IF (J.EQ.4) ZTEMP4=ZTEMP
c SUBTRACT SQUARED PARTS FOR INTERFERENCE ONLY
              QQ_2V_FF_OD= QQ_2V_FF_OD - ZTEMP*EIGEN_VAL(I)*CONJG(ZTEMP)
            ENDDO
c SUBTRACT OD INTERFERENCE AS WELL
            QQ_2V_FF_OD= QQ_2V_FF_OD - ZTEMP3*EIGEN_VAL(I)*CONJG(ZTEMP4)
     &                               - ZTEMP4*EIGEN_VAL(I)*CONJG(ZTEMP3)
c SUBTRACT INTERFERENCE ONLY FROM TOTAL FOR SQUARED ONLY
            IF(CONT.eq.2)then
            QQ_2V_FF_OD = QQ_2V_FF_OD_TEMP - QQ_2V_FF_OD
            ENDIF
          ENDIF
ccccccccccccccccccccccccccccccccccccccc             
      ENDDO
C ADDITIONAL COLOR FOR QUARKS
      IF((JF.LE.2).OR.(JF.GE.5))QQ_2V_FF_OD=QQ_2V_FF_OD*3.D0
      END
