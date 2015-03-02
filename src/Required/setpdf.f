      subroutine setpdf(x1, x2, QQ, istructure, icoll,
     &                             *,fxn,fx1,fx2)
      implicit none
c Arguments
      integer istructure, icoll, i
      real*8 x1,x2,QQ,Q2,fxn
      real*8 fx1(13),fx2(13)
c Locals
      real*8 u1, usea1, d1, dsea1, str1, chm1, btm1, g1
      real*8 u2, usea2, d2, dsea2, str2, chm2, btm2, g2
c externals
      real*8 ctq6pdf
c construct polarised hadronic structure functions.
      IF(ISTRUCTURE.LE.4)THEN
        Q2=QQ*QQ
        IF((X1.LE.1.D-6).OR.(X1.GE.1.D0))THEN
          FXN=0.D0
          RETURN 1
        END IF
        IF((X2.LE.1.D-6).OR.(X2.GE.1.D0))THEN
          FXN=0.D0
          RETURN 1
        END IF
        IF((QQ.LE.1.3D0).OR.(QQ.GE.1.D4))THEN
          FXN=0.D0
          RETURN 1
        END IF
        U1=X1*Ctq6Pdf(1,X1,QQ)
        D1=X1*Ctq6Pdf(2,X1,QQ)
        USEA1=X1*Ctq6Pdf(-1,X1,QQ)
        DSEA1=X1*Ctq6Pdf(-2,X1,QQ)
        STR1=X1*Ctq6Pdf(3,X1,QQ)
        CHM1=X1*Ctq6Pdf(4,X1,QQ)
        BTM1=X1*Ctq6Pdf(5,X1,QQ)
        G1=X1*Ctq6Pdf(0,X1,QQ)
        U2=X2*Ctq6Pdf(1,X2,QQ)
        D2=X2*Ctq6Pdf(2,X2,QQ)
        USEA2=X2*Ctq6Pdf(-1,X2,QQ)
        DSEA2=X2*Ctq6Pdf(-2,X2,QQ)
        STR2=X2*Ctq6Pdf(3,X2,QQ)
        CHM2=X2*Ctq6Pdf(4,X2,QQ)
        BTM2=X2*Ctq6Pdf(5,X2,QQ)
        G2=X2*Ctq6Pdf(0,X2,QQ)
      END IF
c actual PDFs.
      IF(ISTRUCTURE.LE.4)THEN 
        fx1(1)=D1
        fx1(2)=U1
        fx1(3)=STR1
        fx1(4)=CHM1
        fx1(5)=BTM1
        fx1(6)=0.D0
        fx1(7)=DSEA1
        fx1(8)=USEA1
        fx1(9)=fx1(3)
        fx1(10)=fx1(4)
        fx1(11)=fx1(5)
        fx1(12)=fx1(6)
        fx1(13)=G1
        do i=1,13
          fx1(i)=fx1(i)/x1
        end do
        fx2(1)=D2*(1-icoll)+DSEA2*icoll
        fx2(2)=U2*(1-icoll)+USEA2*icoll
        fx2(3)=STR2
        fx2(4)=CHM2
        fx2(5)=BTM2
        fx2(6)=0.D0
        fx2(7)=D2*icoll+DSEA2*(1-icoll)
        fx2(8)=U2*icoll+USEA2*(1-icoll)
        fx2(9)=fx2(3)
        fx2(10)=fx2(4)
        fx2(11)=fx2(5)
        fx2(12)=fx2(6)
        fx2(13)=G2
        do i=1,13
           fx2(i)=fx2(i)/x2
        end do
      END IF
      END