      SUBROUTINE INITIALIZE(RMT,GAMT)
C******************************************************************
C     SETS UP MASSES AND COUPLING CONSTANTS 
C******************************************************************
      IMPLICIT NONE
C
      REAL*8 RMT,GAMT
C
C     CONSTANTS
C
      REAL*8     SW2
      PARAMETER (SW2 = .2320D0)
      integer igw

C     MASSES AND WIDTHS OF FERMIONS

      REAL*8     TMASS,      BMASS,    CMASS,    SMASS,    UMASS
      PARAMETER (TMASS=175.D0,BMASS=4.5D0,CMASS=0.D0,
     &                        SMASS=0.D0,UMASS=0.D0)
      REAL*8     TWIDTH,    BWIDTH,    CWIDTH,    SWIDTH,    UWIDTH
      PARAMETER (TWIDTH=0.D0,BWIDTH=0.D0,CWIDTH=0.D0,
     &                         SWIDTH=0.D0,UWIDTH=0.D0)
      REAL*8     DMASS,    EMASS,    MUMASS,    TAUMASS
      PARAMETER (DMASS=0.D0,EMASS=0.D0,MUMASS=0.D0,
     &                                 TAUMASS=1.78D0)
      REAL*8     DWIDTH,    EWIDTH,    MUWIDTH,    TAUWIDTH
      PARAMETER (DWIDTH=0D0,EWIDTH=0D0,MUWIDTH=0D0,TAUWIDTH=0D0)

C     MASSES AND WIDTHS OF BOSONS

      REAL*8     WMASS,      ZMASS,      WWIDTH,     ZWIDTH
      PARAMETER (WMASS=80.23D0, ZMASS=91.19D0, 
     .           WWIDTH=2.08D0, ZWIDTH=2.50D0)
      REAL*8     AMASS,     AWIDTH,     HMASS,        HWIDTH
      PARAMETER (AMASS=0D0, AWIDTH=0D0, HMASS=150.D0, HWIDTH=0.1635D-01)

C
C     LOCAL
C
      INTEGER I
C
C     GLOBAL
C
      REAL*8          GW, GWWA, GWWZ
      COMMON /COUP1/ GW, GWWA, GWWZ
      REAL*8         GAL(2),GAU(2),GAD(2),GWF(2)
      COMMON /COUP2A/GAL,   GAU,   GAD,   GWF
      REAL*8         GZN(2),GZL(2),GZU(2),GZD(2),G1(2)
      COMMON /COUP2B/GZN,   GZL,   GZU,   GZD,   G1
      REAL*8        GWWH,GZZH,GHHH,GWWHH,GZZHH,GHHHH
      COMMON /COUP3/ GWWH,GZZH,GHHH,GWWHH,GZZHH,GHHHH
      COMPLEX*16     GCHF(2,12)
      COMMON /COUP4/ GCHF
      REAL*8         WMASS1,WWIDTH1,ZMASS1,ZWIDTH1
      COMMON /VMASS1/WMASS1,WWIDTH1,ZMASS1,ZWIDTH1
      REAL*8         AMASS1,AWIDTH1,HMASS1,HWIDTH1
      COMMON /VMASS2/AMASS1,AWIDTH1,HMASS1,HWIDTH1
      REAL*8            FMASS(12), FWIDTH(12)
      COMMON /FERMIONS/ FMASS,     FWIDTH
      REAL*8           GG(2), G
      COMMON /COUPQCD/ GG,    G

C-----
C  BEGIN CODE
C-----
      FMASS(1) = EMASS
      FMASS(2) = 0D0
      FMASS(3) = UMASS
      FMASS(4) = DMASS
      FMASS(5) = MUMASS
      FMASS(6) = 0D0
      FMASS(7) = CMASS
      FMASS(8) = SMASS
      FMASS(9) = TAUMASS
      FMASS(10)= 0D0
      FMASS(11)= RMT
      FMASS(12)= BMASS

      FWIDTH(1) = EWIDTH
      FWIDTH(2) = 0D0
      FWIDTH(3) = UWIDTH
      FWIDTH(4) = DWIDTH
      FWIDTH(5) = MUWIDTH
      FWIDTH(6) = 0D0
      FWIDTH(7) = CWIDTH
      FWIDTH(8) = SWIDTH
      FWIDTH(9) = TAUWIDTH
      FWIDTH(10)= 0D0
      FWIDTH(11)= GAMT
      FWIDTH(12)= BWIDTH

C      WMASS1=WMASS
      WMASS1=ZMASS*DSQRT(1.D0-SW2)
      ZMASS1=ZMASS
      AMASS1=AMASS
      HMASS1=HMASS
      WWIDTH1=WWIDTH
      ZWIDTH1=ZWIDTH
      AWIDTH1=AWIDTH
      HWIDTH1=HWIDTH
      CALL COUP1X(SW2,GW,GWWA,GWWZ)
      CALL COUP2X(SW2,GAL,GAU,GAD,GWF,GZN,GZL,GZU,GZD,G1)
      CALL COUP3X(SW2,ZMASS,HMASS,GWWH,GZZH,GHHH,GWWHH,GZZHH,GHHHH)
      DO I=1,12
         CALL COUP4X(SW2,ZMASS,FMASS(I),GCHF(1,I))
      ENDDO

C     QCD COUPLINGS

      G = 1D0
      GG(1)=-G
      GG(2)=-G
      igw=0          !don't include W width effects
c      CALL TOPWID(fmass(11),Wmass,fmass(12),Wwidth,IGW,fwidth(11))
c      call printconstants
      RETURN
      END

      Subroutine printconstants
c*************************************************************************
c     Prints out all masses, widths, and couplings in common blocks
c*************************************************************************
      implicit none
c
c     Local
c
c      integer i
C
C     GLOBAL
C
      REAL*8          GW, GWWA, GWWZ
      COMMON /COUP1/ GW, GWWA, GWWZ
      REAL*8         GAL(2),GAU(2),GAD(2),GWF(2)
      COMMON /COUP2A/GAL,   GAU,   GAD,   GWF
      REAL*8         GZN(2),GZL(2),GZU(2),GZD(2),G1(2)
      COMMON /COUP2B/GZN,   GZL,   GZU,   GZD,   G1
      REAL*8         GWWH,GZZH,GHHH,GWWHH,GZZHH,GHHHH
      COMMON /COUP3/ GWWH,GZZH,GHHH,GWWHH,GZZHH,GHHHH
      COMPLEX*16     GCHF(2,12)
      COMMON /COUP4/ GCHF
      REAL*8         WMASS,WWIDTH,ZMASS,ZWIDTH
      COMMON /VMASS1/WMASS,WWIDTH,ZMASS,ZWIDTH
      REAL*8         AMASS,AWIDTH,HMASS,HWIDTH
      COMMON /VMASS2/AMASS,AWIDTH,HMASS,HWIDTH
      REAL*8            FMASS(12), FWIDTH(12)
      COMMON /FERMIONS/ FMASS,     FWIDTH
      REAL*8           GG(2), G
      COMMON /COUPQCD/ GG,    G

C-----
C  BEGIN CODE
C-----
      write(*,'(a)') 'Boson masses and Widths:'
      write(*,10) 'W',wmass,'W',wwidth,'Z',zmass,'Z',zwidth
      write(*,10) 'A',amass,'A',awidth,'H',hmass,'H',hwidth
      write(*,'(a)') 'Quark masses and Widths:'
      write(*,10) 't',fmass(11),'t',fwidth(11),
     &            'b',fmass(12),'b',fwidth(12)


 10   format(a,1x,5Hmass=,f6.2,1H:,2x,a,1x,6Hwidth=,f6.2,1H:,2x,
     &       a,1x,5Hmass=,f6.2,1H:,2x,a,1x,6Hwidth=,f6.2)

      end

      SUBROUTINE TOPWID(RMT,RMW,RMB,RGW,IGW,RGT)
c*************************************************************************
c     THE TOTAL WEAK DECAY WIDTH OF THE TOP QUARK, INCLUDING
c     THE EFFECTS OF BOTTOM MASS AND, IF IGW=1,  A FINITE W WIDTH.
c     From James Stirling 6-10-94
c*************************************************************************
      IMPLICIT COMPLEX*16(A-H,O-Z)
      REAL*8 RMT,RMB,RMW,XW,XB,RGW,RGT
*
      PI=4.*DATAN(1.D0)
      GF=1.16637D-05
      GW=CDSQRT(RMW**2*GF/DSQRT(2.D0))
*                            FLAVOUR & COLOUR
      XB=RMB/RMT
      XW=RMW/RMT
      IF(IGW.EQ.1) GOTO 10
      IF(RMT.LE.(RMW+RMB)) THEN
          WRITE(6,*)'WARNING: mt < mw + mb !!!!'
          STOP
          ENDIF
      RGT = GF*RMT**3/8D0/PI/DSQRT(2D0) 
     .   * DSQRT( (1D0-(XW+XB)**2)*(1D0-(XW-XB)**2) )
     .   * ( (1D0-XB**2)**2 + (1D0+XB**2)*XW**2 - 2D0*XW**4 )
      RETURN
  10  CONTINUE
      RM=XB**2
      OM=1.+RM-DCMPLX(RMW**2,RMW*RGW)/RMT**2
      Y1=OM+CDSQRT(OM*OM-4.*RM)
      Y0=OM-CDSQRT(OM*OM-4.*RM)
      Z1=2.
      Z0=2.*CDSQRT(RM)
*
      D0=(-Y0**8+3.*Y0**7*RM+3.*Y0**7-8.*Y0**6*RM-12.*Y0**5*RM**
     . 2-12.*Y0**5*RM+96.*Y0**4*RM**2-48.*Y0**3*RM**3-48.*Y0**3*
     . RM**2-128.*Y0**2*RM**3+192.*Y0*RM**4+192.*Y0*RM**3-256.*
     . RM**4)/(24.*Y0**4*(Y1-Y0))
      D1=(-Y1**8+3.*Y1**7*RM+3.*Y1**7-8.*Y1**6*RM-12.*Y1**5*RM**
     . 2-12.*Y1**5*RM+96.*Y1**4*RM**2-48.*Y1**3*RM**3-48.*Y1**3*
     . RM**2-128.*Y1**2*RM**3+192.*Y1*RM**4+192.*Y1*RM**3-256.*
     . RM**4)/(24.*Y1**4*(Y1-Y0))
      A4=(32.*RM**4*(Y1-Y0))/(3.*Y1*Y0*(Y1-Y0))
      A3=(8.*RM**3*(-3.*Y1**2*Y0*RM-3.*Y1**2*Y0+4.*Y1**2*RM+3.*
     . Y1*Y0**2*RM+3.*Y1*Y0**2-4.*Y0**2*RM))/(3.*Y1**2*Y0**2*(Y1
     . -Y0))
      A2=(8.*RM**3*(2.*Y1**3*Y0**2-3.*Y1**3*Y0*RM-3.*Y1**3*Y0+4.
     . *Y1**3*RM-2.*Y1**2*Y0**3+3.*Y1*Y0**3*RM+3.*Y1*Y0**3-4.*Y0
     . **3*RM))/(3.*Y1**3*Y0**3*(Y1-Y0))
      A1=(2.*RM**2*(3.*Y1**4*Y0**3*RM+3.*Y1**4*Y0**3+8.*Y1**4*Y0
     . **2*RM-12.*Y1**4*Y0*RM**2-12.*Y1**4*Y0*RM+16.*Y1**4*RM**2
     . -3.*Y1**3*Y0**4*RM-3.*Y1**3*Y0**4-8.*Y1**2*Y0**4*RM+12.*
     . Y1*Y0**4*RM**2+12.*Y1*Y0**4*RM-16.*Y0**4*RM**2))/(3.*Y1**
     . 4*Y0**4*(Y1-Y0))
      B0=(Y1**3-3.*Y1**2*RM-3.*Y1**2+8.*Y1*RM-Y0**3+3.*Y0**2*RM+
     . 3.*Y0**2-8.*Y0*RM)/(24.*(Y1-Y0))
      B1=(Y1+Y0-3.*RM-3.)/24.
      B2=1./24.
*
      RINT=D0*CDLOG((Z1-Y0)/(Z0-Y0))
     .    -D1*CDLOG((Y1-Z1)/(Y1-Z0))
     .    -A4/3.*(1./Z1**3-1./Z0**3)
     .    -A3/2.*(1./Z1**2-1./Z0**2)
     .    -A2   *(1./Z1   -1./Z0   )
     .    +A1*CDLOG(Z1/Z0)
     .    +B0   *(Z1   -Z0   )
     .    +B1/2.*(Z1**2-Z0**2)
     .    +B2/3.*(Z1**3-Z0**3)
*
      GW4=GW**4
*
* TOTAL WIDTH INCLUDES FLAVOUR & COLOUR FACTORS
      RGT=RMT**3/(RMW*RGW)*GW4/(8.*PI**3)*DIMAG(RINT)
      RGT=9.*RGT
      RETURN
      END
