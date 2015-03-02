      real*8 rmW,gamW
      common/W/rmW,gamW
      real*8 rmZ,gamZ
      common/Z/rmZ,gamZ
      real*8 rmt,gamt
      common/t/rmt,gamt
      real*8 rmb
      common/b/rmb
      real*8 rmh, gamh
      common/H/rmH,gamH
      real*8 a_em,s2w,zc_norm
      common/EW/a_em,s2w,zc_norm
c 'madgraph' sm mass declarations as in initialize.f
      real*8         WMASS,WWIDTH,ZMASS,ZWIDTH
      common /VMASS1/WMASS,WWIDTH,ZMASS,ZWIDTH
      real*8         AMASS,AWIDTH,HMASS,HWIDTH
      common /VMASS2/AMASS,AWIDTH,HMASS,HWIDTH
      real*8            FMASS(12), FWIDTH(12)
      common /FERMIONS/ FMASS,     FWIDTH


