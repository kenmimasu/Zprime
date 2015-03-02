c 'madgraph' sm coupling declarations as in initialize.f
      REAL*8         GW, GWWA, GWWZ
      COMMON /COUP1/ GW, GWWA, GWWZ
      REAL*8         GAL(2),GAU(2),GAD(2),GWF(2)
      COMMON /COUP2A/GAL,   GAU,   GAD,   GWF
      REAL*8         GZN(2),GZL(2),GZU(2),GZD(2),G1(2)
      COMMON /COUP2B/GZN,   GZL,   GZU,   GZD,   G1
      REAL*8         GWWH,GZZH,GHHH,GWWHH,GZZHH,GHHHH
      COMMON /COUP3/ GWWH,GZZH,GHHH,GWWHH,GZZHH,GHHHH
      COMPLEX*16     GH(2,12)
      COMMON /COUP4/ GH
