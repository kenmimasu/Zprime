C ffbar system z-momentum
      real*8 Pzffmax,Pzffmin,Pzffw
      common/ext_Pzff/Pzffmax,Pzffmin,Pzffw
      real*8 xPzff(500),fxPzff(500,20),fxPzfftot(500)
      common/dist_Pzff/xPzff,fxPzff,fxPzfftot
      integer m_Pzff,ndiv_Pzff 
      common/inp_Pzff/m_Pzff
      common/div_Pzff/ndiv_Pzff 
