C ffbar invariant mass
      real*8 rmassmax,rmassmin,rmassw
      common/ext_rmass/rmassmax,rmassmin,rmassw
      real*8 xrmass(500),fxrmass(500,20),fxrmasstot(500)
      common/dist_rmass/xrmass,fxrmass,fxrmasstot
      integer m_rmass,ndiv_rmass 
      common/inp_rmass/m_rmass
      common/div_rmass/ndiv_rmass 
