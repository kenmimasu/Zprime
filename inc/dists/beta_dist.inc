C Lorentz boost
      real*8 betamax,betamin,betaw
      common/ext_beta/betamax,betamin,betaw
      real*8 xbeta(500),fxbeta(500,20),fxbetatot(500)
      common/dist_beta/xbeta,fxbeta,fxbetatot
      integer m_beta,ndiv_beta 
      common/inp_beta/m_beta
      common/div_beta/ndiv_beta 
