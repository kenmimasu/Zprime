      real*8 sigpmax,sigpmin,sigpw
      common/ext_sigp/sigpmax,sigpmin,sigpw
      real*8 xsigp(1000),fxsigp(20,1000,20),fxsigptot(20,1000)
      common/dist_sigp/xsigp,fxsigp,fxsigptot
      integer m_sigp, ndiv_sigp
      common/inp_sigp/m_sigp
      common/div_sigp/ndiv_sigp
      real*8 sigmmax,sigmmin,sigmw
      common/ext_sigm/sigmmax,sigmmin,sigmw
      real*8 xsigm(1000),fxsigm(20,1000,20),fxsigmtot(20,1000)
      common/dist_sigm/xsigm,fxsigm,fxsigmtot
      integer m_sigm, ndiv_sigm
      common/inp_sigm/m_sigm
      common/div_sigm/ndiv_sigm
      real*8 polcross(20,-1:1,-1:1),polerror(20,-1:1,-1:1)
      common/polarised/polcross,polerror
      real*8 asycross(20,-1:1,20),asyerror(20,-1:1,20)
      common/ASYX/asycross,asyerror
      real*8 asy(20),asyerr(20)
      common/ASYV/asy,asyerr  
      real*8 yc,yffcut,pzffcut,yffs(8)
      common/ACENTcut/yc,yffcut,pzffcut
      common/AFBSTcut/yffs
      integer nasy
      common/asyparam/nasy