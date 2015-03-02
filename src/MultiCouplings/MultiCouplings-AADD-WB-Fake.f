      SUBROUTINE ZPCOUP(rmZZ,gp,gZZu,gZZd,gZZc,gZZs,gZZe,gZZn,
     & gZZt,gZZb,gZZta,gZZnt,gamZZ,neumass,ineu)
      implicit none
c takes vector and axial Zp couplings as input and returns left and right arrays
      real*8 param(10)
      integer nzp
      common/pars/param(10),nzp
      real*8 NeuMass
      INTEGER iNeu,i
      real*8 Rinv
c sm coupling declarations
      include 'ewcoups.inc'
      include 'ewparams.inc'
      dimension rmZZ(10),gamZZ(10),gp(10)
      dimension gZZu(2,10),gZZd(2,10),gZZe(2,10),gZZn(2,10),
     &          gZZc(2,10),gZZs(2,10),gZZt(2,10),gZZb(2,10),
     &          gZZta(2,10),gZZnt(2,10)
      real*8    ggvu(10),ggau(10),ggvd(10),ggad(10),
     &          ggvt(10),ggat(10),ggvb(10),ggab(10),
     &          ggve(10),ggae(10),ggvn(10),ggan(10),
     &          ggvta(10),ggata(10),ggvnt(10),ggant(10)
      REAL*8         GT3N(2),GT3L(2),GT3U(2),GT3D(2)
      REAL*8         GYN(2), GYL(2), GYU(2), GYD(2)
C       external rand
C WEAK ISOSPIN QUANTUM NUMBERS
      data GT3N/0.5d0,0.d0/,GT3L/-0.5d0,0.d0/,
     &     GT3U/0.5d0,0.d0/,GT3D/-0.5d0,0.d0/
C HYPERCHARGE QUANTUM NUMBERS
      data GYN/-0.5d0,0.d0/,            GYL/-0.5d0,-1.d0/,
     &     GYU/0.1666667d0,0.6666667d0/,GYD/0.1666667d0,-0.3333333d0/
      parameter (pi=3.14159265358979323846d0)

c cumstom rmZZ assignments for AADD where param = R^-1
c custom VA coupling assignments for quarks in AADD (sqrt(2)*SM) 
C assumption that radiative corrections have suppressed all mass mixing between W3 and B
c successive pair of indices correspond to ith level KK B' followed by KK W3'
c EW params
      em_ch=sqrt(4.d0*pi*a_em)
      sw=sqrt(s2w)
      cw=sqrt(1-s2w)
      g_w=em_ch/sw
      g_y=em_ch/cw
      Rinv = param(1)
      rmZZ(1)=Rinv
      ggvu(1)=g_y*sqrt(2.d0)*(gyu(1)+gyu(2))/2.d0
      ggau(1)=g_y*sqrt(2.d0)*(gyu(1)-gyu(2))/2.d0
      ggvu(2)=g_w*sqrt(2.d0)*(gt3u(1)+gt3u(2))/2.d0
      ggau(2)=g_w*sqrt(2.d0)*(gt3u(1)-gt3u(2))/2.d0
      ggvd(1)=g_y*sqrt(2.d0)*(gyd(1)+gyd(2))/2.d0
      ggad(1)=g_y*sqrt(2.d0)*(gyd(1)-gyd(2))/2.d0
      ggvd(2)=g_w*sqrt(2.d0)*(gt3d(1)+gt3d(2))/2.d0
      ggad(2)=g_w*sqrt(2.d0)*(gt3d(1)-gt3d(2))/2.d0

c no heavy neutrinos
      NeuMass=1.d6
      iNeu=0
      SCL=1.24d0
c conversion to LR
c toy couplings, up only uL^2+uR^2=0.556^2:
      seed=MOD(TIME(),12345)
      print*, seed
      umax=rand(TIME())
      print*, umax
      dmax=0.d0
      gZZu(1,1)= rand(0)*umax
      gZZu(2,1)= sqrt(umax**2-gZZu(1,1)**2)
      gZZd(1,1)= rand(0)*dmax
      gZZd(2,1)= sqrt(dmax**2-gZZd(1,1)**2)
      gZZe(1,1)= 0d0
      gZZe(2,1)= 0d0
      gZZn(1,1)= 0d0
      gZZn(2,1)= 0d0
      gZZc(1,1)= gZZu(1,1)
      gZZc(2,1)= gZZu(2,1)
      gZZs(1,1)= gZZd(1,1)
      gZZs(2,1)= gZZd(2,1)
      gZZt(1,1)= gZZu(1,1)
      gZZt(2,1)= gZZu(2,1)
      gZZb(1,1)= gZZd(1,1)
      gZZb(2,1)= gZZd(2,1)
      gZZta(1,1)= gZZe(1,1)
      gZZta(2,1)= gZZe(2,1)
      gZZnt(1,1)= gZZn(1,1)
      gZZnt(2,1)= gZZn(2,1)
c calculate Z width from couplings assuming all massless except top
      gamZZ(1)=90.d0
ccccccccccccccccccccc
      RETURN
      END
      
      
      
      
