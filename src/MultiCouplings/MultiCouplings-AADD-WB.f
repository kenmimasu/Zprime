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
      Rinv=param(1)
      do i=1,nzp/2
            rmZZ(2*i-1)= float(i)*Rinv
            rmZZ(2*i)= float(i)*Rinv
            ggvu(2*i-1)=g_y*sqrt(2.d0)*(gyu(1)+gyu(2))/2.d0
            ggau(2*i-1)=g_y*sqrt(2.d0)*(gyu(1)-gyu(2))/2.d0
            ggvu(2*i)  =g_w*sqrt(2.d0)*(gt3u(1)+gt3u(2))/2.d0
            ggau(2*i)  =g_w*sqrt(2.d0)*(gt3u(1)-gt3u(2))/2.d0
            ggvd(2*i-1)=g_y*sqrt(2.d0)*(gyd(1)+gyd(2))/2.d0
            ggad(2*i-1)=g_y*sqrt(2.d0)*(gyd(1)-gyd(2))/2.d0
            ggvd(2*i)  =g_w*sqrt(2.d0)*(gt3d(1)+gt3d(2))/2.d0
            ggad(2*i)  =g_w*sqrt(2.d0)*(gt3d(1)-gt3d(2))/2.d0
      end do
c no heavy neutrinos
      NeuMass=1.d6
      iNeu=0
c conversion to LR
      do i=1,nzp
            gZZu(1,i)= (ggvu(i)+ggau(i))*gp(i)
            gZZu(2,i)= (ggvu(i)-ggau(i))*gp(i)
            gZZd(1,i)= (ggvd(i)+ggad(i))*gp(i)
            gZZd(2,i)= (ggvd(i)-ggad(i))*gp(i)
            gZZe(1,i)= (ggve(i)+ggae(i))*gp(i)
            gZZe(2,i)= (ggve(i)-ggae(i))*gp(i)
            gZZn(1,i)= (ggvn(i)+ggan(i))*gp(i)
            gZZn(2,i)= (ggvn(i)-ggan(i))*gp(i)
            gZZc(1,i)= gZZu(1,i)
            gZZc(2,i)= gZZu(2,i)
            gZZs(1,i)= gZZd(1,i)
            gZZs(2,i)= gZZd(2,i)
            gZZt(1,i)= gZZu(1,i)
            gZZt(2,i)= gZZu(2,i)
            gZZb(1,i)= gZZd(1,i)
            gZZb(2,i)= gZZd(2,i)
            gZZta(1,i)= gZZe(1,i)
            gZZta(2,i)= gZZe(2,i)
            gZZnt(1,i)= gZZn(1,i)
            gZZnt(2,i)= gZZn(2,i)
      end do
c calculate Z width from couplings assuming all massless except top
      call ZZwidthMulti(nzp,rmZZ,gamZZ,neumass,ineu)
      
ccccccccccccccccccccc
      RETURN
      END
      
      
      
      
