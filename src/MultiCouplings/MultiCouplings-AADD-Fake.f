c xsection matching AADD-2000 1600<Mtt<2400 : (AADD)0.5569 pb:0.5573 pb(FAKE) ~ dsigma < 1%
c obtained by defining effectiove L/R couplings by adding those of the gamma' and Z' in quad. 
c and then scaling by and overall factor to retain normalisation (account for difference from interference terms)
c SCL=1.15
c uL =   0.29311986
c uR =   0.406641767
c dL =   0.14655993
c dR =   0.406641767
c uL^2+uR^2 =   0.251276779
c dL^2+dR^2 =   0.18683734
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

c cumstom rmZZ assignments for AADD where param = R^-1
c custom VA coupling assignments for quarks in AADD (sqrt(2)*SM) 
c successive pair of indices correspond to ith level KK photon followed by KK Z
      Rinv=param(1)
      rmZZ(1)=Rinv
      ggvu(1)=sqrt(2.d0)*(gau(1)+gau(2))/2.d0
      ggau(1)=sqrt(2.d0)*(gau(1)-gau(2))/2.d0
      ggvu(2)=sqrt(2.d0)*(gzu(1)+gzu(2))/2.d0
      ggau(2)=sqrt(2.d0)*(gzu(1)-gzu(2))/2.d0
      ggvd(1)=sqrt(2.d0)*(gad(1)+gad(2))/2.d0
      ggad(1)=sqrt(2.d0)*(gad(1)-gad(2))/2.d0
      ggvd(2)=sqrt(2.d0)*(gzd(1)+gzd(2))/2.d0
      ggad(2)=sqrt(2.d0)*(gzd(1)-gzd(2))/2.d0


c no heavy neutrinos
      NeuMass=1.d6
      iNeu=0
      SCL=1.18
c conversion to LR
      gZZu(1,1)= dsqrt((ggvu(1)+ggau(1))**2 + (ggvu(2)+ggau(2))**2)/SCL
      gZZu(2,1)= dsqrt((ggvu(1)-ggau(1))**2 + (ggvu(2)-ggau(2))**2)/SCL
      gZZd(1,1)= dsqrt((ggvd(1)+ggad(1))**2 + (ggvd(2)+ggad(2))**2)/SCL
      gZZd(2,1)= dsqrt((ggvd(1)-ggad(1))**2 + (ggvd(2)-ggad(2))**2)/SCL
      gZZe(1,1)= dsqrt((ggve(1)+ggae(1))**2 + (ggve(2)+ggae(2))**2)/SCL
      gZZe(2,1)= dsqrt((ggve(1)-ggae(1))**2 + (ggve(2)-ggae(2))**2)/SCL
      gZZn(1,1)= dsqrt((ggvn(1)+ggan(1))**2 + (ggvn(2)+ggan(2))**2)/SCL
      gZZn(2,1)= dsqrt((ggvn(1)-ggan(1))**2 + (ggvn(2)-ggan(2))**2)/SCL
      gZZc(1,1)= gZZu(1,1)/SCL
      gZZc(2,1)= gZZu(2,1)/SCL
      gZZs(1,1)= gZZd(1,1)/SCL
      gZZs(2,1)= gZZd(2,1)/SCL
      gZZt(1,1)= gZZu(1,1)/SCL
      gZZt(2,1)= gZZu(2,1)/SCL
      gZZb(1,1)= gZZd(1,1)/SCL
      gZZb(2,1)= gZZd(2,1)/SCL
      gZZta(1,1)= gZZe(1,1)/SCL
      gZZta(2,1)= gZZe(2,1)/SCL
      gZZnt(1,1)= gZZn(1,1)/SCL
      gZZnt(2,1)= gZZn(2,1)/SCL
c calculate Z width from couplings assuming all massless except top
c      call ZZwidthMulti(1,rmZZ,rlambdaQCD4,nloop,gamZZ,neumass,ineu)
      gamZZ(1)=95.d0
ccccccccccccccccccccc
ccccccccccccccccccccc
      RETURN
      END
      
      
      
      
