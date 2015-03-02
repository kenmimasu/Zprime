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
      Rinv = param(1)
      do i=1,nzp/2
c            rmZZ(2*i-1)=gp(i)*float(i)*Rinv
            rmZZ(2*i)=gp(i)*sqrt(rmZ*rmZ + (float(i)*Rinv)**2.d0 )
c            ggvu(2*i-1)=gp(i)*sqrt(2.d0)*(gau(1)+gau(2))/2.d0
c            ggau(2*i-1)=gp(i)*sqrt(2.d0)*(gau(1)-gau(2))/2.d0
            ggvu(2*i)=gp(i)*sqrt(2.d0)*(gzu(1)+gzu(2))/2.d0
            ggau(2*i)=gp(i)*sqrt(2.d0)*(gzu(1)-gzu(2))/2.d0
c            ggvd(2*i-1)=gp(i)*sqrt(2.d0)*(gad(1)+gad(2))/2.d0
c            ggad(2*i-1)=gp(i)*sqrt(2.d0)*(gad(1)-gad(2))/2.d0
            ggvd(2*i)=gp(i)*sqrt(2.d0)*(gzd(1)+gzd(2))/2.d0
            ggad(2*i)=gp(i)*sqrt(2.d0)*(gzd(1)-gzd(2))/2.d0
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
      RETURN
      END