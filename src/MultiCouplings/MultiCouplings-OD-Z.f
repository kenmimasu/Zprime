      SUBROUTINE ZPCOUP(rmZZ,gp,gZZu,gZZd,gZZc,gZZs,gZZe,gZZn,
     & gZZt,gZZb,gZZta,gZZnt,gamZZ,neumass,ineu)
      implicit none
c takes vector and axial Zp couplings as input and returns left and right arrays
c gZZx(1,i)=Left ; gZZx(1,i)=Right
      real*8 NeuMass,pm
      INTEGER iNeu,nzp2,i,j
      real*8 rmZZ(10),gamZZ(10),gp(10)
      real*8 gZZu(2,10),gZZd(2,10),gZZe(2,10),gZZn(2,10),
     &       gZZc(2,10),gZZs(2,10),gZZt(2,10),gZZb(2,10),
     &       gZZta(2,10),gZZnt(2,10)
      include 'VAcoups.inc'
      include 'ewcoups.inc'
      include 'ewparams.inc'
c conversion to LR

      nzp2=1
      do j=1,2 ! set first Z' couplings to Z couplings
        gZZu(j,1) = GZU(j)
        gZZd(j,1) = GZD(j)
        gZZe(j,1) = GZL(j)
        gZZn(j,1) = GZN(j)
        gZZc(j,1) = GZU(j)
        gZZs(j,1) = GZD(j)
        gZZt(j,1) = GZU(j)
        gZZb(j,1) = GZD(j)
        gZZta(j,1)= GZL(j)
        gZZnt(j,1)= GZN(j)
      enddo
      ! set second Z' couplings to Z' couplings
      do j=1,2
        pm=3d0-2d0*float(j)
        gZZu(j,2) = (ggvu(1)+pm*ggau(1))*gp(1)/2.d0
        gZZd(j,2) = (ggvd(1)+pm*ggad(1))*gp(1)/2.d0
        gZZe(j,2) = (ggve(1)+pm*ggae(1))*gp(1)/2.d0
        gZZn(j,2) = (ggvn(1)+pm*ggan(1))*gp(1)/2.d0
        gZZc(j,2) = gZZu(1,2)
        gZZs(j,2) = gZZd(1,2)
        gZZt(j,2) = (ggvt(1)+pm*ggat(1))*gp(1)/2.d0
        gZZb(j,2) = (ggvb(1)+pm*ggab(1))*gp(1)/2.d0
        gZZta(j,2)= (ggvta(1)+pm*ggata(1))*gp(1)/2.d0
        gZZnt(j,2)= (ggvnt(1)+pm*ggant(1))*gp(1)/2.d0
      enddo
      
c no heavy neutrinos
      NeuMass=1.d6
      iNeu=0      
      call ZZwidthMulti(nzp2,rmZZ,gamZZ,neumass,ineu)
C hack Z' masses & widths to include SM Z
      rmZZ(2)=rmZZ(1)
      rmZZ(1)=ZMASS
      gamZZ(2)=gamZZ(1)
      gamZZ(1)=ZWIDTH

      RETURN
      END