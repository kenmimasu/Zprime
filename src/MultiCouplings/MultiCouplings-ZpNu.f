      SUBROUTINE ZPCOUP(rmZZ,gp,gZZu,gZZd,gZZc,gZZs,gZZe,gZZn,
     & gZZt,gZZb,gZZta,gZZnt,gamZZ,neumass,ineu)
      implicit none
c takes vector and axial Zp couplings as input and returns left and right arrays
      real*8 NeuMass,QQ
      INTEGER iNeu
      real*8 rmZZ(10),gamZZ(10),gp(10)
      real*8 gZZu(2,10),gZZd(2,10),gZZe(2,10),gZZn(2,10),
     &          gZZc(2,10),gZZs(2,10),gZZt(2,10),gZZb(2,10),
     &          gZZta(2,10),gZZnt(2,10)
      include 'ewparams.inc'
      include 'VAcoups.inc'

c ZpNu model: gp(1)=g1p, gp(2)=gtilde ggvn(1)=neumass ggvn(2)=iNeu switch
c conversion to LR
      gZZu(1,1)= (ggvu(1)+ggau(1))*gp(1)+1.d0/6.d0*gp(2)
      gZZu(2,1)= (ggvu(1)-ggau(1))*gp(1)+2.d0/3.d0*gp(2)
      gZZd(1,1)= (ggvd(1)+ggad(1))*gp(1)+1.d0/6.d0*gp(2)
      gZZd(2,1)= (ggvd(1)-ggad(1))*gp(1)-1.d0/3.d0*gp(2)
      gZZe(1,1)= (ggve(1)+ggae(1))*gp(1)-1.d0/2.d0*gp(2)
      gZZe(2,1)= (ggve(1)-ggae(1))*gp(1)-1.d0*gp(2)
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
c Majorana heavy neutrinos can only have axial couplings
      gZZn(1,1)= (ggvn(1)+ggan(1))*gp(1)-1.d0/2.d0*gp(2)
      gZZn(2,1)= (ggvn(1)+ggan(1))*gp(1)
      gZZnt(1,1)= gZZn(1,1)
      gZZnt(2,1)= gZZn(2,1)
      NeuMass=ggan(2)
      iNeu=int(ggan(3))
c      print*,'gamZZ =',gamZZ,'    (GeV)'
      call ZZwidthMulti(1,rmZZ,gamZZ,neumass,ineu)
      call ZZwidthOther(1,rmZZ,gp,gamZZ)
c extra info
c      print*,'rmZZ =',rmZZ,'    (GeV)'
c      print*,'gamZZ =',gamZZ,'    (GeV)'
      print*,'##### ZpNu params: #####'
      print*,'g1p:',gp(1)
      print*,'gtilde:',gp(2)

      
      RETURN
      END