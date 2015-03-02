      SUBROUTINE ZPCOUP(rmZZ,gp,gZZu,gZZd,gZZc,gZZs,gZZe,gZZn,
     & gZZt,gZZb,gZZta,gZZnt,gamZZ,neumass,ineu)
      implicit none
c takes vector and axial Zp couplings as input and returns left and right arrays
      include 'VAcoups.inc'
      real*8 NeuMass
      INTEGER iNeu,nzp,i
      real*8 rmZZ(10),gamZZ(10),gp(10)
      real*8 gZZu(2,10),gZZd(2,10),gZZe(2,10),gZZn(2,10),
     &          gZZc(2,10),gZZs(2,10),gZZt(2,10),gZZb(2,10),
     &          gZZta(2,10),gZZnt(2,10)
c conversion to LR
      nzp=5
      do i=1,5
            if (rmZZ(i).ne.0.d0) then
                  nzp=i
c                  gp(i)=zc_norm
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
                  gZZta(1,i)= (ggve(i)+ggae(i))*gp(i)
                  gZZta(2,i)= (ggve(i)-ggae(i))*gp(i)
                  gZZnt(1,i)= (ggvn(i)+ggan(i))*gp(i)
                  gZZnt(2,i)= (ggvn(i)-ggan(i))*gp(i)
                  gZZt(1,i)= (ggvu(i+5)+ggau(i+5))*gp(i)
                  gZZt(2,i)= (ggvu(i+5)-ggau(i+5))*gp(i)
                  gZZb(1,i)= (ggvd(i+5)+ggad(i+5))*gp(i)
                  gZZb(2,i)= (ggvd(i+5)-ggad(i+5))*gp(i)
            else 
                  gZZu(1,i)= 0.d0
                  gZZu(2,i)= 0.d0
                  gZZd(1,i)= 0.d0
                  gZZd(2,i)= 0.d0
                  gZZe(1,i)= 0.d0
                  gZZe(2,i)= 0.d0
                  gZZn(1,i)= 0.d0
                  gZZn(2,i)= 0.d0
                  gZZc(1,i)= 0.d0
                  gZZc(2,i)= 0.d0
                  gZZs(1,i)= 0.d0
                  gZZs(2,i)= 0.d0
                  gZZta(1,i)= 0.d0
                  gZZta(2,i)= 0.d0
                  gZZnt(1,i)= 0.d0
                  gZZnt(2,i)= 0.d0
                  gZZt(1,i)= 0.d0
                  gZZt(2,i)= 0.d0
                  gZZb(1,i)= 0.d0
                  gZZb(2,i)= 0.d0
            end if
            
      end do
c no heavy neutrinos
      NeuMass=1.d6
      iNeu=0
      call ZZwidthMulti(nzp,rmZZ,gamZZ,neumass,ineu)
      RETURN
      END