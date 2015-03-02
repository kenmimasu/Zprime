      SUBROUTINE ZPCOUP(rmZZ,gp,gZZu,gZZd,gZZc,gZZs,gZZe,gZZn,
     & gZZt,gZZb,gZZta,gZZnt,gamZZ,neumass,ineu)
      implicit none
c takes vector and axial Zp couplings as input and returns left and right arrays
c gZZx(1,i)=Left ; gZZx(1,i)=Right
      real*8 NeuMass,pm
      INTEGER iNeu,nzp2,i,j
      include 'VAcoups.inc'
      real*8 rmZZ(10),gamZZ(10),gp(10)
      real*8 gZZu(2,10),gZZd(2,10),gZZe(2,10),gZZn(2,10),
     &          gZZc(2,10),gZZs(2,10),gZZt(2,10),gZZb(2,10),
     &          gZZta(2,10),gZZnt(2,10)
c conversion to LR
      nzp2=0
      do i=1,10
            if (rmZZ(i).ne.0.d0) then
                nzp2=i
                do j=1,2
                  pm=3d0-2d0*float(j)
                  gZZu(j,i)= (ggvu(i)+pm*ggau(i))*gp(i)/2.d0
                  gZZd(j,i)= (ggvd(i)+pm*ggad(i))*gp(i)/2.d0
                  gZZe(j,i)= (ggve(i)+pm*ggae(i))*gp(i)/2.d0
                  gZZn(j,i)= (ggvn(i)+pm*ggan(i))*gp(i)/2.d0
                  gZZc(j,i)= gZZu(1,i)
                  gZZs(j,i)= gZZd(1,i)
                  gZZt(j,i)= (ggvt(i)+pm*ggat(i))*gp(i)/2.d0
                  gZZb(j,i)= (ggvb(i)+pm*ggab(i))*gp(i)/2.d0
                  gZZta(j,i)= (ggvta(i)+pm*ggata(i))*gp(i)/2.d0
                  gZZnt(j,i)= (ggvnt(i)+pm*ggant(i))*gp(i)/2.d0
                enddo
            else 
                do j=1,2
                  gZZu(j,i)= 0.d0
                  gZZd(j,i)= 0.d0
                  gZZe(j,i)= 0.d0
                  gZZn(j,i)= 0.d0
                  gZZc(j,i)= 0.d0
                  gZZs(j,i)= 0.d0
                  gZZt(j,i)= 0.d0
                  gZZb(j,i)= 0.d0
                  gZZta(j,i)= 0.d0
                  gZZnt(j,i)= 0.d0
                enddo
            end if
            
      end do
c no heavy neutrinos
      NeuMass=1.d6
      iNeu=0      
      call ZZwidthMulti(nzp2,rmZZ,gamZZ,neumass,ineu)
c PDF scales
      RETURN
      END