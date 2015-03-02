      SUBROUTINE ZPCOUP(rmZZ,gp,gZZu,gZZd,gZZc,gZZs,gZZe,gZZn,
     & gZZt,gZZb,gZZta,gZZnt,gamZZ,neumass,ineu)
      implicit real*8 (a-h,o-z)
c takes vector and axial Zp(i) couplings as input and returns left and right arrays
c gZZx(1,i)=Left ; gZZx(2,i)=Right
      common/pars/param(10),nzp
      include 'VAcoups.inc'
      real*8 NeuMass,pm
      INTEGER iNeu
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
                pm=3d0-2d0*float(j) ! +1 if j=1, -1 if j=2
                gZZu(j,i)= (ggvu(i)+pm*ggau(i))*gp(i)/2.d0
                gZZd(j,i)= (ggvd(i)+pm*ggad(i))*gp(i)/2.d0
                gZZe(j,i)= (ggve(i)+pm*ggae(i))*gp(i)/2.d0
                gZZn(j,i)= (ggvn(i)+pm*ggan(i))*gp(i)/2.d0
                gZZc(j,i)= gZZu(j,i)
                gZZs(j,i)= gZZd(j,i)
                gZZt(j,i)= gZZu(j,i)
                gZZb(j,i)= gZZd(j,i)
                gZZta(j,i)= gZZe(j,i)
                gZZnt(j,i)= gZZn(j,i)
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
C       QQ=0.d0
      RETURN
      END