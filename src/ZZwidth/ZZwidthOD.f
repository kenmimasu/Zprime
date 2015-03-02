      SUBROUTINE ZZwidthOD(p1,p2,mass,width,widthmatrix)
c calculates off-diagonal Zp width contributions from light SM fermions 
c in a nearly degenerate 2 gauge boson system at scale psq.
      implicit real*8 (a-h,o-z)
      real*8 mass(10),width(10)
      real*8 p1(0:3),p2(0:3),psq
      common/Z/rmZ,gamZ
      common/t/rmt,gamt
      common/b/rmb
      common/LR/gZZu(2,10),gZZd(2,10),gZZe(2,10),gZZn(2,10),
     &          gZZc(2,10),gZZs(2,10),gZZt(2,10),gZZb(2,10),
     &          gZZta(2,10),gZZnt(2,10)
c Locals
      real*8 gL(2),gR(2),lambda
c Returns
      real*8 widthmatrix(2,2)
ccccccccccccccccccccc
c couplings.
C       print*,"Couplings into width:"
C       print*,"gZZu: (",gZZu,")"
C       print*,"gZZd: (",gZZd,")"
C       print*,"gZZt: (",gZZt,")"
C       print*,"gZZb: (",gZZb,")"
C       print*,"gZZe: (",gZZe,")"
C       print*,"gZZn: (",gZZn,")"
C       print*,"rmZZ: (",mass,")"
C       print*,"gamZZ: (",width,")"
ccccccccccccccccccccc
      psq = (p1(0)+p2(0))**2 - (p1(1)+p2(1))**2 
     &    - (p1(2)+p2(2))**2 - (p1(3)+p2(3))**2 
      pi=dacos(-1.d0) 
c ZZ width.
      ZZwidth=0d0
      temptt=0.d0 
      tempqq=0.d0
      do i=1,6
        tempZZ=0.d0
        do j=1,2
            if((i.eq.2).or.(i.eq.4))then 
c u-quark. 
              rmq=0.d0
              gL(j) = gzzu(1,j)
              gR(j) = gzzu(2,j)
            else if((i.eq.1).or.(i.eq.3))then 
c d-quark. 
              rmq=0.d0
              gL(j) = gzzd(1,j)
              gR(j) = gzzd(2,j)
            else if(i.eq.5)then
c b-quark. 
              rmq=rmb
C               rmq=0d0
              gL(j) = gzzb(1,j)
              gR(j) = gzzb(2,j)
            else if(i.eq.6)then
c t-quark. 
              rmq=rmt
C                 rmq=0d0
              gL(j) = gzzt(1,j)
              gR(j) = gzzt(2,j)
            end if
        enddo
        if(sqrt(psq).le.2.d0*rmq)goto 123
C         if(mass(1).le.2.d0*rmq)goto 123
        
        lambda=(psq**2 - 4d0*psq*rmq**2)/psq**2
        tempZZ=3d0*lambda/24d0/pi*(gL(1)*gL(2)+gR(1)*gR(2))*(psq-rmq**2)
cccccccccccccccccccccc
      if(i.eq.6) then
            temptt=tempZZ
      endif
      if((i.eq.5).or.(i.eq.4).or.(i.eq.3).or.(i.eq.2).or.(i.eq.1)) then
            tempqq=tempqq+tempZZ
      endif
      ZZwidth=ZZwidth+tempZZ
cccccccccccccccccccccc
 123  continue
      end do 
ccccccccccccccccccccc
      temp=0.d0 
      temp1=0.d0
      temp2=0.d0
      temp3=0.d0
cccccccccccccccccccccc
      do i=1,6 
        tempZZ=0.d0
        rmq=0.d0
        do j=1,2
            if((i.eq.2).or.(i.eq.4).or.(i.eq.6))then 
c neutrino.
                gL(j) = gzzn(1,j)
                gR(j) = gzzn(2,j)
            else if(i.eq.6)then ! Tau neutrino
                gL(j) = gzznt(1,j)
                gR(j) = gzznt(2,j)
            else if((i.eq.1).or.(i.eq.3))then 
c lepton. 
                gL(j) = gzze(1,j)
                gR(j) = gzze(2,j)
            else if(i.eq.5)then ! Tau
                gL(j) = gzzta(1,j)
                gR(j) = gzzta(2,j)
            end if 
        enddo
        if(mass(1).le.2.d0*rmq)goto 456
        lambda=(psq-4d0*rmq**2)/psq
        tempZZ=sqrt(lambda)/24d0/pi*(gL(1)*gL(2)+gR(1)*gR(2))*
     &        (psq-rmq**2)
        ZZwidth=ZZwidth+tempZZ
cccccccccccccccccccccc
        temp=temp+tempZZ
        if((i.eq.1).or.(i.eq.3).or.(i.eq.5))then
        temp2=temp2+tempZZ
        else if((i.eq.2).or.(i.eq.4).or.(i.eq.6))then
        temp1=temp1+tempZZ
        endif
        tempZZ=0.d0
      end do
cccccccccccccccccccccc
 456  continue
      do i=1,2
        do j=1,2
          if (i.eq.j) then
              widthmatrix(i,j)=width(i)
          else
              widthmatrix(i,j)=ZZwidth/mass(i) ! Sigma/M
          endif
        enddo
      enddo 
cccccccccccccccccccccc
C       print*,'########### Off-diag widths ###########'
C       print *,'GammaOD(Zp->ff)=',ZZwidth/mass(1),' GeV'
C       print *,'GammaOD(Zp->tt)=',temptt/mass(1),' GeV'
C       print *,'GammaOD(Zp->qq)=',tempqq/mass(1),' GeV'
C       print *,'GammaOD(Zp->ll)=',temp/mass(1),' GeV'
C       print *,'GammaOD(Zp->e/mu/tau)=',temp2/mass(1),' GeV'
C       print *,'GammaOD(Zp->nu)=',temp1/mass(1),' GeV'
cccccccccccccccccccccc
      return
      end
