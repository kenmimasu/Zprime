      SUBROUTINE ZZwidthOD_GB(p1,p2,GZZW,mass,width,widthmatrix)
c calculates off-diagonal Zp width contributions from W bosons
c in a 2 gauge boson system at scale psq.
      implicit real*8 (a-h,o-z)
      real*8 mass(10),width(10)
      real*8 p1(0:3),p2(0:3),psq
      real*8 GZZW(2)
      include 'ewparams.inc'
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
        if(psq.le.2.d0*rmq)goto 456
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
          tempZZ = psq**3*/192d0/pi/GZZW(i)*GZZW(j)/rmW**2    
            
            
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
