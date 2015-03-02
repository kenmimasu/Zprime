      SUBROUTINE ZZwidthOD(p1,p2,mass,width,widthmatrix)
c reads in (diagonal) off diagonal width entries from file Models/{$MODEL}_ODwidths.txt
      implicit real*8 (a-h,o-z)
      real*8 mass(10),width(10)
      real*8 p1(0:3),p2(0:3)
c Locals
      real*8 sigma12
      include 'names.inc'
      
c Returns
      real*8 widthmatrix(2,2)

      open(unit=55,file='Models/'//model(1:imodel)//
     &                                    '_ODwidths.txt',status='old')
      read(55,*) sigma12
      close(unit=55)
      do i=1,2
        do j=1,2
          if (i.eq.j) then
              widthmatrix(i,j)=width(i)
          else
              widthmatrix(i,j)=sigma12
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
