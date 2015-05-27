      SUBROUTINE ZZwidthOD(p1,p2,mass,width,widthmatrix)
c reads in  parametrisation for p^2 dependent OD widths entries from file Models/{$MODEL}_sigma12s.txt
      implicit real*8 (a-h,o-z)
      real*8 mass(10),width(10)
      real*8 p1(0:3),p2(0:3)
c Locals
      real*8 sigma12, psq, fqq, fww, ftt, mbar
      include 'names.inc'
      include 'ewparams.inc'
      
c Returns
      real*8 widthmatrix(2,2)

      open(unit=55,file='Models/'//model(1:imodel)//
     &                          '_ODwidths_pdep.txt',status='old')
      read(55,*) fqq,fww,ftt
      close(unit=55)
      
      psq = (p1(0)+p2(0))**2 - (p1(1)+p2(1))**2 
     &    - (p1(2)+p2(2))**2 - (p1(3)+p2(3))**2 
      mbar = sqrt((mass(1**2)+mass(2)**2)/2d0) 
      
      
      sigma12 = psq*fqq
      if (sqrt(psq).gt.(2d0*WMASS)) then
          sigma12= sigma12 + psq**3/WMASS**4*fww
      endif
      if (sqrt(psq).gt.(2d0*rmt)) then
          sigma12= sigma12 + psq*ftt
      endif
      
      print*, sqrt(psq),sigma12/mbar
      do i=1,2
        do j=1,2
          if (i.eq.j) then
              widthmatrix(i,j)=width(i)
          else
              widthmatrix(i,j)=sigma12/mbar
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
