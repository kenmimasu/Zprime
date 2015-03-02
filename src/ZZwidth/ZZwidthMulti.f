      SUBROUTINE ZZwidthMulti(nzprime,mass,width,neumass,ineu)
c calculates Zp width contributions from light fermions
      implicit real*8 (a-h,o-z)
      include 'ewparams.inc'
      include 'zpparams.inc'
      include 'LRcoups.inc'
      COMMON/ALFASTRONG/rlambdaQCD4,nloops
c arguments
      real*8 mass(10),width(10)
      real*8 NeuMass
      integer nzprime,ineu
ccccccccccccccccccccc
c initialise
      do n=1,10
            width(n)=0.d0
      enddo
cccccccccccccccccccccc
c couplings.
      pi=dacos(-1.d0) 
      do n = 1,nzprime
c ZZ width.
      ZZwidth=0.d0
      temptt=0.d0 
      tempqq=0.d0
      if (mass(n).gt.0.d0) then 
            a_s=alfas(mass(n),rlambdaQCD4,nloop)
      else 
            goto 124
      endif
      do i=1,6
        tempZZ=0.d0
        if((i.eq.2).or.(i.eq.4))then 
c u-quark. 
          rmq=0.d0
          gv = (gzzu(1,n)+gzzu(2,n))
          ga = (gzzu(1,n)-gzzu(2,n))
        else if((i.eq.1).or.(i.eq.3))then 
c d-quark. 
          rmq=0.d0
          gv = (gzzd(1,n)+gzzd(2,n))
          ga = (gzzd(1,n)-gzzd(2,n))
        else if(i.eq.5)then
          rmq=rmb
          gv = (gzzb(1,n)+gzzb(2,n))
          ga = (gzzb(1,n)-gzzb(2,n))
        else if(i.eq.6)then
          rmq=rmt
          gv = (gzzt(1,n)+gzzt(2,n))
          ga = (gzzt(1,n)-gzzt(2,n))
        end if
        if(mass(n).le.2.d0*rmq)goto 123
c with QCD kfactor
        tempZZ=3.d0/48.d0/pi*mass(n)
     &               *sqrt(1.d0-4.d0*rmq**2/mass(n)**2)
     &               *(gv**2*(1.d0+2.d0*rmq**2/mass(n)**2)
     &                +ga**2*(1.d0-4.d0*rmq**2/mass(n)**2))
     &               *(1.d0+1.045d0*a_s/pi)
            ZZwidth=ZZwidth+tempZZ
c without QCD kfactor
C        tempZZ = 3.d0/12.d0/pi*mass(n)
C     &               *sqrt(1.d0-4.d0*rmq**2/mass(n)**2)
C     &               *(gv**2*(1.d0+2.d0*rmq**2/mass(n)**2)
C     &                +ga**2*(1.d0-4.d0*rmq**2/mass(n)**2))
C        ZZwidth=ZZwidth+tempZZ
C             print*,i,tempZZ
cccccccccccccccccccccc
      if(i.eq.6) then
            temptt=tempZZ
      endif
      if((i.eq.5).or.(i.eq.4).or.(i.eq.3).or.(i.eq.2).or.(i.eq.1)) then
            tempqq=tempqq+tempZZ
      endif
      
cccccccccccccccccccccc
 123  continue
      end do 

c      print *, 'ineu = ', ineu
C       print*,'Z'' width due to quarks:',ZZwidth,' GeV'
ccccccccccccccccccccc
      temp=0.d0 
      temp1=0.d0
      temp2=0.d0
      temp3=0.d0
cccccccccccccccccccccc
      do i=1,9 
        tempZZ=0.d0
        rmq=0.d0
        if((i.eq.2).or.(i.eq.4))then 
c neutrino.
          if(iNeu.eq.0) then
          gV = (gzzn(1,n)+gzzn(2,n))
          gA = (gzzn(1,n)-gzzn(2,n))
          dof=1.d0
          else if(iNeu.eq.1) then
          gV=0.d0
          gA = gzzn(1,n)
          dof=0.5d0
          endif
        else if(i.eq.6)then ! Tau neutrino
          if(iNeu.eq.0) then
          gV = (gzznt(1,n)+gzznt(2,n))
          gA = (gzznt(1,n)-gzznt(2,n))
          dof=1.d0
          else if(iNeu.eq.1) then
          gV=0.d0
          gA = gzznt(1,n)
          dof=0.5d0
          endif
        else if((i.eq.1).or.(i.eq.3))then 
c lepton. 
          gv = (gzze(1,n)+gzze(2,n))
          ga = (gzze(1,n)-gzze(2,n))
          dof=1.d0
        else if(i.eq.5)then ! Tau
          gv = (gzzta(1,n)+gzzta(2,n))
          ga = (gzzta(1,n)-gzzta(2,n))
          dof=1.d0
        else if((i.eq.7).or.(i.eq.8).or.(i.eq.9))then
c heavy neutrinos 
          rmq=NeuMass
          if(iNeu.eq.0) then
          gV = (gzzn(1,n)+gzzn(2,n))
          gA = (gzzn(1,n)-gzzn(2,n))
          dof=1.d0
          else if(iNeu.eq.1) then
          gV=0.d0
          gA = gzzn(1,n)
          dof=0.5d0
          endif
        end if 
        if(mass(n).le.2.d0*rmq)goto 456
        tempZZ = dof/48.d0/pi*mass(n)
     &               *sqrt(1.d0-4.d0*rmq**2/mass(n)**2)
     &               *(gv**2*(1.d0+2.d0*rmq**2/mass(n)**2)
     &                +ga**2*(1.d0-4.d0*rmq**2/mass(n)**2))
        ZZwidth=ZZwidth+tempZZ
        
cccccccccccccccccccccc
        temp=temp+tempZZ
        if((i.eq.1).or.(i.eq.3).or.(i.eq.5))then
        temp2=temp2+tempZZ
        else if((i.eq.2).or.(i.eq.4).or.(i.eq.6))then
        temp1=temp1+tempZZ
        else if((i.eq.7).or.(i.eq.8).or.(i.eq.9))then
        temp3=temp3+tempZZ
        endif
        tempZZ=0.d0
cccccccccccccccccccccc
      end do
 456  continue
      width(n)=ZZwidth
 124  continue      
cccccccccccccccccccccc
      print *,'########### Width and BRs ###########'
      print *,'Gamma(Zp(',n,')->ff)=',width(n),' GeV'
      print *,'Gamma(Zp(',n,')->tt)=',temptt,' GeV'
      print *,'Gamma(Zp(',n,')->qq)=',tempqq,' GeV'
      print *,'Gamma(Zp(',n,')->ll)=',temp,' GeV'
      print *,'Gamma(Zp(',n,')->e/mu/tau)=',temp2,' GeV'
      print *,'Gamma(Zp(',n,')->nu)=',temp1,' GeV'
      print *,'Gamma(Zp(',n,')->nuH)=',temp3,' GeV'
cccccccccccccccccccccc
      end do
      return
      end
