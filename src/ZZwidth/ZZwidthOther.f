      SUBROUTINE ZZwidthOther(nzp,mass,gp,width)
c calculates Zp width contributions from light fermions
      implicit real*8 (a-h,o-z)
      real*8 gp(10),mass(10),width(10)
c Lorenzos extra vars
      include 'ewparams.inc'

      pi=dacos(-1.d0) 
      ctw=sqrt(1.d0-s2w)
      stw=sqrt(s2w) 
      e=sqrt(4.d0*pi*a_em) 
      g=e/sqrt(s2w) 
      GF=1.16639D-5
      v=sqrt(1.d0/(sqrt(2.d0)*GF))
      Sp=0.d0
      Cp=0.d0
c      if (gp(1).eq.0.d0) then
c            gp(1)=1.d-5 
c      endif
      do n = 1,nzp
      rmZZ=mass(n)

      Sp2num=2.d0*gp(2)*sqrt((e/stw)**2+(e/ctw)**2)

c removed g1p dependence
      x=rmZZ/(2.d0)*sqrt(1-gp(2)**2*v**2/(4.d0*rmZZ**2-
     &      v**2*((e/stw)**2+(e/ctw)**2)))

      Cp2num=gp(2)**2+16.d0*(x/v)**2-(e/stw)**2-(e/ctw)**2

      Sp=sin(dasin(Sp2num/sqrt(Sp2num**2+Cp2num**2))/2.d0)  
      Cp=sqrt(1-Sp**2)
cccccccccccccccccccccc

      
c    Other decays: Zp->Zh
      
      ww1=s2w*(1.d0-s2w)*Sp*Cp*rmW*gp(2)**2
      ww2=-(1.d0-2.d0*Sp**2)*stw*ctw*e*rmW*gp(2)
      ww3=-Sp*Cp*e**2*rmW
c      ww4=-stw*(1-s2w)*8*Sa*Sp*Cp*e*x*gp(2)**2
      COUPZpZh = 1.d0/(ctw**2)/e/stw*(ww1+ww2+ww3)
      
      
      pZ=sqrt((rmZZ**2-(rmZ+rmH)**2)
     &       *(rmZZ**2-(rmZ-rmH)**2))
     &  /2.d0/rmZZ     
      EZ=sqrt(pZ**2+rmZ**2)
      wdtZpZh=COUPZpZh**2/(24.d0*dacos(-1.d0))
     &    *pZ*((EZ/rmZ)**2+2.d0)
      wdtZpZh=wdtZpZh/rmZZ**2


      
c    Other decays:  Zp->WW

      COUPZpWW = ctw*e*Sp/sqrt(s2w)
      wdtZpWW=COUPZpWW**2/(48.d0*dacos(-1.d0))*rmZZ
     &    *(1.d0-4.d0*rmW**2/rmZZ**2)**(1.d0/2.d0)
     &    *(0.25d0*rmZZ**4/rmW**4
     &     +4.00d0*rmZZ**2/rmW**2
     &     -17.0d0
     &     -12.d0*rmW**2/rmZZ**2)

c      print *,'Zp-WW coupling', COUPZpWW,ww1,ww2,ww3
   
      Zpwidth=Zpwidth+wdtZpZh+wdtZpWW
c      if (gp(1).eq.0.d0) then
c            Zpwidth=0.d0
c            wdtZpZh=0.d0
c            wdtZpWW=0.d0
c      endif
      width(n)=width(n)+Zpwidth

      print *,'Gamma(Zp(',n,')->Zh)=',wdtZpZh,' GeV'
      print *,'Gamma(Zp(',n,')->WW)=',wdtZpWW,' GeV'
      end do
      return
      end
      
      
      
      