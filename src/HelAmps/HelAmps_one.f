c Calls relevant MADGRAPH helicity amplitude functions for all initial states
c QCD = switch for including QCD ttbar,EW = switch for gamma/Z,
c  GS=strong coupling, IF = fermion number, MASSF = fermion mass
c p1,p2,p3,p4 = MADGRAPH momenta
c CONT = 0: ALL; 1: INTERFERENCE only; 2: NO INTERFERENCE; 3: NONE
      SUBROUTINE ME(GS,MASSF,p1,p2,p3,p4,resgg,
     &      resqq,resqqalt,resuu,resuualt,resdd,resddalt,
     &      rescc,resccalt,resss,resssalt,resbb,resbbalt)
      implicit none
      include 'runparams.inc'
      include 'zpparams.inc'
      include 'ewparams.inc'
C Arguments
      real*8 MASSF,gs
      real*8 p1(0:3),p2(0:3),p3(0:3),p4(0:3)
      real*8 resgg(-1:1,-1:1),resqq(-1:1,-1:1),resqqalt(-1:1,-1:1)
      real*8 resdd(-1:1,-1:1),resuu(-1:1,-1:1),resss(-1:1,-1:1),
     &          rescc(-1:1,-1:1),resbb(-1:1,-1:1)
      real*8 resddalt(-1:1,-1:1),resuualt(-1:1,-1:1),
     &       resssalt(-1:1,-1:1),resccalt(-1:1,-1:1),
     &       resbbalt(-1:1,-1:1)
C Locals
      integer khel,lhel
CC
      do khel=-1,1,2
        do lhel=-1,1,2
cc gg QCD.
          if(QCD.eq.1.AND.jf.ne.3.AND.jf.ne.4) then 
          resgg(khel,lhel)= 0.25d0
c qq-bar QCD
          resqq(khel,lhel)= 0.25d0
          resqqalt(khel,lhel)= 0.25d0
          resbb(khel,lhel)= 0.25d0
          resbbalt(khel,lhel)= 0.25d0
          endif
c EW/Z' CONTRIBUTIONS
          if(jf.ne.4.and.CONT.ne.3)then
c dd-bar 
         resdd(khel,lhel)= 0.25d0
         resddalt(khel,lhel)= 0.25d0
c uu-bar 
         resuu(khel,lhel)= 0.25d0
         resuualt(khel,lhel)= 0.25d0
c ss-bar 
         resss(khel,lhel)= 0.25d0
         resssalt(khel,lhel)= 0.25d0
c cc-bar 
         rescc(khel,lhel)= 0.25d0
         resccalt(khel,lhel)= 0.25d0
c bb-bar
         resbb(khel,lhel)= 0.25d0
         resbbalt(khel,lhel)= 0.25d0
          else 
c DIJET CASE !!!NOT IMPLEMENTED FOR NON UNIVERSAL C,S COUPLINGS!!!
c dd-bar -> ddbar,ssbar
         resdd(khel,lhel)= 0.25d0
         resddalt(khel,lhel)= 0.25d0
c dd-bar -> uubar,ccbar
         resdd(khel,lhel)= 0.25d0
         resddalt(khel,lhel)= 0.25d0
c dd-bar -> bbbar
         resdd(khel,lhel)= 0.25d0
         resddalt(khel,lhel)= 0.25d0
c uu-bar -> ddbar,ssbar
         resuu(khel,lhel)= 0.25d0
         resuualt(khel,lhel)= 0.25d0
c uu-bar -> uubar,ccbar
         resuu(khel,lhel)= 0.25d0
         resuualt(khel,lhel)= 0.25d0
c uu-bar -> bbbar
         resuu(khel,lhel)= 0.25d0
         resuualt(khel,lhel)= 0.25d0
c ss-bar and cc-bar
         rescc(khel,lhel)=resuu(khel,lhel)
         resccalt(khel,lhel)=resuualt(khel,lhel)
         resss(khel,lhel)=resdd(khel,lhel)
         resssalt(khel,lhel)=resddalt(khel,lhel)
c bb-bar -> ddbar,ssbar
         resbb(khel,lhel)= 0.25d0
         resbbalt(khel,lhel)= 0.25d0
c bb-bar -> uubar,ccbar
         resbb(khel,lhel)= 0.25d0
         resbbalt(khel,lhel)= 0.25d0 
c bb-bar -> bbbar
         resbb(khel,lhel)= 0.25d0
         resbbalt(khel,lhel)= 0.25d0
          end if
        end do
      end do
      RETURN
      END