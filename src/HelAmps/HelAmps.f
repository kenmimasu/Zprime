c Calls relevant MADGRAPH helicity amplitude functions for all initial states
c QCD = switch for including QCD ttbar,EW = switch for gamma/Z,
c  GS=strong coupling, IF = fermion number, MASSF = fermion mass
c p1,p2,p3,p4 = MADGRAPH momenta
c CONT = 0: ALL; 1: INTERFERENCE only; 2: NO INTERFERENCE; 3: NONE
      SUBROUTINE ME(GS,MASSF,p1,p2,p3,p4,resgg,
     &      resqq,resqqalt,resuu,resuualt,resdd,resddalt,
     &      rescc,resccalt,resss,resssalt,resbb,resbbalt)
C       implicit real*8 (a-h,o-z)
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
      real*8 sgg_ttb,sqqb_ttb,sbbb_ttb,sqqb_ffb_multi_INT
CC
      do khel=-1,1,2
        do lhel=-1,1,2
cc gg QCD.
          if(QCD.eq.1.AND.jf.ne.3.AND.jf.ne.4) then 
          resgg(khel,lhel)= 
     &    sgg_ttb(JF,MASSF,khel,lhel,p1,p2,p3,p4)*gs**4
c qq-bar QCD
          resqq(khel,lhel)=
     &    sqqb_ttb(JF,MASSF,khel,lhel,p1,p2,p3,p4)*gs**4
          resqqalt(khel,lhel)=
     &    sqqb_ttb(JF,MASSF,khel,lhel,p2,p1,p3,p4)*gs**4
          resbb(khel,lhel)=
     &    sbbb_ttb(JF,MASSF,khel,lhel,p1,p2,p3,p4)*gs**4
          resbbalt(khel,lhel)=
     &    sbbb_ttb(JF,MASSF,khel,lhel,p2,p1,p3,p4)*gs**4
          endif
c EW/Z' CONTRIBUTIONS
          if(jf.ne.4.and.CONT.ne.3)then
c dd-bar 
         resdd(khel,lhel)=sqqb_ffb_multi_INT(khel,lhel,1,JF,
     &                  p1,p2,p3,p4,MASSF,rmzz,gamzz,CONT,EW,BSM,0)
         resddalt(khel,lhel)=sqqb_ffb_multi_INT(khel,lhel,1,JF,
     &                  p2,p1,p3,p4,MASSF,rmzz,gamzz,CONT,EW,BSM,0)
c uu-bar 
         resuu(khel,lhel)=sqqb_ffb_multi_INT(khel,lhel,2,JF,
     &                  p1,p2,p3,p4,MASSF,rmzz,gamzz,CONT,EW,BSM,0)
         resuualt(khel,lhel)=sqqb_ffb_multi_INT(khel,lhel,2,JF,
     &                  p2,p1,p3,p4,MASSF,rmzz,gamzz,CONT,EW,BSM,0)
c ss-bar 
         resss(khel,lhel)=sqqb_ffb_multi_INT(khel,lhel,3,JF,
     &                  p1,p2,p3,p4,MASSF,rmzz,gamzz,CONT,EW,BSM,0)
         resssalt(khel,lhel)=sqqb_ffb_multi_INT(khel,lhel,3,JF,
     &                  p2,p1,p3,p4,MASSF,rmzz,gamzz,CONT,EW,BSM,0)
c cc-bar 
         rescc(khel,lhel)=sqqb_ffb_multi_INT(khel,lhel,4,JF,
     &                  p1,p2,p3,p4,MASSF,rmzz,gamzz,CONT,EW,BSM,0)
         resccalt(khel,lhel)=sqqb_ffb_multi_INT(khel,lhel,4,JF,
     &                  p2,p1,p3,p4,MASSF,rmzz,gamzz,CONT,EW,BSM,0)
c bb-bar
         resbb(khel,lhel)=resbb(khel,lhel) + sqqb_ffb_multi_INT(
     &      khel,lhel,5,JF,p1,p2,p3,p4,MASSF,rmzz,
     &      gamzz,CONT,EW,BSM,1)
         resbbalt(khel,lhel)=resbbalt(khel,lhel) + sqqb_ffb_multi_INT(
     &      khel,lhel,5,JF,p2,p1,p3,p4,MASSF,rmzz,
     &      gamzz,CONT,EW,BSM,1)
          else 
c DIJET CASE !!!NOT IMPLEMENTED FOR NON UNIVERSAL C,S COUPLINGS!!!
c dd-bar -> ddbar,ssbar
         resdd(khel,lhel)=2.d0*sqqb_ffb_multi_INT(khel,lhel,1,1,
     &                  p1,p2,p3,p4,0.d0,rmzz,gamzz,CONT,EW,BSM,0)
         resddalt(khel,lhel)=2.d0*sqqb_ffb_multi_INT(khel,lhel,1,1,
     &                  p2,p1,p3,p4,0.d0,rmzz,gamzz,CONT,EW,BSM,0)
c dd-bar -> uubar,ccbar
         resdd(khel,lhel)=resdd(khel,lhel)+
     &                  2.d0*sqqb_ffb_multi_INT(khel,lhel,1,2,
     &                  p1,p2,p3,p4,0.d0,rmzz,gamzz,CONT,EW,BSM,0)
         resddalt(khel,lhel)=resddalt(khel,lhel)+
     &                  2.d0*sqqb_ffb_multi_INT(khel,lhel,1,2,
     &                  p2,p1,p3,p4,0.d0,rmzz,gamzz,CONT,EW,BSM,0)
c dd-bar -> bbbar
         resdd(khel,lhel)=resdd(khel,lhel)+
     &                  sqqb_ffb_multi_INT(khel,lhel,1,5,
     &                  p1,p2,p3,p4,rmb,rmzz,gamzz,CONT,EW,BSM,0)
         resddalt(khel,lhel)=resddalt(khel,lhel)+
     &                  sqqb_ffb_multi_INT(khel,lhel,1,5,
     &                  p2,p1,p3,p4,rmb,rmzz,gamzz,CONT,EW,BSM,0)
c uu-bar -> ddbar,ssbar
         resuu(khel,lhel)=2.d0*sqqb_ffb_multi_INT(khel,lhel,2,1,
     &                  p1,p2,p3,p4,0.d0,rmzz,gamzz,CONT,EW,BSM,0)
         resuualt(khel,lhel)=2.d0*sqqb_ffb_multi_INT(khel,lhel,2,1,
     &                  p2,p1,p3,p4,0.d0,rmzz,gamzz,CONT,EW,BSM,0)
c uu-bar -> uubar,ccbar
         resuu(khel,lhel)=resuu(khel,lhel)+
     &                  2.d0*sqqb_ffb_multi_INT(khel,lhel,2,2,
     &                  p1,p2,p3,p4,0.d0,rmzz,gamzz,CONT,EW,BSM,0)
         resuualt(khel,lhel)=resuualt(khel,lhel)+
     &                  2.d0*sqqb_ffb_multi_INT(khel,lhel,2,2,
     &                  p2,p1,p3,p4,0.d0,rmzz,gamzz,CONT,EW,BSM,0) 
c uu-bar -> bbbar
         resuu(khel,lhel)=resuu(khel,lhel)+
     &                  sqqb_ffb_multi_INT(khel,lhel,2,5,
     &                  p1,p2,p3,p4,rmb,rmzz,gamzz,CONT,EW,BSM,0)
         resuualt(khel,lhel)=resuualt(khel,lhel)+
     &                  sqqb_ffb_multi_INT(khel,lhel,2,5,
     &                  p2,p1,p3,p4,rmb,rmzz,gamzz,CONT,EW,BSM,0)
c ss-bar and cc-bar
         rescc(khel,lhel)=resuu(khel,lhel)
         resccalt(khel,lhel)=resuualt(khel,lhel)
         resss(khel,lhel)=resdd(khel,lhel)
         resssalt(khel,lhel)=resddalt(khel,lhel)
c bb-bar -> ddbar,ssbar
         resbb(khel,lhel)=2.d0*sqqb_ffb_multi_INT(khel,lhel,5,1,
     &                  p1,p2,p3,p4,0.d0,rmzz,gamzz,CONT,EW,BSM,1)
         resbbalt(khel,lhel)=2.d0*sqqb_ffb_multi_INT(khel,lhel,5,1,
     &                  p2,p1,p3,p4,0.d0,rmzz,gamzz,CONT,EW,BSM,1)
c bb-bar -> uubar,ccbar
         resbb(khel,lhel)=resbb(khel,lhel)+
     &                  2.d0*sqqb_ffb_multi_INT(khel,lhel,5,2,
     &                  p1,p2,p3,p4,0.d0,rmzz,gamzz,CONT,EW,BSM,1)
         resbbalt(khel,lhel)=resbbalt(khel,lhel)+
     &                  2.d0*sqqb_ffb_multi_INT(khel,lhel,5,2,
     &                  p2,p1,p3,p4,0.d0,rmzz,gamzz,CONT,EW,BSM,1) 
c bb-bar -> bbbar
         resbb(khel,lhel)=resbb(khel,lhel)+
     &                  sqqb_ffb_multi_INT(khel,lhel,5,5,
     &                  p1,p2,p3,p4,rmb,rmzz,gamzz,CONT,EW,BSM,1)
         resbbalt(khel,lhel)=resbbalt(khel,lhel)+
     &                  sqqb_ffb_multi_INT(khel,lhel,5,5,
     &                  p2,p1,p3,p4,rmb,rmzz,gamzz,CONT,EW,BSM,1)
          end if
        end do
      end do
      RETURN
      END