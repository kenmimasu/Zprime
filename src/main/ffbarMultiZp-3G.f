c      program d_pptt_LO
      implicit real*8 (a-h,o-z)
      include 'vegas.inc'
      include 'ewparams.inc'
      include 'zpparams.inc'
      include 'runparams.inc'
      include 'VAcoups.inc'
      include 'LRcoups.inc'
c Distributions.
      include 'dists_common.inc'
cc Asymmetries
      include 'asy_common.inc'
c few extra variables for calculating CL, invariant mass cutt
      common/lum/alumpb
      common/cut/ainvmasscut
c gauge boson masses and widths.
      data gamW/2.08d0/
      data rmb/4.5d0/
      data rmt/175.d0/,gamt/1.55d0/
      data rmZ/91.19d0/,gamZ/2.50d0/
      data rmH/150.d0/,gamH/0.1635D-01/
      data ierr/0/
      data nzp/10/
      data iNeu/0/
      data NeuMass/0.d0/
      include 'plotparams.inc'
      dimension cnorm(20)
      dimension snorm(20),ave(20)
      dimension poltot(-1:1,-1:1),polchi(-1:1,-1:1)
      dimension asytot(-1:1,20),asychi(-1:1,20)
      parameter (pi=3.14159265358979323846d0)
      parameter (conv=.38937966d9)
C       external fxn,fxncosx
      external fxn
c     filenames
      character*4 name
      character*1 intstr
      character*6 mff
c some EW parameters.
      a_em=1.d0/128.d0
      s2w=.2320d0
      rmW=rmZ*sqrt(1.d0-s2w)
      sw=sqrt(s2w)
      c2w=1.d0-s2w
      cw=sqrt(c2w)
      em_ch=sqrt(4.d0*pi*a_em)
      zc_norm=em_ch/sw/cw/2.d0
C     Read in .com file
      read(*,*) ecm_coll,icoll
      read(*,*) ISTRUCTURE
      read(*,*) y
      read(*,*) ptCut
      read(*,*) iseed,jdummy,kdummy
      read(*,*) ncall
      read(*,*) itmx
      read(*,*) acc
      read(*,*) idummy
      read(*,*) alumpb
      read(*,*) ainvmasscut
      read(*,*) ainvmassup
      read(*,*) ndist
      read(*,*) isumcosx
      read(*,*) cont
      read(*,*) QCD
      read(*,*) EW
      read(*,*) BSM
      read(*,*) mf
      read(*,*) name
      read(*,*) model
      read(*,*) QQ
      imodel = len(model)
      do while(model(imodel:imodel).eq.'') 
            imodel = imodel-1
      end do
c itmx no more than 20 + 1 or 1 for ncutmtt,ndist
      if(itmx.gt.20)then
        write(*,*)'itmx does not have to exceed 20 !'
        stop
      end if
      if(ndist.ne.1.AND.ndist.ne.0)then
        write(*,*)'Invalid value for ndist! (1 for yes and 0 for no)'
        stop
      end if
c files containing 95% CL total cross section as a function of M_s 
c distributions setup.
c UPPER BOUNDS FOR Mtt, Et, Pt: either total available energy OR up to the cutoff (0.8*M_S) (ncutmtt = 0 OR 1)
      if(ainvmassup.eq.0.d0)then  
            rmassmax=ecm_coll
      else 
            rmassmax=ainvmassup
      endif
      rmassmin=ainvmasscut
      if (pTmax.eq.-999d0) pTmax=rmassmax/2.d0     
      if (Etmax.eq.-999d0) Etmax=rmassmax/2.d0
      if (Pzffmax.eq.-999d0) then
      Pzffmax=ecm_coll/2.d0*(1.d0 - rmassmin*rmassmin/ecm_coll/ecm_coll)
      endif
      if (ndist.ne.-1) then
        m_cost3 = ndist
        m_cost4 = ndist
        m_cost3col = ndist
        m_cost4col = ndist
        m_Y3 = ndist
        m_Y4 = ndist
        m_Y3col = ndist
        m_Y4col = ndist
        m_eta3 = ndist
        m_eta4 = ndist
        m_eta3col = ndist
        m_eta4col = ndist
        m_Yff = ndist
        m_Pzff = ndist
        m_rmass = ndist
        m_beta = ndist
        m_delY = ndist
        m_Et = ndist
        m_pT = ndist
      endif
      m_sigp=m_rmass
      sigpmax=rmassmax
      sigpmin=rmassmin
      ndiv_sigp=ndiv_rmass
      m_sigm=m_rmass
      sigmmax=rmassmax
      sigmmin=rmassmin
      ndiv_sigm=ndiv_rmass
c model filename
      imodel = len(model)
      do while(model(imodel:imodel).eq.'') 
            imodel = imodel-1
      end do
      open(unit=42,file='Models/'//model(1:imodel)//'.mdl',status='old')
c read in model file and set couplings
      read(42,*) rmZZ
      read(42,*) gp
      read(42,*) param  
      read(42,*) ggvu
      read(42,*) ggau
      read(42,*) ggvd
      read(42,*) ggad
      read(42,*) ggvt
      read(42,*) ggat
      read(42,*) ggvb
      read(42,*) ggab
      read(42,*) ggve
      read(42,*) ggae
      read(42,*) ggvn
      read(42,*) ggan
      read(42,*) ggvta
      read(42,*) ggata
      read(42,*) ggvnt
      read(42,*) ggant
c iuniv selects universal (0) or non-universal (1) couplings
      iuniv=0
      call initialize(rmt,gamt)
      call SMCORR(rmt,gamt,rmb,rmW,rmH,
     &            DgZu,DgZd,DgZe,DgZn,DgZt,DgZb,DgZta,DgZnt,iuniv)
c initialize MadGraph for QCD MEs.     
ccccccccccccccccccccccccccc
c SET FINAL STATE
      dycut=1.d3
c ttbar
      if(mf.eq.0)then
            rmf=rmt
            if (iuniv.eq.0) then
                  jf=2
            else if (iuniv.eq.1) then
                  jf=6
            endif
            f="t"
            ff="tt"
            ptCut=0.d0
c bbbar
      else if(mf.eq.1)then
            rmf=rmb
            if (iuniv.eq.0) then
                  jf=1
            else if (iuniv.eq.1) then
                  jf=5
            endif            
            f="b"
            ff="bb"
            
            if(ptCut.eq.0.d0) ptCut=300.d0
c epem
      else if(mf.eq.2)then
            rmf=0.d0
            jf=3
            f="e"
            ff="ee"
            if(ptCut.eq.0.d0)ptCut=25.d0
c dijet ONLY USE FOR SIGNAL! QCD BKG MISSING LOADS OF BITS FOR jj
      else if(mf.eq.3)then
            rmf=0.d0
            jf=4
            f="j"
            ff="jj"
C ETA CUT FOR DIJET RESONANCE SEARCH MATCHING
            y=2.5
            dycut=1.3d0
            ptCut=30.d0
      endif
ccccccccccccccccccccccccccc
c unphysical invariant mass cut
      if(ainvmasscut.lt.2.d0*rmf.AND.ainvmasscut.ne.0.d0)then
            write(*,*)'unphysical invariant mass cut ( < 2*Mt ) !'
        stop
      end if
C QCDL4 IS QCD LAMBDA4 (to match PDF fits).
      IF(ISTRUCTURE.EQ.1)THEN
      QCDL4=0.326D0
      nloops=2
      ELSE IF(ISTRUCTURE.EQ.2)THEN
      QCDL4=0.326D0
      nloops=2
      ELSE IF(ISTRUCTURE.EQ.3)THEN
      QCDL4=0.326D0
      nloops=2
      ELSE IF(ISTRUCTURE.EQ.4)THEN
      QCDL4=0.215D0
      nloops=1
      ENDIF
      rlambdaQCD4=QCDL4
c initialise CTEQ grids.
      IF(ISTRUCTURE.LE.4)THEN
        ICTEQ=ISTRUCTURE
        call SetCtq6(ICTEQ)
      END IF 
cccccccccccccccccccccccccccccccc
      CALL ZPCOUP(rmZZ,gp,
     & gZZu,gZZd,gZZc,gZZs,gZZe,gZZn,gZZt,gZZb,gZZta,gZZnt,gamZZ,
     & neumass,ineu,QQ) 
c collider CM energy squared.      
      s=ecm_coll*ecm_coll
c output label (undirected).
      jh=6
c fac=constant factor outside integration: conversion GeV^2 -> pb.
      fac=conv       
c azimuthal angle integrated out (no initial transverse polarisation).
      fac=fac*2.d0*pi
c string of cuts for AC,ARFB,AOFB,AFBSTAR
      yc = 0.5d0
      yffcut = 0.5d0
      do i=1,8
         yffs(i)=0.2d0*i
      end do
      pzffcut = 700.d0
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c switch of whether or not to manually sum over cos and x
      if(ISUMCOSX.eq.1)then
c VEGAS parameters.
      ndim=2
c if nprn<0 no print-out.
      nprn=0
c integrates on:
c x(1)=(x1-tau)/(1-tau),
c x(2)=(ecm-rm3-rm4)/(ecm_max-rm3-rm4),
c x(3)=ct NOT integrated over in vegas but generated symmetrically for nicer asymmetries!!
c limits:
      do i=1,2
        xl(i)=0.d0
        xu(i)=1.d0
      end do
c set number of cos theta points to sum over
      nctpoints = 200
      else 
c VEGAS parameters.
      ndim=3
c if nprn<0 no print-out.
      nprn=0
c integrates on:
c x(1)=(x1-tau)/(1-tau),
c x(2)=(ecm-rm3-rm4)/(ecm_max-rm3-rm4),
c x(3)=ct
c limits:
      do i=1,2
        xl(i)=0.d0
        xu(i)=1.d0
      end do
      do i=3,3
        xl(i)=-1.d0
        xu(i)=1.d0
      end do
c set number of cos theta points to sum over
      nctpoints = 0
      end if
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c binning distributions
      include 'dists_bin.inc'
      if(m_sigp.eq.1)then
c sigp bin size.
        sigpw=(sigpmax-sigpmin)/ndiv_sigp
c generate bin in sigp.
        do i=1,ndiv_sigp
          xsigp(i)=sigpmin+sigpw*(i-1)+sigpw/2.d0
        end do
      end if
      if(m_sigm.eq.1)then
c bin sigm bin size.
        sigmw=(sigmmax-sigmmin)/ndiv_sigm
c generate bin in sigm.
        do i=1,ndiv_sigm
          xsigm(i)=sigmmin+sigmw*(i-1)+sigmw/2.d0
        end do
      end if
c output information.
  633 FORMAT(A8,F10.2,A10) 
  634 FORMAT(A8,10F10.2,A10) 
      print*,'#########################################################'
      write(jh,*)'Cross section of the process:'
      if(icoll.eq.0)write(jh,*)'pp     -> ',ff,'-bar'
      if(icoll.eq.1)write(jh,*)'pp-bar -> ',ff,'-bar'
      write(jh,*)'with unpolarised beam(s)'
      write(jh,*)'(quarks except top all strictly massless !)'
      write(jh,633)'rmb =',rmb,'    (GeV)'
      write(jh,633)'rmt =',rmt,'    (GeV)'
      write(jh,633)'gamt =',gamt,'    (GeV)'
      write(jh,633)'rmZ =',rmZ,'    (GeV)'
      write(jh,633)'gamZ =',gamZ,'    (GeV)'
      write(jh,633)'rmW =',rmW,'    (GeV)'
      write(jh,633)'gamW =',gamW,'    (GeV)'
      write(jh,633)'rmH =',rmH,'    (GeV)'
      write(jh,633)'gamH =',gamH,'    (GeV)'
      write(jh,634)'rmZZ =',rmZZ,'    (GeV)'
      write(jh,634)'gamZZ =',gamZZ,'    (GeV)'
      write(jh,633)'rmnuH =',NeuMass,'    (GeV)'
      write(jh,"(A13,F10.2,A10)")'at sqrt(s) =',ecm_coll,'    (GeV)'
      write(jh,"(A9,10F10.2)")'param  =',param
      IF(ISTRUCTURE.EQ.1)WRITE(jh,*)'PDFs from CTEQ6M'
      IF(ISTRUCTURE.EQ.2)WRITE(jh,*)'PDFs from CTEQ6D'
      IF(ISTRUCTURE.EQ.3)WRITE(jh,*)'PDFs from CTEQ6L'
      IF(ISTRUCTURE.EQ.4)WRITE(jh,*)'PDFs from CTEQ6L1'
      IF(QQ.gt.0.d0)THEN
            WRITE(jh,633)'QQ=',QQ,' (GeV)'
      ELSE IF(QQ.eq.-1.d0)THEN
            WRITE(jh,*)'QQ= ecm'
      ELSE 
            IF(ff.eq.'tt')THEN
            WRITE(jh,*)'QQ= 2*rmt'
            ELSE 
            WRITE(jh,*)'QQ= rmZ'
            ENDIF
      ENDIF
      WRITE(jh,"(A16,F6.2,A6)")'LAMBDA_QCD(4)= ',rlambdaQCD4,' (GeV)'
      write(jh,"(A23,I1,A8)")'with a_s evaluated at ',nloops,' loop(s)'   
      write(jh,"(A11,F7.5)")'a_s(MZ) = ',alfas(rmZ,rlambdaQCD4,nloops)
      write(jh,"(A6,F5.1)")'|y| <',abs(y)
      write(jh,"(A7,F7.1)")'|dy| <',dycut
      write(jh,"(A5,F5.0)")'pT >',ptCut
      write(jh,*)'OTHER SETTINGS'
      write(jh,*)'.mdl file = ', model
677   FORMAT(' dsigma integrated over: ',F6.1,' < M',A2,' < ',F6.1)
      write(jh,677)rmassmin,ff,rmassmax
      write(jh,"(A13,I5)")'nctpoints = ', nctpoints
      write(jh,*)'Contributions:'
      if(QCD.eq.1)then 
            write(jh,*)'    QCD: YES'
      else if(QCD.eq.0)then
            write(jh,*)'    QCD: NO'
      endif
      if(EW.eq.1)then 
            write(jh,*)'    EW: YES'
      else if(EW.eq.0)then
            write(jh,*)'    EW: NO'
      endif
      if(BSM.eq.1)then 
            write(jh,*)'    BSM: YES'
      else if(BSM.eq.0)then
            write(jh,*)'    BSM: NO'
      endif

      if(CONT.eq.0)then
            write(jh,*)'    EW/BSM:ALL'      
      else if(CONT.eq.1)then
            write(jh,*)'    EW/BSM:INT ONLY'
      else if(CONT.eq.2)then
            write(jh,*)'    EW/BSM:SQ ONLY'
      else if(CONT.eq.3)then
            write(jh,*)'    EW/BSM:NONE'
      endif
      print*,'################### START INTEGRATION ###################'  
c reset counter.
      npoints=0
c      kpts=0
c reset various iterative quantities.
      do i=1,20
        resl(i)=0.d0
        standdevl(i)=0.d0
        cnorm(i)=0.d0
        do iphel=-1,+1,2
          do jphel=-1,+1,2
            polcross(i,iphel,jphel)=0.d0
            polerror(i,iphel,jphel)=0.d0
          end do
        end do
        do iasy=-1,+1,2
          do j = 1,nasy
            asycross(i,iasy,j)=0.d0
            asyerror(i,iasy,j)=0.d0
          end do
        end do 
      end do
c starts integrations.
      it=0
c calls different function depending on cos x sum switch
      if(ISUMCOSX.eq.1)then 
c            call vegas(ndim,fxncosx,avgi,sd,chi2a)
            call vegas(ndim,fxn,avgi,sd,chi2a)
      else 
            call vegas(ndim,fxn,avgi,sd,chi2a)
      endif
      print*,'#################### END INTEGRATION ####################'
c collect total cross-section.
      cross=avgi
      error=sd
      xevents = alumpb*cross
      clup = (xevents + dsqrt(xevents))/alumpb
      cldn = (xevents - dsqrt(xevents))/alumpb         
c re-weight distributions for different it.
      stantot=0.d0
      do i=1,it
        stantot=stantot+1.d0/standdevl(i)/standdevl(i)
      end do
      do i=1,it
        standdevl(i)=standdevl(i)*standdevl(i)*stantot
      end do
      do i=1,it
        cnorm(i)=resl(i)*standdevl(i)
      end do
c collect polarised cross sections.
      do iphel=-1,+1,2
        do jphel=-1,+1,2
          do i=1,it
            polcross(i,iphel,jphel)=polcross(i,iphel,jphel)
     &                             *avgi/cnorm(i)
            polerror(i,iphel,jphel)=polcross(i,iphel,jphel)
     &                             *sd/cnorm(i)
          end do
          poltot(iphel,jphel)=0.d0
          polchi(iphel,jphel)=0.d0
          do i=1,it
            poltot(iphel,jphel)=poltot(iphel,jphel)
     &                     +polcross(i,iphel,jphel)
            polchi(iphel,jphel)=polchi(iphel,jphel)
     &                     +polerror(i,iphel,jphel)
          end do
          polchi(iphel,jphel)=polchi(iphel,jphel)
     &                       /poltot(iphel,jphel)
        end do
      end do
c collect unpolarised FB asymmetry.
      do iasy=-1,+1,2
        do i=1,it
          do j=1,nasy
            asycross(i,iasy,j)=asycross(i,iasy,j)*avgi/cnorm(i)
            asyerror(i,iasy,j)=asycross(i,iasy,j)*sd/cnorm(i)
          end do
        end do
        do j=1,nasy
          asytot(iasy,j)=0.d0
          asychi(iasy,j)=0.d0
        end do
        do i=1,it
          do j=1,nasy
            asytot(iasy,j)=asytot(iasy,j)+asycross(i,iasy,j)
            asychi(iasy,j)=asychi(iasy,j)+asyerror(i,iasy,j)
          enddo
        end do
        do j=1,nasy
          asychi(iasy,j)=asychi(iasy,j)/asytot(iasy,j)
        enddo
        
      end do
c define asymmetries.
      asy(1) = ! ALL
     &       +(poltot(+1,+1)-poltot(+1,-1)
     &       - poltot(-1,+1)+poltot(-1,-1))
     &        /cross
      asyerr(1) =
     &          +(polchi(+1,+1)+polchi(+1,-1)
     &          + polchi(-1,+1)+polchi(-1,-1))
     &           /4.d0*asy(1)
      asy(2) = ! AL 
     &       +(poltot(-1,-1)-poltot(+1,-1) 
     &       + poltot(-1,+1)-poltot(+1,+1))
     &        /cross
      asyerr(2) =
     &          +(polchi(-1,-1)+polchi(+1,-1) 
     &          + polchi(-1,+1)+polchi(+1,+1))
     &          /4.d0*asy(2)
      asy(3) = ! APV
     &       +(poltot(-1,-1)-poltot(+1,+1))
     &       /cross/2.d0
      asyerr(3) = 
     &          +(polchi(-1,-1)+polchi(+1,+1))
     &          /2.d0*asy(3)
      do j=4,nasy
        asy(j)= (asytot(+1,j)-asytot(-1,j))
     &              /(asytot(+1,j)+asytot(-1,j))
        asyerr(j)= sd/avgi*asy(j)
      enddo
C 
 964  FORMAT(F10.2,E15.5)
      mff = 'M_{'//ff//'}'
c plotting distributions
      include 'dists_plot.inc'
C plotting asymmetries
      include 'asy_plot.inc'
 1508 FORMAT((F10.4,F25.2))
c write all out.
      write(jh,*)'(using ',npoints,' points)'
      write(jh,*)'                                       '
      write(jh,*)'sigma (pb):           error (same units):'
      write(jh,'(2E15.5)')cross,error 
      write(jh,*)'ALL:                  error (same units):'
      write(jh,1508)asy(1),asyerr(1) 
      write(jh,*)'AL:                   error (same units):'
      write(jh,1508)asy(2),asyerr(2)
      write(jh,*)'APV:                  error (same units):'
      write(jh,1508)asy(3),asyerr(3)
      write(jh,*)'AFB:                  error (same units):'
      write(jh,1508)asy(4),asyerr(4)
      write(jh,*)'AFB2:                 error (same units):'
      write(jh,1508)asy(7),asyerr(7)
      write(jh,"(A4,F3.1,A2,A33)")'AC(',yc,'):','error (same units):'
      write(jh,1508)asy(5),asyerr(5)
      write(jh,"(A4,F3.1,A2,A33)")'AF(',yc,'):','error (same units):'
      write(jh,1508)asy(6),asyerr(6)
      write(jh,"(A6,F3.1,A2,A31)")'ARFB(',yffcut,'):',
     &                                          'error (same units):'
      write(jh,1508)asy(8),asyerr(8)
      write(jh,"(A6,F4.0,A2,A30)")'AOFB(',pzffcut,'):',
     &                                          'error (same units):'
      write(jh,1508)asy(9),asyerr(9)
      write(jh,*)'AFBSTAR:              error (same units):'
      write(jh,1508)asy(10),asyerr(10)
      do j=1,8
        write(jh,"(A9,F3.1,A2,A28)")'AFBSTAR(',yffs(j),'):'
     &                                          ,'error (same units):'
        write(jh,1508)asy(10+j),asyerr(10+j)
      enddo
cc testing:
      stop
      end