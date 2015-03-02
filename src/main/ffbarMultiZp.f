c      program d_pptt_LO
      implicit real*8 (a-h,o-z)
      common/bveg1/ncall,itmx,nprn,ndev,xl(100),xu(100),acc
      common/bveg2/it,ndo,si,swgt,schi,xi(50,100)
      common/rndm/iseed
      common/reslocal/resl(20),standdevl(20)
      common/par/s
      common/limfac/fac
      common/EW/a_em,s2w,zc_norm
      common/W/rmW,gamW  
      common/Z/rmZ,gamZ  
      common/t/rmt,gamt
      common/b/rmb
      common/H/rmH,gamH
      common/stat/npoints,nctpoints
      common/coll/ecm_coll
      common/ptr/y,ptCut,dycut
      COMMON/PARTDIST/ISTRUCTURE
      COMMON/PDF/QQ
      COMMON/ALFASTRONG/rlambdaQCD4,nloops
      common/collider/icoll
      common/errors/ierr
c test dummies
      common/test/acbinp(20),acbinm(20),antiacbinp(20),antiacbinm(20)
c sm coupling declarations
      REAL*8         GAL(2),GAU(2),GAD(2),GWF(2)
      COMMON /COUP2A/GAL,   GAU,   GAD,   GWF
      REAL*8         GZN(2),GZL(2),GZU(2),GZD(2),G1(2)
      COMMON /COUP2B/GZN,   GZL,   GZU,   GZD,   G1
c zprime model params
      common/pars/param(10),nzp
      common/ZZ/rmZZ(10),gamZZ(10),gp(10)
      common/VA/ggvu(10),ggau(10),ggvd(10),ggad(10),
     &          ggvt(10),ggat(10),ggvb(10),ggab(10),
     &          ggve(10),ggae(10),ggvn(10),ggan(10),
     &          ggvta(10),ggata(10),ggvnt(10),ggant(10)
      common/LR/gZZu(2,10),gZZd(2,10),gZZe(2,10),gZZn(2,10),
     &          gZZc(2,10),gZZs(2,10),gZZt(2,10),gZZb(2,10),
     &          gZZta(2,10),gZZnt(2,10)
c Possible corrections to Z couplings
      common/ZCORR/DgZu(2),DgZd(2),DgZe(2),DgZn(2),
     &             DgZt(2),DgZb(2),DgZta(2),DgZnt(2) 
c Lorenzos extra vars
      real*8 NeuMass
      common/neu/NeuMass
      common/neutype/iNeu 
c few extra variables for calculating CL, invariant mass cutt
      common/lum/alumpb
      common/cut/ainvmasscut
      common/msloop2/ndist
c contribution and QCD switches
      INTEGER CONT,QCD,EW,BSM
      common/switch/CONT,QCD,EW,BSM
c final state switches, parameters
      CHARACTER*2 ff
      CHARACTER*1 f
      common/fstr/ff,f
      common/fermion/rmf
      common/fints/jf,mf
c distributions.
      common/ext_pT/pTmax,pTmin,pTw
      common/dist_pT/xpT(500),fxpT(500,20),fxpTtot(500)
      common/inp_pT/m_pT
      common/div_pT/ndiv_pT
      common/ext_rmass/rmassmax,rmassmin,rmassw
      common/dist_rmass/xrmass(500),fxrmass(500,20),fxrmasstot(500)
      common/inp_rmass/m_rmass
      common/div_rmass/ndiv_rmass
      common/ext_eta/etamax,etamin,etaw
      common/dist_eta3/xeta3(500),fxeta3(500,20),fxeta3tot(500)
      common/dist_eta4/xeta4(500),fxeta4(500,20),fxeta4tot(500)
      common/dist_eta3col/xeta3col(500),fxeta3col(500,20)
     &       ,fxeta3coltot(500)
      common/dist_eta4col/xeta4col(500),fxeta4col(500,20)
     &       ,fxeta4coltot(500)
      common/inp_eta/m_eta
      common/div_eta/ndiv_eta
      common/ext_beta/betamax,betamin,betaw
      common/dist_beta/xbeta(500),fxbeta(500,20),fxbetatot(500)
      common/inp_beta/m_beta
      common/div_beta/ndiv_beta
      common/ext_cost/costmax,costmin,costw
      common/dist_cost/xcost3(500),fxcost3(500,20),fxcost3tot(500)
      common/dist_cost/xcost4(500),fxcost4(500,20),fxcost4tot(500)
      common/dist_cost/xcost3col(500),fxcost3col(500,20)
     &       ,fxcost3coltot(500)
      common/dist_cost/xcost4col(500),fxcost4col(500,20)
     &       ,fxcost4coltot(500)
      common/inp_cost/m_cost
      common/div_cost/ndiv_cost
      common/ext_Et/Etmax,Etmin,Etw
      common/dist_Et/xEt(500),fxEt(500,20),fxEttot(500)
      common/inp_Et/m_Et
      common/div_Et/ndiv_Et
cc Chi distributions      
      common/ext_Chi/Chimax,Chimin,Chiw
      common/dist_Chi/xChi(500),fxChi(500,20),fxChitot(500)
      common/inp_Chi/m_Chi
      common/div_Chi/ndiv_Chi   
      common/ext_Chi2/Chi2max,Chi2min,Chi2w
      common/dist_Chi2/xChi2(500),fxChi2(500,20),fxChi2tot(500)
      common/inp_Chi2/m_Chi2
      common/div_Chi2/ndiv_Chi2 
cc Rapidity distributions
      common/ext_Y/Ymax,Ymin,Yw
      common/dist_Y3/xY3(500),fxY3(500,20),fxY3tot(500)
      common/dist_Y4/xY4(500),fxY4(500,20),fxY4tot(500)
      common/inp_Y/m_Y
      common/div_Y/ndiv_Y 
      common/ext_Ycol/Ycolmax,Ycolmin,Ycolw
      common/dist_Ycol3/xYcol3(500),fxYcol3(500,20),fxYcol3tot(500)
      common/dist_Ycol4/xYcol4(500),fxYcol4(500,20),fxYcol4tot(500)
      common/inp_Ycol/m_Ycol
      common/div_Ycol/ndiv_Ycol 
      common/ext_DelY/DelYmax,DelYmin,DelYw
      common/dist_DelY/xDelY(500),fxDelY(500,20),fxDelYtot(500)
      common/inp_DelY/m_DelY
      common/div_DelY/ndiv_DelY 
      common/ext_Ytt/Yttmax,Yttmin,Yttw
      common/dist_Ytt/xYtt(500),fxYtt(500,20),fxYtttot(500)
      common/inp_Ytt/m_Ytt
      common/div_Ytt/ndiv_Ytt  
      common/ext_Pztt/Pzttmax,Pzttmin,Pzttw
      common/dist_Pztt/xPztt(500),fxPztt(500,20),fxPztttot(500)
      common/inp_Pztt/m_Pztt
      common/div_Pztt/ndiv_Pztt 
cc x1/x2 distribution
      common/ext_x1x2/x1x2max,x1x2min,x1x2w
      common/dist_x1x2/xx1x2(500),fxx1x2(500,20),fxx1x2tot(500)
      common/inp_x1x2/m_x1x2
      common/div_x1x2/ndiv_x1x2 
cc Asymmetries      
      common/ext_sigp/sigpmax,sigpmin,sigpw
      common/dist_sigp/xsigp(1000),fxsigp(20,1000,20),fxsigptot(20,1000)
      common/inp_sigp/m_sigp
      common/div_sigp/ndiv_sigp
      common/ext_sigm/sigmmax,sigmmin,sigmw
      common/dist_sigm/xsigm(1000),fxsigm(20,1000,20),fxsigmtot(20,1000)
      common/inp_sigm/m_sigm
      common/div_sigm/ndiv_sigm
      common/nfnb/fb(1000,-1:1,5)
      common/polarised/polcross(20,-1:1,-1:1),polerror(20,-1:1,-1:1)
      common/FB/asycross(20,-1:1),asyerror(20,-1:1)
      common/FB2/asy2cross(20,-1:1),asy2error(20,-1:1)
      common/ACENT/accross(20,-1:1),acerror(20,-1:1)
      common/AFORW/afcross(20,-1:1),aferror(20,-1:1)
      common/ARFB/arfbcross(20,-1:1),arfberror(20,-1:1)
      common/AOFB/aofbcross(20,-1:1),aofberror(20,-1:1)
      common/AFBST/afbstcross(20,-1:1),afbsterror(20,-1:1)
      common/AFBST1/afbst1cross(20,-1:1,8),afbst1error(20,-1:1,8)    
      common/ACENTcut/yc,yff,pff
      common/AFBSTcut/yffs(8)
      dimension cnorm(20)
      dimension snorm(20),ave(20)
      dimension poltot(-1:1,-1:1),polchi(-1:1,-1:1)
      dimension asytot(-1:1),asychi(-1:1)
      dimension asy2tot(-1:1),asy2chi(-1:1)
      dimension actot(-1:1),acchi(-1:1)
      dimension aftot(-1:1),afchi(-1:1)
      dimension arfbtot(-1:1),arfbchi(-1:1)
      dimension aofbtot(-1:1),aofbchi(-1:1)
      dimension afbsttot(-1:1),afbstchi(-1:1)
      dimension afbst1tot(-1:1,8),afbst1chi(-1:1,8)
      dimension AFBSTCUT(8),AFBSTCUTerr(8)
      real*8 widthmatrix(2,2)
      parameter (pi=3.14159265358979323846d0)
      parameter (conv=.38937966d9)
c     filenames
      CHARACTER*4 name
      CHARACTER*1 intstr
c      CHARACTER*4 mzp
      CHARACTER*50 model
      common/mod/model,imodel
      external fxn,fxncosx
c gauge boson masses and widths.
      data gamW/2.08d0/
      data rmb/4.5d0/
      data rmZ/91.19d0/,gamZ/2.50d0/
      data rmH/150.d0/,gamH/0.1635D-01/
      data ierr/0/
      data nzp/10/
      data iNeu/0/
      data NeuMass/0.d0/
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
      pTmax=rmassmax/2.d0     
      Etmax=rmassmax/2.d0
      m_pT=ndist
      pTmin=0.d0
      ndiv_pT=175
      m_rmass=ndist
      rmassmin=ainvmasscut
      ndiv_rmass=200
      m_eta=ndist
      etamax=+5.d0
      etamin=-5.d0
      ndiv_eta=50
      m_beta=ndist
      betamax=+1.d+3
      betamin=+1.d-3
      ndiv_beta=100
      m_cost=ndist
      costmax=+1.d0
      costmin=-1.d0
      ndiv_cost=50
      m_Et=ndist
      Etmin=0.d0
      ndiv_Et=175    
      
      m_Chi=ndist
      Chimin=1.d0
      Chimax=30.d0
      ndiv_Chi=100
      m_Chi2=ndist
      Chi2min=1.d0
      Chi2max=Chimax
      ndiv_Chi2=100

      m_Y=ndist
      Ymin=-4.d0
      Ymax=4.d0
      ndiv_Y=100
      m_Ycol=ndist
      Ycolmin=-4.d0
      Ycolmax=4.d0
      ndiv_Ycol=100
      m_DelY=ndist
      DelYmin=-4.d0
      DelYmax=4.d0
      ndiv_DelY=100
    
      m_x1x2=ndist
      x1x2max=20d0
      x1x2min=0d0
      ndiv_x1x2=100
    
      m_Ytt=ndist
      Yttmin=-4.d0
      Yttmax=4.d0
      ndiv_Ytt=100
      
      m_Pztt=ndist
      Pzttmin=0.d0
      Pzttmax=ecm_coll/2.d0*(1.d0 - rmassmin*rmassmin/ecm_coll/ecm_coll)
      ndiv_Pztt=100

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
      read(42,*) ggve
      read(42,*) ggae
      read(42,*) ggvn
      read(42,*) ggan
c some EW parameters.
      a_em=1.d0/128.d0
      s2w=.2320d0
      rmW=rmZ*sqrt(1.d0-s2w)
      sw=sqrt(s2w)
      c2w=1.d0-s2w
      cw=sqrt(c2w)
      em_ch=sqrt(4.d0*pi*a_em)
      zc_norm=em_ch/sw/cw/2.d0

c default top mass/width
      rmt=175.d0
      gamt=1.55d0
c iuniv selects universal (0) or non-universal (1) couplings
      iuniv=0      
      CALL SMCORR(rmt,gamt,rmb,rmW,rmH,
     &            DgZu,DgZd,DgZe,DgZn,DgZt,DgZb,DgZta,DgZnt,iuniv)
c initialize MadGraph for QCD MEs.     
ccccccccccccccccccccccccccc
c SET FINAL STATE
      dycut=1000000.d0
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
cc Chi limits for ttbar
      qmax2 = ((rmassmax*rmassmax-2.d0*rmf*rmf)**2-(2.d0*rmf*rmf)**2)/
     &    (4.d0*rmassmax*rmassmax)
      qmax = dsqrt(qmax2)
      Emax = dsqrt(qmax2+rmf*rmf)
      upperChi = int((Emax+qmax)/(Emax-qmax))+1
      if(upperChi.lt.Chimax)then
            Chimax=upperChi
            Chi2max=upperChi
      end if
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
      CALL ZPCOUP(rmZZ,gp,ggvu,ggau,ggvd,ggad,
     & ggve,ggae,ggvn,ggan,ggvt,ggat,ggvb,ggab,ggvta,ggata,ggvnt,ggant,
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
      yff = 0.5d0
      do i=1,8
         yffs(i)=0.2d0*i
      end do
      pff = 700.d0
c compare with LB code
c      neumass=1.d5
c      ineu = 0
c      gamZZalt = Zpwidth(rmW,rmZ,rmZZ,a_em,s2w,rlambdaQCD4,nloop) 
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
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      
c binning distributions
      if(m_pT.eq.1)then
c pT bin size.
        pTw=(pTmax-pTmin)/ndiv_pT
c generate bin in pT.
        do i=1,ndiv_pT
          xpT(i)=pTmin+pTw*(i-1)+pTw/2.d0
        end do
      end if
      if(m_rmass.eq.1)then
c rmass bin size.
        rmassw=(rmassmax-rmassmin)/ndiv_rmass
c generate bin in rmass.
        do i=1,ndiv_rmass
          xrmass(i)=rmassmin+rmassw*(i-1)+rmassw/2.d0
        end do
      end if
      if(m_eta.eq.1)then
c eta bin size.
        etaw=(etamax-etamin)/ndiv_eta
c generate bin in eta.
        do i=1,ndiv_eta
          xeta3(i)=etamin+etaw*(i-1)+etaw/2.d0
          xeta4(i)=xeta3(i)
          xeta3col(i)=xeta3(i)
          xeta4col(i)=xeta3(i)
        end do
      end if
      if(m_beta.eq.1)then
c beta bin size.
        betaw=(betamax-betamin)/ndiv_beta
c generate bin in beta.
        do i=1,ndiv_beta
          xbeta(i)=betamin+betaw*(i-1)+betaw/2.d0
        end do
      end if
      if(m_cost.eq.1)then
c cost bin size.
        costw=(costmax-costmin)/ndiv_cost
c generate bin in cost.
        do i=1,ndiv_cost
          xcost3(i)=costmin+costw*(i-1)+costw/2.d0
          xcost4(i)=xcost3(i)
          xcost3col(i)=xcost3(i)
          xcost4col(i)=xcost3(i)
        end do
      end if
      if(m_Et.eq.1)then
c Et bin size.
        Etw=(Etmax-Etmin)/ndiv_Et
c generate bin in Et.
        do i=1,ndiv_Et
          xEt(i)=Etmin+Etw*(i-1)+Etw/2.d0
        end do
      end if
      
      if(m_Chi.eq.1)then
c Chi bin size.
        Chiw=(Chimax-Chimin)/ndiv_Chi
c generate bin in Chi.
        do i=1,ndiv_Chi
          xChi(i)=Chimin+Chiw*(i-1)+Chiw/2.d0
        end do
      end if
      if(m_Chi2.eq.1)then
c Chi2 bin size.
        Chi2w=(Chi2max-Chi2min)/ndiv_Chi2
c generate bin in Chi2.
        do i=1,ndiv_Chi2
          xChi2(i)=Chi2min+Chi2w*(i-1)+Chi2w/2.d0
        end do
      end if

      if(m_Y.eq.1)then
c Y3,4 bin size.
        Yw=(Ymax-Ymin)/ndiv_Y
c generate bin in Y3.
        do i=1,ndiv_Y
          xY3(i)=Ymin+Yw*(i-1)+Yw/2.d0
          xY4(i)=xY3(i)
        end do
      end if
      if(m_Ycol.eq.1)then
c Ycol3,4 bin size.
        Ycolw=(Ycolmax-Ycolmin)/ndiv_Ycol
c generate bin in Ycol3.
        do i=1,ndiv_Ycol
          xYcol3(i)=Ycolmin+Ycolw*(i-1)+Ycolw/2.d0
          xYcol4(i)=xYcol3(i)
        end do
      end if
      if(m_DelY.eq.1)then
c DelY bin size.
        DelYw=(DelYmax-DelYmin)/ndiv_DelY
c generate bin in DelY.
        do i=1,ndiv_DelY
          xDelY(i)=DelYmin+DelYw*(i-1)+DelYw/2.d0
        end do
      end if

      if(m_Ytt.eq.1)then
c Ytt bin size.
        Yttw=(Yttmax-Yttmin)/ndiv_Ytt
c generate bin in Ytt.
        do i=1,ndiv_Ytt
          xYtt(i)=Yttmin+Yttw*(i-1)+Yttw/2.d0
        end do
      end if
      
      if(m_Pztt.eq.1)then
c Pztt bin size.
        Pzttw=(Pzttmax-Pzttmin)/ndiv_Pztt
c generate bin in Pztt.
        do i=1,ndiv_Pztt
          xPztt(i)=Pzttmin+Pzttw*(i-1)+Pzttw/2.d0
        end do
      end if
      
      if(m_x1x2.eq.1)then
c x1x2 bin size. MAKE LOGSCALE
        x1x2w=(x1x2max-x1x2min)/ndiv_x1x2
c generate bin in x1x2.
        do i=1,ndiv_x1x2
          xx1x2(i)=x1x2min+x1x2w*(i-1)+x1x2w/2.d0
        end do
      end if
      
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
      print*,'######################'
      write(jh,*)'Cross section of the process:'
      if(icoll.eq.0)write(jh,*)'pp     -> ',ff,'-bar'
      if(icoll.eq.1)write(jh,*)'pp-bar -> ',ff,'-bar'
      write(jh,*)'with unpolarised beam(s)'
      write(jh,*)'(quarks except top all strictly massless !)'
      write(jh,*)'rmb =',rmb,'    (GeV)'
      write(jh,*)'rmt =',rmt,'    (GeV)'
      write(jh,*)'gamt =',gamt,'    (GeV)'
      write(jh,*)'rmZ =',rmZ,'    (GeV)'
      write(jh,*)'gamZ =',gamZ,'    (GeV)'
      write(jh,*)'rmW =',rmW,'    (GeV)'
      write(jh,*)'gamW =',gamW,'    (GeV)'
      write(jh,*)'rmH =',rmH,'    (GeV)'
      write(jh,*)'gamH =',gamH,'    (GeV)'
      write(jh,*)'rmZZ =',rmZZ,'    (GeV)'
      write(jh,*)'gamZZ =',gamZZ,'    (GeV)'
      write(jh,*)'rmnuH =',NeuMass,'    (GeV)'
      write(jh,*)'at sqrt(s)  =',ecm_coll,'    (GeV)'
      write(jh,*)'param  =',param
      IF(ISTRUCTURE.EQ.1)WRITE(jh,*)'PDFs from CTEQ6M'
      IF(ISTRUCTURE.EQ.2)WRITE(jh,*)'PDFs from CTEQ6D'
      IF(ISTRUCTURE.EQ.3)WRITE(jh,*)'PDFs from CTEQ6L'
      IF(ISTRUCTURE.EQ.4)WRITE(jh,*)'PDFs from CTEQ6L1'
      IF(QQ.gt.0.d0)THEN
            WRITE(jh,*)'QQ=',QQ,' (GeV)'
      ELSE IF(QQ.eq.-1.d0)THEN
            WRITE(jh,*)'QQ= ecm'
      ELSE 
            IF(ff.eq.'tt')THEN
            WRITE(jh,*)'QQ= 2*rmt'
            ELSE 
            WRITE(jh,*)'QQ= rmZ'
            ENDIF
      ENDIF
      WRITE(jh,*)'LAMBDA_QCD(4)=',rlambdaQCD4,' (GeV)'
      write(jh,*)'with a_s evaluated at',nloops,' loop(s)'   
      write(jh,*)'a_s(MZ) =',alfas(rmZ,rlambdaQCD4,nloops)
      write(jh,*)'|y| <',abs(y)
      write(jh,*)'|dy| <',dycut
      write(jh,*)'pT >',ptCut
      write(jh,*)'OTHER SETTINGS'
      write(jh,*)'.mdl file = ', model
      write(jh,*)'dsigma integrated over: ',rmassmin,' < M',ff,' < '
     & ,rmassmax
      write(jh,*)'nctpoints = ', nctpoints
      write(jh,*)'Contributions:'
      if(QCD.eq.1)then 
            write(jh,*)'QCD: YES'
      else if(QCD.eq.0)then
            write(jh,*)'QCD: NO'
      endif
      if(EW.eq.1)then 
            write(jh,*)'EW: YES'
      else if(EW.eq.0)then
            write(jh,*)'EW: NO'
      endif
      if(BSM.eq.1)then 
            write(jh,*)'BSM: YES'
      else if(BSM.eq.0)then
            write(jh,*)'BSM: NO'
      endif

      if(CONT.eq.0)then
            write(jh,*)'EW/BSM:ALL'      
      else if(CONT.eq.1)then
            write(jh,*)'EW/BSM:INT ONLY'
      else if(CONT.eq.2)then
            write(jh,*)'EW/BSM:SQ ONLY'
      else if(CONT.eq.3)then
            write(jh,*)'EW/BSM:NONE'
      endif
      print*,'########### START INTEGRATION ###########'  
c reset counter.
      npoints=0
c      kpts=0
      do j=1,1000
            do i=1,it
                  fb(j,0,i)=0.d0
            end do
      end do
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
          asycross(i,iasy)=0.d0
          asyerror(i,iasy)=0.d0
          asy2cross(i,iasy)=0.d0
          asy2error(i,iasy)=0.d0
          accross(i,iasy)=0.d0
          acerror(i,iasy)=0.d0
          afcross(i,iasy)=0.d0
          aferror(i,iasy)=0.d0
          arfbcross(i,iasy)=0.d0
          arfberror(i,iasy)=0.d0
          aofbcross(i,iasy)=0.d0
          aofberror(i,iasy)=0.d0
          afbstcross(i,iasy)=0.d0
          afbsterror(i,iasy)=0.d0
          do j = 1,8
            afbst1cross(i,iasy,j)=0.d0
            afbst1error(i,iasy,j)=0.d0
          end do
        end do 
      end do
cc initialize test dummmies
c      do j = 1,20
c      acbinp(j) = 0.d0
c      acbinm(j) = 0.d0
c      antiacbinp(j) = 0.d0
c      antiacbinm(j) = 0.d0
c      end do
c starts integrations.
      it=0
c      print*, 'calling vegas'
c calls different function depending on cos x sum switch
      if(ISUMCOSX.eq.1)then 
C             print*,"calling VEGAS : ndim,avgi,sd,chi2a = ",
C      &      ndim,avgi,sd,chi2a
            call vegas(ndim,fxncosx,avgi,sd,chi2a)
      else 
C             print*,"calling VEGAS : ndim,avgi,sd,chi2a = ",
C      &     ndim,avgi,sd,chi2a
            call vegas(ndim,fxn,avgi,sd,chi2a)
      endif
      print*,'########### END INTEGRATION ###########'
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
c collect unpolarised asymmetries.
      do iasy=-1,+1,2
        do i=1,it
          asycross(i,iasy)=asycross(i,iasy)*avgi/cnorm(i)
          asyerror(i,iasy)=asycross(i,iasy)*sd/cnorm(i)
          asy2cross(i,iasy)=asy2cross(i,iasy)*avgi/cnorm(i)
          asy2error(i,iasy)=asy2cross(i,iasy)*sd/cnorm(i)
          accross(i,iasy)=accross(i,iasy)*avgi/cnorm(i)
          acerror(i,iasy)=accross(i,iasy)*sd/cnorm(i)
          afcross(i,iasy)=afcross(i,iasy)*avgi/cnorm(i)
          aferror(i,iasy)=afcross(i,iasy)*sd/cnorm(i)
          arfbcross(i,iasy)=arfbcross(i,iasy)*avgi/cnorm(i)
          arfberror(i,iasy)=arfbcross(i,iasy)*sd/cnorm(i)
          aofbcross(i,iasy)=aofbcross(i,iasy)*avgi/cnorm(i)
          aofberror(i,iasy)=aofbcross(i,iasy) *sd/cnorm(i)
          afbstcross(i,iasy)=afbstcross(i,iasy)*avgi/cnorm(i)
          afbsterror(i,iasy)=afbstcross(i,iasy)*sd/cnorm(i)
          do j=1,8
            afbst1cross(i,iasy,j)=afbst1cross(i,iasy,j)*avgi/cnorm(i)
            afbst1error(i,iasy,j)=afbst1cross(i,iasy,j)*sd/cnorm(i)
          end do
        end do
        asytot(iasy)=0.d0
        asychi(iasy)=0.d0
        asy2tot(iasy)=0.d0
        asy2chi(iasy)=0.d0
        actot(iasy)=0.d0
        acchi(iasy)=0.d0
        aftot(iasy)=0.d0
        afchi(iasy)=0.d0
        arfbtot(iasy)=0.d0
        arfbchi(iasy)=0.d0
        aofbtot(iasy)=0.d0
        aofbchi(iasy)=0.d0
        afbsttot(iasy)=0.d0
        afbstchi(iasy)=0.d0
        do j=1,8
          afbst1tot(iasy,j)=0.d0
          afbst1chi(iasy,j)=0.d0
        end do
        do i=1,it
          asytot(iasy)=asytot(iasy)+asycross(i,iasy)
          asychi(iasy)=asychi(iasy)+asyerror(i,iasy)
          asy2tot(iasy)=asy2tot(iasy)+asy2cross(i,iasy)
          asy2chi(iasy)=asy2chi(iasy)+asy2error(i,iasy)
          actot(iasy)=actot(iasy)+accross(i,iasy)
          acchi(iasy)=acchi(iasy)+acerror(i,iasy)
          aftot(iasy)=aftot(iasy)+afcross(i,iasy)
          afchi(iasy)=afchi(iasy)+aferror(i,iasy)
          arfbtot(iasy)=arfbtot(iasy)+arfbcross(i,iasy)
          arfbchi(iasy)=arfbchi(iasy)+arfberror(i,iasy)
          aofbtot(iasy)=aofbtot(iasy)+aofbcross(i,iasy)
          aofbchi(iasy)=aofbchi(iasy)+aofberror(i,iasy)
          afbsttot(iasy)=afbsttot(iasy)+afbstcross(i,iasy)
          afbstchi(iasy)=afbstchi(iasy)+afbsterror(i,iasy)
          do j=1,8
            afbst1tot(iasy,j)=afbst1tot(iasy,j)+afbst1cross(i,iasy,j)
            afbst1chi(iasy,j)=afbst1chi(iasy,j)+afbst1error(i,iasy,j)
          enddo
        end do
        asychi(iasy)=asychi(iasy)/asytot(iasy)
        asy2chi(iasy)=asy2chi(iasy)/asy2tot(iasy)
        acchi(iasy)=acchi(iasy)/actot(iasy)
        afchi(iasy)=afchi(iasy)/aftot(iasy)
        arfbchi(iasy)=arfbchi(iasy)/arfbtot(iasy)
        aofbchi(iasy)=aofbchi(iasy)/aofbtot(iasy)
        afbstchi(iasy)=afbstchi(iasy)/afbsttot(iasy)
        do j=1,8
          afbst1chi(iasy,j)=afbst1chi(iasy,j)/afbst1tot(iasy,j)
        enddo
      end do
c define integrated asymmetries.
      ALL=     (poltot(+1,+1)-poltot(+1,-1)
     &         -poltot(-1,+1)+poltot(-1,-1))
     &        /cross
      ALLerr=  (polchi(+1,+1)+polchi(+1,-1)
     &         +polchi(-1,+1)+polchi(-1,-1))
     &           /4.d0*ALL
      AL=      (poltot(-1,-1)-poltot(+1,-1) 
     &         +poltot(-1,+1)-poltot(+1,+1))
     &        /cross
      ALerr=     (polchi(-1,-1)+polchi(+1,-1) 
     &           +polchi(-1,+1)+polchi(+1,+1))
     &          /4.d0*AL
      APV=    (poltot(-1,-1)-poltot(+1,+1))
     &        /cross/2.d0
      APVerr= (polchi(-1,-1)+polchi(+1,+1))
     &           /2.d0*APV
      AFB=    (asytot(+1)-asytot(-1))
     &        /cross
      AFBerr= sd/avgi*AFB
      AFB2= (asy2tot(+1)-asy2tot(-1))
     &        /cross
      AFB2err= sd/avgi*AFB2
cc corrected to not normalise to total xsect
      AC=      (actot(+1)-actot(-1))
     &        /(actot(+1)+actot(-1))
      ACerr=  sd/avgi*AC
      AF=      (aftot(+1)-aftot(-1))
     &        /(aftot(+1)+aftot(-1))
      AFerr=  sd/avgi*AF
      ARFB=    (arfbtot(+1)-arfbtot(-1))
     &        /(arfbtot(+1)+arfbtot(-1))
      ARFBerr=  sd/avgi*ARFB
      AOFB=    (aofbtot(+1)-aofbtot(-1))
     &        /(aofbtot(+1)+aofbtot(-1))
      AOFBerr=  sd/avgi*AOFB
      AFBST=   (afbsttot(+1)-afbsttot(-1))
     &        /(afbsttot(+1)+afbsttot(-1))
      AFBSTerr= sd/avgi*AFBST
      do j=1,8
        AFBSTCUT(j)= (afbst1tot(+1,j)-afbst1tot(-1,j))
     &              /(afbst1tot(+1,j)+afbst1tot(-1,j))
        AFBSTCUTerr(j)= sd/avgi*AFBSTCUT(j)
      enddo
c plotting distributions
      if(m_pT.eq.1)then
c plot distribution in pT.
        do j=1,ndiv_pT
          do i=1,it
            fxpT(j,i)=fxpT(j,i)*avgi/cnorm(i)/pTw
          end do
        end do
        do j=1,ndiv_pT
          fxpTtot(j)=0.d0 
          do i=1,it
            fxpTtot(j)=fxpTtot(j)+fxpT(j,i)
          end do
        end do
        write(jh,*)'DISTRIBUTION:'
        write(jh,*)'DSTRING: "pT"'
        write(jh,*)'TITLE: "p_{T} Distribution"'
        write(jh,*)'AXES: "p_{T} (GeV)"',
     &' "d{/Symbol s}/dp_{T} (pb/GeV)"'
        do i=1,ndiv_pT
          write(jh,*)xpT(i),fxpTtot(i)
        end do
        write(jh,*)'END'
      end if  
      
      if(m_rmass.eq.1)then
c plot distribution in rmass.
        do j=1,ndiv_rmass
          do i=1,it
            fxrmass(j,i)=fxrmass(j,i)*avgi/cnorm(i)/rmassw
          end do
        end do
        do j=1,ndiv_rmass
          fxrmasstot(j)=0.d0 
          do i=1,it
            fxrmasstot(j)=fxrmasstot(j)+fxrmass(j,i)
          end do
        end do
        write(jh,*)'DISTRIBUTION:'
        write(jh,*)'DSTRING: "M',ff,'"'
        write(jh,*)'TITLE: "M_{',ff,'} Distribution"'
        write(jh,*)'AXES: "M_{',ff,'} (GeV)"',
     &' "d{/Symbol s}/dM_{',ff,'} (pb/GeV)"'
        do i=1,ndiv_rmass
          write(jh,*)xrmass(i),fxrmasstot(i)
        end do
        write(jh,*)'END'        
      end if  
      
      if(m_Et.eq.1)then
c plot distribution in Et.
        do j=1,ndiv_Et
          do i=1,it
            fxEt(j,i)=fxEt(j,i)*avgi/cnorm(i)/Etw
          end do
        end do
        do j=1,ndiv_Et
          fxEttot(j)=0.d0 
          do i=1,it
            fxEttot(j)=fxEttot(j)+fxEt(j,i)
          end do
        end do
        write(jh,*)'DISTRIBUTION:'
        write(jh,*)'DSTRING: "E',f,'"'
        write(jh,*)'TITLE: "E_{',f,'} Distribution"'
        write(jh,*)'AXES: "E_{',f,'} (GeV)"', 
     &'"d{/Symbol s}/dE_{t} (pb/GeV)"'
        do i=1,ndiv_Et
          write(jh,*)xEt(i),fxEttot(i)
        end do
        write(jh,*)'END'
      end if
      
      if(m_eta.eq.1)then
c plot distribution in eta.
        do j=1,ndiv_eta
          do i=1,it
            fxeta3(j,i)=fxeta3(j,i)*avgi/cnorm(i)/etaw
            fxeta4(j,i)=fxeta4(j,i)*avgi/cnorm(i)/etaw
            fxeta3col(j,i)=fxeta3col(j,i)*avgi/cnorm(i)/etaw
            fxeta4col(j,i)=fxeta4col(j,i)*avgi/cnorm(i)/etaw
          end do
        end do
        do j=1,ndiv_eta
          fxeta3tot(j)=0.d0 
          fxeta4tot(j)=0.d0 
          fxeta3coltot(j)=0.d0 
          fxeta4coltot(j)=0.d0 
          do i=1,it
            fxeta3tot(j)=fxeta3tot(j)+fxeta3(j,i)
            fxeta4tot(j)=fxeta4tot(j)+fxeta4(j,i)
            fxeta3coltot(j)=fxeta3coltot(j)+fxeta3col(j,i)
            fxeta4coltot(j)=fxeta4coltot(j)+fxeta4col(j,i)
          end do
        end do
        write(jh,*)'DISTRIBUTION:'
        write(jh,*)'DSTRING: "eta3cm"'
        write(jh,*)'TITLE: "{/Symbol h}_{',f,'bar}(CM) Distribution"'
        write(jh,*)'AXES: "{/Symbol h}"', 
     &'"d{/Symbol s}/d{/Symbol h} (pb)"' 
        do i=1,ndiv_eta
          write(jh,*)xeta3(i),fxeta3tot(i)
        end do
        write(jh,*)'END'
        
        
        write(jh,*)'DISTRIBUTION:'
        write(jh,*)'DSTRING: "eta4cm"'
        write(jh,*)'TITLE: "{/Symbol h}_{',f,'}(CM) Distribution"'
        write(jh,*)'AXES: "{/Symbol h}"',
     &'"d{/Symbol s}/d{/Symbol h} (pb)"'
        do i=1,ndiv_eta
          write(jh,*)xeta4(i),fxeta4tot(i)
        end do
        write(jh,*)'END'
        
        
        write(jh,*)'DISTRIBUTION:'
        write(jh,*)'DSTRING: "eta3"'
        write(jh,*)'TITLE: "{/Symbol h}_{',f,'bar}', 
     &' Distribution"'
        write(jh,*)'AXES: "{/Symbol h}"',
     &' "d{/Symbol s}/d{/Symbol h} (pb)"'
        do i=1,ndiv_eta
          write(jh,*)xeta3col(i),fxeta3coltot(i)
        end do
        write(jh,*)'END'
        
        
        write(jh,*)'DISTRIBUTION:'
        write(jh,*)'DSTRING: "eta4"'
        write(jh,*)'TITLE: "{/Symbol h}_{',f,'} Distribution"'
        write(jh,*)'AXES: "{/Symbol h}"',
     &' "d{/Symbol s}/d{/Symbol h} (pb)"'
        do i=1,ndiv_eta
          write(jh,*)xeta4col(i),fxeta4coltot(i)
        end do
        write(jh,*)'END'
        
      end if  
      if(m_beta.eq.1)then
c plot distribution in beta.
        do j=1,ndiv_beta
          do i=1,it
            fxbeta(j,i)=fxbeta(j,i)*avgi/cnorm(i)/betaw
          end do
        end do
        do j=1,ndiv_beta
          fxbetatot(j)=0.d0 
          do i=1,it
            fxbetatot(j)=fxbetatot(j)+fxbeta(j,i)
          end do
        end do
        write(jh,*)'DISTRIBUTION:'
        write(jh,*)'DSTRING: "beta"'
        write(jh,*)'TITLE: "{/Symbol b} Distribution"'
        write(jh,*)'AXES: "{/Symbol b}"',
     &'"d{/Symbol s}/d{/Symbol b} (pb)"'
        do i=1,ndiv_beta
          write(jh,*)xbeta(i),fxbetatot(i)
        end do
        write(jh,*)'END'
        
      end if  
      if(m_cost.eq.1)then
c plot distribution in cost.
        do j=1,ndiv_cost
          do i=1,it
            fxcost3(j,i)=fxcost3(j,i)*avgi/cnorm(i)/costw
            fxcost4(j,i)=fxcost4(j,i)*avgi/cnorm(i)/costw
            fxcost3col(j,i)=fxcost3col(j,i)*avgi/cnorm(i)/costw
            fxcost4col(j,i)=fxcost4col(j,i)*avgi/cnorm(i)/costw
          end do
        end do
        do j=1,ndiv_cost
          fxcost3tot(j)=0.d0 
          fxcost4tot(j)=0.d0
          fxcost3coltot(j)=0.d0 
          fxcost4coltot(j)=0.d0 
          do i=1,it
            fxcost3tot(j)=fxcost3tot(j)+fxcost3(j,i)
            fxcost4tot(j)=fxcost4tot(j)+fxcost4(j,i)
            fxcost3coltot(j)=fxcost3coltot(j)+fxcost3col(j,i)
            fxcost4coltot(j)=fxcost4coltot(j)+fxcost4col(j,i)
          end do
        end do
        write(jh,*)'DISTRIBUTION:'
        write(jh,*)'DSTRING: "cos3cm"'
        write(jh,*)'TITLE: "cos{/Symbol q}_{',f,'bar}(CM)',
     &' Distribution"'
        write(jh,*)'AXES: "cos{/Symbol q}" 
     &"d{/Symbol s}/dcos{/Symbol q}} (pb)"'
        do i=1,ndiv_cost
          write(jh,*)xcost3(i),fxcost3tot(i)
        end do
        write(jh,*)'END'
        
        write(jh,*)'DISTRIBUTION:'
        write(jh,*)'DSTRING: "cos4cm"'
        write(jh,*)'TITLE: "cos{/Symbol q}_{',f,'}(CM) Distribution"'
        write(jh,*)'AXES: "cos{/Symbol q}"',
     &' "d{/Symbol s}/dcos{/Symbol q}} (pb)"'
        do i=1,ndiv_cost
          write(jh,*)xcost4(i),fxcost4tot(i)
        end do
        write(jh,*)'END'
        
        write(jh,*)'DISTRIBUTION:'
        write(jh,*)'DSTRING: "cos3"'
        write(jh,*)'TITLE: "cos{/Symbol q}_{',f,'bar} Distribution"'
        write(jh,*)'AXES: "cos{/Symbol q}"',
     &' "d{/Symbol s}/dcos{/Symbol q}} (pb)"'
        do i=1,ndiv_cost
          write(jh,*)xcost3col(i),fxcost3coltot(i)
        end do
        write(jh,*)'END'
        
        write(jh,*)'DISTRIBUTION:'
        write(jh,*)'DSTRING: "cos4"'
        write(jh,*)'TITLE: "cos{/Symbol q}_{',f,'} Distribution"'
        write(jh,*)'AXES: "cos{/Symbol q}"',
     &' "d{/Symbol s}/dcos{/Symbol q}} (pb)"'
        do i=1,ndiv_cost
          write(jh,*)xcost4col(i),fxcost4coltot(i)
        end do
        write(jh,*)'END'
        
      end if   
      
      if(m_Chi.eq.1)then
c plot distribution in Chi.
        do j=1,ndiv_Chi
          do i=1,it
            fxChi(j,i)=fxChi(j,i)*avgi/cnorm(i)/Chiw
          end do
        end do
        do j=1,ndiv_Chi
          fxChitot(j)=0.d0 
          do i=1,it
            fxChitot(j)=fxChitot(j)+fxChi(j,i)
          end do
        end do
        write(jh,*)'DISTRIBUTION:'
        write(jh,*)'DSTRING: "chi"'
        write(jh,*)'TITLE: "{/Symbol c} Distribution"'
        write(jh,*)'AXES: "{/Symbol c}"',
     &' "d{/Symbol s}/d{/Symbol c} (pb)"'
        do i=1,ndiv_Chi
          write(jh,*)xChi(i),fxChitot(i)
        end do
        write(jh,*)'END'
        
      end if
      
      if(m_Chi2.eq.1)then                
c plot distribution in Chi2.
        do j=1,ndiv_Chi2
          do i=1,it
            fxChi2(j,i)=fxChi2(j,i)*avgi/cnorm(i)/Chi2w
          end do
        end do
        do j=1,ndiv_Chi2
          fxChi2tot(j)=0.d0 
          do i=1,it
            fxChi2tot(j)=fxChi2tot(j)+fxChi2(j,i)
          end do
        end do
        write(jh,*)'DISTRIBUTION:'
        write(jh,*)'DSTRING: "chi2"'
        write(jh,*)'TITLE: "{/Symbol c} Distribution"'
        write(jh,*)'AXES: "{/Symbol c}"',
     &' "1/{/Symbol s}d{/Symbol s}/d{/Symbol c} (pb)"'
        do i=1,ndiv_Chi2
          write(jh,*)xChi2(i),fxChi2tot(i)/cross
        end do
        write(jh,*)'END'
        
      end if

      if(m_Y.eq.1)then
c plot distribution in Y.
        do j=1,ndiv_Y
          do i=1,it
            fxY3(j,i)=fxY3(j,i)*avgi/cnorm(i)/Yw
            fxY4(j,i)=fxY4(j,i)*avgi/cnorm(i)/Yw
          end do
        end do
        do j=1,ndiv_Y
          fxY3tot(j)=0.d0 
          fxY4tot(j)=0.d0 
          do i=1,it
            fxY3tot(j)=fxY3tot(j)+fxY3(j,i)
            fxY4tot(j)=fxY4tot(j)+fxY4(j,i)
          end do
        end do
        
        write(jh,*)'DISTRIBUTION:'
        write(jh,*)'DSTRING: "y3cm"'
        write(jh,*)'TITLE: "',f,'bar CM Rapidity Distribution"'
        write(jh,*)'AXES: "y" "d{/Symbol s}/dy (pb)"'
        do i=1,ndiv_Y
          write(jh,*)xY3(i),fxY3tot(i)
        end do
        write(jh,*)'END'
        
        
        write(jh,*)'DISTRIBUTION:'
        write(jh,*)'DSTRING: "y4cm"'
        write(jh,*)'TITLE: "',f,' CM Rapidity Distribution"'
        write(jh,*)'AXES: "y" "d{/Symbol s}/dy (pb)"'
        do i=1,ndiv_Y
          write(jh,*)xY4(i),fxY4tot(i)
        end do
        write(jh,*)'END'
        

      end if
      
      if(m_Ycol.eq.1)then
c plot distribution in Ycol.
        do j=1,ndiv_Ycol
          do i=1,it
            fxYcol3(j,i)=fxYcol3(j,i)*avgi/cnorm(i)/Ycolw
            fxYcol4(j,i)=fxYcol4(j,i)*avgi/cnorm(i)/Ycolw
          end do
        end do
        do j=1,ndiv_Ycol
          fxYcol3tot(j)=0.d0 
          fxYcol4tot(j)=0.d0 
          do i=1,it
            fxYcol3tot(j)=fxYcol3tot(j)+fxYcol3(j,i)
            fxYcol4tot(j)=fxYcol4tot(j)+fxYcol4(j,i)
          end do
        end do
        write(jh,*)'DISTRIBUTION:'
        write(jh,*)'DSTRING: "y3"'
        write(jh,*)'TITLE: "',f,'bar Rapidity Distribution"'
        write(jh,*)'AXES: "y" "d{/Symbol s}/dy (pb)"'
        do i=1,ndiv_Ycol
          write(jh,*)xYcol3(i),fxYcol3tot(i)
        end do
        write(jh,*)'END'
        
        write(jh,*)'DISTRIBUTION:'
        write(jh,*)'DSTRING: "y4"'
        write(jh,*)'TITLE: "',f,' Rapidity Distribution"'
        write(jh,*)'AXES: "y" "d{/Symbol s}/dy (pb)"'
        do i=1,ndiv_Ycol
          write(jh,*)xYcol4(i),fxYcol4tot(i)
        end do
        write(jh,*)'END'
        
      end if

      if(m_DelY.eq.1)then
c plot distribution in DelY.
        do j=1,ndiv_DelY
          do i=1,it
            fxDelY(j,i)=fxDelY(j,i)*avgi/cnorm(i)/DelYw
          end do
        end do
        do j=1,ndiv_DelY
          fxDelYtot(j)=0.d0 
          do i=1,it
            fxDelYtot(j)=fxDelYtot(j)+fxDelY(j,i)
          end do
        end do
        write(jh,*)'DISTRIBUTION:'
        write(jh,*)'DSTRING: "dely"'
        write(jh,*)'TITLE: "{/Symbol D}y Distribution"'
        write(jh,*)'AXES: "{/Symbol D}y"', 
     &' "d{/Symbol s}/d{/Symbol D}y (pb)"'
        do i=1,ndiv_DelY
          write(jh,*)xDelY(i),fxDelYtot(i)
        end do
        write(jh,*)'END'
        
      end if
 
      if(m_Ytt.eq.1)then
c plot distribution in Ytt.
        do j=1,ndiv_Ytt
          do i=1,it
            fxYtt(j,i)=fxYtt(j,i)*avgi/cnorm(i)/Yttw
          end do
        end do
        do j=1,ndiv_Ytt
          fxYtttot(j)=0.d0 
          do i=1,it
            fxYtttot(j)=fxYtttot(j)+fxYtt(j,i)
          end do
        end do
        write(jh,*)'DISTRIBUTION:'
        write(jh,*)'DSTRING: "yff"'
        write(jh,*)'TITLE: "y_{',ff,'} Distribution"'
        write(jh,*)'AXES: "y_{',ff,'}"',
     &' "d{/Symbol s}/dy_{',ff,'} (pb)"'
        do i=1,ndiv_Ytt
          write(jh,*)xYtt(i),fxYtttot(i)
        end do
        write(jh,*)'END'
        
      end if 

      if(m_x1x2.eq.1)then
c plot distribution in x1x2.
        do j=1,ndiv_x1x2
          do i=1,it
            fxx1x2(j,i)=fxx1x2(j,i)*avgi/cnorm(i)/x1x2w
          end do
        end do
        do j=1,ndiv_x1x2
          fxx1x2tot(j)=0.d0 
          do i=1,it
            fxx1x2tot(j)=fxx1x2tot(j)+fxx1x2(j,i)
          end do
        end do
        write(jh,*)'DISTRIBUTION:'
        write(jh,*)'DSTRING: "x1x2"'
        write(jh,*)'TITLE: "x1/x2 Distribution"'
        write(jh,*)'AXES: "x1/x2"',  
     &' "d{/Symbol s}/d(x1/x2) (pb/GeV)"'        
        do i=1,ndiv_x1x2
          write(jh,*)xx1x2(i),fxx1x2tot(i)
        end do
        write(jh,*)'END'
        
      end if  
     
      if((m_sigp.eq.1).and.(m_sigm.eq.1))then
c plot distribution in asymmetries.
        do j =1,18
        ave(j)=1.d0
        enddo

        do 1234 jasy=1,18
        snorm(jasy)=0.d0
        do j=1,ndiv_sigp
          do i=1,it
            fxsigp(jasy,j,i)=fxsigp(jasy,j,i)*avgi/cnorm(i)/sigpw
          end do
        end do
        do j=1,ndiv_sigp
          fxsigptot(jasy,j)=0.d0 
          do i=1,it
            fxsigptot(jasy,j)=fxsigptot(jasy,j)+fxsigp(jasy,j,i)
          end do
        end do
        do j=1,ndiv_sigm
          do i=1,it
            fxsigm(jasy,j,i)=fxsigm(jasy,j,i)*avgi/cnorm(i)/sigmw
          end do
        end do
        do j=1,ndiv_sigm
          fxsigmtot(jasy,j)=0.d0 
          do i=1,it
            fxsigmtot(jasy,j)=fxsigmtot(jasy,j)+fxsigm(jasy,j,i)
          end do
        end do
        write(jh,*)'ASYMMETRY:'
        
        if(jasy.eq.1)then
              write(jh,*)'DSTRING: "ALL"'
              write(jh,*)'TITLE: "A_{LL} Distribution"'
              write(jh,*)'AXES: "M_{',ff,'} (GeV)" "A_{LL}"'
        else if(jasy.eq.2)then
              write(jh,*)'DSTRING: "AL"'
              write(jh,*)'TITLE: "A_{L} Distribution"'
              write(jh,*)'AXES: "M_{',ff,'} (GeV)" "A_{L}"'
        else if(jasy.eq.3)then
              write(jh,*)'DSTRING: "APV"'
              write(jh,*)'TITLE: "A_{PV} Distribution"'
              write(jh,*)'AXES: "M_{',ff,'} (GeV)" "A_{PV}"'
        else if(jasy.eq.4)then
              write(jh,*)'DSTRING: "AFB2"'
              write(jh,*)'TITLE: "A_{FB}(CM) Distribution"'
              write(jh,*)'AXES: "M_{',ff,'} (GeV)" "A_{FB}"'
        else if(jasy.eq.5)then
              write(jh,*)'DSTRING: "AC"'
              write(jh,*)'TITLE: "A_{C}(',yc,') Distribution"'
              write(jh,*)'AXES: "M_{',ff,'} (GeV)" "A_{C}"'
        else if(jasy.eq.6)then
              write(jh,*)'DSTRING: "AF"'
              write(jh,*)'TITLE: "A_{F}(',yc,') Distribution"'
              write(jh,*)'AXES: "M_{',ff,'} (GeV)" "A_{F}"'
        else if(jasy.eq.7)then
              write(jh,*)'DSTRING: "AFB"'
              write(jh,*)'TITLE: "A_{FB} Distribution"'
              write(jh,*)'AXES: "M_{',ff,'} (GeV)" "A_{FB}"'
        else if(jasy.eq.8)then
              write(jh,*)'DSTRING: "ARFB"'
              write(jh,*)'TITLE: "A_{RFB}(',yff,') Distribution"'
              write(jh,*)'AXES: "M_{',ff,'} (GeV)" "A_{RFB}"'
        else if(jasy.eq.9)then
              write(jh,*)'DSTRING: "AOFB"'
              write(jh,*)'TITLE: "A_{OFB}(',ipff,') Distribution"'
              write(jh,*)'AXES: "M_{',ff,'} (GeV)" "A_{OFB}"'
        else if(jasy.eq.10)then
              write(jh,*)'DSTRING: "AFBSTAR"'
              write(jh,*)'TITLE: "A^{*}_{FB} Distribution"'
              write(jh,*)'AXES: "M_{',ff,'} (GeV)" "A^{*}_{FB}"'
        endif
        do j=1,8
          if(jasy.eq.(10+j))then
            write(intstr,'(i1)') j 
            write(jh,*)'DSTRING: "AFBSTAR'//intstr//'"'
            write(jh,*)'TITLE: "A^{*}_{FB}(',yffs(j),') Distribution"'
            write(jh,*)'AXES: "M_{',ff,'} (GeV)" "A^{*}_{FB}"'
          endif
        enddo
        ndiv_sig=(ndiv_sigm+ndiv_sigp)/2
        do i=1,ndiv_sig
          if(fxsigptot(jasy,i)+fxsigmtot(jasy,i).eq.0.d0)then
            write(jh,*)(xsigm(i)+xsigp(i))/2.d0,0.d0,
     &                fxsigptot(jasy,i),
     &                fxsigmtot(jasy,i)           
            snorm(jasy)=snorm(jasy)+0.d0
          else  
            write(jh,*)(xsigm(i)+xsigp(i))/2.d0,
     &               (fxsigptot(jasy,i)-fxsigmtot(jasy,i))/
     &               (fxsigptot(jasy,i)+fxsigmtot(jasy,i))/ave(jasy),
     &                fxsigptot(jasy,i),
     &                fxsigmtot(jasy,i)
            snorm(jasy)=snorm(jasy)+
     &               (fxsigptot(jasy,i)-fxsigmtot(jasy,i))/
     &               (fxsigptot(jasy,i)+fxsigmtot(jasy,i))/ave(jasy)
     &               *fxrmasstot(i)*rmassw/avgi
          end if
        end do
        write(jh,*)'INTEGRATED:',snorm(jasy)
        write(jh,*)'END'
        
 1234   continue
      end if

c write all out.
      write(jh,*)'(using ',npoints,' points)'
      write(jh,*)'                                       '
      write(jh,*)'sigma (pb):  error (same units):'
      write(jh,*)cross,error 
      write(jh,*)'ALL:                  error (same units):'
      write(jh,*)ALL,ALLerr 
      write(jh,*)'AL:                   error (same units):'
      write(jh,*)AL,ALerr 
      write(jh,*)'APV:                  error (same units):'
      write(jh,*)APV,APVerr 
      write(jh,*)'AFB:                  error (same units):'
      write(jh,*)AFB,AFBerr 
      write(jh,*)'AFB2:                 error (same units):'
      write(jh,*)AFB2,AFB2err 
      write(jh,*)'AC(',yc,'):           error (same units):'
      write(jh,*)AC,ACerr
      write(jh,*)'AF(',yc,'):           error (same units):'
      write(jh,*)AF,AFerr
      write(jh,*)'ARFB(',yff,'):        error (same units):'
      write(jh,*)ARFB,ARFBerr
      write(jh,*)'AOFB(',pff,'):        error (same units):'
      write(jh,*)AOFB,AOFBerr
      write(jh,*)'AFBSTAR:              error (same units):'
      write(jh,*)AFBST,AFBSTerr
      do j=1,8
        write(jh,*)'AFBSTAR(',yffs(j),'):    error (same units):'
        write(jh,*)AFBSTCUT(j),AFBSTCUTerr(j)
      enddo
cc testing:
      stop
      end
c
c     --------------------------------------------------------------------
c
      
      function fxn(x,wgt)
      implicit real*8 (a-h,o-z)
      common/bveg1/ncall,itmx,nprn,ndev,xl(100),xu(100),acc
      common/bveg2/it,ndo,si,swgt,schi,xi(50,100)
      common/rndm/iseed
      common/reslocal/resl(20),standdevl(20)
      common/par/s
      common/limfac/fac
      common/EW/a_em,s2w,zc_norm
      common/W/rmW,gamW  
      common/Z/rmZ,gamZ  
      common/b/rmb
      common/t/rmt,gamt
      common/stat/npoints,nctpoints
      common/coll/ecm_coll
      common/ptr/y,ptCut,dycut
      COMMON/PARTDIST/ISTRUCTURE
      COMMON/ALFASTRONG/rlambdaQCD4,nloops
      COMMON/PDF/QQ
      common/collider/icoll
      common/errors/ierr
c test dummies
      common/test/acbinp(20),acbinm(20),antiacbinp(20),antiacbinm(20)
c zprime model params
      common/pars/param(10),nzp
      common/ZZ/rmZZ(10),gamZZ(10),gp(10)
      common/angle/theta
      common/VA/ggvu(10),ggau(10),ggvd(10),ggad(10),
     &          ggvt(10),ggat(10),ggvb(10),ggab(10),
     &          ggve(10),ggae(10),ggvn(10),ggan(10),
     &          ggvta(10),ggata(10),ggvnt(10),ggant(10)
      common/LR/gZZu(2,10),gZZd(2,10),gZZe(2,10),gZZn(2,10),
     &          gZZc(2,10),gZZs(2,10),gZZt(2,10),gZZb(2,10),
     &          gZZta(2,10),gZZnt(2,10)
c Possible corrections to Z couplings
      common/ZCORR/DgZu(2),DgZd(2),DgZe(2),DgZn(2),
     &             DgZt(2),DgZb(2),DgZta(2),DgZnt(2) 
c few extra variables for looping over M_s, calculating CL, invariant mass cutt
      common/lum/alumpb
      common/cut/ainvmasscut
      common/msloop2/ndist
c contribution and QCD switches
      INTEGER CONT,QCD,EW,BSM
      common/switch/CONT,QCD,EW,BSM
c final state switches, parameters
      CHARACTER*2 ff
      CHARACTER*1 f
      common/fstr/ff,f
      common/fermion/rmf
      common/fints/jf,mf
c distributions.
      common/ext_pT/pTmax,pTmin,pTw
      common/dist_pT/xpT(500),fxpT(500,20),fxpTtot(500)
      common/inp_pT/m_pT
      common/div_pT/ndiv_pT
      common/ext_rmass/rmassmax,rmassmin,rmassw
      common/dist_rmass/xrmass(500),fxrmass(500,20),fxrmasstot(500)
      common/inp_rmass/m_rmass
      common/div_rmass/ndiv_rmass
      common/ext_eta/etamax,etamin,etaw
      common/dist_eta3/xeta3(500),fxeta3(500,20),fxeta3tot(500)
      common/dist_eta4/xeta4(500),fxeta4(500,20),fxeta4tot(500)
      common/dist_eta3col/xeta3col(500),fxeta3col(500,20)
     &       ,fxeta3coltot(500)
      common/dist_eta4col/xeta4col(500),fxeta4col(500,20)
     &       ,fxeta4coltot(500)
      common/inp_eta/m_eta
      common/div_eta/ndiv_eta
      common/ext_beta/betamax,betamin,betaw
      common/dist_beta/xbeta(500),fxbeta(500,20),fxbetatot(500)
      common/inp_beta/m_beta
      common/div_beta/ndiv_beta
      common/ext_cost/costmax,costmin,costw
      common/dist_cost/xcost3(500),fxcost3(500,20),fxcost3tot(500)
      common/dist_cost/xcost4(500),fxcost4(500,20),fxcost4tot(500)
      common/dist_cost/xcost3col(500),fxcost3col(500,20)
     &       ,fxcost3coltot(500)
      common/dist_cost/xcost4col(500),fxcost4col(500,20)
     &       ,fxcost4coltot(500)
      common/inp_cost/m_cost
      common/div_cost/ndiv_cost
      common/ext_Et/Etmax,Etmin,Etw
      common/dist_Et/xEt(500),fxEt(500,20),fxEttot(500)
      common/inp_Et/m_Et
      common/div_Et/ndiv_Et
cc Chi distributions      
      common/ext_Chi/Chimax,Chimin,Chiw
      common/dist_Chi/xChi(500),fxChi(500,20),fxChitot(500)
      common/inp_Chi/m_Chi
      common/div_Chi/ndiv_Chi   
      common/ext_Chi2/Chi2max,Chi2min,Chi2w
      common/dist_Chi2/xChi2(500),fxChi2(500,20),fxChi2tot(500)
      common/inp_Chi2/m_Chi2
      common/div_Chi2/ndiv_Chi2 
cc Rapidity distributions
      common/ext_Y/Ymax,Ymin,Yw
      common/dist_Y3/xY3(500),fxY3(500,20),fxY3tot(500)
      common/dist_Y4/xY4(500),fxY4(500,20),fxY4tot(500)
      common/inp_Y/m_Y
      common/div_Y/ndiv_Y 
      common/ext_Ycol/Ycolmax,Ycolmin,Ycolw
      common/dist_Ycol3/xYcol3(500),fxYcol3(500,20),fxYcol3tot(500)
      common/dist_Ycol4/xYcol4(500),fxYcol4(500,20),fxYcol4tot(500)
      common/inp_Ycol/m_Ycol
      common/div_Ycol/ndiv_Ycol 
      common/ext_DelY/DelYmax,DelYmin,DelYw
      common/dist_DelY/xDelY(500),fxDelY(500,20),fxDelYtot(500)
      common/inp_DelY/m_DelY
      common/div_DelY/ndiv_DelY 
      common/ext_Ytt/Yttmax,Yttmin,Yttw
      common/dist_Ytt/xYtt(500),fxYtt(500,20),fxYtttot(500)
      common/inp_Ytt/m_Ytt
      common/div_Ytt/ndiv_Ytt  
      common/ext_Pztt/Pzttmax,Pzttmin,Pzttw
      common/dist_Pztt/xPztt(500),fxPztt(500,20),fxPztttot(500)
      common/inp_Pztt/m_Pztt
      common/div_Pztt/ndiv_Pztt 
cc x1/x2 distribution
      common/ext_x1x2/x1x2max,x1x2min,x1x2w
      common/dist_x1x2/xx1x2(500),fxx1x2(500,20),fxx1x2tot(500)
      common/inp_x1x2/m_x1x2
      common/div_x1x2/ndiv_x1x2 
cc Asymmetries 
      common/ext_sigp/sigpmax,sigpmin,sigpw
      common/dist_sigp/xsigp(1000),fxsigp(20,1000,20),fxsigptot(20,1000)
      common/inp_sigp/m_sigp
      common/div_sigp/ndiv_sigp
      common/ext_sigm/sigmmax,sigmmin,sigmw
      common/dist_sigm/xsigm(1000),fxsigm(20,1000,20),fxsigmtot(20,1000)
      common/inp_sigm/m_sigm
      common/div_sigm/ndiv_sigm
      common/nfnb/fb(1000,-1:1,5)
      common/polarised/polcross(20,-1:1,-1:1),polerror(20,-1:1,-1:1)
      common/FB/asycross(20,-1:1),asyerror(20,-1:1)
      common/FB2/asy2cross(20,-1:1),asy2error(20,-1:1)
      common/ACENT/accross(20,-1:1),acerror(20,-1:1)
      common/AFORW/afcross(20,-1:1),aferror(20,-1:1)
      common/ARFB/arfbcross(20,-1:1),arfberror(20,-1:1)
      common/AOFB/aofbcross(20,-1:1),aofberror(20,-1:1)
      common/AFBST/afbstcross(20,-1:1),afbsterror(20,-1:1)
      common/AFBST1/afbst1cross(20,-1:1,8),afbst1error(20,-1:1,8)
      common/ACENTcut/yc,yff,pff
      common/AFBSTcut/yffs(8)
      dimension x(100)
      dimension q(4,4),qcol(4,4)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3)
      dimension resgg(-1:1,-1:1),resqq(-1:1,-1:1),resqqalt(-1:1,-1:1)       
      dimension resdd(-1:1,-1:1),resuu(-1:1,-1:1),resss(-1:1,-1:1),
     &          rescc(-1:1,-1:1),resbb(-1:1,-1:1)      
      dimension resddalt(-1:1,-1:1),resuualt(-1:1,-1:1),
     &          resssalt(-1:1,-1:1),resccalt(-1:1,-1:1),
     &          resbbalt(-1:1,-1:1)
      dimension fx1(13),fx2(13)
      dimension pfx(-1:1,-1:1)
      dimension pfxalt(-1:1,-1:1)
      dimension weight(20,-1:1,-1:1) 
      parameter (pi=3.14159265358979323846d0)
c internal random number seed.
      data jseed/987654321/
c initialise.
      do ii=-1,1,2
        do jj=-1,1,2
          resgg(ii,jj)=0.d0
          resqq(ii,jj)=0.d0
          resqqalt(ii,jj)=0.d0
          resdd(ii,jj)=0.d0
          resddalt(ii,jj)=0.d0
          resuu(ii,jj)=0.d0
          resuualt(ii,jj)=0.d0
          resss(ii,jj)=0.d0
          resssalt(ii,jj)=0.d0
          rescc(ii,jj)=0.d0
          resccalt(ii,jj)=0.d0
          resbb(ii,jj)=0.d0
          resbbalt(ii,jj)=0.d0
          do kk=1,20
            weight(kk,ii,jj)=0.d0
          end do
        end do
      end do
      fxn=0.d0
      fxn1=0.d0
      fxn2=0.d0
ccccccccccccccccccccccccccccccccccccccccccccc
c either total energy available OR cutoff M_S
c rmassmax is either ecm_coll or ams depending on ncutmtt switch (0 or 1)
      ecm_max=rmassmax     
ccccccccccccccccccccccccccccccccccccccccccccc

c final state masses are non-zero (y=/=eta).
      rm3=rmf
      rm4=rmf
c allow for implementation of invariant mass cut
      if (rmassmin.eq.0.d0)then
      ecm_min = rm3+rm4
        else 
        ecm_min = ainvmasscut 
        end if
      ecm=x(2)*(ecm_max-ecm_min) + ecm_min
      qcm2=((ecm*ecm-rm3*rm3-rm4*rm4)**2-(2.d0*rm3*rm4)**2)/
     &    (4.d0*ecm*ecm)
      if (qcm2.lt.0.d0) then
        fxn=0.d0
        return
      else
        qcm=sqrt(qcm2)
      endif
c x1 and x2 of the partons.
      shat=ecm*ecm
      tau=shat/s
      x1=x(1)*(1.d0-tau)+tau
      x2=tau/x1
c recontruct kinematics (event is planar and symmetric around z-axis).
      phi=2.d0*pi*ran(jseed)
c symmetry around ct.
      ct = x(3)
      st=sqrt(1.d0-ct*ct)
      q(4,3)=sqrt(qcm**2+rmf**2)
      q(3,3)=qcm*ct
      q(2,3)=qcm*st*cos(phi)
      q(1,3)=qcm*st*sin(phi)
      q(4,4)=+q(4,3)
      q(3,4)=-q(3,3)
      q(2,4)=-q(2,3)
      q(1,4)=-q(1,3)
c this is the scale for the PDFs.
      if(QQ.eq.0.d0)then
          QQ=2.d0*rmf
      else if((QQ.lt.rmZ).and.(QQ.gt.0.d0))then 
          QQ=rmZ
      else if(QQ.eq.-1d0)then ! Dynamical scale set by QQ=-1.0
          QQ = ecm
      endif
c construct polarised hadronic structure functions.
      IF(ISTRUCTURE.LE.4)THEN
        Q2=QQ*QQ
        IF((X1.LE.1.D-6).OR.(X1.GE.1.D0))THEN
          FXN=0.D0
          RETURN
        END IF
        IF((X2.LE.1.D-6).OR.(X2.GE.1.D0))THEN
          FXN=0.D0
          RETURN
        END IF
        IF((QQ.LE.1.3D0).OR.(QQ.GE.1.D4))THEN
          FXN=0.D0
          RETURN
        END IF
        U1=X1*Ctq6Pdf(1,X1,QQ)
        D1=X1*Ctq6Pdf(2,X1,QQ)
        USEA1=X1*Ctq6Pdf(-1,X1,QQ)
        DSEA1=X1*Ctq6Pdf(-2,X1,QQ)
        STR1=X1*Ctq6Pdf(3,X1,QQ)
        CHM1=X1*Ctq6Pdf(4,X1,QQ)
c        BTM1=0.d0
        BTM1=X1*Ctq6Pdf(5,X1,QQ)
        G1=X1*Ctq6Pdf(0,X1,QQ)
        U2=X2*Ctq6Pdf(1,X2,QQ)
        D2=X2*Ctq6Pdf(2,X2,QQ)
        USEA2=X2*Ctq6Pdf(-1,X2,QQ)
        DSEA2=X2*Ctq6Pdf(-2,X2,QQ)
        STR2=X2*Ctq6Pdf(3,X2,QQ)
        CHM2=X2*Ctq6Pdf(4,X2,QQ)
c        BTM2=0.d0
        BTM2=X2*Ctq6Pdf(5,X2,QQ)
        G2=X2*Ctq6Pdf(0,X2,QQ)
      END IF
c actual PDFs.
      IF(ISTRUCTURE.LE.4)THEN 
        fx1(1)=D1
        fx1(2)=U1
c        fx1(2)=1.d0
        fx1(3)=STR1
        fx1(4)=CHM1
        fx1(5)=BTM1
        fx1(6)=0.D0
        fx1(7)=DSEA1
c        fx1(8)=1.d0
        fx1(8)=USEA1
        fx1(9)=fx1(3)
        fx1(10)=fx1(4)
        fx1(11)=fx1(5)
        fx1(12)=fx1(6)
        fx1(13)=G1
        do i=1,13
          fx1(i)=fx1(i)/x1
C            fx1(i)=0.d0
C            if((i.eq.2).or.(i.eq.8)) fx1(i)=1.d0
        end do
        fx2(1)=D2*(1-icoll)+DSEA2*icoll
        fx2(2)=U2*(1-icoll)+USEA2*icoll
c        fx2(2)=1.d0
        fx2(3)=STR2
        fx2(4)=CHM2
        fx2(5)=BTM2
        fx2(6)=0.D0
        fx2(7)=D2*icoll+DSEA2*(1-icoll)
        fx2(8)=U2*icoll+USEA2*(1-icoll)
c        fx2(8)=1.d0
        fx2(9)=fx2(3)
        fx2(10)=fx2(4)
        fx2(11)=fx2(5)
        fx2(12)=fx2(6)
        fx2(13)=G2
        do i=1,13
           fx2(i)=fx2(i)/x2
C            fx2(i)=0.d0
C            if((i.eq.2).or.(i.eq.8)) fx2(i)=1.d0
        end do
      END IF
c incoming momenta (massless).
      pcm=ecm/2.d0
      q(4,1)=pcm
      q(3,1)=pcm
      q(2,1)=0.d0
      q(1,1)=0.d0
      q(4,2)=pcm
      q(3,2)=-pcm
      q(2,2)=0.d0
      q(1,2)=0.d0
c cuts in pseudorapidity.
      rps3=q(3,3)/sqrt(q(1,3)**2+q(2,3)**2+q(3,3)**2)
      if(rps3.lt.-1.d0)rps3=-1.d0
      if(rps3.gt.+1.d0)rps3=+1.d0
      rpl3=dacos(rps3)
      arg3=tan(rpl3/2.d0)
      if(arg3.le.0.d0)arg3=1.d-9
      eta3=-log(arg3)
      eta3alt=-eta3
      rps4=q(3,4)/sqrt(q(1,4)**2+q(2,4)**2+q(3,4)**2)
      if(rps4.lt.-1.d0)rps4=-1.d0
      if(rps4.gt.+1.d0)rps4=+1.d0
      rpl4=dacos(rps4)
      arg4=tan(rpl4/2.d0)
      if(arg4.le.0.d0)arg4=1.d-9
      eta4=-log(arg4)
      eta4alt=-eta4
c      if(abs(eta).gt.y)then
c        fxn=0.d0
c        return
c      end if
      y3 = 0.5d0*dlog((q(4,3)+q(3,3))/(q(4,3)-q(3,3))) 
      y3alt = -y3
      y4 = 0.5d0*dlog((q(4,4)+q(3,4))/(q(4,4)-q(3,4))) 
      y4alt = -y4 
c MadGraph momenta.
      do i=1,3
        p1(i)=q(i,1)
        p2(i)=q(i,2)
        p3(i)=q(i,3)
        p4(i)=q(i,4)
      end do
      p1(0)=q(4,1)
      p2(0)=q(4,2)
      p3(0)=q(4,3)
      p4(0)=q(4,4)
c MANDELSTAMS
      that =(p1(0)-p3(0))**2.d0
      uhat =(p1(0)-p4(0))**2.d0
      do i=1,3
      that = that - (p1(i)-p3(i))**2.d0   
      uhat = uhat - (p1(i)-p4(i))**2.d0   
      end do
ccccccccccccccccccccccc
      goto 777
            write(*,*) 'root s = ',ecm_coll
            write(*,*) 'root shat = ',ecm
        write(*,*) 'p1  =',p1
        write(*,*) 'p2  =',p2
        write(*,*) 'p3  =',p3
        write(*,*) 'p4  =',p4
        delta_E=p1(0)+p2(0)
     &         -p3(0)-p4(0)
        delta_x=p1(1)+p2(1)
     &         -p3(1)-p4(1)
        delta_y=p1(2)+p2(2)
     &         -p3(2)-p4(2)
        delta_z=p1(3)+p2(3)
     &         -p3(3)-p4(3)
        write(*,*) 'delta_E=',delta_E
        write(*,*) 'delta_x=',delta_x
        write(*,*) 'delta_y=',delta_y
        write(*,*) 'delta_z=',delta_z
        rmassa1=sqrt(abs(p1(0)**2-p1(1)**2-p1(2)**2-p1(3)**2))
        rmassa2=sqrt(abs(p2(0)**2-p2(1)**2-p2(2)**2-p2(3)**2))
        rmassa3=sqrt(abs(p3(0)**2-p3(1)**2-p3(2)**2-p3(3)**2))
        rmassa4=sqrt(abs(p4(0)**2-p4(1)**2-p4(2)**2-p4(3)**2))
        write(*,*) 'rm_1 =',rmassa1
        write(*,*) 'rm_2 =',rmassa2
        write(*,*) 'rm_3 =',rmassa3
        write(*,*) 'rm_4 =',rmassa4
 777    continue
ccccccccccccccccccccccc
c initial and final state momenta in the collider CM.
      vcol=(x1-x2)/(x1+x2)
      gcol=(x1+x2)/2.d0/sqrt(x1*x2)
      do i=1,4
        qcol(4,i)=gcol*(q(4,i)+vcol*q(3,i))
        qcol(3,i)=gcol*(q(3,i)+vcol*q(4,i))
        qcol(2,i)=q(2,i)
        qcol(1,i)=q(1,i)
      end do
      if(qcol(3,3).gt.qcol(4,3)) print *, "Pz > E for particle 3"
      if(qcol(3,4).gt.qcol(4,4)) print *, "Pz > E for particle 4"

c collider frame pseudorapidity 
c particle 3
      rpscol3=qcol(3,3)/sqrt(qcol(1,3)**2+qcol(2,3)**2+qcol(3,3)**2)
      if(rpscol3.lt.-1.d0)rpscol3=-1.d0
      if(rpscol3.gt.+1.d0)rpscol3=+1.d0
      rplcol3=dacos(rpscol3)
      argcol3=tan(rplcol3/2.d0)
      if(argcol3.le.0.d0)argcol3=1.d-9
      eta3col=-log(argcol3)
      ycol3 = 0.5d0*dlog((qcol(4,3)+qcol(3,3))/(qcol(4,3)-qcol(3,3))) 
c particle 4
      rpscol4=qcol(3,4)/sqrt(qcol(1,4)**2+qcol(2,4)**2+qcol(3,4)**2)
      if(rpscol4.lt.-1.d0)rpscol4=-1.d0
      if(rpscol4.gt.+1.d0)rpscol4=+1.d0
      rplcol4=dacos(rpscol4)
      argcol4=tan(rplcol4/2.d0)
      if(argcol4.le.0.d0)argcol4=1.d-9
      eta4col=-log(argcol4)
      ycol4 = 0.5d0*dlog((qcol(4,4)+qcol(3,4))/(qcol(4,4)-qcol(3,4)))
c eta cuts
      if(dabs(ycol3).gt.y.OR.dabs(ycol4).gt.y)then
        fxn=0.d0
        return
      end if
c deta cuts
      if (rmf.gt.0.d0) then
      if(dabs(ycol3-ycol4).gt.dycut)then
        fxn=0.d0
        return
      end if
      else
      if(dabs(eta3col-eta4col).gt.dycut)then
        fxn=0.d0
        return
      end if
      end if
c pT cut
      pT=sqrt(qcol(1,3)**2+qcol(2,3)**2)
      if(pT.lt.ptCut)then
        fxn=0.d0
        return
      end if
c CM boost
      ycm = 0.5d0*dlog(x1/x2)
      pzcm= qcol(3,3)+qcol(3,4)
c calculates alphas.
      a_s=alfas(QQ,rlambdaQCD4,nloops)      
      gs2=4.d0*pi*a_s
      gs=sqrt(gs2)
cccccccccccccccccccc
       CALL ME(CONT,QCD,EW,BSM,gs,rmf,p1,p2,p3,p4,resgg,
     &      resqq,resqqalt,resuu,resuualt,resdd,resddalt,
     &      rescc,resccalt,resss,resssalt,resbb,resbbalt)
ccccccccccccccccccccccc
c initial luminosity for total unpolarised cross section.
      pfxtot=0.d0
      pfxalttot=0.d0
      do khel=-1,1,2
        do lhel=-1,1,2
          pfx(khel,lhel)=
     &  resgg(khel,lhel)*fx1(13)*fx2(13)/2.d0
     &+(resqq(khel,lhel)+resdd(khel,lhel))*fx1( 1)*fx2( 7)
     &+(resqq(khel,lhel)+resuu(khel,lhel))*fx1( 2)*fx2( 8)
     &+(resqq(khel,lhel)+resss(khel,lhel))*fx1( 3)*fx2( 9)
     &+(resqq(khel,lhel)+rescc(khel,lhel))*fx1( 4)*fx2(10)
     &+(resbb(khel,lhel))*fx1( 5)*fx2(11)
          pfxalt(khel,lhel)= 
     &  resgg(khel,lhel)*fx1(13)*fx2(13)/2.d0
     &+(resqqalt(khel,lhel)+resddalt(khel,lhel))*fx1( 7)*fx2( 1)
     &+(resqqalt(khel,lhel)+resuualt(khel,lhel))*fx1( 8)*fx2( 2)
     &+(resqqalt(khel,lhel)+resssalt(khel,lhel))*fx1( 9)*fx2( 3)
     &+(resqqalt(khel,lhel)+resccalt(khel,lhel))*fx1(10)*fx2( 4)
     &+(resbbalt(khel,lhel))*fx1(11)*fx2( 5)
          pfxtot=pfxtot
     &          +pfx(khel,lhel)
          pfxalttot=pfxalttot
     &          +pfxalt(khel,lhel)
        end do
      end do
      if(pfxtot.eq.0.d0.and.pfxalttot.eq.0.d0)then
        fxn=0.d0
        return
      end if
c for distributions, 
      do khel=-1,1,2
        do lhel=-1,1,2
          pfx(khel,lhel)=pfx(khel,lhel)/(pfxtot+pfxalttot)
          pfxalt(khel,lhel)=pfxalt(khel,lhel)/(pfxtot+pfxalttot)
        end do
      end do      
c Jacobians from dx1 dx2 -> dx(1) dEcm.
      pfxtot=pfxtot*(1.d0-tau)*2.d0*ecm/s/x1
     &      *(ecm_max-ecm_min)
      pfxalttot=pfxalttot*(1.d0-tau)*2.d0*ecm/s/x1
     &      *(ecm_max-ecm_min)
c MEs and PDFs.
      fxn1=pfxtot
      fxn2=pfxalttot
c constant factor.
      fxn1=fxn1*fac
      fxn2=fxn2*fac
c jacobians by hand, pi's and flux (1/2s)
cccccccccccccccccccccc
      fxn1=fxn1*qcm/(2.d0*pcm)*2.d0**(4-3*(2))
      fxn1=fxn1/2.d0/ecm/ecm*(2.d0*pi)**(4-3*(2))
      fxn2=fxn2*qcm/(2.d0*pcm)*2.d0**(4-3*(2))
      fxn2=fxn2/2.d0/ecm/ecm*(2.d0*pi)**(4-3*(2))
      fxn=fxn1+fxn2
c      print*,"ecm,ct,x1 = ",ecm,ct,x1
cccccccccccccccccccccc
c polarised cross sections.
      do iphel=-1,+1,2
        do jphel=-1,+1,2
          polcross(it,iphel,jphel)=polcross(it,iphel,jphel)
     &                            +fxn
     &                            *wgt          
     &                            *(pfx(iphel,jphel)         
     &                            +pfxalt(iphel,jphel))
          weight(it,iphel,jphel)=+fxn
     &                           *wgt
     &                           *(pfx(iphel,jphel)
     &                           +pfxalt(iphel,jphel))
          polerror(it,iphel,jphel)=polerror(it,iphel,jphel)
     &                            +polcross(it,iphel,jphel)**2
        end do
      end do
c FB asymmetry. IN CM FRAME!
        cost4=(q(1,4)*q(1,1)
     &       +q(2,4)*q(2,1)
     &       +q(3,4)*q(3,1))
     &  /sqrt(q(1,4)*q(1,4)
     &       +q(2,4)*q(2,4)
     &       +q(3,4)*q(3,4))
     &  /sqrt(q(1,1)*q(1,1)
     &       +q(2,1)*q(2,1)
     &       +q(3,1)*q(3,1))
        cost4alt=-cost4
      if(cost4.eq.0.d0)then
        continue
      else if(cost4.gt.0.d0)then
        asy2cross(it,+1)=asy2cross(it,+1)+fxn1*wgt 
        asy2error(it,+1)=asy2error(it,+1)
     &                 +asy2cross(it,+1)**2
        asy2cross(it,-1)=asy2cross(it,-1)+fxn2*wgt
        asy2error(it,-1)=asy2error(it,-1)
     &                 +asy2cross(it,-1)**2
      else if(cost4.lt.0.d0)then
        asy2cross(it,-1)=asy2cross(it,-1)+fxn1*wgt
        asy2error(it,-1)=asy2error(it,-1)
     &                 +asy2cross(it,-1)**2
        asy2cross(it,+1)=asy2cross(it,+1)+fxn2*wgt
        asy2error(it,+1)=asy2error(it,+1)
     &                 +asy2cross(it,-1)**2
      end if

c FB asymmetry. IN Collider FRAME!
        cost4col=(qcol(1,4)*qcol(1,1)
     &       +qcol(2,4)*qcol(2,1)
     &       +qcol(3,4)*qcol(3,1))
     &  /sqrt(qcol(1,4)*qcol(1,4)
     &       +qcol(2,4)*qcol(2,4)
     &       +qcol(3,4)*qcol(3,4))
     &  /sqrt(qcol(1,1)*qcol(1,1)
     &       +qcol(2,1)*qcol(2,1)
     &       +qcol(3,1)*qcol(3,1))
      if(cost4col.eq.0.d0)then
        continue
      else if(cost4col.gt.0.d0)then
        asycross(it,+1)=asycross(it,+1)+fxn*wgt
        asyerror(it,+1)=asyerror(it,+1)
     &                 +asycross(it,+1)**2
      else if(cost4col.lt.0.d0)then
        asycross(it,-1)=asycross(it,-1)+fxn*wgt
        asyerror(it,-1)=asyerror(it,-1)
     &                 +asycross(it,-1)**2
      end if 
c Central Asymmetry

      if(dabs(ycol4).le.yc)then
        accross(it,+1)=accross(it,+1)+fxn*wgt
        acerror(it,+1)=acerror(it,+1)
     &                 +accross(it,+1)**2
      end if
      if(dabs(ycol3).le.yc)then
        accross(it,-1)=accross(it,-1)+fxn*wgt
        acerror(it,-1)=acerror(it,-1)
     &                 +accross(it,-1)**2
      end if
c Forward Asymmetry
      if(dabs(ycol4).gt.yc)then
        afcross(it,+1)=afcross(it,+1)+fxn*wgt
        aferror(it,+1)=aferror(it,+1)
     &                 +afcross(it,+1)**2
      end if
      if(dabs(ycol3).gt.yc)then
        afcross(it,-1)=afcross(it,-1)+fxn*wgt
        aferror(it,-1)=aferror(it,-1)
     &                 +afcross(it,-1)**2
      end if
c ARFB/AOFB (rapidity dependent and one sided AFB)
      if(dabs(ycol4).gt.dabs(ycol3))then
        if(abs(ycm).gt.yff)then
        arfbcross(it,+1)=arfbcross(it,+1)+fxn*wgt
        arfberror(it,+1)=arfberror(it,+1)
     &                 +arfbcross(it,+1)**2
        end if
        if(abs(pzcm).gt.pff)then
        aofbcross(it,+1)=aofbcross(it,+1)+fxn*wgt
        aofberror(it,+1)=aofberror(it,+1)
     &                 +aofbcross(it,+1)**2
        end if
      end if
      if(dabs(ycol3).gt.dabs(ycol4))then
        if(abs(ycm).gt.yff)then
        arfbcross(it,-1)=arfbcross(it,-1)+fxn*wgt
        arfberror(it,-1)=arfberror(it,-1)
     &                 +arfbcross(it,-1)**2
        end if
        if(abs(pzcm).gt.pff)then
        aofbcross(it,-1)=aofbcross(it,-1)+fxn*wgt
        aofberror(it,-1)=aofberror(it,-1)
     &                 +aofbcross(it,-1)**2
        end if      
      end if
c AFBSTAR
         costst=int(ycm/dabs(ycm))*cost4
      if(costst.gt.0.d0)then
         afbstcross(it,+1)=afbstcross(it,+1)+fxn*wgt
         afbsterror(it,+1)=afbsterror(it,+1)
     &                 +afbstcross(it,+1)**2
          do j=1,8
            if(abs(ycm).gt.yffs(j))then
              afbst1cross(it,+1,j)=afbst1cross(it,+1,j)+fxn*wgt
              afbst1error(it,+1,j)=afbst1error(it,+1,j)
     &                 +afbst1cross(it,+1,j)**2
            end if
          enddo
       else if(costst.lt.0.d0)then
         afbstcross(it,-1)=afbstcross(it,-1)+fxn*wgt
         afbsterror(it,-1)=afbsterror(it,-1)
     &                 +afbstcross(it,-1)**2
          do j=1,8
            if(abs(ycm).gt.yffs(j))then
              afbst1cross(it,-1,j)=afbst1cross(it,-1,j)+fxn*wgt
              afbst1error(it,-1,j)=afbst1error(it,-1,j)
     &                 +afbst1cross(it,-1,j)**2
            end if
          enddo
       end if
c histograms.
      hist1=fxn1*wgt
      hist2=fxn2*wgt
      hist=hist1+hist2
      if(m_pT.eq.1)then
c generate distribution in pT.
        nbin=int((pT-pTmin)/pTw)+1
        if(nbin.ge.(ndiv_pT+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxpT(nbin,it)=fxpT(nbin,it)+hist
        end if
      end if  
      if(m_rmass.eq.1)then
c generate distribution in rmass.
        rmass2=(qcol(4,3)+qcol(4,4))**2
        do i=1,3
          rmass2=rmass2-(qcol(i,3)+qcol(i,4))**2
        end do
        rmass=sqrt(abs(rmass2))
        nbin=int((rmass-rmassmin)/rmassw)+1
        if(nbin.ge.(ndiv_rmass+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxrmass(nbin,it)=fxrmass(nbin,it)+hist
        end if
      end if  
      if(m_eta.eq.1)then
c generate distribution in eta3.
        nbin=int((eta3-etamin)/etaw)+1
        if(nbin.ge.(ndiv_eta+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxeta3(nbin,it)=fxeta3(nbin,it)+hist1
        end if
        nbin=int((eta3alt-etamin)/etaw)+1
        if(nbin.ge.(ndiv_eta+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxeta3(nbin,it)=fxeta3(nbin,it)+hist2
        end if
c Collider eta
        nbin=int((eta3col-etamin)/etaw)+1
        if(nbin.ge.(ndiv_eta+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxeta3col(nbin,it)=fxeta3col(nbin,it)+hist
        end if
      end if  
      if(m_eta.eq.1)then
c generate distribution in eta4.
        nbin=int((eta4-etamin)/etaw)+1
        if(nbin.ge.(ndiv_eta+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxeta4(nbin,it)=fxeta4(nbin,it)+hist1
        end if
        nbin=int((eta4alt-etamin)/etaw)+1
        if(nbin.ge.(ndiv_eta+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxeta4(nbin,it)=fxeta4(nbin,it)+hist2
        end if
c Collider eta
        nbin=int((eta4col-etamin)/etaw)+1
        if(nbin.ge.(ndiv_eta+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxeta4col(nbin,it)=fxeta4col(nbin,it)+hist
        end if
      end if  
      if(m_beta.eq.1)then
c generate distribution in beta.
        beta=shat/4.d0/rmf**2-1.d0
        nbin=int((beta-betamin)/betaw)+1
        if(nbin.ge.(ndiv_beta+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxbeta(nbin,it)=fxbeta(nbin,it)+hist
        end if
      end if  
      if(m_cost.eq.1)then
c generate distribution in cost.
c collider antitop
        cost3col=(qcol(1,3)*qcol(1,1)
     &       +qcol(2,3)*qcol(2,1)
     &       +qcol(3,3)*qcol(3,1))
     &  /sqrt(qcol(1,3)*qcol(1,3)
     &       +qcol(2,3)*qcol(2,3)
     &       +qcol(3,3)*qcol(3,3))
     &  /sqrt(qcol(1,1)*qcol(1,1)
     &       +qcol(2,1)*qcol(2,1)
     &       +qcol(3,1)*qcol(3,1))
        nbin=int((cost3col-costmin)/costw)+1
        if(nbin.ge.(ndiv_cost+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxcost3col(nbin,it)=fxcost3col(nbin,it)+hist
        end if
c collider top
        cost4col=(qcol(1,4)*qcol(1,1)
     &       +qcol(2,4)*qcol(2,1)
     &       +qcol(3,4)*qcol(3,1))
     &  /sqrt(qcol(1,4)*qcol(1,4)
     &       +qcol(2,4)*qcol(2,4)
     &       +qcol(3,4)*qcol(3,4))
     &  /sqrt(qcol(1,1)*qcol(1,1)
     &       +qcol(2,1)*qcol(2,1)
     &       +qcol(3,1)*qcol(3,1))
        nbin=int((cost4col-costmin)/costw)+1
        if(nbin.ge.(ndiv_cost+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxcost4col(nbin,it)=fxcost4col(nbin,it)+hist
        end if
c CM antitop
        cost3=(q(1,3)*q(1,1)
     &       +q(2,3)*q(2,1)
     &       +q(3,3)*q(3,1))
     &  /sqrt(q(1,3)*q(1,3)
     &       +q(2,3)*q(2,3)
     &       +q(3,3)*q(3,3))
     &  /sqrt(q(1,1)*q(1,1)
     &       +q(2,1)*q(2,1)
     &       +q(3,1)*q(3,1))
      cost3alt = -cost3
        nbin=int((cost3-costmin)/costw)+1
        if(nbin.ge.(ndiv_cost+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxcost3(nbin,it)=fxcost3(nbin,it)+hist1
        end if
        nbin=int((cost3alt-costmin)/costw)+1
        if(nbin.ge.(ndiv_cost+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxcost3(nbin,it)=fxcost3(nbin,it)+hist2
        end if
c CM top
        cost4=(q(1,4)*q(1,1)
     &       +q(2,4)*q(2,1)
     &       +q(3,4)*q(3,1))
     &  /sqrt(q(1,4)*q(1,4)
     &       +q(2,4)*q(2,4)
     &       +q(3,4)*q(3,4))
     &  /sqrt(q(1,1)*q(1,1)
     &       +q(2,1)*q(2,1)
     &       +q(3,1)*q(3,1))
        cost4alt=-cost4
        nbin=int((cost4-costmin)/costw)+1
        if(nbin.ge.(ndiv_cost+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxcost4(nbin,it)=fxcost4(nbin,it)+hist1
        end if
        nbin=int((cost4alt-costmin)/costw)+1
        if(nbin.ge.(ndiv_cost+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxcost4(nbin,it)=fxcost4(nbin,it)+hist2
        end if
      end if  
      if(m_Et.eq.1)then
c generate distribution in Et.
        Et=qcol(4,4)
        nbin=int((Et-Etmin)/Etw)+1
        if(nbin.ge.(ndiv_Et+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxEt(nbin,it)=fxEt(nbin,it)+hist
        end if
      end if  
      
      if(m_Chi.eq.1)then
c generate distribution in Chi.
        Chi = dexp(dabs(ycol3-ycol4))
        nbin=int((Chi-Chimin)/Chiw)+1
        if(nbin.ge.(ndiv_Chi+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxChi(nbin,it)=fxChi(nbin,it)+hist
        end if
      end if  
      if(m_Chi2.eq.1)then
c generate distribution in Chi2.
        Chi2 = dexp(dabs(ycol3-ycol4))
        nbin=int((Chi2-Chi2min)/Chi2w)+1
        if(nbin.ge.(ndiv_Chi2+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxChi2(nbin,it)=fxChi2(nbin,it)+hist
        end if
      end if 

      if(m_Y.eq.1)then
c generate distribution in Y3.
        nbin=int((Y3-Ymin)/Yw)+1
        if(nbin.ge.(ndiv_Y+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxY3(nbin,it)=fxY3(nbin,it)+hist1
        end if
        nbin=int((Y3alt-Ymin)/Yw)+1
        if(nbin.ge.(ndiv_Y+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxY3(nbin,it)=fxY3(nbin,it)+hist2
        end if
c generate distribution in Y4.
        nbin=int((Y4-Ymin)/Yw)+1
        if(nbin.ge.(ndiv_Y+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxY4(nbin,it)=fxY4(nbin,it)+hist1
        end if
        nbin=int((Y4alt-Ymin)/Yw)+1
        if(nbin.ge.(ndiv_Y+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxY4(nbin,it)=fxY4(nbin,it)+hist2
        end if
      end if
      if(m_Ycol.eq.1)then
c generate distribution in Ycol3.
        nbin=int((Ycol3-Ycolmin)/Ycolw)+1
        if(nbin.ge.(ndiv_Ycol+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxYcol3(nbin,it)=fxYcol3(nbin,it)+hist
        end if
c generate distribution in Ycol4.
        nbin=int((Ycol4-Ycolmin)/Ycolw)+1
        if(nbin.ge.(ndiv_Ycol+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxYcol4(nbin,it)=fxYcol4(nbin,it)+hist
        end if
      end if
      if(m_DelY.eq.1)then
c generate distribution in DelY.
        DelY = ycol4-ycol3
        nbin=int((DelY-DelYmin)/DelYw)+1
        if(nbin.ge.(ndiv_DelY+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxDelY(nbin,it)=fxDelY(nbin,it)+hist
        end if
      end if  

      if(m_Ytt.eq.1)then
c generate distribution in Ytt.
        nbin=int((ycm-Yttmin)/Yttw)+1
        if(nbin.ge.(ndiv_Ytt+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxYtt(nbin,it)=fxYtt(nbin,it)+hist
        end if
      end if
      x1x2 = X1/X2
      if(m_x1x2.eq.1)then
c generate distribution in x1x2.
        nbin=int((x1x2-x1x2min)/x1x2w)+1
        if(nbin.ge.(ndiv_x1x2+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxx1x2(nbin,it)=fxx1x2(nbin,it)+hist
        end if
      end if

      if(m_Pztt.eq.1)then
c generate distribution in Pztt.
        nbin=int((dabs(pzcm)-Pzttmin)/Pzttw)+1
        if(nbin.ge.(ndiv_Pztt+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxPztt(nbin,it)=fxPztt(nbin,it)+hist
        end if
      end if


      if(m_sigp.eq.1)then
c generate distribution in sigp for ALL.
        sigp=ecm
        nbin=int((sigp-sigpmin)/sigpw)+1
        if(nbin.ge.(ndiv_sigp+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxsigp(1,nbin,it)=fxsigp(1,nbin,it)+
     &    (weight(it,+1,+1)+weight(it,-1,-1))
        end if
      end if  
      if(m_sigm.eq.1)then
c generate distribution in sigm for ALL.
        sigm=ecm
        nbin=int((sigm-sigmmin)/sigmw)+1
        if(nbin.ge.(ndiv_sigm+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxsigm(1,nbin,it)=fxsigm(1,nbin,it)+
     &    (weight(it,+1,-1)+weight(it,-1,+1))
        end if
      end if  
      if(m_sigp.eq.1)then
c generate distribution in sigp for AL.
        sigp=ecm
        nbin=int((sigp-sigpmin)/sigpw)+1
        if(nbin.ge.(ndiv_sigp+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxsigp(2,nbin,it)=fxsigp(2,nbin,it)+
     &    (weight(it,-1,-1)+weight(it,-1,+1))
        end if
      end if  
      if(m_sigm.eq.1)then
c generate distribution in sigm for AL.
        sigm=ecm
        nbin=int((sigm-sigmmin)/sigmw)+1
        if(nbin.ge.(ndiv_sigm+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxsigm(2,nbin,it)=fxsigm(2,nbin,it)+
     &    (weight(it,+1,-1)+weight(it,+1,+1))
        end if
      end if  
      if(m_sigp.eq.1)then
c generate distribution in sigp for APV.
        sigp=ecm
        nbin=int((sigp-sigpmin)/sigpw)+1
        if(nbin.ge.(ndiv_sigp+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxsigp(3,nbin,it)=fxsigp(3,nbin,it)+
     &    (weight(it,-1,-1))
        end if
      end if  
      if(m_sigm.eq.1)then
c generate distribution in sigm for APV.
        sigm=ecm
        nbin=int((sigm-sigmmin)/sigmw)+1
        if(nbin.ge.(ndiv_sigm+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxsigm(3,nbin,it)=fxsigm(3,nbin,it)+
     &    (weight(it,+1,+1))
        end if
      end if  
      if((m_sigp.eq.1).and.(cost4col.gt.0.d0))then
c generate distribution in sigp for AFB.
        sigp=ecm
        nbin=int((sigp-sigpmin)/sigpw)+1
        if(nbin.ge.(ndiv_sigp+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxsigp(4,nbin,it)=fxsigp(4,nbin,it)+hist
        end if
      end if  
      if((m_sigm.eq.1).and.(cost4col.lt.0.d0))then
c generate distribution in sigm for AFB.
        sigm=ecm
        nbin=int((sigm-sigmmin)/sigmw)+1
        if(nbin.ge.(ndiv_sigm+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxsigm(4,nbin,it)=fxsigm(4,nbin,it)+hist
        end if
      end if  
      if((m_sigp.eq.1).and.(cost4.gt.0.d0))then
c generate distribution in sigp for AFB2.
        sigp=ecm
        nbin=int((sigp-sigpmin)/sigpw)+1
        if(nbin.ge.(ndiv_sigp+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxsigp(7,nbin,it)=fxsigp(7,nbin,it)+hist1
          fxsigm(7,nbin,it)=fxsigm(7,nbin,it)+hist2
        end if
      end if  
      if((m_sigm.eq.1).and.(cost4.lt.0.d0))then
c generate distribution in sigm for AFB2.
        sigm=ecm
        nbin=int((sigm-sigmmin)/sigmw)+1
        if(nbin.ge.(ndiv_sigm+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxsigm(7,nbin,it)=fxsigm(7,nbin,it)+hist1
          fxsigp(7,nbin,it)=fxsigp(7,nbin,it)+hist2
        end if
      end if  
      if((m_sigp.eq.1).and.(dabs(ycol4).le.yc))then
c generate distribution in sigp for AC.
        sigp=ecm
        nbin=int((sigp-sigpmin)/sigpw)+1
        if(nbin.ge.(ndiv_sigp+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxsigp(5,nbin,it)=fxsigp(5,nbin,it)+hist
          if(nbin.eq.11)acbinp(it)=acbinp(it)+hist
        end if
      end if
      if((m_sigp.eq.1).and.(dabs(ycol4).gt.yc))then
        sigp=ecm
        nbin=int((sigp-sigpmin)/sigpw)+1
        if(nbin.eq.11)antiacbinp(it)=antiacbinp(it)+hist
      end if
        
      if((m_sigm.eq.1).and.(dabs(ycol3).le.yc))then
c generate distribution in sigm for AC.
        sigm=ecm
        nbin=int((sigm-sigmmin)/sigmw)+1
        if(nbin.ge.(ndiv_sigm+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxsigm(5,nbin,it)=fxsigm(5,nbin,it)+hist
          if(nbin.eq.11)acbinm(it)=acbinm(it)+hist
        end if
      end if
      if((m_sigm.eq.1).and.(dabs(ycol3).gt.yc))then
        sigm=ecm
        nbin=int((sigm-sigmmin)/sigmw)+1
        if(nbin.eq.11)antiacbinm(it)=antiacbinm(it)+hist  
      end if
c generate distribution in sigp for AF.
      if((m_sigp.eq.1).and.(dabs(ycol4).gt.yc))then
        sigp=ecm
        nbin=int((sigp-sigpmin)/sigpw)+1
        if(nbin.ge.(ndiv_sigp+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxsigp(6,nbin,it)=fxsigp(6,nbin,it)+hist
        end if
      end if
        
      if((m_sigm.eq.1).and.(dabs(ycol3).gt.yc))then
c generate distribution in sigm for AF.
        sigm=ecm
        nbin=int((sigm-sigmmin)/sigmw)+1
        if(nbin.ge.(ndiv_sigm+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxsigm(6,nbin,it)=fxsigm(6,nbin,it)+hist
        end if
      end if
c generate distribution in sigp for ARFB/AOFB.
      if((m_sigp.eq.1).and.(dabs(ycol4).gt.dabs(ycol3)))then
        sigp=ecm
        nbin=int((sigp-sigpmin)/sigpw)+1
        if(nbin.ge.(ndiv_sigp+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
            if(dabs(ycm).gt.yff)fxsigp(8,nbin,it)=fxsigp(8,nbin,it)+hist
            if(dabs(pzcm).gt.pff)fxsigp(9,nbin,it)=fxsigp(9,nbin,it)+hist         
        end if
      end if
        
      if((m_sigm.eq.1).and.(dabs(ycol3).gt.dabs(ycol4)))then
c generate distribution in sigm for ARFB/AOFB.
        sigm=ecm
        nbin=int((sigm-sigmmin)/sigmw)+1
        if(nbin.ge.(ndiv_sigm+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
            if(dabs(ycm).gt.yff)fxsigm(8,nbin,it)=fxsigm(8,nbin,it)+hist
            if(dabs(pzcm).gt.pff)fxsigm(9,nbin,it)=fxsigm(9,nbin,it)+hist 
        end if
      end if
c generate distribution in sigp for AFBSTAR.
      if((m_sigp.eq.1).and.(costst.gt.0.d0))then
        sigp=ecm
        nbin=int((sigp-sigpmin)/sigpw)+1
        if(nbin.ge.(ndiv_sigp+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
            fxsigp(10,nbin,it)=fxsigp(10,nbin,it)+hist 
            do j=1,8
              if(dabs(ycm).gt.yffs(j))then
                fxsigp(10+j,nbin,it)=fxsigp(10+j,nbin,it)+hist
              endif
            enddo
        end if
      end if
c generate distribution in sigm for AFBSTAR.
      if((m_sigm.eq.1).and.(costst.lt.0.d0))then
        sigm=ecm
        nbin=int((sigm-sigmmin)/sigmw)+1
        if(nbin.ge.(ndiv_sigm+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
            fxsigm(10,nbin,it)=fxsigm(10,nbin,it)+hist 
            do j=1,8
              if(dabs(ycm).gt.yffs(j))then
                fxsigm(10+j,nbin,it)=fxsigm(10+j,nbin,it)+hist
              endif
            enddo
        end if
      end if
c statistics.
      npoints=npoints+1  
      return
      end
      
c
c     --------------------------------------------------------------
c
c     end program d_pptt_LO

      function fxncosx(x,wgt)
      implicit real*8 (a-h,o-z)
      common/bveg1/ncall,itmx,nprn,ndev,xl(100),xu(100),acc
      common/bveg2/it,ndo,si,swgt,schi,xi(50,100)
      common/rndm/iseed
      common/reslocal/resl(20),standdevl(20)
      common/par/s
      common/limfac/fac
      common/EW/a_em,s2w,zc_norm
      common/W/rmW,gamW  
      common/Z/rmZ,gamZ  
      common/b/rmb
      common/t/rmt,gamt
      common/stat/npoints,nctpoints
      common/coll/ecm_coll
      common/ptr/y,ptCut,dycut
      COMMON/PARTDIST/ISTRUCTURE
      COMMON/PDF/QQ
      COMMON/ALFASTRONG/rlambdaQCD4,nloops
      common/collider/icoll
      common/errors/ierr
c zprime model params
      common/pars/param(10),nzp
      common/ZZ/rmZZ(10),gamZZ(10),gp(10)
      common/angle/theta
      common/VA/ggvu(10),ggau(10),ggvd(10),ggad(10),
     &          ggvt(10),ggat(10),ggvb(10),ggab(10),
     &          ggve(10),ggae(10),ggvn(10),ggan(10),
     &          ggvta(10),ggata(10),ggvnt(10),ggant(10)
      common/LR/gZZu(2,10),gZZd(2,10),gZZe(2,10),gZZn(2,10),
     &          gZZc(2,10),gZZs(2,10),gZZt(2,10),gZZb(2,10),
     &          gZZta(2,10),gZZnt(2,10)
c Possible corrections to Z couplings
      common/ZCORR/DgZu(2),DgZd(2),DgZe(2),DgZn(2),
     &             DgZt(2),DgZb(2),DgZta(2),DgZnt(2) 
c few extra variables for looping over M_s, calculating CL, invariant mass cutt
      common/lum/alumpb
      common/cut/ainvmasscut
      common/msloop2/ndist
c contribution and QCD switches
      INTEGER CONT,QCD,EW,BSM
      common/switch/CONT,QCD,EW,BSM
c final state switches, parameters
      CHARACTER*2 ff
      CHARACTER*1 f
      common/fstr/ff,f
      common/fermion/rmf
      common/fints/jf,mf
c distributions.
      common/ext_pT/pTmax,pTmin,pTw
      common/dist_pT/xpT(500),fxpT(500,20),fxpTtot(500)
      common/inp_pT/m_pT
      common/div_pT/ndiv_pT
      common/ext_rmass/rmassmax,rmassmin,rmassw
      common/dist_rmass/xrmass(500),fxrmass(500,20),fxrmasstot(500)
      common/inp_rmass/m_rmass
      common/div_rmass/ndiv_rmass
      common/ext_eta/etamax,etamin,etaw
      common/dist_eta3/xeta3(500),fxeta3(500,20),fxeta3tot(500)
      common/dist_eta4/xeta4(500),fxeta4(500,20),fxeta4tot(500)
      common/dist_eta3col/xeta3col(500),fxeta3col(500,20)
     &       ,fxeta3coltot(500)
      common/dist_eta4col/xeta4col(500),fxeta4col(500,20)
     &       ,fxeta4coltot(500)
      common/inp_eta/m_eta
      common/div_eta/ndiv_eta
      common/ext_beta/betamax,betamin,betaw
      common/dist_beta/xbeta(500),fxbeta(500,20),fxbetatot(500)
      common/inp_beta/m_beta
      common/div_beta/ndiv_beta
      common/ext_cost/costmax,costmin,costw
      common/dist_cost/xcost3(500),fxcost3(500,20),fxcost3tot(500)
      common/dist_cost/xcost4(500),fxcost4(500,20),fxcost4tot(500)
      common/dist_cost/xcost3col(500),fxcost3col(500,20)
     &       ,fxcost3coltot(500)
      common/dist_cost/xcost4col(500),fxcost4col(500,20)
     &       ,fxcost4coltot(500)
      common/inp_cost/m_cost
      common/div_cost/ndiv_cost
      common/ext_Et/Etmax,Etmin,Etw
      common/dist_Et/xEt(500),fxEt(500,20),fxEttot(500)
      common/inp_Et/m_Et
      common/div_Et/ndiv_Et
cc Chi distributions      
      common/ext_Chi/Chimax,Chimin,Chiw
      common/dist_Chi/xChi(500),fxChi(500,20),fxChitot(500)
      common/inp_Chi/m_Chi
      common/div_Chi/ndiv_Chi   
      common/ext_Chi2/Chi2max,Chi2min,Chi2w
      common/dist_Chi2/xChi2(500),fxChi2(500,20),fxChi2tot(500)
      common/inp_Chi2/m_Chi2
      common/div_Chi2/ndiv_Chi2 
cc Rapidity distributions
      common/ext_Y/Ymax,Ymin,Yw
      common/dist_Y3/xY3(500),fxY3(500,20),fxY3tot(500)
      common/dist_Y4/xY4(500),fxY4(500,20),fxY4tot(500)
      common/inp_Y/m_Y
      common/div_Y/ndiv_Y 
      common/ext_Ycol/Ycolmax,Ycolmin,Ycolw
      common/dist_Ycol3/xYcol3(500),fxYcol3(500,20),fxYcol3tot(500)
      common/dist_Ycol4/xYcol4(500),fxYcol4(500,20),fxYcol4tot(500)
      common/inp_Ycol/m_Ycol
      common/div_Ycol/ndiv_Ycol 
      common/ext_DelY/DelYmax,DelYmin,DelYw
      common/dist_DelY/xDelY(500),fxDelY(500,20),fxDelYtot(500)
      common/inp_DelY/m_DelY
      common/div_DelY/ndiv_DelY 
      common/ext_Ytt/Yttmax,Yttmin,Yttw
      common/dist_Ytt/xYtt(500),fxYtt(500,20),fxYtttot(500)
      common/inp_Ytt/m_Ytt
      common/div_Ytt/ndiv_Ytt  
      common/ext_Pztt/Pzttmax,Pzttmin,Pzttw
      common/dist_Pztt/xPztt(500),fxPztt(500,20),fxPztttot(500)
      common/inp_Pztt/m_Pztt
      common/div_Pztt/ndiv_Pztt 
cc x1/x2 distribution
      common/ext_x1x2/x1x2max,x1x2min,x1x2w
      common/dist_x1x2/xx1x2(500),fxx1x2(500,20),fxx1x2tot(500)
      common/inp_x1x2/m_x1x2
      common/div_x1x2/ndiv_x1x2       
cc Asymmetries 
      common/ext_sigp/sigpmax,sigpmin,sigpw
      common/dist_sigp/xsigp(1000),fxsigp(20,1000,20),fxsigptot(20,1000)
      common/inp_sigp/m_sigp
      common/div_sigp/ndiv_sigp
      common/ext_sigm/sigmmax,sigmmin,sigmw
      common/dist_sigm/xsigm(1000),fxsigm(20,1000,20),fxsigmtot(20,1000)
      common/inp_sigm/m_sigm
      common/div_sigm/ndiv_sigm
      common/nfnb/fb(1000,-1:1,5)
      common/polarised/polcross(20,-1:1,-1:1),polerror(20,-1:1,-1:1)
      common/FB/asycross(20,-1:1),asyerror(20,-1:1)
      common/FB2/asy2cross(20,-1:1),asy2error(20,-1:1)
      common/ACENT/accross(20,-1:1),acerror(20,-1:1)
      common/AFORW/afcross(20,-1:1),aferror(20,-1:1)
      common/ARFB/arfbcross(20,-1:1),arfberror(20,-1:1)
      common/AOFB/aofbcross(20,-1:1),aofberror(20,-1:1)
      common/AFBST/afbstcross(20,-1:1),afbsterror(20,-1:1)
      common/AFBST1/afbst1cross(20,-1:1,8),afbst1error(20,-1:1,8)
      common/ACENTcut/yc,yff,pff
      common/AFBSTcut/yffs(8)
cc Array for looping over x1 <-> x2
      dimension ex(2,2)
      dimension x(100)
      dimension q(4,4),qcol(4,4)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3)
      dimension resgg(-1:1,-1:1),resqq(-1:1,-1:1),resqqalt(-1:1,-1:1)   
      dimension resdd(-1:1,-1:1),resuu(-1:1,-1:1),resss(-1:1,-1:1),
     &          rescc(-1:1,-1:1),resbb(-1:1,-1:1)      
      dimension resddalt(-1:1,-1:1),resuualt(-1:1,-1:1),
     &          resssalt(-1:1,-1:1),resccalt(-1:1,-1:1),
     &          resbbalt(-1:1,-1:1)     
      dimension fx1(13),fx2(13)
      dimension pfx(-1:1,-1:1)
      dimension pfxalt(-1:1,-1:1)
      dimension weight(20,-1:1,-1:1) 
      dimension fasycross(20,-1:1),fasyerror(20,-1:1),
     &          fasy2cross(20,-1:1),fasy2error(20,-1:1),
     &          faccross(20,-1:1),facerror(20,-1:1),
     &          fafcross(20,-1:1),faferror(20,-1:1),
     &          farfbcross(20,-1:1),farfberror(20,-1:1),
     &          faofbcross(20,-1:1),faofberror(20,-1:1),
     &          fafbstcross(20,-1:1),fafbsterror(20,-1:1),
     &          fafbst1cross(20,-1:1,8),fafbst1error(20,-1:1,8),
     &          fpolcross(20,-1:1,-1:1),fpolerror(20,-1:1,-1:1)
      parameter (pi=3.14159265358979323846d0)
c internal random number seed.
      data jseed/987654321/
      fxncosx=0.d0

ccccccccccccccccccccccccccccccccccccccccccccc
c either total energy available OR cutoff M_S
c rmassmax is either ecm_coll or ams depending on ncutmtt switch (0 or 1)
      ecm_max=rmassmax     
ccccccccccccccccccccccccccccccccccccccccccccc

c final state masses are non-zero (y=/=eta).
      rm3=rmf
      rm4=rmf
c allow for implementation of invariant mass cut
      if (ainvmasscut.eq.0.d0)then
      ecm_min = rm3+rm4            
        else 
        ecm_min = ainvmasscut 
        end if
      ecm=x(2)*(ecm_max-ecm_min) + ecm_min
      qcm2=((ecm*ecm-rm3*rm3-rm4*rm4)**2-(2.d0*rm3*rm4)**2)/
     &    (4.d0*ecm*ecm)
      if (qcm2.lt.0.d0) then
        fxncosx=0.d0
        return
      else
        qcm=sqrt(qcm2)
      endif
c x1 and x2 of the partons.
      shat=ecm*ecm
      tau=shat/s
      xx1=x(1)*(1.d0-tau)+tau
      xx2=tau/xx1
      ex(1,1)=xx1
      ex(1,2)=xx2
      ex(2,1)=xx2
      ex(2,2)=xx1      
      do ix=1,2
      x1 = ex(ix,1)
      x2 = ex(ix,2)
      fffxn=0.d0
c recontruct kinematics (event is planar and symmetric around z-axis).
      phi=2.d0*pi*ran(jseed)
c symmetry around ct.
c generates a symmetric set of ct points
cccccc START COS THETA LOOP cccccc
      do ictloop = -nctpoints,nctpoints
ccccccccccccccccccccccccccccccccccccccccccccc
c reset internal iterative quantities
      fxn1=0.d0
      fxn2=0.d0
      ffxn = 0.d0
      do i=1,20
        do iphel=-1,+1,2
          do jphel=-1,+1,2
            fpolcross(i,iphel,jphel)=0.d0
            fpolerror(i,iphel,jphel)=0.d0
          end do
        end do
        do iasy=-1,+1,2
          fasycross(i,iasy)=0.d0
          fasyerror(i,iasy)=0.d0
          fasy2cross(i,iasy)=0.d0
          fasy2error(i,iasy)=0.d0
          faccross(i,iasy)=0.d0
          facerror(i,iasy)=0.d0
          fafcross(i,iasy)=0.d0
          faferror(i,iasy)=0.d0
          fafbstcross(i,iasy)=0.d0
          fafbsterror(i,iasy)=0.d0
          do j=1,8
            fafbst1cross(i,iasy,j)=0.d0
            fafbst1error(i,iasy,j)=0.d0
          enddo
        end do 
      end do
c initialise.
      do ii=-1,1,2
        do jj=-1,1,2
          resgg(ii,jj)=0.d0
          resqq(ii,jj)=0.d0
          resqqalt(ii,jj)=0.d0
          resdd(ii,jj)=0.d0
          resddalt(ii,jj)=0.d0
          resuu(ii,jj)=0.d0
          resuualt(ii,jj)=0.d0
          resss(ii,jj)=0.d0
          resssalt(ii,jj)=0.d0
          rescc(ii,jj)=0.d0
          resccalt(ii,jj)=0.d0
          resbb(ii,jj)=0.d0
          resbbalt(ii,jj)=0.d0
          do kk=1,20
            weight(kk,ii,jj)=0.d0
          end do
        end do
      end do
cccccccccccccccccccccccccccccccccccccccccccccc
      ct= float(ictloop)/float(nctpoints)
      st=sqrt(1.d0-ct*ct)
      q(4,3)=sqrt(qcm**2+rmf**2)
      q(3,3)=qcm*ct
      q(2,3)=qcm*st*cos(phi)
      q(1,3)=qcm*st*sin(phi)
      q(4,4)=+q(4,3)
      q(3,4)=-q(3,3)
      q(2,4)=-q(2,3)
      q(1,4)=-q(1,3)
c this is the scale for the PDFs.
      if(QQ.eq.0.d0)then
          QQ=2.d0*rmf
      else if((QQ.lt.rmZ).and.(QQ.gt.0.d0))then 
          QQ=rmZ
      else if(QQ.eq.-1d0)then ! Dynamical scale set by QQ=-1.0
          QQ = ecm
      endif
c construct polarised hadronic structure functions.
      IF(ISTRUCTURE.LE.4)THEN
        Q2=QQ*QQ
        IF((X1.LE.1.D-6).OR.(X1.GE.1.D0))THEN
          FXNCOSX=0.D0
          RETURN
        END IF
        IF((X2.LE.1.D-6).OR.(X2.GE.1.D0))THEN
          FXNCOSX=0.D0
          RETURN
        END IF
        IF((QQ.LE.1.3D0).OR.(QQ.GE.1.D4))THEN
          FXNCOSX=0.D0
          RETURN
        END IF
        U1=X1*Ctq6Pdf(1,X1,QQ)
        D1=X1*Ctq6Pdf(2,X1,QQ)
        USEA1=X1*Ctq6Pdf(-1,X1,QQ)
        DSEA1=X1*Ctq6Pdf(-2,X1,QQ)
        STR1=X1*Ctq6Pdf(3,X1,QQ)
        CHM1=X1*Ctq6Pdf(4,X1,QQ)
c        BTM1=0.d0
        BTM1=X1*Ctq6Pdf(5,X1,QQ)
        G1=X1*Ctq6Pdf(0,X1,QQ)
        U2=X2*Ctq6Pdf(1,X2,QQ)
        D2=X2*Ctq6Pdf(2,X2,QQ)
        USEA2=X2*Ctq6Pdf(-1,X2,QQ)
        DSEA2=X2*Ctq6Pdf(-2,X2,QQ)
        STR2=X2*Ctq6Pdf(3,X2,QQ)
        CHM2=X2*Ctq6Pdf(4,X2,QQ)
c        BTM2=0.d0
        BTM2=X2*Ctq6Pdf(5,X2,QQ)
        G2=X2*Ctq6Pdf(0,X2,QQ)
      END IF
c actual PDFs.
      IF(ISTRUCTURE.LE.4)THEN 
        fx1(1)=D1
        fx1(2)=U1
        fx1(3)=STR1
        fx1(4)=CHM1
        fx1(5)=BTM1
        fx1(6)=0.D0
        fx1(7)=DSEA1
        fx1(8)=USEA1
        fx1(9)=fx1(3)
        fx1(10)=fx1(4)
        fx1(11)=fx1(5)
        fx1(12)=fx1(6)
        fx1(13)=G1
        do i=1,13
          fx1(i)=fx1(i)/x1
        end do
        fx2(1)=D2*(1-icoll)+DSEA2*icoll
        fx2(2)=U2*(1-icoll)+USEA2*icoll
        fx2(3)=STR2
        fx2(4)=CHM2
        fx2(5)=BTM2
        fx2(6)=0.D0
        fx2(7)=D2*icoll+DSEA2*(1-icoll)
        fx2(8)=U2*icoll+USEA2*(1-icoll)
        fx2(9)=fx2(3)
        fx2(10)=fx2(4)
        fx2(11)=fx2(5)
        fx2(12)=fx2(6)
        fx2(13)=G2
        do i=1,13
          fx2(i)=fx2(i)/x2
        end do
      END IF
c incoming momenta (massless).
      pcm=ecm/2.d0
      q(4,1)=pcm
      q(3,1)=pcm
      q(2,1)=0.d0
      q(1,1)=0.d0
      q(4,2)=pcm
      q(3,2)=-pcm
      q(2,2)=0.d0
      q(1,2)=0.d0
c cuts in pseudorapidity.
      rps3=q(3,3)/sqrt(q(1,3)**2+q(2,3)**2+q(3,3)**2)
      if(rps3.lt.-1.d0)rps3=-1.d0
      if(rps3.gt.+1.d0)rps3=+1.d0
      rpl3=dacos(rps3)
      arg3=tan(rpl3/2.d0)
      if(arg3.le.0.d0)arg3=1.d-9
      eta3=-log(arg3)
      eta3alt=-eta3
      rps4=q(3,4)/sqrt(q(1,4)**2+q(2,4)**2+q(3,4)**2)
      if(rps4.lt.-1.d0)rps4=-1.d0
      if(rps4.gt.+1.d0)rps4=+1.d0
      rpl4=dacos(rps4)
      arg4=tan(rpl4/2.d0)
      if(arg4.le.0.d0)arg4=1.d-9
      eta4=-log(arg4)
      eta4alt=-eta4
c      if(abs(eta).gt.y)then
c        fxncosx=0.d0
c        return
c      end if
      y3 = 0.5d0*dlog((q(4,3)+q(3,3))/(q(4,3)-q(3,3))) 
      y3alt = -y3
      y4 = 0.5d0*dlog((q(4,4)+q(3,4))/(q(4,4)-q(3,4))) 
      y4alt = -y4 
c MadGraph momenta.
      do i=1,3
        p1(i)=q(i,1)
        p2(i)=q(i,2)
        p3(i)=q(i,3)
        p4(i)=q(i,4)
      end do
      p1(0)=q(4,1)
      p2(0)=q(4,2)
      p3(0)=q(4,3)
      p4(0)=q(4,4)
c MANDELSTAMS
      that =(p1(0)-p3(0))**2.d0
      uhat =(p1(0)-p4(0))**2.d0
      do i=1,3
      that = that - (p1(i)-p3(i))**2.d0   
      uhat = uhat - (p1(i)-p4(i))**2.d0   
      end do
ccccccccccccccccccccccc
      goto 777
            write(*,*) 'root s = ',ecm_coll
            write(*,*) 'root shat = ',ecm
        write(*,*) 'p1  =',p1
        write(*,*) 'p2  =',p2
        write(*,*) 'p3  =',p3
        write(*,*) 'p4  =',p4
        delta_E=p1(0)+p2(0)
     &         -p3(0)-p4(0)
        delta_x=p1(1)+p2(1)
     &         -p3(1)-p4(1)
        delta_y=p1(2)+p2(2)
     &         -p3(2)-p4(2)
        delta_z=p1(3)+p2(3)
     &         -p3(3)-p4(3)
        write(*,*) 'delta_E=',delta_E
        write(*,*) 'delta_x=',delta_x
        write(*,*) 'delta_y=',delta_y
        write(*,*) 'delta_z=',delta_z
        rmassa1=sqrt(abs(p1(0)**2-p1(1)**2-p1(2)**2-p1(3)**2))
        rmassa2=sqrt(abs(p2(0)**2-p2(1)**2-p2(2)**2-p2(3)**2))
        rmassa3=sqrt(abs(p3(0)**2-p3(1)**2-p3(2)**2-p3(3)**2))
        rmassa4=sqrt(abs(p4(0)**2-p4(1)**2-p4(2)**2-p4(3)**2))
        write(*,*) 'rm_1 =',rmassa1
        write(*,*) 'rm_2 =',rmassa2
        write(*,*) 'rm_3 =',rmassa3
        write(*,*) 'rm_4 =',rmassa4
 777    continue
ccccccccccccccccccccccc
c initial and final state momenta in the collider CM.
      vcol=(x1-x2)/(x1+x2)
      gcol=(x1+x2)/2.d0/sqrt(x1*x2)
      do i=1,4
        qcol(4,i)=gcol*(q(4,i)+vcol*q(3,i))
        qcol(3,i)=gcol*(q(3,i)+vcol*q(4,i))
        qcol(2,i)=q(2,i)
        qcol(1,i)=q(1,i)
      end do
      
      if(qcol(3,3).gt.qcol(4,3)) print *, "Pz > E for particle 3"
      if(qcol(3,4).gt.qcol(4,4)) print *, "Pz > E for particle 4"
      
c collider frame pseudorapidity 
c particle 3
      rpscol3=qcol(3,3)/sqrt(qcol(1,3)**2+qcol(2,3)**2+qcol(3,3)**2)
      if(rpscol3.lt.-1.d0)rpscol3=-1.d0
      if(rpscol3.gt.+1.d0)rpscol3=+1.d0
      rplcol3=dacos(rpscol3)
      argcol3=tan(rplcol3/2.d0)
      if(argcol3.le.0.d0)argcol3=1.d-9
      eta3col=-log(argcol3)
      ycol3 = 0.5d0*dlog((qcol(4,3)+qcol(3,3))/(qcol(4,3)-qcol(3,3))) 
c particle 4
      rpscol4=qcol(3,4)/sqrt(qcol(1,4)**2+qcol(2,4)**2+qcol(3,4)**2)
      if(rpscol4.lt.-1.d0)rpscol4=-1.d0
      if(rpscol4.gt.+1.d0)rpscol4=+1.d0
      rplcol4=dacos(rpscol4)
      argcol4=tan(rplcol4/2.d0)
      if(argcol4.le.0.d0)argcol4=1.d-9
      eta4col=-log(argcol4)
      ycol4 = 0.5d0*dlog((qcol(4,4)+qcol(3,4))/(qcol(4,4)-qcol(3,4)))
c eta cuts
      if(dabs(ycol3).gt.y.OR.dabs(ycol4).gt.y)then
         ffxn=0.d0
         goto 333
c        return
      end if
c deta cuts
      if (rmf.gt.0.d0) then
      if(dabs(ycol3-ycol4).gt.dycut)then
        ffxn=0.d0
        return
      end if
      else
      if(dabs(eta3col-eta4col).gt.dycut)then
        ffxn=0.d0
        return
      end if
      end if
c pT cut
      pT=sqrt(qcol(1,3)**2+qcol(2,3)**2)
      if(pT.lt.ptCut)then
        ffxn=0.d0
        goto 333
c        return
      end if
c CM boost
      ycm = 0.5d0*dlog(x1/x2)
      pzcm= qcol(3,3)+qcol(3,4)
c calculates alphas.
      a_s=alfas(QQ,rlambdaQCD4,nloops)      
      gs2=4.d0*pi*a_s
      gs=sqrt(gs2)
cccccccccccccccccccc
      CALL ME(CONT,QCD,EW,BSM,gs,rmf,p1,p2,p3,p4,resgg,
     &      resqq,resqqalt,resuu,resuualt,resdd,resddalt,
     &      rescc,resccalt,resss,resssalt,resbb,resbbalt)
ccc MEs.
c initial luminosity for total unpolarised cross section.
      pfxtot=0.d0
      pfxalttot = 0.d0
      do khel=-1,1,2
        do lhel=-1,1,2
          pfx(khel,lhel)=resgg(khel,lhel) *fx1(13)*fx2(13)
     &+(resqq(khel,lhel)+resdd(khel,lhel))*fx1( 1)*fx2( 7)
     &+(resqq(khel,lhel)+resuu(khel,lhel))*fx1( 2)*fx2( 8)
     &+(resqq(khel,lhel)+resss(khel,lhel))*fx1( 3)*fx2( 9)
     &+(resqq(khel,lhel)+rescc(khel,lhel))*fx1( 4)*fx2(10)
     &+(resbb(khel,lhel))*fx1( 5)*fx2(11)
          pfxalt(khel,lhel)= 
     &+(resqqalt(khel,lhel)+resddalt(khel,lhel))*fx1( 7)*fx2( 1)
     &+(resqqalt(khel,lhel)+resuualt(khel,lhel))*fx1( 8)*fx2( 2)
     &+(resqqalt(khel,lhel)+resssalt(khel,lhel))*fx1( 9)*fx2( 3)
     &+(resqqalt(khel,lhel)+resccalt(khel,lhel))*fx1(10)*fx2( 4)
     &+(resbbalt(khel,lhel))*fx1(11)*fx2( 5)
          pfxtot=pfxtot
     &          +pfx(khel,lhel)
          pfxalttot=pfxalttot
     &          +pfxalt(khel,lhel)
         end do
       end do
      if(pfxtot.eq.0.d0.and.pfxalttot.eq.0.d0)then
        ffxn=0.d0
        goto 333
c        return
      end if
c for distributions, 
      do khel=-1,1,2
        do lhel=-1,1,2
          pfx(khel,lhel)=pfx(khel,lhel)/(pfxtot+pfxalttot)
          pfxalt(khel,lhel)=pfxalt(khel,lhel)/(pfxtot+pfxalttot)
        end do
      end do      
c Jacobians from dx1 dx2 -> dx(1) dEcm.
      if (ix.eq.1) then
      pfxtot=pfxtot*(1.d0-tau)*2.d0*ecm/s/x1
     &      *(ecm_max-ecm_min)
      pfxalttot=pfxalttot*(1.d0-tau)*2.d0*ecm/s/x1
     &      *(ecm_max-ecm_min)
      else
      pfxtot=pfxtot*(1.d0-tau)*2.d0*ecm/s/x2
     &      *(ecm_max-ecm_min)
      pfxalttot=pfxalttot*(1.d0-tau)*2.d0*ecm/s/x2
     &      *(ecm_max-ecm_min)
      endif
c MEs and PDFs.
      fxn1=pfxtot
      fxn2=pfxalttot
c constant factor.
      fxn1=fxn1*fac
      fxn2=fxn2*fac
c jacobians by hand, pi's and flux (1/2s)
cccccccccccccccccccccc
      fxn1=fxn1*qcm/(2.d0*pcm)*2.d0**(4-3*(2))
      fxn1=fxn1/2.d0/ecm/ecm*(2.d0*pi)**(4-3*(2))
      fxn2=fxn2*qcm/(2.d0*pcm)*2.d0**(4-3*(2))
      fxn2=fxn2/2.d0/ecm/ecm*(2.d0*pi)**(4-3*(2))
      fxn1=fxn1/(2.d0*float(nctpoints)+1.d0)
      fxn2=fxn2/(2.d0*float(nctpoints)+1.d0)
      ffxn = fxn1+fxn2
cccccccccccccccccccccc
c polarised cross sections.
      do iphel=-1,+1,2
        do jphel=-1,+1,2
          fpolcross(it,iphel,jphel)=fpolcross(it,iphel,jphel)
     &                            +ffxn
     &                            *wgt          
     &                            *(pfx(iphel,jphel)         
     &                            +pfxalt(iphel,jphel))
          weight(it,iphel,jphel)=+ffxn
     &                           *wgt
     &                           *(pfx(iphel,jphel)
     &                           +pfxalt(iphel,jphel))
          fpolerror(it,iphel,jphel)=fpolerror(it,iphel,jphel)
     &                            +fpolcross(it,iphel,jphel)**2
        end do
      end do
c FB asymmetry. IN CM FRAME!
        cost4=(q(1,4)*q(1,1)
     &       +q(2,4)*q(2,1)
     &       +q(3,4)*q(3,1))
     &  /sqrt(q(1,4)*q(1,4)
     &       +q(2,4)*q(2,4)
     &       +q(3,4)*q(3,4))
     &  /sqrt(q(1,1)*q(1,1)
     &       +q(2,1)*q(2,1)
     &       +q(3,1)*q(3,1))
        cost4alt=-cost4
      if(cost4.eq.0.d0)then
        continue
      else if(cost4.gt.0.d0)then
        fasy2cross(it,+1)=fasy2cross(it,+1)
     &                            +fxn1
     &                            *wgt          
        fasy2error(it,+1)=fasy2error(it,+1)
     &                 +fasy2cross(it,+1)**2
        fasy2cross(it,-1)=fasy2cross(it,-1)
     &                            +fxn2
     &                            *wgt          
        fasy2error(it,-1)=fasy2error(it,-1)
     &                 +fasy2cross(it,-1)**2
      else if(cost4.lt.0.d0)then
        fasy2cross(it,-1)=fasy2cross(it,-1)
     &                            +fxn1
     &                            *wgt          
        fasy2error(it,-1)=fasy2error(it,-1)
     &                 +fasy2cross(it,-1)**2
        fasy2cross(it,+1)=fasy2cross(it,+1)
     &                            +fxn2
     &                            *wgt          
        fasy2error(it,+1)=fasy2error(it,+1)
     &                 +fasy2cross(it,-1)**2
      end if 
c FB asymmetry. IN Collider FRAME!
        cost4col=(qcol(1,4)*qcol(1,1)
     &       +qcol(2,4)*qcol(2,1)
     &       +qcol(3,4)*qcol(3,1))
     &  /sqrt(qcol(1,4)*qcol(1,4)
     &       +qcol(2,4)*qcol(2,4)
     &       +qcol(3,4)*qcol(3,4))
     &  /sqrt(qcol(1,1)*qcol(1,1)
     &       +qcol(2,1)*qcol(2,1)
     &       +qcol(3,1)*qcol(3,1))
      if(cost4col.eq.0.d0)then
        continue
      else if(cost4col.gt.0.d0)then
        fasycross(it,+1)=fasycross(it,+1)
     &                            +ffxn
     &                            *wgt          
        fasyerror(it,+1)=fasyerror(it,+1)
     &                 +fasycross(it,+1)**2
      else if(cost4col.lt.0.d0)then
        fasycross(it,-1)=fasycross(it,-1)
     &                            +ffxn
     &                            *wgt          
        fasyerror(it,-1)=fasyerror(it,-1)
     &                 +fasycross(it,-1)**2
      end if 
      
c Central Asymmetry
      if(dabs(ycol4).le.yc)then
        faccross(it,+1)=faccross(it,+1)
     &                            +ffxn
     &                            *wgt          
        facerror(it,+1)=facerror(it,+1)
     &                 +faccross(it,+1)**2
      end if
      if(dabs(ycol3).le.yc)then
        faccross(it,-1)=faccross(it,-1)
     &                            +ffxn
     &                            *wgt          
        facerror(it,-1)=facerror(it,-1)
     &                 +faccross(it,-1)**2
      end if
c Forward Asymmetry
      if(dabs(ycol4).gt.yc)then
        fafcross(it,+1)=fafcross(it,+1)
     &                            +ffxn
     &                            *wgt          
        faferror(it,+1)=faferror(it,+1)
     &                 +fafcross(it,+1)**2
      end if
      if(dabs(ycol3).gt.yc)then
        fafcross(it,-1)=fafcross(it,-1)
     &                            +ffxn
     &                            *wgt          
        faferror(it,-1)=faferror(it,-1)
     &                 +fafcross(it,-1)**2
      end if
c ARFB/AOFB (rapidity dependent and one sided AFB)
      if(dabs(ycol4).gt.dabs(ycol3))then
        if(abs(ycm).gt.yff)then
        farfbcross(it,+1)=farfbcross(it,+1)
     &                            +ffxn
     &                            *wgt          
        farfberror(it,+1)=farfberror(it,+1)
     &                 +farfbcross(it,+1)**2
        end if
        if(abs(pzcm).gt.pff)then
        faofbcross(it,+1)=faofbcross(it,+1)
     &                            +ffxn
     &                            *wgt          
        faofberror(it,+1)=faofberror(it,+1)
     &                 +faofbcross(it,+1)**2
        end if
      end if
      if(dabs(ycol3).gt.dabs(ycol4))then
        if(abs(ycm).gt.yff)then
        farfbcross(it,-1)=farfbcross(it,-1)
     &                            +ffxn
     &                            *wgt          
        farfberror(it,-1)=farfberror(it,-1)
     &                 +farfbcross(it,-1)**2
        end if
        if(abs(pzcm).gt.pff)then
        faofbcross(it,-1)=faofbcross(it,-1)
     &                            +ffxn
     &                            *wgt          
        faofberror(it,-1)=faofberror(it,-1)
     &                 +faofbcross(it,-1)**2
        end if      
      end if
c AFBSTAR
         costst=int(ycm/dabs(ycm))*cost4
      if(costst.gt.0.d0)then
         fafbstcross(it,+1)=fafbstcross(it,+1)
     &                            +ffxn
     &                            *wgt          
         fafbsterror(it,+1)=fafbsterror(it,+1)
     &                 +fafbstcross(it,+1)**2
         do j=1,8
           if (dabs(ycm).gt.yffs(j)) then
             fafbst1cross(it,+1,j)=fafbst1cross(it,+1,j)
     &                              +ffxn
     &                              *wgt          
             fafbst1error(it,+1,j)=fafbst1error(it,+1,j)
     &                   +fafbst1cross(it,+1,j)**2
           endif
         enddo
       else if(costst.lt.0.d0)then
         fafbstcross(it,-1)=fafbstcross(it,-1)
     &                            +ffxn
     &                            *wgt          
         fafbsterror(it,-1)=fafbsterror(it,-1)
     &                 +fafbstcross(it,-1)**2
         do j=1,8
             if (dabs(ycm).gt.yffs(j))then
               fafbst1cross(it,-1,j)=fafbst1cross(it,-1,j)
     &                            +ffxn
     &                            *wgt          
               fafbst1error(it,-1,j)=fafbst1error(it,-1,j)
     &                 +fafbst1cross(it,-1,j)**2
             endif
         enddo
       end if
c histograms.
      hist1=fxn1*wgt
      hist2=fxn2*wgt
      hist=hist1+hist2
      if(m_pT.eq.1)then
c generate distribution in pT.
        pT=sqrt(qcol(1,3)**2+qcol(2,3)**2)
        nbin=int((pT-pTmin)/pTw)+1
        if(nbin.ge.(ndiv_pT+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxpT(nbin,it)=fxpT(nbin,it)+hist
        end if
      end if  
      if(m_rmass.eq.1)then
c generate distribution in rmass.
        rmass2=(qcol(4,3)+qcol(4,4))**2
        do i=1,3
          rmass2=rmass2-(qcol(i,3)+qcol(i,4))**2
        end do
        rmass=sqrt(abs(rmass2))
        nbin=int((rmass-rmassmin)/rmassw)+1
        if(nbin.ge.(ndiv_rmass+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxrmass(nbin,it)=fxrmass(nbin,it)+hist
        end if
      end if  
      if(m_eta.eq.1)then
c generate distribution in eta3.
        nbin=int((eta3-etamin)/etaw)+1
        if(nbin.ge.(ndiv_eta+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxeta3(nbin,it)=fxeta3(nbin,it)+hist1
        end if
        nbin=int((eta3alt-etamin)/etaw)+1
        if(nbin.ge.(ndiv_eta+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxeta3(nbin,it)=fxeta3(nbin,it)+hist2
        end if
c Collider eta
        nbin=int((eta3col-etamin)/etaw)+1
        if(nbin.ge.(ndiv_eta+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxeta3col(nbin,it)=fxeta3col(nbin,it)+hist
        end if
      end if  
      if(m_eta.eq.1)then
c generate distribution in eta4.
        nbin=int((eta4-etamin)/etaw)+1
        if(nbin.ge.(ndiv_eta+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxeta4(nbin,it)=fxeta4(nbin,it)+hist1
        end if
        nbin=int((eta4alt-etamin)/etaw)+1
        if(nbin.ge.(ndiv_eta+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxeta4(nbin,it)=fxeta4(nbin,it)+hist2
        end if
c Collider eta
        nbin=int((eta4col-etamin)/etaw)+1
        if(nbin.ge.(ndiv_eta+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxeta4col(nbin,it)=fxeta4col(nbin,it)+hist
        end if
      end if  
      if(m_beta.eq.1)then
c generate distribution in beta.
        beta=shat/4.d0/rmf**2-1.d0
        nbin=int((beta-betamin)/betaw)+1
        if(nbin.ge.(ndiv_beta+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxbeta(nbin,it)=fxbeta(nbin,it)+hist
        end if
      end if  
      if(m_cost.eq.1)then
c generate distribution in cost.
c collider antitop
        cost3col=(qcol(1,3)*qcol(1,1)
     &       +qcol(2,3)*qcol(2,1)
     &       +qcol(3,3)*qcol(3,1))
     &  /sqrt(qcol(1,3)*qcol(1,3)
     &       +qcol(2,3)*qcol(2,3)
     &       +qcol(3,3)*qcol(3,3))
     &  /sqrt(qcol(1,1)*qcol(1,1)
     &       +qcol(2,1)*qcol(2,1)
     &       +qcol(3,1)*qcol(3,1))
        nbin=int((cost3col-costmin)/costw)+1
        if(nbin.ge.(ndiv_cost+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxcost3col(nbin,it)=fxcost3col(nbin,it)+hist
        end if
c collider top
        cost4col=(qcol(1,4)*qcol(1,1)
     &       +qcol(2,4)*qcol(2,1)
     &       +qcol(3,4)*qcol(3,1))
     &  /sqrt(qcol(1,4)*qcol(1,4)
     &       +qcol(2,4)*qcol(2,4)
     &       +qcol(3,4)*qcol(3,4))
     &  /sqrt(qcol(1,1)*qcol(1,1)
     &       +qcol(2,1)*qcol(2,1)
     &       +qcol(3,1)*qcol(3,1))
        nbin=int((cost4col-costmin)/costw)+1
        if(nbin.ge.(ndiv_cost+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxcost4col(nbin,it)=fxcost4col(nbin,it)+hist
        end if
c CM antitop
        cost3=(q(1,3)*q(1,1)
     &       +q(2,3)*q(2,1)
     &       +q(3,3)*q(3,1))
     &  /sqrt(q(1,3)*q(1,3)
     &       +q(2,3)*q(2,3)
     &       +q(3,3)*q(3,3))
     &  /sqrt(q(1,1)*q(1,1)
     &       +q(2,1)*q(2,1)
     &       +q(3,1)*q(3,1))
      cost3alt = -cost3
        nbin=int((cost3-costmin)/costw)+1
        if(nbin.ge.(ndiv_cost+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxcost3(nbin,it)=fxcost3(nbin,it)+hist1
        end if
        nbin=int((cost3alt-costmin)/costw)+1
        if(nbin.ge.(ndiv_cost+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxcost3(nbin,it)=fxcost3(nbin,it)+hist2
        end if
c CM top
        cost4=(q(1,4)*q(1,1)
     &       +q(2,4)*q(2,1)
     &       +q(3,4)*q(3,1))
     &  /sqrt(q(1,4)*q(1,4)
     &       +q(2,4)*q(2,4)
     &       +q(3,4)*q(3,4))
     &  /sqrt(q(1,1)*q(1,1)
     &       +q(2,1)*q(2,1)
     &       +q(3,1)*q(3,1))
        cost4alt=-cost4
        nbin=int((cost4-costmin)/costw)+1
        if(nbin.ge.(ndiv_cost+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxcost4(nbin,it)=fxcost4(nbin,it)+hist1
        end if
        nbin=int((cost4alt-costmin)/costw)+1
        if(nbin.ge.(ndiv_cost+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxcost4(nbin,it)=fxcost4(nbin,it)+hist2
        end if
      end if  
      if(m_Et.eq.1)then
c generate distribution in Et.
        Et=qcol(4,4)
        nbin=int((Et-Etmin)/Etw)+1
        if(nbin.ge.(ndiv_Et+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxEt(nbin,it)=fxEt(nbin,it)+hist
        end if
      end if  
      
      if(m_Chi.eq.1)then
c generate distribution in Chi.
        Chi = dexp(dabs(ycol3-ycol4))
        nbin=int((Chi-Chimin)/Chiw)+1
        if(nbin.ge.(ndiv_Chi+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxChi(nbin,it)=fxChi(nbin,it)+hist
        end if
      end if  
      if(m_Chi2.eq.1)then
c generate distribution in Chi2.
        Chi2 = dexp(dabs(ycol3-ycol4))
        nbin=int((Chi2-Chi2min)/Chi2w)+1
        if(nbin.ge.(ndiv_Chi2+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxChi2(nbin,it)=fxChi2(nbin,it)+hist
        end if
      end if 

      if(m_Y.eq.1)then
c generate distribution in Y3.
        nbin=int((Y3-Ymin)/Yw)+1
        if(nbin.ge.(ndiv_Y+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxY3(nbin,it)=fxY3(nbin,it)+hist1
        end if
        nbin=int((Y3alt-Ymin)/Yw)+1
        if(nbin.ge.(ndiv_Y+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxY3(nbin,it)=fxY3(nbin,it)+hist2
        end if
c generate distribution in Y4.
        nbin=int((Y4-Ymin)/Yw)+1
        if(nbin.ge.(ndiv_Y+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxY4(nbin,it)=fxY4(nbin,it)+hist1
        end if
        nbin=int((Y4alt-Ymin)/Yw)+1
        if(nbin.ge.(ndiv_Y+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxY4(nbin,it)=fxY4(nbin,it)+hist2
        end if
      end if
      if(m_Ycol.eq.1)then
c generate distribution in Ycol3.
        nbin=int((Ycol3-Ycolmin)/Ycolw)+1
        if(nbin.ge.(ndiv_Ycol+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxYcol3(nbin,it)=fxYcol3(nbin,it)+hist
        end if
c generate distribution in Ycol4.
        nbin=int((Ycol4-Ycolmin)/Ycolw)+1
        if(nbin.ge.(ndiv_Ycol+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxYcol4(nbin,it)=fxYcol4(nbin,it)+hist
        end if
      end if
      if(m_DelY.eq.1)then
c generate distribution in DelY.
        DelY = ycol4-ycol3
        nbin=int((DelY-DelYmin)/DelYw)+1
        if(nbin.ge.(ndiv_DelY+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxDelY(nbin,it)=fxDelY(nbin,it)+hist
        end if
      end if

      if(m_Ytt.eq.1)then
c generate distribution in Ytt.
        nbin=int((ycm-Yttmin)/Yttw)+1
        if(nbin.ge.(ndiv_Ytt+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxYtt(nbin,it)=fxYtt(nbin,it)+hist
        end if
      end if

      x1x2=X1/X2
      if(m_x1x2.eq.1)then
c generate distribution in x1x2.
        nbin=int(x1x2-x1x2min)/x1x2w+1
        if(nbin.ge.(ndiv_x1x2+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxx1x2(nbin,it)=fxx1x2(nbin,it)+hist
        end if
      end if
     

      if(m_Pztt.eq.1)then
c generate distribution in Pztt.
        nbin=int((dabs(pzcm)-Pzttmin)/Pzttw)+1
        if(nbin.ge.(ndiv_Pztt+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxPztt(nbin,it)=fxPztt(nbin,it)+hist
        end if
      end if
     
      if(m_sigp.eq.1)then
c generate distribution in sigp for ALL.
        sigp=ecm
        nbin=int((sigp-sigpmin)/sigpw)+1
        if(nbin.ge.(ndiv_sigp+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxsigp(1,nbin,it)=fxsigp(1,nbin,it)+
     &    (weight(it,+1,+1)+weight(it,-1,-1))
        end if
      end if  
      if(m_sigm.eq.1)then
c generate distribution in sigm for ALL.
        sigm=ecm
        nbin=int((sigm-sigmmin)/sigmw)+1
        if(nbin.ge.(ndiv_sigm+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxsigm(1,nbin,it)=fxsigm(1,nbin,it)+
     &    (weight(it,+1,-1)+weight(it,-1,+1))
        end if
      end if  
      if(m_sigp.eq.1)then
c generate distribution in sigp for AL.
        sigp=ecm
        nbin=int((sigp-sigpmin)/sigpw)+1
        if(nbin.ge.(ndiv_sigp+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxsigp(2,nbin,it)=fxsigp(2,nbin,it)+
     &    (weight(it,-1,-1)+weight(it,-1,+1))
        end if
      end if  
      if(m_sigm.eq.1)then
c generate distribution in sigm for AL.
        sigm=ecm
        nbin=int((sigm-sigmmin)/sigmw)+1
        if(nbin.ge.(ndiv_sigm+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxsigm(2,nbin,it)=fxsigm(2,nbin,it)+
     &    (weight(it,+1,-1)+weight(it,+1,+1))
        end if
      end if  
      if(m_sigp.eq.1)then
c generate distribution in sigp for APV.
        sigp=ecm
        nbin=int((sigp-sigpmin)/sigpw)+1
        if(nbin.ge.(ndiv_sigp+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxsigp(3,nbin,it)=fxsigp(3,nbin,it)+
     &    (weight(it,-1,-1))
        end if
      end if  
      if(m_sigm.eq.1)then
c generate distribution in sigm for APV.
        sigm=ecm
        nbin=int((sigm-sigmmin)/sigmw)+1
        if(nbin.ge.(ndiv_sigm+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxsigm(3,nbin,it)=fxsigm(3,nbin,it)+
     &    (weight(it,+1,+1))
        end if
      end if  
      if((m_sigp.eq.1).and.(cost4col.gt.0.d0))then
c generate distribution in sigp for AFB.
        sigp=ecm
        nbin=int((sigp-sigpmin)/sigpw)+1
        if(nbin.ge.(ndiv_sigp+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxsigp(4,nbin,it)=fxsigp(4,nbin,it)+hist
        end if
      end if  
      if((m_sigm.eq.1).and.(cost4col.lt.0.d0))then
c generate distribution in sigm for AFB.
        sigm=ecm
        nbin=int((sigm-sigmmin)/sigmw)+1
        if(nbin.ge.(ndiv_sigm+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxsigm(4,nbin,it)=fxsigm(4,nbin,it)+hist
        end if
      end if  
      if((m_sigp.eq.1).and.(cost4.gt.0.d0))then
c generate distribution in sigp for AFB2.
        sigp=ecm
        nbin=int((sigp-sigpmin)/sigpw)+1
        if(nbin.ge.(ndiv_sigp+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxsigp(7,nbin,it)=fxsigp(7,nbin,it)+hist1
          fxsigm(7,nbin,it)=fxsigm(7,nbin,it)+hist2
        end if
      end if  
      if((m_sigm.eq.1).and.(cost4.lt.0.d0))then
c generate distribution in sigm for AFB2.
        sigm=ecm
        nbin=int((sigm-sigmmin)/sigmw)+1
        if(nbin.ge.(ndiv_sigm+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxsigm(7,nbin,it)=fxsigm(7,nbin,it)+hist1
          fxsigp(7,nbin,it)=fxsigp(7,nbin,it)+hist2
        end if
      end if  
      if((m_sigp.eq.1).and.(dabs(ycol4).le.yc))then
c generate distribution in sigp for AC.
        sigp=ecm
        nbin=int((sigp-sigpmin)/sigpw)+1
        if(nbin.ge.(ndiv_sigp+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxsigp(5,nbin,it)=fxsigp(5,nbin,it)+hist
        end if
      end if  
      if((m_sigm.eq.1).and.(dabs(ycol3).le.yc))then
c generate distribution in sigm for AC.
        sigm=ecm
        nbin=int((sigm-sigmmin)/sigmw)+1
        if(nbin.ge.(ndiv_sigm+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxsigm(5,nbin,it)=fxsigm(5,nbin,it)+hist
        end if
      end if
c generate distribution in sigp for AF.
      if((m_sigp.eq.1).and.(dabs(ycol4).gt.yc))then
        sigp=ecm
        nbin=int((sigp-sigpmin)/sigpw)+1
        if(nbin.ge.(ndiv_sigp+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxsigp(6,nbin,it)=fxsigp(6,nbin,it)+hist
        end if
      end if
        
      if((m_sigm.eq.1).and.(dabs(ycol3).gt.yc))then
c generate distribution in sigm for AF.
        sigm=ecm
        nbin=int((sigm-sigmmin)/sigmw)+1
        if(nbin.ge.(ndiv_sigm+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxsigm(6,nbin,it)=fxsigm(6,nbin,it)+hist
        end if
      end if
c generate distribution in sigp for ARFB/AOFB.
      if((m_sigp.eq.1).and.(dabs(ycol4).gt.dabs(ycol3)))then
        sigp=ecm
        nbin=int((sigp-sigpmin)/sigpw)+1
        if(nbin.ge.(ndiv_sigp+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
            if(dabs(ycm).gt.yff)fxsigp(8,nbin,it)=fxsigp(8,nbin,it)+hist
            if(dabs(pzcm).gt.pff)fxsigp(9,nbin,it)=fxsigp(9,nbin,it)+hist         
        end if
      end if
        
      if((m_sigm.eq.1).and.(dabs(ycol3).gt.dabs(ycol4)))then
c generate distribution in sigm for ARFB/AOFB.
        sigm=ecm
        nbin=int((sigm-sigmmin)/sigmw)+1
        if(nbin.ge.(ndiv_sigm+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
            if(dabs(ycm).gt.yff)fxsigm(8,nbin,it)=fxsigm(8,nbin,it)+hist
            if(dabs(pzcm).gt.pff)fxsigm(9,nbin,it)=fxsigm(9,nbin,it)+hist 
        end if
      end if

c generate distribution in sigp for AFBSTAR.
      if((m_sigp.eq.1).and.(costst.gt.0.d0))then
        sigp=ecm
        nbin=int((sigp-sigpmin)/sigpw)+1
        if(nbin.ge.(ndiv_sigp+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
            fxsigp(10,nbin,it)=fxsigp(10,nbin,it)+hist 
            do j=1,8
            if(dabs(ycm).gt.yffs(j))then
                fxsigp(10+j,nbin,it)=fxsigp(10+j,nbin,it)+hist 
            endif
            enddo
        end if
      end if
c generate distribution in sigm for AFBSTAR.
      if((m_sigm.eq.1).and.(costst.lt.0.d0))then
        sigm=ecm
        nbin=int((sigm-sigmmin)/sigmw)+1
        if(nbin.ge.(ndiv_sigm+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
            fxsigm(10,nbin,it)=fxsigm(10,nbin,it)+hist 
            do j=1,8
            if(dabs(ycm).gt.yffs(j))then
                fxsigm(10+j,nbin,it)=fxsigm(10+j,nbin,it)+hist 
            endif
            enddo 
        end if
      end if


      fffxn = fffxn + ffxn
      do i=-1,1,2
        asycross(it,i)=asycross(it,i)+fasycross(it,i)
        asyerror(it,i)=asyerror(it,i)+fasyerror(it,i)
        asy2cross(it,i)=asy2cross(it,i)+fasy2cross(it,i)
        asy2error(it,i)=asy2error(it,i)+fasy2error(it,i)
        accross(it,i)=accross(it,i)+faccross(it,i)
        acerror(it,i)=acerror(it,i)+facerror(it,i)
        afcross(it,i)=afcross(it,i)+fafcross(it,i)
        aferror(it,i)=aferror(it,i)+faferror(it,i)
        arfbcross(it,i)=arfbcross(it,i)+farfbcross(it,i)
        arfberror(it,i)=arfberror(it,i)+farfberror(it,i)
        aofbcross(it,i)=aofbcross(it,i)+faofbcross(it,i)
        aofberror(it,i)=aofberror(it,i)+faofberror(it,i)
        afbstcross(it,i)=afbstcross(it,i)+fafbstcross(it,i)
        afbsterror(it,i)=afbsterror(it,i)+fafbsterror(it,i)
        do j=1,8
          afbst1cross(it,i,j)=afbst1cross(it,i,j)+fafbst1cross(it,i,j)
          afbst1error(it,i,j)=afbst1error(it,i,j)+fafbst1error(it,i,j)
        enddo
        do j=-1,1,2
          polcross(it,i,j)=polcross(it,i,j)+fpolcross(it,i,j)
          polerror(it,i,j)=polerror(it,i,j)+fpolerror(it,i,j)
        end do
      end do
c statistics.
      npoints=npoints+1 
333   continue
c      print*,"ffxn,fffxn= ",ffxn,fffxn
cccccc END COS THETA LOOP cccccc
      end do
      fxncosx=fxncosx+fffxn
cccccc END X LOOP cccccc
      end do
c      print*,"fxn(x,wgt)",fxncosx 
      return
      end
      
c
c     --------------------------------------------------------------
c
c     end program d_pptt_LO

