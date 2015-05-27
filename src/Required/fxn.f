      function fxn(x,wgt)
      implicit real*8 (a-h,o-z)
      include 'vegas.inc'
      include 'ewparams.inc'
      include 'zpparams.inc'
      include 'runparams.inc'
c distributions.
      include 'dists_common.inc'
cc Asymmetries 
      include 'asy_common.inc'
c few extra variables for looping over M_s, calculating CL, invariant mass cut
      common/lum/alumpb
      common/cut/ainvmasscut
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
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c rmassmax is either ecm_coll or ams depending on ncutmtt switch (0 or 1)
      ecm_max=rmassmax     
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
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c factorisation scale
      if(QQ.eq.0.d0)then
          QQ=2.d0*rmf
      else if((QQ.lt.rmZ).and.(QQ.gt.0.d0))then 
          QQ=rmZ
      else if(QQ.eq.-1d0)then ! Dynamical scale set by QQ=-1.0
          QQ = ecm
      endif
c construct polarised hadronic structure functions.
      call setpdf(x1, x2, QQ, istructure, icoll, *626, fxn, fx1, fx2)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c CM frame outgoing momenta
      q(4,3)=sqrt(qcm**2+rmf**2)
      q(3,3)=qcm*ct
      q(2,3)=qcm*st*cos(phi)
      q(1,3)=qcm*st*sin(phi)
      q(4,4)=+q(4,3)
      q(3,4)=-q(3,3)
      q(2,4)=-q(2,3)
      q(1,4)=-q(1,3)
c CM frame incoming momenta (massless).
      pcm=ecm/2.d0
      q(4,1)=pcm
      q(3,1)=pcm
      q(2,1)=0.d0
      q(1,1)=0.d0
      q(4,2)=pcm
      q(3,2)=-pcm
      q(2,2)=0.d0
      q(1,2)=0.d0
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
c initial and final state momenta in the collider (LAB) frame.
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
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C Kinematical quantities for cuts & plots
      rmass = ecm
      pT = sqrt(q(1,3)**2+q(2,3)**2) ! p_T(fbar)
c pT cut
      if(pT.lt.ptCut)then
        fxn=0.d0
        return
      end if
      Et = qcol(4,4)
c collider frame polar angles
      cost3col=qcol(3,3)/sqrt(qcol(1,3)**2+qcol(2,3)**2+qcol(3,3)**2)
      if(cost3col.lt.-1.d0)cost3col=-1.d0
      if(cost3col.gt.+1.d0)cost3col=+1.d0
C 
      cost4col=qcol(3,4)/sqrt(qcol(1,4)**2+qcol(2,4)**2+qcol(3,4)**2)
      if(cost4col.lt.-1.d0)cost4col=-1.d0
      if(cost4col.gt.+1.d0)cost4col=+1.d0
c collider frame pseudorapidity 
      thetacol3=dacos(cost3col)
      argcol3=tan(thetacol3/2.d0)
      if(argcol3.le.0.d0)argcol3=1.d-9
      eta3col=-log(argcol3)
C       
      thetacol4=dacos(cost4col)
      argcol4=tan(thetacol4/2.d0)
      if(argcol4.le.0.d0)argcol4=1.d-9
      eta4col=-log(argcol4)
c collider frame rapidity 
      y3col = 0.5d0*dlog((qcol(4,3)+qcol(3,3))/(qcol(4,3)-qcol(3,3))) 
      y4col = 0.5d0*dlog((qcol(4,4)+qcol(3,4))/(qcol(4,4)-qcol(3,4)))
c collider rapidity cuts
      if(dabs(y3col).gt.y.OR.dabs(y4col).gt.y)then
        fxn=0.d0
        return
      end if
c collider frame rapidity difference
      delY = y4col-y3col
      deltaeta = eta4col-eta3col
c dy/deta cuts
      if (rmf.gt.0.d0) then
        if(dabs(delY).gt.dycut)then
          fxn=0.d0
          return
        end if
      else
        if(dabs(deltaeta).gt.dycut)then
          fxn=0.d0
          return
        end if
      end if
c CM boost
      yff = 0.5d0*dlog(x1/x2) ! y(ffbar)
      pzff= qcol(3,3)+qcol(3,4) ! p_z(ffbar)
c Lorentz factor
      beta=shat/4.d0/rmf**2-1.d0
c CM polar angle.      
      cost3 = ct
      cost3alt = -ct
      cost4 = -ct
      cost4alt = -ct
      costst=int(yff/dabs(yff))*cost4 ! cos thetastar
c CM pseudorapidity.
      theta3=dacos(cost3)
      arg3=tan(theta3/2.d0)
      if(arg3.le.0.d0)arg3=1.d-9
      eta3=-log(arg3)
      eta3alt=-eta3
c  
      theta4=dacos(cost4)
      arg4=tan(theta4/2.d0)
      if(arg4.le.0.d0)arg4=1.d-9
      eta4=-log(arg4)
      eta4alt=-eta4
c CM rapidity.
      y3 = 0.5d0*dlog((q(4,3)+q(3,3))/(q(4,3)-q(3,3))) 
      y3alt = -y3
      y4 = 0.5d0*dlog((q(4,4)+q(3,4))/(q(4,4)-q(3,4))) 
      y4alt = -y4 
ccccccccccccccccccccccc
c calculates alphas.
      a_s=alfas(QQ,rlambdaQCD4,nloops)      
      gs2=4.d0*pi*a_s
      gs=sqrt(gs2)
c compute matrix element
      CALL ME(gs,rmf,p1,p2,p3,p4,resgg,
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
cccccccccccccccccccccc
C histogram weights
      hist1=fxn1*wgt ! weight for q-qbar initial state
      hist2=fxn2*wgt ! weight for qbar-q initial state
      hist=hist1+hist2
C       print*,hist
cccccccccccccccccccccc
c polarised cross sections.
      do iphel=-1,+1,2
        do jphel=-1,+1,2
          polcross(it,iphel,jphel)=polcross(it,iphel,jphel)
     &                            + hist         
     &                            * (pfx(iphel,jphel)         
     &                            + pfxalt(iphel,jphel))
          weight(it,iphel,jphel)=   hist
     &                            * (pfx(iphel,jphel)
     &                            + pfxalt(iphel,jphel))
          polerror(it,iphel,jphel)=polerror(it,iphel,jphel)
     &                            +polcross(it,iphel,jphel)**2
        end do
      end do
cccccccccccccccccccccc
c collect total forward/backward cross sections for asymmetries
      include 'asy_tot.inc'
c fill histograms.
      include 'dists_fill.inc'
c fill binned asymmetries.
      include 'asy_fill.inc'
c statistics.
      npoints=npoints+1  
 626  return
      end
      
c
c     --------------------------------------------------------------
c