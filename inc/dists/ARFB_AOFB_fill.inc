c generate distribution in sigp for ARFB/AOFB.
      if((m_sigp.eq.1).and.(dabs(y4col).gt.dabs(y3col)))then
        sigp=ecm
        nbin=int((sigp-sigpmin)/sigpw)+1
        if(nbin.ge.(ndiv_sigp+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
        if(dabs(yff).gt.yffcut)
     &          fxsigp(iARFB,nbin,it)=fxsigp(iARFB,nbin,it) + hist
        if(dabs(pzff).gt.pzffcut)
     &          fxsigp(iAOFB,nbin,it)=fxsigp(iAOFB,nbin,it) + hist
        end if
      end if
        
      if((m_sigm.eq.1).and.(dabs(y3col).gt.dabs(y4col)))then
c generate distribution in sigm for ARFB/AOFB.
        sigm=ecm
        nbin=int((sigm-sigmmin)/sigmw)+1
        if(nbin.ge.(ndiv_sigm+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
        if(dabs(yff).gt.yffcut)
     &          fxsigm(iARFB,nbin,it)=fxsigm(iARFB,nbin,it) + hist
        if(dabs(pzff).gt.pzffcut)
     &          fxsigm(iAOFB,nbin,it)=fxsigm(iAOFB,nbin,it) + hist
        end if
      end if