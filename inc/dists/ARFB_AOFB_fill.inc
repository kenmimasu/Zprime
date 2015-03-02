c generate distribution in sigp for ARFB/AOFB.
      if((m_sigp.eq.1).and.(dabs(y4col).gt.dabs(y3col)))then
        sigp=ecm
        nbin=int((sigp-sigpmin)/sigpw)+1
        if(nbin.ge.(ndiv_sigp+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
        if(dabs(yff).gt.yffcut)fxsigp(8,nbin,it)=fxsigp(8,nbin,it)+hist
        if(dabs(pzff).gt.pzffcut)fxsigp(9,nbin,it)=fxsigp(9,nbin,it)
     &                                                            +hist
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
        if(dabs(yff).gt.yffcut)fxsigm(8,nbin,it)=fxsigm(8,nbin,it)+hist
        if(dabs(pzff).gt.pzffcut)fxsigm(9,nbin,it)=fxsigm(9,nbin,it)
     &                                                            +hist
        end if
      end if