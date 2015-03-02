
      if((m_sigp.eq.1).and.(dabs(y4col).le.yc))then
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
      if((m_sigp.eq.1).and.(dabs(y4col).gt.yc))then
        sigp=ecm
        nbin=int((sigp-sigpmin)/sigpw)+1
      end if
        
      if((m_sigm.eq.1).and.(dabs(y3col).le.yc))then
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
      if((m_sigm.eq.1).and.(dabs(y3col).gt.yc))then
        sigm=ecm
        nbin=int((sigm-sigmmin)/sigmw)+1
      end if