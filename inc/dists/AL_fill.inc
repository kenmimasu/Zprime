      if(m_sigp.eq.1)then
c generate distribution in sigp for AL.
        sigp=ecm
        nbin=int((sigp-sigpmin)/sigpw)+1
        if(nbin.ge.(ndiv_sigp+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxsigp(iAL,nbin,it)=fxsigp(iAL,nbin,it)+
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
          fxsigm(iAL,nbin,it)=fxsigm(iAL,nbin,it)+
     &    (weight(it,+1,-1)+weight(it,+1,+1))
        end if
      end if  