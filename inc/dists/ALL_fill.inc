      if(m_sigp.eq.1)then
c generate distribution in sigp for ALL.
        sigp=ecm
        nbin=int((sigp-sigpmin)/sigpw)+1
        if(nbin.ge.(ndiv_sigp+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxsigp(iALL,nbin,it)=fxsigp(iALL,nbin,it)+
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
          fxsigm(iALL,nbin,it)=fxsigm(iALL,nbin,it)+
     &    (weight(it,+1,-1)+weight(it,-1,+1))
        end if
      end if  