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
