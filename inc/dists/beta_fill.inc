      if(m_beta.eq.1)then
c generate distribution in beta.
        nbin=int((beta-betamin)/betaw)+1
        if(nbin.ge.(ndiv_beta+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxbeta(nbin,it)=fxbeta(nbin,it)+hist
        end if
      end if
