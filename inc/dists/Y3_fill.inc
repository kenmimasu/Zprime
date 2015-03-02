      if(m_Y3.eq.1)then
c generate distribution in Y3.
        nbin=int((Y3-Y3min)/Y3w)+1
        if(nbin.ge.(ndiv_Y3+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxY3(nbin,it)=fxY3(nbin,it)+hist
        end if
      end if
