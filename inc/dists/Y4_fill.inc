      if(m_Y4.eq.1)then
c generate distribution in Y4.
        nbin=int((Y4-Y4min)/Y4w)+1
        if(nbin.ge.(ndiv_Y4+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxY4(nbin,it)=fxY4(nbin,it)+hist
        end if
      end if
