      if(m_Y3col.eq.1)then
c generate distribution in Y3col.
        nbin=int((Y3col-Y3colmin)/Y3colw)+1
        if(nbin.ge.(ndiv_Y3col+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxY3col(nbin,it)=fxY3col(nbin,it)+hist
        end if
      end if
