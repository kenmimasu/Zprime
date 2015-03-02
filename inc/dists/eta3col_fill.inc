      if(m_eta3col.eq.1)then
c generate distribution in eta3col.
        nbin=int((eta3col-eta3colmin)/eta3colw)+1
        if(nbin.ge.(ndiv_eta3col+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxeta3col(nbin,it)=fxeta3col(nbin,it)+hist
        end if
      end if
