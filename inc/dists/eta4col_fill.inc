      if(m_eta4col.eq.1)then
c generate distribution in eta4col.
        nbin=int((eta4col-eta4colmin)/eta4colw)+1
        if(nbin.ge.(ndiv_eta4col+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxeta4col(nbin,it)=fxeta4col(nbin,it)+hist
        end if
      end if
