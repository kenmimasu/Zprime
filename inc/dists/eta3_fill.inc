      if(m_eta3.eq.1)then
c generate distribution in eta3.
        nbin=int((eta3-eta3min)/eta3w)+1
        if(nbin.ge.(ndiv_eta3+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxeta3(nbin,it)=fxeta3(nbin,it)+hist1
        end if
      nbin=int((eta3alt-eta3min)/eta3w)+1
      if(nbin.ge.(ndiv_eta3+1))then
        continue
      else if(nbin.lt.1)then
        continue
      else
        fxeta3(nbin,it)=fxeta3(nbin,it)+hist2
      end if
      end if
