      if(m_eta4.eq.1)then
c generate distribution in eta4.
        nbin=int((eta4-eta4min)/eta4w)+1
        if(nbin.ge.(ndiv_eta4+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxeta4(nbin,it)=fxeta4(nbin,it)+hist1
        end if
      nbin=int((eta4alt-eta4min)/eta4w)+1
      if(nbin.ge.(ndiv_eta4+1))then
        continue
      else if(nbin.lt.1)then
        continue
      else
        fxeta4(nbin,it)=fxeta4(nbin,it)+hist2
      end if
      end if
