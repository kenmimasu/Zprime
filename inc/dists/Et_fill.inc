      if(m_Et.eq.1)then
c generate distribution in Et.
        nbin=int((Et-Etmin)/Etw)+1
        if(nbin.ge.(ndiv_Et+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxEt(nbin,it)=fxEt(nbin,it)+hist
        end if
      end if
