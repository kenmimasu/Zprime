      if(m_delY.eq.1)then
c generate distribution in delY.
        nbin=int((delY-delYmin)/delYw)+1
        if(nbin.ge.(ndiv_delY+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxdelY(nbin,it)=fxdelY(nbin,it)+hist
        end if
      end if
