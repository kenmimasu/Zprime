      if(m_rmass.eq.1)then
c generate distribution in rmass.
        nbin=int((rmass-rmassmin)/rmassw)+1
        if(nbin.ge.(ndiv_rmass+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxrmass(nbin,it)=fxrmass(nbin,it)+hist
        end if
      end if
