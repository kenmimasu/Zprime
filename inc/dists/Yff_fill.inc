      if(m_Yff.eq.1)then
c generate distribution in Yff.
        nbin=int((Yff-Yffmin)/Yffw)+1
        if(nbin.ge.(ndiv_Yff+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxYff(nbin,it)=fxYff(nbin,it)+hist
        end if
      end if
