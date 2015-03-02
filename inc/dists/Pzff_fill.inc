      if(m_Pzff.eq.1)then
c generate distribution in Pzff.
        nbin=int((Pzff-Pzffmin)/Pzffw)+1
        if(nbin.ge.(ndiv_Pzff+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxPzff(nbin,it)=fxPzff(nbin,it)+hist
        end if
      end if
