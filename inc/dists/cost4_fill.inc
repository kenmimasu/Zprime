      if(m_cost4.eq.1)then
c generate distribution in cost4.
        nbin=int((cost4-cost4min)/cost4w)+1
        if(nbin.ge.(ndiv_cost4+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxcost4(nbin,it)=fxcost4(nbin,it)+hist
        end if
      end if
