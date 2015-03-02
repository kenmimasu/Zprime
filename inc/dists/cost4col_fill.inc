      if(m_cost4col.eq.1)then
c generate distribution in cost4col.
        nbin=int((cost4col-cost4colmin)/cost4colw)+1
        if(nbin.ge.(ndiv_cost4col+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxcost4col(nbin,it)=fxcost4col(nbin,it)+hist
        end if
      end if
