      if(m_cost3col.eq.1)then
c generate distribution in cost3col.
        nbin=int((cost3col-cost3colmin)/cost3colw)+1
        if(nbin.ge.(ndiv_cost3col+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxcost3col(nbin,it)=fxcost3col(nbin,it)+hist
        end if
      end if
