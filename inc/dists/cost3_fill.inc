      if(m_cost3.eq.1)then
c generate distribution in cost3.
        nbin=int((cost3-cost3min)/cost3w)+1
        if(nbin.ge.(ndiv_cost3+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxcost3(nbin,it)=fxcost3(nbin,it)+hist1
        end if
      nbin=int((cost3alt-cost3min)/cost3w)+1
      if(nbin.ge.(ndiv_cost3+1))then
        continue
      else if(nbin.lt.1)then
        continue
      else
        fxcost3(nbin,it)=fxcost3(nbin,it)+hist2
      end if
      end if
