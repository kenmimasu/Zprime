      if(m_Y4col.eq.1)then
c generate distribution in Y4col.
        nbin=int((Y4col-Y4colmin)/Y4colw)+1
        if(nbin.ge.(ndiv_Y4col+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxY4col(nbin,it)=fxY4col(nbin,it)+hist
        end if
      end if
