c generate distribution in sigp for AFBSTAR.
      if((m_sigp.eq.1).and.(costst.gt.0.d0))then
        sigp=ecm
        nbin=int((sigp-sigpmin)/sigpw)+1
        if(nbin.ge.(ndiv_sigp+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
            fxsigp(10,nbin,it)=fxsigp(10,nbin,it)+hist 
            do j=1,8
              if(dabs(yff).gt.yffs(j))then
                fxsigp(10+j,nbin,it)=fxsigp(10+j,nbin,it)+hist
              endif
            enddo 
        end if
      end if
c generate distribution in sigm for AFBSTAR.
      if((m_sigm.eq.1).and.(costst.lt.0.d0))then
        sigm=ecm
        nbin=int((sigm-sigmmin)/sigmw)+1
        if(nbin.ge.(ndiv_sigm+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
            fxsigm(10,nbin,it)=fxsigm(10,nbin,it)+hist 
            do j=1,8
              if(dabs(yff).gt.yffs(j))then
                fxsigm(10+j,nbin,it)=fxsigm(10+j,nbin,it)+hist
              endif
            enddo
        end if
      end if