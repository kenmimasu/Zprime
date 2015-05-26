c generate distribution in sigp for AFBSTAR.
      if((m_sigp.eq.1).and.(costst.gt.0.d0))then
        sigp=ecm
        nbin=int((sigp-sigpmin)/sigpw)+1
        if(nbin.ge.(ndiv_sigp+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
            fxsigp(iAFBM,nbin,it)=fxsigp(iAFBM,nbin,it)+hist 
            do j=1,3
                yecm =
              if(dabs(yff).gt.yffs(j))then
                fxsigp(iAFBM+j,nbin,it)=fxsigp(iAFBM+j,nbin,it)+hist
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
            fxsigm(iAFBM,nbin,it)=fxsigm(iAFBM,nbin,it)+hist 
            do j=1,3
              if(dabs(yff).gt.yffs(j))then
                fxsigm(iAFBM+j,nbin,it)=fxsigm(iAFBM+j,nbin,it)+hist
              endif
            enddo
        end if
      end if