      if((m_sigp.eq.1).and.(cost4col.gt.0.d0))then
c generate distribution in sigp for AFB.
        sigp=ecm
        nbin=int((sigp-sigpmin)/sigpw)+1
        if(nbin.ge.(ndiv_sigp+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxsigp(4,nbin,it)=fxsigp(4,nbin,it)+hist
        end if
      end if  
      if((m_sigm.eq.1).and.(cost4col.lt.0.d0))then
c generate distribution in sigm for AFB.
        sigm=ecm
        nbin=int((sigm-sigmmin)/sigmw)+1
        if(nbin.ge.(ndiv_sigm+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxsigm(4,nbin,it)=fxsigm(4,nbin,it)+hist
        end if
      end if  