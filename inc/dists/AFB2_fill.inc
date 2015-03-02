      if((m_sigp.eq.1).and.(cost4.gt.0.d0))then
c generate distribution in sigp for AFB2.
        sigp=ecm
        nbin=int((sigp-sigpmin)/sigpw)+1
        if(nbin.ge.(ndiv_sigp+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxsigp(7,nbin,it)=fxsigp(7,nbin,it)+hist1
          fxsigm(7,nbin,it)=fxsigm(7,nbin,it)+hist2
        end if
      end if  
      if((m_sigm.eq.1).and.(cost4.lt.0.d0))then
c generate distribution in sigm for AFB2.
        sigm=ecm
        nbin=int((sigm-sigmmin)/sigmw)+1
        if(nbin.ge.(ndiv_sigm+1))then
          continue
        else if(nbin.lt.1)then
          continue
        else
          fxsigm(7,nbin,it)=fxsigm(7,nbin,it)+hist1
          fxsigp(7,nbin,it)=fxsigp(7,nbin,it)+hist2
        end if
      end if  