C fermion tranverse momentum
      if(m_pT.eq.1)then
c pT bin size.
        pTw=(pTmax-pTmin)/ndiv_pT
c generate bin in pT.
        do i=1,ndiv_pT
          xpT(i)=pTmin+pTw*(i-1)+pTw/2.d0
          do it=1,20
            fxpT(i,it)=0d0
          end do
        end do
      end if
