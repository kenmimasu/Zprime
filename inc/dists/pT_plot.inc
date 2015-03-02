C fermion tranverse momentum
      if(m_pT.eq.1)then
c plot distribution in pT.
        do j=1,ndiv_pT
          do i=1,it
            fxpT(j,i)=fxpT(j,i)*avgi/cnorm(i)/pTw
          end do
        end do
        do j=1,ndiv_pT
          fxpTtot(j)=0.d0 
          do i=1,it
            fxpTtot(j)=fxpTtot(j)+fxpT(j,i)
          end do
        end do
        write(jh,*)'DISTRIBUTION:'
        write(jh,*)'DSTRING: "pT"'
        write(jh,*)'TITLE: "p_T^',f,' Distribution"'
        write(jh,*)'AXES: "p_T (GeV)"',
     &' "d{/Symbol s}/p_T (pb/GeV)"'
        do i=1,ndiv_pT
          write(jh,964)xpT(i),fxpTtot(i)
        end do
        write(jh,*)'END'
      end if
