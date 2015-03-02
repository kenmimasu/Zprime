C anti-fermion CM frame pseudorapidity
      if(m_eta3.eq.1)then
c eta3 bin size.
        eta3w=(eta3max-eta3min)/ndiv_eta3
c generate bin in eta3.
        do i=1,ndiv_eta3
          xeta3(i)=eta3min+eta3w*(i-1)+eta3w/2.d0
        end do
      end if
