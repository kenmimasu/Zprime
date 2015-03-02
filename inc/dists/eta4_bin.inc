C fermion CM frame pseudorapidity
      if(m_eta4.eq.1)then
c eta4 bin size.
        eta4w=(eta4max-eta4min)/ndiv_eta4
c generate bin in eta4.
        do i=1,ndiv_eta4
          xeta4(i)=eta4min+eta4w*(i-1)+eta4w/2.d0
        end do
      end if
