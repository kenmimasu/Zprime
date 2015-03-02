C anti-fermion CM frame rapidity
      if(m_Y3.eq.1)then
c Y3 bin size.
        Y3w=(Y3max-Y3min)/ndiv_Y3
c generate bin in Y3.
        do i=1,ndiv_Y3
          xY3(i)=Y3min+Y3w*(i-1)+Y3w/2.d0
        end do
      end if
