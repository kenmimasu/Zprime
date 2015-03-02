C fermion CM frame rapidity
      if(m_Y4.eq.1)then
c Y4 bin size.
        Y4w=(Y4max-Y4min)/ndiv_Y4
c generate bin in Y4.
        do i=1,ndiv_Y4
          xY4(i)=Y4min+Y4w*(i-1)+Y4w/2.d0
        end do
      end if
