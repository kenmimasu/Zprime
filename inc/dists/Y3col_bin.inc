C anti-fermion lab frame rapidity
      if(m_Y3col.eq.1)then
c Y3col bin size.
        Y3colw=(Y3colmax-Y3colmin)/ndiv_Y3col
c generate bin in Y3col.
        do i=1,ndiv_Y3col
          xY3col(i)=Y3colmin+Y3colw*(i-1)+Y3colw/2.d0
        end do
      end if
